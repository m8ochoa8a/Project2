#!/usr/bin/env python3
"""Scan the full AAV9 capsid against virion surface proteins of 50 human viruses.

The query is AAV9 VP1 (736 residues), which includes the VP2 and VP3 regions.
The targets are proteins on the outside of the infectious particle: the outer
capsid of non-enveloped viruses, and the envelope glycoproteins of enveloped
viruses. Internal capsids, tegument, replication proteins, and infected-cell
membrane proteins are left out. Sequence identity is a screen, not evidence
that an antibody or T cell would cross-react.
"""

from __future__ import annotations

import csv
import json
import random
import re
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

from homology import (
    AA,
    AA_INDEX,
    NINE,
    NULL_REPLICATES,
    NULL_SEED,
    RESULTS,
    alignment_stats,
    classify,
    clean_protein,
    composition_shuffle,
    encode_kmer,
    epitope_in_alignment,
    load_aav9,
    load_epitopes,
    load_viruses,
    nine_mer_hit_count,
    nine_mer_variants,
    read_fasta,
    to_array,
    write_tsv,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PROT_DIR = DATA / "proteomes"
FEATURE_TABLE = DATA / "features" / "polyproteins.ft"
VACCINIA_FAA = DATA / "reference" / "vaccinia_copenhagen.faa"
WINDOW_LENGTHS = (9, 11, 15)
MCV_MIN_IDENTITY = 0.30
MCV_MIN_COVERAGE = 0.50

# Copenhagen vaccinia (NC_001559.1) virion-membrane proteins. Protein ids are the
# RefSeq CDS whose prot_desc carries that gene name.
VACCINIA_SURFACE = {
    "A27L": "NP_063815.1",
    "D8L": "NP_063770.1",
    "H3L": "NP_063753.1",
    "L1R": "NP_063739.1",
    "A26L": "NP_063814.1",
    "A17L": "NP_063802.1",
    "A14L": "NP_063799.1",
    "F13L": "NP_063694.1",
    "B5R": "NP_063864.1",
    "A33R": "NP_063823.1",
    "A34R": "NP_063824.1",
    "A56R": "NP_063854.1",
    "A28L": "NP_063816.1",
    "F9L": "NP_063688.1",
    "H2R": "NP_063752.1",
    "G9R": "NP_063738.1",
    "J5L": "NP_063748.1",
    "A21L": "NP_063805.1",
    "L5R": "NP_063743.1",
    "A16L": "NP_063801.1",
    "G3L": "NP_063729.1",
}


def penalized_blosum62():
    """BLOSUM62 with X treated as a mismatch, including X versus X."""
    base = substitution_matrices.load("BLOSUM62")
    matrix = substitution_matrices.Array(alphabet=base.alphabet, dims=2)
    for left in base.alphabet:
        for right in base.alphabet:
            matrix[left, right] = base[left, right]
    for residue in base.alphabet:
        matrix[residue, "X"] = -4
        matrix["X", residue] = -4
    return matrix


def for_alignment(seq: str) -> str:
    return "".join(char if char in AA else "X" for char in seq)


def load_protein_index() -> dict[str, dict]:
    index = {}
    for path in sorted(PROT_DIR.glob("*.faa")):
        virus_id = path.stem
        for header, seq in read_fasta(path):
            match = re.search(r"\[protein_id=([^\]]+)\]", header)
            product_match = re.search(r"\[protein=([^\]]+)\]", header)
            if not match:
                continue
            protein_id = match.group(1)
            index[protein_id] = {
                "virus_id": virus_id,
                "protein_id": protein_id,
                "product": product_match.group(1) if product_match else protein_id,
                "seq": clean_protein(seq),
            }
    return index


def parse_feature_table(path: Path) -> dict[str, list[dict]]:
    records: dict[str, list[dict]] = {}
    accession = None
    features: list[dict] = []
    current = None
    for line in path.read_text().splitlines():
        if line.startswith(">Feature"):
            if accession:
                records[accession] = features
            parts = line.split("|")
            accession = parts[1] if len(parts) >= 2 else line.split()[-1]
            features = []
            current = None
            continue
        if line.startswith("\t\t\t"):
            qual = line.strip().split("\t")
            if current is not None and len(qual) >= 2:
                current["quals"].setdefault(qual[0], []).append(qual[1])
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit():
            current = {
                "start": int(parts[0]),
                "end": int(parts[1]),
                "type": parts[2],
                "quals": {},
            }
            if parts[2] in {"mat_peptide", "Region"}:
                features.append(current)
            else:
                current = None
    if accession:
        records[accession] = features
    return records


def picorna_chain(products: list[str]) -> str | None:
    text = " ".join(products)
    for name in ("VP0", "VP1", "VP2", "VP3"):
        if re.search(rf"\b{name}\b", text):
            return name
    return None


def slice_entry(
    index: dict[str, dict],
    virus_id: str,
    protein_id: str,
    start: int,
    end: int,
    product: str,
    compartment: str,
    note: str,
) -> dict:
    record = index[protein_id]
    if record["virus_id"] != virus_id:
        raise SystemExit(f"{protein_id} belongs to {record['virus_id']}, not {virus_id}")
    if start < 1 or end > len(record["seq"]) or start > end:
        raise SystemExit(
            f"{protein_id} slice {start}-{end} is outside a protein of {len(record['seq'])} aa"
        )
    seq = record["seq"][start - 1 : end]
    return {
        "virus_id": virus_id,
        "role": "panel",
        "protein_id": protein_id,
        "product": product,
        "compartment": compartment,
        "start": start,
        "end": end,
        "length": len(seq),
        "selection": note,
        "seq": seq,
        "label": f"{product} ({protein_id}:{start}-{end})",
    }


def full_entry(
    index: dict[str, dict],
    virus_id: str,
    protein_id: str,
    compartment: str,
    note: str,
    product: str | None = None,
    role: str = "panel",
) -> dict:
    record = index[protein_id]
    if record["virus_id"] != virus_id:
        raise SystemExit(f"{protein_id} belongs to {record['virus_id']}, not {virus_id}")
    name = product or record["product"]
    seq = record["seq"]
    return {
        "virus_id": virus_id,
        "role": role,
        "protein_id": protein_id,
        "product": name,
        "compartment": compartment,
        "start": 1,
        "end": len(seq),
        "length": len(seq),
        "selection": note,
        "seq": seq,
        "label": f"{name} ({protein_id})",
    }


def mature_chain_entries(index: dict[str, dict], features: dict[str, list[dict]]) -> list[dict]:
    """Outer-capsid chains cut from polyproteins at annotated mature-peptide bounds."""
    picorna = {
        "HRV_A": "NP_042288.1",
        "HRV_B": "NP_041009.1",
        "HRV_C": "YP_001552411.1",
        "EVA": "NP_042242.1",
        "EVB": "NP_040958.1",
        "EVD68": "YP_009505605.1",
        "HPEV": "NP_046804.1",
    }
    entries = []
    for virus_id, protein_id in picorna.items():
        found = []
        for feature in features.get(protein_id, []):
            if feature["type"] != "mat_peptide" or feature["start"] > feature["end"]:
                continue
            name = picorna_chain(feature["quals"].get("product", []))
            if name is None:
                continue
            found.append((name, feature["start"], feature["end"]))
        names = [item[0] for item in found]
        expected = ["VP0", "VP3", "VP1"] if virus_id == "HPEV" else ["VP2", "VP3", "VP1"]
        if names != expected:
            raise SystemExit(f"{virus_id} mature chains parsed as {names}, expected {expected}")
        for name, start, end in found:
            compartment = "outer capsid" if name != "VP0" else "outer capsid (VP0; parechoviruses do not cleave a VP4)"
            entries.append(
                slice_entry(
                    index,
                    virus_id,
                    protein_id,
                    start,
                    end,
                    f"capsid protein {name}",
                    compartment,
                    "NCBI protein feature table mat_peptide; VP4 is internal and was not taken",
                )
            )

    dengue = []
    for feature in features.get("NP_059433.1", []):
        if feature["type"] != "mat_peptide":
            continue
        product = " | ".join(feature["quals"].get("product", []))
        if product == "membrane glycoprotein M":
            dengue.append(("M", "membrane protein M", feature["start"], feature["end"]))
        elif product == "envelope protein E":
            dengue.append(("E", "envelope protein E", feature["start"], feature["end"]))
    if [item[0] for item in dengue] != ["M", "E"]:
        raise SystemExit(f"dengue mature surface chains parsed as {dengue}")
    for _key, product, start, end in dengue:
        entries.append(
            slice_entry(
                index,
                "DENV1",
                "NP_059433.1",
                start,
                end,
                product,
                "virion envelope",
                "NCBI mat_peptide on the dengue polyprotein; capsid C, pr, and nonstructural proteins were not taken",
            )
        )
    return entries


def hav_entries(index: dict[str, dict]) -> list[dict]:
    """HAV RefSeq feature table only annotates an alternate VP4. Chains are UniProt P08617."""
    parent = index["NP_041007.1"]
    if parent["virus_id"] != "HAV" or len(parent["seq"]) != 2227 or not parent["seq"].startswith("MNMSRQGIFQTVGSGLDHIL"):
        raise SystemExit("HAV NP_041007.1 does not match the UniProt P08617 polyprotein")
    chains = (
        ("VP2", 24, 245, "outer capsid"),
        ("VP3", 246, 491, "outer capsid"),
        ("VP1", 492, 764, "outer capsid"),
    )
    return [
        slice_entry(
            index,
            "HAV",
            "NP_041007.1",
            start,
            end,
            f"capsid protein {name}",
            compartment,
            "UniProt P08617 chains; the sequence is identical to NP_041007.1. VP4 (1-23) is internal. VP1 stops at 764, before the 2A tail",
        )
        for name, start, end, compartment in chains
    ]


def sapovirus_entry(index: dict[str, dict], features: dict[str, list[dict]]) -> dict:
    coat = [
        feature
        for feature in features.get("YP_077278.1", [])
        if feature["type"] == "Region" and "Calici_coat" in feature["quals"].get("region", [])
    ]
    if len(coat) != 1 or (coat[0]["start"], coat[0]["end"]) != (1738, 2008):
        raise SystemExit("sapovirus Calici_coat domain was not at residues 1738-2008")
    parent = index["YP_077278.1"]
    if len(parent["seq"]) != 2280:
        raise SystemExit(f"sapovirus ORF1 length is {len(parent['seq'])}, expected 2280")
    return slice_entry(
        index,
        "SAPV",
        "YP_077278.1",
        1738,
        2280,
        "capsid protein VP1",
        "outer capsid",
        "ORF1 has no mat_peptide. VP1 is taken from the Calici_coat domain (1738) through the ORF1 C-terminus. ORF2, the minor capsid protein, is internal to the particle and was left out",
    )


def explicit_entries(index: dict[str, dict]) -> list[dict]:
    """Proteins chosen by id because product-name filters pull in non-virion membrane proteins."""
    rows: list[tuple] = []

    def add(virus_id, protein_id, compartment, note, product=None, role="panel"):
        rows.append((virus_id, protein_id, compartment, note, product, role))

    envelope = "virion envelope glycoprotein"
    # EBV. BMRF2 is an UL43-like envelope protein and is not part of the virion surface set.
    for protein_id, product in (
        ("YP_401713.1", "envelope glycoprotein B"),
        ("YP_401700.1", "envelope glycoprotein H"),
        ("YP_401678.1", "envelope glycoprotein L"),
        ("YP_401672.1", "envelope glycoprotein gp42"),
        ("YP_401665.1", "envelope glycoprotein N"),
        ("YP_401685.1", "envelope glycoprotein M"),
        ("YP_401694.1", "envelope glycoprotein gp150"),
        ("YP_401695.1", "envelope glycoprotein gp48"),
        ("YP_401667.1", "glycoprotein 350"),
    ):
        add("EBV", protein_id, envelope, "EBV virion envelope glycoprotein", product)

    # CMV ids are the conserved entry and envelope glycoproteins. US2-US11, RL10, UL4,
    # UL33, and UL37 are infected-cell or tegument proteins and are not included.
    for protein_id, product in (
        ("YP_081514.1", "envelope glycoprotein B"),
        ("YP_081523.1", "envelope glycoprotein H"),
        ("YP_081555.1", "envelope glycoprotein L"),
        ("YP_081522.1", "envelope glycoprotein O"),
        ("YP_081547.1", "envelope glycoprotein M"),
        ("YP_081521.1", "envelope glycoprotein N"),
        ("YP_081565.1", "envelope glycoprotein UL130"),
        ("YP_081566.1", "envelope protein UL131A"),
        ("YP_081567.1", "envelope glycoprotein UL132"),
    ):
        add("CMV", protein_id, envelope, "CMV virion envelope glycoprotein (gB/gH/gL/gO/gM/gN or the pentamer subunits UL130/UL131A/UL132)", product)

    # HHV-6B. U46 gN is the 84-aa ORF immediately upstream of gO; RefSeq does not name it gN.
    # U20-U23 immunomodulatory glycoproteins are excluded.
    for protein_id, product, note in (
        ("NP_050220.1", "glycoprotein B", "HHV-6B U39 gB"),
        ("NP_050229.1", "glycoprotein H", "HHV-6B U48 gH"),
        ("NP_050261.1", "glycoprotein L", "HHV-6B U82 gL"),
        ("NP_050251.1", "glycoprotein M", "HHV-6B U72 gM"),
        ("NP_050227.1", "glycoprotein N", "HHV-6B U46 gN; RefSeq product is 'putative membrane /secreted protein', the 84-aa ORF directly upstream of gO"),
        ("NP_050228.1", "glycoprotein O", "HHV-6B U47 gO"),
        ("NP_050271.1", "envelope glycoprotein Q", "HHV-6B U100 gQ, the spliced envelope glycoprotein"),
        ("NP_050230.1", "putative fusion protein U49", "HHV-6B ORF immediately downstream of gH; RefSeq calls it a putative fusion protein"),
        ("NP_050264.1", "glycoprotein U85", "HHV-6B glycoprotein annotated next to U84; included as an annotated glycoprotein and not as a core entry protein"),
    ):
        add("HHV6B", protein_id, envelope, note, product)

    for protein_id, product, note in (
        ("YP_073779.1", "envelope glycoprotein B", "HHV-7 U39 gB"),
        ("YP_073788.1", "envelope glycoprotein H", "HHV-7 U48 gH"),
        ("YP_073820.1", "envelope glycoprotein L", "HHV-7 U82 gL"),
        ("YP_073812.1", "envelope glycoprotein M", "HHV-7 U72 gM"),
        ("YP_073786.1", "envelope glycoprotein N", "HHV-7 U46 gN"),
        ("YP_073787.1", "envelope glycoprotein O", "HHV-7 U47 gO"),
        ("YP_073827.1", "envelope glycoprotein Q", "HHV-7 U100 gQ"),
    ):
        add("HHV7", protein_id, envelope, note, product)

    for virus_id, protein_id in (
        ("BKV", "YP_717939.1"),
        ("JCV", "NP_043511.1"),
        ("MCPyV", "YP_009111420.1"),
        ("KIPYV", "YP_001111258.1"),
        ("WUPYV", "YP_001285487.1"),
    ):
        note = "Polyomavirus surface is the major capsid VP1. VP2 and VP3 are interior."
        product = "major capsid protein VP1"
        if virus_id == "JCV":
            note += " The FASTA product is 'hypothetical protein'; this accession is the VP1-length capsid protein."
        add(virus_id, protein_id, "outer capsid", note, product)

    for virus_id, pairs in (
        ("B19", (("YP_004928146.1", "capsid protein VP1"), ("YP_004928148.1", "capsid protein VP2"))),
        ("HBOV1", (("YP_338088.1", "capsid protein VP1"), ("YP_338089.1", "capsid protein VP2"))),
    ):
        for protein_id, product in pairs:
            add(virus_id, protein_id, "outer capsid", "Parvovirus capsid surface. NS1 and NP-1 are not capsid proteins.", product)
    add("TTV1", "NP_817122.1", "outer capsid", "Anellovirus surface is ORF1 VP1. VP2 and VP3 were not taken.", "capsid protein VP1")

    for virus_id, ids in (
        ("HADV_C", (("NP_040525.1", "hexon"), ("NP_040521.1", "penton base"), ("NP_040533.1", "fiber"), ("NP_040514.1", "capsid protein IX"))),
        ("HADV_B", (("YP_002213779.1", "hexon"), ("YP_002213774.1", "penton base"), ("YP_002213796.1", "fiber"), ("YP_002213768.1", "capsid protein IX"))),
    ):
        for protein_id, product in ids:
            add(virus_id, protein_id, "outer capsid", "Adenovirus exterior capsid. pIIIa, pVI, pVIII, and E3 membrane proteins were not taken.", product)

    for virus_id, ids, note in (
        ("FLUA_H1N1", (("YP_009118626.1", "hemagglutinin"), ("YP_009118627.1", "neuraminidase"), ("YP_009118622.1", "matrix protein 2")), "Influenza A virion surface. NP, M1, and the polymerase were not taken."),
        ("FLUA_H3N2", (("YP_308839.1", "hemagglutinin"), ("YP_308842.1", "neuraminidase"), ("YP_308840.1", "matrix protein 2")), "Influenza A virion surface. NP, M1, and the polymerase were not taken."),
        ("FLUB", (("NP_056660.1", "hemagglutinin"), ("NP_056663.1", "neuraminidase"), ("NP_056662.1", "NB protein")), "Influenza B virion surface. BM2 is not annotated on this RefSeq. NP and M1 were not taken."),
    ):
        for protein_id, product in ids:
            compartment = "envelope small protein" if product in {"matrix protein 2", "NB protein"} else "virion envelope glycoprotein"
            add(virus_id, protein_id, compartment, note, product)

    for virus_id, ids in (
        ("RSV_A", (("NP_044596.1", "fusion protein F"), ("NP_044595.1", "attachment glycoprotein G"), ("NP_044594.1", "small hydrophobic protein SH"))),
        ("HMPV", (("YP_012608.1", "fusion protein F"), ("YP_012612.1", "attachment glycoprotein G"), ("YP_012611.1", "small hydrophobic protein SH"))),
    ):
        for protein_id, product in ids:
            compartment = "envelope small protein" if product.endswith("SH") else "virion envelope glycoprotein"
            add(virus_id, protein_id, compartment, "Pneumovirus envelope. Matrix and polymerase proteins were not taken.", product)
    for virus_id, fusion_id, hn_id in (
        ("HPIV1", "NP_604440.1", "NP_604441.1"),
        ("HPIV2", "NP_598404.1", "NP_598405.1"),
        ("HPIV3", "NP_067151.1", "NP_067152.1"),
        ("HPIV4", "YP_008378663.1", "YP_008378664.1"),
    ):
        add(virus_id, fusion_id, "virion envelope glycoprotein", "Parainfluenza envelope. The nucleocapsid, matrix, and polymerase were not taken.", "fusion protein F")
        add(virus_id, hn_id, "virion envelope glycoprotein", "Parainfluenza envelope. The nucleocapsid, matrix, and polymerase were not taken.", "hemagglutinin-neuraminidase")

    for virus_id, ids in (
        ("COV_229E", (("NP_073551.1", "surface glycoprotein"), ("NP_073554.1", "envelope protein"), ("NP_073555.1", "membrane protein"))),
        ("COV_NL63", (("YP_003767.1", "spike protein"), ("YP_003769.1", "envelope protein"), ("YP_003770.1", "membrane protein"))),
        ("COV_OC43", (("YP_009555241.1", "spike glycoprotein"), ("YP_009555240.1", "hemagglutinin-esterase"), ("YP_009555243.1", "envelope protein"), ("YP_009555244.1", "membrane protein"))),
        ("COV_HKU1", (("YP_173238.1", "spike glycoprotein"), ("YP_173237.1", "hemagglutinin-esterase"), ("YP_173240.1", "envelope protein"), ("YP_173241.1", "membrane glycoprotein"))),
        ("SARS2", (("YP_009724390.1", "surface glycoprotein"), ("YP_009724392.1", "envelope protein"), ("YP_009724393.1", "membrane glycoprotein"))),
    ):
        for protein_id, product in ids:
            compartment = "virion envelope glycoprotein" if "envelope protein" not in product and not product.startswith("membrane") else "envelope small protein"
            if "hemagglutinin" in product or "spike" in product or "surface" in product:
                compartment = "virion envelope glycoprotein"
            add(virus_id, protein_id, compartment, "Coronavirus virion envelope. Nucleocapsid and the replicase were not taken.", product)

    add("NORO_GI", "NP_056821.2", "outer capsid", "Norovirus surface is VP1. VP2 is a minor interior protein and was not taken.", "capsid protein VP1")
    add("NORO_GII", "YP_009518842.1", "outer capsid", "Norovirus surface is VP1. VP2 is a minor interior protein and was not taken.", "capsid protein VP1")
    add("ROTA", "YP_002302230.1", "outer capsid", "Rotavirus outer capsid only. VP6 and the inner capsid were not taken.", "outer capsid spike protein VP4")
    add("ROTA", "YP_002302222.1", "outer capsid", "Rotavirus outer capsid only. VP6 and the inner capsid were not taken.", "outer capsid protein VP7")
    add("ASTV", "NP_059444.1", "outer capsid", "Astrovirus capsid precursor. The nonstructural protein was not taken.", "capsid precursor")
    add("HEV", "NP_056788.1", "outer capsid", "Hepatitis E capsid. ORF1 is the replicase and was not taken.", "capsid protein")
    add("HBV", "YP_009173869.1", "virion envelope", "Large envelope protein only. It contains preS1, preS2, and S. Middle and small envelope proteins are truncations of this sequence, so they were not scanned again. Core is internal.", "large envelope protein")
    add("HPV16", "NP_041332.2", "outer capsid", "Papillomavirus major capsid.", "major capsid protein L1")
    add("HPV16", "NP_041331.2", "minor capsid", "L2 is partly internal. A portion is exposed on the particle, so it is included and labeled as the minor capsid protein.", "minor capsid protein L2")
    add("HPV6", "NP_040304.1", "outer capsid", "Papillomavirus major capsid.", "major capsid protein L1")
    add("HPV6", "NP_040303.1", "minor capsid", "L2 is partly internal. A portion is exposed on the particle, so it is included and labeled as the minor capsid protein.", "minor capsid protein L2")
    add("AAV2", "YP_680426.1", "outer capsid", "Positive control. AAV2 VP1 only; VP2 and VP3 are truncations of the same sequence. Not one of the 50.", "major coat protein VP1", "control")

    return [full_entry(index, *row) for row in rows]


def named_envelope_entries(index: dict[str, dict]) -> list[dict]:
    """HSV-1, HSV-2, and VZV name every virion envelope glycoprotein as such."""
    entries = []
    excluded_bits = ("UL20", "UL43", "UL45", "UL56", "US8A", "US9")
    for virus_id in ("HSV1", "HSV2", "VZV"):
        matched = []
        for record in index.values():
            if record["virus_id"] != virus_id:
                continue
            if "envelope glycoprotein" not in record["product"].lower():
                continue
            if any(bit in record["product"] for bit in excluded_bits):
                continue
            matched.append(record)
        if virus_id == "VZV" and len(matched) != 9:
            raise SystemExit(f"VZV envelope glycoprotein count is {len(matched)}, expected 9")
        if virus_id in {"HSV1", "HSV2"} and len(matched) != 12:
            raise SystemExit(f"{virus_id} envelope glycoprotein count is {len(matched)}, expected 12")
        for record in matched:
            entries.append(
                full_entry(
                    index,
                    virus_id,
                    record["protein_id"],
                    "virion envelope glycoprotein",
                    "Product name contains 'envelope glycoprotein'. UL20, UL43, UL45, UL56, US8A, and US9 are not in this set.",
                )
            )
    return entries


def load_vaccinia() -> dict[str, str]:
    if not VACCINIA_FAA.exists():
        raise SystemExit(f"missing {VACCINIA_FAA}")
    found = {}
    for header, seq in read_fasta(VACCINIA_FAA):
        match = re.search(r"\[protein_id=([^\]]+)\]", header)
        if match:
            found[match.group(1)] = clean_protein(seq)
    missing = [protein_id for protein_id in VACCINIA_SURFACE.values() if protein_id not in found]
    if missing:
        raise SystemExit(f"vaccinia FASTA is missing {missing}")
    return found


def molluscum_entries(index: dict[str, dict], aligner: PairwiseAligner) -> tuple[list[dict], list[dict]]:
    vaccinia = load_vaccinia()
    queries = [(gene, vaccinia[protein_id], protein_id) for gene, protein_id in VACCINIA_SURFACE.items()]
    mcv = [record for record in index.values() if record["virus_id"] == "MCV"]
    audit = []
    entries = []
    print(f"aligning {len(mcv)} molluscum ORFs to {len(queries)} vaccinia membrane proteins", flush=True)
    for record in mcv:
        best = None
        passed = None
        target_seq = for_alignment(record["seq"])
        for gene, seq, protein_id in queries:
            alignment = next(aligner.align(target_seq, for_alignment(seq)))
            stats = alignment_stats(alignment)
            identity = stats["identity"]
            coverage = stats["aligned_columns"] / len(seq) if seq else 0.0
            hit = {
                "gene": gene,
                "protein_id": protein_id,
                "vaccinia_length": len(seq),
                "score": stats["score"],
                "aligned_columns": stats["aligned_columns"],
                "identity": identity,
                "coverage": coverage,
            }
            if best is None or hit["score"] > best["score"]:
                best = hit
            # 30% identity across half of the vaccinia protein is the main rule.
            # A slightly lower identity is kept when the alignment still covers
            # at least three quarters of the vaccinia protein and the score is
            # in the range of the full-length homologs. That recovers the H3L
            # alignment, which runs the length of both proteins at 29.7%.
            full_length = coverage >= 0.75 and identity >= 0.25 and hit["score"] >= 200
            ok = coverage >= MCV_MIN_COVERAGE and (identity >= MCV_MIN_IDENTITY or full_length)
            if ok and (passed is None or hit["score"] > passed["score"]):
                passed = hit
        if best is None:
            continue
        chosen = passed or best
        audit.append(
            {
                "mcv_protein_id": record["protein_id"],
                "mcv_product": record["product"],
                "mcv_length": len(record["seq"]),
                "vaccinia_gene": chosen["gene"],
                "vaccinia_protein_id": chosen["protein_id"],
                "vaccinia_length": chosen["vaccinia_length"],
                "score": round(chosen["score"], 1),
                "aligned_columns": chosen["aligned_columns"],
                "percent_identity": round(100 * chosen["identity"], 1),
                "coverage_of_vaccinia": round(100 * chosen["coverage"], 1),
                "accepted": "yes" if passed else "no",
                "note": "",
            }
        )
        if passed is None:
            continue
        entries.append(
            full_entry(
                index,
                "MCV",
                record["protein_id"],
                "virion membrane",
                (
                    f"Best BLOSUM62 local hit to vaccinia {passed['gene']} ({passed['protein_id']}): "
                    f"{passed['identity']:.1%} identity over {passed['aligned_columns']} columns "
                    f"({passed['coverage']:.0%} of {passed['gene']}). "
                    "Accepted at >=30% identity and >=50% coverage of the vaccinia protein."
                ),
                product=f"{record['product']} ({passed['gene']} homolog)",
            )
        )
    # A long unrelated ORF can clear 30% on a short slice of a gene that already
    # has a full-length homolog. Keep the stronger alignment.
    best_score = {}
    for row in audit:
        if row["accepted"] == "yes":
            gene = row["vaccinia_gene"]
            best_score[gene] = max(best_score.get(gene, 0), row["score"])
    drop_ids = set()
    for row in audit:
        if row["accepted"] != "yes":
            continue
        if row["score"] < 0.5 * best_score[row["vaccinia_gene"]]:
            row["accepted"] = "no"
            row["note"] = (
                f"partial overlap of vaccinia {row['vaccinia_gene']}; "
                f"a stronger molluscum homolog scores {best_score[row['vaccinia_gene']]}"
            )
            drop_ids.add(row["mcv_protein_id"])
    entries = [entry for entry in entries if entry["protein_id"] not in drop_ids]
    if not entries:
        raise SystemExit("no molluscum ORF passed the vaccinia membrane-protein filter")
    return entries, audit


def assert_panel(entries: list[dict], viruses: list[dict]) -> None:
    panel_ids = {row["virus_id"] for row in viruses if row["role"] == "panel"}
    covered = {entry["virus_id"] for entry in entries if entry["role"] == "panel"}
    missing = sorted(panel_ids - covered)
    extra = sorted(covered - panel_ids)
    if missing or extra or len(panel_ids) != 50:
        raise SystemExit(f"surface catalog coverage error: missing {missing}, extra {extra}")
    if not any(entry["virus_id"] == "AAV2" and entry["role"] == "control" for entry in entries):
        raise SystemExit("AAV2 VP1 control is missing from the surface catalog")


def write_catalog(entries: list[dict]) -> None:
    catalog = DATA / "surface_catalog.tsv"
    fasta_path = DATA / "surface_proteins.faa"
    write_tsv(
        catalog,
        entries,
        [
            "virus_id",
            "role",
            "protein_id",
            "product",
            "compartment",
            "start",
            "end",
            "length",
            "selection",
        ],
    )
    lines = []
    for entry in entries:
        lines.append(f">{entry['virus_id']}|{entry['protein_id']}|{entry['start']}-{entry['end']} {entry['product']}")
        seq = entry["seq"]
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    fasta_path.write_text("\n".join(lines) + "\n")


def best_window(query: np.ndarray, target: np.ndarray, length: int) -> tuple[int, int, int] | None:
    if len(query) < length or len(target) < length:
        return None
    query_windows = np.lib.stride_tricks.sliding_window_view(query, length)
    target_windows = np.lib.stride_tricks.sliding_window_view(target, length)
    query_ok = np.all(query_windows >= 0, axis=1)
    target_ok = np.all(target_windows >= 0, axis=1)
    if not np.any(query_ok) or not np.any(target_ok):
        return None
    matches = (query_windows[:, None, :] == target_windows[None, :, :]).sum(axis=2)
    if not query_ok.all():
        matches[~query_ok, :] = -1
    if not target_ok.all():
        matches[:, ~target_ok] = -1
    flat = int(matches.argmax())
    query_index, target_index = divmod(flat, matches.shape[1])
    identity = int(matches[query_index, target_index])
    if identity < 0:
        return None
    return query_index, target_index, identity


def epitope_overlap(start: int, end: int, epitopes: list[dict]) -> str:
    names = [
        row["epitope_id"]
        for row in epitopes
        if start <= int(row["end"]) and end >= int(row["start"])
    ]
    return ",".join(names)


def surface_kmer_sets(entries: list[dict], k: int) -> tuple[set[str], dict[str, list[tuple[str, str, int]]]]:
    found: set[str] = set()
    locations: dict[str, list[tuple[str, str, int]]] = {}
    for entry in entries:
        if entry["role"] != "panel":
            continue
        seq = entry["seq"]
        for pos in range(0, len(seq) - k + 1):
            peptide = seq[pos : pos + k]
            if any(char not in AA_INDEX for char in peptide):
                continue
            found.add(peptide)
            bucket = locations.get(peptide)
            item = (entry["virus_id"], entry["label"], pos + 1)
            if bucket is None:
                locations[peptide] = [item]
            elif len(bucket) < 8 and item not in bucket:
                bucket.append(item)
    return found, locations


def surface_nine_codes(entries: list[dict]) -> tuple[set[int], dict[int, list[tuple[str, str, int]]]]:
    found: set[int] = set()
    locations: dict[int, list[tuple[str, str, int]]] = {}
    for entry in entries:
        if entry["role"] != "panel":
            continue
        seq = entry["seq"]
        for pos in range(0, len(seq) - NINE + 1):
            code = encode_kmer(seq[pos : pos + NINE])
            if code is None:
                continue
            found.add(code)
            bucket = locations.get(code)
            item = (entry["virus_id"], entry["label"], pos + 1)
            if bucket is None:
                locations[code] = [item]
            elif len(bucket) < 8 and item not in bucket:
                bucket.append(item)
    return found, locations


def best_epitope_in_alignment(alignment, epitopes: list[dict]) -> dict | None:
    best = None
    for epitope in epitopes:
        stats = epitope_in_alignment(alignment, int(epitope["start"]), int(epitope["end"]))
        if stats is None or stats["aligned_residues"] < min(6, len(epitope["sequence"])):
            continue
        rank = (stats["matches"], stats["percent_identity"], stats["aligned_residues"])
        if best is None or rank > best[0]:
            best = (rank, epitope, stats)
    if best is None:
        return None
    _rank, epitope, stats = best
    return {
        "epitope_id": epitope["epitope_id"],
        "matches": stats["matches"],
        "aligned_residues": stats["aligned_residues"],
        "percent_identity": stats["percent_identity"],
    }


def exact_window_count(query: str, k: int, target: set[str]) -> int:
    count = 0
    for pos in range(0, len(query) - k + 1):
        if query[pos : pos + k] in target:
            count += 1
    return count


def global_best_window(query: str, arrays: list[np.ndarray], length: int) -> int:
    query_array = to_array(query)
    best = -1
    for target in arrays:
        hit = best_window(query_array, target, length)
        if hit is not None and hit[2] > best:
            best = hit[2]
    return best


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    viruses = load_viruses()
    virus_by_id = {row["virus_id"]: row for row in viruses}
    aav9 = load_aav9()
    epitopes = load_epitopes(aav9)
    index = load_protein_index()
    if not FEATURE_TABLE.exists():
        raise SystemExit(f"missing cached feature table {FEATURE_TABLE}")
    features = parse_feature_table(FEATURE_TABLE)

    entries = []
    entries.extend(explicit_entries(index))
    entries.extend(named_envelope_entries(index))
    entries.extend(mature_chain_entries(index, features))
    entries.extend(hav_entries(index))
    entries.append(sapovirus_entry(index, features))

    aligner = PairwiseAligner(
        mode="local",
        substitution_matrix=penalized_blosum62(),
        open_gap_score=-11,
        extend_gap_score=-1,
    )
    mcv_entries, mcv_audit = molluscum_entries(index, aligner)
    entries.extend(mcv_entries)
    write_tsv(
        RESULTS / "mcv_surface_homologs.tsv",
        mcv_audit,
        [
            "mcv_protein_id",
            "mcv_product",
            "mcv_length",
            "vaccinia_gene",
            "vaccinia_protein_id",
            "vaccinia_length",
            "score",
            "aligned_columns",
            "percent_identity",
            "coverage_of_vaccinia",
            "accepted",
            "note",
        ],
    )
    accepted = [row for row in mcv_audit if row["accepted"] == "yes"]
    print("molluscum surface homologs:", flush=True)
    for row in accepted:
        print(
            f"  {row['mcv_product']} {row['mcv_protein_id']} -> {row['vaccinia_gene']} "
            f"{row['percent_identity']}% over {row['aligned_columns']} columns",
            flush=True,
        )

    assert_panel(entries, viruses)
    write_catalog(entries)
    panel_entries = [entry for entry in entries if entry["role"] == "panel"]
    panel_residues = sum(entry["length"] for entry in panel_entries)
    print(
        f"surface set: {len(panel_entries)} proteins, {panel_residues} residues, "
        f"{len({entry['virus_id'] for entry in panel_entries})} viruses",
        flush=True,
    )

    print("scoring ungapped windows", flush=True)
    window_rows = []
    for virus_id in [row["virus_id"] for row in viruses]:
        group = [entry for entry in entries if entry["virus_id"] == virus_id]
        if not group:
            continue
        group = sorted(group, key=lambda entry: entry["length"], reverse=True)
        meta = virus_by_id[virus_id]
        for length in WINDOW_LENGTHS:
            best = None
            for entry in group:
                hit = best_window(to_array(aav9), to_array(entry["seq"]), length)
                if hit is None:
                    continue
                if best is None or hit[2] > best[0]:
                    query_index, target_index, identity = hit
                    best = (identity, query_index, target_index, entry)
            if best is None:
                continue
            identity, query_index, target_index, entry = best
            peptide = aav9[query_index : query_index + length]
            viral = entry["seq"][target_index : target_index + length]
            start = query_index + 1
            end = query_index + length
            window_rows.append(
                {
                    "virus_id": virus_id,
                    "virus_name": meta["name"],
                    "role": meta["role"],
                    "group": meta["group"],
                    "window": length,
                    "vp1_start": start,
                    "vp1_end": end,
                    "peptide": peptide,
                    "identity": identity,
                    "percent_identity": round(100 * identity / length, 1),
                    "mismatches": length - identity,
                    "label": classify(identity, length),
                    "viral_peptide": viral,
                    "viral_start": target_index + 1,
                    "viral_protein": entry["label"],
                    "epitope_overlap": epitope_overlap(start, end, epitopes),
                    "surface_proteins": len(group),
                    "surface_residues": sum(item["length"] for item in group),
                }
            )
        print(f"  {virus_id}", flush=True)

    print("indexing surface k-mers", flush=True)
    peptide_rows = []
    seven_set, seven_locations = surface_kmer_sets(panel_entries, 7)
    for pos in range(0, len(aav9) - 6):
        peptide = aav9[pos : pos + 7]
        locations = seven_locations.get(peptide)
        if not locations:
            continue
        for virus_id, label, viral_start in locations:
            peptide_rows.append(
                {
                    "kind": "exact7",
                    "vp1_start": pos + 1,
                    "vp1_end": pos + 7,
                    "peptide": peptide,
                    "identity": 7,
                    "virus_id": virus_id,
                    "role": "panel",
                    "viral_peptide": peptide,
                    "viral_start": viral_start,
                    "viral_protein": label,
                    "epitope_overlap": epitope_overlap(pos + 1, pos + 7, epitopes),
                }
            )

    nine_codes, nine_locations = surface_nine_codes(panel_entries)
    nine_pairs = 0
    nine_windows = 0
    for pos in range(0, len(aav9) - NINE + 1):
        peptide = aav9[pos : pos + NINE]
        variants = nine_mer_variants(peptide)
        if not variants:
            continue
        hits = []
        for code, identity in variants:
            locations = nine_locations.get(code)
            if not locations:
                continue
            for virus_id, label, viral_start in locations:
                hits.append((identity, code, virus_id, label, viral_start))
        if not hits:
            continue
        nine_windows += 1
        nine_pairs += len({(item[1], item[2], item[3]) for item in hits})
        for identity, _code, virus_id, label, viral_start in hits:
            entry = next(item for item in panel_entries if item["virus_id"] == virus_id and item["label"] == label)
            viral_peptide = entry["seq"][viral_start - 1 : viral_start - 1 + NINE]
            peptide_rows.append(
                {
                    "kind": "nine",
                    "vp1_start": pos + 1,
                    "vp1_end": pos + NINE,
                    "peptide": peptide,
                    "identity": identity,
                    "virus_id": virus_id,
                    "role": "panel",
                    "viral_peptide": viral_peptide,
                    "viral_start": viral_start,
                    "viral_protein": label,
                    "epitope_overlap": epitope_overlap(pos + 1, pos + NINE, epitopes),
                }
            )

    print("aligning VP1 to each surface protein", flush=True)
    alignment_rows = []
    for entry in entries:
        alignment = next(aligner.align(for_alignment(aav9), for_alignment(entry["seq"])))
        stats = alignment_stats(alignment)
        epitope_stats = best_epitope_in_alignment(alignment, epitopes)
        meta = virus_by_id[entry["virus_id"]]
        alignment_rows.append(
            {
                "virus_id": entry["virus_id"],
                "virus_name": meta["name"],
                "role": entry["role"],
                "group": meta["group"],
                "protein": entry["label"],
                "protein_length": entry["length"],
                "compartment": entry["compartment"],
                "score": round(stats["score"], 1),
                "aligned_columns": stats["aligned_columns"],
                "matches": stats["matches"],
                "percent_identity": round(100 * stats["identity"], 1),
                "query_start": stats["query_start"],
                "query_end": stats["query_end"],
                "target_start": stats["target_start"],
                "target_end": stats["target_end"],
                "best_epitope": epitope_stats["epitope_id"] if epitope_stats else "",
                "best_epitope_matches": epitope_stats["matches"] if epitope_stats else "",
                "best_epitope_aligned": epitope_stats["aligned_residues"] if epitope_stats else "",
                "best_epitope_percent": epitope_stats["percent_identity"] if epitope_stats else "",
            }
        )

    print("running composition-matched null", flush=True)
    panel_arrays = [to_array(entry["seq"]) for entry in panel_entries if len(entry["seq"]) >= 15]
    rng = random.Random(NULL_SEED)
    null_exact = []
    null_nines = []
    null_best15 = []
    for replicate in range(NULL_REPLICATES):
        shuffled = composition_shuffle(aav9, rng)
        null_exact.append(exact_window_count(shuffled, 7, seven_set))
        windows, _pairs = nine_mer_hit_count(shuffled, nine_codes)
        null_nines.append(windows)
        null_best15.append(global_best_window(shuffled, panel_arrays, 15))
        print(f"  null {replicate + 1}/{NULL_REPLICATES}", flush=True)

    observed_exact = exact_window_count(aav9, 7, seven_set)
    observed_nine_windows, observed_nine_pairs = nine_mer_hit_count(aav9, nine_codes)
    if observed_nine_windows != nine_windows:
        raise SystemExit(
            f"9-mer window count mismatch: scan {nine_windows}, counter {observed_nine_windows}"
        )

    def mean(values: list[int]) -> float:
        return round(sum(values) / len(values), 2) if values else 0.0

    best15_panel = [
        row for row in window_rows if row["window"] == 15 and row["role"] == "panel"
    ]
    best15_panel.sort(key=lambda row: (-row["identity"], row["virus_id"]))
    best_alignment = [
        row for row in alignment_rows if row["role"] == "panel"
    ]
    best_alignment.sort(key=lambda row: (-row["percent_identity"], -row["aligned_columns"]))
    aav2_alignment = next(row for row in alignment_rows if row["virus_id"] == "AAV2")
    if aav2_alignment["percent_identity"] < 75 or aav2_alignment["aligned_columns"] < 700:
        raise SystemExit(f"AAV2 control alignment looks wrong: {aav2_alignment}")

    summary = {
        "query": "AAV9 VP1",
        "query_length": len(aav9),
        "surface_proteins": len(panel_entries),
        "surface_residues": panel_residues,
        "panel_viruses": 50,
        "exact_7mer_windows": observed_exact,
        "nine_mer_windows_at_least_7": observed_nine_windows,
        "nine_mer_pairs_at_least_7": observed_nine_pairs,
        "null_replicates": NULL_REPLICATES,
        "null_seed": NULL_SEED,
        "null_exact_7mer_mean": mean(null_exact),
        "null_exact_7mer_max": max(null_exact),
        "null_exact_7mer_values": null_exact,
        "null_nine_windows_mean": mean(null_nines),
        "null_nine_windows_max": max(null_nines),
        "null_nine_windows_values": null_nines,
        "null_best_15mer_identity_mean": mean(null_best15),
        "null_best_15mer_identity_max": max(null_best15),
        "null_best_15mer_identity_values": null_best15,
        "best_15mer_panel": best15_panel[0],
        "aav2_alignment_percent_identity": aav2_alignment["percent_identity"],
        "aav2_alignment_columns": aav2_alignment["aligned_columns"],
        "viruses": [],
    }
    for virus in viruses:
        virus_windows = [row for row in window_rows if row["virus_id"] == virus["virus_id"] and row["window"] == 15]
        virus_alignments = [row for row in alignment_rows if row["virus_id"] == virus["virus_id"]]
        if not virus_windows or not virus_alignments:
            continue
        top_alignment = max(
            virus_alignments,
            key=lambda row: (float(row["score"]), int(row["aligned_columns"])),
        )
        top_window = virus_windows[0]
        summary["viruses"].append(
            {
                "virus_id": virus["virus_id"],
                "role": virus["role"],
                "surface_proteins": top_window["surface_proteins"],
                "surface_residues": top_window["surface_residues"],
                "best_15_identity": top_window["identity"],
                "best_15_percent": top_window["percent_identity"],
                "best_15_peptide": top_window["peptide"],
                "best_15_viral_peptide": top_window["viral_peptide"],
                "best_15_protein": top_window["viral_protein"],
                "best_15_vp1_start": top_window["vp1_start"],
                "best_15_label": top_window["label"],
                "best_15_epitope_overlap": top_window["epitope_overlap"],
                "best_alignment_score": top_alignment["score"],
                "best_alignment_percent": top_alignment["percent_identity"],
                "best_alignment_columns": top_alignment["aligned_columns"],
                "best_alignment_protein": top_alignment["protein"],
                "best_alignment_query": f"{top_alignment['query_start']}-{top_alignment['query_end']}",
            }
        )

    outside = [
        row for row in alignment_rows if row["virus_id"] not in {"B19", "HBOV1", "AAV2"}
    ]
    summary["max_alignment_score_outside_parvoviruses"] = max(float(row["score"]) for row in outside)
    summary["max_15mer_identity_outside_parvoviruses"] = max(
        int(row["identity"]) for row in best15_panel if row["virus_id"] not in {"B19", "HBOV1"}
    )

    write_tsv(
        RESULTS / "surface_best_windows.tsv",
        window_rows,
        [
            "virus_id",
            "virus_name",
            "role",
            "group",
            "window",
            "vp1_start",
            "vp1_end",
            "peptide",
            "identity",
            "percent_identity",
            "mismatches",
            "label",
            "viral_peptide",
            "viral_start",
            "viral_protein",
            "epitope_overlap",
            "surface_proteins",
            "surface_residues",
        ],
    )
    write_tsv(
        RESULTS / "surface_peptide_hits.tsv",
        peptide_rows,
        [
            "kind",
            "vp1_start",
            "vp1_end",
            "peptide",
            "identity",
            "virus_id",
            "role",
            "viral_peptide",
            "viral_start",
            "viral_protein",
            "epitope_overlap",
        ],
    )
    write_tsv(
        RESULTS / "surface_alignments.tsv",
        alignment_rows,
        [
            "virus_id",
            "virus_name",
            "role",
            "group",
            "protein",
            "protein_length",
            "compartment",
            "score",
            "aligned_columns",
            "matches",
            "percent_identity",
            "query_start",
            "query_end",
            "target_start",
            "target_end",
            "best_epitope",
            "best_epitope_matches",
            "best_epitope_aligned",
            "best_epitope_percent",
        ],
    )
    (RESULTS / "surface_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in summary if key != "viruses"}, indent=2))


if __name__ == "__main__":
    main()
