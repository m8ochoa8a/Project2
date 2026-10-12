#!/usr/bin/env python3
"""Compare AAV9 capsid epitopes with proteomes of common human viruses.

The screen is ungapped peptide identity. It asks whether a published AAV9
epitope, or any capsid 6-/7-/9-mer, has a close counterpart in each viral
proteome. Sequence identity is not evidence that T or B cells cross-react.
"""

from __future__ import annotations

import csv
import json
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PROT_DIR = DATA / "proteomes"
RESULTS = ROOT / "results"
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA)}
KMER_LENGTHS = (6, 7)
NINE = 9
NINE_MIN_IDENTITY = 7
NULL_REPLICATES = 25
NULL_SEED = 20261012


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip().upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def clean_protein(seq: str) -> str:
    return seq.replace("*", "").replace(" ", "")


def load_aav9() -> str:
    records = read_fasta(DATA / "aav9_vp1.fasta")
    if len(records) != 1:
        raise SystemExit("expected one AAV9 VP1 record")
    seq = clean_protein(records[0][1])
    if len(seq) != 736:
        raise SystemExit(f"AAV9 VP1 length is {len(seq)}, expected 736")
    return seq


def load_epitopes(aav9: str) -> list[dict]:
    epitopes = []
    with (DATA / "epitopes.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            start = int(row["start"])
            end = int(row["end"])
            sliced = aav9[start - 1 : end]
            if sliced != row["sequence"]:
                raise SystemExit(
                    f"{row['epitope_id']} sequence does not match VP1 {start}-{end}: {sliced}"
                )
            row["start"] = start
            row["end"] = end
            epitopes.append(row)
    return epitopes


def load_viruses() -> list[dict]:
    with (DATA / "viruses.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    panel = [row for row in rows if row["role"] == "panel"]
    if len(panel) != 50:
        raise SystemExit(f"expected 50 panel viruses, found {len(panel)}")
    return rows


def protein_label(header: str) -> str:
    protein = ""
    protein_id = ""
    if "[protein=" in header:
        protein = header.split("[protein=", 1)[1].split("]", 1)[0]
    if "[protein_id=" in header:
        protein_id = header.split("[protein_id=", 1)[1].split("]", 1)[0]
    if protein and protein_id:
        return f"{protein} ({protein_id})"
    if protein:
        return protein
    return header.split(" [", 1)[0][:140]


def iter_proteins(virus_id: str):
    for header, seq in read_fasta(PROT_DIR / f"{virus_id}.faa"):
        seq = clean_protein(seq)
        if len(seq) < 6:
            continue
        yield protein_label(header), seq


def to_array(seq: str) -> np.ndarray:
    arr = np.empty(len(seq), dtype=np.int16)
    for i, char in enumerate(seq):
        arr[i] = AA_INDEX.get(char, -1)
    return arr


def best_ungapped(epi: np.ndarray, prot: np.ndarray) -> tuple[int, int] | None:
    length = len(epi)
    if len(prot) < length or np.any(epi < 0):
        return None
    windows = np.lib.stride_tricks.sliding_window_view(prot, length)
    valid = np.all(windows >= 0, axis=1)
    if not np.any(valid):
        return None
    matches = (windows == epi).sum(axis=1)
    matches = np.where(valid, matches, -1)
    index = int(matches.argmax())
    identity = int(matches[index])
    if identity < 0:
        return None
    return index, identity


def classify(identity: int, length: int) -> str:
    if identity == length and length >= 7:
        return "exact"
    fraction = identity / length
    mismatches = length - identity
    if length >= 8 and fraction >= 0.8 and mismatches <= 2:
        return "high"
    if length >= 8 and fraction >= 0.7 and mismatches <= 3:
        return "moderate"
    return ""


def index_kmers(virus_rows: list[dict], k: int) -> dict[str, list[tuple[str, str, int]]]:
    index: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    for row in virus_rows:
        for name, seq in iter_proteins(row["virus_id"]):
            for pos in range(0, len(seq) - k + 1):
                mer = seq[pos : pos + k]
                if any(char not in AA_INDEX for char in mer):
                    continue
                hits = index[mer]
                if len(hits) < 8:
                    hits.append((row["virus_id"], name, pos + 1))
    return index


def viral_kmer_set(virus_rows: list[dict], k: int) -> set[str]:
    found = set()
    for row in virus_rows:
        for _name, seq in iter_proteins(row["virus_id"]):
            for pos in range(0, len(seq) - k + 1):
                mer = seq[pos : pos + k]
                if all(char in AA_INDEX for char in mer):
                    found.add(mer)
    return found


def viral_encoded_nines(virus_rows: list[dict]) -> tuple[set[int], dict[int, list[tuple[str, str, int]]]]:
    found: set[int] = set()
    locations: dict[int, list[tuple[str, str, int]]] = {}
    for row in virus_rows:
        for name, seq in iter_proteins(row["virus_id"]):
            for pos in range(0, len(seq) - NINE + 1):
                code = encode_kmer(seq[pos : pos + NINE])
                if code is None:
                    continue
                found.add(code)
                bucket = locations.get(code)
                if bucket is None:
                    locations[code] = [(row["virus_id"], name, pos + 1)]
                elif len(bucket) < 8:
                    bucket.append((row["virus_id"], name, pos + 1))
    return found, locations


def aav_kmers(seq: str, k: int) -> list[tuple[int, str]]:
    return [
        (i + 1, seq[i : i + k])
        for i in range(0, len(seq) - k + 1)
        if all(char in AA_INDEX for char in seq[i : i + k])
    ]


POW9 = tuple(20 ** i for i in range(NINE - 1, -1, -1))


def encode_kmer(peptide: str) -> int | None:
    value = 0
    for char in peptide:
        index = AA_INDEX.get(char)
        if index is None:
            return None
        value = value * 20 + index
    return value


def decode_kmer(value: int, length: int = NINE) -> str:
    chars = []
    for _ in range(length):
        value, index = divmod(value, 20)
        chars.append(AA[index])
    return "".join(reversed(chars))


def nine_mer_variants(peptide: str) -> list[tuple[int, int]]:
    """Encoded 9-mers within 2 substitutions, as (code, identity)."""
    code = encode_kmer(peptide)
    if code is None:
        return []
    digits = [AA_INDEX[char] for char in peptide]
    found = [(code, NINE)]
    for i, digit in enumerate(digits):
        weight = POW9[i]
        for new_digit in range(20):
            if new_digit == digit:
                continue
            one = code + (new_digit - digit) * weight
            found.append((one, NINE - 1))
            for j in range(i + 1, NINE):
                other = digits[j]
                other_weight = POW9[j]
                for second in range(20):
                    if second == other:
                        continue
                    found.append((one + (second - other) * other_weight, NINE - 2))
    return found


def nine_mer_hit_count(query: str, viral_nines: set[int]) -> tuple[int, int]:
    """Return (query 9-mers with a hit, total hit pairs) at >=7/9 identity."""
    peptides_with_hit = 0
    pairs = 0
    for i in range(0, len(query) - NINE + 1):
        variants = nine_mer_variants(query[i : i + NINE])
        if not variants:
            continue
        hit_pairs = 0
        for code, _identity in variants:
            if code in viral_nines:
                hit_pairs += 1
        if hit_pairs:
            peptides_with_hit += 1
            pairs += hit_pairs
    return peptides_with_hit, pairs


def composition_shuffle(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def alignment_stats(alignment) -> dict:
    columns = 0
    matches = 0
    for q_char, t_char in zip(alignment[0], alignment[1]):
        if q_char == "-" or t_char == "-":
            continue
        columns += 1
        if q_char == t_char:
            matches += 1
    return {
        "score": float(alignment.score),
        "aligned_columns": columns,
        "matches": matches,
        "identity": matches / columns if columns else 0.0,
        "query_start": int(alignment.coordinates[0][0]) + 1,
        "query_end": int(alignment.coordinates[0][-1]),
        "target_start": int(alignment.coordinates[1][0]) + 1,
        "target_end": int(alignment.coordinates[1][-1]),
    }


def epitope_in_alignment(alignment, start: int, end: int) -> dict | None:
    """Identity of one AAV9 interval inside a local alignment. Coordinates are 1-based inclusive."""
    qpos = int(alignment.coordinates[0][0])
    tpos = int(alignment.coordinates[1][0])
    matched = 0
    compared = 0
    target_chars: list[str] = []
    query_chars: list[str] = []
    target_start = None
    for q_char, t_char in zip(alignment[0], alignment[1]):
        this_q = None
        if q_char != "-":
            qpos += 1
            this_q = qpos
        if t_char != "-":
            tpos += 1
        if this_q is None or not (start <= this_q <= end):
            continue
        query_chars.append(q_char)
        target_chars.append(t_char)
        if t_char == "-":
            continue
        compared += 1
        if target_start is None:
            target_start = tpos
        if q_char == t_char:
            matched += 1
    if compared == 0:
        return None
    return {
        "aligned_residues": compared,
        "matches": matched,
        "percent_identity": round(100 * matched / compared, 1),
        "target_start": target_start or "",
        "query_aligned": "".join(query_chars),
        "target_aligned": "".join(target_chars),
    }


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    aav9 = load_aav9()
    epitopes = load_epitopes(aav9)
    viruses = load_viruses()
    panel = [row for row in viruses if row["role"] == "panel"]
    controls = [row for row in viruses if row["role"] == "control"]

    epi_arrays = {row["epitope_id"]: to_array(row["sequence"]) for row in epitopes}
    best_rows = []
    notable_rows = []

    print("scoring published epitopes", flush=True)
    for virus in panel + controls:
        protein_count = 0
        residue_count = 0
        best_for_virus = {
            row["epitope_id"]: {"identity": -1, "start": None, "peptide": "", "protein": ""}
            for row in epitopes
        }
        for protein_name, seq in iter_proteins(virus["virus_id"]):
            protein_count += 1
            residue_count += len(seq)
            prot = to_array(seq)
            for epi in epitopes:
                hit = best_ungapped(epi_arrays[epi["epitope_id"]], prot)
                if hit is None:
                    continue
                start0, identity = hit
                current = best_for_virus[epi["epitope_id"]]
                if identity > current["identity"]:
                    length = len(epi["sequence"])
                    peptide = seq[start0 : start0 + length]
                    current.update(
                        {
                            "identity": identity,
                            "start": start0 + 1,
                            "peptide": peptide,
                            "protein": protein_name,
                        }
                    )
        for epi in epitopes:
            current = best_for_virus[epi["epitope_id"]]
            length = len(epi["sequence"])
            identity = current["identity"]
            best_rows.append(
                {
                    "virus_id": virus["virus_id"],
                    "virus_name": virus["name"],
                    "role": virus["role"],
                    "group": virus["group"],
                    "epitope_id": epi["epitope_id"],
                    "epitope_class": epi["epitope_class"],
                    "epitope": epi["sequence"],
                    "length": length,
                    "identity": identity,
                    "percent_identity": round(100 * identity / length, 1) if identity >= 0 else "",
                    "mismatches": length - identity if identity >= 0 else "",
                    "label": classify(identity, length) if identity >= 0 else "",
                    "viral_peptide": current["peptide"],
                    "viral_start": current["start"] or "",
                    "viral_protein": current["protein"],
                    "proteins": protein_count,
                    "residues": residue_count,
                }
            )
            if best_rows[-1]["label"]:
                notable_rows.append(
                    {
                        "virus_id": virus["virus_id"],
                        "virus_name": virus["name"],
                        "role": virus["role"],
                        "epitope_id": epi["epitope_id"],
                        "epitope": epi["sequence"],
                        "epitope_class": epi["epitope_class"],
                        "vp1_start": epi["start"],
                        "vp1_end": epi["end"],
                        "identity": identity,
                        "length": length,
                        "percent_identity": best_rows[-1]["percent_identity"],
                        "mismatches": best_rows[-1]["mismatches"],
                        "label": best_rows[-1]["label"],
                        "viral_peptide": current["peptide"],
                        "viral_start": current["start"] or "",
                        "viral_protein": current["protein"],
                    }
                )
        print(f"  {virus['virus_id']}: {protein_count} proteins, {residue_count} aa", flush=True)

    print("indexing exact k-mers", flush=True)
    kmer_rows = []
    kmer_summary = {}
    for k in KMER_LENGTHS:
        index = index_kmers(panel, k)
        hits = []
        for start, mer in aav_kmers(aav9, k):
            locations = index.get(mer)
            if not locations:
                continue
            viruses_hit = sorted({item[0] for item in locations})
            hits.append(
                {
                    "k": k,
                    "vp1_start": start,
                    "vp1_end": start + k - 1,
                    "peptide": mer,
                    "n_viruses": len(viruses_hit),
                    "viruses": ",".join(viruses_hit),
                    "example_virus": locations[0][0],
                    "example_protein": locations[0][1],
                    "example_start": locations[0][2],
                }
            )
        kmer_rows.extend(hits)
        kmer_summary[str(k)] = {
            "shared_windows": len(hits),
            "unique_peptides": len({row["peptide"] for row in hits}),
        }
        print(f"  exact {k}-mers shared with the panel: {len(hits)}", flush=True)

    print("null model", flush=True)
    viral_nines, nine_locations = viral_encoded_nines(panel)
    viral_six = viral_kmer_set(panel, 6)
    viral_seven = viral_kmer_set(panel, 7)
    nine_hit_rows = []
    for start in range(0, len(aav9) - NINE + 1):
        peptide = aav9[start : start + NINE]
        for code, identity in nine_mer_variants(peptide):
            locations = nine_locations.get(code)
            if not locations:
                continue
            viral_peptide = decode_kmer(code)
            viruses_hit = sorted({item[0] for item in locations})
            nine_hit_rows.append(
                {
                    "vp1_start": start + 1,
                    "vp1_end": start + NINE,
                    "aav9_peptide": peptide,
                    "identity": identity,
                    "viral_peptide": viral_peptide,
                    "n_viruses": len(viruses_hit),
                    "viruses": ",".join(viruses_hit),
                    "example_virus": locations[0][0],
                    "example_protein": locations[0][1],
                    "example_start": locations[0][2],
                }
            )
    observed_9_peptides = len({row["vp1_start"] for row in nine_hit_rows})
    observed_9_pairs = len(nine_hit_rows)
    observed_6 = sum(1 for _start, mer in aav_kmers(aav9, 6) if mer in viral_six)
    observed_7 = sum(1 for _start, mer in aav_kmers(aav9, 7) if mer in viral_seven)
    rng = random.Random(NULL_SEED)
    null_9 = []
    null_6 = []
    null_7 = []
    for replicate in range(NULL_REPLICATES):
        shuffled = composition_shuffle(aav9, rng)
        peptides, pairs = nine_mer_hit_count(shuffled, viral_nines)
        null_9.append({"peptides": peptides, "pairs": pairs})
        null_6.append(sum(1 for _start, mer in aav_kmers(shuffled, 6) if mer in viral_six))
        null_7.append(sum(1 for _start, mer in aav_kmers(shuffled, 7) if mer in viral_seven))
        print(
            f"  shuffle {replicate + 1}/{NULL_REPLICATES}: "
            f"9-mer peptides {peptides}, 7-mers {null_7[-1]}, 6-mers {null_6[-1]}",
            flush=True,
        )

    def mean(values: list[float]) -> float:
        return float(np.mean(values)) if values else 0.0

    null_summary = {
        "replicates": NULL_REPLICATES,
        "seed": NULL_SEED,
        "observed_9mer_peptides_with_ge_7_of_9": observed_9_peptides,
        "observed_9mer_pairs_ge_7_of_9": observed_9_pairs,
        "null_9mer_peptides_mean": mean([row["peptides"] for row in null_9]),
        "null_9mer_peptides_max": max(row["peptides"] for row in null_9),
        "null_9mer_pairs_mean": mean([row["pairs"] for row in null_9]),
        "observed_exact_7mers": observed_7,
        "null_exact_7mers_mean": mean(null_7),
        "null_exact_7mers_max": max(null_7),
        "observed_exact_6mers": observed_6,
        "null_exact_6mers_mean": mean(null_6),
        "null_exact_6mers_max": max(null_6),
        "viral_unique_9mers": len(viral_nines),
        "viral_unique_7mers": len(viral_seven),
        "viral_unique_6mers": len(viral_six),
    }

    print("parvovirus local alignments", flush=True)
    aligner = PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )
    alignment_rows = []
    conserved_rows = []
    for virus in panel + controls:
        if virus["virus_id"] not in {"B19", "HBOV1", "AAV2"}:
            continue
        for protein_name, seq in iter_proteins(virus["virus_id"]):
            if len(seq) > 2500:
                continue
            alignment = aligner.align(aav9, seq)[0]
            result = alignment_stats(alignment)
            if result["aligned_columns"] < 80:
                continue
            result.update(
                {
                    "virus_id": virus["virus_id"],
                    "virus_name": virus["name"],
                    "role": virus["role"],
                    "viral_protein": protein_name,
                    "viral_length": len(seq),
                    "percent_identity": round(100 * result["identity"], 1),
                }
            )
            result.pop("alignment", None)
            alignment_rows.append(result)
            print(
                f"  {virus['virus_id']} {protein_name[:70]} "
                f"{result['percent_identity']}% over {result['aligned_columns']} aa",
                flush=True,
            )
            if "VP1" in protein_name or "capsid protein 1" in protein_name or "major coat protein VP1" in protein_name:
                for epi in epitopes:
                    mapped = epitope_in_alignment(alignment, epi["start"], epi["end"])
                    if mapped is None:
                        continue
                    conserved_rows.append(
                        {
                            "virus_id": virus["virus_id"],
                            "viral_protein": protein_name,
                            "epitope_id": epi["epitope_id"],
                            "epitope_class": epi["epitope_class"],
                            "vp1_start": epi["start"],
                            "vp1_end": epi["end"],
                            "epitope": epi["sequence"],
                            **mapped,
                        }
                    )

    # Best panel hit per epitope, for the narrative table.
    overall = []
    for epi in epitopes:
        candidates = [
            row
            for row in best_rows
            if row["epitope_id"] == epi["epitope_id"] and row["role"] == "panel" and row["identity"] >= 0
        ]
        candidates.sort(key=lambda row: (row["identity"], row["percent_identity"]), reverse=True)
        best = candidates[0]
        control_hits = [
            row
            for row in best_rows
            if row["epitope_id"] == epi["epitope_id"] and row["role"] == "control"
        ]
        control = control_hits[0] if control_hits else None
        overall.append(
            {
                **best,
                "aav2_identity": control["identity"] if control else "",
                "aav2_percent_identity": control["percent_identity"] if control else "",
                "aav2_peptide": control["viral_peptide"] if control else "",
            }
        )

    write_tsv(
        RESULTS / "nine_mer_hits.tsv",
        nine_hit_rows,
        [
            "vp1_start",
            "vp1_end",
            "aav9_peptide",
            "identity",
            "viral_peptide",
            "n_viruses",
            "viruses",
            "example_virus",
            "example_protein",
            "example_start",
        ],
    )
    write_tsv(
        RESULTS / "epitope_by_virus.tsv",
        best_rows,
        [
            "virus_id",
            "virus_name",
            "role",
            "group",
            "epitope_id",
            "epitope_class",
            "epitope",
            "length",
            "identity",
            "percent_identity",
            "mismatches",
            "label",
            "viral_peptide",
            "viral_start",
            "viral_protein",
            "proteins",
            "residues",
        ],
    )
    write_tsv(
        RESULTS / "epitope_best_hits.tsv",
        overall,
        [
            "epitope_id",
            "epitope_class",
            "epitope",
            "length",
            "virus_id",
            "virus_name",
            "identity",
            "percent_identity",
            "mismatches",
            "label",
            "viral_peptide",
            "viral_start",
            "viral_protein",
            "aav2_identity",
            "aav2_percent_identity",
            "aav2_peptide",
        ],
    )
    write_tsv(
        RESULTS / "notable_epitope_hits.tsv",
        notable_rows,
        [
            "virus_id",
            "virus_name",
            "role",
            "epitope_id",
            "epitope",
            "epitope_class",
            "vp1_start",
            "vp1_end",
            "identity",
            "length",
            "percent_identity",
            "mismatches",
            "label",
            "viral_peptide",
            "viral_start",
            "viral_protein",
        ],
    )
    write_tsv(
        RESULTS / "exact_kmers.tsv",
        kmer_rows,
        [
            "k",
            "vp1_start",
            "vp1_end",
            "peptide",
            "n_viruses",
            "viruses",
            "example_virus",
            "example_protein",
            "example_start",
        ],
    )
    write_tsv(
        RESULTS / "domain_alignments.tsv",
        alignment_rows,
        [
            "virus_id",
            "virus_name",
            "role",
            "viral_protein",
            "viral_length",
            "score",
            "aligned_columns",
            "matches",
            "percent_identity",
            "query_start",
            "query_end",
            "target_start",
            "target_end",
        ],
    )
    write_tsv(
        RESULTS / "parvovirus_epitope_identity.tsv",
        conserved_rows,
        [
            "virus_id",
            "viral_protein",
            "epitope_id",
            "epitope_class",
            "vp1_start",
            "vp1_end",
            "epitope",
            "aligned_residues",
            "matches",
            "percent_identity",
            "target_start",
            "query_aligned",
            "target_aligned",
        ],
    )

    panel_notable = [row for row in notable_rows if row["role"] == "panel"]
    control_notable = [row for row in notable_rows if row["role"] == "control"]
    summary = {
        "aav9_length": len(aav9),
        "n_epitopes": len(epitopes),
        "n_panel_viruses": len(panel),
        "panel_virus_ids": [row["virus_id"] for row in panel],
        "notable_definition": {
            "exact": "identical peptide, length >= 7",
            "high": "length >= 8, identity >= 80%, at most 2 mismatches",
            "moderate": "length >= 8, identity >= 70%, at most 3 mismatches",
        },
        "panel_notable_hits": len(panel_notable),
        "control_notable_hits": len(control_notable),
        "exact_kmer_summary": kmer_summary,
        "nine_mer_hits_ge_7_of_9": len(nine_hit_rows),
        "null": null_summary,
        "best_panel_percent_by_epitope": [
            {
                "epitope_id": row["epitope_id"],
                "percent_identity": row["percent_identity"],
                "virus_id": row["virus_id"],
            }
            for row in overall
        ],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"panel_notable_hits": len(panel_notable), "null": null_summary}, indent=2))


if __name__ == "__main__":
    main()
