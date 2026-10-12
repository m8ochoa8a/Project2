#!/usr/bin/env python3
"""Download reference viral CDS translations for the homology panel."""

from __future__ import annotations

import csv
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
VIRUSES = ROOT / "data" / "viruses.tsv"
OUTDIR = ROOT / "data" / "proteomes"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT = "https://rest.uniprot.org/uniprotkb/Q6JC40.fasta"
AA = set("ACDEFGHIKLMNPQRSTVWY")


def fetch_cds(accession: str) -> str:
    last_error = None
    for attempt in range(4):
        try:
            response = requests.get(
                EUTILS,
                params={
                    "db": "nuccore",
                    "id": accession,
                    "rettype": "fasta_cds_aa",
                    "retmode": "text",
                    "tool": "aav9-epitope-homology",
                    "email": "research@example.com",
                },
                timeout=180,
            )
            response.raise_for_status()
            text = response.text.strip()
            if not text.startswith(">"):
                raise RuntimeError(f"non-FASTA response for {accession}: {text[:200]}")
            return text + "\n"
        except Exception as exc:  # noqa: BLE001 - retry network and NCBI errors
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to fetch {accession}: {last_error}")


def protein_stats(fasta: str) -> tuple[int, int]:
    proteins = 0
    residues = 0
    seq_chars: list[str] = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if seq_chars:
                residues += sum(char in AA for char in seq_chars)
                seq_chars = []
            proteins += 1
        else:
            seq_chars.extend(line.strip().upper())
    if seq_chars:
        residues += sum(char in AA for char in seq_chars)
    if proteins == 0 or residues == 0:
        raise RuntimeError("FASTA contained no protein sequence")
    return proteins, residues


def confirm_aav9() -> None:
    response = requests.get(UNIPROT, timeout=60)
    response.raise_for_status()
    remote = "".join(
        line.strip() for line in response.text.splitlines() if line and not line.startswith(">")
    )
    local_path = ROOT / "data" / "aav9_vp1.fasta"
    local = "".join(
        line.strip() for line in local_path.read_text().splitlines() if line and not line.startswith(">")
    )
    if remote != local:
        raise SystemExit(
            f"Local AAV9 VP1 does not match UniProt Q6JC40 (local {len(local)} aa, UniProt {len(remote)} aa)"
        )
    print(f"AAV9 VP1 matches UniProt Q6JC40 ({len(local)} aa)")


def main() -> None:
    confirm_aav9()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTDIR / "manifest.tsv"
    rows = list(csv.DictReader(VIRUSES.open(), delimiter="\t"))
    manifest = []
    for row in rows:
        virus_id = row["virus_id"]
        accessions = [item for item in row["accessions"].split(";") if item]
        parts = []
        for accession in accessions:
            print(f"fetch {virus_id} {accession}", flush=True)
            parts.append(fetch_cds(accession))
            time.sleep(0.4)
        fasta = "".join(parts)
        proteins, residues = protein_stats(fasta)
        out = OUTDIR / f"{virus_id}.faa"
        out.write_text(fasta)
        print(f"  {proteins} proteins, {residues} aa", flush=True)
        manifest.append(
            {
                "virus_id": virus_id,
                "proteins": proteins,
                "residues": residues,
                "accessions": row["accessions"],
            }
        )
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["virus_id", "proteins", "residues", "accessions"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(manifest)
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
