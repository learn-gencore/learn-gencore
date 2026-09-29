#!/usr/bin/env python3
"""Download GEO supplementary data and assemble the selected D0 gene counts."""
from __future__ import annotations

import csv
import gzip
import io
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "metadata" / "GSE327304_D0_sample_manifest.csv"
RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "processed"
ARCHIVE = RAW_DIR / "GSE327304_RAW.tar"
COUNTS_OUT = OUT_DIR / "GSE327304_D0_raw_counts.tsv"
URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE327304&format=file"


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open(newline="") as handle:
        samples = list(csv.DictReader(handle))
    if not ARCHIVE.exists():
        print(f"Downloading {URL}")
        urllib.request.urlretrieve(URL, ARCHIVE)
    gene_ids: list[str] | None = None
    columns: list[list[str]] = []
    with tarfile.open(ARCHIVE, "r") as archive:
        members = set(archive.getnames())
        missing = [row["file"] for row in samples if row["file"] not in members]
        if missing:
            raise FileNotFoundError(f"Missing expected GEO archive members: {missing}")
        for row in samples:
            member = archive.extractfile(row["file"])
            if member is None:
                raise OSError(f"Could not read archive member {row['file']}")
            with gzip.open(io.BytesIO(member.read()), "rt") as stream:
                pairs = [line.rstrip("\n").split("\t") for line in stream]
            ids = [pair[0] for pair in pairs]
            if gene_ids is None:
                gene_ids = ids
            elif ids != gene_ids:
                raise ValueError(f"Gene ID order differs in {row['file']}")
            columns.append([pair[1] for pair in pairs])
    if gene_ids is None:
        raise ValueError("The sample manifest is empty")
    with COUNTS_OUT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *[row["sample_id"] for row in samples]])
        for gene_index, gene_id in enumerate(gene_ids):
            writer.writerow([gene_id, *[column[gene_index] for column in columns]])
    print(f"Wrote {len(gene_ids):,} genes x {len(samples)} samples: {COUNTS_OUT}")


if __name__ == "__main__":
    main()
