#!/usr/bin/env python3
"""
Fetch current GO annotations for Aedes aegypti (taxid 7159) NCBI Gene IDs.

Strategy: rather than hitting NCBI E-utilities 18,727 times (slow, rate-limited,
and easy to get throttled/banned), this pulls NCBI's bulk gene2go file, which is
updated daily and contains every gene->GO annotation for every organism NCBI
tracks. We filter it down to Aedes aegypti and then to your gene ID list.
gene2go already includes GO term names, evidence codes, and category
(Function/Process/Component), so no extra API calls are needed.

Usage:
    python get_go_annotations.py --ids my_gene_ids.txt --out go_annotations.csv

    --ids can be:
      - a .txt file with one NCBI GeneID per line
      - a .csv/.tsv file with a column named 'GeneID' (or pass --id-column)
"""

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import pandas as pd

GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
AEDES_AEGYPTI_TAXID = 7159


def download_gene2go(dest: Path) -> Path:
    if dest.exists():
        print(f"Using cached {dest} (delete it to force a fresh download)")
        return dest
    print(f"Downloading {GENE2GO_URL} (this file is large, ~300MB, please wait)...")
    urllib.request.urlretrieve(GENE2GO_URL, dest)
    print("Download complete.")
    return dest


def normalize_id(raw: str) -> str:
    """Strip a leading 'LOC' from NCBI-style placeholder symbols.

    When a gene has no official symbol, NCBI names it LOC<GeneID>
    (e.g. LOC110673977 -> GeneID 110673977). gene2go is keyed on the
    bare numeric GeneID, so we need the digits only.
    """
    raw = raw.strip()
    if raw.upper().startswith("LOC") and raw[3:].isdigit():
        return raw[3:]
    return raw


def load_gene_ids(path: Path, id_column: str) -> tuple:
    if path.suffix.lower() == ".txt":
        raw_ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    else:
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(path, sep=sep, dtype=str)
        if id_column not in df.columns:
            sys.exit(f"Column '{id_column}' not found in {path}. Columns are: {list(df.columns)}")
        raw_ids = df[id_column].dropna().tolist()
    raw_ids = [str(i).strip() for i in raw_ids if str(i).strip()]

    # map normalized GeneID -> original string the user supplied, so we can
    # report back in their original naming (e.g. LOC110673977) as well.
    id_map = {normalize_id(r): r for r in raw_ids}
    ids = set(id_map.keys())
    n_stripped = sum(1 for r in raw_ids if r.upper().startswith("LOC") and r[3:].isdigit())
    if n_stripped:
        print(f"Stripped 'LOC' prefix from {n_stripped} IDs to get numeric GeneIDs")
    print(f"Loaded {len(ids)} unique gene IDs from {path}")
    return ids, id_map


def extract_annotations(gene2go_gz: Path, gene_ids: set, id_map: dict) -> pd.DataFrame:
    print(f"Scanning gene2go for tax_id {AEDES_AEGYPTI_TAXID} (Aedes aegypti)...")
    cols = ["tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier",
            "GO_term", "PubMed", "Category"]
    chunks = []
    with gzip.open(gene2go_gz, "rt") as fh:
        reader = pd.read_csv(fh, sep="\t", names=cols, header=0, dtype=str,
                              chunksize=500_000, na_values=["-"])
        for chunk in reader:
            chunk = chunk[chunk["tax_id"] == str(AEDES_AEGYPTI_TAXID)]
            if not chunk.empty:
                chunks.append(chunk)

    if not chunks:
        sys.exit("No Aedes aegypti rows found in gene2go — check the download.")

    aa = pd.concat(chunks, ignore_index=True)
    print(f"Found {aa['GeneID'].nunique()} Aedes aegypti genes with GO annotations in gene2go total.")

    matched = aa[aa["GeneID"].isin(gene_ids)].copy()
    found_ids = set(matched["GeneID"].unique())
    missing = gene_ids - found_ids

    print(f"Matched {len(found_ids)} / {len(gene_ids)} of your gene IDs.")
    if missing:
        print(f"{len(missing)} gene IDs had no GO annotation in gene2go "
              f"(either no annotations exist yet, or the ID is wrong/retired).")

    matched = matched.drop(columns=["tax_id"])
    matched.insert(0, "InputID", matched["GeneID"].map(id_map))
    return matched, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, type=Path, help="File with your NCBI Gene IDs")
    ap.add_argument("--id-column", default="GeneID", help="Column name if --ids is csv/tsv")
    ap.add_argument("--out", default="go_annotations.csv", type=Path, help="Output CSV path")
    ap.add_argument("--missing-out", default="genes_without_go.txt", type=Path,
                     help="Where to write gene IDs that had no GO annotation")
    ap.add_argument("--cache", default="gene2go.gz", type=Path,
                     help="Local path to cache the downloaded gene2go file")
    args = ap.parse_args()

    gene_ids, id_map = load_gene_ids(args.ids, args.id_column)
    gene2go_path = download_gene2go(args.cache)
    annotations, missing = extract_annotations(gene2go_path, gene_ids, id_map)

    annotations.to_csv(args.out, index=False)
    print(f"Wrote {len(annotations)} annotation rows to {args.out}")

    if missing:
        missing_original = sorted(id_map[m] for m in missing)
        args.missing_out.write_text("\n".join(missing_original) + "\n")
        print(f"Wrote {len(missing)} unannotated gene IDs to {args.missing_out}")


if __name__ == "__main__":
    main()
