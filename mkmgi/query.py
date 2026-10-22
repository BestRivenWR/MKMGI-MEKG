"""Read an archived prediction row without changing or recalculating its values."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import re


def read_gene(table, gene):
    if not re.fullmatch(r"Zm00001eb\d{6}", gene):
        raise ValueError("Use an exact B73 v5 gene ID (Zm00001eb plus six digits).")
    table = Path(table)
    opener = gzip.open if table.suffix == ".gz" else open
    found = []
    with opener(table, "rt", encoding="utf-8-sig", newline="") as f:
        rows = csv.DictReader(f, delimiter="\t")
        required = {"gene_id", "priority_score", "priority_percentile"}
        if not required.issubset(rows.fieldnames or []):
            raise ValueError("Missing required prediction columns.")
        for row in rows:
            if row["gene_id"] == gene:
                found.append(row)
    if len(found) > 1:
        raise ValueError("Duplicate gene ID in prediction table.")
    if not found:
        return {"gene_id": gene, "supported_by_table": False}
    row = found[0]
    for field in ("priority_score", "priority_percentile"):
        value = float(row[field])
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{field} must be a finite value between zero and one.")
    return {
        "gene_id": gene,
        "supported_by_table": True,
        "score_semantics": "candidate priority, not phenotype probability",
        "prediction": row,  # strings preserve the original stored precision
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    args = parser.parse_args()
    print(json.dumps(read_gene(args.table, args.gene), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
