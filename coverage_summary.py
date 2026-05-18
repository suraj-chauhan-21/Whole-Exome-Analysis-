#!/usr/bin/env python3
"""
coverage_summary.py
===================
Parse a samtools depth output file and produce a coverage summary table
at multiple depth thresholds.

Usage:
    python scripts/coverage_summary.py \
        --depth   results/metrics/SRR22317682.depth.txt \
        --sample  SRR22317682 \
        --out     results/metrics/SRR22317682.coverage_summary.tsv \
        --thresholds 1 10 20 30 50 100

Output TSV columns:
    sample | threshold | covered_bases | mean_depth | pct_covered
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth",       required=True,  help="samtools depth output (3-col: chrom pos depth)")
    p.add_argument("--sample",      required=True,  help="Sample ID for output label")
    p.add_argument("--out",         required=True,  help="Output TSV file path")
    p.add_argument("--thresholds",  nargs="+", type=int,
                   default=[1, 10, 20, 30, 50, 100],
                   help="Depth thresholds to evaluate (default: 1 10 20 30 50 100)")
    return p.parse_args()


def load_depth(path: str) -> list[int]:
    """Read the third column (depth) from a samtools depth file."""
    depths = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            depths.append(int(parts[2]))
    return depths


def summarise(depths: list[int], thresholds: list[int], sample: str) -> list[dict]:
    total_bases = len(depths)
    rows = []
    for t in sorted(thresholds):
        passing = [d for d in depths if d >= t]
        covered = len(passing)
        mean_d  = sum(passing) / covered if covered else 0.0
        pct     = 100.0 * covered / total_bases if total_bases else 0.0
        rows.append({
            "sample":        sample,
            "threshold_x":   f"{t}x",
            "covered_bases": covered,
            "total_bases":   total_bases,
            "pct_covered":   round(pct, 2),
            "mean_depth":    round(mean_d, 2),
        })
    return rows


def write_tsv(rows: list[dict], out_path: str):
    headers = ["sample", "threshold_x", "covered_bases",
               "total_bases", "pct_covered", "mean_depth"]
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[h]) for h in headers) + "\n")


def main():
    args = parse_args()

    print(f"[coverage_summary] Loading depth file: {args.depth}", file=sys.stderr)
    depths = load_depth(args.depth)

    if not depths:
        print("[coverage_summary] ERROR: No depth data found.", file=sys.stderr)
        sys.exit(1)

    print(f"[coverage_summary] {len(depths):,} positions loaded. "
          f"Computing thresholds: {args.thresholds}", file=sys.stderr)

    rows = summarise(depths, args.thresholds, args.sample)
    write_tsv(rows, args.out)

    # Also print to stdout for quick inspection
    print(f"{'Sample':<20} {'Threshold':<12} {'Covered':<15} "
          f"{'Total':<15} {'% Covered':<12} {'Mean Depth'}")
    print("-" * 80)
    for r in rows:
        print(f"{r['sample']:<20} {r['threshold_x']:<12} {r['covered_bases']:<15,} "
              f"{r['total_bases']:<15,} {r['pct_covered']:<12} {r['mean_depth']}")

    print(f"\n[coverage_summary] Summary written to: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
