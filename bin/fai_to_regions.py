#! /usr/bin/env python3
from pathlib import Path

from cyclopts import run


def main(fai: Path, min_windows: int, out: Path, primary_assembly: bool = False):
    import re

    import polars as pl
    import pyranges as pr

    if fai.suffix != ".fai":
        raise ValueError(f"Fai file must have .fai suffix: {fai}")
    if not re.search(r"\.bed(\.gz)?$", str(out.name)):
        raise ValueError(f"Output file must have .bed or .bed.gz suffix: {out}")

    bed = pl.read_csv(
        fai,
        separator="\t",
        has_header=False,
        columns=range(2),
        new_columns=["Chromosome", "End"],
    ).with_columns(Start=0)

    if primary_assembly:
        primary_chroms = list(map(str, range(1, 23))) + ["X", "Y", "M", "MT"]
        chr_prefixed = bed["Chromosome"].str.contains("^chr").any()
        if chr_prefixed:
            primary_chroms = ["chr" + chrom for chrom in primary_chroms]
        bed = bed.filter(pl.col("Chromosome").is_in(primary_chroms))

    total_length = int(bed["End"].sum())
    window_size = -(-total_length // min_windows)

    (
        pr.PyRanges(bed.to_pandas())
        .window_ranges(window_size=window_size)
        .sort_ranges()
        .to_bed(str(out))
    )


if __name__ == "__main__":
    run(main)
