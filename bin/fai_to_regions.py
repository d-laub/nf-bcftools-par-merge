#! /usr/bin/env python3
from pathlib import Path

from cyclopts import run


def main(fai: Path, min_windows: int, out: Path, primary_assembly: bool = False):
    import re

    import polars as pl
    from natsort import natsorted

    if fai.suffix != ".fai":
        raise ValueError(f"Fai file must have .fai suffix: {fai}")
    if not re.search(r"\.bed(\.gz)?$", str(out.name)):
        raise ValueError(f"Output file must have .bed or .bed.gz suffix: {out}")

    bed = pl.read_csv(
        fai,
        separator="\t",
        has_header=False,
        columns=range(2),
        new_columns=["chrom", "end"],
    ).with_columns(start=0)

    if primary_assembly:
        primary_chroms = list(map(str, range(1, 23))) + ["X", "Y", "M", "MT"]
        chr_prefixed = bed["chrom"].str.contains("^chr").any()
        if chr_prefixed:
            primary_chroms = ["chr" + chrom for chrom in primary_chroms]
        bed = bed.filter(pl.col("chrom").is_in(primary_chroms))

    total_length = int(bed["end"].sum())
    window_size = -(-total_length // min_windows)

    chroms = natsorted(bed["chrom"].unique())

    (
        bed.with_columns(
            borders=pl.int_ranges(
                end=pl.col("end") + 1,
                step=pl.col("end") // -(-pl.col("end") // window_size),
            )
        )
        .with_columns(
            start=pl.col("borders").list.head(pl.col("borders").list.len() - 1),
            end=pl.col("borders").list.slice(1),
        )
        .explode("start", "end")
        .select(pl.col("chrom").cast(pl.Enum(chroms)), "start", "end")
        .sort("chrom", "start")
        .write_csv(out, separator="\t", include_header=False)
    )


if __name__ == "__main__":
    run(main)
