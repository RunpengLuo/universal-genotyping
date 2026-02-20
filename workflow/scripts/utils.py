import os
import sys
import logging
import subprocess
from io import StringIO
from collections import OrderedDict

import pandas as pd
import numpy as np

SPOT_ASSAYS = {"VISIUM"}
CELL_ASSAYS = {"scRNA", "scATAC", "multiome"}
GEX_ASSAYS = {"scRNA", "multiome", "VISIUM"}
ATAC_ASSAYS = {"scATAC", "multiome"}
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}
ALL_ASSAYS = ["scRNA", "scATAC", "multiome", "VISIUM"]

def setup_logging(log):
    """Configure root logger to write INFO-level messages with timestamps to a file.

    Parameters
    ----------
    log : str
        Path to the log file.
    """
    logging.basicConfig(
        filename=log,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )


def symlink_force(src, dst):
    """Create a symlink from *src* to *dst*, removing any existing file at *dst*.

    Parameters
    ----------
    src : str
        Source path (symlink target).
    dst : str
        Destination path (symlink location).
    """
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    os.symlink(os.path.abspath(src), os.path.abspath(dst))


def maybe_path(x):
    """Return None if *x* is an empty list or None, otherwise return *x* unchanged.

    Useful for coercing Snakemake optional inputs to ``None``.
    """
    if x == [] or x is None:
        return None
    return x


def sort_chroms(chromosomes: list):
    """Sort chromosome names in standard genomic order (1-22, X, Y, M).

    Handles both ``chr``-prefixed and plain chromosome names.

    Parameters
    ----------
    chromosomes : list of str
        Chromosome name strings.

    Returns
    -------
    list of str
        Sorted chromosome names.
    """
    assert len(chromosomes) != 0
    chromosomes = [str(c) for c in chromosomes]
    ch = "chr" if str(chromosomes[0]).startswith("chr") else ""
    chr2ord = {}
    for i in range(1, 23):
        chr2ord[f"{ch}{i}"] = i
    chr2ord[f"{ch}X"] = 23
    chr2ord[f"{ch}Y"] = 24
    chr2ord[f"{ch}M"] = 25
    return sorted(chromosomes, key=lambda x: chr2ord[x])


def sort_df_chr(df: pd.DataFrame, ch="#CHR", pos="POS"):
    """Sort a DataFrame by chromosome (genomic order) then by position, in-place.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with chromosome and position columns.
    ch : str
        Name of the chromosome column.
    pos : str
        Name of the position column.

    Returns
    -------
    pd.DataFrame
        The same DataFrame, sorted in-place.
    """
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df
