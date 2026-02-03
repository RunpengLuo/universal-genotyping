import os
import sys
import subprocess
from io import StringIO
from collections import OrderedDict

import pandas as pd
import numpy as np


def symlink_force(src, dst):
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    os.symlink(os.path.abspath(src), os.path.abspath(dst))


def maybe_path(x):
    if x == [] or x is None:
        return None
    return x


def get_chr2ord(ch):
    chr2ord = {}
    for i in range(1, 23):
        chr2ord[f"{ch}{i}"] = i
    chr2ord[f"{ch}X"] = 23
    chr2ord[f"{ch}Y"] = 24
    chr2ord[f"{ch}M"] = 25
    return chr2ord


def sort_chroms(chromosomes: list):
    assert len(chromosomes) != 0
    chromosomes = [str(c) for c in chromosomes]
    ch = "chr" if str(chromosomes[0]).startswith("chr") else ""
    chr2ord = get_chr2ord(ch)
    return sorted(chromosomes, key=lambda x: chr2ord[x])


def sort_df_chr(df: pd.DataFrame, ch="#CHR", pos="POS"):
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df
