import os, sys, gzip, logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import read_VCF, sort_df_chr

##################################################
"""
Parse genetic map files
"""

logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

snps = read_VCF(sm.input["phased_vcf"])
gmap_files = list(sm.input["gmap_files"])
chrnames = list(sm.params["chrnames"])
phaser = sm.params["phaser"]

logging.info(f"parse genetic map files, phaser={phaser}")

required_columns = ["#CHR", "POS", "cM"]
if phaser == "eagle":
    if len(set(gmap_files)) == 1:
        # commonly eagle genetic map is single file covers all chroms.
        genetic_map = pd.read_table(gmap_files[0], sep=" ", index_col=None).rename(
            columns={
                "chr": "#CHR",
                "position": "POS",
                "COMBINED_rate(cM/Mb)": "recomb_rate",
                "Genetic_Map(cM)": "cM",
            }
        )
    else:
        genetic_maps = []
        for chrname, gmap_file in zip(chrnames, gmap_files):
            genetic_map = pd.read_table(gmap_file, sep=" ", index_col=None).rename(
                columns={
                    "chr": "#CHR",
                    "position": "POS",
                    "COMBINED_rate(cM/Mb)": "recomb_rate",
                    "Genetic_Map(cM)": "cM",
                }
            )
            genetic_maps.append(genetic_map)
        genetic_map = pd.concat(genetic_maps, ignore_index=True)

    genetic_map["#CHR"] = genetic_map["#CHR"].astype(str)
    genetic_map.loc[genetic_map["#CHR"] == "23", "#CHR"] = "X"
    if not str(genetic_map["#CHR"].loc[0]).startswith("chr"):
        if str(snps["#CHR"].iloc[0]).startswith("chr"):
            genetic_map["#CHR"] = "chr" + genetic_map["#CHR"].astype(str)
    genetic_map = genetic_map[
        genetic_map["#CHR"].isin({f"chr{c}" for c in chrnames})
    ].reset_index(drop=True)
    genetic_map = sort_df_chr(genetic_map, ch="#CHR", pos="POS")
    genetic_map[required_columns].to_csv(
        sm.output["gmap_tsv"], sep="\t", header=True, index=False
    )

if phaser == "shapeit":
    genetic_maps = []
    for chrname, gmap_file in zip(chrnames, gmap_files):
        genetic_map = pd.read_csv(
            gmap_file,
            sep="\t",
            comment="#",
        )
        assert "pos" in genetic_map.columns and "cM" in genetic_map.columns, (
            f"gmap.gz file is invalid"
        )
        genetic_map["#CHR"] = f"chr{chrname}"
        genetic_map["POS"] = genetic_map["pos"]

        genetic_maps.append(genetic_map[["#CHR", "POS", "cM"]].reset_index(drop=True))

    genetic_map = pd.concat(genetic_maps, ignore_index=True)
    genetic_map = sort_df_chr(genetic_map, ch="#CHR", pos="POS")
    genetic_map[required_columns].to_csv(
        sm.output["gmap_tsv"], sep="\t", header=True, index=False
    )
