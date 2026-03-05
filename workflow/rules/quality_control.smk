##################################################
# Bin filtering rules (bulk only)
# Removes bins with NaN RDR, extreme GC, or low mappability
##################################################

if workflow_mode == "bulk_genotyping":

    _rdr_method = config.get("rdr_method", "bin")
    _filter_cfg = config.get("params_filter_bins", {})

    rule qc_filter_bulk_bb:
        input:
            bb_file=config["bb_dir"] + "/{assay_type}/raw/bb.raw.tsv.gz",
            rdr_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.rdr.npz",
            dp_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.depth.npz",
            tot_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.Tallele.npz",
            a_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.Aallele.npz",
            b_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.Ballele.npz",
            baf_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.baf.npz",
            corr_factors=config["bb_dir"] + "/{assay_type}/raw/bb.raw.corr_factors.tsv.gz",
            sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
        output:
            bb_file=config["bb_dir"] + "/{assay_type}/bb.tsv.gz",
            rdr_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.rdr.npz",
            dp_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.depth.npz",
            tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
            a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
            b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
            baf_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.baf.npz",
        wildcard_constraints:
            assay_type="(bulkDNA|bulkWGS|bulkWES)",
        params:
            gc_lo=_filter_cfg.get("gc_lo", 0.01),
            gc_hi=_filter_cfg.get("gc_hi", 0.99),
            map_tol=_filter_cfg.get("map_tol", 0.2),
            rdr_hi=_filter_cfg.get("rdr_hi", 0.999),
        log:
            config["log_dir"] + "/qc_filter_bulk_bb.{assay_type}.log",
        run:
            import numpy as np
            import pandas as pd
            import logging
            from scripts.utils import setup_logging

            setup_logging(log[0])

            bb = pd.read_csv(input.bb_file, sep="\t")
            rdr = np.load(input.rdr_mtx_bb)["mat"]
            dp = np.load(input.dp_mtx_bb)["mat"]
            tot = np.load(input.tot_mtx_bb)["mat"]
            a_mat = np.load(input.a_mtx_bb)["mat"]
            b_mat = np.load(input.b_mtx_bb)["mat"]
            baf = np.load(input.baf_mtx_bb)["mat"]
            corr = pd.read_csv(input.corr_factors, sep="\t")

            n = len(bb)
            mask = np.ones(n, dtype=bool)

            nan_rdr = np.any(np.isnan(rdr), axis=1)
            mask &= ~nan_rdr
            logging.info(f"NaN RDR filter: {nan_rdr.sum()}/{n}")

            max_rdr = np.nanmax(rdr, axis=1)
            rdr_hi_val = np.nanquantile(max_rdr[mask], params.rdr_hi)
            bad_rdr = max_rdr > rdr_hi_val
            mask &= ~bad_rdr
            logging.info(f"outlier RDR filter (>{params.rdr_hi} quantile = {rdr_hi_val:.3f}): {bad_rdr.sum()}/{n}")

            gc = corr["GC"].to_numpy()
            gc_lo_val = np.nanquantile(gc[mask], params.gc_lo)
            gc_hi_val = np.nanquantile(gc[mask], params.gc_hi)
            bad_gc = (gc < gc_lo_val) | (gc > gc_hi_val)
            mask &= ~bad_gc
            logging.info(f"extreme GC filter: {bad_gc.sum()}/{n}")

            if "MAP" in corr.columns:
                mapp = corr["MAP"].to_numpy()
                for t in [0.1, 0.2, 0.3, 0.5]:
                    logging.info(f"  MAP < {t}: {(mapp < t).sum()}/{n}")
                bad_map = mapp < params.map_tol
                mask &= ~bad_map
                logging.info(f"low MAP filter (tol={params.map_tol}): {bad_map.sum()}/{n}")

            logging.info(f"kept {mask.sum()}/{n} bins after all filters")

            bb[mask].to_csv(output.bb_file, sep="\t", index=False, compression="gzip")
            np.savez_compressed(output.rdr_mtx_bb, mat=rdr[mask])
            np.savez_compressed(output.dp_mtx_bb, mat=dp[mask])
            np.savez_compressed(output.tot_mtx_bb, mat=tot[mask])
            np.savez_compressed(output.a_mtx_bb, mat=a_mat[mask])
            np.savez_compressed(output.b_mtx_bb, mat=b_mat[mask])
            np.savez_compressed(output.baf_mtx_bb, mat=baf[mask])
