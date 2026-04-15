"""Assay-type, data-type, and reference constants shared by Snakefile and scripts."""

import os

NONBULK_ASSAYS = {"scATAC", "scRNA", "VISIUM", "VISIUM3prime"}
BULK_ASSAYS = {"bulkWGS", "bulkWGS-lr", "bulkWES"}
ALLOWED_ASSAY_TYPES = list(BULK_ASSAYS) + list(NONBULK_ASSAYS)
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}
CHROM_ORDER = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

ASSAY_TYPE2MODALITY = {
    "bulkWGS": "DNA",
    "bulkWGS-lr": "DNA",
    "bulkWES": "DNA",
    "scATAC": "DNA",
    "scRNA": "RNA",
    "VISIUM": "RNA",
    "VISIUM3prime": "RNA",
}

ASSAY_TYPE2FEATURE_TYPE = {
    "bulkWGS": "dna",
    "bulkWGS-lr": "dna",
    "bulkWES": "exon",
    "scATAC": "tile",
    "scRNA": "gene",
    "VISIUM": "gene",
    "VISIUM3prime": "gene",
}

ALLOWED_REFVERS = ["hg19", "hg38", "chm13v2"]

GTF_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes",
]


def get_eagle_gmap_path(phaser_dir, refvers):
    """Return the genetic map path for Eagle2."""
    gmap_dir = os.path.join(phaser_dir, "tables")
    return gmap_dir, lambda chrname: os.path.join(
        gmap_dir, f"genetic_map_{refvers}_withX.txt.gz"
    )


def get_shapeit_gmap_path(phaser_dir, refvers):
    """Return the genetic map directory and per-chromosome path function for SHAPEIT5."""
    _gmap_patterns = {
        "hg19": ("b37", "b37/chr{chrname}.b37.gmap.gz"),
        "hg38": ("b38", "b38/chr{chrname}.b38.gmap.gz"),
        "chm13v2": ("chm13v2", "chm13v2/chr{chrname}.t2t.scaled.gmap.gz"),
    }
    subdir, pattern = _gmap_patterns[refvers]
    gmap_dir = os.path.join(phaser_dir, "resources/maps", subdir)
    return gmap_dir, lambda chrname: os.path.join(
        phaser_dir, "resources/maps", pattern.format(chrname=chrname)
    )


def get_phasing_panel_path(phasing_panel):
    """Return a per-chromosome phasing panel path function."""
    return lambda chrname: os.path.join(phasing_panel, f"chr{chrname}.genotypes.bcf")
