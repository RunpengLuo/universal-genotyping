"""Assay-type, data-type, and reference constants shared by Snakefile and scripts."""

import os

NONBULK_ASSAYS = {"scATAC", "scRNA", "VISIUM", "VISIUM3prime"}
BULK_ASSAYS = {"bulkWGS", "bulkWGS-lr", "bulkWES"}
ALLOWED_ASSAY_TYPES = list(BULK_ASSAYS) + list(NONBULK_ASSAYS)
SPATIAL_ASSAYS = {"VISIUM", "VISIUM3prime"}

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


def get_eagle_gmap_path(phaser_dir, refvers):
    """Return the genetic map path for Eagle2."""
    gmap_dir = os.path.join(phaser_dir, "tables")
    return gmap_dir, lambda chrname: os.path.join(
        gmap_dir, f"genetic_map_{refvers}_withX.txt.gz"
    )


def get_shapeit_gmap_path(phaser_dir, refvers, gmap_dir=None):
    """Return the genetic map directory and per-chromosome path function for SHAPEIT5."""
    if refvers == "chm13v2":
        assert gmap_dir is not None, "gmap_dir required for chm13v2 shapeit phasing"
        return gmap_dir, lambda chrname: os.path.join(
            gmap_dir, f"chr{chrname}.t2t.scaled.gmap.gz"
        )
    _refvers2 = {"hg19": "b37", "hg38": "b38"}
    refvers2 = _refvers2[refvers]
    gmap_dir = os.path.join(phaser_dir, "resources/maps")
    return gmap_dir, lambda chrname: os.path.join(
        gmap_dir, f"{refvers2}/chr{chrname}.{refvers2}.gmap.gz"
    )


def get_phasing_panel_path(phasing_panel):
    """Return a per-chromosome phasing panel path function."""
    return lambda chrname: os.path.join(
        phasing_panel, f"chr{chrname}.genotypes.bcf"
    )
