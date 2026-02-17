import os
import sys

import pandas as pd

"""
Take a barcode/label visium annotation file, convert to correct format.
"""

_, ann_file, out_file = sys.argv

def simplify_label(v):
    if v[0] == "T":
        if "_" not in v:
            return "tumor"
        else:
            return v[str(v).find("_") + 1 :]
    else:
        return v

# change accordingly
barcode_col = "Barcode"
raw_label = "Microregion_annotation"
target_label = "path_label"

anns = pd.read_table(ann_file, sep="\t", keep_default_na=True)
print(anns.head())

anns = anns.rename(columns={barcode_col: "BARCODE"})

# NA labels are non-tumor
anns[raw_label].fillna("normal", inplace=True)
anns[target_label] = anns.apply(
    func=lambda r: simplify_label(r[raw_label]), axis=1
).astype("str")

print(anns.head())

anns.to_csv(out_file, sep="\t", header=True, index=False)
