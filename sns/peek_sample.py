import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "tcga-brca" / "raw"
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "TCGA-A8-A09K-01"

LONG = {
    "clinical.tsv": ("sampleID", [
        "sample_type",
        "gender",
        "Age_at_Initial_Pathologic_Diagnosis_nature2012",
        "breast_carcinoma_estrogen_receptor_status",
        "breast_carcinoma_progesterone_receptor_status",
        "PAM50Call_RNAseq",
    ]),
    "survival.tsv": ("sample", ["_PATIENT", "OS", "OS.time"]),
}
WIDE = {
    "protein_rppa.tsv.gz": ["ERALPHA", "PR", "HER2", "AKT"],
    "rna_hiseq.tsv.gz": ["ESR1", "PGR", "ERBB2", "GATA3"],
    "mutation_gene_level.tsv.gz": ["TP53", "PIK3CA", "CDH1", "GATA3"],
    "copy_number_gistic.tsv.gz": ["ERBB2", "MYC", "CCND1"],
}

print(f"sample: {SAMPLE}")
print(f"patient: {SAMPLE[:12]}")
print(f"sample_type_code: {SAMPLE[13:15]}")
print()

for name, (id_col, cols) in LONG.items():
    df = pd.read_csv(RAW / name, sep="\t", low_memory=False).set_index(id_col)
    print(f"[{name}] samples={df.shape[0]} columns={df.shape[1]}")
    if SAMPLE in df.index:
        keep = [c for c in cols if c in df.columns]
        for c in keep:
            print(f"  {c}: {df.at[SAMPLE, c]}")
    else:
        print("  absent")
    print()

for name, feats in WIDE.items():
    header = pd.read_csv(RAW / name, sep="\t", nrows=0)
    first = header.columns[0]
    n_samples = sum(str(c).startswith("TCGA-") for c in header.columns)
    print(f"[{name}] samples={n_samples}")
    if SAMPLE not in header.columns:
        print("  absent")
        print()
        continue
    df = pd.read_csv(RAW / name, sep="\t", usecols=[first, SAMPLE]).set_index(first)
    for f in feats:
        if f in df.index:
            v = df.loc[f, SAMPLE]
            v = v.iloc[0] if hasattr(v, "iloc") else v
            print(f"  {f}: {v}")
    print()