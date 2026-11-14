import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tcga-brca"
A = DATA / "raw"
B = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA / "releases" / "practice-B"
SAMPLES = sys.argv[2:] or ["TCGA-A8-A09K-01", "TCGA-B6-A0WX-01"]
PROTEINS = ["ERALPHA", "PR", "HER2", "AKT", "BCL2"]

clin = pd.read_csv(A / "clinical.tsv", sep="\t", low_memory=False).set_index("sampleID")
prot_a = pd.read_csv(A / "protein_rppa.tsv.gz", sep="\t", index_col=0)
prot_b = pd.read_csv(B / "protein_rppa.tsv.gz", sep="\t", index_col=0)

print(f"protein samples: A={prot_a.shape[1]} B={prot_b.shape[1]}")
print()
for s in SAMPLES:
    er = clin.at[s, "breast_carcinoma_estrogen_receptor_status"] if s in clin.index else None
    pam = clin.at[s, "PAM50Call_RNAseq"] if s in clin.index else None
    print(f"{s}  ER={er}  PAM50={pam}")
    for p in PROTEINS:
        va = prot_a.at[p, s] if s in prot_a.columns and p in prot_a.index else None
        vb = prot_b.at[p, s] if s in prot_b.columns and p in prot_b.index else None
        print(f"  {p:<8} A={va}  B={vb}")
    print()