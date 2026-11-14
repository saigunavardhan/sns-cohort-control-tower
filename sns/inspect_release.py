import gzip
import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "tcga-brca"
RELEASES = {"A": BASE / "raw", "B": BASE / "releases" / "practice-B"}


def content_hash(path):
    opener = gzip.open if path.suffix == ".gz" else open
    h = hashlib.sha256()
    with opener(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def profile(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    tcga_cols = [c for c in df.columns if str(c).startswith("TCGA-")]
    if len(tcga_cols) > 10:
        layout = "wide"
        id_field = "columns"
        ids = pd.Series(tcga_cols, dtype=str)
    else:
        layout = "long"
        id_field = next(
            (c for c in df.columns
             if df[c].astype(str).str.startswith("TCGA-").mean() > 0.9),
            None,
        )
        ids = df[id_field].astype(str) if id_field else pd.Series(dtype=str)
    lengths = ids.str.len().value_counts().sort_index()
    return {
        "rows": df.shape[0],
        "cols": df.shape[1],
        "layout": layout,
        "id_field": id_field,
        "n_ids": int(ids.size),
        "n_unique_ids": int(ids.nunique()),
        "id_lengths": {int(k): int(v) for k, v in lengths.items()},
        "first_col": str(df.columns[0]),
        "first_values": df.iloc[:3, 0].astype(str).tolist(),
        "head_cols": [str(c) for c in df.columns[:6]],
        "hash": content_hash(path),
    }


def main():
    names = sorted(
        {p.name for p in RELEASES["A"].glob("*.tsv*")}
        | {p.name for p in RELEASES["B"].glob("*.tsv*")}
    )
    results = {}
    for tag, folder in RELEASES.items():
        for name in names:
            path = folder / name
            if path.exists():
                results[(tag, name)] = profile(path)

    for (tag, name), r in results.items():
        print(f"[{tag}] {name}")
        for k, v in r.items():
            if k != "hash":
                print(f"  {k}: {v}")

    print("[compare]")
    for name in names:
        a = results.get(("A", name))
        b = results.get(("B", name))
        if a is None or b is None:
            print(f"  {name}: missing_in={'A' if a is None else 'B'}")
            continue
        status = "same" if a["hash"] == b["hash"] else "different"
        print(
            f"  {name}: content={status} "
            f"shape_A={a['rows']}x{a['cols']} shape_B={b['rows']}x{b['cols']} "
            f"unique_ids_A={a['n_unique_ids']} unique_ids_B={b['n_unique_ids']}"
        )


if __name__ == "__main__":
    main()