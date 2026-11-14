import gzip
import hashlib
from dataclasses import dataclass, replace
from pathlib import Path

import pandas as pd

TCGA_PREFIX = "TCGA-"
DATA_SUFFIXES = (".tsv", ".tsv.gz", ".csv", ".csv.gz")


def content_sha256(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    h = hashlib.sha256()
    with opener(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _mostly_tcga(values, threshold=0.9):
    s = pd.Series(values).astype(str)
    return len(s) > 0 and s.str.startswith(TCGA_PREFIX).mean() >= threshold

def _pick_id_column(raw):
    best, best_score = None, None
    for c in raw.columns:
        s = raw[c].astype(str)
        if len(s) == 0 or s.str.startswith(TCGA_PREFIX).mean() < 0.9:
            continue
        name = str(c).lower()
        score = (
            int((s.str.len() >= 15).mean() >= 0.9),
            s.nunique() / len(s),
            int(name in ("sampleid", "sample", "sample_id", "sample_barcode")),
        )
        if best_score is None or score > best_score:
            best, best_score = c, score
    return best

@dataclass
class Layer:
    name: str
    path: Path
    layout: str
    id_field: str
    feature_field: str
    table: pd.DataFrame
    sha256: str


def read_table(path):
    sep = "," if ".csv" in path.name.lower() else "\t"
    try:
        return pd.read_csv(path, sep=sep, engine="pyarrow")
    except Exception:
        return pd.read_csv(path, sep=sep, low_memory=False)


def load_layer(path, sha=None):
    path = Path(path)
    raw = read_table(path)
    tcga_cols = [c for c in raw.columns if str(c).startswith(TCGA_PREFIX)]
    if len(tcga_cols) > 10 and len(tcga_cols) >= 0.5 * (raw.shape[1] - 1):
        feature_field = str(raw.columns[0])
        table = raw.set_index(raw.columns[0])[tcga_cols].T
        layout, id_field = "wide", "columns"
    else:
        id_field = _pick_id_column(raw)
        if id_field is None:
            return None
        table = raw.set_index(id_field)
        layout, feature_field = "long", "columns"
    table.index = table.index.astype(str)
    table.index.name = "sample_id"
    table.columns = [str(c) for c in table.columns]
    return Layer(
        name=path.name,
        path=path,
        layout=layout,
        id_field=str(id_field),
        feature_field=feature_field,
        table=table,
        sha256=sha or content_sha256(path),
    )


def load_release(folder, reuse=None, skip_files=None):
    folder = Path(folder)
    reuse = reuse or {}
    skip = set(skip_files or [])
    layers = {}
    for path in sorted(folder.iterdir()):
        if not (path.is_file() and path.name.lower().endswith(DATA_SUFFIXES)):
            continue
        if path.name in skip:
            continue
        sha = content_sha256(path)
        prev = reuse.get(path.name)
        if prev is not None and prev.sha256 == sha:
            layers[path.name] = replace(prev, path=path)
            continue
        layer = load_layer(path, sha)
        if layer is not None:
            layers[path.name] = layer
    return layers