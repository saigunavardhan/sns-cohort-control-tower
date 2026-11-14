import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from engine.release import load_layer  # noqa: E402

TEMPLATE = HERE / "dashboard_template.html"
RAW = ROOT / "data" / "tcga-brca" / "raw"

KIND_COLOR = {
    "sample_swap": "#ff5d73",
    "duplicate": "#ffb347",
    "membership_change": "#b18cff",
    "missing_data": "#5fb3ff",
    "label_conflict": "#3ddc97",
    "other": "#9aa5b1",
}
KIND_LABEL = {
    "sample_swap": "swapped",
    "duplicate": "copies",
    "membership_change": "lost or renamed",
    "missing_data": "values missing",
    "label_conflict": "labels changed",
    "other": "values changed",
}

_cache = {}


def table(folder, name):
    p = Path(folder) / name
    key = str(p)
    if key in _cache:
        return _cache[key]
    t = None
    if p.exists():
        layer = load_layer(p)
        if layer is not None:
            t = layer.table[~layer.table.index.duplicated()]
    _cache[key] = t
    return t


def numeric(t):
    if all(pd.api.types.is_numeric_dtype(d) for d in t.dtypes):
        return t.astype(float)
    return t.apply(pd.to_numeric, errors="coerce")


def r(v, n=4):
    return [None if (x is None or not np.isfinite(x)) else round(float(x), n) for x in v]


def evidence_by_file(doc):
    out = {}
    for f in doc.get("findings", []):
        for e in f.get("evidence", []):
            if e.get("kind") == "sample" and e.get("file"):
                out.setdefault(e["file"], {}).setdefault(f["kind"], set()).add(e["id"])
    return out


def release_files(doc):
    rels = doc.get("data", {}).get("releases", [])
    a = {x["file"]: x["sha256"] for x in rels[0]["files"]} if rels else {}
    b = {x["file"]: x["sha256"] for x in rels[1]["files"]} if len(rels) > 1 else {}
    return a, b


def clinical_name(files):
    names = [n for n in files if "clinical" in n.lower()]
    return names[0] if names else None


def main_file(doc):
    ev = evidence_by_file(doc)
    if not ev:
        return None
    return max(ev, key=lambda f: sum(len(s) for s in ev[f].values()))


def build_pipeline(doc, a_dir):
    fa, fb = release_files(doc)
    changed = [n for n in fb if n in fa and fa[n] != fb[n]]
    rows = 0
    for n in changed:
        t = table(a_dir, n)
        rows += 0 if t is None else len(t)
    samples = set()
    for kinds in evidence_by_file(doc).values():
        for s in kinds.values():
            samples |= s
    findings = doc.get("findings", [])
    subs = doc.get("sns_run", {}).get("subgroup_checks", [])
    analyses = doc.get("analyses", [])
    affected = [a["id"] for a in analyses if a.get("affected_by")]
    safe = [a["id"] for a in analyses if not a.get("affected_by")]
    restored = [a["id"] for a in analyses if a.get("affected_by") and a.get("repair_restores_before")]
    group = None
    if subs:
        h = subs[0]
        group = {"value": h.get("value"), "column": h.get("column"),
                 "a": h.get("with_data_in_a"), "b": h.get("with_data_in_b")}
    if not affected:
        proof = "n/a"
    elif len(restored) == len(affected):
        proof = "yes"
    else:
        proof = "partial"
    return {
        "files_total": len(fb),
        "files_changed": changed,
        "rows_checked": rows,
        "samples_changed": len(samples),
        "n_findings": len(findings),
        "kinds": sorted({KIND_LABEL.get(f.get("kind"), f.get("kind", "")) for f in findings}),
        "group": group,
        "affected": affected,
        "safe": safe,
        "restored": restored,
        "proof": proof,
    }


def build_flow(doc, a_dir, b_dir):
    f = main_file(doc)
    if not f:
        return None
    ta, tb = table(a_dir, f), table(b_dir, f)
    if ta is None or tb is None:
        return None
    ids_a, ids_b = set(ta.index), set(tb.index)
    ev = evidence_by_file(doc).get(f, {})
    lost, new, both = ids_a - ids_b, ids_b - ids_a, ids_a & ids_b
    groups = {k: v & both for k, v in ev.items() if v & both}
    flagged = set().union(*groups.values()) if groups else set()
    unchanged = both - flagged
    labels = [f"Release A: {len(ids_a)}", f"Release B: {len(ids_b)}", f"unchanged: {len(unchanged)}"]
    colors = ["#2f81f7", "#2f81f7", "#6e7681"]
    src, tgt, val = [0, 2], [2, 1], [len(unchanged), len(unchanged)]
    lcol = ["rgba(110,118,129,0.35)", "rgba(110,118,129,0.35)"]
    for k, v in groups.items():
        labels.append(f"{KIND_LABEL.get(k, k)}: {len(v)}")
        colors.append(KIND_COLOR.get(k, "#9aa5b1"))
        i = len(labels) - 1
        src += [0, i]
        tgt += [i, 1]
        val += [len(v), len(v)]
        lcol += ["rgba(255,93,115,0.55)", "rgba(255,93,115,0.55)"]
    if lost:
        labels.append(f"lost: {len(lost)}")
        colors.append(KIND_COLOR["membership_change"])
        src.append(0)
        tgt.append(len(labels) - 1)
        val.append(len(lost))
        lcol.append("rgba(177,140,255,0.6)")
    if new:
        labels.append(f"new IDs: {len(new)}")
        colors.append(KIND_COLOR["duplicate"])
        src.append(len(labels) - 1)
        tgt.append(1)
        val.append(len(new))
        lcol.append("rgba(255,179,71,0.65)")
    return {"file": f, "n_a": len(ids_a), "n_b": len(ids_b), "lost": len(lost), "new": len(new),
            "changed_same_id": len(flagged), "unchanged": len(unchanged),
            "labels": labels, "colors": colors, "source": src, "target": tgt, "value": val, "link_colors": lcol}


def color_column(doc, f, clin):
    for a in doc.get("analyses", []):
        if f in a.get("inputs", []):
            cols = a.get("uses_columns", {}).get(clin, [])
            if cols:
                return cols[0]
    return None


def build_pca(doc, a_dir, b_dir, max_features=2000):
    f = main_file(doc)
    if not f:
        return None
    ta, tb = table(a_dir, f), table(b_dir, f)
    if ta is None or tb is None:
        return None
    xa = numeric(ta)
    xa = xa.loc[:, xa.notna().mean() > 0.5]
    if xa.shape[1] < 3:
        return None
    if xa.shape[1] > max_features:
        xa = xa[xa.var().sort_values(ascending=False).index[:max_features]]
    mu = xa.mean()
    sd = xa.std().replace(0, 1).fillna(1)
    za = ((xa - mu) / sd).fillna(0).to_numpy()
    _, s, vt = np.linalg.svd(za, full_matrices=False)
    comps = vt[:3].T
    explained = (s[:3] ** 2 / (s ** 2).sum() * 100).round(1).tolist()
    xb = numeric(tb.reindex(columns=xa.columns))
    zb = ((xb - mu) / sd).fillna(0).to_numpy()
    pa = pd.DataFrame(za @ comps, index=xa.index, columns=["x", "y", "z"])
    pb = pd.DataFrame(zb @ comps, index=xb.index, columns=["x", "y", "z"])

    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    clin = table(a_dir, clin_name) if clin_name else None
    col = color_column(doc, f, clin_name) if clin is not None else None
    groups = []
    if col and col in clin.columns:
        lab = clin[col].reindex(pb.index).fillna("unknown").astype(str)
        for v in lab.value_counts().index[:6]:
            idx = lab.index[lab == v]
            groups.append({"label": f"{col}: {v}", "x": r(pb.loc[idx, "x"], 3), "y": r(pb.loc[idx, "y"], 3),
                           "z": r(pb.loc[idx, "z"], 3), "ids": list(idx)})
    else:
        groups.append({"label": "samples", "x": r(pb["x"], 3), "y": r(pb["y"], 3), "z": r(pb["z"], 3),
                       "ids": list(pb.index)})

    def pts(frame, ids):
        return {"x": r(frame.loc[ids, "x"], 3), "y": r(frame.loc[ids, "y"], 3),
                "z": r(frame.loc[ids, "z"], 3), "ids": list(ids)}

    highlights = []
    for kind, ids in evidence_by_file(doc).get(f, {}).items():
        in_b = [s for s in ids if s in pb.index]
        moved = [s for s in in_b if s in pa.index]
        gone = [s for s in ids if s in pa.index and s not in pb.index]
        highlights.append({"kind": kind, "label": KIND_LABEL.get(kind, kind), "color": KIND_COLOR.get(kind, "#ffffff"),
                           "new": pts(pb, in_b), "old": pts(pa, moved), "gone": pts(pa, gone)})
    return {"file": f, "color_column": col, "explained": explained, "groups": groups, "highlights": highlights}


def build_analysis_views(doc, a_dir, b_dir):
    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    if not clin_name:
        return []
    ev = evidence_by_file(doc)
    views = []
    for a in doc.get("analyses", []):
        uc = a.get("uses_columns", {})
        gcols = uc.get(clin_name, [])
        others = [(fn, feats) for fn, feats in uc.items() if fn != clin_name and feats]
        if len(gcols) != 1 or not others:
            continue
        fn, feats = others[0]
        feat, gcol = feats[0], gcols[0]
        clin_a, clin_b = table(a_dir, clin_name), table(b_dir, clin_name)
        ta, tb = table(a_dir, fn), table(b_dir, fn)
        if any(x is None for x in (clin_a, clin_b, ta, tb)) or feat not in ta.columns:
            continue
        hits = set().union(*ev.get(fn, {}).values()) if ev.get(fn) else set()
        series = []
        for rel_name, clin, t in (("old release", clin_a, ta), ("new release", clin_b, tb)):
            if feat not in t.columns:
                continue
            ids = [s for s in t.index if str(s)[13:15] == "01"]
            g = clin[gcol].reindex(ids)
            v = pd.to_numeric(t.loc[ids, feat], errors="coerce")
            counts = g.value_counts()
            for grp in counts.index[counts >= 10][:3]:
                sel = [s for s in g.index[g == grp] if pd.notna(v.get(s))]
                hs = [s for s in sel if s in hits]
                series.append({"name": f"{rel_name} | {grp}", "release": rel_name, "group": str(grp),
                               "values": r(v.loc[sel]), "mean": round(float(v.loc[sel].mean()), 3),
                               "n": len(sel), "hit_ids": hs, "hit_values": r(v.loc[hs])})
        if series:
            views.append({"id": a["id"], "feature": feat, "group_col": gcol, "file": fn, "series": series,
                          "before": a.get("before", {}), "after": a.get("after", {}),
                          "after_repair": a.get("after_repair", {})})
    return views


def build_swap(doc, a_dir, b_dir):
    swaps = [f for f in doc.get("findings", []) if f.get("kind") == "sample_swap"]
    if not swaps:
        return None
    samples = [e for e in swaps[0]["evidence"] if e.get("kind") == "sample"]
    if not samples:
        return None
    file, s1 = samples[0]["file"], samples[0]["id"]
    note = samples[0].get("note", "")
    if "profile of " in note:
        s2 = note.split("profile of ")[-1].split(",")[0].strip()
    elif len(samples) > 1:
        s2 = samples[1]["id"]
    else:
        return None
    ta, tb = table(a_dir, file), table(b_dir, file)
    if ta is None or tb is None or s1 not in tb.index or s1 not in ta.index or s2 not in ta.index:
        return None
    cols = [c for c in ta.columns if c in tb.columns]
    b1 = pd.to_numeric(tb.loc[s1, cols], errors="coerce")
    a2 = pd.to_numeric(ta.loc[s2, cols], errors="coerce")
    a1 = pd.to_numeric(ta.loc[s1, cols], errors="coerce")
    ok = b1.notna() & a2.notna() & a1.notna()
    b1, a2, a1 = b1[ok], a2[ok], a1[ok]
    return {"file": file, "s1": s1, "s2": s2, "features": list(b1.index),
            "new_s1": b1.round(6).tolist(), "old_s2": a2.round(6).tolist(), "old_s1": a1.round(6).tolist(),
            "maxdiff_match": float((b1 - a2).abs().max()), "maxdiff_own": float((b1 - a1).abs().max()),
            "n": int(len(b1))}


def safe(fn, *args):
    try:
        return fn(*args)
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--findings", required=True)
    p.add_argument("--b", required=True)
    p.add_argument("--a", default=str(RAW))
    p.add_argument("--out", default=None)
    args = p.parse_args()

    doc = json.loads(Path(args.findings).read_text(encoding="utf-8"))
    data = {
        "meta": {"a": str(args.a), "b": str(args.b), "findings": str(args.findings),
                 "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")},
        "doc": doc,
        "pipeline": safe(build_pipeline, doc, args.a),
        "flow": safe(build_flow, doc, args.a, args.b),
        "pca": safe(build_pca, doc, args.a, args.b),
        "analysis_views": safe(build_analysis_views, doc, args.a, args.b),
        "swap": safe(build_swap, doc, args.a, args.b),
    }
    payload = json.dumps(data, default=str).replace("</", "<\\/")
    html = TEMPLATE.read_text(encoding="utf-8").replace("__SNS_DATA__", payload)
    out = Path(args.out) if args.out else Path(args.findings).with_suffix(".html")
    out.write_text(html, encoding="utf-8")
    print(json.dumps({"dashboard": str(out), "bytes": out.stat().st_size}))


if __name__ == "__main__":
    main()