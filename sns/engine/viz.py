import hashlib

import numpy as np
import pandas as pd

from .release import load_layer

KIND_LABEL = {
    "sample_swap": "swapped",
    "duplicate": "copies",
    "membership_change": "lost or renamed",
    "missing_data": "values missing",
    "label_conflict": "labels changed",
    "other": "values changed",
}


class Tables:
    def __init__(self):
        self._cache = {}

    def get(self, folder, name):
        if not folder or not name:
            return None
        from pathlib import Path
        p = Path(folder) / name
        key = str(p)
        if key not in self._cache:
            t = None
            if p.exists():
                layer = load_layer(p)
                if layer is not None:
                    t = layer.table[~layer.table.index.duplicated()]
            self._cache[key] = t
        return self._cache[key]


def _num(t):
    if all(pd.api.types.is_numeric_dtype(d) for d in t.dtypes):
        return t.astype(float)
    return t.apply(pd.to_numeric, errors="coerce")


def _r(v, n=4):
    out = []
    for x in v:
        try:
            x = float(x)
        except (TypeError, ValueError):
            out.append(None)
            continue
        out.append(round(x, n) if np.isfinite(x) else None)
    return out


def _jitter(sid):
    return (int(hashlib.md5(sid.encode()).hexdigest()[:6], 16) / 0xFFFFFF - 0.5) * 0.6


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


def pipeline(doc, tables, a_dir):
    fa, fb = release_files(doc)
    changed = [n for n in fb if n in fa and fa[n] != fb[n]]
    rows = sum(len(t) for t in (tables.get(a_dir, n) for n in changed) if t is not None)
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
    proof = "n/a" if not affected else ("yes" if len(restored) == len(affected) else "partial")
    return {"files_total": len(fb), "files_changed": changed, "rows_checked": rows,
            "samples_changed": len(samples), "n_findings": len(findings),
            "kinds": sorted({KIND_LABEL.get(f.get("kind"), f.get("kind", "")) for f in findings}),
            "group": group, "affected": affected, "safe": safe, "restored": restored, "proof": proof}


def ribbon(doc, tables, a_dir, b_dir):
    f = main_file(doc)
    if not f:
        return None
    ta, tb = tables.get(a_dir, f), tables.get(b_dir, f)
    if ta is None or tb is None:
        return None
    ids_a, ids_b = set(ta.index), set(tb.index)
    allids = sorted(ids_a | ids_b)
    ev = evidence_by_file(doc).get(f, {})
    status = {}
    for kind, ids in ev.items():
        for s in ids:
            status.setdefault(s, kind)
    marks = []
    for i, s in enumerate(allids):
        if s in status:
            where = "both" if (s in ids_a and s in ids_b) else ("only_a" if s in ids_a else "only_b")
            marks.append({"i": i, "id": s, "kind": status[s], "where": where})
    return {"file": f, "n": len(allids), "n_a": len(ids_a), "n_b": len(ids_b), "marks": marks}


def flow(doc, tables, a_dir, b_dir):
    f = main_file(doc)
    if not f:
        return None
    ta, tb = tables.get(a_dir, f), tables.get(b_dir, f)
    if ta is None or tb is None:
        return None
    ids_a, ids_b = set(ta.index), set(tb.index)
    ev = evidence_by_file(doc).get(f, {})
    lost, new, both = ids_a - ids_b, ids_b - ids_a, ids_a & ids_b
    groups = {k: v & both for k, v in ev.items() if v & both}
    flagged = set().union(*groups.values()) if groups else set()
    nodes = [{"label": f"Release A, {len(ids_a)}", "kind": "release"},
             {"label": f"Release B, {len(ids_b)}", "kind": "release"},
             {"label": f"unchanged, {len(both - flagged)}", "kind": "unchanged"}]
    links = [{"s": 0, "t": 2, "v": len(both - flagged), "kind": "unchanged"},
             {"s": 2, "t": 1, "v": len(both - flagged), "kind": "unchanged"}]
    for k, v in groups.items():
        nodes.append({"label": f"{KIND_LABEL.get(k, k)}, {len(v)}", "kind": k})
        i = len(nodes) - 1
        links += [{"s": 0, "t": i, "v": len(v), "kind": k}, {"s": i, "t": 1, "v": len(v), "kind": k}]
    if lost:
        nodes.append({"label": f"lost, {len(lost)}", "kind": "membership_change"})
        links.append({"s": 0, "t": len(nodes) - 1, "v": len(lost), "kind": "membership_change"})
    if new:
        nodes.append({"label": f"new IDs, {len(new)}", "kind": "duplicate"})
        links.append({"s": len(nodes) - 1, "t": 1, "v": len(new), "kind": "duplicate"})
    return {"file": f, "n_a": len(ids_a), "n_b": len(ids_b), "lost": len(lost), "new": len(new),
            "changed_same_id": len(flagged), "nodes": nodes, "links": links}


def _color_column(doc, f, clin):
    for a in doc.get("analyses", []):
        if f in a.get("inputs", []):
            cols = a.get("uses_columns", {}).get(clin, [])
            if cols:
                return cols[0]
    return None


def pca(doc, tables, a_dir, b_dir, max_features=2000):
    f = main_file(doc)
    if not f:
        return None
    ta, tb = tables.get(a_dir, f), tables.get(b_dir, f)
    if ta is None or tb is None:
        return None
    xa = _num(ta)
    xa = xa.loc[:, xa.notna().mean() > 0.5]
    if xa.shape[1] < 3:
        return None
    if xa.shape[1] > max_features:
        xa = xa[xa.var().sort_values(ascending=False).index[:max_features]]
    mu, sd = xa.mean(), xa.std().replace(0, 1).fillna(1)
    za = ((xa - mu) / sd).fillna(0).to_numpy()
    _, s, vt = np.linalg.svd(za, full_matrices=False)
    comps = vt[:3].T
    explained = (s[:3] ** 2 / (s ** 2).sum() * 100).round(1).tolist()
    xb = _num(tb.reindex(columns=xa.columns))
    zb = ((xb - mu) / sd).fillna(0).to_numpy()
    pa = pd.DataFrame(za @ comps, index=xa.index, columns=list("xyz"))
    pb = pd.DataFrame(zb @ comps, index=xb.index, columns=list("xyz"))
    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    clin = tables.get(a_dir, clin_name)
    col = _color_column(doc, f, clin_name) if clin is not None else None

    def pts(frame, ids):
        ids = list(ids)
        return {"x": _r(frame.loc[ids, "x"], 3), "y": _r(frame.loc[ids, "y"], 3),
                "z": _r(frame.loc[ids, "z"], 3), "ids": ids}

    groups = []
    if col and col in clin.columns:
        lab = clin[col].reindex(pb.index).fillna("unknown").astype(str)
        for v in lab.value_counts().index[:5]:
            groups.append({"label": v, **pts(pb, lab.index[lab == v])})
    else:
        groups.append({"label": "samples", **pts(pb, pb.index)})
    highlights = []
    for kind, ids in evidence_by_file(doc).get(f, {}).items():
        in_b = [s for s in ids if s in pb.index]
        moved = [s for s in in_b if s in pa.index]
        gone = [s for s in ids if s in pa.index and s not in pb.index]
        highlights.append({"kind": kind, "label": KIND_LABEL.get(kind, kind),
                           "new": pts(pb, in_b), "old": pts(pa, moved), "gone": pts(pa, gone)})
    return {"file": f, "color_column": col, "explained": explained, "groups": groups, "highlights": highlights}


def _repaired_position(s, kind, note, in_a, pa, pb):
    if kind is None:
        return pb
    if kind == "duplicate":
        return None if note.startswith("identical to") and not in_a else pb
    if kind == "membership_change" and not in_a:
        return None if note.startswith("same data as") else pb
    if in_a:
        return pa
    return pb


def replays(doc, tables, a_dir, b_dir):
    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    if not clin_name:
        return []
    out = []
    for a in doc.get("analyses", []):
        uc = a.get("uses_columns", {})
        gcols = uc.get(clin_name, [])
        others = [(fn, feats) for fn, feats in uc.items() if fn != clin_name and feats]
        if len(gcols) != 1 or not others:
            continue
        fn, feats = others[0]
        feat, gcol = feats[0], gcols[0]
        ca, cb = tables.get(a_dir, clin_name), tables.get(b_dir, clin_name)
        ta, tb = tables.get(a_dir, fn), tables.get(b_dir, fn)
        if any(x is None for x in (ca, cb, ta, tb)) or feat not in ta.columns:
            continue
        prim_a = [s for s in ta.index if str(s)[13:15] == "01"]
        g_a = ca[gcol].reindex(prim_a)
        v_a = pd.to_numeric(ta.loc[prim_a, feat], errors="coerce")
        counts = g_a[v_a.notna()].dropna().value_counts()
        groups = [str(g) for g in counts.index[counts >= 10][:4]]
        if not groups:
            continue
        kind_of, note_of = {}, {}
        for f in doc.get("findings", []):
            for e in f.get("evidence", []):
                if e.get("kind") == "sample" and e.get("file") == fn:
                    kind_of.setdefault(e["id"], f["kind"])
                    note_of.setdefault(e["id"], e.get("note") or "")

        def side(clin, t, s):
            if feat not in t.columns or s not in t.index or str(s)[13:15] != "01" or s not in clin.index:
                return None
            g = clin.at[s, gcol]
            v = pd.to_numeric(t.at[s, feat], errors="coerce")
            if pd.isna(g) or str(g) not in groups or pd.isna(v):
                return None
            return [groups.index(str(g)), round(float(v), 4)]

        samples = []
        for s in sorted(set(ta.index) | set(tb.index)):
            pa, pb = side(ca, ta, s), side(cb, tb, s)
            kind, note = kind_of.get(s), note_of.get(s, "")
            pr = _repaired_position(s, kind, note, s in ta.index, pa, pb)
            if pa is None and pb is None and pr is None:
                continue
            samples.append({"id": s, "j": round(_jitter(s), 3), "a": pa, "b": pb, "r": pr,
                            "kind": kind, "note": note})
        stats = {}
        for key in ("a", "b", "r"):
            means, ns = [], []
            for gi in range(len(groups)):
                vals = [x[key][1] for x in samples if x[key] and x[key][0] == gi]
                means.append(round(float(np.mean(vals)), 3) if vals else None)
                ns.append(len(vals))
            stats[key] = {"means": means, "n": ns}
        out.append({"id": a["id"], "description": a.get("description", ""), "feature": feat, "group_col": gcol,
                    "groups": groups, "samples": samples, "stats": stats,
                    "before": a.get("before", {}), "after": a.get("after", {}),
                    "after_repair": a.get("after_repair", {}),
                    "repair_restores_before": a.get("repair_restores_before")})
    return out


def swap(doc, tables, a_dir, b_dir):
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
    ta, tb = tables.get(a_dir, file), tables.get(b_dir, file)
    if ta is None or tb is None or s1 not in tb.index or s1 not in ta.index or s2 not in ta.index:
        return None
    cols = [c for c in ta.columns if c in tb.columns]
    b1 = pd.to_numeric(tb.loc[s1, cols], errors="coerce")
    a2 = pd.to_numeric(ta.loc[s2, cols], errors="coerce")
    a1 = pd.to_numeric(ta.loc[s1, cols], errors="coerce")
    ok = b1.notna() & a2.notna() & a1.notna()
    b1, a2, a1 = b1[ok], a2[ok], a1[ok]
    return {"file": file, "s1": s1, "s2": s2, "features": list(b1.index), "n": int(len(b1)),
            "new_s1": b1.round(6).tolist(), "old_s2": a2.round(6).tolist(), "old_s1": a1.round(6).tolist(),
            "maxdiff_match": float((b1 - a2).abs().max()), "maxdiff_own": float((b1 - a1).abs().max())}


def _safe(fn, *args):
    try:
        return fn(*args)
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def build_viz(doc, a_dir, b_dir):
    tables = Tables()
    return {
        "pipeline": _safe(pipeline, doc, tables, a_dir),
        "ribbon": _safe(ribbon, doc, tables, a_dir, b_dir),
        "replays": _safe(replays, doc, tables, a_dir, b_dir),
        "swap": _safe(swap, doc, tables, a_dir, b_dir),
    }