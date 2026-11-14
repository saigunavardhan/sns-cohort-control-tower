import importlib.util
import inspect
import re
import shutil
import tempfile
from pathlib import Path

import pandas as pd

TOKEN_RE = re.compile(r"[A-Za-z0-9_.\-]+")
RANK = {"Ready": 0, "Conditional": 1, "Insufficient Evidence": 2, "Blocked": 3}
HIDDEN_KEYS = ("analysis", "inputs", "sample_ids")


def load_analysis_module(path):
    path = Path(path)
    spec = importlib.util.spec_from_file_location("sns_reference_analyses", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    funcs = {
        name.upper(): fn
        for name, fn in vars(module).items()
        if re.fullmatch(r"a\d+", name) and callable(fn)
    }
    ordered = dict(sorted(funcs.items(), key=lambda kv: int(kv[0][1:])))
    return module, ordered


def describe(fn):
    lines = (fn.__doc__ or "").strip().splitlines()
    return lines[0].strip() if lines else fn.__name__


def used_tokens(fn):
    try:
        return set(TOKEN_RE.findall(inspect.getsource(fn)))
    except (OSError, TypeError):
        return set()


def run_analyses(module, funcs, release_dir):
    try:
        loaded = module.load(str(release_dir))
    except Exception as exc:
        err = f"could not load release: {type(exc).__name__}: {exc}"
        return {aid: {"analysis": aid, "error": err} for aid in funcs}
    if not isinstance(loaded, tuple):
        loaded = (loaded,)
    results = {}
    for aid, fn in funcs.items():
        n = len(inspect.signature(fn).parameters)
        try:
            results[aid] = fn(*loaded[:n])
        except Exception as exc:
            results[aid] = {"analysis": aid, "error": f"{type(exc).__name__}: {exc}"}
    return results


def public(result):
    return {k: v for k, v in (result or {}).items() if k not in HIDDEN_KEYS}


def _restore_rows(table, source, ids):
    ids = [s for s in ids if s in source.index]
    if not ids:
        return table
    rows = source.loc[ids]
    if list(rows.columns) != list(table.columns):
        rows = rows.loc[:, ~rows.columns.duplicated()].reindex(
            columns=table.columns[~table.columns.duplicated()]
        )
        table = table.loc[:, ~table.columns.duplicated()]
    table = table.drop(index=[s for s in ids if s in table.index])
    return pd.concat([table, rows])


def repaired_table(layer_a, layer_b, rep):
    a_tab = layer_a.table[~layer_a.table.index.duplicated()]
    table = layer_b.table[~layer_b.table.index.duplicated()].copy()
    drop = {c for d in rep.get("duplicates", []) for c in d["new_copies"]}
    drop |= {r["to"] for r in rep.get("renamed", [])}
    restore = {s for cycle in rep.get("swaps", []) for s in cycle}
    restore |= {r["sample"] for r in rep.get("relabeled", [])}
    restore |= {v["sample"] for v in rep.get("value_changes", [])}
    restore |= set(rep.get("dropped", []))
    restore |= {r["from"] for r in rep.get("renamed", [])}
    table = table.drop(index=[s for s in drop if s in table.index])
    return _restore_rows(table, a_tab, sorted(restore))


def write_layer(layer, table, path):
    if layer.layout == "wide":
        out = table.T
        out.index.name = layer.feature_field
    else:
        out = table.copy()
        out.index.name = layer.id_field
    out.to_csv(path, sep="," if ".csv" in path.name.lower() else "\t")


def build_repaired_release(rel_a, rel_b, diff_report):
    out_dir = Path(tempfile.mkdtemp(prefix="sns_repaired_"))
    reports = {r["file"]: r for r in diff_report["layers"]}
    for name, layer_b in rel_b.items():
        rep = reports.get(name)
        target = out_dir / name
        if rep is None or rep["content_same"] or name not in rel_a:
            shutil.copy2(layer_b.path, target)
        else:
            write_layer(layer_b, repaired_table(rel_a[name], layer_b, rep), target)
    for name, layer_a in rel_a.items():
        if name not in rel_b:
            shutil.copy2(layer_a.path, out_dir / name)
    return out_dir


def find_clinical(release):
    longs = [layer for layer in release.values() if layer.layout == "long"]
    named = [layer for layer in longs if "clinical" in layer.name.lower()]
    pool = named or longs
    return max(pool, key=lambda layer: layer.table.shape[1]) if pool else None


def concentration(clin_tab, samples, preferred=(), max_levels=12, limit=3, min_share=0.8):
    samples = [s for s in samples if s in clin_tab.index]
    if len(samples) < 2:
        return []
    hits = []
    for col in clin_tab.columns:
        series = clin_tab[col]
        if pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
            continue
        n_levels = series.nunique(dropna=True)
        if n_levels < 2 or n_levels > max_levels:
            continue
        vals = series.loc[samples].dropna()
        if len(vals) < max(2, min_share * len(samples)):
            continue
        top = vals.mode().iloc[0]
        n_top = int((vals == top).sum())
        share = n_top / len(samples)
        if share < min_share:
            continue
        base = float((series.dropna() == top).mean())
        if base >= 0.5 or base == 0:
            continue
        hits.append({
            "column": col,
            "value": str(top),
            "n_samples": len(samples),
            "n_with_value": n_top,
            "share_of_lost": round(share, 3),
            "baseline_share": round(base, 3),
            "enrichment": round(share / base, 2),
            "used_by_an_analysis": col in preferred,
        })
    hits.sort(key=lambda h: (not h["used_by_an_analysis"], -h["enrichment"]))
    return hits[:limit]


def _present(layer):
    t = layer.table
    return set(t.index[~t.isna().all(axis=1)])


def group_coverage(clin_tab, column, value, layer_a, layer_b):
    members = set(clin_tab.index[clin_tab[column] == value])
    return {
        "in_a": len(members & _present(layer_a)),
        "in_b": len(members & _present(layer_b)),
    }


def numeric_fields(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            out[key] = float(v)
        elif isinstance(v, dict):
            out.update(numeric_fields(v, key + "."))
    return out


def verdict_for(affected_by, before, after, min_group=10):
    if "error" in (after or {}):
        return "Insufficient Evidence", f"analysis failed on the new release: {after['error']}"
    if not affected_by:
        return "Ready", "no finding touches its input files, columns or samples"
    nb, na = numeric_fields(public(before)), numeric_fields(public(after))
    small = {
        k: v for k, v in na.items()
        if k.split(".")[-1].startswith("n_") and v < min_group
    }
    if small:
        listed = ", ".join(f"{k}={int(v)}" for k, v in small.items())
        return "Insufficient Evidence", f"too few samples after the change: {listed} (minimum {min_group})"
    shifts = []
    for k in sorted(nb.keys() & na.keys()):
        if k.split(".")[-1].startswith("n_"):
            continue
        b, a = nb[k], na[k]
        if b == a:
            continue
        if b != 0 and a != 0 and (b > 0) != (a > 0):
            return "Blocked", f"{k} changed sign: {b:g} in the old release, {a:g} in the new release"
        rel = abs(a - b) / abs(b) if b else float("inf")
        shifts.append((rel, k, b, a))
    touched = ", ".join(affected_by)
    if not shifts:
        return "Conditional", f"touched by {touched}, but its results did not move"
    rel, k, b, a = max(shifts)
    text = (
        f"largest shift is {k}: {b:g} in the old release, {a:g} in the new release, "
        f"a change of {abs(a - b):.4g}, which is {rel:.1%} of {abs(b):g}"
    )
    if rel >= 0.10:
        return "Blocked", f"touched by {touched}; {text} (10% or more)"
    return "Conditional", f"touched by {touched}; {text} (under 10%, but the data errors must be fixed first)"


def assess(analysis_path, rel_a, rel_b, folder_a, folder_b, diff_report):
    module, funcs = load_analysis_module(analysis_path)
    context = {}
    for aid, fn in funcs.items():
        tokens = used_tokens(fn)
        used = {
            name: sorted(tokens & set(layer.table.columns))
            for name, layer in rel_a.items()
        }
        context[aid] = {
            "description": describe(fn),
            "used_features": {k: v for k, v in used.items() if v},
        }
    before = run_analyses(module, funcs, folder_a)
    after = run_analyses(module, funcs, folder_b)
    repaired_dir = build_repaired_release(rel_a, rel_b, diff_report)
    try:
        repaired = run_analyses(module, funcs, repaired_dir)
    finally:
        shutil.rmtree(repaired_dir, ignore_errors=True)
    for aid, ctx in context.items():
        res_a, res_b = before.get(aid, {}), after.get(aid, {})
        ctx["inputs"] = res_a.get("inputs") or res_b.get("inputs") or sorted(ctx["used_features"])
        ctx["samples"] = set(res_a.get("sample_ids", [])) | set(res_b.get("sample_ids", []))
    return {
        "analysis_path": str(analysis_path),
        "context": context,
        "before": before,
        "after": after,
        "repaired": repaired,
    }