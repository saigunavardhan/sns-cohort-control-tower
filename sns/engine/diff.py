import warnings

import numpy as np
import pandas as pd

NA_TOKEN = "<NA>"


def is_numeric_table(table):
    return all(pd.api.types.is_numeric_dtype(t) for t in table.dtypes)


def _norm_value(x):
    if x is None or (not isinstance(x, str) and pd.isna(x)):
        return NA_TOKEN
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    if isinstance(x, (int, float, np.integer, np.floating)):
        return repr(float(x))
    return str(x)


def _col_str(s):
    if pd.api.types.is_bool_dtype(s):
        s = s.astype(object)
    elif pd.api.types.is_numeric_dtype(s):
        s = s.astype("float64")
    return s.map(_norm_value)


def _as_str(obj):
    if isinstance(obj, pd.Series):
        return _col_str(obj)
    parts = [_col_str(obj.iloc[:, i]) for i in range(obj.shape[1])]
    out = pd.concat(parts, axis=1)
    out.columns = obj.columns
    return out


def row_hashes(table):
    prepared = table.astype("float64") if is_numeric_table(table) else _as_str(table)
    return pd.util.hash_pandas_object(prepared, index=False)


def align(ta, tb):
    if list(ta.columns) == list(tb.columns):
        return ta, tb, [], []
    ca = list(dict.fromkeys(ta.columns))
    cb = list(dict.fromkeys(tb.columns))
    sa, sb = set(ca), set(cb)
    common = [c for c in ca if c in sb]
    xa = ta.loc[:, ~ta.columns.duplicated()][common]
    xb = tb.loc[:, ~tb.columns.duplicated()][common]
    return xa, xb, [c for c in cb if c not in sa], [c for c in ca if c not in sb]


def hash_groups(hashes, empty):
    groups = {}
    for sid, h in hashes.items():
        if not empty[sid]:
            groups.setdefault(h, []).append(sid)
    return groups


def find_cycles(edges):
    cycles, seen = [], set()
    for start in edges:
        if start in seen:
            continue
        path, node = [], start
        while node in edges and node not in path and node not in seen:
            path.append(node)
            node = edges[node]
        if node in path:
            cycle = path[path.index(node):]
            if len(cycle) >= 2:
                cycles.append(cycle)
        seen.update(path)
    return cycles


def _val(x):
    if pd.isna(x):
        return None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return float(x)
    return str(x)


def cell_diff(ra, rb, numeric, max_cols=15):
    na_a = ra.isna().to_numpy()
    na_b = rb.isna().to_numpy()
    if numeric:
        va = ra.to_numpy(dtype=float)
        vb = rb.to_numpy(dtype=float)
        mask = ~((va == vb) | (na_a & na_b))
    else:
        mask = (_as_str(ra).to_numpy() != _as_str(rb).to_numpy())
    out = {
        "n_changed": int(mask.sum()),
        "n_features": int(len(ra)),
        "newly_missing": int((~na_a & na_b).sum()),
        "newly_present": int((na_a & ~na_b).sum()),
        "all_missing_in_b": bool(na_b.all() and not na_a.all()),
    }
    if numeric:
        d = np.abs(va - vb)
        out["max_abs_diff"] = float(np.nanmax(d)) if np.isfinite(d).any() else None
    idx = np.flatnonzero(mask)[:max_cols]
    out["changed_columns"] = [
        {"column": str(ra.index[i]), "a": _val(ra.iloc[i]), "b": _val(rb.iloc[i])}
        for i in idx
    ]
    return out


def _unit_rows(values, mu):
    v = np.nan_to_num(values.astype(np.float32) - mu)
    n = np.linalg.norm(v, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return v / n


def nearest_profiles(xa, xb, queries):
    a_vals = xa.to_numpy(dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mu = np.nan_to_num(np.nanmean(a_vals, axis=0))
    a_unit = _unit_rows(a_vals, mu)
    q_unit = _unit_rows(xb.loc[queries].to_numpy(dtype=np.float32), mu)
    sims = q_unit @ a_unit.T
    ids = list(xa.index)
    pos = {s: i for i, s in enumerate(ids)}
    result = {}
    for i, s in enumerate(queries):
        row = sims[i]
        order = np.argsort(-row)
        best = next((j for j in order if ids[j] != s), None)
        result[s] = {
            "self_corr": round(float(row[pos[s]]), 4) if s in pos else None,
            "best_other": ids[best] if best is not None else None,
            "best_other_corr": round(float(row[best]), 4) if best is not None else None,
        }
    return result


def compare_layers(a, b, near_threshold=0.95, margin=0.05):
    report = {
        "file": a.name,
        "layout": a.layout,
        "content_same": a.sha256 == b.sha256,
        "n_samples_a": int(a.table.shape[0]),
        "n_samples_b": int(b.table.shape[0]),
    }
    if report["content_same"]:
        return report

    xa, xb, f_add, f_rem = align(a.table, b.table)
    report["n_features_added"] = len(f_add)
    report["n_features_removed"] = len(f_rem)
    report["features_added"] = f_add[:50]
    report["features_removed"] = f_rem[:50]
    report["duplicate_ids_a"] = sorted(set(xa.index[xa.index.duplicated()]))
    report["duplicate_ids_b"] = sorted(set(xb.index[xb.index.duplicated()]))
    xa = xa[~xa.index.duplicated()]
    xb = xb[~xb.index.duplicated()]

    numeric = is_numeric_table(xa) and is_numeric_table(xb)
    ha, hb = row_hashes(xa), row_hashes(xb)
    empty_a = xa.isna().all(axis=1)
    empty_b = xb.isna().all(axis=1)
    ga, gb = hash_groups(ha, empty_a), hash_groups(hb, empty_b)

    ids_a, ids_b = set(xa.index), set(xb.index)
    removed = sorted(ids_a - ids_b)
    added = sorted(ids_b - ids_a)
    common = sorted(ids_a & ids_b)
    modified = [s for s in common if ha[s] != hb[s]]

    origin = {}
    for s in modified + added:
        src = [t for t in ga.get(hb[s], []) if t != s]
        if src:
            origin[s] = src

    edges = {s: src[0] for s, src in origin.items() if len(src) == 1 and s in ids_a}
    cycles = find_cycles(edges)
    in_cycle = {s for c in cycles for s in c}

    duplicates = []
    for h, ids in gb.items():
        if len(ids) > 1:
            a_ids = ga.get(h, [])
            if set(ids) != set(a_ids):
                duplicates.append({
                    "data_owner_in_a": sorted(a_ids),
                    "ids_sharing_profile_in_b": sorted(ids),
                    "new_copies": sorted(set(ids) - set(a_ids)),
                })
    in_dup = {s for d in duplicates for s in d["new_copies"]}

    removed_set = set(removed)
    renamed = []
    for x in added:
        gone = [r for r in origin.get(x, []) if r in removed_set]
        if gone and x not in in_dup:
            renamed.append({"from": gone[0], "to": x})
    renamed_from = {r["from"] for r in renamed}

    relabeled = [
        {"sample": s, "holds_data_of": origin[s]}
        for s in modified
        if s in origin and s not in in_cycle and s not in in_dup
    ]

    unexplained = [s for s in modified if s not in origin]
    new_ids = [x for x in added if x not in origin]
    dropped = [r for r in removed if r not in renamed_from]

    near = {}
    if numeric and (unexplained or new_ids):
        near = nearest_profiles(xa, xb, unexplained + new_ids)

    value_changes = []
    for s in unexplained:
        item = {"sample": s}
        item.update(cell_diff(xa.loc[s], xb.loc[s], numeric))
        if s in near:
            item.update(near[s])
            item["likely_mislabel"] = bool(
                near[s]["best_other_corr"] is not None
                and near[s]["best_other_corr"] >= near_threshold
                and near[s]["best_other_corr"] - (near[s]["self_corr"] or 0) >= margin
            )
        value_changes.append(item)

    new_info = []
    for x in new_ids:
        item = {"sample": x}
        if x in near:
            item.update(near[x])
            item["near_duplicate_of"] = (
                near[x]["best_other"]
                if near[x]["best_other_corr"] is not None
                and near[x]["best_other_corr"] >= near_threshold
                else None
            )
        new_info.append(item)

    report.update({
        "numeric": numeric,
        "n_common": len(common),
        "n_modified": len(modified),
        "removed": removed,
        "added": added,
        "swaps": cycles,
        "duplicates": duplicates,
        "renamed": renamed,
        "relabeled": relabeled,
        "value_changes": value_changes,
        "dropped": dropped,
        "new": new_info,
    })
    return report


def _present(release):
    present = {}
    for name, layer in release.items():
        t = layer.table
        for s in t.index[~t.isna().all(axis=1)]:
            present.setdefault(s, set()).add(name)
    return present


def coverage_changes(rel_a, rel_b, limit=200):
    # Only compare files present in BOTH releases to avoid phantom coverage losses
    rel_a_shared = {k: v for k, v in rel_a.items() if k in rel_b}
    pa, pb = _present(rel_a_shared), _present(rel_b)
    changes = []
    for s in sorted(set(pa) | set(pb)):
        lost = sorted(pa.get(s, set()) - pb.get(s, set()))
        gained = sorted(pb.get(s, set()) - pa.get(s, set()))
        if lost or gained:
            changes.append({"sample": s, "lost_files": lost, "gained_files": gained})
    return {"n_samples_with_coverage_change": len(changes), "changes": changes[:limit]}


def compare_releases(rel_a, rel_b):
    files_a, files_b = set(rel_a), set(rel_b)
    return {
        "files_only_in_a": sorted(files_a - files_b),
        "files_only_in_b": sorted(files_b - files_a),
        "layers": [compare_layers(rel_a[n], rel_b[n]) for n in sorted(files_a & files_b)],
        "coverage": coverage_changes(rel_a, rel_b),
    }
