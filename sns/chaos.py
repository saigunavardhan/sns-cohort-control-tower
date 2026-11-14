import argparse
import json
import random
import re
import shutil
import sys
import tempfile
import time
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from engine.diff import compare_releases, is_numeric_table  # noqa: E402
from engine.findings import build_document  # noqa: E402
from engine.impact import assess, find_clinical, write_layer  # noqa: E402
from engine.release import load_layer, load_release  # noqa: E402

DATA = ROOT / "data" / "tcga-brca"
ANALYSES = DATA / "reference_analyses.py"
SCHEMA = ROOT / "format" / "findings.schema.json"
COPY_LIMIT = 10 * 1024 * 1024
SUBGROUP_COL = "PAM50Call_RNAseq"
FILE_WEIGHTS = {"protein": 4, "clinical": 3, "survival": 1, "mutation": 1, "copy": 1, "rna": 1}
OPS_ANY = ["swap", "dup_replace", "drop", "subgroup_drop", "rename",
           "blank", "partial_blank", "edit", "dup_id"]
OPS_NUMERIC = ["near_dup", "mislabel", "drop_feature"]
DROPPING = {"dup_replace", "drop", "subgroup_drop"}


class Skip(Exception):
    pass


def expect(kinds, file, samples=None, features=None, subgroup=None):
    return {"kinds": list(kinds), "file": file, "samples": sorted(samples or []),
            "features": sorted(features or []), "subgroup": subgroup}


class Base:
    def __init__(self, rel_a):
        self.rel_a = rel_a
        self.clin_name = find_clinical(rel_a).name
        self.nonempty = {n: list(l.table.index[~l.table.isna().all(axis=1)]) for n, l in rel_a.items()}
        self.numeric = {n: is_numeric_table(l.table) for n, l in rel_a.items()}
        self.sparse = {}
        for n, l in rel_a.items():
            if self.numeric[n]:
                v = np.nan_to_num(l.table.to_numpy(dtype=float))
                self.sparse[n] = float(np.mean(v != 0)) < 0.05
            else:
                self.sparse[n] = False
        self.all_ids = sorted(set().union(*(set(l.table.index) for l in rel_a.values())))
        self.weights = {n: next((w for k, w in FILE_WEIGHTS.items() if k in n.lower()), 1) for n in rel_a}
        self.cat_cols = {}
        for n, l in rel_a.items():
            if self.numeric[n]:
                self.cat_cols[n] = []
                continue
            t = l.table
            self.cat_cols[n] = [
                c for c, d in t.dtypes.items()
                if not pd.api.types.is_numeric_dtype(d) and not pd.api.types.is_bool_dtype(d)
                and 2 <= t[c].nunique(dropna=True) <= 12
            ]
        self._std = {}
        analysis_text = ANALYSES.read_text(encoding="utf-8")
        self.analysis_inputs = set(re.findall(r'\b[\w-]+\.tsv(?:\.gz)?\b', analysis_text))

    def std(self, name):
        if name not in self._std:
            t = self.rel_a[name].table.astype("float64")
            s = t.std(axis=0, skipna=True).to_numpy()
            s = np.where(np.isfinite(s) & (s > 0), s, 1.0)
            self._std[name] = s
        return self._std[name]


class Chaos:
    def __init__(self, base, seed):
        self.base = base
        self.rng = random.Random(seed)
        self.nrng = np.random.default_rng(seed)
        self.tables = {}
        self.used = defaultdict(set)
        self.changes = []
        self.dropped_files = []
        self.dropping = set()

    def t(self, f):
        if f not in self.tables:
            t = self.base.rel_a[f].table.copy()
            self.tables[f] = t.astype("float64") if self.base.numeric[f] else t
        return self.tables[f]

    def put(self, f, table):
        self.tables[f] = table

    def floatify(self, f):
        t = self.t(f)
        if self.base.numeric[f]:
            return
        conv = {c: "float64" for c, d in t.dtypes.items()
                if pd.api.types.is_integer_dtype(d) and not pd.api.types.is_bool_dtype(d)}
        if conv:
            self.put(f, t.astype(conv))

    def pick(self, f, k, pool=None):
        pool = self.base.nonempty[f] if pool is None else pool
        current = set(self.t(f).index)
        cands = [s for s in pool if s in current and s not in self.used[f]]
        if len(cands) < k:
            raise Skip
        chosen = self.rng.sample(cands, k)
        self.used[f].update(chosen)
        return chosen

    def outside(self, f, k):
        current = set(self.t(f).index)
        cands = [s for s in self.base.all_ids if s not in current and s not in self.used[f]]
        chosen = self.rng.sample(cands, k) if len(cands) >= k else []
        while len(chosen) < k:
            sid = f"TCGA-ZZ-{self.rng.randint(0, 9999):04d}-01"
            if sid not in current and sid not in self.used[f] and sid not in chosen:
                chosen.append(sid)
        self.used[f].update(chosen)
        return chosen

    def relabel_rows(self, f, sources, targets, noise=0.0, drop=()):
        t = self.t(f)
        sub = t.loc[sources].copy()
        if noise:
            sub = sub + self.nrng.normal(0, 1, sub.shape) * noise * self.base.std(f)
        sub.index = list(targets)
        gone = [s for s in set(targets) | set(drop) if s in t.index]
        self.put(f, pd.concat([t.drop(index=gone), sub]))

    def op_swap(self, f):
        ids = self.pick(f, self.rng.choice([2, 2, 2, 3]))
        self.relabel_rows(f, ids[1:] + ids[:1], ids)
        return [expect(["sample_swap"], f, samples=ids)]

    def op_dup_replace(self, f):
        k = self.rng.randint(1, 5)
        removed, sources = self.pick(f, k), self.pick(f, k)
        new = self.outside(f, k)
        self.relabel_rows(f, sources, new, drop=removed)
        return [expect(["duplicate"], f, samples=new),
                expect(["membership_change"], f, samples=removed)]

    def op_drop(self, f):
        ids = self.pick(f, self.rng.randint(1, 5))
        self.put(f, self.t(f).drop(index=ids))
        return [expect(["membership_change"], f, samples=ids)]

    def op_subgroup_drop(self, f):
        clin = self.base.rel_a[self.base.clin_name].table
        if f == self.base.clin_name or SUBGROUP_COL not in clin.columns:
            raise Skip
        clin = clin[~clin.index.duplicated()]
        vals = clin[SUBGROUP_COL].reindex(self.base.nonempty[f])
        share = clin[SUBGROUP_COL].dropna().value_counts(normalize=True)
        options = [v for v, s in share.items() if s < 0.4 and int((vals == v).sum()) >= 8]
        if not options:
            raise Skip
        value = self.rng.choice(options)
        pool = list(vals.index[vals == value])
        ids = self.pick(f, self.rng.randint(3, 6), pool=pool)
        self.put(f, self.t(f).drop(index=ids))
        return [expect(["membership_change"], f, samples=ids, subgroup=str(value))]

    def op_rename(self, f):
        k = self.rng.randint(1, 3)
        old = self.pick(f, k)
        new = self.outside(f, k)
        self.relabel_rows(f, old, new, drop=old)
        return [expect(["membership_change"], f, samples=new)]

    def op_blank(self, f):
        ids = self.pick(f, self.rng.randint(1, 3))
        self.floatify(f)
        t = self.t(f)
        cols = [c for c, d in t.dtypes.items() if not pd.api.types.is_bool_dtype(d)]
        if self.base.numeric[f]:
            t.loc[ids, :] = np.nan
        else:
            t.loc[ids, cols] = np.nan
        return [expect(["missing_data"], f, samples=ids)]

    def op_partial_blank(self, f):
        ids = self.pick(f, self.rng.randint(1, 3))
        self.floatify(f)
        t = self.t(f)
        for s in ids:
            if self.base.numeric[f]:
                vals = t.loc[s].to_numpy(dtype=float).copy()
                present = np.flatnonzero(~np.isnan(vals))
                m = max(1, int(0.3 * len(present)))
                vals[self.nrng.choice(present, m, replace=False)] = np.nan
                t.loc[s] = vals
            else:
                cols = [c for c, d in t.dtypes.items()
                        if not pd.api.types.is_bool_dtype(d) and pd.notna(t.at[s, c])]
                m = max(1, int(0.3 * len(cols)))
                t.loc[s, self.rng.sample(cols, m)] = np.nan
        return [expect(["missing_data"], f, samples=ids)]

    def op_edit(self, f):
        ids = self.pick(f, self.rng.randint(1, 3))
        if self.base.numeric[f]:
            t = self.t(f)
            std = self.base.std(f)
            for s in ids:
                vals = t.loc[s].to_numpy(dtype=float).copy()
                present = np.flatnonzero(~np.isnan(vals))
                m = max(1, len(present) // 10)
                pos = self.nrng.choice(present, m, replace=False)
                vals[pos] = vals[pos] + self.nrng.normal(0, 1, m) * std[pos]
                t.loc[s] = vals
            return [expect(["other"], f, samples=ids)]
        self.floatify(f)
        t = self.t(f)
        num_cols = [c for c, d in t.dtypes.items()
                    if pd.api.types.is_numeric_dtype(d) and not pd.api.types.is_bool_dtype(d)]
        for s in ids:
            cats = [c for c in self.base.cat_cols[f] if pd.notna(t.at[s, c])]
            if cats:
                c = self.rng.choice(cats)
                levels = [v for v in t[c].dropna().unique() if v != t.at[s, c]]
                t.at[s, c] = self.rng.choice(levels)
            else:
                nums = [c for c in num_cols if pd.notna(t.at[s, c])]
                if not nums:
                    raise Skip
                c = self.rng.choice(nums)
                t.at[s, c] = float(t.at[s, c]) + 1.0
        return [expect(["label_conflict"], f, samples=ids)]

    def op_near_dup(self, f):
        if self.base.sparse[f]:
            raise Skip
        k = self.rng.randint(1, 3)
        sources = self.pick(f, k)
        new = self.outside(f, k)
        self.relabel_rows(f, sources, new, noise=0.05)
        return [expect(["duplicate"], f, samples=new)]

    def op_mislabel(self, f):
        k = self.rng.randint(1, 2)
        targets, sources = self.pick(f, k), self.pick(f, k)
        self.relabel_rows(f, sources, targets, noise=0.05)
        return [expect(["sample_swap"], f, samples=targets)]

    def op_drop_feature(self, f):
        t = self.t(f)
        cols = list(dict.fromkeys(t.columns))
        if len(cols) <= 20:
            raise Skip
        chosen = self.rng.sample(cols, self.rng.randint(1, 3))
        if "ERALPHA" in cols and "ERALPHA" not in chosen and self.rng.random() < 0.3:
            chosen[0] = "ERALPHA"
        self.put(f, t.drop(columns=chosen))
        return [expect(["missing_data"], f, features=chosen)]

    def op_dup_id(self, f):
        ids = self.pick(f, 1)
        t = self.t(f)
        self.put(f, pd.concat([t, t.loc[ids]]))
        return [expect(["duplicate"], f, samples=ids)]

    def op_drop_file(self, f):
        self.dropped_files.append(f)
        if f in self.base.analysis_inputs:
            return [expect(["missing_data"], f)]
        # No finding expected; file should appear in sns_run.files_not_in_release
        return [expect([], f)]

    def scenario(self, max_ops=4):
        n_ops = self.rng.randint(1, max_ops)
        files = list(self.base.rel_a)
        weights = [self.base.weights[f] for f in files]
        tries = 0
        while len(self.changes) < n_ops and tries < 60:
            tries += 1
            f = self.rng.choices(files, weights)[0]
            if f in self.dropped_files:
                continue
            if f not in self.tables and self.rng.random() < 0.04:
                op = "drop_file"
            else:
                op = self.rng.choice(OPS_ANY + (OPS_NUMERIC if self.base.numeric[f] else []))
            if op in DROPPING and f in self.dropping:
                continue
            try:
                exps = getattr(self, "op_" + op)(f)
            except Skip:
                continue
            if op in DROPPING:
                self.dropping.add(f)
            self.changes.append({"type": op, "file": f, "expect": exps})

    def build(self, out_dir):
        rel_b = {}
        for name, layer in self.base.rel_a.items():
            if name in self.dropped_files:
                continue
            path = out_dir / name
            if name in self.tables:
                write_layer(layer, self.tables[name], path)
                rel_b[name] = load_layer(path)
            else:
                if layer.path.stat().st_size <= COPY_LIMIT:
                    shutil.copy2(layer.path, path)
                rel_b[name] = layer
        return rel_b


def norm_id(i):
    return re.sub(r"\.\d+$", "", str(i))


def match(e, findings):
    need = set(e["samples"]) or set(e["features"]) or {e["file"]}
    for f in findings:
        if f["kind"] not in e["kinds"]:
            continue
        ids = {
            norm_id(r["id"]) for r in f["evidence"]
            if r.get("file") == e["file"] or (r["kind"] == "file" and r["id"] == e["file"])
        }
        if need <= ids:
            return f
    return None


def score(changes, doc):
    rows, matched = [], set()
    expected_ids = defaultdict(set)
    for ch in changes:
        for e in ch["expect"]:
            expected_ids[e["file"]] |= set(e["samples"])
            if not e["kinds"]:
                # Empty kinds means: expect no finding; file should be in files_not_in_release
                any_finding = any(
                    r.get("file") == e["file"]
                    or (r["kind"] == "file" and r["id"] == e["file"])
                    for fd in doc["findings"] for r in fd["evidence"]
                )
                in_run = e["file"] in doc.get("sns_run", {}).get("files_not_in_release", [])
                row = {"type": ch["type"], "file": ch["file"], "caught": not any_finding and in_run}
            else:
                f = match(e, doc["findings"])
                row = {"type": ch["type"], "file": ch["file"], "caught": f is not None}
                if f is not None:
                    matched.add(f["id"])
                if e["subgroup"]:
                    v = e["subgroup"]
                    row["subgroup_named"] = bool(f) and (
                        f"{v} samples" in f["summary"] or f"of them {v}" in f["summary"]
                    )
            rows.append(row)
    false_findings = [fd["id"] for fd in doc["findings"] if fd["id"] not in matched]
    extra_ids = sorted({
        norm_id(r["id"])
        for fd in doc["findings"] for r in fd["evidence"]
        if r["kind"] == "sample" and norm_id(r["id"]) not in expected_ids.get(r.get("file"), set())
    })
    return rows, false_findings, extra_ids


def validate(doc, validator):
    errors = [e.message for e in validator.iter_errors(doc)]
    fids = {f["id"] for f in doc["findings"]}
    errors += [f"repair {r['id']} points at a missing finding"
               for r in doc.get("repairs", []) if not set(r["addresses"]) <= fids]
    return errors


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=30)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--max-ops", type=int, default=4)
    p.add_argument("--out", default=str(HERE / "chaos_report.json"))
    args = p.parse_args()

    t_start = time.time()
    rel_a = load_release(DATA / "raw")
    base = Base(rel_a)
    validator = Draft202012Validator(json.loads(SCHEMA.read_text(encoding="utf-8")))

    per_type = defaultdict(lambda: [0, 0])
    sub_ok = sub_n = 0
    report, crashes, invalid, restore_fail = [], [], 0, []
    n_false = n_extra = 0

    for i in range(args.n):
        seed = args.seed * 1000 + i
        chaos = Chaos(base, seed)
        chaos.scenario(args.max_ops)
        tmp = Path(tempfile.mkdtemp(prefix="sns_chaos_"))
        t0 = time.time()
        try:
            rel_b = chaos.build(tmp)
            diff = compare_releases(rel_a, rel_b)
            impact = assess(ANALYSES, rel_a, rel_b, DATA / "raw", tmp, diff)
            doc = build_document(diff, rel_a, rel_b, "raw", f"chaos-{seed}", impact)
        except Exception as exc:
            crashes.append({"seed": seed, "changes": [f"{c['type']}:{c['file']}" for c in chaos.changes],
                            "error": traceback.format_exc(limit=4)})
            print(json.dumps({"seed": seed, "crash": f"{type(exc).__name__}: {exc}",
                              "changes": [f"{c['type']}:{c['file']}" for c in chaos.changes]}))
            continue
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        rows, false_findings, extra_ids = score(chaos.changes, doc)
        errors = validate(doc, validator)
        invalid += bool(errors)
        n_false += len(false_findings)
        n_extra += len(extra_ids)
        for r in rows:
            per_type[r["type"]][1] += 1
            per_type[r["type"]][0] += r["caught"]
            if "subgroup_named" in r:
                sub_n += 1
                sub_ok += r["subgroup_named"]
        not_restored = [
            a["id"] for a in doc["analyses"]
            if a["affected_by"] and "error" not in a["after"] and not a["repair_restores_before"]
        ]
        if not_restored:
            restore_fail.append({"seed": seed, "analyses": not_restored,
                                 "changes": [f"{c['type']}:{c['file']}" for c in chaos.changes]})

        line = {
            "seed": seed,
            "changes": [f"{c['type']}:{c['file']}" for c in chaos.changes],
            "caught": f"{sum(r['caught'] for r in rows)}/{len(rows)}",
            "missed": [f"{r['type']}:{r['file']}" for r in rows if not r["caught"]],
            "false_findings": false_findings,
            "extra_sample_ids": len(extra_ids),
            "verdict": doc["verdict"],
            "valid": not errors,
            "seconds": round(time.time() - t0, 1),
        }
        print(json.dumps(line))
        report.append({**line, "answer_key": chaos.changes, "schema_errors": errors[:5],
                       "extra_ids": extra_ids[:20], "findings": [
                           {k: f[k] for k in ("id", "kind", "severity", "summary")} for f in doc["findings"]]})

    caught = sum(v[0] for v in per_type.values())
    total = sum(v[1] for v in per_type.values())
    summary = {
        "releases": args.n,
        "crashes": len(crashes),
        "invalid_files": invalid,
        "changes_caught": f"{caught}/{total}",
        "catch_rate": round(caught / total, 4) if total else None,
        "by_type": {k: f"{v[0]}/{v[1]}" for k, v in sorted(per_type.items())},
        "false_findings": n_false,
        "extra_sample_ids": n_extra,
        "subgroup_named": f"{sub_ok}/{sub_n}",
        "repair_not_restoring": len(restore_fail),
        "seconds_total": round(time.time() - t_start, 1),
    }
    print(json.dumps({"summary": summary}, indent=2))
    Path(args.out).write_text(json.dumps({"summary": summary, "crashes": crashes,
                                          "repair_not_restoring": restore_fail, "releases": report},
                                         indent=2, default=str), encoding="utf-8")


if __name__ == "__main__":
    main()