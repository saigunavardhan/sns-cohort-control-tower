from datetime import datetime, timezone

from .impact import RANK, concentration, find_clinical, group_coverage, public, verdict_for

ENGINE_VERSION = "0.2"


def ref(kind, id_, release=None, file=None, note=None):
    r = {"kind": kind, "id": str(id_)}
    if release:
        r["release"] = release
    if file:
        r["file"] = file
    if note:
        r["note"] = note
    return r


def feature_kind(file):
    n = file.lower()
    if "protein" in n or "rppa" in n:
        return "protein"
    if any(k in n for k in ("rna", "mutation", "copy", "gene", "cnv", "expr")):
        return "gene"
    return "row"


def _plural(n, word):
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def _v(n, one, many):
    return one if n == 1 else many


def _join(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _feature_set(items):
    feats = set()
    for v in items:
        if v["n_changed"] > len(v["changed_columns"]):
            return None
        feats |= {c["column"] for c in v["changed_columns"]}
    return feats


class Builder:
    def __init__(self, diff, rel_a, rel_b, a_name, b_name, impact):
        self.diff, self.rel_a, self.rel_b = diff, rel_a, rel_b
        self.a, self.b = a_name, b_name
        self.ctx = impact["context"]
        clin = find_clinical(rel_a) or find_clinical(rel_b)
        self.clin_name = clin.name if clin else None
        self.clin = clin.table[~clin.table.index.duplicated()] if clin else None
        self.used_clinical = {
            col
            for c in self.ctx.values()
            for col in c["used_features"].get(self.clin_name, [])
        } if clin else set()
        self.cands = []
        self.subgroup_checks = []
        self.files_not_in_release = []
        self._all_analysis_inputs = set()
        for c in self.ctx.values():
            self._all_analysis_inputs |= set(c.get("inputs", []))

    def label(self, sample, file):
        if self.clin is None or sample not in self.clin.index:
            return sample
        cols = sorted({
            col
            for c in self.ctx.values()
            if file in c["inputs"]
            for col in c["used_features"].get(self.clin_name, [])
        })
        parts = [f"{col}={self.clin.at[sample, col]}" for col in cols if col in self.clin.columns]
        return f"{sample} ({', '.join(parts)})" if parts else sample

    def used_feature_refs(self, file):
        feats = sorted({
            f for c in self.ctx.values() if file in c["inputs"]
            for f in c["used_features"].get(file, [])
        })
        if file == self.clin_name:
            return []
        return [ref(feature_kind(file), f, file=file, note="used by an analysis") for f in feats]

    def add(self, **cand):
        cand.setdefault("samples", set())
        cand.setdefault("features", None)
        cand.setdefault("whole_file", False)
        self.cands.append(cand)

    def file_level(self):
        for f in self.diff["files_only_in_a"]:
            if f not in self._all_analysis_inputs:
                self.files_not_in_release.append(f)
                continue
            self.add(
                kind="missing_data", file=f, severity="high", confidence="high", whole_file=True,
                summary=f"File {f} exists in {self.a} but is missing from {self.b}.",
                evidence=[ref("file", f, self.a)],
                limitation="The tool cannot tell whether the file was dropped on purpose.",
                detected_by="file listing comparison",
                repair=("restore_from_previous_release", [ref("file", f, self.a)],
                        f"Restore {f} from {self.a} unless the data provider confirms it was removed on purpose."),
            )
        for f in self.diff["files_only_in_b"]:
            self.add(
                kind="other", file=f, severity="low", confidence="high",
                summary=f"New file {f} appears in {self.b} and was not in {self.a}.",
                evidence=[ref("file", f, self.b)],
                limitation="New files are not checked against the previous release because there is nothing to compare.",
                detected_by="file listing comparison",
                repair=("flag_for_review", [ref("file", f, self.b)],
                        "Confirm what the new file contains before any analysis uses it."),
            )

    def layer_level(self, rep):
        f = rep["file"]
        same_count = rep["n_samples_a"] == rep["n_samples_b"]
        feat_refs = self.used_feature_refs(f)

        for cycle in rep.get("swaps", []):
            n = len(cycle)
            what = "appear exchanged" if n == 2 else f"appear rotated among {n} samples"
            labelled = _join(self.label(s, f) for s in cycle)
            self.add(
                kind="sample_swap", file=f, severity="high", confidence="high", samples=set(cycle),
                summary=f"Profiles in {f} of {labelled} {what} between {self.a} and {self.b}.",
                evidence=[
                    ref("sample", s, self.b, f, note=f"now holds the {self.a} profile of {cycle[(i + 1) % n]}")
                    for i, s in enumerate(cycle)
                ] + feat_refs,
                limitation="The swap is inferred from exact profile matches trading places; the release notes do not say why.",
                detected_by="exact profile fingerprint match in both directions",
                repair=("restore_from_previous_release", [ref("sample", s, self.a, f) for s in cycle],
                        f"Restore these profiles as they were in {self.a} until the swap is explained."),
            )

        copies = []
        for d in rep.get("duplicates", []):
            owners = d["data_owner_in_a"] or [x for x in d["ids_sharing_profile_in_b"] if x not in d["new_copies"]]
            for c in d["new_copies"]:
                copies.append((c, owners))
        for sid in rep.get("duplicate_ids_b", []):
            if sid not in rep.get("duplicate_ids_a", []):
                copies.append((sid, ["the same ID listed twice"]))
        if copies:
            tail = ""
            if same_count and rep.get("dropped"):
                tail = f", while {_plural(len(rep['dropped']), 'real sample')} left the file, so the total count looks unchanged"
            self.add(
                kind="duplicate", file=f, severity="high", confidence="high",
                samples={c for c, _ in copies},
                summary=(f"{_plural(len(copies), 'sample ID')} in {f} {_v(len(copies), 'carries a profile', 'carry profiles')} "
                         f"identical to existing samples{tail}."),
                evidence=[
                    ref("sample", c, self.b, f, note=f"identical to {', '.join(o)}") for c, o in copies
                ] + feat_refs,
                limitation="Identical profiles strongly suggest copies, not new measurements.",
                detected_by="exact profile fingerprint shared by more IDs than in the previous release",
                repair=("exclude", [ref("sample", c, self.b, f) for c, _ in copies],
                        f"Exclude the duplicated profiles from analyses that use {f}."),
            )

        renamed = rep.get("renamed", [])
        if renamed:
            self.add(
                kind="membership_change", file=f, severity="medium", confidence="high",
                samples={r["from"] for r in renamed} | {r["to"] for r in renamed},
                summary=f"{_plural(len(renamed), 'sample')} in {f} changed ID while keeping identical data.",
                evidence=[ref("sample", r["to"], self.b, f, note=f"same data as {r['from']} in {self.a}") for r in renamed],
                limitation="An ID change with identical data may be a correction or an error; the data alone cannot tell which.",
                detected_by="exact profile fingerprint moved from a removed ID to a new ID",
                repair=("relabel", [ref("sample", r["to"], self.b, f, note=f"rename back to {r['from']}") for r in renamed],
                        f"Rename the IDs back to their {self.a} names unless the provider confirms the new names."),
            )

        relabeled = rep.get("relabeled", [])
        if relabeled:
            self.add(
                kind="sample_swap", file=f, severity="high", confidence="high",
                samples={r["sample"] for r in relabeled},
                summary=(f"{_plural(len(relabeled), 'sample')} in {f} now {_v(len(relabeled), 'holds a profile', 'hold profiles')} "
                         f"that belonged to other samples in {self.a}."),
                evidence=[
                    ref("sample", r["sample"], self.b, f, note=f"holds the {self.a} profile of {', '.join(r['holds_data_of'])}")
                    for r in relabeled
                ] + feat_refs,
                limitation="Inferred from exact profile matches; the original owner may also still hold its data.",
                detected_by="exact profile fingerprint match to a different sample",
                repair=("restore_from_previous_release", [ref("sample", r["sample"], self.a, f) for r in relabeled],
                        f"Restore these profiles from {self.a} until the change is explained."),
            )

        self.value_changes(rep, f, feat_refs)

        dropped = rep.get("dropped", [])
        if dropped:
            hits = concentration(self.clin, dropped, self.used_clinical) if self.clin is not None else []
            lead, extra = _plural(len(dropped), "sample"), ""
            whose = _v(len(dropped), "its", "their")
            if hits:
                h = hits[0]
                lead = f"{len(dropped)} {h['value']} samples ({h['column']})"
                if f in self.rel_a and f in self.rel_b:
                    cov = group_coverage(self.clin, h["column"], h["value"], self.rel_a[f], self.rel_b[f])
                    h.update({"file": f, "with_data_in_a": cov["in_a"], "with_data_in_b": cov["in_b"]})
                    extra = (f" {h['value']} samples with data in {f}: {cov['in_a']} in {self.a}, "
                             f"{cov['in_b']} in {self.b}.")
                self.subgroup_checks.append(h)
            tail = ", though the total sample count in this file is unchanged" if same_count else ""
            self.add(
                kind="membership_change", file=f, severity="medium", confidence="high", samples=set(dropped),
                summary=f"{lead} lost {whose} data in {f}{tail}.{extra}",
                evidence=[ref("sample", s, self.a, f, note=f"present in {self.a}, absent in {self.b}") for s in dropped],
                limitation=("Coverage for this subgroup drops; any per-subgroup result that uses this file is affected."
                            if hits else "The data does not say why these samples were removed."),
                detected_by="sample ID present in the old release and absent in the new one, with no matching profile elsewhere",
                repair=("flag_for_review", [ref("sample", s, self.a, f) for s in dropped],
                        f"Ask the data provider why these samples lost data in {f}; restore them from {self.a} if the loss was unintended."),
            )

        new = rep.get("new", [])
        near = [x for x in new if x.get("near_duplicate_of")]

        def fresh_note(x):
            if x.get("best_other") is None:
                return f"not in {self.a}"
            return f"not in {self.a}; closest existing sample {x['best_other']}, similarity {x['best_other_corr']:.3f}"

        fresh = [x for x in new if not x.get("near_duplicate_of")]
        if near:
            self.add(
                kind="duplicate", file=f, severity="medium", confidence="medium",
                samples={x["sample"] for x in near},
                summary=f"{_plural(len(near), 'new sample')} in {f} {_v(len(near), 'looks', 'look')} almost identical to existing samples.",
                evidence=[
                    ref("sample", x["sample"], self.b, f,
                        note=f"profile similarity {x['best_other_corr']:.3f} with {x['near_duplicate_of']}")
                    for x in near
                ],
                limitation="Based on profile similarity, not an exact match; a lab replicate can also look this similar.",
                detected_by="centered profile correlation against every sample in the old release",
                repair=("flag_for_review", [ref("sample", x["sample"], self.b, f) for x in near],
                        "Confirm these are new measurements and not edited copies before using them."),
            )
        if fresh:
            self.add(
                kind="membership_change", file=f, severity="low", confidence="high",
                samples={x["sample"] for x in fresh},
                summary=(f"{_plural(len(fresh), 'new sample')} appeared in {f} with "
                         f"{_v(len(fresh), 'a profile that matches', 'profiles that match')} no existing sample."),
                evidence=[ref("sample", x["sample"], self.b, f, note=fresh_note(x)) for x in fresh],
                limitation=("New samples cannot be checked against the previous release. In files where most values are 0 "
                            "(like mutation calls), a near copy of a sample with very few non-zero values cannot be told "
                            "apart from a real sample."),
                detected_by="sample ID absent in the old release, profile matches nothing",
                repair=("flag_for_review", [ref("sample", x["sample"], self.b, f) for x in fresh],
                        "Confirm the new samples are expected in this release."),
            )

        if rep.get("n_features_removed"):
            feats = rep["features_removed"]
            self.add(
                kind="missing_data", file=f, severity="medium", confidence="high",
                features=set(feats),
                summary=f"{_plural(rep['n_features_removed'], 'column')} of {f} disappeared in {self.b}.",
                evidence=[ref(feature_kind(f), x, self.a, f, note=f"absent in {self.b}") for x in feats],
                limitation="Any analysis that reads these columns will fail or change.",
                detected_by="feature list comparison",
                repair=("restore_from_previous_release", [ref(feature_kind(f), x, self.a, f) for x in feats],
                        f"Restore the missing columns from {self.a}."),
            )
        if rep.get("n_features_added"):
            feats = rep["features_added"]
            self.add(
                kind="other", file=f, severity="low", confidence="high",
                features=set(feats),
                summary=f"{_plural(rep['n_features_added'], 'new column')} appeared in {f}.",
                evidence=[ref(feature_kind(f), x, self.b, f, note=f"not in {self.a}") for x in feats],
                limitation="New columns are not used by the current analyses.",
                detected_by="feature list comparison",
                repair=("flag_for_review", [ref(feature_kind(f), x, self.b, f) for x in feats],
                        "Document the new columns before analyses use them."),
            )

    def value_changes(self, rep, f, feat_refs):
        vcs = rep.get("value_changes", [])
        lost_all = [v for v in vcs if v["all_missing_in_b"]]
        mislabel = [v for v in vcs if v not in lost_all and v.get("likely_mislabel")]
        partial = [
            v for v in vcs
            if v not in lost_all and v not in mislabel
            and v["newly_missing"] > 0 and v["n_changed"] == v["newly_missing"]
        ]
        edited = [v for v in vcs if v not in lost_all and v not in mislabel and v not in partial]

        if lost_all:
            self.add(
                kind="missing_data", file=f, severity="medium", confidence="high",
                samples={v["sample"] for v in lost_all},
                summary=(f"{_plural(len(lost_all), 'sample')} in {f} lost all {_v(len(lost_all), 'its', 'their')} values, "
                         f"though {_v(len(lost_all), 'its ID is', 'their IDs are')} still listed."),
                evidence=[ref("sample", v["sample"], self.b, f, note=f"{v['newly_missing']} values became empty") for v in lost_all],
                limitation="Empty rows can look present in a simple count.",
                detected_by="cell-by-cell comparison with the old release",
                repair=("restore_from_previous_release", [ref("sample", v["sample"], self.a, f) for v in lost_all],
                        f"Restore these rows from {self.a}."),
            )
        if partial:
            total = sum(v["newly_missing"] for v in partial)
            self.add(
                kind="missing_data", file=f, severity="medium", confidence="high",
                samples={v["sample"] for v in partial}, features=_feature_set(partial),
                summary=f"{_plural(len(partial), 'sample')} in {f} lost some values ({_plural(total, 'value')} in total became empty).",
                evidence=[
                    ref("sample", v["sample"], self.b, f,
                        note=f"{v['newly_missing']} of {v['n_features']} values became empty")
                    for v in partial
                ],
                limitation="Partial loss may only matter for analyses that read the emptied columns.",
                detected_by="cell-by-cell comparison with the old release",
                repair=("restore_from_previous_release", [ref("sample", v["sample"], self.a, f) for v in partial],
                        f"Restore the emptied values from {self.a}."),
            )
        if mislabel:
            self.add(
                kind="sample_swap", file=f, severity="high", confidence="medium",
                samples={v["sample"] for v in mislabel},
                summary=(f"{_plural(len(mislabel), 'sample')} in {f} now {_v(len(mislabel), 'looks', 'look')} like a different sample "
                         f"than before (likely mislabel with small edits)."),
                evidence=[
                    ref("sample", v["sample"], self.b, f,
                        note=(f"new profile similarity {v['best_other_corr']:.3f} with {v['best_other']}, "
                              f"only {v['self_corr']:.3f} with its own {self.a} profile"))
                    for v in mislabel
                ] + feat_refs,
                limitation="Inferred from profile similarity, not an exact match; needs human review.",
                detected_by="centered profile correlation against every sample in the old release",
                repair=("restore_from_previous_release", [ref("sample", v["sample"], self.a, f) for v in mislabel],
                        f"Restore these profiles from {self.a} until the change is explained."),
            )
        if edited:
            cells = sum(v["n_changed"] for v in edited)
            numeric = rep.get("numeric", False)

            def note(v):
                shown = "; ".join(f"{c['column']}: {c['a']} -> {c['b']}" for c in v["changed_columns"][:3])
                more = "" if v["n_changed"] <= 3 else f"; {v['n_changed'] - 3} more"
                return shown + more

            self.add(
                kind="other" if numeric else "label_conflict", file=f, severity="medium", confidence="high",
                samples={v["sample"] for v in edited}, features=_feature_set(edited),
                summary=(f"{_plural(len(edited), 'sample')} in {f} {_v(len(edited), 'has', 'have')} changed values "
                         f"({_plural(cells, 'value')} in total) that match no other sample."),
                evidence=[ref("sample", v["sample"], self.b, f, note=note(v)) for v in edited],
                limitation="Changed values may be a correction or reprocessing; the data alone cannot tell which.",
                detected_by="cell-by-cell comparison with the old release",
                repair=("flag_for_review", [ref("sample", v["sample"], self.b, f) for v in edited],
                        f"Confirm whether the new values are corrections; restore them from {self.a} if not."),
            )

    def affects(self, cand, c):
        if cand["file"] not in set(c["inputs"]):
            return False
        if cand["whole_file"]:
            return True
        used = set(c["used_features"].get(cand["file"], []))
        feats = cand["features"]
        if feats is not None and used and not (used & feats):
            return False
        if not cand["samples"]:
            return feats is not None and bool(used & feats)
        if not c["samples"]:
            return True
        return bool(cand["samples"] & c["samples"])


def build_document(diff, rel_a, rel_b, a_name, b_name, impact,
                   team="SNS", product="Cohort Control Tower", dataset="TCGA-BRCA",
                   extra_files_not_in_release=None):
    b = Builder(diff, rel_a, rel_b, a_name, b_name, impact)
    if extra_files_not_in_release:
        b.files_not_in_release.extend(extra_files_not_in_release)
    b.file_level()
    for rep in diff["layers"]:
        if not rep["content_same"]:
            b.layer_level(rep)

    findings, repairs = [], []
    for i, cand in enumerate(b.cands, 1):
        fid = f"F{i}"
        aff = [aid for aid, c in b.ctx.items() if b.affects(cand, c)]
        sev = cand["severity"]
        if not aff and sev == "high":
            sev = "medium"
        findings.append({
            "id": fid,
            "kind": cand["kind"],
            "summary": cand["summary"],
            "severity": sev,
            "confidence": cand["confidence"],
            "evidence": cand["evidence"],
            "affects_analyses": aff,
            "not_affected_analyses": [aid for aid in b.ctx if aid not in aff],
            "limitation": cand["limitation"],
            "needs_human_review": sev in ("high", "medium"),
            "detected_by": cand["detected_by"],
        })
        action, targets, why = cand["repair"]
        repairs.append({
            "id": f"R{len(repairs) + 1}",
            "action": action,
            "targets": targets,
            "justification": why,
            "addresses": [fid],
        })

    analyses, verdicts = [], []
    for aid, c in b.ctx.items():
        by = [f["id"] for f in findings if aid in f["affects_analyses"]]
        before = impact["before"].get(aid, {})
        after = impact["after"].get(aid, {})
        rep = impact["repaired"].get(aid, {})
        verdict, reason = verdict_for(by, before, after)
        verdicts.append(verdict)
        entry = {
            "id": aid,
            "status": "changed" if public(before) != public(after) else "unchanged",
            "before": public(before),
            "after": public(after),
            "after_repair": public(rep),
            "repair_restores_before": public(rep) == public(before),
            "description": c["description"],
            "inputs": c["inputs"],
            "uses_columns": c["used_features"],
            "affected_by": by,
            "verdict": verdict,
            "reason": reason,
        }
        if not by:
            entry["note"] = "not affected by any finding; rerun only to confirm it did not move"
        analyses.append(entry)

    if verdicts:
        overall = max(verdicts, key=RANK.get)
    else:
        overall = "Conditional" if findings else "Ready"

    def files(rel):
        return [{"file": n, "sha256": layer.sha256, "hash_of": "content"} for n, layer in sorted(rel.items())]

    asked = _join(f"{aid} ({c['description'].rstrip('.')})" for aid, c in b.ctx.items())
    return {
        "format_version": "0.2",
        "producer": {"team": team, "product": product, "track": "cohort-firewall"},
        "data": {
            "dataset": dataset,
            "releases": [
                {"name": a_name, "files": files(rel_a)},
                {"name": b_name, "files": files(rel_b)},
            ],
        },
        "question": f"Did release {b_name} change anything that affects analyses {asked}?",
        "verdict": overall,
        "findings": findings,
        "analyses": analyses,
        "repairs": repairs,
        "sns_run": {
            "engine_version": ENGINE_VERSION,
            "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "subgroup_checks": b.subgroup_checks,
            "files_not_in_release": b.files_not_in_release,
            "files_not_in_release_note": (
                "not provided in this release; no analysis reads it; not compared"
                if b.files_not_in_release else None
            ),
        },
    }