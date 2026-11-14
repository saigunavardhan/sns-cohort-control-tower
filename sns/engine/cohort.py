"""Cohort export engine: reads findings.json, applies repairs, runs spec steps."""
import json
from pathlib import Path

import pandas as pd

from .release import TCGA_PREFIX, load_layer

BUILTIN_SPECS = {
    "A1": [
        {"step": "require_file", "file": "protein_rppa.tsv.gz"},
        {"step": "filter_sample_type", "type": "01"},
        {"step": "require_value", "column": "breast_carcinoma_estrogen_receptor_status",
         "allowed": ["Positive", "Negative"]},
    ],
    "A2": [
        {"step": "require_file", "file": "clinical.tsv"},
        {"step": "filter_sample_type", "type": "01"},
        {"step": "require_value", "column": "PAM50Call_RNAseq", "not_empty": True},
    ],
}


def _patient(sid):
    return str(sid)[:12]


def _file_samples(folder, name):
    path = Path(folder) / name
    if not path.exists():
        return set()
    layer = load_layer(path)
    return set(layer.table.index) if layer else set()


def _load_clinical(folder):
    for p in sorted(Path(folder).iterdir()):
        if "clinical" in p.name.lower() and p.suffix in (".tsv", ".csv"):
            layer = load_layer(p)
            if layer is not None:
                return layer.table
    return None


def parse_repairs(findings_doc):
    """Return dict: file -> {excluded, flagged, restored, relabeled}."""
    result = {}
    for rep in findings_doc.get("repairs", []):
        action = rep["action"]
        for t in rep["targets"]:
            if t["kind"] != "sample":
                continue
            sid = str(t["id"])
            file = t.get("file") or ""
            if file not in result:
                result[file] = {"excluded": set(), "flagged": set(),
                                "restored": set(), "relabeled": {}}
            if action == "exclude":
                result[file]["excluded"].add(sid)
            elif action == "restore_from_previous_release":
                result[file]["restored"].add(sid)
            elif action == "flag_for_review":
                result[file]["flagged"].add(sid)
            elif action == "relabel":
                note = t.get("note", "")
                if "rename back to" in note:
                    old_id = note.split("rename back to")[-1].strip()
                    result[file]["relabeled"][sid] = old_id
    return result


def _apply_file_repairs(samples_b, folder_a, file_name, repairs):
    """Return (repaired_samples: set, flagged: set, applied: list)."""
    fr = repairs.get(file_name, {})
    excluded = fr.get("excluded", set())
    flagged_targets = fr.get("flagged", set())
    restored = fr.get("restored", set())
    relabeled = fr.get("relabeled", {})

    working = set(samples_b)
    applied = []
    flagged_out = set()

    for sid in sorted(excluded):
        if sid in working:
            working.discard(sid)
            applied.append({"action": "exclude", "sample": sid, "file": file_name})

    # Flag-for-review: if sample is still in working, mark it; if dropped from B, add from A
    samples_a = _file_samples(folder_a, file_name) if folder_a else set()
    for sid in sorted(flagged_targets):
        if sid in working:
            flagged_out.add(sid)
            applied.append({"action": "flag_for_review", "sample": sid, "file": file_name})
        elif sid in samples_a:
            working.add(sid)
            flagged_out.add(sid)
            applied.append({"action": "flag_for_review_restored", "sample": sid, "file": file_name})

    for sid in sorted(restored):
        if sid in working:
            applied.append({"action": "restore_from_previous_release", "sample": sid, "file": file_name})

    for new_id, old_id in sorted(relabeled.items()):
        if new_id in working:
            working.discard(new_id)
            working.add(old_id)
            if new_id in flagged_out:
                flagged_out.discard(new_id)
                flagged_out.add(old_id)
            applied.append({"action": "relabel", "sample": new_id, "renamed_to": old_id})

    return working, flagged_out, applied


def _step_counts(sample_set):
    ids = sorted(sample_set)
    n_distinct = len(set(ids))
    n_patients = len({_patient(s) for s in ids})
    return {"n_rows": len(ids), "n_distinct_samples": n_distinct, "n_distinct_patients": n_patients}


def run_cohort(spec, folder_b, folder_a, findings_doc, slides_path=None):
    """Execute spec steps. Return (samples_df, steps_log, applied_repairs)."""
    repairs = parse_repairs(findings_doc)
    all_applied = []

    # Collect all samples from B to start
    all_b = set()
    for p in sorted(Path(folder_b).iterdir()):
        if p.is_file() and p.name.lower().endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz")):
            layer = load_layer(p)
            if layer:
                all_b.update(layer.table.index)

    working = set(all_b)
    flagged = set()
    steps = []

    clinical = _load_clinical(folder_b)
    if clinical is None and folder_a:
        clinical = _load_clinical(folder_a)

    for step_def in spec:
        step_name = step_def["step"]
        step_rec = {"step": step_name, "params": {k: v for k, v in step_def.items() if k != "step"}}

        if step_name == "require_file":
            fname = step_def["file"]
            step_rec["params"]["file"] = fname
            file_b = _file_samples(folder_b, fname)
            repaired, new_flagged, applied = _apply_file_repairs(file_b, folder_a, fname, repairs)
            all_applied.extend(applied)
            flagged |= new_flagged
            working &= repaired
            # Flagged samples restored from A may not be in working after the intersect
            # Add them back if they're in repaired (which they are)
            for sid in new_flagged:
                if sid in repaired:
                    working.add(sid)

        elif step_name == "filter_sample_type":
            typ = step_def.get("type", "01")
            working = {s for s in working if str(s)[13:15] == typ}

        elif step_name == "require_value":
            col = step_def["column"]
            allowed = step_def.get("allowed")
            not_empty = step_def.get("not_empty", False)
            if clinical is not None and col in clinical.columns:
                clin_idx = clinical[~clinical.index.duplicated()]
                if not_empty:
                    valid = {s for s in working if s in clin_idx.index and pd.notna(clin_idx.at[s, col])}
                else:
                    valid = {s for s in working if s in clin_idx.index and clin_idx.at[s, col] in allowed}
                working = valid

        elif step_name == "join_slides":
            if slides_path and Path(slides_path).exists():
                slides = pd.read_csv(slides_path)
                slide_kind = step_def.get("slide_kind")
                if slide_kind:
                    slides = slides[slides["slide_kind"].str.contains(slide_kind, case=False, na=False)]
                slides = slides[slides["sample_id"].isin(working)]
                n_before = len(working)
                n_after_rows = len(slides)
                n_after_distinct = slides["sample_id"].nunique()
                inflated = n_after_rows > n_after_distinct
                if inflated:
                    top = slides["sample_id"].value_counts().head(5).to_dict()
                    step_rec["inflation"] = {"inflated": True, "top_samples": top}
                else:
                    step_rec["inflation"] = {"inflated": False}
                step_rec["n_slides_matched"] = n_after_rows
                # Collapse back to one row per sample
                working = set(slides["sample_id"].tolist())

        counts = _step_counts(working)
        step_rec.update(counts)
        steps.append(step_rec)

    return working, flagged, steps, all_applied


def build_cohort_output(working, flagged):
    """Build the cohort DataFrame."""
    rows = []
    for sid in sorted(working):
        rows.append({
            "sample_id": sid,
            "patient_id": _patient(sid),
            "flagged": sid in flagged,
        })
    return pd.DataFrame(rows)
