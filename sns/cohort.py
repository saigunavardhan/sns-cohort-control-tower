"""Clean Cohort Export: reads findings.json and builds a safe, repaired sample list."""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from engine.cohort import BUILTIN_SPECS, build_cohort_output, run_cohort  # noqa: E402

DATA = ROOT / "data" / "tcga-brca"
SLIDES = DATA / "slide_manifest.csv"


def _make_reproduce_script(out_dir, folder_b, folder_a, spec, applied_repairs):
    """Write a standalone pandas script that rebuilds cohort_samples.csv."""
    lines = [
        '"""Standalone script: rebuilds cohort_samples.csv from the release folder."""',
        "import pandas as pd",
        "from pathlib import Path",
        "",
        f'RELEASE_B = Path({repr(str(folder_b))})',
        f'RELEASE_A = Path({repr(str(folder_a))})',
        "",
        "SPEC = " + json.dumps(spec, indent=4),
        "",
        "REPAIRS = " + json.dumps(applied_repairs, indent=4),
        "",
        "# --- helper ---",
        "def _load(folder, name):",
        "    p = Path(folder) / name",
        "    sep = ',' if '.csv' in p.name.lower() else '\\t'",
        "    df = pd.read_csv(p, sep=sep, index_col=0)",
        "    tcga_cols = [c for c in df.columns if str(c).startswith('TCGA-')]",
        "    if len(tcga_cols) > 10:",
        "        return df[tcga_cols].T",
        "    id_col = next((c for c in df.reset_index().columns if",
        "                   df.reset_index()[c].astype(str).str.startswith('TCGA-').mean() > 0.9), None)",
        "    return df",
        "",
        "# --- build sample set ---",
        "working = set()",
        "for p in sorted(RELEASE_B.iterdir()):",
        "    if p.suffix in ('.tsv', '.gz', '.csv'):",
        "        try:",
        "            sep = ',' if '.csv' in p.name.lower() else '\\t'",
        "            df = pd.read_csv(p, sep=sep, nrows=0)",
        "            tcga_cols = [c for c in df.columns if str(c).startswith('TCGA-')]",
        "            if tcga_cols:",
        "                full = pd.read_csv(p, sep=sep, index_col=0)",
        "                working.update(str(c) for c in full.columns if str(c).startswith('TCGA-'))",
        "            else:",
        "                for col in df.columns:",
        "                    full = pd.read_csv(p, sep=sep)",
        "                    if full[col].astype(str).str.startswith('TCGA-').mean() > 0.9:",
        "                        working.update(str(v) for v in full[col].dropna())",
        "                        break",
        "        except Exception:",
        "            pass",
        "",
        "flagged = set()",
        "",
        "# Apply repairs",
        "for r in REPAIRS:",
        "    if r['action'] == 'exclude':",
        "        working.discard(r['sample'])",
        "    elif r['action'] in ('flag_for_review', 'flag_for_review_restored'):",
        "        if r.get('action') == 'flag_for_review_restored':",
        "            working.add(r['sample'])",
        "        flagged.add(r['sample'])",
        "    elif r['action'] == 'relabel':",
        "        working.discard(r['sample'])",
        "        working.add(r['renamed_to'])",
        "",
        "# Apply spec filters",
        "clinical_path = next((RELEASE_B / f for f in ['clinical.tsv'] if (RELEASE_B / f).exists()), None)",
        "clinical = pd.read_csv(clinical_path, sep='\\t') if clinical_path else None",
        "if clinical is not None:",
        "    clinical = clinical.set_index('sampleID') if 'sampleID' in clinical.columns else clinical",
        "",
        "for step in SPEC:",
        "    if step['step'] == 'require_file':",
        "        fname = step['file']",
        "        p = RELEASE_B / fname",
        "        if p.exists():",
        "            sep = ',' if '.csv' in fname.lower() else '\\t'",
        "            df = pd.read_csv(p, sep=sep, nrows=0)",
        "            tcga_cols = [c for c in df.columns if str(c).startswith('TCGA-')]",
        "            if tcga_cols:",
        "                full = pd.read_csv(p, sep=sep, index_col=0)",
        "                file_sids = set(str(c) for c in full.columns if str(c).startswith('TCGA-'))",
        "            else:",
        "                full = pd.read_csv(p, sep=sep)",
        "                id_col = next((c for c in full.columns if",
        "                               full[c].astype(str).str.startswith('TCGA-').mean() > 0.9), None)",
        "                file_sids = set(str(v) for v in (full[id_col].dropna() if id_col else []))",
        "            # Add back flagged samples from A that were dropped",
        "            for r in REPAIRS:",
        "                if r['action'] == 'flag_for_review_restored' and r.get('file') == fname:",
        "                    file_sids.add(r['sample'])",
        "            # Remove excluded from file set",
        "            for r in REPAIRS:",
        "                if r['action'] == 'exclude' and r.get('file') == fname:",
        "                    file_sids.discard(r['sample'])",
        "            working &= file_sids",
        "            # Re-add flagged that were added",
        "            for r in REPAIRS:",
        "                if r['action'] == 'flag_for_review_restored' and r.get('file') == fname:",
        "                    if r['sample'] in file_sids:",
        "                        working.add(r['sample'])",
        "    elif step['step'] == 'filter_sample_type':",
        "        t = step.get('type', '01')",
        "        working = {s for s in working if str(s)[13:15] == t}",
        "    elif step['step'] == 'require_value':",
        "        col = step['column']",
        "        if clinical is not None and col in clinical.columns:",
        "            if step.get('not_empty'):",
        "                working = {s for s in working if s in clinical.index and pd.notna(clinical.at[s, col])}",
        "            else:",
        "                allowed = set(step['allowed'])",
        "                working = {s for s in working if s in clinical.index and clinical.at[s, col] in allowed}",
        "",
        "rows = [{'sample_id': s, 'patient_id': s[:12], 'flagged': s in flagged} for s in sorted(working)]",
        "out = pd.DataFrame(rows)",
        "out.to_csv('cohort_samples.csv', index=False)",
        f"print(f'cohort_samples.csv: {{len(out)}} samples')",
    ]
    path = Path(out_dir) / "reproduce_cohort.py"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    p = argparse.ArgumentParser(description="Build a clean cohort from a Control Tower findings.json.")
    p.add_argument("--release", required=True, help="Release B folder")
    p.add_argument("--prev", default=str(DATA / "raw"), help="Release A folder (default: raw)")
    p.add_argument("--findings", required=True, help="findings.json from run_tower.py")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--analysis", choices=list(BUILTIN_SPECS.keys()),
                       help="Built-in analysis spec (A1 or A2)")
    group.add_argument("--spec", help="Custom spec JSON file")
    p.add_argument("--slides", default=str(SLIDES) if SLIDES.exists() else None,
                   help="Slide manifest CSV (optional)")
    p.add_argument("--out-dir", required=True, help="Output folder")
    args = p.parse_args()

    folder_b = Path(args.release)
    folder_a = Path(args.prev)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    findings_doc = json.loads(Path(args.findings).read_text(encoding="utf-8"))

    if args.analysis:
        spec = BUILTIN_SPECS[args.analysis]
        analysis_id = args.analysis
    else:
        spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
        analysis_id = "custom"

    working, flagged, steps, applied_repairs = run_cohort(
        spec, folder_b, folder_a, findings_doc,
        slides_path=args.slides,
    )

    df = build_cohort_output(working, flagged)

    # Duplicate check
    if df["sample_id"].duplicated().any():
        dupes = df[df["sample_id"].duplicated(keep=False)]["sample_id"].tolist()
        print(f"WARNING: duplicate sample_ids detected: {dupes[:5]}", flush=True)

    csv_path = out_dir / "cohort_samples.csv"
    df.to_csv(csv_path, index=False)

    steps_data = {
        "analysis": analysis_id,
        "release_b": str(folder_b),
        "release_a": str(folder_a),
        "steps": steps,
        "repairs_applied": applied_repairs,
        "n_samples": len(df),
        "n_flagged": int(df["flagged"].sum()),
        "n_patients": int(df["patient_id"].nunique()),
    }
    steps_path = out_dir / "cohort_steps.json"
    steps_path.write_text(json.dumps(steps_data, indent=2, default=str), encoding="utf-8")

    reproduce_path = _make_reproduce_script(
        out_dir, folder_b, folder_a, spec, applied_repairs
    )

    print(json.dumps({
        "analysis": analysis_id,
        "n_samples": len(df),
        "n_flagged": int(df["flagged"].sum()),
        "n_patients": int(df["patient_id"].nunique()),
        "out_dir": str(out_dir),
        "files": [str(csv_path.name), str(steps_path.name), str(reproduce_path.name)],
    }, indent=2))


if __name__ == "__main__":
    main()
