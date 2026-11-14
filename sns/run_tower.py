import argparse
import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from engine.diff import compare_releases  # noqa: E402
from engine.findings import build_document  # noqa: E402
from engine.impact import assess  # noqa: E402
from engine.release import load_release  # noqa: E402

DATA = ROOT / "data" / "tcga-brca"


def quick_scan_analysis_files(analysis_path):
    text = Path(analysis_path).read_text(encoding="utf-8")
    return set(re.findall(r'\b[\w-]+\.tsv(?:\.gz)?\b', text))


def main():
    t0 = time.time()
    p = argparse.ArgumentParser()
    p.add_argument("--a", default=str(DATA / "raw"))
    p.add_argument("--b", default=str(DATA / "releases" / "practice-B"))
    p.add_argument("--analyses", default=str(DATA / "reference_analyses.py"))
    p.add_argument("--out", default="findings.json")
    p.add_argument("--diff-out", default=None)
    p.add_argument("--team", default="SNS")
    p.add_argument("--a-name", default=None)
    p.add_argument("--b-name", default=None)
    p.add_argument("--pretty", action="store_true",
                   help="Print a human-readable summary after writing the file.")
    args = p.parse_args()

    folder_a, folder_b = Path(args.a), Path(args.b)

    # Determine b_name: if folder is named "release", use parent name (envelope convention)
    if args.b_name:
        b_name = args.b_name
    elif folder_b.name == "release":
        b_name = folder_b.parent.name
    else:
        b_name = folder_b.name

    # Skip files from A that are absent in B and not referenced by any analysis
    analysis_files = quick_scan_analysis_files(args.analyses)
    files_in_b = {p.name for p in folder_b.iterdir() if p.is_file()}
    files_in_a_folder = {p.name for p in folder_a.iterdir() if p.is_file()}
    skip_from_a = {f for f in files_in_a_folder - files_in_b if f not in analysis_files}

    rel_a = load_release(folder_a, skip_files=skip_from_a)
    rel_b = load_release(folder_b, reuse=rel_a)
    diff = compare_releases(rel_a, rel_b)
    if args.diff_out:
        Path(args.diff_out).write_text(json.dumps(diff, indent=2, default=str), encoding="utf-8")

    impact = assess(args.analyses, rel_a, rel_b, folder_a, folder_b, diff)
    doc = build_document(
        diff, rel_a, rel_b,
        args.a_name or folder_a.name,
        b_name,
        impact, team=args.team,
        extra_files_not_in_release=sorted(skip_from_a),
    )
    out = Path(args.out)
    out.write_text(json.dumps(doc, indent=2, default=str), encoding="utf-8")

    elapsed = round(time.time() - t0, 1)

    if args.pretty:
        _print_pretty(doc, out, elapsed)
    else:
        print(json.dumps({
            "out": str(out),
            "verdict": doc["verdict"],
            "findings": [
                {"id": f["id"], "kind": f["kind"], "severity": f["severity"],
                 "n_evidence": len(f["evidence"]), "affects": f["affects_analyses"]}
                for f in doc["findings"]
            ],
            "analyses": [
                {"id": a["id"], "status": a["status"], "verdict": a["verdict"],
                 "repair_restores_before": a["repair_restores_before"]}
                for a in doc["analyses"]
            ],
        }, indent=2))


def _print_pretty(doc, out, elapsed):
    print(f"verdict: {doc['verdict']}")
    for f in doc["findings"]:
        aff = ", ".join(f["affects_analyses"]) or "none"
        naf = ", ".join(f["not_affected_analyses"]) or "none"
        print(f"  {f['id']} [{f['severity']}] {f['kind']}: {f['summary']}")
        print(f"    affects: {aff}  not affected: {naf}")
        review = "yes" if f["needs_human_review"] else "no"
        print(f"    needs human review: {review}")
    repairs_map = {}
    for a in doc["analyses"]:
        repairs_map[a["id"]] = a["repair_restores_before"]
    repair_line = ", ".join(
        f"{a['id']} {'yes' if a['repair_restores_before'] else 'no'}"
        for a in doc["analyses"]
    )
    print(f"repair restores original results: {repair_line}")
    for a in doc["analyses"]:
        status = a["status"]
        verdict = a["verdict"]
        reason = a.get("reason", "")
        print(f"  {a['id']} {status} | verdict {verdict} | {reason}")
    not_in_release = doc.get("sns_run", {}).get("files_not_in_release", [])
    if not_in_release:
        print(f"files not in release: {', '.join(not_in_release)}")
    print(f"output: {out}  ({elapsed}s)")


if __name__ == "__main__":
    main()
