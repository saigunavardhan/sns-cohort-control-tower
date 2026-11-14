import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / "data" / "tcga-brca" / "raw"
RUNS = HERE / "runs"
UPLOADS = HERE / "uploads"
PY = sys.executable
ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}

VERDICT_STYLE = {
    "Ready": st.success,
    "Conditional": st.warning,
    "Blocked": st.error,
    "Insufficient Evidence": st.info,
}

st.set_page_config(page_title="SNS Cohort Control Tower", layout="wide")


def run_cmd(cmd):
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=ENV)


def find_releases():
    found = {}
    opened = ROOT / "envelopes" / "opened"
    if opened.exists():
        for p in sorted(opened.glob("*/release")):
            found[f"envelope: {p.parent.name}"] = p
    rel = ROOT / "data" / "tcga-brca" / "releases"
    if rel.exists():
        for p in sorted(rel.iterdir()):
            if p.is_dir():
                found[f"release: {p.name}"] = p
    if UPLOADS.exists():
        for p in sorted(UPLOADS.iterdir()):
            if p.is_dir():
                found[f"uploaded: {p.name}"] = p
    return found


def flat(d, prefix=""):
    out = {}
    for k, v in (d or {}).items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flat(v, key + "."))
        elif not isinstance(v, list):
            out[key] = v
    return out


st.title("SNS Cohort Control Tower")
st.caption("Compares a trusted release with a new release, finds silent changes, "
           "shows which analyses they affect, and proposes repairs. Research use only.")

with st.sidebar:
    st.header("1. Get a new release")
    tab_env, tab_up = st.tabs(["Sealed envelope", "Upload files"])

    with tab_env:
        env_no = st.text_input("Envelope number", placeholder="07 or practice")
        env_words = st.text_input("Four words", placeholder="amber falcon river quartz")
        if st.button("Unseal"):
            words = env_words.split()
            if not env_no.strip() or len(words) != 4:
                st.error("Enter an envelope number and exactly four words.")
            else:
                r = run_cmd([PY, "envelopes/unseal.py", env_no.strip(), *words])
                if r.returncode == 0:
                    st.success("Envelope opened.")
                else:
                    st.error("Unseal failed.")
                st.code((r.stdout + r.stderr)[-3000:])

    with tab_up:
        up_name = st.text_input("Name for this release", placeholder="my-release")
        files = st.file_uploader("Release files (.tsv, .tsv.gz, .csv)",
                                 accept_multiple_files=True)
        if st.button("Save upload"):
            if not up_name.strip() or not files:
                st.error("Give a name and at least one file.")
            else:
                target = UPLOADS / up_name.strip()
                target.mkdir(parents=True, exist_ok=True)
                for f in files:
                    (target / f.name).write_bytes(f.getbuffer())
                st.success(f"Saved {len(files)} files to {target.relative_to(ROOT)}")

    st.header("2. Choose and run")
    releases = find_releases()
    if not releases:
        st.info("No releases found yet. Unseal an envelope or upload files.")
    choice = st.selectbox("New release (B)", list(releases)) if releases else None
    old = st.text_input("Trusted release (A)", value=str(RAW))
    run_clicked = st.button("Run check", type="primary", disabled=choice is None)

if run_clicked and choice:
    b_path = releases[choice]
    RUNS.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = choice.split(": ", 1)[-1].replace(" ", "_")
    out = RUNS / f"findings_{label}_{stamp}.json"
    t0 = time.time()
    with st.spinner("Comparing releases and rerunning analyses..."):
        r = run_cmd([PY, "sns/run_tower.py", "--a", old, "--b", str(b_path), "--out", str(out)])
    if r.returncode != 0 or not out.exists():
        st.error("The run failed.")
        st.code((r.stdout + r.stderr)[-4000:])
    else:
        st.session_state.update({"out": str(out), "b": str(b_path), "a": old,
                                 "label": label, "seconds": round(time.time() - t0, 1)})

if "out" not in st.session_state:
    st.info("Pick a release in the sidebar and press Run check.")
    st.stop()

out = Path(st.session_state["out"])
doc = json.loads(out.read_text(encoding="utf-8"))
sns_run = doc.get("sns_run", {})

verdict = doc.get("verdict", "Unknown")
VERDICT_STYLE.get(verdict, st.info)(f"### Verdict: {verdict}")
st.write(doc.get("question", ""))
st.caption(f"Release A: {st.session_state['a']}  |  Release B: {st.session_state['b']}  |  "
           f"{st.session_state['seconds']} s  |  saved as {out.relative_to(ROOT)}")

st.header("Analyses")
for a in doc.get("analyses", []):
    st.subheader(f"{a['id']}: {a.get('description', '')}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Status", a.get("status", ""))
    c2.metric("Verdict", a.get("verdict", ""))
    c3.metric("Repair restores original", "yes" if a.get("repair_restores_before") else "no")
    st.write(a.get("reason", ""))
    keys = sorted(set(flat(a.get("before"))) | set(flat(a.get("after"))) | set(flat(a.get("after_repair"))))
    if keys:
        fb, fa, fr = flat(a.get("before")), flat(a.get("after")), flat(a.get("after_repair"))
        st.dataframe(pd.DataFrame({
            "value": keys,
            "old release": [fb.get(k) for k in keys],
            "new release": [fa.get(k) for k in keys],
            "after repair": [fr.get(k) for k in keys],
        }), hide_index=True, use_container_width=True)
    if a.get("affected_by"):
        st.caption("Affected by: " + ", ".join(a["affected_by"]))

st.header(f"Findings ({len(doc.get('findings', []))})")
if not doc.get("findings"):
    st.success("No changes found that need review.")
for f in doc.get("findings", []):
    title = f"{f['id']} [{f['severity']}] {f.get('kind', '')}: {f['summary']}"
    with st.expander(title, expanded=True):
        st.write(f"**Affects:** {', '.join(f.get('affects_analyses', [])) or 'none'}  |  "
                 f"**Not affected:** {', '.join(f.get('not_affected_analyses', [])) or 'none'}  |  "
                 f"**Confidence:** {f.get('confidence', '')}  |  "
                 f"**Needs human review:** {'yes' if f.get('needs_human_review') else 'no'}")
        if f.get("limitation"):
            st.write(f"**Limitation:** {f['limitation']}")
        if f.get("detected_by"):
            st.caption(f"Detected by: {f['detected_by']}")
        st.dataframe(pd.DataFrame(f.get("evidence", [])), hide_index=True, use_container_width=True)

st.header("Proposed repairs")
reps = doc.get("repairs", [])
if reps:
    st.dataframe(pd.DataFrame([{
        "id": r["id"],
        "action": r["action"],
        "fixes": ", ".join(r.get("addresses", [])),
        "targets": len(r.get("targets", [])),
        "justification": r.get("justification", ""),
    } for r in reps]), hide_index=True, use_container_width=True)
    st.caption("Repairs are proposals. A person must approve or reject each one.")
else:
    st.write("No repairs proposed.")

missing = sns_run.get("files_not_in_release") or []
if missing:
    st.info("Files not in this release (no analysis reads them, so they were not compared): "
            + ", ".join(missing))

st.header("Limitations")
st.markdown(
    "- The tool sees what changed between releases, not why it changed.\n"
    "- It treats release A as the trusted version. Every finding needs human review.\n"
    "- Tables only: slide images are not checked.\n"
    "- Sample IDs must start with TCGA-.\n"
    "- In sparse files (like mutation calls), a near copy of a sample with very few "
    "non-zero values cannot be told apart from a real sample.\n"
    "- Research use only. No diagnosis, treatment advice or clinical claims."
)

st.header("Downloads and clean cohort")
st.download_button("Download findings.json", out.read_bytes(), file_name=f"findings_{st.session_state['label']}.json",
                   mime="application/json")

analysis = st.selectbox("Clean cohort for analysis", [a["id"] for a in doc.get("analyses", [])] or ["A1"])
if st.button("Build clean cohort"):
    cdir = RUNS / f"cohort_{st.session_state['label']}_{analysis}_{datetime.now().strftime('%H%M%S')}"
    r = run_cmd([PY, "sns/cohort.py", "--release", st.session_state["b"], "--prev", st.session_state["a"],
                 "--findings", str(out), "--analysis", analysis, "--out-dir", str(cdir)])
    csv = cdir / "cohort_samples.csv"
    if r.returncode == 0 and csv.exists():
        cohort = pd.read_csv(csv)
        st.success(f"{len(cohort)} samples, {cohort['sample_id'].duplicated().sum()} duplicates")
        st.dataframe(cohort.head(50), hide_index=True, use_container_width=True)
        st.download_button("Download cohort_samples.csv", csv.read_bytes(),
                           file_name=f"cohort_samples_{analysis}.csv", mime="text/csv")
    else:
        st.error("Cohort export failed.")
        st.code((r.stdout + r.stderr)[-3000:])