import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from engine.release import load_layer  # noqa: E402

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
KIND_COLOR = {
    "sample_swap": "#e63946",
    "duplicate": "#f4a261",
    "membership_change": "#7b2cbf",
    "missing_data": "#457b9d",
    "label_conflict": "#2a9d8f",
    "other": "#8d99ae",
}
KIND_LABEL = {
    "sample_swap": "swapped",
    "duplicate": "copies",
    "membership_change": "lost or renamed",
    "missing_data": "values missing",
    "label_conflict": "labels changed",
    "other": "values changed",
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


@st.cache_data(show_spinner=False)
def layer_table(path_str, mtime):
    layer = load_layer(Path(path_str))
    if layer is None:
        return None
    t = layer.table
    return t[~t.index.duplicated()]


def get_table(folder, name):
    p = Path(folder) / name
    if not p.exists():
        return None
    return layer_table(str(p), p.stat().st_mtime)


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


def numeric(t):
    return t.apply(pd.to_numeric, errors="coerce")


# ---------- visual tabs ----------

def tab_pipeline(doc, a_dir, b_dir):
    fa, fb = release_files(doc)
    changed = [n for n in fb if n in fa and fa[n] != fb[n]]
    n_rows = 0
    for n in changed:
        t = get_table(a_dir, n)
        n_rows += 0 if t is None else len(t)
    ev = evidence_by_file(doc)
    changed_samples = set()
    for kinds in ev.values():
        for s in kinds.values():
            changed_samples |= s
    findings = doc.get("findings", [])
    subs = doc.get("sns_run", {}).get("subgroup_checks", [])
    analyses = doc.get("analyses", [])
    affected = [a["id"] for a in analyses if a.get("affected_by")]
    safe = [a["id"] for a in analyses if not a.get("affected_by")]
    restored = [a["id"] for a in analyses if a.get("affected_by") and a.get("repair_restores_before")]

    c = st.columns(6)
    c[0].metric("1. Files fingerprinted", f"{len(fb)}", f"{len(changed)} changed", delta_color="inverse")
    c[1].metric("2. Rows fingerprinted", f"{n_rows}", f"{len(changed_samples)} samples changed", delta_color="inverse")
    c[2].metric("3. Problems named", f"{len(findings)}",
                ", ".join(sorted({KIND_LABEL.get(f.get('kind'), f.get('kind', '')) for f in findings})) or "none",
                delta_color="off")
    if subs:
        h = subs[0]
        c[3].metric("4. Group hit", f"{h['value']}",
                    f"{h.get('with_data_in_a', '?')} to {h.get('with_data_in_b', '?')}", delta_color="inverse")
    else:
        c[3].metric("4. Group hit", "none", "no single group", delta_color="off")
    c[4].metric("5. Analyses affected", ", ".join(affected) or "none",
                f"safe: {', '.join(safe) or 'none'}", delta_color="off")
    c[5].metric("6. Repair proof", "yes" if affected and len(restored) == len(affected) else ("n/a" if not affected else "partial"),
                f"restored: {', '.join(restored) or 'none'}", delta_color="off")
    st.caption("1: a sha256 code per file. 2: a sha256 code per sample row. 3: codes matched between releases. "
               "4: lost samples looked up in clinical data. 5: analysis code read to see which files and columns it uses. "
               "6: analyses rerun on the old, new and repaired data.")


def main_file(doc):
    ev = evidence_by_file(doc)
    if not ev:
        return None
    return max(ev, key=lambda f: sum(len(s) for s in ev[f].values()))


def tab_flow(doc, a_dir, b_dir):
    f = main_file(doc)
    if not f:
        st.info("No sample-level changes to show.")
        return
    ta, tb = get_table(a_dir, f), get_table(b_dir, f)
    if ta is None or tb is None:
        st.info(f"{f} is not in both releases.")
        return
    ids_a, ids_b = set(ta.index), set(tb.index)
    ev = evidence_by_file(doc).get(f, {})
    lost = ids_a - ids_b
    new = ids_b - ids_a
    both = ids_a & ids_b
    groups = {k: v & both for k, v in ev.items() if v & both}
    flagged = set().union(*groups.values()) if groups else set()
    unchanged = both - flagged

    labels = [f"Release A: {len(ids_a)} samples", f"Release B: {len(ids_b)} samples",
              f"unchanged: {len(unchanged)}"]
    colors = ["#264653", "#264653", "#adb5bd"]
    src, tgt, val, lcol = [0, 2], [2, 1], [len(unchanged), len(unchanged)], ["rgba(173,181,189,0.4)"] * 2
    for k, v in groups.items():
        labels.append(f"{KIND_LABEL.get(k, k)}: {len(v)}")
        colors.append(KIND_COLOR.get(k, "#8d99ae"))
        i = len(labels) - 1
        src += [0, i]
        tgt += [i, 1]
        val += [len(v), len(v)]
        lcol += ["rgba(230,57,70,0.5)"] * 2
    if lost:
        labels.append(f"lost: {len(lost)}")
        colors.append(KIND_COLOR["membership_change"])
        src.append(0)
        tgt.append(len(labels) - 1)
        val.append(len(lost))
        lcol.append("rgba(123,44,191,0.5)")
    if new:
        labels.append(f"new IDs: {len(new)}")
        colors.append(KIND_COLOR["duplicate"])
        src.append(len(labels) - 1)
        tgt.append(1)
        val.append(len(new))
        lcol.append("rgba(244,162,97,0.6)")
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(label=labels, color=colors, pad=18, thickness=22),
        link=dict(source=src, target=tgt, value=val, color=lcol),
    ))
    fig.update_layout(height=460, margin=dict(l=10, r=10, t=30, b=10),
                      title=f"Where every sample of {f} went between the two releases")
    st.plotly_chart(fig, use_container_width=True)
    st.caption(f"Release A has {len(ids_a)} samples and release B has {len(ids_b)}. "
               f"Same count or not, {len(lost)} left, {len(new)} new IDs came in, and {len(flagged)} kept their ID but changed.")


def pca3(ta, tb, max_features=2000):
    xa = numeric(ta)
    xa = xa.loc[:, xa.notna().mean() > 0.5]
    if xa.shape[1] > max_features:
        xa = xa[xa.var().sort_values(ascending=False).index[:max_features]]
    mu = xa.mean()
    sd = xa.std().replace(0, 1).fillna(1)
    za = ((xa - mu) / sd).fillna(0).to_numpy()
    _, _, vt = np.linalg.svd(za, full_matrices=False)
    comps = vt[:3].T
    xb = numeric(tb.reindex(columns=xa.columns))
    zb = ((xb - mu) / sd).fillna(0).to_numpy()
    pa = pd.DataFrame(za @ comps, index=xa.index, columns=["x", "y", "z"])
    pb = pd.DataFrame(zb @ comps, index=xb.index, columns=["x", "y", "z"])
    return pa, pb


def color_column(doc, f, clin):
    for a in doc.get("analyses", []):
        if f in a.get("inputs", []):
            cols = a.get("uses_columns", {}).get(clin, [])
            if cols:
                return cols[0]
    return None


def tab_3d(doc, a_dir, b_dir):
    f = main_file(doc)
    if not f:
        st.info("No sample-level changes to show.")
        return
    ta, tb = get_table(a_dir, f), get_table(b_dir, f)
    if ta is None or tb is None:
        st.info(f"{f} is not in both releases.")
        return
    if numeric(ta).notna().mean().mean() < 0.5:
        st.info(f"{f} is not a number table, so it cannot be drawn in 3D.")
        return
    with st.spinner("Building the 3D map..."):
        pa, pb = pca3(ta, tb)
    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    clin = get_table(a_dir, clin_name) if clin_name else None
    col = color_column(doc, f, clin_name) if clin is not None else None

    fig = go.Figure()
    if col and col in clin.columns:
        lab = clin[col].reindex(pb.index).fillna("unknown").astype(str)
        for val_ in lab.value_counts().index[:6]:
            idx = lab.index[lab == val_]
            fig.add_trace(go.Scatter3d(x=pb.loc[idx, "x"], y=pb.loc[idx, "y"], z=pb.loc[idx, "z"],
                                       mode="markers", name=f"{col}: {val_}",
                                       marker=dict(size=2.5, opacity=0.35)))
    else:
        fig.add_trace(go.Scatter3d(x=pb["x"], y=pb["y"], z=pb["z"], mode="markers", name="samples",
                                   marker=dict(size=2.5, opacity=0.35, color="#adb5bd")))

    ev = evidence_by_file(doc).get(f, {})
    for kind, ids in ev.items():
        color = KIND_COLOR.get(kind, "#000000")
        in_b = [s for s in ids if s in pb.index]
        in_a_only = [s for s in ids if s in pa.index and s not in pb.index]
        if in_b:
            fig.add_trace(go.Scatter3d(
                x=pb.loc[in_b, "x"], y=pb.loc[in_b, "y"], z=pb.loc[in_b, "z"],
                mode="markers+text", text=[s[5:12] for s in in_b], textposition="top center",
                name=f"{KIND_LABEL.get(kind, kind)} (new position)",
                marker=dict(size=7, color=color, symbol="diamond", line=dict(width=1, color="black"))))
        moved = [s for s in in_b if s in pa.index]
        for s in moved:
            fig.add_trace(go.Scatter3d(
                x=[pa.at[s, "x"], pb.at[s, "x"]], y=[pa.at[s, "y"], pb.at[s, "y"]], z=[pa.at[s, "z"], pb.at[s, "z"]],
                mode="lines", line=dict(color=color, width=5, dash="dash"), showlegend=False,
                hovertext=f"{s}: old to new position"))
        if moved:
            fig.add_trace(go.Scatter3d(
                x=pa.loc[moved, "x"], y=pa.loc[moved, "y"], z=pa.loc[moved, "z"], mode="markers",
                name=f"{KIND_LABEL.get(kind, kind)} (old position)",
                marker=dict(size=6, color="white", line=dict(width=2, color=color))))
        if in_a_only:
            fig.add_trace(go.Scatter3d(
                x=pa.loc[in_a_only, "x"], y=pa.loc[in_a_only, "y"], z=pa.loc[in_a_only, "z"],
                mode="markers", name=f"{KIND_LABEL.get(kind, kind)} (gone from B)",
                marker=dict(size=6, color=color, symbol="x")))
    fig.update_layout(height=650, margin=dict(l=0, r=0, t=30, b=0),
                      legend=dict(itemsizing="constant"),
                      scene=dict(xaxis_title="PC1", yaxis_title="PC2", zaxis_title="PC3"),
                      title=f"Every sample of {f} as one dot (rotate with the mouse)")
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Each dot is one tumor. Its values in this file are reduced to 3 numbers (principal components), "
               "so tumors with similar values sit close together. Big markers are samples with findings. "
               "Dashed lines go from a sample's old position (release A) to its new position (release B). "
               "Copies sit exactly on top of the sample they were copied from.")


def tab_analysis(doc, a_dir, b_dir):
    fa, _ = release_files(doc)
    clin_name = clinical_name(fa)
    if not clin_name:
        st.info("No clinical file found.")
        return
    ev = evidence_by_file(doc)
    shown = False
    for a in doc.get("analyses", []):
        uc = a.get("uses_columns", {})
        group_cols = uc.get(clin_name, [])
        others = [(fn, feats) for fn, feats in uc.items() if fn != clin_name and feats]
        if len(group_cols) != 1 or not others:
            continue
        fn, feats = others[0]
        feat, gcol = feats[0], group_cols[0]
        clin_a, clin_b = get_table(a_dir, clin_name), get_table(b_dir, clin_name)
        ta, tb = get_table(a_dir, fn), get_table(b_dir, fn)
        if any(x is None for x in (clin_a, clin_b, ta, tb)) or feat not in ta.columns:
            continue
        hits = set().union(*ev.get(fn, {}).values()) if ev.get(fn) else set()
        fig = go.Figure()
        for rel_name, clin, t in (("old release", clin_a, ta), ("new release", clin_b, tb)):
            if feat not in t.columns:
                continue
            ids = [s for s in t.index if str(s)[13:15] == "01"]
            g = clin[gcol].reindex(ids)
            v = pd.to_numeric(t.loc[ids, feat], errors="coerce")
            keep = g.value_counts()
            for grp in keep.index[keep >= 10][:3]:
                sel = g.index[g == grp]
                name = f"{rel_name}<br>{grp}"
                fig.add_trace(go.Box(y=v.loc[sel], name=name, boxpoints="all", jitter=0.5, pointpos=0,
                                     marker=dict(size=3, opacity=0.4), line=dict(width=1), showlegend=False,
                                     boxmean=True))
                hs = [s for s in sel if s in hits]
                if hs:
                    fig.add_trace(go.Scatter(x=[name] * len(hs), y=v.loc[hs], mode="markers+text",
                                             text=[s[5:12] for s in hs], textposition="middle right",
                                             marker=dict(size=13, color="#e63946", symbol="diamond",
                                                         line=dict(width=1, color="black")),
                                             showlegend=False))
        fig.update_layout(height=520, margin=dict(l=10, r=10, t=40, b=10),
                          title=f"{a['id']}: {feat} by {gcol} (dashed line = mean, red = samples with findings)",
                          yaxis_title=feat)
        st.plotly_chart(fig, use_container_width=True)
        before, after = a.get("before", {}), a.get("after", {})
        if "difference" in before and "difference" in after:
            st.caption(f"{a['id']} difference between the groups: {before['difference']} in the old release, "
                       f"{after['difference']} in the new release, "
                       f"{a.get('after_repair', {}).get('difference', '?')} after repair.")
        shown = True
    if not shown:
        st.info("No analysis compares a number across clinical groups in this run.")


def tab_swap(doc, a_dir, b_dir):
    swaps = [f for f in doc.get("findings", []) if f.get("kind") == "sample_swap"]
    if not swaps:
        st.info("No swaps in this release.")
        return
    f = swaps[0]
    samples = [e for e in f["evidence"] if e.get("kind") == "sample"]
    if not samples:
        st.info("Swap has no sample evidence.")
        return
    file = samples[0]["file"]
    s1 = samples[0]["id"]
    note = samples[0].get("note", "")
    s2 = note.split("profile of ")[-1].strip() if "profile of " in note else (samples[1]["id"] if len(samples) > 1 else None)
    ta, tb = get_table(a_dir, file), get_table(b_dir, file)
    if ta is None or tb is None or s2 is None or s1 not in tb.index or s2 not in ta.index or s1 not in ta.index:
        st.info("Cannot draw this swap.")
        return
    cols = [c for c in ta.columns if c in tb.columns]
    b1 = pd.to_numeric(tb.loc[s1, cols], errors="coerce")
    a2 = pd.to_numeric(ta.loc[s2, cols], errors="coerce")
    a1 = pd.to_numeric(ta.loc[s1, cols], errors="coerce")
    ok = b1.notna() & a2.notna() & a1.notna()
    b1, a2, a1 = b1[ok], a2[ok], a1[ok]
    lo, hi = float(min(b1.min(), a2.min(), a1.min())), float(max(b1.max(), a2.max(), a1.max()))
    c1, c2 = st.columns(2)
    for col_, x, title in ((c1, a2, f"new {s1} vs old {s2}"), (c2, a1, f"new {s1} vs old {s1}")):
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(dash="dash", color="#adb5bd"),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=x, y=b1, mode="markers", text=list(b1.index), showlegend=False,
                                 marker=dict(size=7, color="#e63946")))
        diff = float((b1 - x).abs().max())
        fig.update_layout(height=420, margin=dict(l=10, r=10, t=50, b=10), title=title,
                          xaxis_title="old release value", yaxis_title=f"new release value of {s1}")
        col_.plotly_chart(fig, use_container_width=True)
        col_.caption(f"{len(b1)} values compared. Largest difference: {diff:.12g}. "
                     + ("Identical: every dot is on the line." if diff == 0 else "Not the same profile."))
    st.caption(f"Left: the new values of {s1} are exactly the old values of {s2}. "
               f"Right: they do not match {s1}'s own old values. So the data was moved, not measured again.")


# ---------- page ----------

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
                (st.success if r.returncode == 0 else st.error)(
                    "Envelope opened." if r.returncode == 0 else "Unseal failed.")
                st.code((r.stdout + r.stderr)[-3000:])
    with tab_up:
        up_name = st.text_input("Name for this release", placeholder="my-release")
        files = st.file_uploader("Release files (.tsv, .tsv.gz, .csv)", accept_multiple_files=True)
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
        st.session_state.pop("dash", None)

if "out" not in st.session_state:
    st.info("Pick a release in the sidebar and press Run check.")
    st.stop()

out = Path(st.session_state["out"])
doc = json.loads(out.read_text(encoding="utf-8"))
sns_run = doc.get("sns_run", {})
a_dir, b_dir = st.session_state["a"], st.session_state["b"]

verdict = doc.get("verdict", "Unknown")
VERDICT_STYLE.get(verdict, st.info)(f"### Verdict: {verdict}")
st.write(doc.get("question", ""))
st.caption(f"Release A: {a_dir}  |  Release B: {b_dir}  |  "
           f"{st.session_state['seconds']} s  |  saved as {out.relative_to(ROOT)}")

if st.button("Open React dashboard", type="primary"):
    dash = out.with_suffix(".html")
    r = run_cmd([PY, "sns/dashboard.py", "--findings", str(out), "--a", a_dir, "--b", b_dir, "--out", str(dash)])
    if r.returncode == 0 and dash.exists():
        st.session_state["dash"] = str(dash)
    else:
        st.error("Dashboard build failed.")
        st.code((r.stdout + r.stderr)[-3000:])
if st.session_state.get("dash") and Path(st.session_state["dash"]).exists():
    import streamlit.components.v1 as components
    dash_path = Path(st.session_state["dash"])
    st.caption(f"React dashboard file: {dash_path.relative_to(ROOT)} (you can also open it directly in a browser)")
    components.html(dash_path.read_text(encoding="utf-8"), height=1500, scrolling=True)
    st.download_button("Download React dashboard", dash_path.read_bytes(), file_name=dash_path.name, mime="text/html")

st.header("Visual story")
tabs = st.tabs(["Pipeline", "Sample flow", "3D map", "Analysis view", "Swap proof"])
for tab, fn in zip(tabs, (tab_pipeline, tab_flow, tab_3d, tab_analysis, tab_swap)):
    with tab:
        try:
            fn(doc, a_dir, b_dir)
        except Exception as exc:
            st.info(f"This view is not available for this release ({type(exc).__name__}).")

st.header("Analyses")
for a in doc.get("analyses", []):
    st.subheader(f"{a['id']}: {a.get('description', '')}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Status", a.get("status", ""))
    c2.metric("Verdict", a.get("verdict", ""))
    c3.metric("Repair restores original", "yes" if a.get("repair_restores_before") else "no")
    st.write(a.get("reason", ""))
    fb_, fa_, fr_ = flat(a.get("before")), flat(a.get("after")), flat(a.get("after_repair"))
    keys = sorted(set(fb_) | set(fa_) | set(fr_))
    if keys:
        st.dataframe(pd.DataFrame({
            "value": keys,
            "old release": [fb_.get(k) for k in keys],
            "new release": [fa_.get(k) for k in keys],
            "after repair": [fr_.get(k) for k in keys],
        }), hide_index=True, use_container_width=True)
    if a.get("affected_by"):
        st.caption("Affected by: " + ", ".join(a["affected_by"]))

st.header(f"Findings ({len(doc.get('findings', []))})")
if not doc.get("findings"):
    st.success("No changes found that need review.")
for f in doc.get("findings", []):
    with st.expander(f"{f['id']} [{f['severity']}] {f.get('kind', '')}: {f['summary']}", expanded=True):
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
        "id": r["id"], "action": r["action"], "fixes": ", ".join(r.get("addresses", [])),
        "targets": len(r.get("targets", [])), "justification": r.get("justification", ""),
    } for r in reps]), hide_index=True, use_container_width=True)
    st.caption("Repairs are proposals. A person must approve or reject each one.")
else:
    st.write("No repairs proposed.")

missing = sns_run.get("files_not_in_release") or []
if missing:
    st.info("Files not in this release (no analysis reads them, so they were not compared): " + ", ".join(missing))

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
st.download_button("Download findings.json", out.read_bytes(),
                   file_name=f"findings_{st.session_state['label']}.json", mime="application/json")
analysis = st.selectbox("Clean cohort for analysis", [a["id"] for a in doc.get("analyses", [])] or ["A1"])
if st.button("Build clean cohort"):
    cdir = RUNS / f"cohort_{st.session_state['label']}_{analysis}_{datetime.now().strftime('%H%M%S')}"
    r = run_cmd([PY, "sns/cohort.py", "--release", b_dir, "--prev", a_dir,
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