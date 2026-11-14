"""Animated A1 replay (release A -> release B -> B repaired) for the Streamlit app.

Usage in app.py:
    from swap_replay import build_html
    html = build_html(st.session_state["a"], st.session_state["b"], doc)
    if html:
        components.html(html, height=900, scrolling=True)
"""
import json
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "swap_replay_template.html"
ER = "breast_carcinoma_estrogen_receptor_status"
GROUPS = ("Positive", "Negative")


def _a1_points(release):
    """ER label and ERALPHA for every sample A1 uses (same rule as reference_analyses.a1)."""
    r = Path(release)
    if not (r / "clinical.tsv").exists() or not (r / "protein_rppa.tsv.gz").exists():
        return None
    clinical = pd.read_csv(r / "clinical.tsv", sep="\t", low_memory=False)
    protein = pd.read_csv(r / "protein_rppa.tsv.gz", sep="\t", index_col=0)
    if "ERALPHA" not in protein.index or ER not in clinical.columns:
        return None
    er = clinical.drop_duplicates("sampleID").set_index("sampleID")[ER]
    row = protein.loc["ERALPHA"]
    return {s: (er[s], round(float(row[s]), 3)) for s in protein.columns
            if str(s)[13:15] == "01" and er.get(s) in GROUPS and pd.notna(row[s])}


def _statuses(doc):
    """Map sample id -> (status, source) from the findings."""
    st = {}
    for f in doc.get("findings", []):
        kind = f.get("kind")
        for e in f.get("evidence", []):
            if e.get("kind") != "sample":
                continue
            sid, note = e.get("id"), e.get("note") or ""
            if kind == "sample_swap":
                st[sid] = ("swap", None)
            elif kind == "duplicate":
                m = re.search(r"identical to (TCGA-[\w-]+)", note)
                st.setdefault(sid, ("copy", m.group(1) if m else None))
            elif kind == "label_conflict" and ER in note:
                st.setdefault(sid, ("relabel", None))
            elif kind == "membership_change" and "absent" in note:
                st.setdefault(sid, ("removed", None))
    return st


def _a1_numbers(doc, key):
    a1 = next((a for a in doc.get("analyses", []) if a.get("id") == "A1"), None)
    return (a1 or {}).get(key)


def build_html(release_a, release_b, doc):
    """Return the page as a string, or None if A1 cannot be drawn for these releases."""
    pa, pb = _a1_points(release_a), _a1_points(release_b)
    before, after = _a1_numbers(doc, "before"), _a1_numbers(doc, "after")
    if not pa or not pb or not before or not after:
        return None
    repaired = _a1_numbers(doc, "after_repair") or before
    status = _statuses(doc)
    points = []
    for s in sorted(set(pa) | set(pb)):
        st, src = status.get(s, (None, None))
        if st is None:
            if s in pa and s not in pb:
                st = "removed"
            elif s in pa and s in pb and pa[s][0] != pb[s][0]:
                st = "relabel"
            elif s in pa and s in pb and pa[s][1] != pb[s][1]:
                st = "swap"
            else:
                st = "same"
        points.append({"id": s, "st": st, "src": src,
                       "erA": pa[s][0] if s in pa else None, "a": pa[s][1] if s in pa else None,
                       "erB": pb[s][0] if s in pb else None, "b": pb[s][1] if s in pb else None})
    data = {"points": points, "A": before, "B": after, "R": repaired,
            "restores": bool(_a1_numbers(doc, "repair_restores_before"))}
    html = TEMPLATE.read_text(encoding="utf-8")
    return html.replace("__DATA__", json.dumps(data).replace("</", "<\\/"))
