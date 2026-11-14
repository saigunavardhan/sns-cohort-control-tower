import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from pydantic import BaseModel

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from engine.viz import build_viz  # noqa: E402

RAW = ROOT / "data" / "tcga-brca" / "raw"
RUNS = HERE / "runs"
UPLOADS = HERE / "uploads"
WEB = HERE / "web"
PY = sys.executable
ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}
SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]{1,60}$")

app = FastAPI(title="SNS Cohort Control Tower", version="1.0")


def run_cmd(cmd, timeout=600):
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
                          errors="replace", env=ENV, timeout=timeout)


def releases():
    found = []
    opened = ROOT / "envelopes" / "opened"
    if opened.exists():
        for p in sorted(opened.glob("*/release")):
            found.append({"id": f"envelope:{p.parent.name}", "label": f"Envelope {p.parent.name}",
                          "kind": "envelope", "path": str(p)})
    rel = ROOT / "data" / "tcga-brca" / "releases"
    if rel.exists():
        for p in sorted(rel.iterdir()):
            if p.is_dir():
                found.append({"id": f"release:{p.name}", "label": p.name, "kind": "release", "path": str(p)})
    if UPLOADS.exists():
        for p in sorted(UPLOADS.iterdir()):
            if p.is_dir():
                found.append({"id": f"upload:{p.name}", "label": p.name, "kind": "upload", "path": str(p)})
    return found


def release_path(rid):
    for r in releases():
        if r["id"] == rid:
            return Path(r["path"])
    raise HTTPException(404, f"Release {rid} not found. Refresh the list.")


def run_dir(run_id):
    if not SAFE_NAME.match(run_id):
        raise HTTPException(400, "Bad run id.")
    d = RUNS / run_id
    if not (d / "findings.json").exists():
        raise HTTPException(404, f"Run {run_id} not found.")
    return d


@app.get("/api/health")
def health():
    return {"ok": True, "release_a": RAW.exists()}


@app.get("/api/releases")
def list_releases():
    return {"releases": releases(), "default_a": str(RAW)}


class UnsealBody(BaseModel):
    number: str
    words: str


@app.post("/api/unseal")
def unseal(body: UnsealBody):
    words = body.words.split()
    number = body.number.strip()
    if not number or len(words) != 4 or not SAFE_NAME.match(number):
        raise HTTPException(400, "Give the envelope number and exactly four words.")
    r = run_cmd([PY, "envelopes/unseal.py", number, *words], timeout=300)
    log = (r.stdout + r.stderr)[-4000:]
    if r.returncode != 0:
        return JSONResponse({"ok": False, "log": log}, status_code=400)
    name = "practice" if number.lower() == "practice" else f"E{number.zfill(2)}" if number.isdigit() else number
    rid = f"envelope:{name}"
    known = [x["id"] for x in releases()]
    return {"ok": True, "log": log, "release_id": rid if rid in known else None}


@app.post("/api/upload")
async def upload(name: str = Form(...), files: List[UploadFile] = File(...)):
    name = name.strip()
    if not SAFE_NAME.match(name):
        raise HTTPException(400, "Use letters, numbers, dots, dashes or underscores for the name.")
    target = UPLOADS / name
    target.mkdir(parents=True, exist_ok=True)
    saved = []
    for f in files:
        fname = Path(f.filename).name
        if not fname.lower().endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz")):
            continue
        (target / fname).write_bytes(await f.read())
        saved.append(fname)
    if not saved:
        raise HTTPException(400, "No .tsv, .tsv.gz, .csv or .csv.gz files were given.")
    return {"ok": True, "release_id": f"upload:{name}", "files": saved}


class RunBody(BaseModel):
    release_id: str
    a_path: Optional[str] = None


@app.post("/api/run")
def run(body: RunBody):
    b = release_path(body.release_id)
    a = Path(body.a_path) if body.a_path else RAW
    if not a.exists():
        raise HTTPException(400, f"Release A folder not found: {a}")
    label = re.sub(r"[^A-Za-z0-9_-]", "_", body.release_id.split(":", 1)[-1])
    run_id = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}-{label}"
    d = RUNS / run_id
    d.mkdir(parents=True, exist_ok=True)
    out = d / "findings.json"
    t0 = time.time()
    r = run_cmd([PY, "sns/run_tower.py", "--a", str(a), "--b", str(b), "--out", str(out)])
    seconds = round(time.time() - t0, 1)
    if r.returncode != 0 or not out.exists():
        return JSONResponse({"ok": False, "log": (r.stdout + r.stderr)[-4000:]}, status_code=500)
    doc = json.loads(out.read_text(encoding="utf-8"))
    viz = build_viz(doc, str(a), str(b))
    meta = {"run_id": run_id, "release_id": body.release_id, "a": str(a), "b": str(b),
            "seconds": seconds, "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
    (d / "viz.json").write_text(json.dumps(viz, default=str), encoding="utf-8")
    (d / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    return {"ok": True, "meta": meta, "doc": doc, "viz": viz}


@app.get("/api/runs")
def list_runs():
    out = []
    if RUNS.exists():
        for d in sorted(RUNS.iterdir(), reverse=True):
            m = d / "meta.json"
            f = d / "findings.json"
            if m.exists() and f.exists():
                meta = json.loads(m.read_text(encoding="utf-8"))
                doc = json.loads(f.read_text(encoding="utf-8"))
                meta["verdict"] = doc.get("verdict")
                meta["n_findings"] = len(doc.get("findings", []))
                out.append(meta)
    return {"runs": out[:20]}


@app.get("/api/runs/{run_id}")
def get_run(run_id: str):
    d = run_dir(run_id)
    doc = json.loads((d / "findings.json").read_text(encoding="utf-8"))
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8")) if (d / "meta.json").exists() else {"run_id": run_id}
    if (d / "viz.json").exists():
        viz = json.loads((d / "viz.json").read_text(encoding="utf-8"))
    else:
        viz = build_viz(doc, meta.get("a", str(RAW)), meta.get("b", ""))
    return {"ok": True, "meta": meta, "doc": doc, "viz": viz}


@app.get("/api/runs/{run_id}/findings.json")
def download_findings(run_id: str):
    d = run_dir(run_id)
    return FileResponse(d / "findings.json", media_type="application/json", filename=f"findings_{run_id}.json")


class CohortBody(BaseModel):
    analysis: str


@app.post("/api/runs/{run_id}/cohort")
def cohort(run_id: str, body: CohortBody):
    d = run_dir(run_id)
    if not SAFE_NAME.match(body.analysis):
        raise HTTPException(400, "Bad analysis id.")
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    cdir = d / f"cohort_{body.analysis}"
    r = run_cmd([PY, "sns/cohort.py", "--release", meta["b"], "--prev", meta["a"],
                 "--findings", str(d / "findings.json"), "--analysis", body.analysis, "--out-dir", str(cdir)])
    csv = cdir / "cohort_samples.csv"
    if r.returncode != 0 or not csv.exists():
        return JSONResponse({"ok": False, "log": (r.stdout + r.stderr)[-3000:]}, status_code=500)
    import pandas as pd
    c = pd.read_csv(csv)
    flagged = int(c["flagged"].astype(str).str.lower().isin(["true", "1"]).sum()) if "flagged" in c.columns else 0
    return {"ok": True, "analysis": body.analysis, "n_samples": int(len(c)),
            "n_patients": int(c["patient_id"].nunique()) if "patient_id" in c.columns else None,
            "n_duplicates": int(c["sample_id"].duplicated().sum()), "n_flagged": flagged,
            "download": f"/api/runs/{run_id}/cohort/{body.analysis}.csv"}


@app.get("/api/runs/{run_id}/cohort/{analysis}.csv")
def cohort_csv(run_id: str, analysis: str):
    d = run_dir(run_id)
    if not SAFE_NAME.match(analysis):
        raise HTTPException(400, "Bad analysis id.")
    csv = d / f"cohort_{analysis}" / "cohort_samples.csv"
    if not csv.exists():
        raise HTTPException(404, "Build the cohort first.")
    return FileResponse(csv, media_type="text/csv", filename=f"cohort_samples_{analysis}_{run_id}.csv")


@app.get("/", response_class=HTMLResponse)
def index():
    return HTMLResponse((WEB / "index.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("SNS_PORT", "8000"))
    print(json.dumps({"url": f"http://127.0.0.1:{port}"}))
    uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")