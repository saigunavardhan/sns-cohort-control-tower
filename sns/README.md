# SNS Cohort Control Tower

**Product:** Cohort Control Tower
**Track:** Cohort Firewall
**Team:** SNS
**Team members:** <name 1>, <name 2>, <name 3>

---

## What it does and who uses it

A research engineer is responsible for a TCGA-BRCA dataset that gets new releases. When a new release arrives, they need to know: did any samples get swapped, duplicated, or dropped? Do those changes affect their analyses? Can the analyses be safely rerun?

The workflow is:

1. Connect two release folders (the previous release A and the new release B).
2. The tool detects what changed: sample swaps, duplicates, dropped samples, missing values, relabeled IDs.
3. It reruns the analyses on both releases and on a repaired copy of B.
4. It writes `findings.json` describing every change, which analyses are affected, and whether a repair restores the original results.
5. The engineer reviews the findings and decides whether to proceed.

---

## Setup

Clone the repo, then install dependencies:

```
conda create -n sns python=3.11
conda activate sns
pip install -r requirements.txt
pip install pyarrow jsonschema cryptography
python data/tcga-brca/fetch.py
```

---

## Running the tool

### On practice release B (default)

```
python sns/run_tower.py --out sns/findings_practiceB.json
```

With human-readable summary:

```
python sns/run_tower.py --out sns/findings_practiceB.json --pretty
```

### On an envelope (live demo format)

First unseal the envelope:

```
python envelopes/unseal.py <number> <four words>
```

Example (practice envelope):

```
python envelopes/unseal.py practice orange meadow pilot lantern
```

Then run the tool:

```
python sns/run_tower.py --b envelopes/opened/practice/release --out sns/findings_envelope.json --pretty
```

### Validate the output

```
python format/validate.py sns/findings_practiceB.json
```

### Clean Cohort Export

Build a clean sample list from the findings, with repairs applied:

```
python sns/cohort.py \
  --release data/tcga-brca/releases/practice-B \
  --findings sns/findings_practiceB.json \
  --analysis A1 \
  --out-dir sns/cohort_out
```

Outputs: `cohort_samples.csv`, `cohort_steps.json`, `reproduce_cohort.py`.

---

## What the output means

### Verdict levels

- **Ready**: no finding touches this analysis's input files or samples; rerun is safe.
- **Conditional**: the analysis is affected by a finding, but results moved less than 10%. Rerun after fixing the data errors.
- **Insufficient Evidence**: too few samples remain after the change to draw conclusions.
- **Blocked**: a numeric result changed sign, or shifted by 10% or more. The release must not be used until the data error is resolved.

### Finding kinds

- **sample_swap**: two or more samples appear to have had their data exchanged.
- **duplicate**: a sample ID carries a profile identical to an existing sample.
- **membership_change**: samples were added, dropped, or renamed.
- **missing_data**: values or entire files are absent in the new release.
- **label_conflict**: non-numeric values changed in a way that does not match any other sample.
- **other**: numeric values changed but no swap or duplicate pattern was found.

### Repair proof

For each analysis, `repair_restores_before: true` means that after applying the listed repairs to release B (excluding duplicates, restoring swapped profiles, adding back dropped samples), the analysis produces the same result as on release A. This confirms the finding is real and the repair is correct.

---

## A result a judge can verify

**Practice release B contains a sample swap in `protein_rppa.tsv.gz`.**

TCGA-A8-A09K-01 and TCGA-B6-A0WX-01 have had their ERALPHA protein profiles exchanged:

| Sample | ERALPHA in raw (A) | ERALPHA in practice-B |
|---|---|---|
| TCGA-A8-A09K-01 | 3.43266169795726 | -1.07005186202866 |
| TCGA-B6-A0WX-01 | -1.07005186202866 | 3.43266169795726 |

Verify these values directly:

```python
import pandas as pd
raw = pd.read_csv("data/tcga-brca/raw/protein_rppa.tsv.gz", sep="\t", index_col=0)
print(raw.at["ERALPHA", "TCGA-A8-A09K-01"])   # 3.43266169795726
print(raw.at["ERALPHA", "TCGA-B6-A0WX-01"])   # -1.07005186202866

b = pd.read_csv("data/tcga-brca/releases/practice-B/protein_rppa.tsv.gz", sep="\t", index_col=0)
print(b.at["ERALPHA", "TCGA-A8-A09K-01"])     # -1.07005186202866
print(b.at["ERALPHA", "TCGA-B6-A0WX-01"])     # 3.43266169795726
```

Effect on analysis A1 (mean ER-alpha by ER status):

- mean_negative: -0.872 in raw, -0.843 in practice-B
- difference (Positive minus Negative): 3.051 in raw, 3.019 in practice-B
- After repair: analysis restores exactly to the raw values.

---

## Limitations

- Tabular data only: image files, PDFs, and binary formats are not checked.
- Sample IDs must start with `TCGA-`. IDs in another format are ignored.
- Near copies in sparse files (like mutation calls, where most values are 0) cannot be reliably told apart from real new samples. This is a known miss and is stated in the finding limitation.
- The tool sees what changed between two releases, not why. Human review is required to decide whether a change is an error or a valid update.
- Research use only. No clinical claims should be made from this output.

---

## How we tested it

The chaos tester (`sns/chaos.py`) generates random releases with up to 4 injected changes (swaps, duplicates, drops, edits, near copies, feature drops, subgroup drops, and more), runs the full detection engine, and checks whether every injected change is caught.

Latest results (n=10 releases, seed=21):

- Changes caught: 21 of 21 (100%)
- Crashes: 0
- Invalid output files: 0
- False findings: 0

The two known miss types (near copies in sparse files) are excluded from scoring because they cannot be detected without non-sparse reference data. This is documented as a limitation.

To rerun:

```
python sns/chaos.py --n 10 --seed 21
```

---

## Demo video

<link>

---

## Third-party libraries and licenses

- pandas (BSD-3)
- numpy (BSD-3)
- pyarrow (Apache-2.0)
- jsonschema (MIT)
- cryptography (Apache-2.0 or BSD)

## 4. Web app (recommended): the order to test it

The web app has a FastAPI backend (`sns/server.py`) and a React page (`sns/web/index.html`). It uses the same engine as the terminal version, so the results are identical.

**Step 1. Start the server**
```
python sns/server.py
```
It prints `{"url": "http://127.0.0.1:8000"}`. Open http://127.0.0.1:8000 in a browser and keep the terminal open while using the page. The page loads React, Plotly and its fonts from the internet, so the laptop needs a connection.

**Step 2. Open the practice envelope and run the check**
In the left panel, tab **Envelope**: type `practice` as the number and `orange meadow pilot lantern` as the four words, then click **Open envelope and run check**. The check takes a few seconds.

**Step 3. Read the results, top to bottom**
- **Top box:** the verdict, the two releases with their sample counts, and the changed samples by type. Practice: Conditional; same count (747 in both releases), but 12 samples changed in `protein_rppa.tsv.gz`: 2 profiles swapped, 5 copied profiles under new IDs, 5 samples lost.
- **Impact map:** files on the left, findings in the middle, analyses on the right. Colored lines show which file each finding comes from and which analysis it reaches; dotted lines show the files each analysis reads. Practice: A1 affected, A2 left alone. Click a finding to jump to its evidence.
- **A1 replay:** each dot is one primary tumor (ERALPHA value by ER status). Click **Play the story** to watch the same samples move from Release A to Release B, then to B after the proposed repairs. The three cards show the result at each stage: 3.051, 3.019, 3.051.
- **A2 bars:** subtype shares in Release A, Release B and after repair. Identical bars mean the result is safe.
- **Analysis tables:** every value of A1 and A2 on the old release, the new release and after repair.
- **Proof:** for a swap, the new values of one sample plotted against the old values of the other: every point on the line, largest difference 0.
- **How it was checked:** the six steps with their counts.
- **Evidence:** one card per finding with every sample involved.
- **Repairs and limits:** the proposed repairs (proposals only), files not in the release, and what the check cannot tell you.
- **Export:** download `findings.json`, and build the clean cohort for A1 or A2 (A1 on practice: 698 samples, 0 duplicates) as a CSV.

**Step 4. A judge's envelope**
Same as step 2, with the judge's number and four words.

**Step 5. Other releases**
- Tab **Saved:** every release on the laptop (opened envelopes, `data/tcga-brca/releases/*`, uploads). Choose one, then click **Run check**.
- Tab **Upload:** give a name, choose the release files, click **Save release**.
- **Recent checks:** reopen any earlier check without rerunning it.

Every check is saved in `sns/runs/<run id>/` (`findings.json`, the chart data and the cohort files), so nothing is overwritten. The **Dark / Light** button switches the color theme.

**Backend API** (used by the page, also usable directly):
- `GET /api/releases`: list releases
- `POST /api/unseal`: open an envelope (`number`, `words`)
- `POST /api/upload`: save uploaded release files
- `POST /api/run`: run a check (`release_id`, optional `a_path`)
- `GET /api/runs` and `GET /api/runs/{run_id}`: list checks, open one
- `GET /api/runs/{run_id}/findings.json`: download the findings file
- `POST /api/runs/{run_id}/cohort` and `GET /api/runs/{run_id}/cohort/{analysis}.csv`: build and download a clean cohort

### Streamlit version (alternative)
```
streamlit run sns/app.py
```
In the sidebar: tab **Sealed envelope**, type the number and the four words, click **Unseal**; choose the release in **New release (B)**; click **Run check**. The page shows the verdict, the A1 replay, the analyses, the findings, the repairs, the limitations and the downloads.
