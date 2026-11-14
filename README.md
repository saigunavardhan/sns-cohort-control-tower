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
