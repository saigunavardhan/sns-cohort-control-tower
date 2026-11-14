# SNS Cohort Control Tower

The organizers' original challenge description is in `README_CHALLENGE.md`.

**Product:** Cohort Control Tower
**Track:** Cohort Firewall
**Team:** SNS
**Team members:** Sai Guna Vardhan Vudatala, Nishant Naravarajula, Shravya Ushake
**Demo video:** [Watch the demo (Google Drive)](https://drive.google.com/drive/folders/1eTEWI1hvT8cqYyvnRcI8ZyK73mfqb9SQ?usp=sharing)

Research use only. No diagnosis, no treatment advice, no clinical claims.

---

## 1. What it does and who uses it

**User:** a research engineer who looks after the TCGA breast cancer dataset, which gets new releases.

**Problem:** a new release can silently break research. Sample counts can stay the same while samples are swapped, copied, lost or edited inside the files. A researcher who reruns their analysis on the new release gets a different result and never knows why.

**Workflow (one run):**

1. Give the tool two release folders: the trusted old release (A) and the new release (B).
2. It makes a sha256 fingerprint of every file, to see which files changed. Unchanged files are skipped.
3. It makes a sha256 fingerprint of every sample's row, to see which samples changed, even when counts stay the same.
4. It matches the fingerprints between A and B and names each change: swap, duplicate (copy), lost sample, renamed ID, emptied values, edited values, near copy, removed column, missing file.
5. It checks if the lost samples all belong to one group (for example one PAM50 subtype).
6. It reads the code of the two analyses (A1 and A2) to see which files and columns each one uses, so it knows which analysis each change affects.
7. It reruns A1 and A2 three times: on A, on B, and on B after its own proposed repairs.
8. It gives a verdict for each analysis and overall, and writes everything to `findings.json` in the shared format v0.2.

**The two analyses (given by the organizers in `data/tcga-brca/reference_analyses.py`):**

- **A1:** mean ERALPHA protein in ER positive tumors minus mean ERALPHA in ER negative tumors. Reads `clinical.tsv` and `protein_rppa.tsv.gz`.
- **A2:** percent of tumors in each PAM50 subtype. Reads `clinical.tsv` only.

---

## 2. Setup (do once)

Windows, macOS and Linux all work. The commands below use forward slashes, which work in all of them.

```
git clone https://github.com/saigunavardhan/sns-cohort-control-tower.git
cd sns-cohort-control-tower
conda create -n sns python=3.11 -y
conda activate sns
pip install -r requirements.txt -r sns/requirements.txt
python data/tcga-brca/fetch.py
python data/tcga-brca/make_release.py
python data/tcga-brca/make_release.py --verify
```

- `fetch.py` downloads release A into `data/tcga-brca/raw/`.
- `make_release.py` builds practice release B into `data/tcga-brca/releases/practice-B/`. The `--verify` line must print `practice-B matches the published release`.

No conda? Use `python -m venv .venv`, activate it, then run the same `pip install` line.

---

## 3. Terminal version: the order to test it

Run these from the repo folder, with the environment active.

**Step 1. Practice release B**
```
python sns/run_tower.py --out sns/findings_practiceB.json --pretty
python format/validate.py sns/findings_practiceB.json
```
Expected: verdict `Conditional`; F1 sample_swap (high), F2 duplicate (high), F3 membership_change (medium, 5 Her2 samples, 61 in raw and 56 in practice-B); A1 changed, A2 unchanged; repair restores original results for A1 and A2; the validator prints `valid`. This matches the organizers' reference answer in `format/examples/control-tower-practice-B.json`.

**Step 2. The practice envelope (exactly like the live judging)**
```
python envelopes/unseal.py practice orange meadow pilot lantern
python sns/run_tower.py --b envelopes/opened/practice/release --out sns/findings_envelope_practice.json --pretty
python format/validate.py sns/findings_envelope_practice.json
```
Expected: the same 3 findings, no finding about RNA, and `files not in release: rna_hiseq.tsv.gz` (envelopes leave the RNA file out on purpose, and no analysis reads it). Runs in about 2 seconds.

**Step 3. A judge's envelope**
The judge gives an envelope number and four words. Example for envelope 07:
```
python envelopes/unseal.py 07 amber falcon river quartz
python sns/run_tower.py --b envelopes/opened/E07/release --out sns/findings_E07.json --pretty
python format/validate.py sns/findings_E07.json
```
Use a different `--out` name for every release, otherwise the previous file is overwritten.

**Step 4. Clean cohort export (uses the findings from step 1)**
```
python sns/cohort.py --release data/tcga-brca/releases/practice-B --findings sns/findings_practiceB.json --analysis A1 --out-dir sns/cohort_out
```
Expected: `sns/cohort_out/cohort_samples.csv` with 698 samples and no duplicates, the same samples A1 uses on release A. Also written: `cohort_steps.json` (counts after every step) and `reproduce_cohort.py` (a script that rebuilds the same list).

**Step 5. Test on a release the tool has never seen, with a hidden answer key**
```
python sns/make_test_release.py --seed 101 --out data/tcga-brca/releases/test-1
python sns/run_tower.py --b data/tcga-brca/releases/test-1 --out sns/findings_test1.json --pretty
python sns/score_release.py --key data/tcga-brca/releases/test-1.ANSWER_KEY.json --findings sns/findings_test1.json
```
Do not open the `ANSWER_KEY.json` file before running the tool. Any seed number gives a different random release.

**Step 6 (optional, slow). Many random releases at once**
```
python sns/chaos.py --n 10 --seed 21
```

---

## 4. Web app version: the order to test it

The web app is a FastAPI server (`sns/server.py`) with a browser page (`sns/web/index.html`). The page loads React, Plotly and fonts from a CDN, so it needs an internet connection. If there is no internet, use the Streamlit or terminal version.

**Step 1. Start the server**
```
python -m uvicorn server:app --app-dir sns --port 8000
```
Open **http://localhost:8000** in a browser. Stop the server with **Ctrl + C**.

**Step 2. Open the practice envelope and run**
Type `practice` as the envelope number and `orange meadow pilot lantern` as the four words, then click **Open envelope and run check**. (Or pick any listed release and click **Run check**.)

**Step 3. Read the results, top to bottom**
- The **verdict**, then **Impact on the analyses** (A1 and A2 on release A, release B and after repair).
- **Proof:** press **Play the story** to watch the swapped, copied and removed tumors move, and A1 return to its release A value after repair.
- **How it was checked**, **Evidence** (every sample involved), **Repairs and limits**.
- **Export:** download `findings.json`, or build the clean cohort for A1 or A2.

**Step 4. A judge's envelope**
Same as step 2 with the judge's number and words.

---

## 5. Streamlit version: the order to test it

**Step 1. Start the page**
```
streamlit run sns/app.py
```
A browser tab opens.

**Step 2. Open the practice envelope**
In the sidebar, section **1. Get a new release**, tab **Sealed envelope**: type `practice` as the number and `orange meadow pilot lantern` as the four words, then click **Unseal**. It shows `Envelope opened`.

**Step 3. Run the check**
In section **2. Choose and run**: choose `envelope: practice` in **New release (B)**. Leave **Trusted release (A)** as it is. Click **Run check**.

**Step 4. Read the results, top to bottom**
- The colored **verdict** box (green Ready, orange Conditional, red Blocked, blue Insufficient Evidence).
- **Watch A1 change:** press **▶ Play the story** to watch the swapped, copied and removed tumors move, and A1 return to its release A value after repair.
- **Analyses:** for A1 and A2, the status, the verdict, whether the repair restores the original, the reason, and a table with the old release, the new release and the after-repair values.
- **Findings:** one box per problem, with what it affects, confidence, whether it needs human review, its limitation, and a table of every sample involved.
- **Proposed repairs**, the **files not in this release**, and the **limitations**.

**Step 5. Downloads and clean cohort**
Click **Download findings.json**. Choose `A1` and click **Build clean cohort**: expect 698 samples and 0 duplicates, then download `cohort_samples.csv`.

**Step 6. A judge's envelope**
Same as step 2 with the judge's number and words, then choose `envelope: E07` (or the number given) and click **Run check**.

Every run is saved with its own name and time in `sns/runs/`, so nothing is overwritten.

---

## 6. How to add a new release

A release is a folder with the same file names as release A: `clinical.tsv`, `survival.tsv`, `protein_rppa.tsv.gz`, `mutation_gene_level.tsv.gz`, `copy_number_gistic.tsv.gz`, and optionally `rna_hiseq.tsv.gz`. Sample IDs must start with `TCGA-`.

**Way 1. A sealed envelope** (judging day)
```
python envelopes/unseal.py <number> <four words>
```
The release appears in `envelopes/opened/E<number>/release/`. Streamlit lists it automatically.

**Way 2. Your own folder** (terminal)
Put the files in any folder and pass it with `--b`:
```
python sns/run_tower.py --b path/to/my-release --out sns/findings_my-release.json --pretty
```
Folders inside `data/tcga-brca/releases/` are also listed automatically in Streamlit.

**Way 3. Upload in Streamlit**
Sidebar, tab **Upload files**: type a name, choose the files, click **Save upload**. The files are saved in `sns/uploads/<name>/`, and the release appears in the list as `uploaded: <name>`.

**Way 4. A random test release with a hidden answer key**
```
python sns/make_test_release.py --seed <any number> --out data/tcga-brca/releases/<name>
```

**Compare any two releases:** A does not have to be the raw release. Example, compare practice-B with an envelope:
```
python sns/run_tower.py --a data/tcga-brca/releases/practice-B --b envelopes/opened/E07/release --out sns/findings_B_vs_E07.json --pretty
```

---

## 7. What the output means

**Findings (F1, F2, ...)** are the problems found. The number changes from release to release. Kinds:
- `sample_swap`: two or more samples traded their data, or a sample now holds another sample's data.
- `duplicate`: a sample ID carries a profile identical (or almost identical) to an existing sample.
- `membership_change`: samples were lost, added or renamed.
- `missing_data`: values, columns or a file needed by an analysis are missing.
- `label_conflict`: values changed in a text table, such as clinical, and match no other sample.
- `other`: numbers changed and match no other sample.

**Analyses (A1, A2)** are the two research calculations we protect. Each one shows its result on release A, on release B and after repair.

**Repairs (R1, R2, ...)** are proposed fixes. They are proposals only: a person must approve or reject each one.

**Verdict rules, for each analysis:**
- **Ready:** no finding touches its input files, columns or samples.
- **Conditional:** it is touched, and its results moved by less than 10%. Fix the data first, then rerun.
- **Blocked:** a result changed sign, or moved by 10% or more.
- **Insufficient Evidence:** a group has fewer than 10 samples, or the analysis cannot run.
The overall verdict is the worst of the analysis verdicts.

Worked example: in practice B, A1 `mean_negative` moved from -0.872 to -0.843. The change is 0.029, and 0.029 divided by 0.872 is 3.3%. That is under 10%, so A1 is Conditional.

**Repair proof:** `repair_restores_before: true` means that after applying our repairs to release B, the analysis gives exactly the same result as on release A. It proves the repair undoes the change. Whether the change was an error or a valid update is a human decision.

---

## 8. A result a judge can verify

In practice B, `TCGA-A8-A09K-01` (ER positive) and `TCGA-B6-A0WX-01` (ER negative) traded all their protein values. ERALPHA:

```
TCGA-A8-A09K-01   raw: 3.43266169795726    practice-B: -1.07005186202866
TCGA-B6-A0WX-01   raw: -1.07005186202866   practice-B: 3.43266169795726
```

Check it yourself:
```
python sns/show_change.py
```

All 131 protein values of both samples match exactly in the other sample, to every decimal. Two real lab measurements never match like that, so this is a swap, not a new measurement.

Effect on A1: difference (ER positive mean minus ER negative mean) 3.051 on raw, 3.019 on practice-B, and 3.051 again after repair.

---

## 9. How we tested it

- **Organizers' reference answer:** our practice B output has the same findings, kinds, severities and A1 numbers as `format/examples/control-tower-practice-B.json`.
- **Organizers' validator:** every output passes `format/validate.py`.
- **Random releases (chaos tester):** 50 random releases with up to 4 hidden changes each, across 13 kinds of change and all files: 124 of 126 changes caught, 0 crashes, 0 invalid files. The 2 misses were near copies in the sparse mutation file, now a stated limitation.
- **Hidden answer keys:** test releases made with `make_test_release.py`, scored after running: all changes caught, no false findings.
- **Independent test:** a separate set of 17 tests designed outside our build: 44 of 44 changes caught, no false findings. It found one bug (reordered clinical columns), which we fixed.
- **Official envelope:** the practice envelope opens with the organizers' `unseal.py`, and our tool gives the correct findings in about 2 seconds.

---

## 10. Limitations

- It sees what changed between two releases, not why it changed.
- It treats release A as the trusted version. Every finding needs human review.
- Tables only. Slide images and other files are not checked.
- Sample IDs must start with `TCGA-`. Other ID formats are not recognized.
- In sparse files (like mutation calls, where most values are 0), a near copy of a sample with very few non-zero values cannot be told apart from a real sample.
- The 10% line between Conditional and Blocked is our own review rule, stated openly.
- Research use only. No clinical claims.

---

## 11. Files

- `sns/run_tower.py`: the tool, terminal version.
- `sns/server.py`: the web app's FastAPI server.
- `sns/web/index.html`: the web app page.
- `sns/app.py`: the Streamlit page.
- `sns/swap_replay.py`, `sns/swap_replay_template.html`: the animated A1 replay in the Streamlit page.
- `sns/engine/`: `release.py` (loading), `diff.py` (finding changes), `impact.py` (analyses, repairs, verdicts), `findings.py` (writing findings.json), `cohort.py` (clean cohort), `viz.py` (data for the web app's charts).
- `sns/cohort.py`: clean cohort export.
- `sns/chaos.py`, `sns/make_test_release.py`, `sns/score_release.py`: testing tools.
- `sns/peek_sample.py`, `sns/show_change.py`: look at real samples.

---

## 12. Third-party libraries and licenses

- pandas (BSD-3-Clause)
- numpy (BSD-3-Clause)
- pyarrow (Apache-2.0)
- jsonschema (MIT)
- cryptography (Apache-2.0 or BSD-3-Clause)
- streamlit (Apache-2.0)
- plotly (MIT)
- fastapi (MIT)
- uvicorn (BSD-3-Clause)
- python-multipart (Apache-2.0)
- Web page, loaded from the jsDelivr CDN: React and ReactDOM (MIT), htm (Apache-2.0), plotly.js (MIT)
- Fonts: IBM Plex Sans and IBM Plex Mono via Google Fonts (SIL Open Font License 1.1)

Data: TCGA via the UCSC Xena TCGA hub and NCI Imaging Data Commons; see `DATA_LICENSES.md`.