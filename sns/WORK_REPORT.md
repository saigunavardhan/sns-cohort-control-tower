# SNS Cohort Control Tower: Work Report

## Backup

Created before any changes: `C:\Users\gunav\Sem 4\hackathon\sns_backup_20261114_0001.zip`
Size: 2.9 MB, 151 files. Excludes: data/tcga-brca/raw, data/tcga-brca/releases, envelopes/opened, __pycache__, and cohort-control-tower-main.zip.

**File that must not be committed:** `cohort-control-tower-main.zip` (the teammate zip used to source the envelopes folder).

---

## Task 1: Envelopes

The `envelopes/` folder already existed in the repo with sealed files and an older `unseal.py`. The zip `cohort-control-tower-main.zip` contained a newer version.

We extracted the zip to `C:\Users\gunav\Sem 4\hackathon\sns_tmp_extract` (outside the repo), then copied only its `envelopes/` subfolder (unseal.py, SHA256SUMS, README.md, sealed/) into the repo root, overwriting the older files.

The `cryptography` package was already installed (version 50.0.2). No pip install was needed.

We unsealed the practice envelope with:
```
PYTHONIOENCODING=utf-8 python envelopes/unseal.py practice orange meadow pilot lantern
```

This produced 5 files in `envelopes/opened/practice/release/`: clinical.tsv, copy_number_gistic.tsv.gz, mutation_gene_level.tsv.gz, protein_rppa.tsv.gz, survival.tsv. No RNA file, as expected.

---

## Task 2: RNA Rule

**Problem observed:** running on the envelope release produced a spurious F1 `missing_data` finding for `rna_hiseq.tsv.gz`, and coverage changes showed "lost rna_hiseq.tsv.gz" for hundreds of samples.

**Changes made:**

`sns/engine/diff.py` `coverage_changes()`: restricted the A-side presence table to only files also present in B, so missing files from B do not appear as coverage losses.

`sns/engine/release.py` `load_release()`: added a `skip_files` parameter. Files in this set are not read from disk at all.

`sns/engine/findings.py` `Builder`: added `files_not_in_release` list and `_all_analysis_inputs` set. In `file_level()`, files in A but not in B are only reported as findings if they appear in some analysis's `inputs` list. Files that no analysis reads go into `files_not_in_release` instead.

`sns/engine/findings.py` `build_document()`: added `extra_files_not_in_release` parameter for files skipped before loading. Added `files_not_in_release` and `files_not_in_release_note` to `sns_run`.

`sns/run_tower.py`: before loading release A, we scan the analysis source text for file names (regex `\b[\w-]+\.tsv(?:\.gz)?\b`), list files in folder A and folder B, and compute `skip_from_a` = files in A's folder but not in B's folder and not in the analysis source scan. These are skipped when loading A and passed as `extra_files_not_in_release` to `build_document`. This avoids reading the 64 MB RNA file during envelope runs.

`sns/chaos.py` `Base`: reads the analysis source and stores `analysis_inputs` (file names). `op_drop_file`: if the dropped file is in `analysis_inputs`, expects `missing_data` finding as before. If not, expects empty kinds (no finding, file should be in `files_not_in_release`). `score()`: handles empty-kinds expectations by checking that no finding covers the dropped file and that it appears in `sns_run.files_not_in_release`.

**Verification:**
- Practice B: 3 findings unchanged, validator valid, 5.7 seconds.
- Practice envelope: 3 findings, no RNA finding, `files_not_in_release: ['rna_hiseq.tsv.gz']`, validator valid, 2.1 seconds (RNA not loaded).

---

## Task 3: Readable Output

Added `--pretty` flag to `run_tower.py`. With the flag, it prints: verdict, one line per finding with severity, summary, affects, and needs-human-review; then per-analysis status and verdict; then files not in release; then total seconds. Without the flag, the original JSON summary is printed.

Added `--b-name` default: if the B folder is named "release", the parent folder name is used (e.g., "practice" for `envelopes/opened/practice/release`). This gives the correct release name in findings.json for envelope runs.

---

## Task 4: Clean Cohort Export

Created `sns/engine/cohort.py` and `sns/cohort.py`.

The engine defines `BUILTIN_SPECS` for A1 and A2. It parses repairs from `findings.json` per file: excluded samples are removed, flag_for_review samples are kept (and added from A if they were dropped from B), restore_from_previous_release samples stay in the working set (data is restored, IDs unchanged), relabeled samples are renamed to their old IDs. Each spec step records n_rows, n_distinct_samples, and n_distinct_patients. The `join_slides` step checks for row inflation and collapses to one row per sample.

**A1 on practice-B check:** cohort gives exactly 698 samples with 0 duplicates, matching the A1 sample list on release A exactly (verified by set equality).

**A2 on practice-B:** 844 samples, 0 duplicates.

Outputs: `cohort_samples.csv` (sample_id, patient_id, flagged), `cohort_steps.json` (step counts and repairs), `reproduce_cohort.py` (standalone pandas script).

---

## Task 5: SNS README

Created `sns/README.md` with: product and team info, user workflow, setup instructions, run commands (practice-B, envelope, validator, cohort export), verdict levels and finding kinds, the sample swap verification with exact ERALPHA values, limitations, chaos tester results, demo video placeholder, and third-party licenses.

ERALPHA values verified directly from the files:
- TCGA-A8-A09K-01: 3.43266169795726 in raw, -1.07005186202866 in practice-B.
- TCGA-B6-A0WX-01: -1.07005186202866 in raw, 3.43266169795726 in practice-B.

Chaos results (n=10, seed=21): 21 of 21 changes caught (100%), 0 crashes, 0 invalid files, 0 false findings.

Created `sns/requirements.txt` listing pandas, numpy, pyarrow, jsonschema, cryptography.

---

## Acceptance Test Results

1. Practice B: 3 findings (F1 sample_swap high, F2 duplicate high, F3 membership_change medium), verdict Conditional, validator valid. Runtime: 5.7 seconds.
2. Practice envelope: same 3 findings, no RNA finding, rna_hiseq.tsv.gz in files_not_in_release, validator valid. Runtime: 2.1 seconds.
3. Test release seed 101: 2 of 2 changes caught, 0 false findings.
4. Test release seed 303 (RNA deleted): 2 of 2 changes caught, 0 false findings, no RNA finding.
5. Chaos n=10, seed=21: 21 of 21 caught (100%), 0 crashes, 0 invalid files.
6. Cohort A1: 698 samples (exact match with A1 raw), 0 duplicates. Cohort A2: 844 samples, 0 duplicates.
7. git status: no commits made. Changed files: envelopes/README.md, envelopes/SHA256SUMS, envelopes/unseal.py. New files: sns/ (all contents), make_backup.py at repo root.

---

## Decisions Made Without Asking

The task says to use an `--out` path for acceptance test 3 but does not specify the path. We used `sns/findings_test1.json`. We found the same behavior when specifying `sns/findings_test303.json` for test 4.

The `files_not_in_release_note` field is `null` in findings.json when no files are absent (e.g., for practice-B). We chose this over an empty string to keep the output clean.

The `reproduce_cohort.py` script uses hardcoded paths from the `--release` and `--prev` arguments. This is intentional: the script reproduces a specific run, not an abstract workflow.

---

## Open Problems

None. All acceptance tests pass.
