# Reproducibility — X4 v1.1.0 (published 28 October 2026)

This document is the practical companion to the [README](README.md) and [Supplementary Information](supplement/Supplementary_Information_v1.1.0.pdf). It explains how to **independently verify numerical statements from archived coordinates, arrays and logs**, rather than how to regenerate the BioEmu ensemble, HPacker stages, trained NN responses or MD trajectories. The published version DOI is [10.5281/zenodo.23717073](https://doi.org/10.5281/zenodo.23717073).

## 1. Choose the correct data distribution

| Distribution | Contains | Appropriate use |
| --- | --- | --- |
| [Zenodo v1.1.0](https://zenodo.org/records/23717073), file `X4_v1.1.0_Zenodo_DOI23717073_PREPUBLICATION.zip` | All **1,200 original stage PDBs**, E2 saved arrays/interventions, E3 run outputs, E4 pair tables, code, SHA receipts | **Required** for the complete three-command read-only verification |
| This GitHub repository | Readable manuscript/SI, figures, scripts, selected lightweight tables/structures, historical v1.0.0 files | Browse the study and inspect scripts; **not** sufficient on its own for full checks |
| [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.23616075) | Historical previous release | Historical comparison, **not** a replacement for the revised v1.1.0 verification inputs |

The word `PREPUBLICATION` in the **immutable uploaded ZIP filename** describes the build-time packaging state; the Zenodo **record is published**. No rerun requires renaming the archive.

## 2. Install and run: clean path

Use **Python >=3.11** with the pinned **NumPy 2.3.5** installation shown below; the verification scripts themselves are written for Python >=3.10, but NumPy 2.3.5 does not support all those Python versions; the previously recorded validation environment was **Linux / Python 3.12.14 / NumPy 2.3.5** (`verification/current/production_environment.json`). Matching those versions is a sensible starting point, **not a claim of a new independent clean-room run on all platforms**. These three scripts have no model-weight, OpenMM, MDTraj or GPU requirement. An archive of ~69 MB expands to substantially more; ensure free disk space for the data and outputs.

1. Download the named ZIP from Zenodo v1.1.0 and unpack it into an empty directory using any standard ZIP tool (including Windows Explorer). **Do not flatten the directory tree.**
2. In a terminal, change to the directory containing the extracted `SHA256SUMS.txt`, `code/`, `data/` and `provenance/`. Commands below assume this is the **working directory**.
3. For an isolated environment, create and activate a virtual environment; put it **outside the archive directory** if you prefer a clean `unlisted_extra_files` field in the hash check. The platform-specific activation step is optional for an already isolated Python environment.
4. Verify hashes **before creating output files**. Then run E1 and the supporting verification.

```bash
python -m venv ../x4-verify-env
# Activate ../x4-verify-env using your OS (e.g. source ../x4-verify-env/bin/activate on Linux/macOS,
# or ..\x4-verify-env\Scripts\Activate.ps1 in Windows PowerShell).
python -m pip install numpy==2.3.5

python code/verify_package_hashes.py .
python code/verify_e1_from_pdb.py --root . --out verification/recheck
python code/verify_supporting_records.py --root . --out verification/recheck
```

All three checks were implemented as **read-only checks of the original evidence**. The two numerical checkers write **only** their `--out` directory; choose a fresh path if you want to retain an earlier rerun. For a Windows system with `py -3`, `python` above can be replaced by the active environment's Python executable. No command here runs learned inference or MD.

## 3. Expected outputs, tests and failure criteria

| Command | Successful indication | Output |
| --- | --- | --- |
| `verify_package_hashes.py` | JSON `status: PASS`; `mismatches: []`; lists any unmanifested newly created files separately | Console JSON only (no output directory) |
| `verify_e1_from_pdb.py` | JSON `status: PASS`, `structures: 1200`, `atoms_each: 494`, strict recontact **147**, PHE recontact **43**; asserts all expected numbers | `verification/recheck/e1_verification.json`, `e1_stages.csv`, `e1_frames.csv`, `e1_transitions.csv` |
| `verify_supporting_records.py` | JSON `status: PASS`, `e2_paths: 10`, `e2_block_cases: 5`, `e3_runs: 4`, `e4_pairwise_total: 464560` | `verification/recheck/supporting_verification.json`, `e2_paths.csv`, `e2_edits.csv`, `e2_blocks.csv`, `e3_runs.csv`, `e3_separate_m2.csv`, `e4_lj.csv`, `e4_tracked_pairs.csv` |

The E1 expected values are pinned in [`data/expected_e1.json`](data/expected_e1.json). They include stagewise general counts **183/148/140/128/130/134**, PHE counts **39/22/25/23/29/21**, unions **196** and **67**, **147** strict same-atom ordered recontacts and **43** PHE residue-pair recontacts, and the 25 alternating-frame IDs. Synthetic checks in the script reject persistent contact and replacement by a different atom pair as strict same-atom recontact. The historical `153` field is intentionally retained only as historical data.

E1 **fails** when the 200×6 grid is incomplete, source hashes change, atom identities are duplicated or nonfinite, historical stage values diverge, or any expected summary/assertion fails. The supporting script likewise stops on missing evidence or numeric/log/identity disagreement. The hash-only checker is **not** an independent validation of scientific conclusions; a hash pass alone is insufficient.

The archived `verification/current/` outputs are *prior verification records*, not automatically the outcome of your new run. Compare your new outputs to them and the manuscript/SI. Numeric outputs may reflect specified floating-point paths; precision differences are documented in SI S6 rather than silently equated.

## 4. What is and is not independently checked

| Evidence tier | Current recheck | Boundary |
| --- | --- | --- |
| **E1** | Reparse all 1,200 saved PDBs, verify SHAs, calculate heavy-atom geometry, stage prevalence and strict recontact; test definitions | Does **not** regenerate BioEmu sampling or run HPacker |
| **E2** | Recalculate ten saved-array differences, confirm preserved χ records and SHA receipts; inspect one saved full-batch forward receipt | Does **not** run ten new NN forward passes or establish unique cause |
| **E3** | Verify four archived execution-log sequences, per-run hashes, final PDB proximity, and separately saved WT0 minimized PDB geometry | Does **not** reproduce MD trajectories or infer missing WT1 minimized contacts; no robustness inference |
| **E4** | Recalculate **464,560** pairwise LJ values and saved-coordinate / PDB-double comparisons from archived σ/ε | Does **not** independently construct or retype a force field; no full potential energy or ΔG |

Relevant evidence crosswalk: SI Table S5; `data/authoritative/`; `provenance/PDB_SHA256_Manifest.csv`; `provenance/Rights_Manifest.csv`; and the frozen source receipts. E3 original and perturbed systems have different solvent counts and were run under restraints, so final total PE values are **not** comparative stability estimates. E4 values are highly sensitive to very short distances and coordinate precision.

## 5. Figures, manuscripts and optional document rebuild

The editable materials reside under `source/`; the published PDFs are under `manuscript/` and `supplement/`. The historical source and rendering commands are:

```bash
python source/freeze_tables.py
python source/make_figures.py
python source/author_documents.py
python source/build_pdfs.py
```

These **are not part of the standard three-command verification**. They need Matplotlib, Pillow, ReportLab and accessible DejaVu TTF fonts; the builder currently uses `/usr/share/fonts/truetype/dejavu`, so Windows users must adapt font paths. Rendered PDF visual review also uses Poppler. Rebuild differences due to fonts/layout do not change the frozen scientific result; do not claim bitwise PDF reproduction across platforms. Figures link to preserved source coordinates and table inputs; no new experiment is required.

## 6. Integrity, safety and historical files

`SHA256SUMS.txt` at the top of the extracted **Zenodo ZIP** belongs to that **frozen Zenodo distribution**. The root GitHub `SHA256SUMS.txt` covers the GitHub tree and may change after GitHub documentation-only updates; **do not cross-apply the manifests** to different distributions. `verify_package_hashes.py` excludes the manifest itself from the list of digested files and reports extra outputs separately. The README, manuscript and SI inside the published Zenodo ZIP remain exactly as deposited; subsequent GitHub documentation corrections are not a retroactive Zenodo change.

Original code and source records are preserved for transparency, but archived execution/intervention scripts are **not** the recommended verification commands. Review scripts before running any historical execution path. Original sham-array provenance, full model parity, force-field retyping and unrecorded MD trajectories remain outside this verification claim. See [`GITHUB_PACKAGE_SCOPE.md`](GITHUB_PACKAGE_SCOPE.md), the Technical Note and SI for scientific and distribution boundaries.
