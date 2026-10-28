# GitHub vs Zenodo distribution scope — X4 v1.1.0

**Status:** GitHub v1.1.0 materials and [Zenodo v1.1.0](https://doi.org/10.5281/zenodo.23717073) are publicly available (28 October 2026). This document describes why the two distributions have different file sets; the uploaded Zenodo archive itself is frozen.

| Scope | This GitHub repository | Published Zenodo ZIP |
| --- | --- | --- |
| Current Technical Note, SI and figures | Yes | Yes |
| Verification scripts and authoritative summary tables | Yes | Yes |
| All 200 × 6 = 1,200 stage PDBs | **No** | **Yes** (`data/pdb/`) |
| Complete archived E2 arrays and original interventions | Selected lightweight files only | Yes |
| E3 selected OpenMM logs/continuation structures | Selected historical/summary files | Full saved evidence |
| E4 large pairwise LJ records | Lightweight outputs only | Complete archived pair CSVs |
| Full read-only verification of E1–E4 | **Not from GitHub alone** | **Yes**, within limits in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) |

**Start here:** Download `X4_v1.1.0_Zenodo_DOI23717073_PREPUBLICATION.zip` from the [published Zenodo file list](https://zenodo.org/records/23717073), extract it without flattening paths, and follow the three commands in [README](README.md) / [REPRODUCIBILITY](REPRODUCIBILITY.md). The word `PREPUBLICATION` in the filename describes the build-time packaging state, **not** the publication status of the Zenodo record. The [historical v1.0.0 record](https://doi.org/10.5281/zenodo.23616075) is not interchangeable with the corrected v1.1.0 inputs.

The GitHub update intentionally preserves the existing v1.0.0 Technical Note and SI at `manuscript/Technical_Note.pdf` and `supplement/Supplementary_Information.pdf`. Current-version PDFs have `_v1.1.0` in their filenames. Historical receipt dates and historical status narratives remain historical documents rather than the authoritative current release status.

After this three-document GitHub-only usability update, the **GitHub** checksum manifest is refreshed for the amended tree. The **Zenodo** ZIP and its internal checksum manifest are **not** altered. GitHub documentation improvements must not be described as new scientific execution, new model inference, or an amended Zenodo deposit.
