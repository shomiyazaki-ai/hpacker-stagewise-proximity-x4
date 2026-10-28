# Stage-wise Extreme Atomic Proximity in HPacker Side-chain Reconstruction

**Published v1.1.0 · 28 October 2026**  
**Citable version DOI:** [10.5281/zenodo.23717073](https://doi.org/10.5281/zenodo.23717073) · **Historical v1.0.0 DOI:** [10.5281/zenodo.23616075](https://doi.org/10.5281/zenodo.23616075)

**Sho Miyazaki** — Independent researcher, Tokyo, Japan

## At a glance

This fixed-archive case study examines **200 stored WT 2K2A BioEmu frames × six HPacker reconstruction stages = 1,200 PDB structures**. An extreme proximity is a **nonadjacent interresidue heavy-atom distance strictly below 1.6 Å**. Reconstruction stages are *algorithmic updates*, not physical time. This is not a cross-protein benchmark or a general HPacker failure-rate estimate.

| Fixed-archive E1 endpoint | Frames / 200 |
| --- | ---: |
| Any general extreme proximity at any stage | 196 |
| Ordered recontact of the **same atom pair** (present → absent → present) | **147** |
| Ordered PHE6–PHE55 **residue-pair** recontact (atoms may differ) | **43** |
| Complete PHE6–PHE55 six-stage alternation | **25** (010101: 10; 101010: 15) |

Stagewise general counts (stages 0–5): **183, 148, 140, 128, 130, 134**. The historical count **153** did not require intervening absence and **must not** be substituted for the corrected strict 147. Results beyond E1 are **selected exploratory cases**: E2 saved-array NN responses and interventions, E3 four restrained short OpenMM continuations, and E4 a restricted protein heavy-atom LJ proxy. They do **not** prove a general packing improvement, a unique NN causal mechanism, higher MD success probability, or thermodynamic stabilization.

## Read and inspect the evidence

| What you need | Where to find it |
| --- | --- |
| Main claims, figures, limitations | [Technical Note v1.1.0 (PDF)](manuscript/Technical_Note_v1.1.0.pdf) |
| Complete methods, full tables, evidence boundaries | [Supplementary Information v1.1.0 (PDF)](supplement/Supplementary_Information_v1.1.0.pdf) |
| **Complete data and verification inputs** (including all 1,200 PDBs) | [Published Zenodo v1.1.0](https://doi.org/10.5281/zenodo.23717073) — download `X4_v1.1.0_Zenodo_DOI23717073_PREPUBLICATION.zip` (the filename records its *packaging* state; the record itself is published) |
| Read-only verification scripts | [E1 coordinate check](code/verify_e1_from_pdb.py) · [E2–E4 evidence checks](code/verify_supporting_records.py) · [package hashes](code/verify_package_hashes.py) |
| Frozen targets and numerical outputs | [E1 expected values](data/expected_e1.json) · [authoritative tables](data/authoritative/) · [saved verification receipts](verification/current/) |
| Claim-to-evidence paths, input identities and limitations | [Provenance](PROVENANCE.md) · [SI Table S5](supplement/Supplementary_Information_v1.1.0.pdf) · [full walkthrough](REPRODUCIBILITY.md) |
| Rights, previous version and change history | [Rights scope](LICENSES/LICENSE_SCOPE.md) · [NOTICE](NOTICE.md) · [Changes](CHANGELOG.md) · [GitHub-vs-Zenodo scope](GITHUB_PACKAGE_SCOPE.md) |

The original v1.0.0 manuscript (`manuscript/Technical_Note.pdf`) and SI (`supplement/Supplementary_Information.pdf`) remain under their historical paths; use the **versioned v1.1.0 PDFs** for the current results.

## Verify the published results: three steps

**Required for the pinned example below:** Python **3.11+** and NumPy **2.3.5** (the verification scripts are written for Python 3.10+, but this specific NumPy release requires a newer Python), approximately 70 MB to download the archive plus space for extraction and outputs. **Environment previously used:** Python **3.12.14**, NumPy **2.3.5**, Linux (see [`verification/current/production_environment.json`](verification/current/production_environment.json)). A new clean-environment replay on other operating systems is **not claimed**. No GPU, HPacker weights, learned model inference, OpenMM, MDTraj or MD simulation is needed for these **read-only saved-evidence checks**.

1. **Download** the file `X4_v1.1.0_Zenodo_DOI23717073_PREPUBLICATION.zip` from the [v1.1.0 Zenodo file list](https://zenodo.org/records/23717073). **Extract the ZIP into an empty directory**, keeping its internal folders. Open a terminal **in the extracted directory**: it should contain `README.md`, `SHA256SUMS.txt`, `code/` and `data/pdb/`.
2. **Install NumPy in a separate environment** (example below), then run the three verification commands **from that extracted directory**. Use an unused output directory. `verification/recheck/` must be new or clear of earlier generated rechecks; this does not modify the archived input files.
3. **Inspect the results**: the package hash script should report `"status": "PASS"` with no `mismatches`; the E1 script should report 1,200 structures, 147 strict same-atom recontacts, 43 PHE recontacts, and the stage counts above; the supporting script should report `"status": "PASS"`, 10 E2 saved-array comparisons, four E3 runs and 464,560 E4 pair rows. Find JSON/CSV outputs in `verification/recheck/`.

```bash
python -m venv .venv-x4
# Activate .venv-x4 using your operating system's standard command.
python -m pip install numpy==2.3.5

python code/verify_package_hashes.py .
python code/verify_e1_from_pdb.py --root . --out verification/recheck
python code/verify_supporting_records.py --root . --out verification/recheck
```

**Important:** The first command compares *file hashes*, not scientific correctness. The E1/E2–E4 commands independently recompute **defined numerical endpoints from archived evidence**, not the original BioEmu samples, trained HPacker/NN forward passes, or molecular-dynamics trajectories. The hash checker can report unlisted extras such as a newly created virtual environment without failing if original files match; for a clean listing run it before creating the virtual environment or put the virtual environment outside the extracted directory. The E4 script reuses archived force-field parameters and does not retype them. See [detailed interpretation, output inventory, and known limits](REPRODUCIBILITY.md).

**Do not run the three full-data commands against a GitHub clone alone.** GitHub intentionally excludes the complete 1,200-PDB corpus and large supporting data. Obtain the matching Zenodo v1.1.0 ZIP; the v1.0.0 archive is *not* a substitute.

## Scope of this revision

v1.1.0 corrects the definition/reporting of strict same-atom recontact (147, not historical 153), improves provenance, Methods, figures and full saved-evidence verification, while **preserving historical research inputs and execution records**. It introduced no new BioEmu generation, HPacker inference, NN forward, minimization or MD. Source documents and plotting/build scripts are provided; reproducing the rendered PDF also needs plotting/document dependencies and fonts described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Citation, rights, acknowledgments

Cite **Miyazaki, S. (2026). *Stage-wise Extreme Atomic Proximity in HPacker Side-chain Reconstruction* (Version 1.1.0). Zenodo. https://doi.org/10.5281/zenodo.23717073**. The root [`CITATION.cff`](CITATION.cff) identifies this version; [v1.0.0](https://doi.org/10.5281/zenodo.23616075) remains separately citable.

Original author-written prose/figures: **CC BY 4.0**; eligible original author code: **MIT**; upstream BioEmu coordinates, code, and other third-party materials retain their respective scoped terms. **The entire archive is not uniformly CC BY or MIT.** See [LICENSE_SCOPE](LICENSES/LICENSE_SCOPE.md), [NOTICE](NOTICE.md), and the [rights manifest](provenance/Rights_Manifest.csv). BioEmu, H-Packer and OpenMM developers are credited in the Technical Note. Generative AI assistance and author responsibility are disclosed there; no upstream endorsement is implied.
