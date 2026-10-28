# Read-only reproducibility

The full candidate includes the fixed data needed for the checks below. Use Python >=3.10 and NumPy. P1 was executed with Python and NumPy versions recorded in `verification/current/production_environment.json`. No trained model or OpenMM dependency is required for these commands.

```bash
python code/verify_package_hashes.py .
python code/verify_e1_from_pdb.py --root . --out verification/recheck
python code/verify_supporting_records.py --root . --out verification/recheck
```

Hash checking is independent of scientific correctness. Run it before writing recheck outputs. The manifest excludes itself; newly produced files are reported as extras, not silently included in the original checksum set. Input files are never modified by the verifiers.

## Geometry

`verify_e1_from_pdb.py` requires the complete 200×6 grid, per-file SHA/size match, consistent 494-heavy-atom identities, finite coordinates, receipt hashes where present, and the original stage CSV fields. It computes stage prevalence, unions, strict atom-pair recontact, PHE residue-pair recontact and exact alternation IDs. Frozen expectations are in `data/expected_e1.json`. Persistence and different-pair substitution controls distinguish the strict endpoint from the old 153 counter.

## Supporting records

`verify_supporting_records.py` recalculates ten saved-array diagnostics and forced χ records, validates receipt hashes, inspects the one historical full-batch receipt, checks four Phase C log sequences and final PDBs, keeps separate WT0 minimized geometry separate, and recomputes all 464,560 LJ pair rows. P1 uses the archived force-field σ/ε; it does not repeat the historical OpenMM parameter-typing audit. Archived summaries, pair-CSV summation and PDB-double recomputation are reported separately.

## Figures and editable documents

```bash
python source/freeze_tables.py
python source/make_figures.py
python source/author_documents.py
python source/build_pdfs.py
```

The authoring scripts consume `verification/current/` and original saved inputs; they introduce no new scientific results. Figure 2 uses six actual stage PDBs and a fixed orthographic camera. PDF production also requires Matplotlib, Pillow, ReportLab and DejaVu fonts (the builder uses `/usr/share/fonts/truetype/dejavu`, which can be changed for a local installation). Render PDFs with Poppler for visual review. JSON/Markdown are the editable manuscript sources; authoritative table CSV/JSON and figure-source CSV are included.

## Evidence levels and limits

- E1: fresh P1 independent coordinate and strict-sequence calculations.
- E2: fresh P1 saved-array arithmetic, receipt hashes and χ-record consistency. Historical selected one-point full-batch forward evidence; no new P1 forward or full-feedback regeneration.
- E3: fresh P1 log/hash/final-geometry inspection; no new MD and no WT1 minimized coordinates.
- E4: fresh P1 pairwise arithmetic/double-coordinate distances, using archived parameters; historical independent parameter typing retained.

No universal model or MD replay guarantee is made. Archived execution/intervention scripts are supplied for historical transparency and are **not** part of the read-only commands above. Original NN sham provenance and exact internal-coordinate/model parity remain limited as described in SI.

## GitHub-sized candidate

The GitHub candidate contains manuscripts, figures, source code, authoritative tables, expected values, current verification receipts and lightweight original tables/selected structures. It omits the 1,200-PDB corpus, most E2 records and large LJ pair tables. Full commands must be run against the **matching full Zenodo candidate**, not against the existing v1.0.0 DOI package alone. Both candidates are unpublished. The common files match by SHA. No automatic GitHub-to-Zenodo archive action was performed.
