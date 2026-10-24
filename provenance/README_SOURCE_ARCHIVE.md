# X4: Stage-wise extreme proximity in HPacker side-chain reconstruction

**Release candidate (not yet deposited).** Independent, exploratory, fixed-input case record.

## Results

- 200 source frames × 6 HPacker stages = 1,200 archived structures.
- Stage 5 general contact below 1.6 Å: **134/200**.
- PHE6–PHE55 ordered recontact: **43/200**; full six-stage alternation **25/200** (10 `010101`, 15 `101010`).
- Exploratory local NN/chi response results are selected-case evidence, **not** a validated unique causal pathway.
- Four selected WT0/WT1 baseline/rescue structures completed a fixed short BioEmu relaxation.
- A bounded heavy-atom protein LJ proxy attributes extreme baseline repulsion to PHE6–PHE55 overlapping atom pairs, **not** to proven whole-solvent-system PE contributions.

## Manuscript and SI

- `manuscript/Technical_Note_RELEASE_CANDIDATE.pdf`
- `supplement/Supplementary_Information_RELEASE_CANDIDATE.pdf`
- `figures/` includes a six-panel WT0 structural depiction **derived from actual stored PDBs**, as well as source plots.
- `data/` includes 200-frame/1,200-stage diagnostics, selected archived WT0 stage PDBs, selected v8.1–v8.4 comparisons, Phase C summaries if saved, and Stage A pair-resolved results.

## Verification

Run `python code/verify_archive.py` from the package root. This checks fixed primary stage counts, patterns, and packaged SHA integrity. This verification is **not** a re-run of HPacker, OpenMM, or Phase 2B patching.

## Source and chronology

Fixed HPacker commit `6cd5ec33961bb2662b51c0a15d11d6fd9162f2a7`; later BioEmu commit `42c551133f9825f4ec01812075b0d50e268c2209`. See `provenance/` for source identification and limits; originals from BioEmu and HPacker are not claimed as our own. The initial X4 core and the later selected relaxation/LJ continuation are distinct. No independent frozen standalone minimum-protocol document was established.

## Status

Research STOP: October 23, 2026. No new MD, Stage B/C, or comprehensive Phase 2B audit. The scientific endpoints are frozen; public deposition requires human verification of final author identity, third-party redistribution rights, acknowledgments, funding/COI, AI declaration, and metadata.


## Licensing, declarations, citation — 2026-10-24 administrative package

The approved scope is specified in `licenses/LICENSE_SCOPE.md` (NTBE) or `LICENSES/LICENSE_SCOPE.md` (X4). Author-owned prose and original figures are CC BY 4.0; original author-owned code is MIT. Upstream benchmark data and other third-party-derived assets are not covered by a blanket license. `AUTHOR_DECLARATIONS_FOR_APPROVAL.md`, `REFERENCES_VERIFIED.md` and `CITATION.cff` provide author-confirmation metadata. **Author declarations and licensing preferences have been affirmed in the chat; legal ownership of third-party-derived inputs remains subject to upstream terms. No public deposition has occurred.**

## Editorial revision (24 October 2026)
The manuscript now gives quantitative local-perturbation and selected NN/chi response results, defines intervention-derived "rescue" before its downstream use, and records that high/low response paths were selected post hoc. SI S4 includes `supplement/Table_S4_selected_response_evidence.csv`. Lag/masking were used in the original exploratory interventions; the archived v8.1 implementation specifies side-chain masking and a two-stage (s−2) lag; see `code/archived_interventions/`. Do not treat the summary tables as a runnable reconstruction of unprovided source scripts. The WT0 residue-17 to residue-55/56 story remains a proposed chain rather than a uniquely certified causal mechanism.


## Consolidated editorial release (24 October 2026)
Technical Note and SI were typeset anew with explicit Discussion/Conclusions, substantive Results and correctly formatted SI tables. See `provenance/FINAL_EDITORIAL_CHANGELOG.md`. For X4, recorded mask/lag intervention rules and fixed BioEmu Phase C AMBER99SB/TIP3P, 300 K/1 bar protocol are described. No new scientific execution or unverified Phase 2B causal certificate is claimed. Publication URLs remain to be assigned.

## Release-ready document and archived intervention sources
The Technical Note now includes complete references, figure captions, author/AI declarations, and explicitly describes the restricted Stage A LJ definition. Four original v8.1–v8.4 intervention scripts are preserved under `code/archived_interventions/` as research-source evidence, **not independently rerun**. SI S4 includes the selection-history caveat.
