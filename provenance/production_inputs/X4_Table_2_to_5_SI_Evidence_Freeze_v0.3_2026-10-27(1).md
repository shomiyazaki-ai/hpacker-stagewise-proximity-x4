# X4 Table 2–5 / SI evidence freeze (proposed v0.3, 2026-10-27)

**Status:** numerical/source audit, NOT a new scientific run and NOT final manuscript proof. Values are derived from saved source records plus the Phase 1A independent audit. Discrepancies and absent checkpoints remain explicit.

## Table 2 — WT source frame 0, local perturbation (E2)

| Metric | Original stage 5 | Intervention-derived stage 5 | Source |
|---|---:|---:|---|
| PHE6–PHE55 min heavy-atom distance (Å) | 0.277274 | 1.881472 | original PDB / independently inspected Phase1A E3-M2 inputs and E2 outputs |
| PHE6–PHE55 atomic pairs <1.6 Å | 5 | 0 | stage PDB / independent geometry |
| General nonadjacent heavy-atom pairs <1.6 Å | 6 | 1 | stage PDB / independent geometry |

IMPORTANT: This uses the pre-relaxation selected stage-5 structures. In E3, the intervention-derived structure has **global** minimum distance 1.550220 Å, even though PHE6–PHE55 distance is 1.881472 Å. Thus do not confuse PHE separation with removal of every <1.6Å contact.

## Table 3 — selected 10 NN-response paths (E2)

| Source | Stage | Donor→receiver | Post hoc group | Maximum |raw 8-vector component difference| | Maximum circular χ difference (°) |
|---|---:|---|---|---:|---:|
| 0 | 4 | 53→6 | Selected weak | 0.035291 | 0.239754 |
| 0 | 4 | 53→25 | Selected weak | 0.024820 | 0.053169 |
| 0 | 4 | 53→50 | Selected weak | 0.005751 | 0.082947 |
| 0 | 4 | 53→55 | Selected strong | 0.141019 | 54.088766 |
| 1 | 3 | 55→3 | Selected strong | 1.369136 | 148.689018 |
| 1 | 3 | 55→11 | Selected weak | 0.038130 | 0.305740 |
| 1 | 3 | 55→52 | Selected weak | 0.013456 | 0.118828 |
| 1 | 3 | 55→53 | Selected weak | 0.080839 | 0.079979 |
| 1 | 3 | 55→56 | Selected strong | 1.671166 | 178.328796 |
| 1 | 3 | 55→59 | Selected strong | 0.752588 | 52.082260 |

**Table 3 limitation:** Groups are retrospectively selected examples, not prospectively thresholded classes. Differences derive from saved model arrays; the Phase1A audit did not independently execute the full trained NN forward pass. Raw-output differences are dimensionless, not probabilities. χ circular differences are not chemically symmetry-corrected.

## Table 4 — selected OpenMM geometry checkpoints (E3)

| Source | Structure | Before min (pairs <1.6 Å) | After min (pairs <1.6 Å) | Final MD (pairs <1.6 Å) | Completion |
|---|---|---:|---:|---:|---|
| 0 | Original | 6 | 0 | 0 | COMPLETE |
| 0 | Locally perturbed | 1 | 0 | 0 | COMPLETE |
| 1 | Original | 3 | **NOT VERIFIED** | 0 | COMPLETE |
| 1 | Locally perturbed | 1 | **NOT VERIFIED** | 0 | COMPLETE |

**Checkpoint caveat:** WT0 after-min values verified from separate M2 minimized PDBs, whereas Phase C is a distinct continuation. Do not imply a verified three-checkpoint continuous trajectory in a single run until Phase C per-run minimized coordinates and identity are verified. Input 3/1 for WT1 was computed from the Stage A selected PDBs; **all four input PDBs in the Phase C kit and Stage A v2 kit are byte-identical by independently calculated SHA-256**. Matching these kits does not alone prove exact run-time loader inputs absent a per-run file access trace. Final near-pair zeros are inherited from the archived final-pdb geometry verification; archive paths/checkpoint must be collated before QC.

## Table 5 — restricted protein-heavy-atom Lennard–Jones proxy (E4)

| Source | Structure | Archived net LJ proxy (kJ/mol) | Recomputed from original PDB double-coordinates (kJ/mol) | Heavy-atom pairs <1.6 Å |
|---|---|---:|---:|---:|
| 0 | Original | 16616827328239.929688 | 16614078362116.181641 | 6 |
| 0 | Locally perturbed | 27071.418583 | 27071.324155 | 1 |
| 1 | Original | 9951309.663635 | 9951099.620821 | 3 |
| 1 | Locally perturbed | 27582.949266 | 27584.089713 | 1 |

Baseline-dominant same atomic pair tracked in source 0: PHE6 CD1–PHE55 CE2, 0.277→5.514 Å; ~99.9999991% of the defined proxy reduction. Source 1: PHE6 CG–PHE55 CE2, 0.916→6.336 Å; ~99.1753%. **These distances are the same tracked pair, not necessarily the closest residue pair after perturbation.**

Archived LJ values follow an nm-float32 coordinate pathway; direct PDB double-coordinate calculation is not bit-identical, especially for 0.277 Å contacts. No energy stability/free energy conclusions.

## SI S1–S7 content and evidence freeze

### S1 — Origin/lineage
Original source 2K2A XTC and topology are identical by full SHA256 to internal BioEmu archive entries delta_G//2K2A/samples.xtc and topology.pdb; official ZIP full MD5 matching remains a separate receipt to retain. XTC 258900B sha256 0f42d6d669f5e1c765d0be9299aafc3f127676800d72bfee0c3c4e511ff790b9; topology 25163B sha256 ad00dcbeec51b8941430bba94eadb4226a1e3e495e6644bce9492a193bf10798. Upstream coordinate rights require appropriately scoped dataset license/NOTICE. 200 saved frames, 310 input atoms, 252 backbone atoms; 494 heavy atoms in reconstructed stage output.

### S2 — Table S1: all stage counts
| Stage | General frames /200 | PHE6–55 frames /200 |
|---:|---:|---:|
| 0 | 183 | 39 |
| 1 | 148 | 22 |
| 2 | 140 | 25 |
| 3 | 128 | 23 |
| 4 | 130 | 29 |
| 5 | 134 | 21 |

Frame union 196/200 general, 67/200 PHE; strict same-atom ordered recontact **147/200**; PHE residue-pair ordered recontact **43/200**; six-stage alternation 25/200 (010101 10, 101010 15). Historical same-atom counter 153/200 counted persistence and is superseded for this endpoint.

### S3 — Figure S1 and full alternation
WT source 0 PHE6–55 stage distances (Å): 4.291, 0.633, 4.037, 0.560, 4.030, 0.277. Threshold <1.6 Å. Six-stage pattern IDs: 010101 = 0,10,34,53,73,83,88,167,180,191; 101010 = 14,28,49,60,66,77,113,126,134,137,141,145,160,165,190. Include 3D min-pair atom names per stage from data/all_1200_stages.csv at final QC.

### S4 — Intervention and NN response evidence
Split S2A as **S2A-1 masking/lag & pre-Cartesian chi rescue**, **S2A-2 downstream residue chi blocking/propagation**; use separate stage, operation, affected residue and role columns, avoid misleading donor/mediator equivalence. Original SI: source 0 v8.1 stage 2 target6/intervened55/third58/mediator17 lag branch; source1 stage1 target55/intervened6/mediator56 mask branch. Downstream v8.2 original receipt cases (rescue-background targeted block): source0 at stage4: residue 53; source1 at stage3: residue 3, 55, 59, and combined 55+59+3. Each receipt specifies `block_residues`, PASS output hashes, no MD, and explicitly warns that this is path-specific, not unique necessity/sufficiency. The stage is visible in case ID and stored stage PDB; establish whether the intervention operation is pre-Cartesian χ override from the fixed runner before finalizing method wording. Table S2B has 10 paths and more columns than Table 3: Zernikegram L2, raw-max, chi-max, archived result link, evidence tier.

| Source | Stage | Donor→receiver | Zernikegram L2 | raw-max | χ-max (°) | Evidence scope |
|---|---:|---|---:|---:|---:|---|
| 0 | 4 | 53→6 | 0.056752 | 0.035291 | 0.239754 | Saved-array independent consistency; NOT full independent forward |
| 0 | 4 | 53→25 | 0.041523 | 0.024820 | 0.053169 | Saved-array independent consistency; NOT full independent forward |
| 0 | 4 | 53→50 | 0.036039 | 0.005751 | 0.082947 | Saved-array independent consistency; NOT full independent forward |
| 0 | 4 | 53→55 | 0.067109 | 0.141019 | 54.088766 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→3 | 0.211317 | 1.369136 | 148.689018 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→11 | 0.100229 | 0.038130 | 0.305740 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→52 | 0.024546 | 0.013456 | 0.118828 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→53 | 0.137227 | 0.080839 | 0.079979 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→56 | 0.198319 | 1.671166 | 178.328796 | Saved-array independent consistency; NOT full independent forward |
| 1 | 3 | 55→59 | 0.231874 | 0.752588 | 52.082260 | Saved-array independent consistency; NOT full independent forward |

Phase 2B saved 16/16 patched donor-output comparisons exact, but source raw arrays for identity/read-only sham controls unavailable for independent reconstruction. New selected one-point CPU replay must be described with its exact receipt; do not extend globally.

### S5 — Table S3: 4 Phase C records

| Case | Particles | Waters | Steps | Energy states | Final PE (kJ/mol) | Platform |
|---|---:|---:|---:|---:|---:|---|
| wt0_baseline | 11164 | 3395 | 611100 | 226 | -143055.181703 | OpenCL/single |
| wt0_v8_1_rescue | 11461 | 3494 | 611100 | 226 | -146432.789062 | OpenCL/single |
| wt1_baseline | 11404 | 3475 | 611100 | 226 | -145726.650193 | OpenCL/single |
| wt1_v8_1_rescue | 11383 | 3468 | 611100 | 226 | -146144.451858 | OpenCL/single |

Each run: exit0, 204 begin/end calls, nonfinite flags0, 611100 steps = 0.5004ns (not 0.6111ns); NVT 0.1ns + NPT 0.4ns + initial four 0.1ps stages. Note backbone-named atom restraint in fixed code also selects water oxygen O; cannot call free MD. Differing water counts preclude direct final PE comparison.

### S6 — Table S4, LJ precision
Analysis 494 common original heavy atoms; 116140 pairs/case and 464560 total, per-forcefield typing verified with 952 particles and 5144 exceptions per geometry. OXT added temporarily for typing only. Pair-contribution delta is independent proof of restricted calculation, not solvated PE. Include archived and double-coordinate values plus epsilon/sigma mixing and pair exclusions in SI. Figure S2 graph log10 of strictly positive proxy.

### S7 — Claim-evidence-file-verification
Preserve file-level mapping: E1 stage PDB + verify_e1_frames.py, E2 original archives + e2_direct.csv/e2_coordinate_chi.csv, E3 Phase C 4 ZIP + e3_phase_c.csv, E4 stage_a_results.zip + e4_summary.csv. Distinguish primary saved records, retrospective independent calculations, editorial displays, and unresolved original-run receipts.

## Remaining release-blocking validation
1. **COMPLETED for input kits:** all four Phase C kit PDBs are SHA-256 identical to Stage A v2 kit PDBs. Remaining: map the kit files to the actual inputs loaded in each Phase C run (run-time provenance), when available.
2. Locate WT1 minimized structures if present; otherwise set S5/Table4 checkpoint cells to ND and caption to “verified separately for WT0; unavailable for WT1” (do not claim triple-series completeness).
3. Verify Figure2 WT0 stage geometry and Figure S1 atom IDs/positions directly from PDB.
4. **RESOLVED at receipt level:** five v8.2 conditions confirmed (source0 stage4 block53; source1 stage3 block3, block55, block59, combined block55+59+3). Remaining: compare code branch details/terminology before final text.
5. Ensure Table 3 raw8 labels are not confused with probabilities, and confirm one-point CPU forward replay uses exact same arrays/mapping.
6. Final license notice and archived ZIP MD5 record, output file hashes, and table manuscript numeric parity.
