# X4 — Technical Note / SI pre-production numerical, hash and caption audit

**Audit date:** 2026-10-27
**Scope:** existing PDFs, Phase1A independently revalidated archive, WT upstream match JSON, E1 result JSON, previously frozen v0.3 table plan. No new MD/HPacker computation.

## Executive gate

- **Historical source consistency:** PASS for tested Phase1A ZIP ledger (60/60), internal ZIP CRC, saved E1 fixed metrics, E2 archived-response numeric inputs, E3 completion records, E4 pairwise arithmetic, and official WT upstream XTC/topology SHA.
- **New manuscript / new public ZIP QC:** NOT YET TESTABLE. Draft v0.2 English Tables 2–5 and new SI do not yet exist in the material audited. Do not assign release PASS.
- **Table 4:** preserve two matched timepoints, initial and final MD; separately archived WT0 minimized geometry belongs in SI S5.4, not inserted as continuous Phase C checkpoint.

## Executed checks

| Check | Verdict | Direct evidence |
|---|---|---|
| Phase1A内部SHA256台帳 | PASS | 60 entries checked; mismatches=[] |
| Phase1A ZIP CRC | PASS | zipfile.testzip() |
| E1全1200 PDB固定再計算 | PASS | {'structures': 1200, 'general_stage_counts': [183, 148, 140, 128, 130, 134], 'phe655_stage_counts': [39, 22, 25, 23, 29, 21], 'strict_same_atom_recontact': 147, 'phe655_residue_pair_recontact': 43} |
| E3 四実行完走 | PASS | wt0_baseline: 611100 steps/OpenCL; wt1_v8_1_rescue: 611100 steps/OpenCL; wt0_v8_1_rescue: 611100 steps/OpenCL; wt1_baseline: 611100 steps/OpenCL |
| E4 4件のpairwise集計 | PASS | wt0_baseline: 16616827328250.18; wt0_v8_1_rescue: 27071.418583273637; wt1_baseline: 9951309.663635513; wt1_v8_1_rescue: 27582.94926599252 |
| E2 10件保存配列 | PASS | Cases=[('WT_source001', '55', '53'), ('WT_source000', '53', '25'), ('WT_source000', '53', '50'), ('WT_source001', '55', '11'), ('WT_source000', '53', '6'), ('WT_source001', '55', '52'), ('WT_source000', '53', '55'), ('WT_source001', '55', '59'), ('WT_source001', '55', '56'), ('WT_source001', '55', '3')] |
| E2 blocking 5件 | PASS | Cases=[('WT_source001', '[55, 59, 3]', '3'), ('WT_source001', '[55]', '3'), ('WT_source001', '[3]', '3'), ('WT_source000', '[53]', '4'), ('WT_source001', '[59]', '3')] |
| 公式archive WT XTC topology一致 | PASS | [('delta_G//2K2A/samples.xtc', '0f42d6d669f5e1c765d0be9299aafc3f127676800d72bfee0c3c4e511ff790b9'), ('delta_G//2K2A/topology.pdb', 'ad00dcbeec51b8941430bba94eadb4226a1e3e495e6644bce9492a193bf10798')] |
| 旧PDFがv0.2新版表に更新済み | NOT PASS / PENDING | 現行PDF: Table4=False, SI strict147=False |

## Frozen expected scientific displays and location

| Display | Authoritative data | Required caption / precision control |
|---|---|---|
| Fig. 1A / 1B | E1 six stage counts general=183,148,140,128,130,134; PHE6–55=39,22,25,23,29,21 | same 200 source frames, same y axis 0–200; algorithm stages, not physical time; <1.6Å criterion |
| Table 1 / SI S2 | general frame union 196/200; strict exact-atom 147/200; PHE union 67/200; residue-pair 43/200; complete alternation 25/200 | exact-atom recontact versus residue-pair recontact; old 153 is superseded and must not be treated as ordered recontact |
| Fig. 2 / Figure S1 | WT0 six archived stage PDBs, six interresidue distances 4.291,0.633,4.037,0.560,4.030,0.277 Å | independent same atom identifiers per panel, consistent camera and PDB source, 3D distances; Figure S1 horizontal 1.6Å criterion |
| Table 2 | WT0 original vs locally perturbed PHE min 0.277274 vs 1.881472Å; PHE close pairs 5 vs 0; general close pairs 6 vs 1 | 1.6Å strict, source0 only, original vs locally perturbed is not an MD result |
| Table 3 / SI S2B | archived E2 10 path array consistency; four strong selected displays: 0:s4:53→55, 1:s3:55→56/3/59 | raw component difference dimensionless (not a probability); circular χ difference degrees; selected post hoc; SI links to exact result files and replay evidence tier |
| Table 4 | source0 original 6→0; perturbed 1→0; source1 original 3→0; perturbed 1→0 | initial input PDB versus final short MD heavy-atom PDB, not a clean minimization-only study; four selected runs only |
| SI S5.4 | separate WT0 minimized source0 original/rescue 0 pairs; source1 after min unknown | separate M2 run; not Phase C intermediate frame; energy-only AFTER_MIN cannot supply close-pair count |
| Table 5 / SI S6 | 4 archived LJ proxy outputs; dominant same-atom distance comparisons | baseline-dominant atom pair measured in perturbed structure, not necessarily perturbed closest atom; protein heavy-atom proxy ≠ solvated PE or ΔG; float32-nm saved vs double coordinate recomputation |
| SI Table S2A-1/2 | v8.1 lag/mask/pre-Cartesian chi and v8.2 five blocks | distinguish target/perturbed/third/chi override; block-stage table 0:s4:53, 1:s3:3/55/59/combined |
| SI Table S3 | four E3 step counts 611100, 226 energy states, OpenCL single, different waters | 0.5004ns total; restrained O atom-name selection includes water O; differing solvent sizes preclude energy rank inference |
| SI Table S5 | input files & E1 PDB, E2 arrays/scripts/receipts, E3 phase C zipped logs, E4 pairwise & verification code | verification tiers independent computation vs reading archived receipt vs one-point trained forward, never universal replay |

## Mandatory release blockers

1. Render every new table and figure only from a frozen structured row source. Distinguish 2K2A WT frame 0/1 from distinct systems; avoid unintroduced internal labels.
2. Regenerate all values, captions, table indices and references in the **new** Technical Note and SI; check against CSV/JSON with exact units and pre-agreed display rounding.
3. Compute SHA256SUMS for **final** 1200 PDB ZIP contents and verify every entry readably and against underlying frozen PDB outputs. Check stage/file identity and count 200×6, not only ZIP-level hash.
4. Preserve earlier released version; update version/DOI links, scope-specific licenses/NOTICE and exact source Zenodo record; record archive MD5 output explicitly if making that claim.
5. Re-render PDF and inspect figures / table headers / footnotes, including no truncated signs, dates, overlap, or clipping. Do not claim a new PDF has passed these checks before it exists.
6. No replacement of missing WT1 minimization coordinates with a recreated or inferred value; no silent relabeling of separately executed WT0 M2 checkpoint.

## Source locations
- `/mnt/data/X4_Phase1A_Independent_Audit_2026-10-27.zip` / `audit/SHA256SUMS.txt`, `audit/outputs/e2_direct.csv`, `e2_blocks.csv`, `e3_phase_c.csv`, `e4_summary.csv`, `e1_summary.json`.
- `/mnt/data/X4_1200PDB_Verification_Result.json`.
- `/mnt/data/X4_BioEmu_Official_Archive_Match.json`.
- `/mnt/data/X4_Technical_Note.pdf` and `/mnt/data/X4_Supplementary_Information.pdf` are historical v1.0.0, not revised PDF.

## Verdict
**PRE-PRODUCTION NUMERICAL BASIS: CONDITIONALLY VERIFIED. PRE-RELEASE DOCUMENT/PACKAGE QC: PENDING NEW OUTPUTS.**
