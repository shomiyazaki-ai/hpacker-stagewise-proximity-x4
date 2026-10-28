# Provenance

This is an unpublished revised candidate, version 1.1.0. Existing v1.0.0 is at https://doi.org/10.5281/zenodo.23616075. No new DOI is asserted.

## Input audit

The named corrected outer ZIP and its README were not attached. All 14 constituent files were available individually and exactly matched the supplied `ATTACHMENT_MANIFEST.csv` by size and SHA-256. Name suffixes and renamed complete Phase C/verification files were resolved by hash, without substituting content. The 4 reference PDFs informed editorial design only. Full input inventory, recursive archive index and source hashes are in `provenance/`.

Recursive extraction treated Windows backslash members as path separators in working copies, rejected traversal/absolute paths, checked CRC and grouped duplicates by SHA. Original source archive bytes were not rewritten. Duplicate files or nested copies are not independent scientific evidence.

## E1 lineage

BioEmu supplementary archive `samples_delta_G.zip`, WT entry `delta_G//2K2A/`, is identified by the historical official-entry-match JSON. The saved XTC/topology hashes match that receipt; P1 independently verified byte hashes and XTC frame headers (200 frames, 310 atoms), topology (63 residues) and 252 N/CA/C/O backbone atoms. No XTC coordinates were decompressed or new structures generated in P1.

HPacker commit: `6cd5ec33961bb2662b51c0a15d11d6fd9162f2a7`. Candidate PDBs match all 1,200 original stage PDBs **and** their original per-source receipts. The candidate's former NOT_FOR_PUBLICATION filename was not treated as authentication; byte comparisons and independent geometry were required. `PDB_SHA256_Manifest.csv` records the exact source/stage grid and hashes.

## Selected supporting records

- E2 original intervention/diagnostic archives remain byte-identical under `provenance/e2_original_archives/`. Normalized extracted records are in `data/e2/`. Historical mutant examples inside untouched source archives are outside the WT200 analysis. The one selected full-batch receipt was inspected, not rerun.
- E3 full Phase C kit and all four output records are in `data/e3/phase_c/`. Runtime parent `environment.json` input hashes match both Phase C and Stage A inputs. The worker START field is `verified_by_parent`, not a separate OS file-access trace. Separate WT0 M2 minimization-only PDBs remain separately labeled.
- E4 original Stage A inputs, pair CSV, differences, summary and hashes are retained in `data/e4/`. P1 remeasured distances and recalculated the archived pair formulas; prior independent force-field typing is historical audit evidence, not a new P1 operation.

## Current-source access

P1 read the current public GitHub README, LICENSE_SCOPE and CITATION using the connected read-only API. Their Git blob identifiers are recorded in `provenance/Live_Source_Checks.json`. The entire live repository was not cloned or exhaustively byte-compared. The supplied v1.0.0 archive is the byte-level baseline for `X4_Release_Delta_Manifest.csv`.

The supplied `Graylab/HPacker` location returned 404 through the available route. The primary paper, archived provenance and `gvisani/hpacker` README identify the relevant upstream repository; that source is cited in the candidate. BioEmu `md_utils.py` at the fixed commit was read and retained with its Microsoft MIT notice.

Zenodo record/API retrieval failed (including rate limiting for upstream record 15672282). P1 therefore does not claim a fresh live Zenodo DOI/file-list/metadata check. The whole 1.96-GB upstream archive was not downloaded and its MD5 was not independently recalculated. This does not alter the saved-input SHA comparison; it remains a disclosed prepublication metadata/source-access check.
