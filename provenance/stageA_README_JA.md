# X4 Continuation Stage A — Protein-heavy-atom LJ audit

## Scope
WT0/WT1 baseline and v8.1 rescue, fixed Stage-5 heavy-atom PDB. **No production MD, whole-protein minimization, solvation, or HPacker**. This is an exploratory mechanistic energy diagnostic, not a Phase C PE reproduction or free-energy calculation.

AMBER99SB OpenMM parameters are assigned after `Modeller.addHydrogens` at pH 7, to type the full protein, but energies are computed only for pairs of original 494 heavy atoms on different residues separated by at least one residue. `NonbondedForce` particle mixing and 1–4/other exception parameters are respected. Each pair contributes repulsive `4ε(σ/r)^12`, attractive `−4ε(σ/r)^6` and net LJ. The protein-only `NoCutoff` energy diagnostic excludes PME, solvent, ion, hydrogen LJ, bonded and close-sequence terms; **do not compare numeric sum directly with solvated BEFORE_MIN PE**. This is a strictly limited, mechanistically interpretable subtotal.

## How to run on existing Windows installation
Expand this ZIP under `C:\Users\user\NTBE_X4` as its own folder. In PowerShell:

```powershell
cd C:\Users\user\NTBE_X4\X4_StageA_Protein_LJ_2026-10-23
& "..\.venv_md\Scripts\python.exe" .\stage_a_lj.py --bioemu "..\bioemu"
& "..\.venv_md\Scripts\python.exe" .\stage_a_lj.py --bioemu "..\bioemu" --execute
```

Dry run checks frozen input hashes and BioEmu git commit without constructing an OpenMM system. Execute computes static pair energies in seconds/minutes, without running a BioEmu MD; OpenMM addHydrogens may internally place/optimize added hydrogen positions. Share the `stage_a_results` directory as ZIP. Do not re-run if error; share stderr/log.

## Acceptance criteria
Fixed input hash + BioEmu clean commit, 494 heavy atoms per input, unchanged heavy coords, complete pair-key correspondence between paired structures, exact LJ type and exception parameter equality between paired cases; per-pair CSV + grouped JSON + SHA checksums. If mapping or parameter mismatches arise, script fails closed. No claims about causation of whole solvated PE.

## STOP
Stage A only, do not execute Stage B/C, do not tune force field or implement MD. One run only.


## v2修正（Stage A）
初版は末端THRのOXT欠損によりaddHydrogensで停止。v2では元PDBを変更せず、力場テンプレート照合専用の一時的OXTをメモリ上で追加します。OXTは原子対LJ集計から除外し、固定元494重原子のみを解析します。OXT位置は末端CA/C/Oから構築し、研究上の実測座標とは扱いません。Dry runはテンプレート生成自体を検証しないため、実行時の別のエラーはなお起こり得ます。入力SHA-256およびBioEmu commitは変更しません。
