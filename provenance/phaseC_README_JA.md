# X4 Phase C — 固定BioEmu緩和MD 4構造（逐次）

対象は WT0 baseline / WT0 v8.1 rescue / WT1 baseline / WT1 v8.1 rescue。各ケース2400秒（40分）で終了し、**自動的に次のケースには進まない**。途中結果を確認し、別コマンドで次ケースを実行。WT1 v8.2は含まない。

## 重要
- 固定 BioEmu commit `42c551133f9825f4ec01812075b0d50e268c2209`、clean tree必須。入力PDB各494原子・63残基・SHA-256固定。HPacker再推論なし。
- `run_one_md(..., only_energy_minimization=False, simtime_ns=0)`で、初期安定化＋0.1ns NVT＋0.4ns NPTを固定ソースのまま実行。MD後重原子PDBを保存。
- OpenMM Simulation APIを実行時だけ観測するラッパーを挿入。Simulation.stepを1000刻みに分割しチェックポイントを取るため、**数値軌跡・性能の厳密同一性は保証しない**。力場・積分法・時間刻み・拘束・温度・圧力条件は変更しない。
- エネルギーは系全体の瞬時値であり、異なる溶媒系の間でエネルギーの大小比較をしない。
- 時間上限2400秒は計算を途中停止する。TIMEOUT時の最後の完了ステップ・PEはevents.jsonlに残る。
- 目的は完走・NaN/Inf・途中段階・最終PEの保存。自由エネルギー・介入の優劣は評価しない。

## 手順（PowerShell）
展開先：`C:\Users\user\NTBE_X4\X4_PhaseC_4Cases_Sequential_2026-10-23`

```powershell
cd C:\Users\user\NTBE_X4\X4_PhaseC_4Cases_Sequential_2026-10-23
& "..\.venv_md\Scripts\python.exe" .\run_phase_c.py --bioemu "..\bioemu" --case wt0_baseline
& "..\.venv_md\Scripts\python.exe" .\run_phase_c.py --bioemu "..\bioemu" --case wt0_baseline --execute
```
次ケースも同じ形で `wt0_v8_1_rescue`、`wt1_baseline`、`wt1_v8_1_rescue` の順に**1件ずつ**実行。自動一括実行しない。

各結果 `outputs_phase_c\<case>_<timestamp>` をZIP化してチャットに添付。run summaryでPASSでも、ここで検証した後に次ケースへ進む。
