from pathlib import Path
import csv,json
import sys
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'code'))
from verify_e1_from_pdb import geometry,write_csv
V=R/'verification/current';D=R/'data/authoritative';D.mkdir(exist_ok=True)
def read(n):return list(csv.DictReader((V/n).open()))
e1=json.loads((V/'e1_verification.json').read_text());stages=read('e1_stages.csv');paths=read('e2_paths.csv');runs=read('e3_runs.csv');lj=read('e4_lj.csv');tracked=read('e4_tracked_pairs.csv')
tables={}
def table(id,headers,rows,caption,footnote=''):
 tables[id]=dict(headers=headers,rows=rows,caption=caption,footnote=footnote)
 with (D/(id+'.csv')).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(headers);w.writerows(rows)
table('Table_1',['Endpoint','Definition across stages 0-5','Frames / 200'],[
 ['General union','Any qualifying atom pair in at least one stage','196'],
 ['Strict same-atom recontact','The same unordered atom pair: present, absent, present at i < j < k','147'],
 ['PHE6-PHE55 union','Residue pair close in at least one stage','67'],
 ['PHE6-PHE55 recontact','Residue-pair binary sequence: contact, absence, contact','43'],
 ['Complete PHE alternation','010101 (10) or 101010 (15)','25']],
 'Table 1. Fixed-ensemble sequence endpoints.', 'All denominators are the same 200 source frames. Strict atom-pair recontact and residue-pair recontact are different endpoints. The historical 153 counter is superseded for strict recontact.')
g=[]
for case in ['wt0_baseline','wt0_v8_1_rescue']:
 m=geometry(R/'data/e4/inputs'/(case+'.pdb'))[0];g.append(m)
table('Table_2',['Metric (source 0, stage 5)','Original','Locally perturbed'],[
 ['PHE6-PHE55 minimum distance (Å)',*[f"{r['phe655_min_distance_A']:.6f}" for r in g]],
 ['PHE6-PHE55 atom pairs <1.6 Å',*[str(r['phe655_n_pairs_lt_1p6']) for r in g]],
 ['General atom pairs <1.6 Å',*[str(r['general_pairs_lt_1p6']) for r in g]]],
 'Table 2. Selected pre-relaxation local geometry.', 'The perturbed structure still has one general close pair (minimum 1.550220 Å); separating PHE6-PHE55 does not remove every extreme proximity.')
strong=[p for p in paths if (int(p['source']),int(p['receiver'])) in [(0,55),(1,3),(1,56),(1,59)]]
table('Table_3',['Source','Stage','Donor → receiver','Max raw-output |Δ|','Max circular χ |Δ| (°)'],[[p['source'],p['stage'],p['donor']+' → '+p['receiver'],f"{float(p['raw_max']):.6f}",f"{float(p['chi_max_deg']):.6f}"] for p in strong],
 'Table 3. Four selected stronger NN-response examples.', 'Raw eight-component output differences are dimensionless, not probabilities. Selection is post hoc; circular χ differences have no chemical symmetry correction. See Table S2B for all ten saved-array comparisons and the separate selected full-batch receipt.')
table('Table_4',['Source','Structure','Initial pairs\n<1.6 Å','Final MD pairs\n<1.6 Å','Completion'],[[r['source'],r['structure'],r['initial_pairs'],r['final_pairs'],r['completion']] for r in runs],
 'Table 4. Initial inputs and final short-MD structures.', 'These are two matched timepoints per selected Phase C case. No minimization-only checkpoint is inserted. The separate WT0 M2 results appear in S5.4; WT1 minimized geometry is unavailable.')
tr={int(p['source']):p for p in tracked}
table('Table_5',['Source','Structure','Archived LJ proxy\n(kJ/mol)','Tracked baseline-dominant pair','Same-pair\ndistance (Å)'],[[r['source'],r['structure'],f"+{float(r['archived_summary_kj_mol']):.6e}" if r['case']=='wt0_baseline' else f"+{float(r['archived_summary_kj_mol']):,.2f}",tr[int(r['source'])]['pair'].replace('|',' - '),f"{float(tr[int(r['source'])]['original_distance_A' if r['structure']=='Original' else 'perturbed_same_pair_distance_A']):.3f}"] for r in lj],
 'Table 5. Restricted protein heavy-atom LJ proxy and tracked pairs.', 'Values are archived summary outputs. The perturbed distance tracks the baseline-dominant atom pair, which need not remain the nearest pair. Coordinate precision comparisons are in Table S4. These are neither solvated potential energies nor free energies.')
table('Table_S1',['Stage / endpoint','General frames / 200','PHE6-PHE55 frames / 200'],[[str(s),str(e1['general_stage_counts'][s]),str(e1['phe_stage_counts'][s])] for s in range(6)]+[['Union','196','67'],['Ordered recontact','147 (strict same atom)','43 (residue pair)'],['Full alternation','Not defined','25 = 10 + 15']],
 'Table S1. Complete fixed-ensemble counts.', 'Stage 0 is the initial HPacker prediction; stages 1-5 are refinement updates, not physical time. Contact requires a distance strictly below 1.6 Å.')
plan=json.loads(next((R/'data/e2/v81').rglob('v8_plan.json')).read_text())
table('Table_S2A_1',['Source','Pulse stage','Operation','Target residue','Perturbed residue','Measured partner','χ edit stage / residue'],[[x['case'][-1],str(x['pulse']['step']),x['pulse']['mode'].replace('_',' '),str(x['pulse']['target']),str(x['pulse']['intervened']),str(x['pulse']['measured_partner']),str(x['step'])+' / '+str(x['mediator'])] for x in plan['cases'] if x['case'].startswith('WT_')],
 'Table S2A-1. Original local-input pulses and subsequent pre-Cartesian χ interventions.', 'The source-0 third_lag branch targets residue 6, perturbs residue 55 and measures partner 58. χ rescue uses baseline χ on the pulse background; the paired archived sufficiency probe uses branch χ on the original background. These operational names do not establish unique causal necessity or sufficiency.')
blocks=read('e2_blocks.csv')
table('Table_S2A_2',['Source','Block stage','Blocked residue(s)','Operation'],[[x['source'],x['stage'],x['residues'],'Original-branch χ substituted before reconstruction on rescue background'] for x in blocks],
 'Table S2A-2. Five downstream χ blocks.', 'The preceding rescue χ edit is retained: source 0 stage 3 residue 17; source 1 stage 2 residue 56. The blocks are path-specific selected perturbations.')
for p in paths:p['path_key']=f"N{int(p['source'])}{int(p['receiver']):02d}"
table('Table_S2B',['Key','Source / stage','Donor → receiver','Zernikegram L2','Raw max |Δ|','Circular χ max |Δ| (°)','Tier'],[[p['path_key'],p['source']+' / '+p['stage'],p['donor']+' → '+p['receiver'],f"{float(p['zernike_l2']):.6f}",f"{float(p['raw_max']):.6f}",f"{float(p['chi_max_deg']):.6f}",('A; R' if p['path_key']=='N055' else 'A')] for p in paths],
 'Table S2B. Ten selected saved-array NN comparisons.', 'A = P1 independent array arithmetic, no trained forward. R = an archived one-point full-batch receipt exists only for N055 (source 0, stage 4, donor 53 → receiver 55); P1 inspected it without rerunning it. Exact array/result paths are in NN_File_Crosswalk.csv. Strong/weak selection is retrospective.')
write_csv(D/'NN_File_Crosswalk.csv',paths)
table('Table_S3',['Case','Particles','Waters','Steps','Energy states','Final PE (kJ/mol)','Platform / precision'],[[r['case'].replace('wt','s').replace('_baseline',' original').replace('_v8_1_rescue',' perturbed'),r['particles'],r['waters'],r['steps'],r['energy_states'],f"{float(r['final_PE']):.6f}",r['platform']+' / '+r['precision']] for r in runs],
 'Table S3. Four archived Phase C continuation records.', 'Each case has exit code 0, 204 matched begin/end calls, one COMPLETE record and no recorded nonfinite energy/position/velocity/force flags. Integration duration is 0.5004 ns, including four 0.1-ps stabilization stages. Solvent counts differ; final PE values cannot be ranked as stability.')
table('Table_S4A',['Source / structure','Archived summary (kJ/mol)','Archived CSV sum (kJ/mol)','PDB-double recomputation (kJ/mol)'],[[r['source']+' / '+r['structure'],f"{float(r['archived_summary_kj_mol']):.6f}",f"{float(r['archived_pair_sum_kj_mol']):.6f}",f"{float(r['double_PDB_kj_mol']):.6f}"] for r in lj],
 'Table S4A. LJ coordinate-pathway and accumulation comparisons.', 'Saved summary and summation of printed pair CSV values differ slightly (source-0 original: about 10.25 kJ/mol, relative 6.17×10⁻¹³). Direct PDB-double coordinates differ from the original nm-float32 pathway. Values are recorded separately; none replaces the archived scientific result.')
table('Table_S4B',['Source','Baseline-dominant pair','Original distance (Å)','Perturbed same-pair distance (Å)','Share of proxy reduction (%)'],[[r['source'],r['pair'].replace('|',' - '),f"{float(r['original_distance_A']):.6f}",f"{float(r['perturbed_same_pair_distance_A']):.6f}",f"{float(r['percent_proxy_reduction']):.7f}"] for r in tracked],
 'Table S4B. Same-pair distance and contribution tracking.', 'Shares use the archived pair contribution and archived summary reduction. The measured pair is fixed by its baseline dominance, not by the nearest pair in the perturbed structure.')
claims=[['E1 stage counts / recontact','1,200 original stage PDBs + receipts','data/pdb; provenance/PDB_SHA256_Manifest.csv','P1 independent geometry and SHA'],['Source-0 local geometry','Original and perturbed stage-5 PDBs','data/e4/inputs; verification/current/e1_stages.csv','P1 independent geometry'],['Ten NN-response comparisons','Saved NPZ vectors + result JSON','data/e2/v83; data/e2/v84; NN_File_Crosswalk.csv','P1 saved-array arithmetic'],['One selected full-batch diagnostic','54-residue batch receipt, arrays, two PDBs','data/e2/selected_full_batch','Archived forward receipt inspected; no P1 forward'],['Selected continuation','Four run environments, logs, summaries and final PDBs','data/e3/phase_c/outputs_phase_c','P1 log, SHA and final geometry'],['Separate WT0 minimized geometry','Two M2 minimized PDBs','data/e3/separate_m2','P1 geometry; separate run'],['Restricted LJ proxy','4×116,140 pair rows + selected input PDBs','data/e4/results; data/e4/inputs','P1 pairwise arithmetic / PDB doubles'],['Force-field typing','Archived code and Phase1A parameter audit','code/archived_execution; provenance/phase1a_historical','Historical independent typing; no new typing']]
table('Table_S5',['Claim','Evidence','Package file / directory','Verification tier'],claims,'Table S5. Claim-evidence-file-verification map.', 'Full file hashes are in SHA256SUMS.txt. Read-only verifiers are code/verify_e1_from_pdb.py and code/verify_supporting_records.py. No universal NN, HPacker or MD replay is claimed.')
# Compact source-0 figure manifest directly from the verified geometry table.
f0=[r for r in stages if r['frame']=='0'];write_csv(D/'Figure_2_Source_0.csv',f0)
(D/'tables.json').write_text(json.dumps(tables,indent=2,ensure_ascii=False)+'\n')
print('Frozen',len(tables),'display tables; figure-source and exact-file crosswalk saved')
