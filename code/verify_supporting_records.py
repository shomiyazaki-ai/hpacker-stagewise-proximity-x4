#!/usr/bin/env python3
"""Recalculate saved E2 arrays, E3 checkpoints and E4 pairwise arithmetic only.
Does not run a trained model, construct an OpenMM system, minimize, or integrate.
"""
from pathlib import Path
import argparse,csv,hashlib,json,re,sys
import numpy as np
from verify_e1_from_pdb import parse_pdb,geometry,write_csv

def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def circ(a,b):return float(np.max(np.abs((np.asarray(a)-np.asarray(b)+180)%360-180)))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();R=a.root;O=a.out;O.mkdir(parents=True,exist_ok=True)
 # E2: direct numerical differences in the ten archived arrays.
 paths=[]
 for v in ['v83','v84']:
  for p in (R/'data/e2'/v).rglob('model_arrays.npz'):
   z=np.load(p,allow_pickle=False);q=p.with_name('result.json');rec=json.loads(q.read_text());old=rec['blocked_vs_restored']
   dg=z['zernike_blocked'].astype(float)-z['zernike_restored'].astype(float)
   row=dict(source=int(rec['case'][-3:]),stage=rec['source_stage'],donor=rec['donor'],receiver=rec['receiver'],zernike_l2=float(np.linalg.norm(dg)),raw_max=float(np.max(np.abs(z['raw_blocked']-z['raw_restored']))),chi_max_deg=circ(z['chi_blocked'],z['chi_restored']),array_file=str(p.relative_to(R)),result_file=str(q.relative_to(R)),tier='Saved-array recomputation; no new forward')
   assert abs(row['zernike_l2']-old['zernike_l2'])<1e-7 and abs(row['raw_max']-old['raw_yhat_max_abs'])<1e-7 and abs(row['chi_max_deg']-old['chi_max_circular_deg'])<1e-6
   assert np.array_equal(z['raw_restored'],z['raw_rescue']) and np.array_equal(z['zernike_restored'],z['zernike_rescue']) and circ(z['chi_restored'],z['chi_rescue'])==0
   paths.append(row)
 paths.sort(key=lambda x:(x['source'],x['receiver']));assert len(paths)==10
 write_csv(O/'e2_paths.csv',paths)
 # Receipt hashes and pre-Cartesian chi edits, limited to the two WT cases.
 receipt_files=0;edits=[];blocks=[]
 for v in ['v81','v82']:
  for p in (R/'data/e2'/v).rglob('receipt.json'):
   rec=json.loads(p.read_text())
   if not rec['case'].startswith('WT_'):continue
   for name,h in rec['files_sha256'].items():assert sh(p.parent/name.replace('\\','/'))==h;receipt_files+=1
   if 'mediator_edit' in rec:
    ed=rec['mediator_edit'];ss=ed['step'];res=ed['residue'];chi=json.loads((p.parent/f'stages/stage_{ss:02d}_chis.json').read_text());assert circ(ed['forced_chi'],chi[str(res)])<1e-6
    assert ed['injection_timing']=='BEFORE_SIDECHAIN_RECONSTRUCTION';edits.append(dict(source=int(rec['case'][-3:]),kind=rec['kind'],stage=ss,residue=res,receipt=str(p.relative_to(R))))
   else:
    obs=json.loads(p.with_name('observations.json').read_text());ss=obs['block_step'];res=rec['block_residues'];blocks.append(dict(source=int(rec['case'][-3:]),stage=ss,residues='+'.join(map(str,res)),operation='Rescue background; pre-Cartesian chi from original branch',receipt=str(p.relative_to(R))))
    for ed in obs['edits']:
     chi=json.loads((p.parent/f"stages/stage_{ed['step']:02d}_chis.json").read_text());assert circ(ed['forced_chi'],chi[str(ed['residue'])])<1e-6
 blocks.sort(key=lambda x:(x['source'],x['residues']));assert len(blocks)==5
 write_csv(O/'e2_edits.csv',edits);write_csv(O/'e2_blocks.csv',blocks)
 # One historical selected full-batch receipt: inspect only, no execution.
 fb=next((R/'data/e2/selected_full_batch').rglob('original_result.json'));fbr=json.loads(fb.read_text());fp=fb.parent
 assert sh(fp/'blocked_stage04.pdb')==fbr['source_blocked_pdb_sha256'];assert sh(fp/'rescue_stage04.pdb')==fbr['source_rescue_pdb_sha256'];assert fbr['n_batch']==54 and fbr['target_batch_index']==46
 # E3: original per-run log continuity, metadata, hashes and geometry.
 runs=[]
 for envp in (R/'data/e3/phase_c/outputs_phase_c').glob('*/environment.json'):
  d=envp.parent;env=json.loads(envp.read_text());case=env['case'];inp=R/'data/e3/phase_c'/(case+'.pdb');assert sh(inp)==env['input_sha256']==sh(R/'data/e4/inputs'/(case+'.pdb'))
  hm=json.loads((d/'sha256_manifest.json').read_text());assert all(sh(d/k)==v for k,v in hm.items())
  events=[json.loads(l) for l in (d/'events.jsonl').read_text().splitlines()];begin=[e for e in events if e['event']=='STEP_CALL_BEGIN'];end=[e for e in events if e['event']=='STEP_CALL_END'];assert len(begin)==len(end)==204
  total=0
  for b,e in zip(begin,end):assert b['total_before']==total;total+=b['requested_steps'];assert total==e['total_completed']
  assert total==611100
  energy=[e for e in events if e['event']=='ENERGY_STATE'];assert len(energy)==226 and all(np.isfinite(e['potential_kj_mol']) and np.isfinite(e['kinetic_kj_mol']) and not any(e['nonfinite'].values()) for e in energy)
  prep=next(e for e in events if e['event']=='PREPARE_END');platform=next(e for e in events if e['event']=='MIN_BEGIN');complete=next(e for e in events if e['event']=='COMPLETE');summary=json.loads((d/'summary.json').read_text());assert summary['exit_code']==0
  initial=geometry(inp)[0];final=geometry(d/'final_heavy_atoms.pdb')[0];assert final['general_pairs_lt_1p6']==0
  runs.append(dict(case=case,source=int(case[2]),structure='Original' if case.endswith('baseline') else 'Locally perturbed',initial_pairs=initial['general_pairs_lt_1p6'],final_pairs=final['general_pairs_lt_1p6'],particles=prep['n_particles'],waters=prep['water_residues'],steps=total,energy_states=len(energy),platform=platform['platform'],precision=platform['platform_properties']['Precision'],final_PE=complete['final_potential_kj_mol'],completion='Completed',input_sha256=sh(inp),events_file=str((d/'events.jsonl').relative_to(R))))
 runs.sort(key=lambda x:(x['source'],x['structure']!='Original'));assert [r['initial_pairs'] for r in runs]==[6,1,3,1];write_csv(O/'e3_runs.csv',runs)
 m2=[]
 for p in (R/'data/e3/separate_m2').rglob('*_minimized_heavy.pdb'):
  g=geometry(p)[0];assert g['general_pairs_lt_1p6']==0;m2.append(dict(file=str(p.relative_to(R)),general_pairs_lt_1p6=0,run='Separate WT0 M2, not Phase C checkpoint'))
 assert len(m2)==2;write_csv(O/'e3_separate_m2.csv',m2)
 # E4: all original pair distances and energies. Parameters come from archived pair CSV,
 # whose upstream force-field typing was separately checked in Phase 1A (not rerun here).
 summary=json.loads((R/'data/e4/results/stage_a_summary.json').read_text());lj=[];dfs={}
 for r in runs:
  case=r['case'];p=R/'data/e4/inputs'/(case+'.pdb');keys,xyz=parse_pdb(p);coords={(k[1],k[-1]):v for k,v in zip(keys,xyz)}
  rows=list(csv.DictReader((R/'data/e4/results'/(case+'_pairs.csv')).open()));assert len(rows)==116140 and len({z['pair_id'] for z in rows})==116140
  ds=np.array([np.linalg.norm(coords[int(z['res_i']),z['atom_i']]-coords[int(z['res_j']),z['atom_j']]) for z in rows]);eps=np.array([float(z['epsilon_kj_mol']) for z in rows]);sig=np.array([float(z['sigma_nm']) for z in rows]);saved_dist=np.array([float(z['distance_A']) for z in rows]);saved=np.array([float(z['lj_net_kj_mol']) for z in rows]);formula=4*eps*((sig/(saved_dist/10))**12-(sig/(saved_dist/10))**6);double=4*eps*((sig/(ds/10))**12-(sig/(ds/10))**6)
  assert np.max(np.abs(formula-saved)/np.maximum(1,np.abs(saved)))<1e-12
  assert all(abs(int(z['res_i'])-int(z['res_j']))>1 and 'OXT' not in z['pair_id'] for z in rows)
  # Reproduce the original nm-float32 quantization separately.
  cc={k:(v/10).astype(np.float32).astype(float) for k,v in coords.items()};fd=np.array([np.linalg.norm(cc[int(z['res_i']),z['atom_i']]-cc[int(z['res_j']),z['atom_j']]) for z in rows]);assert np.max(abs(fd*10-saved_dist))<1e-12
  old=summary['cases'][case]['lj_net_kj_mol'];pair_sum=float(saved.sum());assert abs(pair_sum-old)/old<1e-10
  lj.append(dict(case=case,source=r['source'],structure=r['structure'],pairs=len(rows),archived_summary_kj_mol=old,archived_pair_sum_kj_mol=pair_sum,double_PDB_kj_mol=float(double.sum()),pairs_lt1p6=int(sum(ds<1.6)),max_saved_vs_double_distance_A=float(max(abs(ds-saved_dist))),max_formula_scaled_error=float(max(abs(formula-saved)/np.maximum(1,abs(saved))))));dfs[case]={z['pair_id']:z for z in rows}
 tracked=[]
 for s in [0,1]:
  akey=f'wt{s}_baseline';bkey=f'wt{s}_v8_1_rescue';top=summary['paired_differences'][f'wt{s}']['top_20_absolute_pair_contributions'][0];key=top['pair_id'];a0=dfs[akey][key];b0=dfs[bkey][key];reduction=summary['cases'][akey]['lj_net_kj_mol']-summary['cases'][bkey]['lj_net_kj_mol'];contribution=float(a0['lj_net_kj_mol'])-float(b0['lj_net_kj_mol']);tracked.append(dict(source=s,pair=key,original_distance_A=float(a0['distance_A']),perturbed_same_pair_distance_A=float(b0['distance_A']),percent_proxy_reduction=100*contribution/reduction))
 write_csv(O/'e4_lj.csv',lj);write_csv(O/'e4_tracked_pairs.csv',tracked)
 result=dict(status='PASS',e2_paths=len(paths),e2_block_cases=len(blocks),e2_receipt_files_hashed=receipt_files,e2_forward_executed=False,e3_runs=len(runs),e3_input_runtime_environment_sha_matches=4,e3_final_zero_pairs=4,e3_separate_m2_structures=2,e4_pairs_per_structure=116140,e4_pairwise_total=464560,e4_new_forcefield_typing=False,md_executed=False)
 (O/'supporting_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
