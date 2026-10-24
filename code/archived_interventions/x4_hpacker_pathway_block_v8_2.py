#!/usr/bin/env python3
"""NT-BE-X4 full-network five-refinement, one-shot-input-perturbation forensic test.
Only existing 8 ineligible diagnostic frames. No OpenMM, MD, acceptance or publication.
"""
import argparse, hashlib, json, subprocess, sys, tempfile, traceback, zipfile
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.spatial.distance import pdist, squareform

COMMIT='6cd5ec33961bb2662b51c0a15d11d6fd9162f2a7'
HPARAMS={'initial_guess':'cded11d51b3aeb3d1b2107707959f7e9b695a556dd649e326068d84dd07caed5','refinement':'888dfc1df9e8a1b5e83a392197cb71a0ad0e4f14b865524c202120857f844f6a','initial_guess_conditioned':'d6cc54a8c45289b26c2d39017522d4ff75d8fc60c9da29538dca8222734f012f'}
PAIRS={'WT_source000':[(6,55),(6,58),(22,58)],'WT_source001':[(6,55),(6,59)],'WT_source002':[(6,58),(16,52),(23,39)],'WT_source003':[(23,39)],'mutant_source000':[(22,58),(6,22)],'mutant_source001':[(6,55),(6,59)],'mutant_source002':[(53,58),(4,7),(1,58)],'mutant_source003':[(23,35),(48,52)]}
BB={'N','CA','C','O','OXT'}

def require(ok,msg):
 if not ok:raise RuntimeError('STOP: '+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def put_json(p,obj):
 p=Path(p);temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8');temp.replace(p)
def dec(v):return v.decode().strip() if isinstance(v,bytes) else str(v).strip()
def circ(a,b):return abs((float(a)-float(b)+180)%360-180)

def pdb_atoms(b):
 d={}
 for line in b.decode('ascii').splitlines():
  if not line.startswith(('ATOM  ','HETATM')):continue
  key=(line[21:22],int(line[22:26]),line[26:27],line[12:16].strip())
  require(key not in d,'duplicate atom '+str(key))
  coord=np.array([float(line[i:i+8]) for i in (30,38,46)])
  element=line[76:78].strip() or key[-1][:1]
  d[key]=(coord,element)
 return d

def compare_pdb(a,b,tol=.05):
 aa=pdb_atoms(a);bb=pdb_atoms(b)
 if set(aa)!=set(bb):return {'pass':False,'missing':len(set(bb)-set(aa)),'extra':len(set(aa)-set(bb))}
 k=max(aa,key=lambda k:float(np.linalg.norm(aa[k][0]-bb[k][0])))
 md=float(np.linalg.norm(aa[k][0]-bb[k][0]))
 return {'pass':md<=tol,'max_A':md,'worst_key':str(k),'tol_A':tol}

def parse_trace(zips):
 traces={};inputs={}
 for zp in zips:
  with zipfile.ZipFile(zp) as z:
   require(z.testzip() is None,'trace ZIP invalid')
   for name in z.namelist():
    if not name.endswith('/trace/trace_manifest.json'):continue
    case=name.split('/')[-3];require(case not in traces,'duplicate case')
    meta=json.loads(z.read(name));require(meta['commit']==COMMIT,'trace git commit')
    p=name[:-len('trace/trace_manifest.json')]
    stage=[]
    for k in range(6):
     n=p+'trace/'+('stage_00_initial.pdb' if k==0 else f'stage_{k:02d}_refinement.pdb')
     b=z.read(n);require(sha(b)==meta['trace'][k]['sha256'],'trace stage hash')
     stage.append(b)
    require(stage[-1]==z.read(p+'trace/official_final.pdb'),'final PDB not official')
    traces[case]=stage;inputs[case]=z.read(p+'official_backbone_input.pdb')
 require(set(traces)==set(PAIRS),'8-case fixed coverage mismatch')
 return traces,inputs

def load_v2(path):
 records={}
 with zipfile.ZipFile(path) as z:
  require(z.testzip() is None,'v2 zip corrupt')
  for name in z.namelist():
   if not name.replace('\\','/').endswith('/observations.json'):continue
   x=json.loads(z.read(name));key=(x['case'],x['input_stage'])
   require(key not in records and not x['failures'],'v2 key duplicates or errors')
   records[key]=x['records']
 require(len(records)==40,'v2 40 source states not found')
 return records

def find(v2,case,stage,target,partner,role,variant):
 matches=[x for x in v2[case,stage] if x['target']==target and x['partner']==partner and x['role']==role and x['variant']==variant]
 require(len(matches)<=1,'duplicate v2 record');return matches[0] if matches else None

def scenario(case,step,target,intervened,measured,kind,mode,v2record):
 return {'case':case,'step':step,'target':target,'intervened':intervened,'measured_partner':measured,
         'kind':kind,'mode':mode,'v2_chi':[float(x) for x in v2record['chi_deg']],
         'v2_yhat':[float(x) for x in v2record['raw_yhat']]}

def plan(v2,phase):
 out=[];cases=sorted(PAIRS) if phase=='all' else ['WT_source000','mutant_source001']
 for case in cases:
  for a,b in PAIRS[case]:
   if phase=='pilot' and (a,b)!=(6,55):continue
   for target,partner in [(a,b),(b,a)]:
    for step in (range(1,6) if phase=='all' else [2,4,5]):
     state=step-1
     for mode in ('partner_mask','partner_lag'):
      if mode.endswith('lag') and step==1:continue
      rec=find(v2,case,state,target,partner,'pair',mode)
      require(rec is not None,'missing primary v2 record')
      out.append(scenario(case,step,target,partner,partner,'pair',mode,rec))
     baseline=find(v2,case,state,target,partner,'pair','baseline')
     require(baseline is not None,'missing v2 baseline')
     thirds=[]
     for rec in v2[case,state]:
      if rec['target']==target and rec['role']==f'third_for_{partner}' and rec['variant']=='third_mask':
       strength=max(circ(u,v) for u,v in zip(baseline['chi_deg'],rec['chi_deg']))
       thirds.append((strength,int(rec['partner'])))
     for strength,third in sorted(thirds,reverse=True)[:(1 if phase=='pilot' else 2)]:
      for mode in ('third_mask','third_lag'):
       if mode.endswith('lag') and step==1:continue
       rr=find(v2,case,state,target,third,f'third_for_{partner}',mode)
       require(rr is not None,'missing tested third v2 record')
       out.append(scenario(case,step,target,third,partner,'third',mode,rr))
 # Placebo intervention: intercept exact native batch, pass identical neighborhood unchanged.
 for case in cases:
  a,b=PAIRS[case][0]
  rec=find(v2,case,2,a,b,'pair','baseline')
  require(rec is not None,'missing placebo v2 baseline')
  out.append(scenario(case,3,a,b,b,'placebo','noop',rec))
 out.sort(key=lambda v:(v['case'],v['step'],v['target'],v['measured_partner'],v['intervened'],v['mode']))
 require(len(out)==len(set(task_id(x) for x in out)),'duplicate task IDs')
 return out

def task_id(x):
 return f"{x['case']}_s{x['step']}_t{x['target']}_i{x['intervened']}_p{x['measured_partner']}_{x['kind']}_{x['mode']}"

def metrics(b,pairs):
 d=pdb_atoms(b)
 keys=[k for k in d if d[k][1].upper() not in ('H','D')]
 xyz=np.asarray([d[k][0] for k in keys]);r=np.asarray([k[1] for k in keys]);dist=squareform(pdist(xyz))
 nonlocal_mask=np.triu(np.abs(r[:,None]-r[None,:])>1,1)
 selected=dist[nonlocal_mask];which=np.argwhere(nonlocal_mask)[int(np.argmin(selected))]
 measures={}
 for a,b in pairs:
  ia=np.where(r==a)[0];ib=np.where(r==b)[0]
  require(len(ia)>0 and len(ib)>0,'missing pair residue')
  dd=dist[np.ix_(ia,ib)];idx=np.unravel_index(int(np.argmin(dd)),dd.shape)
  measures[f'{a}-{b}']={'min_A':float(dd[idx]),'pairs_below_1p6A':int(np.sum(dd<1.6)),
                         'closest_atoms':[keys[int(ia[idx[0]])][-1],keys[int(ib[idx[1]])][-1]]}
 return {'nonlocal_atom_pairs_below_1p6A':int(np.sum(selected<1.6)),
         'global_min_nonlocal_A':float(np.min(selected)),
         'global_closest_atoms':[[keys[int(i)][1],keys[int(i)][-1]] for i in which],
         'watched_pairs':measures}

def propagation(a,b,target):
 aa=pdb_atoms(a);bb=pdb_atoms(b)
 require(set(aa)==set(bb),'atom keys changed between branch and control')
 result={}
 for resid in sorted(set(k[1] for k in aa)):
  keys=[k for k in aa if k[1]==resid and k[-1] not in BB and aa[k][1].upper() not in ('H','D')]
  if not keys:continue
  ds=[float(np.linalg.norm(aa[k][0]-bb[k][0])) for k in keys]
  result[str(resid)]={'max_displacement_A':max(ds),'rms_displacement_A':float(np.sqrt(np.mean(np.square(ds))))}
 other=[(n,v['max_displacement_A']) for n,v in result.items() if int(n)!=target]
 return {'n_other_res_gt_0p1A':sum(d>.1 for _,d in other),'n_other_res_gt_1A':sum(d>1 for _,d in other),
         'top_10_other_residues':sorted(other,key=lambda t:-t[1])[:10], 'per_residue':result}

def network_chi(branch,base,step):
 a=json.loads((branch/'stages'/f'stage_{step:02d}_chis.json').read_text())
 b=json.loads((base/'stages'/f'stage_{step:02d}_chis.json').read_text())
 require(set(a)==set(b),'native prediction residue IDs changed')
 ds={}
 for rid in a:
  require(len(a[rid])==len(b[rid]),'native prediction count changed')
  if a[rid]:ds[rid]=max(circ(u,v) for u,v in zip(a[rid],b[rid]))
 return {'n_res_gt_0p5deg':sum(x>.5 for x in ds.values()),
         'n_res_gt_30deg':sum(x>30 for x in ds.values()),
         'top_10':sorted(ds.items(),key=lambda x:-x[1])[:10],'all':ds}

def env(root):
 require((root/'hpacker'/'.git').exists(),'missing HPacker Git checkout')
 actual=subprocess.check_output(['git','-C',str(root/'hpacker'),'rev-parse','HEAD'],text=True).strip()
 require(actual==COMMIT,'wrong HPacker git SHA')
 require(Path(sys.executable).resolve()==(root/'.venv'/'Scripts'/'python.exe').resolve(),'not frozen .venv Python')
 import importlib
 loaded=Path(importlib.import_module('hpacker.hpacker').__file__).resolve()
 source=root/'hpacker'/'hpacker'/'hpacker.py'
 require(sha(loaded.read_bytes())==sha(source.read_bytes()),'loaded HPacker differs from repo source')
 hashes={'commit':COMMIT,'python':str(Path(sys.executable).resolve()),'loaded_hpacker':str(loaded),
         'hpacker_py_sha256':sha(source.read_bytes()),'hparams':{},'weights':{}}
 for model,digest in HPARAMS.items():
  loc=root/'hpacker'/'hpacker'/'pretrained_models'/model
  inst=loaded.parent/'pretrained_models'/model
  for fname,expected in [('hparams.json',digest)]:
   require(sha((loc/fname).read_bytes())==expected,'config SHA '+model)
   require(sha((inst/fname).read_bytes())==expected,'loaded config SHA '+model)
  hashes['hparams'][model]=digest
  modelname='lowest_valid_loss_model.pt'
  w=sha((loc/modelname).read_bytes())
  require(sha((inst/modelname).read_bytes())==w,'weights SHA '+model)
  hashes['weights'][model]=w
 return hashes

class FullReplay:
 def __init__(self,case,backbone,intervention,out):
  self.case=case;self.backbone=backbone;self.intervention=intervention;self.out=out
  self.step=-1;self.event_count=0;self.event=None;self.states={};self.mediator_edit=None;self.mediator_edit_log=None;self.mediator_edit_logs=[]
 def modified_nb(self,h,native,r,mi):
  sc=self.intervention;idx=sc['intervened'];kind=sc['mode'];note={'mode':kind,'intervened_residue':idx}
  if kind=='noop':
   nb,bb=native(r,mi);note['no_change']=True
  elif kind.endswith('lag'):
   require(self.step>=2,'lag needs earlier constructed stage')
   past=pdb_atoms(self.states[self.step-2])
   old={k[3]:val[0] for k,val in past.items() if k[1]==idx and k[3] not in BB}
   residues=[rr for rr in h.structure.get_residues() if rr.id[1]==idx]
   require(len(residues)==1,'partner residue not unique')
   atoms={a.name:a for a in residues[0] if a.name not in BB}
   require(set(old)==set(atoms),'lag atom lists differ')
   save={n:a.coord.copy() for n,a in atoms.items()}
   try:
    for n,a in atoms.items():a.coord=np.asarray(old[n],dtype=a.coord.dtype)
    nb,bb=native(r,mi)
   finally:
    for n,a in atoms.items():a.coord=save[n]
   note['moved_sidechain_atoms']=len(atoms)
  elif kind.endswith('mask'):
   nb,bb=native(r,mi)
   mask=np.asarray([(int(dec(rr[4]))==idx and dec(an) not in BB) for rr,an in zip(nb['res_ids'],nb['atom_names'])])
   for key in ('res_ids','atom_names','elements','coords','charges'):
    nb[key]=nb[key][~mask].copy()
   note['removed_sidechain_atoms']=int(sum(mask))
  else:raise ValueError('unexpected intervention '+kind)
  note['neighborhood']=[{'resid':int(dec(rr[4])),'atom':dec(an),'element':dec(e),'xyz_CA_A':[float(q) for q in co],'charge':float(ch)} for rr,an,e,co,ch in zip(nb['res_ids'],nb['atom_names'],nb['elements'],nb['coords'],nb['charges'])]
  return nb,bb,note
 def execute(self):
  from hpacker import HPacker
  with tempfile.TemporaryDirectory(prefix='ntbe_x4_feedback_') as temp:
   inp=Path(temp)/'backbone.pdb';inp.write_bytes(self.backbone)
   hp=HPacker(str(inp))
   native_initial=hp.initial_guess;native_refine=hp.refinement
   native_batch=hp._get_predicted_chi_angles_batch
   native_neighborhood=hp._get_neighborhood
   native_zgram=hp._get_zernikegram
   stages=self.out/'stages';stages.mkdir()
   def record(step,angles):
    dest=stages/f'stage_{step:02d}.pdb';hp.write_pdb(str(dest));self.states[step]=dest.read_bytes()
    put_json(stages/f'stage_{step:02d}_chis.json',{str(int(k[1])):[float(v) for v in vals] for k,vals in angles.items()})
    print('  '+self.case+' stage '+str(step),flush=True)
   def initial_wrapper(*args,**kwargs):
    self.step=0;vals=native_initial(*args,**kwargs);record(0,vals);return vals
   def refine_wrapper(*args,**kwargs):
    self.step+=1;vals=native_refine(*args,**kwargs)
    record(self.step,vals);return vals
   def force_before_rebuild(res_ids, model_idx, angles):
    edits=self.mediator_edit or []
    if isinstance(edits,dict):edits=[edits]
    if model_idx!=1:return angles
    result=list(angles)
    for edit in edits:
     if self.step!=edit['step']:continue
     candidates=[i for i,r in enumerate(res_ids) if int(r[1])==edit['residue']]
     require(len(candidates)==1,'edit residue absent or ambiguous')
     i=candidates[0]
     previous=[float(v) for v in result[i]]
     forced=[float(v) for v in edit['chi']]
     require(len(previous)==len(forced) and len(forced)>0,'edit chi dimensional mismatch')
     require(not any(e['step']==self.step and e['residue']==edit['residue'] for e in self.mediator_edit_logs),'duplicate edit')
     result[i]=forced
     info={'step':self.step,'residue':edit['residue'],'native_predicted_chi':previous,'forced_chi':forced,
       'delta_from_native_deg':[circ(x,y) for x,y in zip(forced,previous)],
       'injection_timing':'BEFORE_SIDECHAIN_RECONSTRUCTION','prediction_index':i}
     self.mediator_edit_logs.append(info)
     self.mediator_edit_log=info
    return result
   def batch_wrapper(res_ids,model_idx):
    sc=self.intervention
    active=sc is not None and model_idx==1 and self.step==sc['step'] and any(int(x[1])==sc['target'] for x in res_ids)
    if not active:
     return force_before_rebuild(res_ids,model_idx,native_batch(res_ids,model_idx))
    self.event_count+=1;require(self.event_count==1,'repeated injection')
    index=next(i for i,r in enumerate(res_ids) if r[1]==sc['target'])
    record_nb={};tensors={};y=[]
    def patched_nb(res,mi):
     if mi==1 and res[1]==sc['target']:
      nb,bb,info=self.modified_nb(hp,native_neighborhood,res,mi)
      record_nb.update(info);return nb,bb
     return native_neighborhood(res,mi)
    def patched_zg(nb,mi):
     result=native_zgram(nb,mi)
     if mi==1 and int(dec(nb['res_id'][4]))==sc['target']:
      tensors['zernike']=result.detach().cpu().numpy()
     return result
    def hook(m,i,o):
     out=o[0] if isinstance(o,(tuple,list)) else o
     y.append(out.detach().cpu().numpy().astype(float))
    hp._get_neighborhood=patched_nb;hp._get_zernikegram=patched_zg
    handle=hp.model[1].register_forward_hook(hook)
    try:angles=native_batch(res_ids,model_idx)
    finally:
     handle.remove();hp._get_neighborhood=native_neighborhood;hp._get_zernikegram=native_zgram
    require(len(y)==1 and y[0].shape==(len(res_ids),8),'raw model output not Nx8')
    require('zernike' in tensors,'target Zernikegram not captured')
    chi=[float(v) for v in angles[index]];raw=y[0][index].tolist()
    require(len(chi)==len(sc['v2_chi']),'target χ count mismatch')
    ch=max(circ(a,b) for a,b in zip(chi,sc['v2_chi']))
    dy=max(abs(a-b) for a,b in zip(raw,sc['v2_yhat']))
    self.event={'chi_deg':chi,'raw_yhat':raw,'v2_chi_max_diff_deg':ch,'v2_yhat_max_abs_diff':dy,
      'v2_match':ch<=1 and dy<=.05,'batch_size':len(res_ids),'target_batch_index':index,**record_nb}
    np.savez_compressed(self.out/'intervention_tensor.npz',zernike=tensors['zernike'],raw_yhat=np.array(raw))
    return force_before_rebuild(res_ids,model_idx,angles)
   hp.initial_guess=initial_wrapper;hp.refinement=refine_wrapper
   hp._get_predicted_chi_angles_batch=batch_wrapper
   hp.reconstruct_sidechains(num_refinement_iterations=5)
   require(set(self.states)==set(range(6)),'incomplete six structural outputs')
   if self.intervention is not None:require(self.event_count==1,'input intervention never applied')
  return {'sha256':{str(k):sha(v) for k,v in self.states.items()},'event':self.event,'n_events':self.event_count}

def baseline(args,traces,backbone):
 dest=args.output/'baseline';dest.mkdir(exist_ok=True)
 for case in sorted(traces):
  target=dest/case;receipt=target/'receipt.json'
  if receipt.is_file():
   r=json.loads(receipt.read_text());require(r['status']=='PASS','nonpass prior baseline')
   for k in range(6):
    b=(target/'stages'/f'stage_{k:02d}.pdb').read_bytes()
    require(sha(b)==r['stage_sha'][str(k)] and compare_pdb(b,traces[case][k])['pass'],'prior baseline changed')
   continue
  require(not target.exists(),'incomplete previous baseline '+case)
  target.mkdir()
  try:
   run=FullReplay(case,backbone[case],None,target).execute()
   check={str(k):compare_pdb((target/'stages'/f'stage_{k:02d}.pdb').read_bytes(),traces[case][k]) for k in range(6)}
   status='PASS' if all(x['pass'] for x in check.values()) else 'FAIL_NATIVE_REPLAY'
   put_json(receipt,{'case':case,'status':status,'stage_sha':run['sha256'],'trace_comparison':check,'input_sha':sha(backbone[case]),'MD_steps_executed':False})
   require(status=='PASS','native original five-step replay failed '+case)
  except Exception:
   (target/'ERROR.txt').write_text(traceback.format_exc());raise

def baseline_gate(args,traces,backbones):
 for case,stages in traces.items():
  dest=args.output/'baseline'/case
  p=dest/'receipt.json'
  require(p.is_file(),'baseline not yet completed '+case)
  meta=json.loads(p.read_text())
  require(meta['status']=='PASS' and meta['input_sha']==sha(backbones[case]),'baseline invalid '+case)
  for k in range(6):
   data=(dest/'stages'/f'stage_{k:02d}.pdb').read_bytes()
   require(sha(data)==meta['stage_sha'][str(k)] and compare_pdb(data,stages[k])['pass'],'baseline drift '+case)

def branch(args,traces,backbones,v2):
 tasks=plan(v2,args.phase);put_json(args.output/f'plan_{args.phase}.json',{'phase':args.phase,'n':len(tasks),'tasks':tasks})
 root=args.output/'branches';root.mkdir(exist_ok=True)
 executed=0;counts=Counter()
 for task in tasks:
  name=task_id(task);dest=root/name;receipt=dest/'receipt.json'
  if receipt.is_file():
   data=json.loads(receipt.read_text())
   require(data['task']==task,'existing receipt identity conflict')
   if data['status']=='PASS':
    for rel,digest in data['files_sha256'].items():
     file=dest/rel;require(file.is_file() and sha(file.read_bytes())==digest,'saved output changed '+name+' '+rel)
   counts[data['status']]+=1;continue
  require(not dest.exists(),'incomplete branch requires manual inspection: '+name)
  if executed>=args.max_branches:counts['PENDING']+=1;continue
  dest.mkdir();print('EXPERIMENT',name,flush=True)
  try:
   replay=FullReplay(task['case'],backbones[task['case']],task,dest)
   res=replay.execute();outcomes=[]
   base=args.output/'baseline'/task['case']
   for k in range(6):
    b=(dest/'stages'/f'stage_{k:02d}.pdb').read_bytes()
    ref=traces[task['case']][k]
    vals={'step':k,'geometry':metrics(b,PAIRS[task['case']]),
          'vs_frozen':compare_pdb(b,ref),
          'all_residue_chi_propagation':network_chi(dest,base,k),
          'sidechain_coordinate_propagation':propagation(b,ref,task['target'])}
    outcomes.append(vals)
   placebo=(task['mode']!='noop' or all(x['vs_frozen']['pass'] for x in outcomes))
   prefix_same=all(x['vs_frozen']['pass'] for x in outcomes[:task['step']])
   same_step_non_target_chi_max=max([float(v) for resid,v in outcomes[task['step']]['all_residue_chi_propagation']['all'].items() if int(resid)!=task['target']] or [0])
   # In the official batch inference, all other predictions at the pulse step
   # should be unchanged. Any greater deviation invalidates clean one-target attribution.
   same_step_isolated=same_step_non_target_chi_max<=.5
   if not prefix_same:status='FAILED_PRE_INTERVENTION_REPLAY'
   elif not same_step_isolated:status='FAILED_SAME_STEP_ISOLATION'
   elif not placebo:status='FAILED_NOOP_NEGATIVE_CONTROL'
   elif not res['event']['v2_match']:status='FAILED_V2_REPRO'
   else:status='PASS'
   created={str(p.relative_to(dest)):sha(p.read_bytes()) for p in dest.rglob('*') if p.is_file()}
   put_json(receipt,{'task':task,'status':status,'event':res['event'],'stage_sha':res['sha256'],
     'geometry_by_stage':outcomes,'files_sha256':created,
     'gates':{'same_pre_intervention_states':prefix_same,'same_step_other_residue_max_chi_delta_deg':same_step_non_target_chi_max,
              'same_step_target_only':same_step_isolated,'placebo_no_effect':placebo,'v2_perturbation_reproduced':res['event']['v2_match']},
     'interpretation':'one-shot feature intervention; all other residue predictions at this step and all subsequent full updates native; not MD',
     'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True})
  except Exception as e:
   (dest/'ERROR.txt').write_text(traceback.format_exc(),encoding='utf-8')
   put_json(receipt,{'task':task,'status':'ERROR','error':str(e),'MD_steps_executed':False})
  executed+=1
  counts[json.loads(receipt.read_text())['status']]+=1
  put_json(args.output/'progress.json',{'phase':args.phase,'completed_this_call':executed,'last':name,'counts_so_far':dict(counts)})
 # Recount in case PENDING increment and progress was partial
 status=Counter()
 for task in tasks:
  p=root/task_id(task)/'receipt.json'
  status[json.loads(p.read_text())['status'] if p.is_file() else 'PENDING']+=1
 result={'phase':args.phase,'total':len(tasks),'status_counts':dict(status),
         'phase_complete_pass':status['PASS']==len(tasks),
         'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
 put_json(args.output/f'status_{args.phase}.json',result)
 print(json.dumps(result,indent=2),flush=True)

def test():
 require(circ(179,-179)==2,'circular angle metric')
 require(sum(map(len,PAIRS.values()))==18 and len(PAIRS)==8,'coverage')
 a=b'ATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00  0.00           C\n'
 require(np.allclose(next(iter(pdb_atoms(a).values()))[0],[1,2,3]),'PDB parse')
 print('SELFTEST PASS; no inference executed')

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--selftest',action='store_true');ap.add_argument('--root',type=Path)
 ap.add_argument('--traces',type=Path,nargs=2);ap.add_argument('--v2',type=Path)
 ap.add_argument('--output',type=Path);ap.add_argument('--phase',choices=['baseline','pilot','all'])
 ap.add_argument('--max-branches',type=int,default=8)
 ap.add_argument('--confirm',default='');args=ap.parse_args()
 if args.selftest:test();return
 require(args.confirm=='I_APPROVE_X4_FULL_FEEDBACK_DIAGNOSTIC_ONLY','diagnostic-only confirmation missing')
 require(args.root and args.traces and args.v2 and args.output and args.phase,'missing arguments')
 require(1<=args.max_branches<=200,'batch cap outside 1..200')
 source=env(args.root);traces,backbones=parse_trace(args.traces);v2=load_v2(args.v2)
 args.output.mkdir(exist_ok=True,parents=True)
 manifest={'schema':'NT-BE-X4-full-feedback-v4','model':source,'trace_zip_sha256':[sha(p.read_bytes()) for p in args.traces],
   'v2_zip_sha256':sha(args.v2.read_bytes()),'inputs_sha256':{k:sha(v) for k,v in backbones.items()},
   'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
 out=args.output/'environment_receipt.json'
 if out.is_file():require(json.loads(out.read_text())==manifest,'changing executable or input environment')
 else:put_json(out,manifest)
 if args.phase=='baseline':baseline(args,traces,backbones)
 else:baseline_gate(args,traces,backbones);branch(args,traces,backbones,v2)



# --- X4 v8: supplementary single-stage mediator-edit experiments ---
# No changes to existing baseline / intervention semantics above.
MEDIATION_CASES = [
 {'branch':'WT_source000_s2_t6_i55_p58_third_third_lag','case':'WT_source000','mediator':17},
 {'branch':'WT_source001_s1_t55_i6_p6_pair_partner_mask','case':'WT_source001','mediator':56},
 {'branch':'mutant_source001_s2_t59_i6_p6_pair_partner_lag','case':'mutant_source001','mediator':4},
]

def verified_receipt(directory):
 p=directory/'receipt.json';require(p.is_file(),'missing branch receipt: '+str(p))
 rec=json.loads(p.read_text(encoding='utf-8'))
 require(rec['status']=='PASS','only PASS branches are eligible')
 for rel,digest in rec.get('files_sha256',{}).items():
  f=directory/rel;require(f.is_file() and sha(f.read_bytes())==digest,'branch file SHA differs: '+str(f))
 return rec

def read_chi_file(folder,step,mediator):
 obj=json.loads((folder/'stages'/f'stage_{step:02d}_chis.json').read_text(encoding='utf-8'))
 require(str(mediator) in obj,'mediator not in chi file')
 return obj[str(mediator)]

def run_edit(case,backbone,task,mediator,step,chi,dest):
 require(not dest.exists(),'refuse overwrite '+str(dest))
 dest.mkdir(parents=True)
 run=FullReplay(case,backbone,task,dest)
 run.mediator_edit={'step':step,'residue':mediator,'chi':chi}
 try:
  result=run.execute()
  require(run.mediator_edit_log is not None,'mediator intervention not applied')
  return {'result':result,'edit':run.mediator_edit_log}
 except Exception:
  (dest/'ERROR.txt').write_text(traceback.format_exc(),encoding='utf-8');raise

def mediation_main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--mediation',action='store_true');ap.add_argument('--root',type=Path,required=True)
 ap.add_argument('--traces',type=Path,nargs=2,required=True)
 ap.add_argument('--v2',type=Path,required=True)
 ap.add_argument('--v4-results',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--confirm',required=True);ap.add_argument('--selftest',action='store_true')
 a=ap.parse_args()
 require(a.confirm=='I_APPROVE_X4_MEDIATION_V8_1_DIAGNOSTIC_ONLY','confirmation mismatch')
 if a.selftest:
  test();require(len(MEDIATION_CASES)==3,'plan integrity');print('V8.1 MEDIATION SELFTEST PASS (no inference)');return
 src=env(a.root);traces,backbones=parse_trace(a.traces);v2=load_v2(a.v2)
 base=a.v4_results
 require(base.is_dir(),'v4 output folder missing')
 gate=base/'environment_receipt.json';require(gate.exists(),'v4 environment receipt missing')
 original=json.loads(gate.read_text(encoding='utf-8'))
 require(original['model']==src,'HPacker executable/hparams/weights diverged from v4')
 require(original['v2_zip_sha256']==sha(a.v2.read_bytes()),'v2 data changed')
 require(original['trace_zip_sha256']==[sha(z.read_bytes()) for z in a.traces],'trace ZIPs changed')
 baseline_gate(type('A',(),{'output':base})(),traces,backbones)
 cases=[]
 for entry in MEDIATION_CASES:
  case=entry['case'];name=entry['branch'];med=entry['mediator']
  receipt=verified_receipt(base/'branches'/name)
  task=receipt['task'];require(task['case']==case and task_id(task)==name,'task identity mismatch')
  step=task['step']+1;require(1<=step<=5,'rescue needs a subsequent update step')
  bc=read_chi_file(base/'baseline'/case,step,med)
  changed=read_chi_file(base/'branches'/name,step,med)
  delta=max(circ(x,y) for x,y in zip(bc,changed));require(delta>1.,'mediator prediction did not change >1 deg')
  cases.append({'case':case,'branch':name,'mediator':med,'step':step,'pulse':task,
                'baseline_chi':bc,'branch_chi':changed,'delta_deg':delta})
 a.output.mkdir(parents=True,exist_ok=True)
 plan=a.output/'v8_plan.json'
 identity={'schema':'X4-v8p1-mediation-pre-reconstruction','cases':cases,'model':src,'trace_sha':[sha(z.read_bytes()) for z in a.traces],
           'v2_sha':sha(a.v2.read_bytes()),'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
 if plan.exists():require(json.loads(plan.read_text(encoding='utf-8'))==identity,'plan mismatch; STOP')
 else:put_json(plan,identity)
 all_results=[]
 for entry in cases:
  case=entry['case']; name=entry['branch'];med=entry['mediator'];step=entry['step']
  outcomes={}
  for kind in ('rescue','sufficiency'):
   dest=a.output/case/kind
   final_receipt=dest/'receipt.json'
   if final_receipt.exists():
    obj=json.loads(final_receipt.read_text(encoding='utf-8'))
    require(obj['status']=='PASS' and obj['case']==case and obj['branch']==name and obj.get('structural_gate',{}).get('forced_chi_recorded_before_rebuild') is True,'previous output not validated V8.1')
    for rel,dig in obj['files_sha256'].items():require(sha((dest/rel).read_bytes())==dig,'saved result changed')
   else:
    task=entry['pulse'] if kind=='rescue' else None
    chis=entry['baseline_chi'] if kind=='rescue' else entry['branch_chi']
    obj=run_edit(case,backbones[case],task,med,step,chis,dest)
    if kind=='rescue':require(obj['result']['event']['v2_match'],'rescue original pulse not reproduced')
    b_control=base/'branches'/name
    before=base/'baseline'/case
    max_pre=0.
    for k in range(step):
     comparator=b_control if kind=='rescue' else before
     diff=compare_pdb((dest/'stages'/f'stage_{k:02d}.pdb').read_bytes(),(comparator/'stages'/f'stage_{k:02d}.pdb').read_bytes())
     max_pre=max(max_pre,diff.get('max_A',float('inf')))
     require(diff['max_A']==0,'pre-mediator state differs from parent')
    # Hard evidence gate: mediator sidechain must be physically displaced at edit step,
    # with every other atom identical to its correct parent at the same step.
    parent=b_control if kind=='rescue' else before
    edited_atoms=pdb_atoms((dest/'stages'/f'stage_{step:02d}.pdb').read_bytes())
    parent_atoms=pdb_atoms((parent/'stages'/f'stage_{step:02d}.pdb').read_bytes())
    require(set(edited_atoms)==set(parent_atoms),'atom identities differ after mediator force')
    mediator_dist=[];other_dist=[]
    for atom in edited_atoms:
     d=float(np.linalg.norm(edited_atoms[atom][0]-parent_atoms[atom][0]))
     if atom[1]==med and atom[-1] not in BB:mediator_dist.append(d)
     else:other_dist.append(d)
    require(mediator_dist and max(mediator_dist)>0.01,
      'MEDIATION_INVALID: chi override not reflected in sidechain coordinates')
    require(not other_dist or max(other_dist)<=0.001,
      'MEDIATION_INVALID: another residue or backbone changed in the edit stage')
    require(max(obj['edit']['delta_from_native_deg'])>1.,
      'MEDIATION_INVALID: edit too small for a substantive intervention')
    structural_gate={'mediator_max_displacement_A':max(mediator_dist),
      'other_atoms_max_displacement_A':max(other_dist,default=0.),
      'forced_chi_recorded_before_rebuild':True}
    geom=[]
    for k in range(6):
     struct=(dest/'stages'/f'stage_{k:02d}.pdb').read_bytes()
     geom.append({'stage':k,'geometry':metrics(struct,PAIRS[case]),
                  'versus_baseline':compare_pdb(struct,(before/'stages'/f'stage_{k:02d}.pdb').read_bytes()),
                  'versus_original_branch':compare_pdb(struct,(b_control/'stages'/f'stage_{k:02d}.pdb').read_bytes())})
    files={str(p.relative_to(dest)):sha(p.read_bytes()) for p in dest.rglob('*') if p.is_file() and p.name!='receipt.json'}
    obj={'status':'PASS','case':case,'branch':name,'mediator':med,'kind':kind,'step':step,
         'mediator_edit':obj['edit'],'geometry':geom,'parent_prefix_max_A':max_pre,'structural_gate':structural_gate,
         'files_sha256':files,'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
    put_json(final_receipt,obj)
   outcomes[kind]={'final':obj['geometry'][-1]['geometry'],'forced':obj['mediator_edit']}
  all_results.append({'case':case,'mediator':med,'original_branch':name,
        'baseline_final':metrics((base/'baseline'/case/'stages'/'stage_05.pdb').read_bytes(),PAIRS[case]),
        'branch_final':metrics((base/'branches'/name/'stages'/'stage_05.pdb').read_bytes(),PAIRS[case]),
        'experiments':outcomes})
  put_json(a.output/'progress.json',{'completed':len(all_results),'total':len(cases)})
 put_json(a.output/'v8_summary.json',{'status':'PASS','cases':all_results,'interpretation':'single-step one-residue do(chi) mediation probes, not general necessity/sufficiency',
                                     'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True})
 print('V8.1 COMPLETE PASS',len(all_results),'cases, 6 controlled experiments')



# V8.2 restricted stage-wise pathway-blocking experimental driver.
# WT0 17 rescue (stage3) -> residue53 (stage4);
# WT1 56 rescue (stage2) -> residue55, 59, 3 (stage3).
BLOCK_PLAN = [
 ('WT_source000', [53]),
 ('WT_source001', [55]),
 ('WT_source001', [59]),
 ('WT_source001', [3]),
 ('WT_source001', [55,59,3]),
]

def all_keys_distance(p1,p2):
 a=pdb_atoms(p1);b=pdb_atoms(p2)
 require(set(a)==set(b),'atom list mismatch')
 return {k:float(np.linalg.norm(a[k][0]-b[k][0])) for k in a}

def block_main():
 ap=argparse.ArgumentParser(description='X4 V8.2 restricted rescue pathway block')
 ap.add_argument('--block',action='store_true')
 ap.add_argument('--selftest',action='store_true')
 ap.add_argument('--root',type=Path,required=True)
 ap.add_argument('--traces',type=Path,nargs=2,required=True)
 ap.add_argument('--v2',type=Path,required=True)
 ap.add_argument('--v4-results',type=Path,required=True)
 ap.add_argument('--v8-results',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--confirm',required=True)
 a=ap.parse_args()
 require(a.confirm=='I_APPROVE_X4_V8_2_PATHWAY_BLOCK_DIAGNOSTIC_ONLY','confirmation missing')
 if a.selftest:
  test(); require(len(BLOCK_PLAN)==5,'plan count')
  require(len(set((c,tuple(r)) for c,r in BLOCK_PLAN))==len(BLOCK_PLAN),'duplicate')
  print('V8.2 SELFTEST PASS; no inference');return
 src=env(a.root)
 traces,backbones=parse_trace(a.traces);v2=load_v2(a.v2)
 v4=a.v4_results;v8=a.v8_results
 require(v4.is_dir() and v8.is_dir(),'source folder unavailable')
 e=json.loads((v4/'environment_receipt.json').read_text(encoding='utf-8'))
 require(e['model']==src,'model diverged from v4')
 require(e['v2_zip_sha256']==sha(a.v2.read_bytes()),'v2 changed')
 require(e['trace_zip_sha256']==[sha(q.read_bytes()) for q in a.traces],'trace changed')
 baseline_gate(type('Tmp',(),{'output':v4})(),traces,backbones)
 vp=json.loads((v8/'v8_plan.json').read_text(encoding='utf-8'))
 require(vp['schema']=='X4-v8p1-mediation-pre-reconstruction','v8.1 provenance mismatch')
 metadata={x['case']:x for x in vp['cases']}
 require(all(c in metadata for c,_ in BLOCK_PLAN),'missing case in v8 plan')
 plan={'schema':'NT-BE-X4-v8.2-pathway-block','cases':BLOCK_PLAN,'v8_plan_sha':sha((v8/'v8_plan.json').read_bytes()),'source':src,
       'v2_sha':sha(a.v2.read_bytes()),'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
 a.output.mkdir(parents=True,exist_ok=True)
 pf=a.output/'v8_2_plan.json'
 if pf.exists():require(json.loads(pf.read_text(encoding='utf-8'))==plan,'existing plan mismatch')
 else:put_json(pf,plan)
 summary=[]
 for case,residues in BLOCK_PLAN:
  info=metadata[case];name=info['branch'];pulse=info['pulse'];step=info['step'];block_step=step+1
  require(block_step<=5,'block step outside frozen refinement')
  med=info['mediator']
  original=v4/'branches'/name
  verified_receipt(original)
  rescue=v8/case/'rescue'
  rescue_receipt=verified_receipt(rescue)
  require(rescue_receipt.get('structural_gate',{}).get('forced_chi_recorded_before_rebuild') is True,'invalid rescue parent')
  require(rescue_receipt['step']==step and rescue_receipt['mediator']==med,'rescue identity mismatch')
  forced=[{'step':step,'residue':med,'chi':info['baseline_chi']}]
  before_chi={}
  for r in residues:
   target=read_chi_file(original,block_step,r)
   rescue_chi=read_chi_file(rescue,block_step,r)
   delta=max(circ(x,y) for x,y in zip(target,rescue_chi))
   require(delta>1.,f'candidate lacks >1 degree contrast: {case} {r} {delta}')
   before_chi[str(r)]={'rescue_chi':rescue_chi,'original_branch_chi':target,'pre_edit_delta_deg':delta}
   forced.append({'step':block_step,'residue':r,'chi':target})
  taskname=case+'_block_s'+str(block_step)+'_'+('_'.join(map(str,residues)))
  dest=a.output/taskname;rc=dest/'receipt.json'
  if rc.exists():
   obj=json.loads(rc.read_text(encoding='utf-8'))
   require(obj['status']=='PASS' and obj['case']==case and obj['block_residues']==residues,'existing receipt mismatch')
   for rel,digest in obj['files_sha256'].items():require(sha((dest/rel).read_bytes())==digest,'output corruption')
  else:
   require(not dest.exists(),'partial output exists, preserve for forensic review: '+str(dest))
   dest.mkdir()
   try:
    replay=FullReplay(case,backbones[case],pulse,dest);replay.mediator_edit=forced
    run=replay.execute()
    require(run['event']['v2_match'],'pulse no longer matches v2')
    require(len(replay.mediator_edit_logs)==len(forced),'not every forced chi applied')
    refprefix=[]
    for st in range(block_step):
     d=all_keys_distance((dest/'stages'/f'stage_{st:02d}.pdb').read_bytes(),(rescue/'stages'/f'stage_{st:02d}.pdb').read_bytes())
     refprefix.append(max(d.values()))
    require(max(refprefix)==0,'pre-block PDB differs from validated rescue')
    diffs=all_keys_distance((dest/'stages'/f'stage_{block_step:02d}.pdb').read_bytes(),(rescue/'stages'/f'stage_{block_step:02d}.pdb').read_bytes())
    affected={k:v for k,v in diffs.items() if k[1] in residues and k[-1] not in BB}
    unaffected={k:v for k,v in diffs.items() if not (k[1] in residues and k[-1] not in BB)}
    per_res={str(r):max((v for k,v in affected.items() if k[1]==r),default=0) for r in residues}
    require(all(v>.01 for v in per_res.values()),'block chi not reflected in coordinates')
    require(max(unaffected.values(),default=0)<=.001,'off-target structural edit at pulse')
    stages=[]
    for st in range(6):
     b=(dest/'stages'/f'stage_{st:02d}.pdb').read_bytes()
     original_pdb=(original/'stages'/f'stage_{st:02d}.pdb').read_bytes()
     rescue_pdb=(rescue/'stages'/f'stage_{st:02d}.pdb').read_bytes()
     stages.append({'stage':st,'geometry':metrics(b,PAIRS[case]),
                    'vs_rescue':compare_pdb(b,rescue_pdb,tol=0),
                    'vs_original_branch':compare_pdb(b,original_pdb,tol=0),
                    'propagation_vs_rescue':propagation(b,rescue_pdb,med)})
    evidence={'case':case,'branch':name,'rescue_mediator':med,'rescue_step':step,
       'block_step':block_step,'block_residues':residues,'contrasts':before_chi,
       'edits':replay.mediator_edit_logs,'preblock_max_A':max(refprefix),
       'same_step_changed_residue_max_A':per_res,
       'same_step_unchanged_atom_max_A':max(unaffected.values(),default=0),
       'by_stage':stages}
    put_json(dest/'observations.json',evidence)
    files={str(p.relative_to(dest)):sha(p.read_bytes()) for p in dest.rglob('*') if p.is_file() and p.name!='receipt.json'}
    obj={'status':'PASS','case':case,'block_residues':residues,'files_sha256':files,
         'scientific_scope':'path-specific single-stage rescue-background block, not unique necessity or sufficiency',
         'MD_steps_executed':False,'Core_STOP':True,'Publication_STOP':True}
    put_json(rc,obj)
   except Exception as exc:
    (dest/'ERROR.txt').write_text(traceback.format_exc(),encoding='utf-8')
    put_json(rc,{'status':'ERROR','case':case,'block_residues':residues,'error':str(exc)})
    raise
  obs=json.loads((dest/'observations.json').read_text(encoding='utf-8'))
  summary.append({'task':taskname,'case':case,'residues':residues,
     'final':obs['by_stage'][-1]['geometry'],
     'preblock_max_A':obs['preblock_max_A'],
     'unchanged_atom_max_A':obs['same_step_unchanged_atom_max_A']})
  put_json(a.output/'progress.json',{'completed':len(summary),'total':len(BLOCK_PLAN),'results':summary})
 put_json(a.output/'final_summary.json',{'status':'PASS','results':summary,'not_MD':True,'Core_STOP':True,'Publication_STOP':True})
 print('V8.2 COMPLETED: 5 pathway-block experiments; no MD',flush=True)

if __name__=='__main__':
 try:
  if '--block' in sys.argv: block_main()
  elif '--mediation' in sys.argv: mediation_main()
  else:main()
 except Exception:traceback.print_exc();sys.exit(1)
