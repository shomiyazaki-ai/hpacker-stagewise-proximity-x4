"""X4 Phase C: one case per invocation, fixed BioEmu relax MD and persistent diagnostics."""
import argparse,datetime,hashlib,importlib.metadata,json,math,os,pathlib,subprocess,sys,time,traceback
ROOT=pathlib.Path(__file__).resolve().parent
ORDER=['wt0_baseline','wt0_v8_1_rescue','wt1_baseline','wt1_v8_1_rescue']
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def ts():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
def event(p,k,**data):
 with p.open('a',encoding='utf-8') as f:
  f.write(json.dumps({'utc':ts(),'event':k,**data},allow_nan=False,ensure_ascii=False,default=str)+'\n');f.flush();os.fsync(f.fileno())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bioemu',default=str(ROOT.parent/'bioemu'));ap.add_argument('--case',choices=ORDER,required=True);ap.add_argument('--execute',action='store_true');ap.add_argument('--worker',action='store_true');ap.add_argument('--outdir');ap.add_argument('--timeout-sec',type=int,default=2400)
 a=ap.parse_args()
 if not 60<=a.timeout_sec<=2400:raise ValueError('timeout 60..2400s; no automatic extension')
 bio=pathlib.Path(a.bioemu).resolve();m=json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'));inp=ROOT/m['inputs'][a.case]['file']
 if sha(inp)!=m['inputs'][a.case]['sha256']:raise RuntimeError('PDB SHA mismatch')
 git=lambda *args:subprocess.run(['git','-C',str(bio),*args],capture_output=True,text=True,check=True).stdout.strip()
 if git('rev-parse','HEAD')!=m['bioemu_commit'] or git('status','--porcelain'):raise RuntimeError('source commit/clean gate failed')
 import mdtraj,openmm,numpy
 if importlib.metadata.version('OpenMM')!='8.2.0':raise RuntimeError('OpenMM version mismatch')
 pdb=mdtraj.load_pdb(str(inp))
 if (pdb.n_frames,pdb.n_atoms,pdb.n_residues)!=(1,494,63):raise RuntimeError('input topology mismatch')
 meta={'case':a.case,'input_sha256':sha(inp),'bioemu_commit':m['bioemu_commit'],'python':sys.version,'openmm':importlib.metadata.version('OpenMM'),'mdtraj':mdtraj.__version__,'numpy':numpy.__version__,'time_limit_seconds':a.timeout_sec,'protocol':'run_one_md(traj,only_energy_minimization=False,simtime_ns=0), fixed source default','diagnostic':'Simulation.step split in <=1000 steps with read-only getState probes; these introduce synchronization overhead and can alter stochastic trajectories, no force-field/thermostat/integrator parameters changed'}
 if not a.execute and not a.worker:
  print(json.dumps(meta,indent=2,ensure_ascii=False));print('DRY RUN ONLY — NO MD EXECUTED');return
 if a.worker:
  worker(pdb,bio,pathlib.Path(a.outdir));return
 out=ROOT/'outputs_phase_c'/(a.case+'_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));out.mkdir(parents=True,exist_ok=False);dump(out/'environment.json',meta)
 cmd=[sys.executable,'-u',str(ROOT/'run_phase_c.py'),'--bioemu',str(bio),'--case',a.case,'--worker','--outdir',str(out),'--timeout-sec',str(a.timeout_sec)]
 start=time.monotonic();status='UNKNOWN';terminate=None
 with (out/'stdout.txt').open('wb') as f1,(out/'stderr.txt').open('wb') as f2:
  proc=subprocess.Popen(cmd,stdout=f1,stderr=f2,stdin=subprocess.DEVNULL)
  try:status='PASS' if proc.wait(timeout=a.timeout_sec)==0 else 'FAIL'
  except subprocess.TimeoutExpired:
   status='TIMEOUT'
   try:
    if os.name=='nt':
     rr=subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],capture_output=True,timeout=20);terminate={'returncode':rr.returncode,'stdout':rr.stdout.decode(errors='replace'),'stderr':rr.stderr.decode(errors='replace')}
    else:proc.kill()
    proc.wait(timeout=20)
   except Exception as exc:status='TIMEOUT_TERMINATION_UNCONFIRMED';terminate={'error':repr(exc)}
 summary={'case':a.case,'status':status,'elapsed_seconds':round(time.monotonic()-start,3),'exit_code':proc.poll(),'termination':terminate,'output_dir':str(out),'utc':ts()}
 dump(out/'summary.json',summary)
 # Hash saved files, including logs; hash file created last to prevent self-reference.
 hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='sha256_manifest.json'}
 dump(out/'sha256_manifest.json',hashes)
 print(json.dumps(summary,indent=2,ensure_ascii=False));print('RESULT FOLDER:',out)
def worker(traj,bio,out):
 import numpy as np,openmm.app as app,openmm.unit as u
 sys.path.insert(0,str(bio/'src'))
 import bioemu.sidechain_relax as side
 log=out/'events.jsonl';original_prep=side._prepare_system;original_min=app.Simulation.minimizeEnergy;original_step=app.Simulation.step
 calls=[0];total=[0];phase_counts={'NVT':0,'NPT':0};last=[None]
 def ev(k,**d):event(log,k,**d)
 def checkpoint(sim,label):
  try:
   st=sim.context.getState(getEnergy=True,getPositions=True,getVelocities=True,getForces=True)
   pe=float(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole));ke=float(st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole))
   pos=np.asarray(st.getPositions(asNumpy=True).value_in_unit(u.nanometer));vel=np.asarray(st.getVelocities(asNumpy=True).value_in_unit(u.nanometer/u.picosecond));force=np.asarray(st.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole/u.nanometer))
   finite={name:int((~np.isfinite(arr)).sum()) for name,arr in [('position',pos),('velocity',vel),('force',force)]}
   mx=float(np.max(np.linalg.norm(force,axis=1))) if not finite['force'] else 'NONFINITE'
   ev('ENERGY_STATE',label=label,step_total=total[0],potential_kj_mol=pe if math.isfinite(pe) else str(pe),kinetic_kj_mol=ke if math.isfinite(ke) else str(ke),nonfinite=finite,max_force_kj_mol_nm=mx,n_particles=len(pos))
   last[0]={'label':label,'step_total':total[0],'potential_kj_mol':pe if math.isfinite(pe) else str(pe),'kinetic_kj_mol':ke if math.isfinite(ke) else str(ke),'nonfinite':finite}
  except BaseException as ex:ev('ENERGY_STATE_ERROR',label=label,error=repr(ex))
 def prep(*args,**kw):
  ev('PREPARE_BEGIN');ret=original_prep(*args,**kw)
  sysmod,mod=ret
  box=mod.topology.getPeriodicBoxVectors();vec=[[float(v.value_in_unit(u.nanometer)) for v in row] for row in box] if box else None
  ev('PREPARE_END',n_particles=mod.topology.getNumAtoms(),water_residues=sum(1 for r in mod.topology.residues() if r.name in ('HOH','WAT')),box_vectors_nm=vec);return ret
 def minim(self,*args,**kw):
  platform=self.context.getPlatform();props={}
  for key in platform.getPropertyNames():
   try:props[key]=platform.getPropertyValue(self.context,key)
   except Exception as ex:props[key]='UNAVAILABLE:'+repr(ex)
  ev('MIN_BEGIN',platform=platform.getName(),platform_properties=props);checkpoint(self,'BEFORE_MIN');ans=original_min(self,*args,**kw);ev('MIN_END');checkpoint(self,'AFTER_MIN');return ans
 def step(self,n):
  calls[0]+=1;cid=calls[0];stage=('INITIAL_'+str(cid) if cid<=4 else 'NVT' if cid<=104 else 'NPT')
  # Official method makes first 4 calls for stabilization, 100 NVT calls x1000, 100 NPT calls x4000.
  ev('STEP_CALL_BEGIN',call=cid,stage=stage,requested_steps=int(n),total_before=total[0]);remain=int(n)
  if cid==1:checkpoint(self,'BEFORE_FIRST_MD_STEP')
  while remain:
   chunk=min(1000,remain)
   try:original_step(self,chunk)
   except BaseException as ex:
    ev('STEP_ERROR',call=cid,stage=stage,total_completed=total[0],error=repr(ex));checkpoint(self,'AT_EXCEPTION');raise
   total[0]+=chunk;remain-=chunk
   # only 1st phase every 1000 steps; later every official call completion.
   if (cid==1 and total[0]%5000==0) or remain==0:checkpoint(self,stage+'_after_'+str(total[0]))
  ev('STEP_CALL_END',call=cid,stage=stage,total_completed=total[0])
 side._prepare_system=prep;app.Simulation.minimizeEnergy=minim;app.Simulation.step=step
 ev('START',sha256_input=sha(ROOT/(out.name.split('_2026')[0]+'.pdb')) if False else 'verified_by_parent')
 try:
  result=side.run_one_md(traj,only_energy_minimization=False,simtime_ns=0)
  ev('MD_FUNCTION_RETURN',n_atoms=result.n_atoms,steps_completed=total[0])
  if total[0]!=611100:raise RuntimeError(f'Unexpected steps {total[0]} != 611100; not COMPLETE')
  # source returns protein heavy atom only, store exact returned result
  p=out/'final_heavy_atoms.pdb';result.save_pdb(str(p));ev('FINAL_PDB_WRITTEN',sha256=sha(p),size_bytes=p.stat().st_size)
  ev('COMPLETE',final_potential_kj_mol=last[0].get('potential_kj_mol') if last[0] else None,latest_state=last[0],steps_completed=total[0])
 except BaseException as ex:
  ev('ERROR',error=repr(ex),step_total=total[0],last_state=last[0]);traceback.print_exc();raise
 finally:ev('EXIT',step_total=total[0])
if __name__=='__main__':main()
