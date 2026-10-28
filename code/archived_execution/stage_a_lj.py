"""Stage A: reproducible, heavy-heavy LJ decomposition, no simulation or solvent.
OpenMM constructs protein-only amber99sb.xml protein system AFTER adding hydrogen atoms.
Use Force exceptions for 1-4 scaling; pairwise terms are exact with NoCutoff
in this protein-only surrogate, NOT PME solvated forcefield full energy.
"""
import argparse,csv,hashlib,importlib.metadata,json,math,pathlib,subprocess,sys,datetime,io
import numpy as np
from collections import defaultdict
HERE=pathlib.Path(__file__).resolve().parent
CASES=['wt0_baseline','wt0_v8_1_rescue','wt1_baseline','wt1_v8_1_rescue']
EXPECTED='42c551133f9825f4ec01812075b0d50e268c2209'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')
def check(bio):
 manifest=json.loads((HERE/'input_manifest.json').read_text(encoding='utf8'))
 for c in CASES:
  if sha(HERE/(c+'.pdb'))!=manifest['files'][c]:raise RuntimeError('INPUT SHA MISMATCH '+c)
 def git(*cmd):return subprocess.run(['git','-C',str(bio),*cmd],text=True,capture_output=True,check=True).stdout.strip()
 if git('rev-parse','HEAD')!=EXPECTED or git('status','--porcelain'):raise RuntimeError('FROZEN BIOEMU GIT MISMATCH/DIRTY')
 import openmm,openmm.app as app,openmm.unit as u
 if importlib.metadata.version('OpenMM')!='8.2.0':raise RuntimeError('Expected OpenMM 8.2.0')
 return openmm,app,u

def typed_pdb_with_transient_oxt(pdb_path,app):
 # Generate *in memory* an OXT for C-terminal residue 63. Source PDB stays untouched.
 lines=pdb_path.read_text(encoding='utf8').splitlines(keepends=True)
 terminal=[line for line in lines if line.startswith(('ATOM  ','HETATM')) and line[22:26].strip()=='63']
 if {line[12:16].strip() for line in terminal} & {'OXT'}:
  raise RuntimeError('Unexpected source OXT; source should have 494 heavy atoms')
 xyz={line[12:16].strip():np.array([float(line[i:i+8]) for i in (30,38,46)]) for line in terminal}
 if not {'CA','C','O'}.issubset(xyz):raise RuntimeError('Missing C-terminal backbone atoms')
 C,CA,O=(xyz[k] for k in ('C','CA','O'))
 uo=(O-C)/np.linalg.norm(O-C);uca=(CA-C)/np.linalg.norm(CA-C)
 direction=-(uo+uca);direction/=np.linalg.norm(direction)
 oxt=C+1.25*direction
 template=next(line for line in terminal if line[12:16].strip()=='O')
 serial=max(int(line[6:11]) for line in lines if line.startswith(('ATOM  ','HETATM')))+1
 # Coordinates of temporary OXT serve only for protonation/template assignment.
 oxtline=(f'ATOM  {serial:5d}  OXT {template[17:26]}   '
          f'{oxt[0]:8.3f}{oxt[1]:8.3f}{oxt[2]:8.3f}'
          f'  1.00  0.00           O  \n')
 insertion=next(i for i,line in enumerate(lines) if line.startswith(('TER','END')))
 lines.insert(insertion,oxtline)
 return app.PDBFile(io.StringIO(''.join(lines)))

def run_one(c,app,u,openmm):
 pdb=app.PDBFile(str(HERE/(c+'.pdb')))
 typing_pdb=typed_pdb_with_transient_oxt(HERE/(c+'.pdb'),app)
 ff=app.ForceField('amber99sb.xml','tip3p.xml')
 model=app.Modeller(typing_pdb.topology,typing_pdb.positions)
 # No solvent. Hydrogens used only for consistent OpenMM typing and exception graph;
 # all measured pairs are heavy atoms with original PDB coordinates.
 model.addHydrogens(ff,pH=7.0)
 system=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds,rigidWater=True)
 nb=[f for f in system.getForces() if isinstance(f,openmm.NonbondedForce)]
 if len(nb)!=1:raise RuntimeError('Expected one NonbondedForce')
 nb=nb[0]
 atoms=list(model.topology.atoms()); heavy=[a for a in atoms if a.element is not None and a.element.symbol!='H' and a.name!='OXT']
 # Direct PDB heavy atom name+residue ordering retained; verify against originals
 orig=[a for a in pdb.topology.atoms() if a.element is not None and a.element.symbol!='H']
 def key(a):return (a.residue.index,a.residue.name,a.name)
 if [key(a) for a in heavy]!=[key(a) for a in orig]:raise RuntimeError('PDB heavy atom mapping mismatch')
 if len(heavy)!=494:raise RuntimeError('Expected 494 heavy atoms')
 xyz=np.asarray(model.positions.value_in_unit(u.nanometer));idx=[a.index for a in heavy]
 if not np.allclose(xyz[idx],np.asarray(pdb.positions.value_in_unit(u.nanometer)),atol=1e-8):raise RuntimeError('heavy coordinates changed by protonation')
 pars=[]
 for atom in atoms:
  q,s,e=nb.getParticleParameters(atom.index)
  pars.append((float(s.value_in_unit(u.nanometer)),float(e.value_in_unit(u.kilojoule_per_mole))))
 exceptions={}
 for k in range(nb.getNumExceptions()):
  a,b,q,s,e=nb.getExceptionParameters(k); i,j=sorted([int(a),int(b)])
  exceptions[(i,j)]=(float(s.value_in_unit(u.nanometer)),float(e.value_in_unit(u.kilojoule_per_mole)))
 rows=[];tot=defaultdict(float);near=[]
 for i,ai in enumerate(heavy):
  for aj in heavy[i+1:]:
   if ai.residue.index==aj.residue.index:continue  # focus nonadjacent and other interresidue heavy atom pairs
   if abs(ai.residue.index-aj.residue.index)<=1:continue
   ii,jj=sorted((ai.index,aj.index))
   if (ii,jj) in exceptions:
    sig,eps=exceptions[(ii,jj)];kind='exception'
   else:
    sig=(pars[ii][0]+pars[jj][0])/2;eps=math.sqrt(pars[ii][1]*pars[jj][1]);kind='standard'
   d=float(np.linalg.norm(xyz[ii]-xyz[jj]))
   if d<=0:raise RuntimeError('coincident heavy atoms')
   sr=sig/d if eps>0 else 0.0
   rep=4*eps*sr**12;attr=-4*eps*sr**6
   lj=rep+attr
   if not all(math.isfinite(x) for x in [rep,attr,lj]):raise RuntimeError('nonfinite LJ')
   tot['rep']+=rep;tot['attr']+=attr;tot['net']+=lj;tot['pair_count']+=1
   if d<0.16:near.append((ai,aj,d,lj))
   rows.append({'pair_id':f'{ai.residue.index+1}:{ai.residue.name}:{ai.name}|{aj.residue.index+1}:{aj.residue.name}:{aj.name}',
    'res_i':ai.residue.index+1,'res_j':aj.residue.index+1,'atom_i':ai.name,'atom_j':aj.name,
    'distance_A':d*10,'epsilon_kj_mol':eps,'sigma_nm':sig,'kind':kind,'lj_repulsive_kj_mol':rep,'lj_attractive_kj_mol':attr,'lj_net_kj_mol':lj})
 return rows,{'case':c,'heavy_atoms':len(heavy),'all_particles':len(atoms),'exceptions_total':len(exceptions),'nonlocal_interresidue_heavy_pairs':int(tot['pair_count']),
 'lj_repulsive_kj_mol':tot['rep'],'lj_attractive_kj_mol':tot['attr'],'lj_net_kj_mol':tot['net'],
 'under_1p6A_count':len(near),'min_interresidue_dist_A':min(r['distance_A'] for r in rows),
 'short_contacts':[{'pair':f'{ai.residue.index+1}:{ai.name}|{aj.residue.index+1}:{aj.name}','distance_A':d*10,'lj_net_kj_mol':lj} for ai,aj,d,lj in near],
 'interpretation':'NoCutoff protein-only, heavy-heavy, nonadjacent-residue pairwise LJ; excludes hydrogen-containing and adjacent/same-residue terms. NOT original solvated PE.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--bioemu',default=str(HERE.parent/'bioemu'));p.add_argument('--execute',action='store_true');args=p.parse_args()
 openmm,app,u=check(pathlib.Path(args.bioemu).resolve())
 if not args.execute:
  print('DRY RUN PASS: SHA, frozen git and OpenMM; no system built / simulation.');return
 out=HERE/'stage_a_results';out.mkdir(exist_ok=True)
 meta={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bioemu_commit':EXPECTED,
       'forcefield':'OpenMM amber99sb.xml, tip3p.xml (protein-only NoCutoff surrogate)',
       'selection':'heavy-heavy, interresidue, residue-index separation >=2',
       'exclusions':'original NonbondedForce exceptions (1-2,1-3,1-4) honored; epsilon=0 skipped',
       'limitations':'No solvent, ions, hydrogens LJ, PME, bonded energy, or whole-system initial PE. Temporary C-terminal OXT is constructed solely for FF atom typing, is excluded from 494-heavy-atom analysis; source inputs never edited. Modeller.addHydrogens may optimize added hydrogens.',
       'input_sha256':json.loads((HERE/'input_manifest.json').read_text())['files']}
 results={};allrows={}
 for c in CASES:
  rows,s=run_one(c,app,u,openmm);results[c]=s;allrows[c]={r['pair_id']:r for r in rows}
  with (out/(c+'_pairs.csv')).open('w',newline='',encoding='utf8') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 diffs={}
 for base,rescue in [('wt0_baseline','wt0_v8_1_rescue'),('wt1_baseline','wt1_v8_1_rescue')]:
  a,b=allrows[base],allrows[rescue]
  if a.keys()!=b.keys():raise RuntimeError('PAIR KEY MISMATCH')
  if any((a[k]['epsilon_kj_mol']!=b[k]['epsilon_kj_mol'] or a[k]['sigma_nm']!=b[k]['sigma_nm'] or a[k]['kind']!=b[k]['kind']) for k in a):raise RuntimeError('LJ PARAMETER MISMATCH between paired structures')
  dr=[]
  for key in a:
   x,y=a[key],b[key];dr.append({'pair_id':key,'distance_baseline_A':x['distance_A'],'distance_rescue_A':y['distance_A'],
    'net_baseline_kj_mol':x['lj_net_kj_mol'],'net_rescue_kj_mol':y['lj_net_kj_mol'],
    'delta_net_rescue_minus_base_kj_mol':y['lj_net_kj_mol']-x['lj_net_kj_mol'],
    'delta_repulsive_rescue_minus_base_kj_mol':y['lj_repulsive_kj_mol']-x['lj_repulsive_kj_mol']})
  dr.sort(key=lambda t:abs(t['delta_net_rescue_minus_base_kj_mol']),reverse=True)
  with (out/(base[:3]+'_delta_pairs.csv')).open('w',newline='',encoding='utf8') as f:
   w=csv.DictWriter(f,fieldnames=list(dr[0]));w.writeheader();w.writerows(dr)
  diffs[base[:3]]={'baseline':base,'rescue':rescue,
   'delta_net_rescue_minus_base_kj_mol':results[rescue]['lj_net_kj_mol']-results[base]['lj_net_kj_mol'],
   'delta_repulsive_rescue_minus_base_kj_mol':results[rescue]['lj_repulsive_kj_mol']-results[base]['lj_repulsive_kj_mol'],
   'top_20_absolute_pair_contributions':dr[:20]}
 dump(out/'stage_a_summary.json',{'method':meta,'cases':results,'paired_differences':diffs})
 files=sorted(p for p in out.iterdir() if p.is_file());dump(out/'results_sha256.json',{p.name:sha(p) for p in files})
 print(json.dumps({'status':'STAGE_A_FINISHED','summary':str(out/'stage_a_summary.json'),'comparisons':{k:{j:v for j,v in x.items() if j!='top_20_absolute_pair_contributions'} for k,x in diffs.items()}},indent=2))
if __name__=='__main__':main()
