#!/usr/bin/env python3
"""NT-BE-X4 v8.3. Four pre-registered single-step direct-input contrasts.
Diagnostic-only; no MD, minimization, new acceptance or publication.
"""
import argparse, hashlib, io, json, sys, tempfile, traceback, zipfile
from pathlib import Path
import numpy as np

PLAN = [
    ('WT_source000', 'WT_source000_block_s4_53', 4, 53, 55),
    ('WT_source001', 'WT_source001_block_s3_55', 3, 55, 56),
    ('WT_source001', 'WT_source001_block_s3_55', 3, 55, 3),
    ('WT_source001', 'WT_source001_block_s3_55', 3, 55, 59),
]
EXPECTED_COMMIT = '6cd5ec33961bb2662b51c0a15d11d6fd9162f2a7'
BB = {'N','CA','C','O','OXT'}

def require(x, message):
    if not x: raise RuntimeError('STOP: '+message)
def sha(x): return hashlib.sha256(x).hexdigest()
def write_json(path, value):
    path=Path(path); tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8');tmp.replace(path)
def cdiff(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    require(a.shape==b.shape,'chi shape mismatch')
    return float(np.max(np.abs((a-b+180)%360-180)))
def atom_dict(pdb_bytes):
    out={}
    for line in pdb_bytes.decode('ascii').splitlines():
        if not line.startswith(('ATOM  ','HETATM')):continue
        key=(line[21:22],int(line[22:26]),line[26:27],line[12:16].strip())
        require(key not in out,'duplicate PDB atom')
        out[key]=np.asarray([float(line[i:i+8]) for i in (30,38,46)],float)
    return out
def zip_entry(zip_path, suffix):
    with zipfile.ZipFile(zip_path) as z:
        require(z.testzip() is None, 'bad ZIP CRC: '+str(zip_path))
        hits=[name for name in z.namelist() if name.replace('\\','/').endswith(suffix)]
        require(len(hits)==1,'missing/nonunique ZIP file suffix '+suffix+': '+str(hits))
        return z.read(hits[0])
def snapshot(z, case, task, stage, source):
    if source=='rescue':
        end=f'/{case}/rescue/stages/stage_{stage:02d}.pdb'
    else:
        end=f'/{task}/stages/stage_{stage:02d}.pdb'
    return zip_entry(z,end)
def next_chi(z,case,task,stage,source,receiver):
    if source=='rescue':end=f'/{case}/rescue/stages/stage_{stage:02d}_chis.json'
    else:end=f'/{task}/stages/stage_{stage:02d}_chis.json'
    return json.loads(zip_entry(z,end))[str(receiver)]
def obs_nb(hp,r):
    nb,bb=hp._get_neighborhood(r,1)
    feat=np.asarray(hp._get_zernikegram(nb,1).detach().cpu().numpy(),dtype=np.float32)
    def cv(v):return v.decode().strip() if isinstance(v,bytes) else str(v).strip()
    atoms=[{'residue':int(cv(res[4])),'atom':cv(an),'element':cv(el),
            'xyz_CA_A':[float(v) for v in xyz], 'charge':float(charge)}
           for res,an,el,xyz,charge in zip(nb['res_ids'],nb['atom_names'],nb['elements'],nb['coords'],nb['charges'])]
    return feat, atoms

def predict(hp,receiver,donor,donor_atoms=None):
    import torch
    ids=list(hp.get_res_ids())
    valid=[]
    from hpacker.src.sidechain_reconstruction.biopython_internal_coords.constants import CHI_ANGLES
    for r in ids:
        name=hp._get_resname_from_res_id(r)
        if name in ('ALA','GLY'): continue
        require(name in CHI_ANGLES,'unknown residue '+name)
        if len(CHI_ANGLES[name]):valid.append(r)
    recv_ids=[r for r in valid if int(r[1])==receiver]
    require(len(recv_ids)==1,'receiver not unique valid-chis')
    recv=recv_ids[0]
    donor_res=[r for r in hp.structure.get_residues() if int(r.id[1])==donor]
    require(len(donor_res)==1,'donor residue mismatch')
    moving={a.name:a for a in donor_res[0] if a.name not in BB}
    require(moving,'no donor sidechain')
    orig={name:a.coord.copy() for name,a in moving.items()}
    if donor_atoms is not None:
        require(set(moving)==set(donor_atoms),'donor sidechain atom set changed')
        for name,atom in moving.items():atom.coord=np.asarray(donor_atoms[name],dtype=atom.coord.dtype)
    raw=[]
    def hook(model,inputs,output):
        y=output[0] if isinstance(output,(tuple,list)) else output
        raw.append(np.asarray(y.detach().cpu().numpy(),float))
    handle=hp.model[1].register_forward_hook(hook)
    try:
        feature,neighbors=obs_nb(hp,recv)
        with torch.no_grad(): angles=hp._get_predicted_chi_angles_batch(valid,1)
        require(len(raw)==1 and raw[0].shape==(len(valid),8),'unexpected raw output shape')
        i=valid.index(recv)
        return {'feature':feature,'neighbors':neighbors,'raw':raw[0][i],
                'chi':np.asarray(angles[i],float), 'n_batch':len(valid), 'index':i,
                'all_chi':{str(int(r[1])):[float(v) for v in angles[k]] for k,r in enumerate(valid)}}
    finally:
        handle.remove()
        for name,atom in moving.items():atom.coord=orig[name]

def run(a):
    require(not a.output.exists(),'output already exists; no overwrites')
    require(a.root.exists() and a.v8.exists() and a.v82.exists(),'source missing')
    import subprocess
    actual=subprocess.check_output(['git','-C',str(a.root/'hpacker'),'rev-parse','HEAD'],text=True).strip()
    require(actual==EXPECTED_COMMIT,'HPacker commit mismatch')
    sys.path.insert(0,str((a.root/'hpacker').resolve()))
    from hpacker import HPacker
    a.output.mkdir(parents=True)
    results=[]
    try:
        for case,task,stage,donor,receiver in PLAN:
            rid=f'{case}_s{stage}_{donor}_to_{receiver}'
            rescue=snapshot(a.v8,case,task,stage,'rescue')
            blocked=snapshot(a.v82,case,task,stage,'block')
            atom_r=atom_dict(rescue);atom_b=atom_dict(blocked)
            require(set(atom_r)==set(atom_b),'source PDB atom IDs differ')
            changed=[key for key in atom_r if np.linalg.norm(atom_r[key]-atom_b[key])>0.0001]
            require(changed and all(key[1]==donor and key[-1] not in BB for key in changed),
                    'source PDB differs outside donor sidechain')
            donors={k[3]:v for k,v in atom_r.items() if k[1]==donor and k[3] not in BB}
            with tempfile.TemporaryDirectory() as td:
                p=Path(td)/'source.pdb';p.write_bytes(blocked)
                hp=HPacker(str(p));hp.model[1].eval()
                B=predict(hp,receiver,donor)
                R=predict(hp,receiver,donor,donors)
                # Entire restored source PDB must match rescue within PDB coordinate precision.
                p.write_bytes(rescue)
                hr=HPacker(str(p));hr.model[1].eval()
                P=predict(hr,receiver,donor)
            chi_next_block=next_chi(a.v82,case,task,stage+1,'block',receiver)
            chi_next_rescue=next_chi(a.v8,case,task,stage+1,'rescue',receiver)
            record={
                'case':case,'source_stage':stage,'predicted_next_stage':stage+1,
                'donor':donor,'receiver':receiver,'blocked_task':task,
                'source_blocked_pdb_sha256':sha(blocked),'source_rescue_pdb_sha256':sha(rescue),
                'donor_changed_atom_count':len(changed),
                'n_batch':B['n_batch'], 'target_batch_index':B['index'],
                'blocked_vs_restored':{
                    'zernike_max_abs':float(np.max(np.abs(B['feature']-R['feature']))),
                    'zernike_l2':float(np.linalg.norm((B['feature']-R['feature']).astype(float))),
                    'raw_yhat_max_abs':float(np.max(np.abs(B['raw']-R['raw']))),
                    'chi_max_circular_deg':cdiff(B['chi'],R['chi']),
                    'neighbor_count_blocked':len(B['neighbors']), 'neighbor_count_restored':len(R['neighbors'])},
                'restore_vs_rescue_direct':{
                    'zernike_max_abs':float(np.max(np.abs(R['feature']-P['feature']))),
                    'raw_yhat_max_abs':float(np.max(np.abs(R['raw']-P['raw']))),
                    'chi_max_circular_deg':cdiff(R['chi'],P['chi'])},
                'next_stage_parity':{
                    'blocked_vs_saved_next_stage_chi_deg':cdiff(B['chi'],chi_next_block),
                    'rescue_vs_saved_next_stage_chi_deg':cdiff(P['chi'],chi_next_rescue)},
                'blocked_chi':B['chi'].tolist(),'restored_chi':R['chi'].tolist(),
                'all_full_batch_chi':{'blocked':B['all_chi'],'restored':R['all_chi']},
                'interpretation':'one-step source-sidechain coordinate intervention; not full-feedback new run',
                'Core_STOP':True,'Publication_STOP':True,'MD_executed':False}
            # Both source-specific predictions must reproduce the archived following stage before inference is accepted.
            record['gate']={
                'restoration_matches_rescue_model':record['restore_vs_rescue_direct']['raw_yhat_max_abs']<=1e-4 and record['restore_vs_rescue_direct']['zernike_max_abs']<=1e-5,
                'next_stage_replay_within_1deg':max(record['next_stage_parity'].values())<=1,
                'feature_changed':record['blocked_vs_restored']['zernike_max_abs']>0,
            }
            d=a.output/rid;d.mkdir()
            np.savez_compressed(d/'model_arrays.npz',zernike_blocked=B['feature'],zernike_restored=R['feature'],
                zernike_rescue=P['feature'],raw_blocked=B['raw'],raw_restored=R['raw'],raw_rescue=P['raw'],
                chi_blocked=B['chi'],chi_restored=R['chi'],chi_rescue=P['chi'])
            write_json(d/'neighborhoods.json',{'blocked':B['neighbors'],'restored':R['neighbors'],'rescue':P['neighbors']})
            write_json(d/'result.json',record)
            results.append({'case':case,'stage':stage,'donor':donor,'receiver':receiver,
                            'record':record['blocked_vs_restored'],'gate':record['gate']})
            write_json(a.output/'progress.json',{'completed':len(results),'total':4,'results':results})
        write_json(a.output/'summary.json',{'status':'PASS' if all(all(r['gate'].values()) for r in results) else 'NEEDS_REVIEW',
            'results':results,'Core_STOP':True,'Publication_STOP':True,'MD_executed':False})
    except Exception:
        (a.output/'ERROR.txt').write_text(traceback.format_exc(),encoding='utf-8')
        raise
    finally:
        manifest={str(p.relative_to(a.output)):sha(p.read_bytes()) for p in a.output.rglob('*') if p.is_file() and p.name!='manifest_sha256.json'}
        write_json(a.output/'manifest_sha256.json',manifest)
    print('V8.3 DONE', len(results),'paths',flush=True)
    print('Summary:', a.output/'summary.json',flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--selftest',action='store_true');ap.add_argument('--root',type=Path)
    ap.add_argument('--v8',type=Path);ap.add_argument('--v82',type=Path);ap.add_argument('--output',type=Path)
    ap.add_argument('--confirm',default='')
    a=ap.parse_args()
    if a.selftest:
        require(len(PLAN)==4 and len(set(PLAN))==4,'4 path fixed')
        require(cdiff([179],[-179])==2,'circular distance')
        p=b'ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n'
        require(len(atom_dict(p))==1,'PDB parser')
        print('SELFTEST PASS (no HPacker inference)');return
    require(a.confirm=='I_APPROVE_X4_V8_3_FOUR_PATH_DIAGNOSTIC_ONLY','confirmation missing')
    require(a.root and a.v8 and a.v82 and a.output,'required arguments missing')
    run(a)
if __name__=='__main__':
    try: main()
    except Exception:traceback.print_exc();sys.exit(1)
