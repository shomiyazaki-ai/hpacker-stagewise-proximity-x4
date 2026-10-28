#!/usr/bin/env python3
"""Read-only X4 geometry verification. No HPacker, NN, OpenMM or MD invocation.
Requires Python >=3.10 and NumPy. Writes only to --out; never modifies input.
"""
from pathlib import Path
import argparse,csv,hashlib,itertools,json,sys
import numpy as np

def parse_pdb(path):
    keys=[]; xyz=[]
    for line in Path(path).read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        element=line[76:78].strip() or line[12:16].strip().lstrip('0123456789')[:1]
        if element in ('H','D'): continue
        keys.append((line[21],int(line[22:26]),line[26],line[17:20].strip(),line[12:16].strip()))
        xyz.append([float(line[i:i+8]) for i in (30,38,46)])
    xyz=np.asarray(xyz,dtype=np.float64)
    if len(keys)!=len(set(keys)) or not np.isfinite(xyz).all():raise ValueError(f'Invalid atoms: {path}')
    return keys,xyz

def geometry(path):
    keys,xyz=parse_pdb(path); n=len(keys);i,j=np.triu_indices(n,1)
    # The fixed single chain has contiguous residue numbers 1..63.
    res=np.asarray([k[1] for k in keys]); chain=np.asarray([k[0] for k in keys]); sel=(chain[i]!=chain[j]) | (np.abs(res[i]-res[j])>1)
    i,j=i[sel],j[sel];dist=np.sqrt(np.sum((xyz[i]-xyz[j])**2,axis=1))
    contacts={tuple(sorted((keys[a],keys[b]))) for a,b,d in zip(i,j,dist) if d<1.6}
    p6=[a for a,k in enumerate(keys) if k[1]==6 and k[3]=='PHE'];p55=[a for a,k in enumerate(keys) if k[1]==55 and k[3]=='PHE']
    pp=sorted((float(np.linalg.norm(xyz[a]-xyz[b])),keys[a][-1],keys[b][-1]) for a in p6 for b in p55)
    ring={'CG','CD1','CD2','CE1','CE2','CZ'}
    row=dict(general_min_distance_A=float(dist.min()),general_pairs_lt_1p6=len(contacts),general_residue_pairs_lt_1p6=len({(a[:4],b[:4]) for a,b in contacts}),phe655_min_distance_A=pp[0][0],phe655_n_pairs_lt_1p6=sum(d<1.6 for d,a,b in pp),phe655_atom6=pp[0][1],phe655_atom55=pp[0][2],phe655_ring_min_A=min(d for d,a,b in pp if a in ring and b in ring))
    return row,contacts,keys

def ordered(bits):
    return any(bits[i] and not bits[j] and bits[k] for i,j,k in itertools.combinations(range(6),3))
def strict_same_atom(contacts):
    # A given unordered atom pair must be present, absent, then present.
    return any(bool((contacts[i]&contacts[k])-contacts[j]) for i,j,k in itertools.combinations(range(6),3))
def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--pdb-root',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);pdbroot=a.pdb_root or a.root/'data/pdb'
    expected=json.loads((a.root/'data/expected_e1.json').read_text());manifest=list(csv.DictReader((a.root/'provenance/PDB_SHA256_Manifest.csv').open()))
    assert len(manifest)==1200 and len({(int(r['source']),int(r['stage'])) for r in manifest})==1200
    actual=set(pdbroot.glob('WT_source*/stage_*.pdb'));want={pdbroot/f'WT_source{f:03d}'/f'stage_{s:02d}.pdb' for f in range(200) for s in range(6)};assert actual==want
    m={(int(r['source']),int(r['stage'])):r for r in manifest};stages=[];frames=[];trans=[];ref=None;receipt_matches=0
    for frame in range(200):
        cs=[];bits=[];rs=[]
        for stage in range(6):
            p=pdbroot/f'WT_source{frame:03d}'/f'stage_{stage:02d}.pdb';h=hashlib.sha256(p.read_bytes()).hexdigest();r=m[frame,stage];assert h==r['sha256'] and p.stat().st_size==int(r['bytes'])
            receipt_path=a.root/r['receipt_path']
            if receipt_path.exists():
                rec=json.loads(receipt_path.read_text());assert rec['frame']==frame and rec['file_sha256'][p.name]==h;receipt_matches+=1
            g,c,keys=geometry(p);ref=keys if ref is None else ref;assert keys==ref and len(keys)==494
            stages.append(dict(frame=frame,stage=stage,**g));rs.append(g);cs.append(c);bits.append(g['phe655_n_pairs_lt_1p6']>0)
        frames.append(dict(frame=frame,phe655_pattern=''.join('1' if b else '0' for b in bits),phe655_recontact=int(ordered(bits)),general_strict_recontact=int(strict_same_atom(cs)),general_any=int(any(cs)),phe_any=int(any(bits))))
        for s in range(1,6):trans.append(dict(frame=frame,stage_from=s-1,stage_to=s,new_contact_atoms=len(cs[s]-cs[s-1]),lost_contact_atoms=len(cs[s-1]-cs[s]),persist_contact_atoms=len(cs[s]&cs[s-1])))
    summary=dict(structures=len(stages),atoms_each=len(ref),general_stage_counts=[sum(r['general_pairs_lt_1p6']>0 for r in stages if r['stage']==s) for s in range(6)],phe_stage_counts=[sum(r['phe655_n_pairs_lt_1p6']>0 for r in stages if r['stage']==s) for s in range(6)],general_union=sum(r['general_any'] for r in frames),phe_union=sum(r['phe_any'] for r in frames),strict_same_atom_recontact=sum(r['general_strict_recontact'] for r in frames),phe_residue_pair_recontact=sum(r['phe655_recontact'] for r in frames),alternation_010101=[r['frame'] for r in frames if r['phe655_pattern']=='010101'],alternation_101010=[r['frame'] for r in frames if r['phe655_pattern']=='101010'])
    assert summary==expected,{'observed':summary,'expected':expected}
    old=list(csv.DictReader((a.root/'data/historical_v1.0.0/all_1200_stages.csv').open()));maxerr=0
    for oldr,new in zip(sorted(old,key=lambda r:(int(r['frame']),int(r['stage']))),stages):
        for k,v in oldr.items():
            if k in ['phe655_atom6','phe655_atom55']:assert new[k]==v
            else:
                err=abs(float(v)-float(new[k]));maxerr=max(maxerr,err);assert err<1e-9,(k,oldr,new)
    # Definition controls exclude persistence and different-pair substitution.
    pair=('a','b');other=('c','d')
    assert not strict_same_atom([{pair}]*6)
    assert strict_same_atom([{pair},set(),{pair},set(),set(),set()])
    assert not strict_same_atom([{pair},set(),{other},set(),set(),set()])
    assert ordered([1,0,1,0,0,0]) and not ordered([1]*6)
    write_csv(a.out/'e1_stages.csv',stages);write_csv(a.out/'e1_frames.csv',frames);write_csv(a.out/'e1_transitions.csv',trans)
    summary.update(status='PASS',sha_matches=1200,receipt_sha_matches=receipt_matches,archived_csv_max_abs_error=maxerr,definition_controls='PASS')
    (a.out/'e1_verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
