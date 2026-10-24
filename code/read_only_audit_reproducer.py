import zipfile,json,hashlib,re,csv,collections,io,pathlib,statistics
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=pathlib.Path('/mnt/data');out=p/'X4_WT_200_Final_Independent_Audit_2026-10-21';out.mkdir(exist_ok=True)
z=zipfile.ZipFile(p/'x4_wt_ensemble_continuation_through199.zip');seed=zipfile.ZipFile(p/'x4_ensemble_wt_first_batch_results.zip');manifest=json.loads((p/'x4_wt_ensemble_through199_delivery_manifest.json').read_text());assert hashlib.sha256((p/'x4_wt_ensemble_continuation_through199.zip').read_bytes()).hexdigest()==manifest['archive_sha256'];assert z.testzip() is None
n=0
for k,v in manifest['file_sha256'].items():
 assert hashlib.sha256(z.read(k)).hexdigest()==v,(k,'hash');n+=1
assert n==len(z.namelist())==3659
files={}
for zz in [seed,z]:
 for k in zz.namelist():
  t=k.replace('\\','/');m=re.search(r'WT_source(\d{3})/(.+)$',t)
  if m:files[(int(m.group(1)),m.group(2))]=zz.read(k)
assert sorted({i for i,f in files})==list(range(200))
# verify seed receipt and 8..199 receipt (seed receipt may format differ)
rv=0
for i in range(200):
 receipt=json.loads(files[(i,'receipt.json')]);assert 'PASS' in receipt['status']
 for fn,h in receipt['file_sha256'].items():
  assert hashlib.sha256(files[i,fn]).hexdigest()==h,(i,fn);rv+=1
# supplementary manifest seed digest
assert hashlib.sha256((p/'x4_ensemble_wt_first_batch_results.zip').read_bytes()).hexdigest()==manifest['seed_sha256']
rows=[]; events=[]; summaries=[]; stagekeys=[]
for i in range(200):
 if i<8:
  old=json.loads(files[i,'stage_diagnostics.json'])
  g=[];phe=[]
  for st, entry in enumerate(old):
   a=entry['geometry']; oldcs=a['contacts_below_1p6A']; cs=[{'residue_1':c['atom_1'][0],'atom_1':c['atom_1'][1],'residue_2':c['atom_2'][0],'atom_2':c['atom_2'][1],'distance_A':c['distance_A']} for c in oldcs]
   g.append({'stage':st,'min_nonadjacent_distance_A':a['min_nonadjacent_residue_distance_A'],'n_pairs_below_1p6A':a['n_pairs_below_1p6A'],'n_unique_residue_pairs_below_1p6A':a['n_unique_residue_pairs_below_1p6A'],'contacts_below_1p6A':cs})
   pdb=files[i,f'stage_{st:02d}.pdb'].decode()
   atoms={6:[],55:[]}
   for line in pdb.splitlines():
    if line.startswith(('ATOM  ','HETATM')):
     resid=int(line[22:26]); atom=line[12:16].strip(); element=line[76:78].strip()
     if resid in atoms and element not in ('H','D'):
      atoms[resid].append((atom,np.array([float(line[k:k+8]) for k in (30,38,46)])))
   pairs=sorted([(float(np.linalg.norm(v1-v2)),an1,an2) for an1,v1 in atoms[6] for an2,v2 in atoms[55]])
   d0,an1,an2=pairs[0]
   phe.append({'min_distance_A':d0,'closest_PHE6_atom':an1,'closest_PHE55_atom':an2,'n_PHE6_PHE55_pairs_lt_1p6':sum(x[0]<1.6 for x in pairs),'min_ring_to_ring_A':min(x[0] for x in pairs if x[1] in {'CG','CD1','CD2','CE1','CE2','CZ'} and x[2] in {'CG','CD1','CD2','CE1','CE2','CZ'})})
 else:
  g=json.loads(files[i,'stage_general_diagnostics.json']);phe=json.loads(files[i,'PHE6_55_stages.json'])
 assert len(g)==len(phe)==6
 contacts=[];d=[];bit=[]
 for s,(a,b) in enumerate(zip(g,phe)):
  if isinstance(a,dict):pass
  else:raise Exception((i,s,a))
  ds=float(b['min_distance_A']); d.append(ds);bit.append(int(ds<1.6))
  ct=a['contacts_below_1p6A'];cs={tuple((v['residue_1'],v['atom_1'],v['residue_2'],v['atom_2'])) for v in ct};contacts.append(cs)
  rows.append({'frame':i,'stage':s,'general_min_distance_A':a['min_nonadjacent_distance_A'],'general_pairs_lt_1p6':a['n_pairs_below_1p6A'],'general_residue_pairs_lt_1p6':a['n_unique_residue_pairs_below_1p6A'],'phe655_min_distance_A':ds,'phe655_n_pairs_lt_1p6':b['n_PHE6_PHE55_pairs_lt_1p6'],'phe655_atom6':b['closest_PHE6_atom'],'phe655_atom55':b['closest_PHE55_atom'],'phe655_ring_min_A':b['min_ring_to_ring_A']})
  if s:
   events.append({'frame':i,'stage_from':s-1,'stage_to':s,'new_contact_atoms':len(cs-contacts[s-1]),'lost_contact_atoms':len(contacts[s-1]-cs),'persist_contact_atoms':len(cs&contacts[s-1])})
 reappear=any(bool(contacts[s] & set().union(*contacts[:s-1])) for s in range(2,6))
 patt=''.join(map(str,bit))
 # robust 1,0,1 subsequence
 phe_reappear=any(bit[a]==1 and bit[b]==0 and bit[c]==1 for a in range(4) for b in range(a+1,5) for c in range(b+1,6))
 summaries.append({'frame':i,'phe655_pattern':patt,'phe655_recontact':int(phe_reappear),'phe655_final_close':bit[-1],'phe655_any_close':int(any(bit)),'phe655_final_distance_A':d[-1],'phe655_min_distance_A':min(d),'phe655_range_A':max(d)-min(d),'general_final_close':int(g[-1]['n_pairs_below_1p6A']>0),'general_any_close':int(any(x['n_pairs_below_1p6A']>0 for x in g)),'general_same_atom_recontact':int(reappear),'general_final_min_distance_A':g[-1]['min_nonadjacent_distance_A'],'general_final_pair_count':g[-1]['n_pairs_below_1p6A']})
def save_csv(name,rs):
 with open(out/name,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
save_csv('all_1200_stages.csv',rows);save_csv('frame_summary_200.csv',summaries);save_csv('general_contact_transitions.csv',events)
arr=np.array([[r['phe655_min_distance_A'] for r in rows[i*6:(i+1)*6]] for i in range(200)])
fig,ax=plt.subplots(figsize=(10,12));im=ax.imshow(np.minimum(arr,5),aspect='auto',vmin=0,vmax=5,cmap='viridis');ax.set_xlabel('Stage');ax.set_ylabel('WT source frame');ax.set_xticks(range(6));fig.colorbar(im,ax=ax,label='PHE 6–55 min heavy-atom distance (Å; clipped at 5)');fig.tight_layout();fig.savefig(out/'phe655_200_heatmap.png',dpi=165);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));ax.plot(range(6),[sum(rows[i*6+s]['general_pairs_lt_1p6']>0 for i in range(200)) for s in range(6)],'o-',label='Any general <1.6 Å');ax.plot(range(6),[sum(rows[i*6+s]['phe655_min_distance_A']<1.6 for i in range(200)) for s in range(6)],'s-',label='PHE 6–55 <1.6 Å');ax.set_xlabel('Stage');ax.set_ylabel('Frames / 200');ax.legend();fig.tight_layout();fig.savefig(out/'stage_prevalence.png',dpi=165);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));ax.hist([r['general_final_min_distance_A'] for r in summaries],bins=28,alpha=.65,label='Global min');ax.hist([r['phe655_final_distance_A'] for r in summaries],bins=28,alpha=.55,label='PHE6–55 min');ax.set_xlabel('Stage 5 distance (Å)');ax.set_ylabel('Frames');ax.legend();fig.tight_layout();fig.savefig(out/'final_distance_distributions.png',dpi=165);plt.close(fig)
patterns=collections.Counter(x['phe655_pattern'] for x in summaries)
result={'zip_entries_verified':n,'receipt_file_sha256_verified':rv,'frames':200,'stage_structures':len(rows),'general_final_close':sum(s['general_final_close'] for s in summaries),'general_any_close':sum(s['general_any_close'] for s in summaries),'general_same_atom_recontact':sum(s['general_same_atom_recontact'] for s in summaries),'phe_final_close':sum(s['phe655_final_close'] for s in summaries),'phe_any_close':sum(s['phe655_any_close'] for s in summaries),'phe_recontact':sum(s['phe655_recontact'] for s in summaries),'phe_010101':patterns['010101'],'phe_101010':patterns['101010'],'phe_recontact_frames':[s['frame'] for s in summaries if s['phe655_recontact']],'stage_general_close':[sum(rows[i*6+j]['general_pairs_lt_1p6']>0 for i in range(200)) for j in range(6)],'stage_phe_close':[sum(rows[i*6+j]['phe655_min_distance_A']<1.6 for i in range(200)) for j in range(6)],'stage_general_min_median':[float(np.median([rows[i*6+j]['general_min_distance_A'] for i in range(200)])) for j in range(6)],'patterns_top':patterns.most_common(12),'min_phe_final':min(summaries,key=lambda s:s['phe655_final_distance_A']),'min_general_final':min(summaries,key=lambda s:s['general_final_min_distance_A'])}
(out/'AUDIT_SUMMARY.json').write_text(json.dumps(result,indent=2,ensure_ascii=False));print(json.dumps(result,ensure_ascii=False,indent=2))
