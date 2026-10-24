from pathlib import Path
import csv,collections,hashlib
r=Path(__file__).resolve().parents[1]
with open(r/'data'/'all_1200_stages.csv') as f:s=list(csv.DictReader(f))
with open(r/'data'/'frame_summary_200.csv') as f:p=list(csv.DictReader(f))
assert len(s)==1200 and len(p)==200
assert [sum(int(x['general_pairs_lt_1p6'])>0 for x in s if int(x['stage'])==i) for i in range(6)]==[183,148,140,128,130,134]
assert [sum(float(x['phe655_min_distance_A'])<1.6 for x in s if int(x['stage'])==i) for i in range(6)]==[39,22,25,23,29,21]
assert sum(int(x['phe655_recontact']) for x in p)==43
c=collections.Counter(x['phe655_pattern'] for x in p)
assert (c['010101'],c['101010'])==(10,15)
for line in (r/'SHA256SUMS.txt').read_text().splitlines():
 if not line.strip():continue
 sha,rel=line.split('  ',1);assert hashlib.sha256((r/rel).read_bytes()).hexdigest()==sha,rel
print('PASS: 200 frames, 1200 stages, 43 recontacts, 25 alternations, SHA256 files')
