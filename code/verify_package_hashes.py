#!/usr/bin/env python3
"""Verify a package SHA256SUMS.txt without executing any scientific code."""
from pathlib import Path,PurePosixPath
import hashlib,sys,json
root=Path(sys.argv[1] if len(sys.argv)>1 else '.').resolve();manifest=root/'SHA256SUMS.txt';count=0;listed=set();fail=[]
for line in manifest.read_text().splitlines():
 h,name=line.split('  ',1);p=PurePosixPath(name)
 if len(h)!=64 or p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name:raise ValueError('Unsafe manifest entry: '+name)
 q=root/p
 if q.is_symlink() or not q.is_file() or hashlib.sha256(q.read_bytes()).hexdigest()!=h:fail.append(name)
 listed.add(name);count+=1
extras=[str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and str(p.relative_to(root)) not in listed and p!=manifest]
print(json.dumps(dict(status='FAIL' if fail else 'PASS',entries=count,mismatches=fail,unlisted_extra_files=extras),indent=2));sys.exit(bool(fail))
