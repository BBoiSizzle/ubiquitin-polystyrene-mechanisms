"""Verify public input fingerprints; optional reproduced numerical facts."""
from pathlib import Path
import argparse,hashlib,json,math
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--reproduced',action='store_true');args=ap.parse_args()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((P/'PACKAGE_MANIFEST.json').read_text())
for row in m['files']:
 p=P/row['path'];assert p.is_file(),row['path'];assert p.stat().st_size==row['bytes'],row['path'];assert sha(p)==row['sha256'],row['path']
# Every path the main plotting code can select must resolve inside the package.
for identifier,name in json.loads((P/'INPUT_MAP.json').read_text()).items():
 p=(P/'inputs'/name).resolve();assert p.is_relative_to(P/'inputs') and p.is_file(),identifier
if args.reproduced:
 def same(a,b,key='root'):
  if isinstance(a,dict):
   assert set(a)==set(b),key
   for k,v in a.items():same(v,b[k],key+'.'+k)
  elif isinstance(a,list):
   assert len(a)==len(b),key
   for i,(x,y) in enumerate(zip(a,b)):same(x,y,key+f'[{i}]')
  elif isinstance(a,(int,float)) and not isinstance(a,bool):
   assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),key
  else:assert a==b,key
 same(json.loads((P/'EXPECTED_NUMERICAL_FACTS.json').read_text()),json.loads((P/'FACTS.json').read_text()))
 print('PASS: rebuilt numerical facts agree with the frozen reference (relative/absolute tolerance 1e-12).')
print(f"PASS: {len(m['files'])} packaged files match their SHA-256 manifest.")
