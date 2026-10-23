#!/usr/bin/env python3
"""Restore exact reviewed files from hash-pinned binary parts; standard library only."""
from pathlib import Path,PurePosixPath
import hashlib,json,os,tempfile
H=Path(__file__).resolve().parent

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def safe(name):
 q=PurePosixPath(name)
 if not name or q.is_absolute() or '..' in q.parts or '\\' in name or ':' in name or str(q)!=name:raise ValueError('Unsafe path: '+name)
 p=H.joinpath(*q.parts)
 for a in [p,*p.parents]:
  if a==H:break
  if a.is_symlink():raise ValueError('Symlink refused: '+name)
 if not p.resolve().is_relative_to(H):raise ValueError('Path escapes repository: '+name)
 return p

def checked(p,row):
 return p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256']

def main():
 m=json.loads((H/'LARGE_FILE_TRANSPORT.json').read_text())
 targets=set(); parts=set()
 # Complete preflight before the first write; changed existing targets are never overwritten.
 for row in m['files']:
  target=safe(row['path'])
  if row['path'] in targets:raise ValueError('Duplicate target')
  targets.add(row['path'])
  if target.exists() and not checked(target,row):raise ValueError('Refusing changed existing file: '+row['path'])
  digest=hashlib.sha256();count=0
  for part in row['parts']:
   p=safe(part['path'])
   if part['path'] in parts:raise ValueError('Duplicate transport part')
   parts.add(part['path'])
   if not checked(p,part):raise ValueError('Part integrity failure: '+part['path'])
   b=p.read_bytes();digest.update(b);count+=len(b)
  if count!=row['bytes'] or digest.hexdigest()!=row['sha256']:raise ValueError('Restored file digest mismatch: '+row['path'])
 if targets & parts:raise ValueError('Target/part collision')
 restored=0;already=0
 for row in m['files']:
  target=safe(row['path'])
  if checked(target,row):already+=1;continue
  target.parent.mkdir(parents=True,exist_ok=True)
  fd,name=tempfile.mkstemp(prefix='.restore-',dir=target.parent);temporary=Path(name)
  try:
   with os.fdopen(fd,'wb') as f:
    for part in row['parts']:f.write(safe(part['path']).read_bytes())
    f.flush();os.fsync(f.fileno())
   if not checked(temporary,row):raise ValueError('Temporary file integrity failure')
   # Atomic creation with no overwrite; an intervening file creation fails safely.
   os.link(temporary,target);restored+=1
  finally:temporary.unlink(missing_ok=True)
 for row in m['reviewed_unchanged_files']:
  if not checked(safe(row['path']),row):raise ValueError('Reviewed payload mismatch: '+row['path'])
 print(f'PASS: restored {restored} files; {already} already exact; all {len(m["reviewed_unchanged_files"])} unchanged reviewed files verified.')
if __name__=='__main__':main()
