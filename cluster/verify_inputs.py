"""Read and hash every transferred model without invoking MQC."""
import argparse
import gzip
import hashlib
from pathlib import Path
from common import load_inputs, sha, under

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--inputs',type=Path,required=True)
args=p.parse_args()
rows,meta=load_inputs(args.inputs)
for i,r in enumerate(rows,1):
    path=under(args.inputs,r['input_relative_path'])
    if sha(path)!=r['gzip_sha256']: raise ValueError('File mismatch: '+r['job_id'])
    h=hashlib.sha256()
    with gzip.open(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
    if h.hexdigest()!=r['xml_sha256']: raise ValueError('XML mismatch: '+r['job_id'])
    if i%1000==0: print('Verified',i,flush=True)
print('Verified all',len(rows),'models; database counts:',meta['database_counts'])
