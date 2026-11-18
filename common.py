from pathlib import Path
import hashlib, json, urllib.request
from importlib.metadata import version

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def download(url,path,expected_sha256=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        temporary=path.with_name(path.name+'.part')
        with urllib.request.urlopen(url,timeout=120) as response,temporary.open('wb') as stream:
            for block in iter(lambda:response.read(1024*1024),b''):stream.write(block)
        temporary.replace(path)
    digest=sha256(path)
    if expected_sha256 and digest!=expected_sha256:raise ValueError(f'SHA256 mismatch: {path.name}')
    return {'file':path.name,'url':url,'sha256':digest,'bytes':path.stat().st_size}

def library(data_dir,name):
    import gseapy as gp
    path=Path(data_dir)/(name+'.json')
    if path.exists():members=json.loads(path.read_text())
    else:
        members=gp.get_library(name=name,organism='Human')
        path.write_text(json.dumps(members,sort_keys=True)+'\n')
    return members,{'library':name,'sha256':sha256(path),'provider':'Enrichr via GSEApy','source':'https://maayanlab.cloud/Enrichr/'}

def run_gsea(rank,members):
    import gseapy as gp
    # Descending rank; mergesort makes equal-rank ordering stable.
    rank=rank.dropna().sort_values(ascending=False,kind='mergesort')
    return gp.prerank(rnk=rank,gene_sets=members,permutation_num=1000,seed=42,
        min_size=15,max_size=500,weight=1,threads=2,outdir=None,no_plot=True,
        verbose=False).res2d.sort_values('FDR q-val').reset_index(drop=True)

def versions(names):return {name:version(name) for name in names}
