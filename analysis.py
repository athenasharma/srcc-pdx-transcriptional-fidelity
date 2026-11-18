"""Descriptive transcriptional comparison of six deposited SRCC samples.

Run: python analysis.py --data-dir data --output-dir results/recomputed
Existing manuscript figures/results are retained separately from this rerun.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from common import download,library,run_gsea,versions,sha256

SAMPLES=['N_1','N_2','T_1','T_2','UT_1','UT_2']
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279979/suppl/GSE279979_raw_counts.txt.gz'
SOURCE_SHA='686a7f6304a7294474e5e6f208cddba7e241d8549adc44ed7f99d120847aedc5'

def run(data_dir='data',output_dir='results/recomputed'):
    data_dir=Path(data_dir);out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    source=download(URL,data_dir/'GSE279979_raw_counts.txt.gz',SOURCE_SHA)
    raw=pd.read_csv(data_dir/'GSE279979_raw_counts.txt.gz',sep='\t')
    assert set(SAMPLES).issubset(raw.columns)
    raw['Geneid']=raw.Geneid.str.split('.').str[0]
    counts=raw.set_index('Geneid')[SAMPLES]
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all()
    assert np.equal(counts.to_numpy(),np.floor(counts.to_numpy())).all()
    meta=pd.DataFrame({'sample':SAMPLES,'group':['Normal']*2+['Tumour']*2+['PDX']*2,
        'GEO_accession':['GSM8584545','GSM8584546','GSM8584547','GSM8584548','GSM8584549','GSM8584550']})
    meta.to_csv(out/'sample_manifest.csv',index=False)
    # Suffixes are sample labels, not evidence of independent matched donors.
    counts=counts.loc[(counts>=10).sum(axis=1)>=2]
    # Low-count filtering removes the duplicated X/Y pseudoautosomal IDs in this source.
    assert counts.index.is_unique,'Retained duplicate gene IDs require an explicit aggregation rule.'
    log_cpm=np.log2(counts.div(counts.sum(axis=0),axis=1)*1e6+1)
    log_cpm.to_csv(out/'log2_CPM_plus1.csv')
    pca=PCA(n_components=2);coordinates=pca.fit_transform(log_cpm.T)
    pd.DataFrame(coordinates,index=SAMPLES,columns=['PC1','PC2']).to_csv(out/'PCA.csv')
    log_cpm.corr(method='spearman').to_csv(out/'sample_correlations.csv')
    normal=log_cpm[['N_1','N_2']].mean(axis=1);tumour=log_cpm[['T_1','T_2']].mean(axis=1)
    pdx=log_cpm[['UT_1','UT_2']].mean(axis=1)
    effects=pd.DataFrame({'Tumor_vs_Normal':tumour-normal,'PDX_vs_Normal':pdx-normal})
    effects.to_csv(out/'gene_effects.csv')
    direct=pdx-tumour;direct.to_csv(out/'direct_PDX_vs_tumour.csv',header=['mean_log2_CPM_plus1_difference'])
    centroid={sample:{'tumour':float(spearmanr(log_cpm[sample],tumour).statistic),
        'normal':float(spearmanr(log_cpm[sample],normal).statistic)} for sample in ['UT_1','UT_2']}
    mapping_path=data_dir/'ensembl_symbol_mapping.csv'
    if not mapping_path.exists():
        import mygene
        mapping=mygene.MyGeneInfo().querymany(effects.index.tolist(),scopes='ensembl.gene',
            fields='symbol',species='human',as_dataframe=True,verbose=False)
        mapping=mapping[['symbol']].dropna();mapping=mapping.loc[~mapping.index.duplicated(keep='first')]
        mapping.to_csv(mapping_path,index_label='Geneid')
    mapping=pd.read_csv(mapping_path,index_col=0)
    mapped=effects.join(mapping[['symbol']],how='inner').dropna(subset=['symbol'])
    # Freeze the final historical notebook's rule: largest absolute tumour-normal
    # effect per symbol. This is a rank-construction choice, not gene-level inference.
    mapped=mapped.sort_values('Tumor_vs_Normal',key=lambda x:x.abs(),ascending=False,kind='mergesort')
    mapped=mapped.drop_duplicates('symbol').set_index('symbol')
    mapped.to_csv(out/'mapped_rank_effects.csv')
    hallmark,library_info=library(data_dir,'MSigDB_Hallmark_2020')
    tg=run_gsea(mapped.Tumor_vs_Normal,hallmark);pg=run_gsea(mapped.PDX_vs_Normal,hallmark)
    dg=run_gsea(mapped.PDX_vs_Normal-mapped.Tumor_vs_Normal,hallmark)
    tg.to_csv(out/'tumour_Hallmark_GSEA.csv',index=False);pg.to_csv(out/'PDX_Hallmark_GSEA.csv',index=False)
    dg.to_csv(out/'direct_Hallmark_GSEA.csv',index=False)
    c=tg[['Term','NES','FDR q-val']].merge(pg[['Term','NES','FDR q-val']],on='Term',suffixes=('_tumor','_pdx'))
    c.to_csv(out/'Hallmark_comparison.csv',index=False)
    joint=c[(c['FDR q-val_tumor']<.05)&(c['FDR q-val_pdx']<.05)&(np.sign(c.NES_tumor)==np.sign(c.NES_pdx))]
    joint.to_csv(out/'jointly_significant_same_direction.csv',index=False)
    c[(c['FDR q-val_tumor']<.05)&(c['FDR q-val_pdx']>=.05)].to_csv(out/'tumour_normal_only_significant.csv',index=False)
    c[(c['FDR q-val_tumor']>=.05)&(c['FDR q-val_pdx']<.05)].to_csv(out/'PDX_normal_only_significant.csv',index=False)
    # An FDR threshold difference is not a differential pathway test or equivalence test.
    summary={'n_samples':6,'retained_genes':len(effects),'gene_effect_spearman':float(spearmanr(effects.Tumor_vs_Normal,effects.PDX_vs_Normal).statistic),
        'pathway_NES_spearman':float(spearmanr(c.NES_tumor,c.NES_pdx).statistic),'jointly_significant_same_direction':len(joint),
        'direct_Hallmark_FDR_hits':int((dg['FDR q-val']<.05).sum()),'median_absolute_direct_difference':float(direct.abs().median()),
        'direct_difference_gt1':int((direct.abs()>1).sum()),'centroid_correlations':centroid,
        'inference':'Descriptive sample-level comparison; GSEA uses gene-set permutations, not a population-level model-fidelity test.'}
    fig,ax=plt.subplots(figsize=(7,5));ax.scatter(effects.Tumor_vs_Normal,effects.PDX_vs_Normal,s=5,alpha=.2)
    ax.set(xlabel='Primary tumour - normal: mean log2(CPM+1)',ylabel='PDX - normal: mean log2(CPM+1)',title=f"Expression-effect concordance (Spearman r = {summary['gene_effect_spearman']:.3f})")
    fig.tight_layout();fig.savefig(out/'gene_effect_concordance.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,5));ax.scatter(c.NES_tumor,c.NES_pdx,s=25)
    ax.set(xlabel='Primary tumour - normal: Hallmark NES',ylabel='PDX - normal: Hallmark NES',title=f"Pathway-score concordance (Spearman r = {summary['pathway_NES_spearman']:.3f})")
    fig.tight_layout();fig.savefig(out/'pathway_score_concordance.png',dpi=180);plt.close(fig)
    provenance={'source':source,'symbol_mapping_sha256':sha256(mapping_path),'gene_sets':library_info,
        'versions':versions(['numpy','pandas','scipy','scikit-learn','gseapy','mygene','matplotlib']),
        'method':'Library-size normalized counts; mean log2(CPM+1) contrasts; 1000 gene-set permutations, seed 42.'}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data-dir',default='data');p.add_argument('--output-dir',default='results/recomputed');a=p.parse_args();run(a.data_dir,a.output_dir)
