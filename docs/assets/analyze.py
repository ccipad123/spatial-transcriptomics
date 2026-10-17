"""Reproducible, CPU-only analysis of two public DLPFC Visium sections.
Run from any directory. Fixed design: 151673 development, 151674 held-out section.
No coordinates, author clusters, or labels enter gene selection or PCA.
"""
from pathlib import Path
import json, sys, platform, hashlib, importlib.metadata, time
import numpy as np
import pandas as pd
import h5py
from scipy import sparse, stats
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score,
    roc_auc_score, average_precision_score, log_loss, confusion_matrix,
    classification_report, adjusted_rand_score, normalized_mutual_info_score)
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; DER=ROOT/'data/derived'
SEED=20261012
LAYERS=['L1','L2','L3','L4','L5','L6','WM']
MARKERS=['RELN','CUX2','RORB','PCP4','NTNG2','MBP','MOBP','SNAP25','AQP4','GFAP','PLP1','SLC17A7']

def dump(obj,path):
    def conv(x):
        if isinstance(x,(np.integer,np.floating)):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        raise TypeError(type(x))
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=conv))

def read_section(sample):
    d=ROOT/'data/raw'/sample
    with h5py.File(d/'counts.h5','r') as h:
        m=h['matrix']
        X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr().astype(np.float32)
        genes=np.array([x.decode() for x in m['features']['name'][:]])
        ids=np.array([x.decode() for x in m['features']['id'][:]])
        barcodes=[x.decode() for x in m['barcodes'][:]]
    assert len(set(barcodes))==len(barcodes)
    positions=pd.read_csv(d/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
    assert positions.index.is_unique
    obs=positions.loc[barcodes].copy()
    labels=pd.read_csv(ROOT/'data/raw/barcode_level_layer_map.tsv',sep='\t',header=None,names=['barcode','sample','layer'])
    labels=labels[labels['sample'].astype(str)==sample].set_index('barcode')
    assert labels.index.is_unique
    obs['layer']=labels['layer'].reindex(obs.index)
    meta=pd.read_csv(d/'metadata.tsv',sep='\t',index_col=0)
    assert meta.index.is_unique
    common=obs.index.intersection(meta.index)
    assert (obs.loc[common,'layer'].fillna('NA')==meta.loc[common,'layer_guess_reordered_short'].fillna('NA')).all(), 'Primary/mirror label mismatch'
    assert np.allclose(np.asarray(X.sum(1)).ravel(),meta.loc[obs.index,'sum_umi']), 'Count and barcode alignment mismatch'
    scale=json.loads((d/'scalefactors_json.json').read_text())
    obs['x']=obs.pixel_col*scale['tissue_hires_scalef'];obs['y']=obs.pixel_row*scale['tissue_hires_scalef']
    # Independent coordinate cross-check against the mirror's 600px image coordinates.
    assert np.allclose(obs.x,meta.loc[obs.index,'imagecol']*scale['tissue_hires_scalef']/scale['tissue_lowres_scalef'],atol=.02)
    assert np.allclose(obs.y,meta.loc[obs.index,'imagerow']*scale['tissue_hires_scalef']/scale['tissue_lowres_scalef'],atol=.02)
    obs['sample']=sample
    obs['umi']=np.asarray(X.sum(1)).ravel()
    obs['n_genes']=np.diff(X.indptr)
    mt=np.char.startswith(genes.astype(str),'MT-')
    obs['pct_mt']=100*np.asarray(X[:,mt].sum(1)).ravel()/np.maximum(obs.umi,1)
    obs['qc_pass']=(obs.in_tissue==1)&(obs.umi>=500)&(obs.n_genes>=250)&(obs.pct_mt<=30)
    obs['eligible']=obs.qc_pass&obs.layer.isin(LAYERS)
    obs['block']=(obs.array_row//10).astype(str)+'_'+(obs.array_col//20).astype(str)
    obs.to_csv(DER/f'{sample}_all_spots.csv')
    info={'sample':sample,'n_spots':X.shape[0],'n_genes':X.shape[1],'n_annotated':int(obs.layer.notna().sum()),
          'n_qc_pass':int(obs.qc_pass.sum()),'n_eligible':int(obs.eligible.sum()),
          'median_umi':float(obs.umi.median()),'median_genes':float(obs.n_genes.median()),'median_pct_mt':float(obs.pct_mt.median()),
          'subject':str(meta.subject.iloc[0]),'position':int(meta.position.iloc[0]),'replicate':int(meta.replicate.iloc[0]),
          'nnz':X.nnz,'sparsity':1-X.nnz/np.prod(X.shape),'scale':scale,'missing_labels':int(obs.layer.isna().sum())}
    return X,obs,genes,ids,info

def lognorm(X):
    total=np.asarray(X.sum(1)).ravel()
    Y=sparse.diags(1e4/np.maximum(total,1))@X
    Y=Y.tocsr();Y.data=np.log1p(Y.data)
    return Y

class GeneProjector(BaseEstimator,TransformerMixin):
    """Fit-only-on-training, mean-binned dispersion -> scaling -> 50-PC pipeline."""
    def __init__(self,genes,n_hvg=2000,n_pcs=50):
        self.genes=genes;self.n_hvg=n_hvg;self.n_pcs=n_pcs
    def fit(self,X,y=None):
        mean=np.asarray(X.mean(0)).ravel();var=np.asarray(X.power(2).mean(0)).ravel()-mean**2
        detection=np.asarray((X>0).sum(0)).ravel()
        allowed=(detection>=20)&(~np.char.startswith(self.genes.astype(str),'MT-'))&(~np.char.startswith(self.genes.astype(str),'RPL'))&(~np.char.startswith(self.genes.astype(str),'RPS'))
        ix=np.flatnonzero(allowed)
        bins=pd.qcut(mean[ix],20,labels=False,duplicates='drop')
        disp=var[ix]/np.maximum(mean[ix],1e-8)
        score=np.zeros(len(ix))
        for b in np.unique(bins):
            hit=np.asarray(bins==b);score[hit]=(disp[hit]-disp[hit].mean())/max(disp[hit].std(),1e-8)
        self.selected_=ix[np.argsort(score,kind='stable')[-self.n_hvg:]]
        self.scaler_=StandardScaler().fit(X[:,self.selected_].toarray())
        Z=np.clip(self.scaler_.transform(X[:,self.selected_].toarray()),-10,10)
        self.pca_=PCA(n_components=self.n_pcs,svd_solver='randomized',random_state=SEED).fit(Z)
        self.allowed_count_=len(ix)
        return self
    def transform(self,X):
        Z=np.clip(self.scaler_.transform(X[:,self.selected_].toarray()),-10,10)
        return self.pca_.transform(Z)

def graph(obs):
    xy=obs[['x','y']].to_numpy()
    tree=cKDTree(xy);dist,_=tree.query(xy,k=2)
    radius=float(np.median(dist[:,1])*1.15)
    pairs=tree.query_pairs(radius,output_type='ndarray')
    row=np.r_[pairs[:,0],pairs[:,1]];col=np.r_[pairs[:,1],pairs[:,0]]
    W=sparse.csr_matrix((np.ones(len(row)),(row,col)),shape=(len(obs),len(obs)))
    return W,radius

def moran(v,W,rng,n_perm=999):
    z=v-v.mean();den=float(z@z)
    if den==0:return np.nan,np.nan
    c=len(v)/W.sum()/den
    observed=float(c*(z@(W@z)))
    sim=np.empty(n_perm)
    for i in range(n_perm):
        zz=rng.permutation(z);sim[i]=c*(zz@(W@zz))
    return observed,float((1+(sim>=observed).sum())/(n_perm+1))

def metrics(y,p,prob):
    Y=label_binarize(y,classes=np.arange(7))
    return {'accuracy':accuracy_score(y,p),'balanced_accuracy':balanced_accuracy_score(y,p),
            'macro_f1':f1_score(y,p,labels=np.arange(7),average='macro',zero_division=0),
            'macro_auroc':roc_auc_score(Y,prob,average='macro'),
            'macro_ap':average_precision_score(Y,prob,average='macro'),
            'log_loss':log_loss(y,prob,labels=np.arange(7)),
            'brier_multiclass':float(np.mean(np.sum((Y-prob)**2,axis=1)))}

def new_models():
    return {'Majority':DummyClassifier(strategy='most_frequent'),
            'Logistic':LogisticRegression(C=1.0,solver='lbfgs',max_iter=2000,class_weight='balanced',random_state=SEED),
            'RandomForest':RandomForestClassifier(n_estimators=300,min_samples_leaf=3,max_features='sqrt',class_weight='balanced_subsample',n_jobs=2,random_state=SEED)}

def main():
    t0=time.time();OUT.mkdir(exist_ok=True);DER.mkdir(exist_ok=True)
    sections={};summ=[]
    for s in ['151673','151674']:
        X,o,g,ids,info=read_section(s);summ.append(info)
        keep=o.eligible.to_numpy();o=o.loc[keep].copy();L=lognorm(X[keep]);W,radius=graph(o)
        sections[s]={'X':X[keep],'obs':o,'L':L,'W':W,'genes':g,'ids':ids}
        info.update(graph_edges=W.nnz//2,graph_radius=radius,isolated_spots=int((np.asarray(W.sum(1)).ravel()==0).sum()))
        print('READ',info,flush=True)
    a,b=sections.values();assert np.array_equal(a['ids'],b['ids'])
    gene_index={g:i for i,g in enumerate(a['genes'])}
    # Quality-control retention table uses all spots, including missing reference labels.
    retention=[]
    for s in sections:
        o=pd.read_csv(DER/f'{s}_all_spots.csv')
        for layer,sub in o.groupby(o.layer.fillna('Unlabelled')):
            retention.append({'sample':s,'layer':layer,'before':len(sub),'qc_pass':int(sub.qc_pass.sum()),'eligible':int(sub.eligible.sum()),'retention':float(sub.qc_pass.mean())})
    pd.DataFrame(retention).to_csv(OUT/'qc_by_layer.csv',index=False)
    # Prespecified marker-only exploratory statistics; each section is a separate family.
    marker_stats=[];rng=np.random.default_rng(SEED)
    for s,d in sections.items():
        o,L,W=d['obs'],d['L'],d['W']
        marker_values=np.column_stack([L[:,gene_index[g]].toarray().ravel() for g in MARKERS])
        pd.DataFrame(marker_values,index=o.index,columns=MARKERS).to_csv(DER/f'{s}_markers.csv')
        rows=[]
        for j,g in enumerate(MARKERS):
            v=marker_values[:,j].astype(np.float64)
            H,p=stats.kruskal(*[v[o.layer.to_numpy()==layer] for layer in LAYERS])
            mi,mp=moran(v,W,rng)
            rows.append({'sample':s,'gene':g,'H':float(H),'p_kw':float(p),'epsilon2':max(0,float((H-7+1)/(len(v)-7))),
                         'moran_i':mi,'p_moran':mp,'mean':float(v.mean()),'detection_rate':float(np.mean(v>0))})
        for field,out in [('p_kw','q_kw'),('p_moran','q_moran')]:
            qs=multipletests([r[field] for r in rows],method='fdr_bh')[1]
            for r,q in zip(rows,qs):r[out]=float(q)
        marker_stats+=rows
        d['markers']=marker_values
    pd.DataFrame(marker_stats).to_csv(OUT/'marker_statistics.csv',index=False)
    # Example statistical model with grid-clustered covariance; descriptive on one section.
    o=a['obs'].copy();o['mbp']=a['markers'][:,MARKERS.index('MBP')];o['log_umi']=np.log1p(o.umi)
    ols=smf.ols('mbp ~ C(layer, Treatment(reference="L1")) + log_umi + pct_mt',data=o).fit(cov_type='cluster',cov_kwds={'groups':o.block})
    (OUT/'mbp_cluster_robust_ols.txt').write_text(ols.summary().as_text())
    pd.DataFrame({'term':ols.params.index,'coef':ols.params.values,'se':ols.bse.values,'p':ols.pvalues.values,'ci_low':ols.conf_int()[0].values,'ci_high':ols.conf_int()[1].values}).to_csv(OUT/'mbp_cluster_robust_ols.csv',index=False)
    # Three fixed validation scenarios. The cross-section result is the primary result.
    y_a=np.array([LAYERS.index(l) for l in a['obs'].layer]);y_b=np.array([LAYERS.index(l) for l in b['obs'].layer])
    idx=np.arange(len(y_a))
    tr,te=train_test_split(idx,test_size=.25,stratify=y_a,random_state=SEED)
    btr,bte=next(GroupShuffleSplit(n_splits=1,test_size=.25,random_state=SEED).split(idx,y_a,groups=a['obs'].block))
    scenarios=[('Random within section',tr,te,a,y_a[te]),('Spatial block holdout',btr,bte,a,y_a[bte]),('Held-out section 151674',idx,np.arange(len(y_b)),b,y_b)]
    results=[];preds={};train_sizes=[]
    for scenario,itr,ite,test_section,truth in scenarios:
        projector=GeneProjector(a['genes']).fit(a['L'][itr]);ztr=projector.transform(a['L'][itr]);zte=projector.transform(test_section['L'][ite])
        for name,model in new_models().items():
            model.fit(ztr,y_a[itr]);p=model.predict(zte);pr=model.predict_proba(zte)
            row={'scenario':scenario,'model':name,'n_train':len(itr),'n_test':len(ite),**metrics(truth,p,pr)};results.append(row)
            print('MODEL',row,flush=True)
            preds[(scenario,name)]=(truth,p,pr)
            if scenario.startswith('Held-out'):
                if name=='Logistic':final_lr=model;final_projector=projector;final_ztr=ztr;final_zte=zte
                pd.DataFrame({'barcode':test_section['obs'].index[ite],'true':truth,'pred':p,**{f'p_{l}':pr[:,i] for i,l in enumerate(LAYERS)}}).to_csv(DER/f'151674_predictions_{name}.csv',index=False)
        pd.DataFrame({'barcode':test_section['obs'].index[ite],'role':'test'}).to_csv(DER/f'split_{scenario.split()[0].lower()}_test.csv',index=False)
    pd.DataFrame(results).to_csv(OUT/'model_metrics.csv',index=False)
    pd.DataFrame({'gene':a['genes'][final_projector.selected_],'ensembl':a['ids'][final_projector.selected_]}).to_csv(OUT/'training_selected_genes.csv',index=False)
    # Bootstrap entire spatial tiles, not individual spots. Conditional on the single test section.
    groups=b['obs'].block.to_numpy();unique=np.unique(groups);boot=[]
    for name in new_models():
        y,p,pr=preds[('Held-out section 151674',name)]
        vals=[]
        for rep in range(500):
            drawn=rng.choice(unique,len(unique),replace=True)
            ii=np.concatenate([np.flatnonzero(groups==g) for g in drawn])
            vals.append(f1_score(y[ii],p[ii],labels=np.arange(7),average='macro',zero_division=0))
        lo,hi=np.quantile(vals,[.025,.975]);boot.append({'model':name,'metric':'macro_f1','lower':lo,'upper':hi,'n_blocks':len(unique),'reps':500})
    pd.DataFrame(boot).to_csv(OUT/'block_bootstrap_ci.csv',index=False)
    # Unsupervised comparisons on development section, K=7 fixed from the benchmark task.
    W=a['W'];degree=np.asarray(W.sum(1)).ravel();Wn=sparse.diags(1/np.maximum(degree,1))@W
    embedding=.5*final_ztr+.5*(Wn@final_ztr)
    no_neighbors=degree==0;embedding[no_neighbors]=final_ztr[no_neighbors]
    clusters={};clstats=[]
    for name,Z in [('Expression KMeans',final_ztr),('Neighbour-smoothed KMeans',embedding)]:
        cl=KMeans(n_clusters=7,n_init=20,random_state=SEED).fit_predict(Z)
        cm=confusion_matrix(y_a,cl,labels=np.arange(7));rr,cc=linear_sum_assignment(-cm)
        mapping=dict(zip(cc,rr));mapped=np.array([mapping[x] for x in cl])
        clusters[name]=mapped
        clstats.append({'method':name,'ari':adjusted_rand_score(y_a,cl),'nmi':normalized_mutual_info_score(y_a,cl)})
    pd.DataFrame(clstats).to_csv(OUT/'clustering_metrics.csv',index=False)
    tsne=TSNE(n_components=2,perplexity=30,init='pca',learning_rate='auto',max_iter=1000,random_state=SEED).fit_transform(final_ztr)
    dd=a['obs'].copy();dd['pc1']=final_ztr[:,0];dd['pc2']=final_ztr[:,1];dd['tsne1']=tsne[:,0];dd['tsne2']=tsne[:,1]
    for name,c in clusters.items():dd[name]=c
    dd.to_csv(DER/'151673_embedding_clusters.csv')
    # Back-project PC coefficients to clipped standardized selected-gene features.
    coefs=final_lr.coef_@final_projector.pca_.components_
    importance=np.max(np.abs(coefs),axis=0);order=np.argsort(importance)[::-1]
    pd.DataFrame({'gene':a['genes'][final_projector.selected_][order],'max_abs_coefficient':importance[order],**{l:coefs[i,order] for i,l in enumerate(LAYERS)}}).to_csv(OUT/'logistic_gene_coefficients.csv',index=False)
    y,p,pr=preds[('Held-out section 151674','Logistic')]
    report=classification_report(y,p,target_names=LAYERS,output_dict=True,zero_division=0)
    dump(report,OUT/'classification_report.json')
    # Error rates at reference boundaries versus interiors; a post-hoc diagnostic only.
    rows,cols=b['W'].nonzero();boundary=np.zeros(len(y),bool);boundary[rows[y[rows]!=y[cols]]]=True
    err={'boundary_n':int(boundary.sum()),'interior_n':int((~boundary).sum()),'boundary_error':float((y[boundary]!=p[boundary]).mean()),'interior_error':float((y[~boundary]!=p[~boundary]).mean())}
    d=b['obs'].copy();d['true']=y;d['pred']=p;d['confidence']=pr.max(1);d['entropy']=-(pr*np.log(np.maximum(pr,1e-12))).sum(1)/np.log(7);d['boundary']=boundary;d.to_csv(DER/'151674_diagnostics.csv')
    summary={'design':{'seed':SEED,'train':'151673','test':'151674','qc':{'min_umi':500,'min_genes':250,'max_pct_mt':30},'hvg':2000,'pcs':50,'moran_permutations':999},
             'sections':summ,'pca_variance_50':float(final_projector.pca_.explained_variance_ratio_.sum()),'pca_variance_ratio':final_projector.pca_.explained_variance_ratio_.tolist(),
             'eligible_genes_training':final_projector.allowed_count_,'ols_r2':ols.rsquared,'ols_n_blocks':o.block.nunique(),'error_diagnostics':err,
             'runtime_seconds':time.time()-t0,'versions':{m:importlib.metadata.version(m) for m in ['numpy','pandas','scipy','scikit-learn','h5py','matplotlib','statsmodels']},'python':sys.version,'platform':platform.platform()}
    dump(summary,OUT/'summary.json')
    print('COMPLETE',summary['runtime_seconds'],err,flush=True)
if __name__=='__main__':main()
