"""Fixed, descriptive extensions; original inputs and benchmark outputs are read-only.

Two adjacent sections from ONE donor. No patient-level inference or test-set tuning.
Run with OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2.
"""
from pathlib import Path
import json, time, hashlib
import h5py
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse, stats
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, f1_score, confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from analyze import ROOT, SEED, LAYERS, MARKERS, lognorm, graph, GeneProjector, metrics, dump

E=ROOT/'results/extended'; D=ROOT/'data/derived/extended'
SAMPLES=['151673','151674']
IMAGE_FEATURES=['red_mean','green_mean','blue_mean','red_sd','green_sd','blue_sd','dark_fraction','gradient_energy']

def load_section(sample):
    """Read existing QC decisions and original counts without rewriting either."""
    obs=pd.read_csv(ROOT/f'data/derived/{sample}_all_spots.csv',index_col=0)
    with h5py.File(ROOT/f'data/raw/{sample}/counts.h5','r') as h:
        m=h['matrix']; X=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr().astype(np.float32)
        genes=np.array([x.decode() for x in m['features/name'][:]])
        ids=np.array([x.decode() for x in m['features/id'][:]])
        assert list(obs.index)==[x.decode() for x in m['barcodes'][:]]
    keep=obs.eligible.to_numpy();o=obs.loc[keep].copy();L=lognorm(X[keep]);W,radius=graph(o)
    return dict(obs=o,L=L,genes=genes,ids=ids,W=W,radius=radius,y=np.array([LAYERS.index(x) for x in o.layer]))

def image_features(sample,obs):
    """Simple colour/texture descriptors, not nuclei counts or stain deconvolution."""
    img=np.asarray(Image.open(ROOT/f'data/raw/{sample}/tissue_hires_image.png').convert('RGB'))/255.
    scale=json.loads((ROOT/f'data/raw/{sample}/scalefactors_json.json').read_text())
    # Same physical capture footprint in each sample; square window, no fake resolution.
    radius=max(2,int(round(scale['spot_diameter_fullres']*scale['tissue_hires_scalef']/2)))
    rows=[]
    for x,y in obs[['x','y']].to_numpy():
        x,y=int(round(x)),int(round(y));p=img[max(0,y-radius):min(img.shape[0],y+radius+1),max(0,x-radius):min(img.shape[1],x+radius+1)]
        v=p.reshape(-1,3);gray=p@np.array([.299,.587,.114]);dx=np.diff(gray,axis=1);dy=np.diff(gray,axis=0)
        rows.append([*v.mean(0),*v.std(0),float((gray<.65).mean()),float((np.mean(dx**2)+np.mean(dy**2))/2)])
    return pd.DataFrame(rows,index=obs.index,columns=IMAGE_FEATURES),radius

def bootstrap_paired(y,p0,p1,groups,rng,reps=500):
    unique=np.unique(groups);index=[np.flatnonzero(groups==g) for g in unique];values=[]
    for _ in range(reps):
        ix=np.concatenate([index[i] for i in rng.integers(0,len(index),len(index))])
        values.append(f1_score(y[ix],p1[ix],labels=np.arange(7),average='macro',zero_division=0)-f1_score(y[ix],p0[ix],labels=np.arange(7),average='macro',zero_division=0))
    return np.asarray(values)

def main():
    start=time.time();E.mkdir(exist_ok=True);D.mkdir(exist_ok=True);rng=np.random.default_rng(SEED+100)
    protected=[ROOT/'results/model_metrics.csv',ROOT/'results/marker_statistics.csv',ROOT/'results/summary.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    data={s:load_section(s) for s in SAMPLES};a,b=data.values();assert np.array_equal(a['ids'],b['ids'])
    selected=pd.read_csv(ROOT/'results/training_selected_genes.csv')
    imap={x:i for i,x in enumerate(a['ids'])};ix=np.array([imap[x] for x in selected.ensembl]);names=a['genes'][ix]
    means={};effects={};detection={};Y={};summary={}
    for s,d in data.items():
        Y[s]=d['L'][:,ix].toarray().astype(np.float64);y=d['y']
        means[s]=np.vstack([Y[s][y==k].mean(0) for k in range(7)])
        effects[s]=np.vstack([Y[s][y==k].mean(0)-Y[s][y!=k].mean(0) for k in range(7)])
        detection[s]=np.vstack([(Y[s][y==k]>0).mean(0) for k in range(7)])
    # Correlate layer-versus-rest contrasts rather than baseline transcript abundance.
    cross=np.array([[stats.spearmanr(effects[SAMPLES[0]][i],effects[SAMPLES[1]][j]).statistic for j in range(7)] for i in range(7)])
    pd.DataFrame(cross,index=LAYERS,columns=LAYERS).to_csv(E/'cross_section_layer_correlation.csv')
    rep=[]
    for k,l in enumerate(LAYERS):
        rep.append({'layer':l,'rho':cross[k,k],'best_match':LAYERS[np.argmax(cross[k])],
                    'sign_agreement':np.mean(np.sign(effects[SAMPLES[0]][k])==np.sign(effects[SAMPLES[1]][k]))})
    pd.DataFrame(rep).to_csv(E/'layer_replication.csv',index=False)
    # Discovery restricted to the 2,000 previously selected TRAINING genes.
    discovery=[];top=[];used=set()
    for k,l in enumerate(LAYERS):
        x=Y[SAMPLES[0]];hit=a['y']==k
        U,p=stats.mannwhitneyu(x[hit],x[~hit],axis=0,alternative='two-sided',method='asymptotic')
        q=multipletests(p,method='fdr_bh')[1];auc=U/(hit.sum()*(~hit).sum())
        order=np.argsort(effects[SAMPLES[0]][k])[::-1]
        for j in range(len(ix)):
            discovery.append({'layer':l,'gene':names[j],'ensembl':a['ids'][ix[j]],'training_delta':effects[SAMPLES[0]][k,j],
               'validation_delta':effects[SAMPLES[1]][k,j],'auc':auc[j],'p_nominal':p[j],'q_within_layer':q[j],
               'training_detected':detection[SAMPLES[0]][k,j],'validation_detected':detection[SAMPLES[1]][k,j]})
        chosen=[]
        for j in order:
            if detection[SAMPLES[0]][k,j]<.1 or effects[SAMPLES[0]][k,j]<=0 or names[j] in used:continue
            chosen.append({'layer':l,'gene':names[j],'column':int(j),'training_delta':effects[SAMPLES[0]][k,j],
              'validation_delta':effects[SAMPLES[1]][k,j],'training_auc':auc[j],'q_nominal':q[j]});used.add(names[j])
            if len(chosen)==3:break
        top+=chosen
    pd.DataFrame(discovery).to_csv(E/'training_gene_screen.csv',index=False)
    pd.DataFrame(top).to_csv(E/'discovered_marker_panel.csv',index=False)
    np.savez_compressed(D/'layer_profiles.npz',genes=names,indices=ix,train_means=means[SAMPLES[0]],test_means=means[SAMPLES[1]],
       train_effects=effects[SAMPLES[0]],test_effects=effects[SAMPLES[1]],train_detection=detection[SAMPLES[0]],test_detection=detection[SAMPLES[1]])
    summary['replication']={'matched_layers':int(np.sum(np.argmax(cross,axis=1)==np.arange(7))),
      'rho_min':float(np.diag(cross).min()),'rho_max':float(np.diag(cross).max()),'positive_discovery_replications':sum(t['validation_delta']>0 for t in top),'panel_size':len(top)}
    print('Layer replication and 2,000-gene screen finished.',flush=True)
    # Reference-guided geometric coordinate; NOT pseudotime and NOT a model feature.
    gradients=[];depthstats=[];imageinfo={};neighborhood=[];scalestats=[]
    for s,d in data.items():
        o=d['obs'];xy=o[['x','y']].to_numpy();y=d['y'];centers=np.vstack([xy[y==k].mean(0) for k in range(7)])
        unit=centers[6]-centers[0];unit/=np.linalg.norm(unit);t=xy@unit;lo,hi=np.quantile(t,[.02,.98]);depth=np.clip((t-lo)/(hi-lo),0,1)
        dd=o.copy();dd['geometric_depth']=depth;dd.to_csv(D/f'{s}_spatial_geometry.csv')
        for k,l in enumerate(LAYERS):depthstats.append({'sample':s,'layer':l,'median_depth':np.median(depth[y==k]),'q25':np.quantile(depth[y==k],.25),'q75':np.quantile(depth[y==k],.75)})
        mv=pd.read_csv(ROOT/f'data/derived/{s}_markers.csv',index_col=0).loc[o.index]
        bins=np.minimum((depth*20).astype(int),19)
        for g in MARKERS:
            v=mv[g].to_numpy()
            for q in range(20):
                hit=bins==q
                gradients.append({'sample':s,'gene':g,'bin':q,'depth':(q+.5)/20,'n':int(hit.sum()),'mean':float(v[hit].mean()) if hit.any() else np.nan})
        feat,radius=image_features(s,o);feat.to_csv(D/f'{s}_image_features.csv');d['image_features']=feat
        imageinfo[s]={'patch_radius_pixels':radius,'axis_vector':unit.tolist(),'axis_p02':float(lo),'axis_p98':float(hi)}
        # Label-shuffle null preserves global class counts, not spatial dependence.
        rows,cols=d['W'].nonzero();edge=rows<cols;rr,cc=rows[edge],cols[edge]
        def counts(v):
            c=np.bincount(v[rr]*7+v[cc],minlength=49).reshape(7,7)
            return c+c.T
        observed=counts(y);null=np.array([counts(rng.permutation(y)) for _ in range(499)])
        expected=null.mean(0);logratio=np.log2((observed+.5)/(expected+.5))
        np.savez_compressed(D/f'{s}_neighborhood.npz',observed=observed,expected=expected,log2_enrichment=logratio,
            z=(observed-expected)/np.maximum(null.std(0),1e-9))
        neighborhood.append({'sample':s,'edges':len(rr),'same_layer_fraction':np.mean(y[rr]==y[cc]),'permutations':499,
          'adjacent_layer_fraction_among_cross_layer':np.mean(np.abs(y[rr][y[rr]!=y[cc]]-y[cc][y[rr]!=y[cc]])==1)})
        # Moran's I over increasing graph radii, no scale selection or extra p-values.
        tree=cKDTree(xy);base=np.median(tree.query(xy,k=2)[0][:,1])
        for scale in [1.15,2.1,3.1,4.1,6.1]:
            pairs=tree.query_pairs(base*scale,output_type='ndarray');r=np.r_[pairs[:,0],pairs[:,1]];c=np.r_[pairs[:,1],pairs[:,0]]
            for g in MARKERS:
                v=mv[g].to_numpy();z=v-v.mean();mi=len(v)/len(r)*np.sum(z[r]*z[c])/np.sum(z*z)
                scalestats.append({'sample':s,'gene':g,'radius_in_spacings':scale,'edges':len(pairs),'moran_i':mi})
    pd.DataFrame(gradients).to_csv(E/'depth_expression_profiles.csv',index=False)
    pd.DataFrame(depthstats).to_csv(E/'depth_by_layer.csv',index=False)
    pd.DataFrame(neighborhood).to_csv(E/'neighborhood_summary.csv',index=False)
    pd.DataFrame(scalestats).to_csv(E/'multiscale_moran.csv',index=False)
    summary['image_and_geometry']=imageinfo
    print('Geometry, H&E descriptors, neighborhood permutations and scales finished.',flush=True)
    # Reproduce the original expression representation; fixed sensitivity analysis on training only.
    projector=GeneProjector(a['genes']).fit(a['L']);za=projector.transform(a['L']);zb=projector.transform(b['L'])
    assert np.array_equal(projector.selected_,ix)
    W=a['W'];degree=np.asarray(W.sum(1)).ravel();wn=sparse.diags(1/np.maximum(degree,1))@W
    neighbor=wn@za;neighbor[degree==0]=za[degree==0];stability=[];pairwise=[];support=None
    for alpha in [0,.25,.5,.75,1.]:
        z=(1-alpha)*za+alpha*neighbor;labels=[]
        for seed in range(8):
            cl=KMeans(n_clusters=7,n_init=5,random_state=SEED+seed).fit_predict(z);labels.append(cl)
            stability.append({'alpha':alpha,'seed':SEED+seed,'ari':adjusted_rand_score(a['y'],cl),'nmi':normalized_mutual_info_score(a['y'],cl)})
        for i in range(8):
            for j in range(i):pairwise.append({'alpha':alpha,'seed_i':i,'seed_j':j,'pairwise_ari':adjusted_rand_score(labels[i],labels[j])})
        if alpha==.5:
            aligned=[]
            for cl in labels:
                cm=confusion_matrix(labels[0],cl,labels=np.arange(7));rr,cc=linear_sum_assignment(-cm);mapping=dict(zip(cc,rr));aligned.append([mapping[x] for x in cl])
            aligned=np.array(aligned);support=np.max(np.stack([(aligned==k).mean(0) for k in range(7)]),axis=0)
    pd.DataFrame(stability).to_csv(E/'clustering_sensitivity.csv',index=False)
    pd.DataFrame(pairwise).to_csv(E/'clustering_seed_agreement.csv',index=False)
    sd=a['obs'].copy();sd['seed_consensus_fraction']=support;sd.to_csv(D/'151673_clustering_stability.csv')
    print('Forty fixed clustering sensitivity fits finished.',flush=True)
    # Existing probabilities are the reference; do not retrain the benchmark RF.
    models={};y=b['y'];groups=b['obs'].block.to_numpy()
    for model in ['Logistic','RandomForest']:
        pr=pd.read_csv(ROOT/f'data/derived/151674_predictions_{model}.csv').set_index('barcode').loc[b['obs'].index]
        prob=pr[[f'p_{l}' for l in LAYERS]].to_numpy(dtype=np.float64)
        prob/=prob.sum(axis=1,keepdims=True)  # remove CSV/float32 round-off only
        models[model]=(pr.pred.to_numpy(),prob)
    # Image-only and a fixed concatenation; no coordinate, label or test-fitted transform.
    scaler=StandardScaler().fit(a['image_features']);ha=scaler.transform(a['image_features']);hb=scaler.transform(b['image_features'])
    for name,xtr,xte in [('Image only',ha,hb),('Expression + image',np.c_[za,ha],np.c_[zb,hb])]:
        model=LogisticRegression(C=1,solver='lbfgs',max_iter=2000,class_weight='balanced',random_state=SEED).fit(xtr,a['y'])
        prob=model.predict_proba(xte);pred=prob.argmax(1);models[name]=(pred,prob)
        pd.DataFrame({'barcode':b['obs'].index,'true':y,'pred':pred,**{f'p_{l}':prob[:,k] for k,l in enumerate(LAYERS)}}).to_csv(D/f'{name.replace(" ","_").replace("+","plus")}_predictions.csv',index=False)
    pd.DataFrame([{'model':name,**metrics(y,p,pr)} for name,(p,pr) in models.items()]).to_csv(E/'modality_model_metrics.csv',index=False)
    # Conditional paired spatial-block bootstrap; the section remains the single replicate.
    paired=[]
    for name in ['RandomForest','Expression + image']:
        vals=bootstrap_paired(y,models['Logistic'][0],models[name][0],groups,rng)
        delta=f1_score(y,models[name][0],average='macro')-f1_score(y,models['Logistic'][0],average='macro')
        paired.append({'comparison':name+' minus Logistic','delta_macro_f1':delta,'lower':np.quantile(vals,.025),'upper':np.quantile(vals,.975),'reps':len(vals)})
    pd.DataFrame(paired).to_csv(E/'paired_model_difference.csv',index=False)
    severity=[];risk=[];risk_summary=[]
    for name,(p,pr) in models.items():
        cortex=y<6
        severity.append({'model':name,'gray_white_accuracy':np.mean((y==6)==(p==6)),
            'cortex_exact_accuracy':np.mean(p[cortex]==y[cortex]),'cortex_within_one_layer':np.mean((p[cortex]<6)&(np.abs(p[cortex]-y[cortex])<=1)),
            'cortex_to_white_error':np.mean(p[cortex]==6)})
        order=np.argsort(-pr.max(1),kind='stable');err=(p[order]!=y[order]).astype(float);cum=np.cumsum(err)/np.arange(1,len(y)+1)
        for coverage in np.linspace(.1,1,19):
            n=int(np.ceil(coverage*len(y)));risk.append({'model':name,'coverage':n/len(y),'n':n,'risk':cum[n-1],
                'min_confidence':pr.max(1)[order[n-1]],'retained_classes':len(np.unique(y[order[:n]]))})
        risk_summary.append({'model':name,'aurc_discrete':cum.mean(),'full_risk':err.mean()})
    pd.DataFrame(severity).to_csv(E/'hierarchical_error.csv',index=False)
    pd.DataFrame(risk).to_csv(E/'selective_risk.csv',index=False)
    pd.DataFrame(risk_summary).to_csv(E/'selective_risk_summary.csv',index=False)
    # Boundary rings, with uncertainty resampled at the SAME spatial block level.
    o=pd.read_csv(ROOT/'data/derived/151674_diagnostics.csv',index_col=0).loc[b['obs'].index]
    distance=dijkstra(b['W'],directed=False,indices=np.flatnonzero(o.boundary),min_only=True,unweighted=True)
    ring=np.minimum(distance,3).astype(int);o['boundary_ring']=ring
    p0=models['Logistic'][0];p1=models['RandomForest'][0];o['model_agreement_group']=np.where((p0==y)&(p1==y),0,np.where(p0==y,1,np.where(p1==y,2,3)))
    o.to_csv(D/'151674_extended_diagnostics.csv')
    unique=np.unique(groups);counts=np.zeros((len(unique),4));errors=counts.copy()
    for i,g in enumerate(unique):
        for k in range(4):
            hit=(groups==g)&(ring==k);counts[i,k]=hit.sum();errors[i,k]=np.sum(p0[hit]!=y[hit])
    weights=rng.multinomial(len(unique),np.ones(len(unique))/len(unique),size=500)
    denom=weights@counts;boot=np.divide(weights@errors,denom,out=np.full_like(denom,np.nan),where=denom>0)
    pd.DataFrame([{'ring':k,'n':int(counts[:,k].sum()),'error':errors[:,k].sum()/max(counts[:,k].sum(),1),
       'lower':np.nanquantile(boot[:,k],.025),'upper':np.nanquantile(boot[:,k],.975)} for k in range(4)]).to_csv(E/'boundary_ring_errors.csv',index=False)
    summary['agreement_counts']={str(k):int(np.sum(o.model_agreement_group==k)) for k in range(4)}
    summary['stability']={'median_support':float(np.median(support)),'fraction_below_075':float(np.mean(support<.75))}
    summary['runtime_seconds']=time.time()-start;summary['design']={'seed':SEED+100,'neighborhood_permutations':499,'paired_bootstrap_reps':500,
      'kmeans_seeds':8,'kmeans_n_init':5,'alphas':[0,.25,.5,.75,1],'image_features':IMAGE_FEATURES,'independent_donors':1}
    for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    dump(summary,E/'summary.json');print('EXTENSION COMPLETE',summary,flush=True)

if __name__=='__main__':main()
