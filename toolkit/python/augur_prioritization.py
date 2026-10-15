#!/usr/bin/env python3
"""Audited native pertpy Augur ranking; never interprets cell CV as donor inference."""
import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
from figure_output import figure_format
from scrna_runtime import nested_get as get, resolved_rscript, sha256
from augur_contract import SKILL,validate_config

def path(value):return Path(os.path.expandvars(value)).expanduser().resolve()
def dump(file,value):file.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def read_inputs(config,output,config_path):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread
    if (get(config,'input.type') or 'seurat')=='seurat':
        r=resolved_rscript(SKILL,config)
        if not r:raise ValueError('Registered Seurat exporter runtime unavailable')
        here=Path(__file__).resolve()
        exporter=next((p for p in [here.parent/'export_cnmf_input.R',here.parents[1]/'R/export_cnmf_input.R'] if p.exists()),None)
        if exporter is None:raise ValueError('Bundled count exporter missing')
        dest=output/'_provenance/input';dest.mkdir(parents=True)
        subprocess.run([str(r),str(exporter),str(config_path),str(dest)],check=True)
        files={k:dest/(k+('.mtx' if k=='counts' else '.tsv')) for k in ['counts','features','barcodes','metadata']}
        orientation='genes_by_cells'
    else:
        files={k:path(get(config,'input.'+k)) for k in ['counts','features','barcodes','metadata']}
        orientation=get(config,'input.orientation') or 'genes_by_cells'
    def identifiers(file):
        table=pd.read_csv(file,sep='\t',header=None,dtype=str,keep_default_na=False)
        if table.shape[1]!=1:raise ValueError('Feature/barcode files must each have exactly one ID column')
        values=pd.Index(table.iloc[:,0])
        if values.has_duplicates or (values.str.strip()=='').any():raise ValueError('Duplicate or empty gene/cell identifiers')
        return values
    genes=identifiers(files['features']);cells=identifiers(files['barcodes'])
    counts=sparse.csr_matrix(mmread(files['counts']),dtype=float)
    if orientation=='genes_by_cells':counts=counts.T.tocsr()
    counts.sum_duplicates()
    if counts.shape!=(len(cells),len(genes)):raise ValueError('Matrix dimensions/orientation do not match identifiers')
    if not np.isfinite(counts.data).all() or (counts.data<0).any() or not np.allclose(counts.data,np.rint(counts.data),rtol=0,atol=1e-8):raise ValueError('Expected finite nonnegative integer raw UMI counts')
    if (np.asarray(counts.sum(axis=1)).ravel()<=0).any():raise ValueError('Zero-count input cells must be reviewed before analysis')
    metadata=pd.read_csv(files['metadata'],sep='\t',index_col=0,dtype=str,keep_default_na=False)
    if metadata.index.has_duplicates or set(metadata.index)!=set(cells):raise ValueError('Metadata must match every cell ID exactly')
    metadata=metadata.loc[cells]
    for role in ['sample','condition','cell_type','batch','donor']:
        column=get(config,'metadata.'+role)
        if column and (column not in metadata or metadata[column].str.strip().eq('').any()):raise ValueError('Missing or empty metadata: '+role)
    records={k:{'path':str(p),'sha256':sha256(p)} for k,p in files.items()}
    if get(config,'input.object'):records['object']={'path':str(path(get(config,'input.object'))),'sha256':sha256(path(get(config,'input.object')))}
    return counts,genes,cells,metadata,records

def audit_design(metadata,config):
    import numpy as np
    import pandas as pd
    sample=get(config,'metadata.sample');condition=get(config,'metadata.condition')
    columns=list(dict.fromkeys([sample,condition]+[get(config,'metadata.'+x) for x in ['batch','donor'] if get(config,'metadata.'+x)]))
    if metadata.groupby(sample)[columns[1:]].nunique().gt(1).any().any():raise ValueError('Each sample must map to one condition, batch and donor')
    samples=metadata[columns].drop_duplicates().set_index(sample)
    labels=[get(config,'contrast.denominator'),get(config,'contrast.numerator')]
    if not set(labels).issubset(set(samples[condition])):raise ValueError('Both contrast conditions must be present')
    samples=samples[samples[condition].isin(labels)]
    batch=get(config,'metadata.batch')
    design=pd.DataFrame({'intercept':1.,'condition':(samples[condition]==labels[1]).astype(float)},index=samples.index)
    if batch:
        design=pd.concat([design,pd.get_dummies(samples[batch],prefix='batch',drop_first=True,dtype=float)],axis=1)
        if np.linalg.matrix_rank(design.to_numpy())<design.shape[1]:raise ValueError('Condition is confounded with batch; Augur cannot isolate the biological response')
    return samples,design

def run(config,config_path):
    import numpy as np
    import pandas as pd
    import anndata
    import matplotlib.pyplot as plt
    output=path(config['output_dir']);technical=output/'_provenance';technical.mkdir(parents=True,exist_ok=True)
    if any(output.glob('auc_ranking.tsv')):raise ValueError('Existing Augur results must be preserved; choose a fresh output directory')
    counts,genes,cells,metadata,records=read_inputs(config,output,config_path)
    samples,design=audit_design(metadata,config)
    import pertpy as pt
    samples.to_csv(output/'sample_design.tsv',sep='\t');design.to_csv(output/'model_design.tsv',sep='\t')
    condition=get(config,'metadata.condition');cell_type=get(config,'metadata.cell_type');sample=get(config,'metadata.sample');replicate=get(config,'metadata.donor') or sample
    labels=[get(config,'contrast.denominator'),get(config,'contrast.numerator')]
    selected=metadata[condition].isin(labels)
    pd.DataFrame({'cell_id':cells,'included_in_contrast':selected.to_numpy()}).to_csv(output/'cell_inclusion.tsv',sep='\t',index=False)
    metadata=metadata[selected].copy();counts=counts[selected.to_numpy()]
    totals=np.asarray(counts.sum(axis=1)).ravel()
    qc=metadata[[sample,condition,cell_type]].copy();qc['umi']=totals;qc['detected_genes']=np.asarray((counts>0).sum(axis=1)).ravel()
    qc.groupby([sample,condition,cell_type],observed=True).agg(n_cells=('umi','size'),mean_umi=('umi','mean'),mean_detected_genes=('detected_genes','mean')).to_csv(output/'sample_cell_type_qc.tsv',sep='\t')
    detected=np.asarray((counts>0).sum(axis=0)).ravel();keep=detected>0
    pd.DataFrame({'gene':genes,'n_detected_cells':detected,'status':np.where(keep,'retained','zero_counts_in_contrast')}).to_csv(output/'feature_status.tsv',sep='\t',index=False)
    counts=counts[:,keep];genes=genes[keep]
    normalization=get(config,'preprocessing.normalization') or 'raw_counts'
    if normalization=='log1p_cp10k':
        from scipy import sparse
        counts=sparse.diags(10000/totals)@counts;counts=counts.tocsr();counts.data=np.log1p(counts.data)
    options={key:get(config,'augur.'+key) if get(config,'augur.'+key) is not None else default for key,default in {'n_subsamples':50,'subsample_size':20,'folds':3,'n_threads':1,'n_estimators':100,'seed':14,'feature_perc':.5,'var_quantile':.5,'span':.75,'min_replicates_per_condition':2}.items()}
    tasks=[];ranks=[];fold_tables=[];importance=[];eligible=[]
    for population in sorted(metadata[cell_type].unique()):
        subset=metadata[metadata[cell_type]==population];reason=[];audit={'cell_type':population}
        for index,label in enumerate(labels):
            obs=subset[subset[condition]==label]
            audit[f'condition_{index}_label']=label;audit[f'condition_{index}_cells']=len(obs);audit[f'condition_{index}_replicates']=obs[replicate].nunique()
            if len(obs)<options['subsample_size']:reason.append('insufficient_cells:'+label)
            if obs[replicate].nunique()<options['min_replicates_per_condition']:reason.append('insufficient_replicates:'+label)
        if not reason and get(config,'metadata.batch'):
            try:audit_design(subset,config)
            except ValueError as exc:reason.append(str(exc))
        audit['eligible']=not reason;audit['reason']=';'.join(reason);eligible.append(audit)
        if reason:
            tasks.append({'cell_type':population,'status':'untestable','reason':audit['reason']});continue
        try:
            native_obs=pd.DataFrame({'label':pd.Categorical(subset[condition],categories=labels),'cell_type':pd.Categorical(['population']*len(subset))},index=subset.index)
            data=anndata.AnnData(counts[metadata[cell_type].eq(population).to_numpy()].copy(),obs=native_obs,var=pd.DataFrame(index=genes))
            model=pt.tl.Augur('random_forest_classifier',n_estimators=options['n_estimators'],random_state=options['seed'])
            loaded=model.load(data,label_col='label',cell_type_col='cell_type')
            _,result=model.predict(loaded,n_subsamples=options['n_subsamples'],subsample_size=options['subsample_size'],folds=options['folds'],min_cells=2*options['subsample_size'],feature_perc=options['feature_perc'],var_quantile=options['var_quantile'],span=options['span'],n_threads=options['n_threads'],random_state=options['seed'])
            folds=result['full_results'].copy();folds['cell_type']=population
            if len(folds)!=options['n_subsamples']*options['folds'] or not np.isfinite(folds['augur_score']).all():raise ValueError('Incomplete or non-finite native CV results')
            score=float(result['summary_metrics'].loc['mean_augur_score','population'])
            if not 0<=score<=1 or not np.isclose(score,folds.augur_score.mean()):raise ValueError('Native summary disagrees with fold AUCs')
            fold_tables.append(folds);imp=result['feature_importances'].copy();imp['cell_type']=population;importance.append(imp)
            ranks.append({'cell_type':population,'auc':score,'cv_unit':'cell','interpretation':'exploratory_transcriptional_separability','n_cells':len(subset)})
            tasks.append({'cell_type':population,'status':'completed','reason':''})
        except Exception as exc:tasks.append({'cell_type':population,'status':'failed','reason':str(exc)})
    pd.DataFrame(eligible).to_csv(output/'cell_type_eligibility.tsv',sep='\t',index=False)
    pd.DataFrame(tasks).to_csv(output/'task_status.tsv',sep='\t',index=False)
    ranked=pd.DataFrame(ranks,columns=['cell_type','auc','cv_unit','interpretation','n_cells']).sort_values(['auc','cell_type'],ascending=[False,True]);ranked.insert(0,'rank',range(1,len(ranked)+1));ranked.to_csv(output/'auc_ranking.tsv',sep='\t',index=False)
    if fold_tables:
        folds=pd.concat(fold_tables,ignore_index=True);folds.to_csv(output/'fold_auc.tsv',sep='\t',index=False)
        means=folds.groupby(['cell_type','idx'],sort=False).augur_score.mean().reset_index(name='auc');means.to_csv(output/'subsample_auc.tsv',sep='\t',index=False)
        pd.concat(importance,ignore_index=True).to_csv(output/'feature_importance.tsv',sep='\t',index=False)
        from matplotlib.backends.backend_pdf import PdfPages
        for family in ['auc_ranking','subsample_auc']:
            pdf=PdfPages(output/(family+'.pdf')) if figure_format(config)=='pdf' else None
            try:
                for page,offset in enumerate(range(0,len(ranked),24),1):
                    part=ranked.iloc[offset:offset+24]
                    fig,ax=plt.subplots(figsize=(8,max(3,.4*len(part)+1)))
                    if family=='auc_ranking':
                        ax.barh(part.cell_type,part.auc,color='#0072B2');ax.invert_yaxis()
                        ax.set(xlabel='Cross-validated ROC AUC (cell-level CV)',title='Cell-type response prioritization')
                    else:
                        ax.boxplot([means.loc[means.cell_type==p,'auc'] for p in part.cell_type],tick_labels=part.cell_type,vert=False)
                        ax.set(xlabel='Subsample AUC; sampling variation is not a donor confidence interval')
                    ax.axvline(.5,color='grey',ls='--');ax.set_xlim(0,1);fig.tight_layout()
                    if pdf:pdf.savefig(fig,dpi=300)
                    else:fig.savefig(output/(family+f'_page{page:02d}.png'),dpi=300)
                    plt.close(fig)
            finally:
                if pdf:pdf.close()
    failed=any(row['status']=='failed' for row in tasks) or not ranks
    warnings=['Native Augur uses pooled cell-level CV, not donor-held-out validation. AUC is not an absolute perturbation effect size or a donor-level significance test.','Feature importance is classifier importance, not causal evidence. No batch correction or covariate regression was performed.']
    if not get(config,'metadata.batch'):warnings.append('No batch field declared: technical confounding cannot be excluded')
    dump(output/'interpretation.json',{'status':'failed' if failed else 'completed','contrast':config['contrast'],'normalization':normalization,'cv_unit':'cell','baseline_auc':.5,'warnings':warnings})
    pd.DataFrame([{'family':'Augur','status':'exploratory','reason':message} for message in warnings]).to_csv(output/'interpretation_audit.tsv',sep='\t',index=False)
    source=Path(inspect.getfile(pt.tl.Augur))
    dump(technical/'run_manifest.json',{'skill':SKILL,'status':'failed' if failed else 'completed','inputs':records,'config_sha256':sha256(config_path),'config':config,'options':options,'implementation':'pertpy.Augur','implementation_source_sha256':sha256(source),'versions':{name:importlib.metadata.version(name) for name in ['pertpy','anndata','numpy','pandas','scikit-learn','scipy','scanpy','scikit-misc','matplotlib','joblib']},'warnings':warnings})
    dump(technical/'workflow_state.json',{'status':'failed' if failed else 'complete','completed_cell_types':len(ranks),'untestable_cell_types':sum(t['status']=='untestable' for t in tasks)})
    return 1 if failed else 0

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('config',type=Path);args=parser.parse_args();config=json.loads(args.config.read_text())
    errors=validate_config(config,get)
    if errors:raise ValueError('; '.join(errors))
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    os.environ['MPLCONFIGDIR']=str(path(config['output_dir'])/'_provenance/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    return run(config,args.config.resolve())
if __name__=='__main__':raise SystemExit(main())
