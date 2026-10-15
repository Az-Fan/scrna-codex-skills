"""Dependency-free binary Augur prioritization contract."""
import math
import os
from pathlib import Path
from figure_output import figure_format
SKILL='18-scrna-prioritize-perturbation'

def validate_config(config,get):
    errors=[]
    try:figure_format(config)
    except ValueError as exc:errors.append(str(exc))
    kind=get(config,'input.type') or 'seurat'
    if kind not in ('seurat','matrix'):errors.append('input.type must be seurat or matrix')
    fields=['input.object'] if kind=='seurat' else ['input.counts','input.features','input.barcodes','input.metadata']
    for field in fields:
        value=get(config,field)
        if not isinstance(value,str) or not Path(os.path.expandvars(value)).expanduser().is_file():errors.append(field+' must name an existing file')
    destination=get(config,'output_dir')
    if isinstance(destination,str) and destination.strip():
        output=Path(os.path.expandvars(destination)).expanduser().resolve()
        for field in fields:
            value=get(config,field)
            if isinstance(value,str) and (Path(os.path.expandvars(value)).expanduser().resolve()==output or output in Path(os.path.expandvars(value)).expanduser().resolve().parents):
                errors.append('Input files must be outside output_dir to preserve source data: '+field)
    if get(config,'input.counts_source')!='raw_umi':errors.append('input.counts_source must be raw_umi')
    assay=get(config,'input.assay') or 'RNA'
    if not isinstance(assay,str) or not assay.strip():errors.append('input.assay must be a nonempty name')
    elif assay.lower() in {'sct','integrated'}:errors.append('Use uncorrected RNA counts, not SCT or integrated assays')
    if (get(config,'input.orientation') or 'genes_by_cells') not in ('genes_by_cells','cells_by_genes'):errors.append('Invalid matrix orientation')
    for role in ['sample','condition','cell_type']:
        value=get(config,'metadata.'+role)
        if not isinstance(value,str) or not value.strip():errors.append('metadata.'+role+' is required')
    for role in ['batch','donor']:
        value=get(config,'metadata.'+role)
        if value is not None and (not isinstance(value,str) or not value.strip()):errors.append('metadata.'+role+' must be a nonempty column name')
    labels=[get(config,'contrast.'+key) for key in ['numerator','denominator']]
    if any(not isinstance(x,str) or not x.strip() for x in labels) or labels[0]==labels[1]:errors.append('contrast requires two distinct nonempty numerator/denominator labels')
    defaults={'n_subsamples':50,'subsample_size':20,'folds':3,'n_threads':1,'n_estimators':100,'seed':14,'min_replicates_per_condition':2}
    minimum={'n_subsamples':2,'subsample_size':3,'folds':2,'n_threads':1,'n_estimators':10,'seed':0,'min_replicates_per_condition':2}
    for field,default in defaults.items():
        value=get(config,'augur.'+field);value=default if value is None else value
        if type(value) is not int or value<minimum[field]:errors.append('augur.'+field+' must be an integer >= '+str(minimum[field]))
    folds=get(config,'augur.folds') or 3;size=get(config,'augur.subsample_size') or 20
    if type(folds) is int and type(size) is int and folds>size:errors.append('folds cannot exceed subsample_size per condition')
    for field,default in [('feature_perc',.5),('var_quantile',.5),('span',.75)]:
        value=get(config,'augur.'+field);value=default if value is None else value
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<value<=(1 if field!='var_quantile' else .999999):errors.append('augur.'+field+' must be finite and in its valid (0,1] range')
    if (get(config,'preprocessing.normalization') or 'raw_counts') not in ('raw_counts','log1p_cp10k'):errors.append('normalization must be raw_counts or log1p_cp10k')
    if get(config,'augur.cv_unit') not in (None,'cell'):errors.append('Native Augur implements cell-level CV; donor-blocked CV is not implemented')
    return errors
