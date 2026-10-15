"""Small balanced multi-sample counts with responsive, null and ineligible types."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

def create(destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(42);means=rng.gamma(2,2,300)+.15
    matrices=[];rows=[]
    for population in ['EC','EC_2','rare']:
        for sample,condition,batch in [('c1','control','b1'),('c2','control','b2'),('d1','case','b1'),('d2','case','b2')]:
            n=40 if population!='rare' else 2
            mu=means.copy()
            if population=='EC' and condition=='case':mu[:80]*=8
            matrices.append(rng.poisson(mu,size=(n,len(mu))))
            offset=len(rows)
            rows.extend({'cell_id':f'cell{offset+i:04d}','sample':sample,'condition':condition,'batch':batch,'cell_type':population} for i in range(n))
    metadata=pd.DataFrame(rows);metadata.loc[0,'cell_id']='001';metadata.loc[1,'cell_id']='NA'
    # IDs are deliberately ordinary strings, including values pandas can coerce.
    mmwrite(destination/'counts.mtx',sparse.csr_matrix(np.vstack(matrices)))
    (destination/'features.tsv').write_text('\n'.join(f'g{i:03d}' for i in range(300))+'\n')
    (destination/'barcodes.tsv').write_text('\n'.join(metadata.cell_id)+'\n')
    metadata.to_csv(destination/'metadata.tsv',sep='\t',index=False)
    return destination
if __name__=='__main__':create(sys.argv[1])
