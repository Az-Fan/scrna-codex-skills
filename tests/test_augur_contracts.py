"""Augur boundary checks without importing the scientific engine at startup."""
import copy,json
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'toolkit/python'))
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy import sparse
from augur_contract import SKILL,validate_config
from augur_prioritization import read_inputs,audit_design
from scrna_runtime import nested_get,default_argv

class AugurContracts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        mmwrite(self.root/'counts.mtx',sparse.csr_matrix([[1,2,3],[2,1,1],[3,1,4],[1,2,1]]))
        (self.root/'features.tsv').write_text('001\nNA\ng3\n');(self.root/'barcodes.tsv').write_text('001\nNA\nc3\nc4\n')
        (self.root/'metadata.tsv').write_text('cell_id\tsample\tcondition\tbatch\tcell_type\nc4\td2\tcase\tb2\tEC_2\n001\tc1\tcontrol\tb1\tEC\nc3\td1\tcase\tb1\tEC_2\nNA\tc2\tcontrol\tb2\tEC\n')
        self.cfg={'project':{'id':'test'},'input':{'type':'matrix','orientation':'cells_by_genes','counts_source':'raw_umi',**{x:str(self.root/(x+('.mtx' if x=='counts' else '.tsv'))) for x in ['counts','features','barcodes','metadata']}},'metadata':{'sample':'sample','condition':'condition','cell_type':'cell_type','batch':'batch'},'contrast':{'numerator':'case','denominator':'control'},'output_dir':str(self.root/'result')}
    def test_alignments_and_identifiers(self):
        counts,genes,cells,metadata,records=read_inputs(self.cfg,self.root,self.root/'cfg.json')
        self.assertEqual(list(cells),['001','NA','c3','c4']);self.assertEqual(list(genes),['001','NA','g3'])
        self.assertEqual(list(metadata.index),list(cells));self.assertEqual(metadata.loc['001','sample'],'c1')
        self.assertEqual(counts.shape,(4,3));self.assertEqual(len(records),4)
    def test_count_and_id_rejections(self):
        for value in [-1,.5,float('nan')]:
            mmwrite(self.root/'counts.mtx',np.array([[value,2,3],[2,1,1],[3,1,4],[1,2,1]]))
            with self.assertRaisesRegex(ValueError,'integer raw UMI'):read_inputs(self.cfg,self.root,self.root/'cfg.json')
        (self.root/'features.tsv').write_text('g1\ng1\ng3\n')
        with self.assertRaisesRegex(ValueError,'Duplicate'):read_inputs(self.cfg,self.root,self.root/'cfg.json')
    def test_biological_design_and_batch_confounding(self):
        _,_,_,metadata,_=read_inputs(self.cfg,self.root,self.root/'cfg.json')
        samples,design=audit_design(metadata,self.cfg);self.assertEqual(len(samples),4);self.assertEqual(design.shape,(4,3))
        metadata['batch']=metadata['condition']
        with self.assertRaisesRegex(ValueError,'confounded'):audit_design(metadata,self.cfg)
    def test_configuration_scope_and_explicit_runtime(self):
        self.assertEqual(validate_config(self.cfg,nested_get),[])
        for field,value in [('folds',30),('n_subsamples',True),('subsample_size',0),('cv_unit','donor')]:
            cfg=copy.deepcopy(self.cfg);cfg['augur']={field:value};self.assertTrue(validate_config(cfg,nested_get))
        cfg=copy.deepcopy(self.cfg);cfg['contrast']['numerator']='control';self.assertTrue(validate_config(cfg,nested_get))
        for section,field in [('input','type'),('input','orientation'),('preprocessing','normalization'),('augur','cv_unit')]:
            cfg=copy.deepcopy(self.cfg);cfg.setdefault(section,{})[field]=['invalid'];self.assertTrue(validate_config(cfg,nested_get))
        cfg=copy.deepcopy(self.cfg);cfg['output_dir']=str(self.root);self.assertTrue(validate_config(cfg,nested_get))
        cfg=copy.deepcopy(self.cfg);cfg['runtime']={'augur_python':'/missing/explicit/python'}
        self.assertIsNone(default_argv(SKILL,self.root/'cfg.json',cfg))
