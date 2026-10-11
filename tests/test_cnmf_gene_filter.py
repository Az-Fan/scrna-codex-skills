"""Verify optional raw-count filtering and descriptive detection annotations."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import warnings
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'toolkit/python'))
import numpy as np
import pandas as pd
from scipy import sparse
from gene_filter import filter_counts,write_top_gene_tables
from gene_filter_contract import validate_thresholds
from cnmf_discovery import preparation_parameters

class GeneFilter(unittest.TestCase):
    def setUp(self):
        self.counts=sparse.csr_matrix([[2,1,4,0,0],[1,0,2,0,0],[3,0,0,0,1],[1,0,0,0,1]])
        self.genes=['common','rare','half','zero','half2']

    def test_default_preserves_all_nonzero_genes_and_all_cells(self):
        filtered,keep,audit=filter_counts(self.counts,self.genes)
        np.testing.assert_array_equal(keep,[True,True,True,False,True])
        np.testing.assert_array_equal(filtered.toarray(),self.counts.toarray()[:,keep])
        self.assertEqual(list(audit.n_detected_cells),[4,1,2,0,2])
        self.assertEqual(list(audit.status),['retained','retained','retained','zero_counts','retained'])

    def test_joint_threshold_and_umi_accounting(self):
        filtered,keep,audit=filter_counts(self.counts,self.genes,2,.75)
        np.testing.assert_array_equal(keep,[True,False,False,False,False])
        self.assertEqual(filtered.shape,(4,1))
        self.assertEqual(filtered.sum(),7)
        self.assertEqual(self.counts.sum(),16)
        self.assertEqual(set(audit.effective_min_cells),{3})
        self.assertEqual(audit.status.tolist(),['retained','low_detection','low_detection','zero_counts','low_detection'])

    def test_no_silent_cell_removal_or_empty_discovery(self):
        for matrix in [sparse.csr_matrix([[1,0],[0,1],[1,0]]),sparse.eye(2,format='csr')]:
            with self.assertRaises(ValueError):filter_counts(matrix,['g1','g2'],2)
        for invalid in [True,0,1.5,'20']:
            with self.assertRaises(ValueError):validate_thresholds(invalid)
        for invalid in [True,-.1,1.1,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):validate_thresholds(1,invalid)

    def test_detection_flags_do_not_change_ranked_scores(self):
        scores=pd.DataFrame({'p1':[5,4,-1],'p2':[1,-2,3]},index=['001','NA','common'])
        tpm=pd.DataFrame({'p1':[1,10,5],'p2':[2,3,9]},index=scores.index)
        metrics=pd.DataFrame({'gene':scores.index,'n_detected_cells':[1,2,4],'detection_fraction':[.25,.5,1.],'total_counts':[1,4,7]})
        with tempfile.TemporaryDirectory() as tmp,warnings.catch_warnings(record=True) as notices:
            table=write_top_gene_tables(scores,tpm,metrics,Path(tmp),3,2)
            self.assertTrue(notices)
            self.assertEqual(table.query("program=='p1'").gene.tolist(),['001','NA','common'])
            self.assertEqual(table.query("program=='p1'").score.tolist(),[5,4,-1])
            self.assertEqual(table.query("program=='p1'").low_detection.tolist(),[True,False,False])
            alternative=pd.read_csv(Path(tmp)/'top_genes_by_tpm.tsv',sep='\t',dtype={'gene':str},keep_default_na=False)
            self.assertEqual(alternative.query("program=='p1'").gene.tolist(),['NA','common','001'])
            self.assertTrue((Path(tmp)/'top_gene_detection_summary.tsv').is_file())

    def test_parameters_bind_explicit_filters_and_preserve_legacy(self):
        cfg={'cnmf':{'components':[2,3]}}
        self.assertNotIn('gene_filter',preparation_parameters(cfg))
        cfg['gene_filter']={'min_cells':2,'min_fraction':.1}
        before=preparation_parameters(cfg)
        cfg['gene_filter']['min_cells']=3
        self.assertNotEqual(before,preparation_parameters(cfg))

    def test_contract_import_does_not_initialize_scientific_libraries(self):
        command="import sys;sys.path.insert(0,sys.argv[1]);import cnmf_discovery;assert 'numpy' not in sys.modules;assert 'scipy' not in sys.modules"
        done=subprocess.run([sys.executable,'-c',command,str(ROOT/'toolkit/python')],capture_output=True,text=True)
        self.assertEqual(done.returncode,0,done.stderr)
