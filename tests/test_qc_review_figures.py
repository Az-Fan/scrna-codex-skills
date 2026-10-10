"""Observable QC display contracts for heterogeneous Seurat metadata."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'toolkit/python'))
from scrna_runtime import resolved_rscript
RSCRIPT = resolved_rscript('03-scrna-review-qc', {})


class QCReviewFigures(unittest.TestCase):
    def test_display_units_missing_metrics_pagination_and_preserved_diagnostics(self):
        if not RSCRIPT:
            self.skipTest('QC R runtime unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            creator = root / 'create.R'
            creator.write_text('''args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(Seurat))
x <- matrix(rep(1:20,80),20,80,dimnames=list(paste0("g",1:20),paste0("c",1:80)))
o <- CreateSeuratObject(x)
o$sample <- rep(c("s1","s2"),each=40)
o$seurat_clusters <- rep(c("2","10"),each=40)
o$cell_type <- rep(c("EC","Fib"),each=40)
o$mito_frac <- factor(rep(c("0.01","0.02"),each=40))
o$doublet_score <- rep(c(NA,.2),each=40)
o$nuclear_frac <- seq(.01,.8,length.out=80)
xy <- cbind(seq(-2,2,length.out=80),sin(seq_len(80)))
rownames(xy)<-colnames(o); colnames(xy)<-c("UMAP_1","UMAP_2")
o[["umap"]]<-CreateDimReducObject(embeddings=xy,key="UMAP_",assay="RNA")
if(args[2]=="percent") {
 o$percent.mt<-as.numeric(as.character(o$mito_frac))*100; o$mito_frac<-NULL
 o$UMAP_1<-xy[,1]; o$UMAP_2<-xy[,2]; o[["umap"]]<-NULL
 o$nCount_RNA[1]<-0; o$nFeature_RNA[2]<-Inf
}
if(args[2]=="secondary") {
 o$sample<-"s1"; o$nCount_RNA<-NULL; o$nFeature_RNA<-NULL
 o$mito_frac<-NULL; o$doublet_score<-NULL; o$cell_type<-NULL
 o$seurat_clusters<-NULL; o[["umap"]]<-NULL
}
if(args[2]=="many") {
 o$sample<-paste0("s",rep(1:9,length.out=80))
 o$seurat_clusters<-as.character(rep(0:24,length.out=80))
}
saveRDS(o,args[1])
''')
            for mode in ['fraction', 'percent', 'secondary', 'many']:
                with self.subTest(mode=mode):
                    obj = root / (mode + '.rds')
                    made = subprocess.run([str(RSCRIPT),str(creator),str(obj),mode],capture_output=True,text=True)
                    self.assertEqual(made.returncode,0,made.stderr)
                    before = hashlib.sha256(obj.read_bytes()).hexdigest()
                    out = root / mode
                    cfg = {'input':{'object':str(obj)},'metadata':{'sample':'sample'},'output_dir':str(out),
                           'output':{'detail_level':'full' if mode=='fraction' else 'compact'},
                           'display':{'sample_labels':{'s1':'A','s2':'B'},'focus_clusters':['2']}}
                    if mode=='fraction':
                        cfg['approved_rules']={'mito_frac':{'upper':.05}}
                        cfg['approved_rules_source']='synthetic approved decision for regression test'
                    path = root / (mode + '.json'); path.write_text(json.dumps(cfg))
                    done = subprocess.run([str(RSCRIPT),str(ROOT/'toolkit/R/review_qc.R'),str(path)],capture_output=True,text=True)
                    self.assertEqual(done.returncode,0,done.stderr)
                    self.assertEqual(before,hashlib.sha256(obj.read_bytes()).hexdigest())
                    methods=json.loads((out/'_provenance/figure_methods.json').read_text())
                    self.assertEqual(methods['cells_removed'],0)
                    self.assertFalse(methods['threshold_candidates_applied'])
                    self.assertEqual(methods['mitochondrial_display_unit'],'percent')
                    self.assertEqual(methods['atlas_pages'],5 if mode=='many' else 3)
                    for preview in methods['previews']:
                        self.assertTrue((out/preview).is_file())
                    self.assertGreater((out/'qc_supplement.pdf').stat().st_size,1000)
                    with (out/'threshold_review.tsv').open() as handle:
                        rows=list(csv.DictReader(handle,delimiter='\t'))
                    self.assertTrue(all(not r['approval'] and not r['decision'] for r in rows))
                    if mode in ['fraction','percent']:
                        self.assertEqual(methods['mitochondrial_source_unit'],mode)
                        with (out/'sample_display_summary.tsv').open() as handle:
                            display=list(csv.DictReader(handle,delimiter='\t'))
                        self.assertEqual([float(r['mito_frac']) for r in display],[1.,2.])
                        mito=next(r for r in rows if r['metric']=='mito_frac')
                        self.assertAlmostEqual(float(mito['candidate_upper']),.02 if mode=='fraction' else 2.)
                        self.assertIn('umap_annotation_qc_context.png',methods['supplemental_diagnostics'])
                        with (out/'cluster_sample_composition.tsv').open() as handle:
                            composition=list(csv.DictReader(handle,delimiter='\t'))
                        for cl in ['2','10']:
                            self.assertAlmostEqual(sum(float(r['fraction']) for r in composition if r['cluster']==cl),1.)
                        if mode=='fraction':
                            self.assertTrue((out/'details/plot_status.tsv').is_file())
                            self.assertTrue((out/'details/qc_distribution_by_celltype.png').is_file())
                        else:
                            self.assertEqual(methods['umap_reduction'],'metadata:UMAP_1/UMAP_2')
                    if mode=='secondary':
                        self.assertFalse((out/'qc_summary_by_cluster.tsv').exists())
                        self.assertFalse((out/'details').exists())
