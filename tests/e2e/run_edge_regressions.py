#!/usr/bin/env python3
"""Verify the second audit's input, design, task and failure-path regressions.

Run after run_fixture_e2e.py and supply that run's directory as --baseline.
All new inputs and outputs are isolated; no source fixture is rewritten.
"""
import argparse
import copy
import csv
import hashlib
import gzip
import shutil
import json
import subprocess
import sys
from pathlib import Path


def tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--env-root", type=Path, default=Path('/home/faz_laptop/projects/scrna_envs'))
    parser.add_argument("--start-at", help="Resume at a case while retaining earlier verified outputs")
    args = parser.parse_args()
    repo, baseline, root = args.repo.resolve(), args.baseline.resolve(), args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=bool(args.start_at))
    installed = root / "installed"
    subprocess.run([sys.executable, str(repo / "scripts/install_skills.py"), "--target", str(installed), "--force"], check=True, stdout=subprocess.DEVNULL)
    rscript = args.env_root / "02-annotation/.pixi/envs/default/bin/Rscript"
    fixture, multi = repo / "tests/fixtures/tiny_scrna.rds", repo / "tests/fixtures/tiny_scrna_multilayer.rds"
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    before = {str(p): sha(p) for p in [fixture, multi]}
    prepare = root / "prepare.R"
    prepare.write_text('''args <- commandArgs(TRUE); root <- args[1]; repo <- args[2]
      suppressPackageStartupMessages(library(Seurat))
      obj <- readRDS(file.path(repo,"tests/fixtures/tiny_scrna.rds"))
      counts <- LayerData(obj,assay="RNA",layer="counts")
      tenx <- file.path(root,"tenx"); dir.create(tenx)
      Matrix::writeMM(counts,file.path(tenx,"matrix.mtx"))
      write.table(data.frame(rownames(counts),rownames(counts),"Gene Expression"),file.path(tenx,"features.tsv"),sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      write.table(colnames(counts),file.path(tenx,"barcodes.tsv"),sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      alt <- obj; alt[["alternative"]] <- CreateAssay5Object(counts=counts[1:10,]); DefaultAssay(alt) <- "alternative"
      saveRDS(alt,file.path(root,"alternate_default.rds"))
      one <- subset(readRDS(file.path(repo,"tests/fixtures/tiny_scrna_multilayer.rds")),cells=colnames(obj)[1:20])
      saveRDS(one,file.path(root,"single_split.rds"))
      collision <- obj; collision$cell_type <- ifelse(collision$cell_type=="Endothelial","T cell","T/cell")
      saveRDS(collision,file.path(root,"collision.rds"))
      obj[["alternate_embedding"]] <- CreateDimReducObject(embeddings=Embeddings(obj,"umap")*3,key="ALT_",assay="RNA")
      saveRDS(obj,file.path(root,"two_embeddings.rds"))
      for (kind in c("na","blank")) {
        d <- data.frame(cluster=c("0","1"),broad=c(if(kind=="na") NA else "   ","Stromal"),fine=c("Endothelial","Fibroblast"),decision="confirmed")
        write.table(d,file.path(root,paste0("annotation_",kind,".tsv")),sep="\\t",quote=FALSE,row.names=FALSE)
      }
      x <- read.delim(file.path(repo,"tests/fixtures/cell_abundance_counts.tsv"))
      ids <- unique(x$sample_label); x$age <- match(x$sample_label,ids)+40
      x$donor <- rep(c(101,102,103,104),2)[match(x$sample_label,ids)]
      write.table(x,file.path(root,"age.tsv"),sep="\\t",quote=FALSE,row.names=FALSE)
    ''')
    if not args.start_at:
        subprocess.run([str(rscript), str(prepare), str(root), str(repo)], check=True, stdout=subprocess.DEVNULL)
        compressed = root / 'tenx_gz'; compressed.mkdir()
        for path in (root/'tenx').iterdir():
            with path.open('rb') as src, gzip.open(compressed/(path.name+'.gz'),'wb') as dest:
                shutil.copyfileobj(src,dest)
    results = []
    started = args.start_at is None

    def should_run(label):
        return started or label == args.start_at

    def config(skill):
        return json.loads((baseline / "configs" / (skill + ".json")).read_text())

    def run(skill, label, cfg, expected_error=None):
        nonlocal started
        if not should_run(label):
            out = root/label
            manifest = json.loads((out/'_provenance/run_manifest.json').read_text())
            results.append(dict(case=label,exit_status=manifest['exit_status'],status='passed',retained_execution=True))
            return out
        started = True
        cfg = copy.deepcopy(cfg); cfg['output_dir'] = str(root / label)
        path = root / (label + '.json'); path.write_text(json.dumps(cfg, indent=2))
        done = subprocess.run([sys.executable, str(installed / skill / "scripts/run.py"), '--config', str(path), '--execute'], text=True, capture_output=True)
        (root / (label + '.log')).write_text(done.stdout + done.stderr)
        out = root / label
        if expected_error:
            assert done.returncode != 0 and expected_error in done.stdout + done.stderr + (out / '_provenance/run.log').read_text(), (label, done.stdout, done.stderr)
        else:
            assert done.returncode == 0, (label, done.stdout, done.stderr)
        manifest_path = out / '_provenance/run_manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            assert manifest['exit_status'] == done.returncode, (label, manifest)
        results.append(dict(case=label, exit_status=done.returncode, status='passed'))
        print('PASS', label, flush=True)
        return out

    skill = '01-scrna-standardize-input'
    cfg = config(skill); cfg['input'] = dict(path=str(root/'tenx'), format='10x_dir', sample_id='s1'); cfg['metadata'] = {'sample':'sample_id'}
    out = run(skill,'tenx_directory',cfg)
    record = json.loads((out/'_provenance/run_manifest.json').read_text())['input']
    assert record['type'] == 'directory' and len(record['members']) == 3
    assert len(tsv(out/'cell_metadata.tsv')) == 80
    cfg['input']['path'] = str(root/'tenx_gz')
    out = run(skill,'tenx_directory_gzip',cfg)
    assert len(tsv(out/'cell_metadata.tsv')) == 80

    skill = '02-scrna-calculate-qc-metrics'
    cfg = config(skill); cfg['metadata']['sample'] = 'sample_typo'; cfg['parallel']['workers'] = 1
    run(skill,'qc_sample_typo',cfg,'Configured metadata column not found: sample_typo')

    skill = '04-scrna-apply-qc-filter'
    cfg = config(skill); cfg['input']['object'] = str(root/'alternate_default.rds')
    out = run(skill,'filter_alternate_assay',cfg)
    subprocess.run([str(rscript), '-e', 'suppressPackageStartupMessages(library(Seurat)); x<-readRDS(commandArgs(TRUE)[1]); stopifnot(nrow(x)==20L,ncol(x)==79L)', str(out/'filtered_object.rds')], check=True)

    skill = '05-scrna-benchmark-integration'
    cfg = config(skill); cfg['input']['object'] = str(root/'two_embeddings.rds'); cfg['benchmark']['dims'] = [1,2]
    cfg['benchmark']['methods'] = [{'name':'precomputed','reduction':'umap'},{'name':'precomputed','reduction':'alternate_embedding'}]
    cfg['metrics'] = {'batch_removal':[],'biological_conservation':[]}; cfg['plots'] = []
    out = run(skill,'precomputed_independent',cfg)
    rows = tsv(out/'exchange/embedding_manifest.tsv'); assert len(rows) == 3
    assert len({row['scenario'] for row in rows}) == len({row['path'] for row in rows}) == 3
    first, second = [row for row in rows if row['method']=='precomputed']
    assert sha(Path(first['path'])) != sha(Path(second['path']))

    skill = '08-scrna-annotate-cells'
    for kind in ['na','blank']:
        cfg = json.loads((baseline/'configs'/ (skill+'.apply.json')).read_text())
        cfg['input']['object'] = str(fixture); cfg['input']['decisions'] = str(root/('annotation_'+kind+'.tsv'))
        cfg['annotation'] = {'broad_column':'broad','fine_column':'fine','decision_column':'decision'}
        run(skill,'annotation_'+kind,cfg,'labels cannot be missing or empty')
        assert not (root/('annotation_'+kind)/'annotated_object.qs').exists()

    skill = '10-scrna-score-programs'
    cfg = config(skill); cfg['input']['object'] = str(root/'single_split.rds')
    out = run(skill,'score_single_layer',cfg)
    task = json.loads((out/'_provenance/task_manifest.json').read_text())['vascular_program']
    assert task['n_cells'] == 20

    skill = '11-scrna-run-differential-analysis'
    cfg = config(skill); cfg['input']['object'] = str(root/'collision.rds')
    cfg['metadata']['population'] = 'cell_type'; cfg['population'] = {'mode':'all'}
    run(skill,'differential_collision',cfg,'task IDs collide after sanitization')
    cfg = config(skill); cfg['analysis']['min_samples_per_group'] = 3
    out = root/'failed_rerun'
    if should_run('failed_rerun'):
        out.mkdir(exist_ok=True); (out/'all_comparisons.tsv').write_text('old result\n')
    out = run(skill,'failed_rerun',cfg,'No differential-analysis task completed')
    assert not (out/'all_comparisons.tsv').exists()
    manifest = json.loads((out/'_provenance/run_manifest.json').read_text())
    assert (Path(manifest['previous_output'])/'all_comparisons.tsv').read_text() == 'old result\n'

    skill = '12-scrna-run-pathway-enrichment'
    cfg = config(skill); cfg['input'] = {'differential_tables':[{'path':str(baseline/'11-scrna-run-differential-analysis/all_comparisons.tsv'),'id':name} for name in ['a','b']]}
    out = run(skill,'enrichment_independent_tables',cfg)
    tasks = tsv(out/'task_status.tsv'); assert len(tasks) == len({x['task_id'] for x in tasks}) == 2
    assert all(task['status'] in {'completed','partial','empty'} for task in tasks)
    assert all((out/'comparisons'/task['task_id']/'standardized_input_table.tsv').exists() for task in tasks)

    skill = '13-scrna-test-cell-abundance'
    cfg = config(skill); cfg['analysis']['methods'] = ['propeller']; cfg['input']['counts_table'] = str(root/'age.tsv'); cfg['metadata']['covariates'] = ['age']
    out = run(skill,'abundance_continuous',cfg)
    assert tsv(out/'covariate_types.tsv')[0]['type'] == 'continuous'
    assert tsv(out/'task_status.tsv')[0]['status'] == 'completed'
    design = json.loads((out/'comparisons/case_vs_control/model_design.json').read_text())
    assert design['rank'] == 3 and design['n_samples'] == 8
    cfg['metadata']['covariates'] = ['donor']; cfg['metadata']['covariate_types'] = {'donor':'categorical'}
    out = run(skill,'abundance_categorical',cfg)
    assert tsv(out/'covariate_types.tsv')[0]['type'] == 'categorical'
    cfg = config(skill); cfg['analysis']['methods'] = ['propeller']; cfg['analysis']['min_cells_per_sample'] = 10000
    for policy, error in [('audit_only',None),('exclude_samples','No cell-abundance task completed'),('stop','No cell-abundance task completed')]:
        cfg['analysis']['min_cells_policy'] = policy
        out = run(skill,'min_cells_'+policy,cfg,error)
        audit = tsv(out/'design_audit.tsv'); assert all(row['passes_min_cells']=='FALSE' and row['min_cells_policy']==policy for row in audit)
        assert (tsv(out/'task_status.tsv')[0]['status']=='completed') == (policy=='audit_only')

    assert all(sha(Path(path)) == value for path,value in before.items())
    assert started, 'Requested start case was not found'
    (root/'regression-report.json').write_text(json.dumps({'status':'passed','source_fixtures_unchanged':True,'cases':results},indent=2))
    print('PASS all second-audit edge regressions')


if __name__ == '__main__':
    main()
