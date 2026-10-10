"""Verify one-format choices and rejection before analysis in both runtimes."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'toolkit/python'))
from figure_output import figure_format
from scrna_runtime import resolved_rscript

CASES=[({},'pdf'),({'output':{'figure_format':'png'}},'png'),({'plots':{'figure_format':'pdf'}},'pdf'),
       ({'output':{'preview_png':True}},'png'),({'output':{'figure_format':'PDF'},'enrichment':{'plot_format':'pdf'}},'pdf'),
       ({'output':{'figure_format':'both'}},'error'),({'output':{'figure_format':['pdf','png']}},'error'),
       ({'output':{'figure_format':['pdf']}},'error'),({'output':{'figure_format':'pdf'},'plots':{'figure_format':'png'}},'error'),
       ({'output':{'figure_format':'pdf','preview_png':True}},'error'),({'output':{'preview_png':1}},'error')]
class FigureOutput(unittest.TestCase):
    def test_python_and_r_agree_on_format_and_conflicts(self):
        for cfg,expected in CASES:
            with self.subTest(config=cfg):
                if expected=='error':
                    with self.assertRaises(ValueError): figure_format(cfg)
                else:self.assertEqual(figure_format(cfg),expected)
        r=resolved_rscript('03-scrna-review-qc',{})
        if not r:self.skipTest('QC R runtime unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'cases.json';p.write_text(json.dumps([c for c,_ in CASES]))
            driver=Path(tmp)/'check.R'
            driver.write_text('args<-commandArgs(TRUE);source(args[1]);x<-jsonlite::fromJSON(args[2],simplifyVector=FALSE);cat(jsonlite::toJSON(vapply(x,function(c)tryCatch(figure_format(c),error=function(e)"error"),character(1))))')
            done=subprocess.run([str(r),str(driver),str(ROOT/'toolkit/R/figure_output.R'),str(p)],capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stderr)
            self.assertEqual(json.loads(done.stdout),[e for _,e in CASES])
