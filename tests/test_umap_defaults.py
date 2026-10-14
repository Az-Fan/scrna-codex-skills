from pathlib import Path
import sys,unittest,warnings
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'toolkit/python'))
from figure_output import umap_marker_area,categorical_palette
class UMAPDefaults(unittest.TestCase):
    def test_area_and_overrides(self):
        self.assertGreater(umap_marker_area(60),umap_marker_area(10000))
        self.assertEqual(umap_marker_area(60,2),2)
        for value in [0,-1,True,float('inf'),float('nan')]:
            with self.assertRaises(ValueError):umap_marker_area(60,value)
    def test_unique_palette_and_alias_validation(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            labels=[str(x) for x in range(1000)]
            palette=categorical_palette(labels)
            self.assertEqual(len(set(palette.values())),1000)
            self.assertEqual(palette,categorical_palette(list(reversed(labels))))
        with self.assertRaises(ValueError):categorical_palette(['a','b'],{'a':'red','b':'#ff0000'})
