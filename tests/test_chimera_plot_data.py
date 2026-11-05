import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_plot_data import export_contacts


class ContactPresentation(unittest.TestCase):
    def test_exact_sparse_round_trip_and_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'junction/hic').mkdir(parents=True)
            self.assertEqual(export_contacts(root,{'id':'C1'}),[])
            matrix=np.array([[5,0,2],[0,0,0],[2,0,8]],dtype=np.uint32)
            path=root/'junction/hic/m.npz'
            np.savez_compressed(path,matrix=matrix,lo=1200,resolution=10000)
            summary={'maps':[dict(event='C1',library='only_library',resolution=10000,lo=1200,file='m.npz')]}
            (path.parent/'summary.json').write_text(json.dumps(summary))
            spec=export_contacts(root,{'id':'C1'})[0]
            values=np.loadtxt(gzip.open(root/'junction'/spec['file'],'rt'),skiprows=1,dtype=int)
            restored=np.zeros((spec['bins'],spec['bins']),dtype=int)
            for row,col,count in values:restored[row,col]=count
            np.testing.assert_array_equal(restored,matrix)
            self.assertEqual((spec['zero_cells'],spec['nonzero_cells'],spec['total']),(5,4,17))
            summary['maps'][0]['lo']=999
            (path.parent/'summary.json').write_text(json.dumps(summary))
            with self.assertRaisesRegex(ValueError,'disagree'):export_contacts(root,{'id':'C1'})


if __name__=='__main__':unittest.main()
