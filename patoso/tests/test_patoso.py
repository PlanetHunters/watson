import os
import shutil
import unittest
from patoso.patoso import Patoso


class TestsPatoso(unittest.TestCase):
    def test_vetting_by_params(self):
        object_dir = "TIC25155310_[1,_2]"
        vetting_dir = object_dir + "/vetting_0"
        try:
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, 1, 0.07571,
                                       cadence=120)
            files_in_dir = os.listdir(vetting_dir)
            assert files_in_dir == 9
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)


if __name__ == '__main__':
    unittest.main()
