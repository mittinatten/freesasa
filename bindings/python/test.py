import freesasa
import unittest
from exceptions import Exception
from freesasa import FreeSASAFail, FreeSASAWarn

def testNonExistentArea(f):
    try:
        print f.totalArea()
    except Exception as e:
        print e

class FreeSASATestCase(unittest.TestCase):
    def testCalcPDB(self):
        f = freesasa.FreeSASA()
        freesasa.setVerbosity(freesasa.silent)
        self.assertRaises(IOError,f.calcPDB,"")
        self.assertRaises(FreeSASAWarn,f.totalArea)
        self.assertRaises(FreeSASAWarn,f.areaOfClass,freesasa.polar)
        self.assertRaises(FreeSASAFail,f.calcPDB,"../../tests/data/empty_model.pdb")
        self.assertRaises(FreeSASAWarn,f.totalArea)
        self.assertRaises(FreeSASAWarn,f.areaOfClass,freesasa.polar)
        self.assertRaises(FreeSASAFail,f.calcPDB,"../../tests/data/dummy.pdb")
        self.assertRaises(FreeSASAWarn,f.totalArea)
        self.assertRaises(FreeSASAWarn,f.areaOfClass,freesasa.polar)
        freesasa.setVerbosity(freesasa.normal)
        f.calcPDB("../../tests/data/1ubq.pdb")
        self.assertTrue(abs(f.totalArea() - 4759.86095865) < 1e-7)
        self.assertTrue(abs(f.areaOfClass(freesasa.polar) - 2232.23038567) < 1e-7)
        self.assertTrue(abs(f.areaOfClass(freesasa.apolar) - 2527.63057298) < 1e-7)
        self.assertTrue(abs(f.areaOfClass(freesasa.nucleic)) < 1e-15)
        self.assertTrue(abs(f.areaOfClass(freesasa.unknown)) < 1e-15)
        
if __name__ == '__main__':
    unittest.main()

