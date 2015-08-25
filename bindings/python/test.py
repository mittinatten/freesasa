import freesasa
import unittest
import math
from exceptions import Exception

class FreeSASATestCase(unittest.TestCase):
    def testParameters(self):
        d = freesasa.defaultParameters
        p = freesasa.Parameters()
        self.assertTrue(p.algorithm() == d['algorithm'])
        self.assertTrue(p.probeRadius() == d['probe-radius'])
        self.assertTrue(p.nPoints() == d['n-points'])
        self.assertTrue(p.delta() == d['delta'])
        self.assertTrue(p.nThreads() == d['n-threads'])

        p.setAlgorithm(freesasa.ShrakeRupley)
        self.assertTrue(p.algorithm() == freesasa.ShrakeRupley)
        p.setAlgorithm(freesasa.LeeRichards)
        self.assertTrue(p.algorithm() == freesasa.LeeRichards)
        self.assertRaises(AssertionError,lambda: p.setAlgorithm(-10))

        p.setProbeRadius(1.5)
        self.assertTrue(p.probeRadius() == 1.5)
        self.assertRaises(AssertionError,lambda: p.setProbeRadius(-1))

        p.setNPoints(20)
        self.assertTrue(p.nPoints() == 20)
        self.assertRaises(AssertionError,lambda: p.setNPoints(21))

        p.setDelta(0.5)
        self.assertTrue(p.delta() == 0.5)
        self.assertRaises(AssertionError,lambda: p.setDelta(-0.5))

        p.setNThreads(2)
        self.assertTrue(p.nThreads() == 2)
        self.assertRaises(AssertionError, lambda: p.setNThreads(0))
        
    def testResult(self):
        r = freesasa.Result()
        self.assertRaises(AssertionError,lambda: r.totalArea())
        self.assertRaises(AssertionError,lambda: r.atomArea(0))

    def testClassifier(self):
        c = freesasa.Classifier()
        self.assertTrue(c.classify("ALA"," CB ") == freesasa.apolar)
        self.assertTrue(c.classify("ARG"," NH1") == freesasa.polar)
        self.assertTrue(c.radius("ALA"," CB ") == 2.00)

        freesasa.setVerbosity(freesasa.silent)
        self.assertRaises(Exception,lambda: freesasa.Classifier("data/err.config"))
        self.assertRaises(IOError,lambda: freesasa.Classifier(""))
        freesasa.setVerbosity(freesasa.normal)

        c = freesasa.Classifier("data/test.config")
        self.assertTrue(c.class2str(c.classify("AA","aa")) == "a")
        self.assertTrue(c.class2str(c.classify("BB","bb")) == "b")
        self.assertTrue(c.radius("AA","aa") == 1.0)
        self.assertTrue(c.radius("BB","bb") == 2.0)

        c = freesasa.Classifier("data/oons.config")
        self.assertTrue(c.radius("ALA"," CB ") == 2.00)

    def testStructure(self):
        self.assertRaises(IOError,lambda: freesasa.Structure("xyz#$%"))
        freesasa.setVerbosity(freesasa.silent)
        # test any file that's not a PDB file
        self.assertRaises(Exception,lambda: freesasa.Structure("data/err.config")) 
        self.assertRaises(Exception,lambda: freesasa.Structure("data/empty.pdb"))
        self.assertRaises(Exception,lambda: freesasa.Structure("data/empty_model.pdb"))
        freesasa.setVerbosity(freesasa.normal)

        s = freesasa.Structure("data/1ubq.pdb")
        self.assertTrue(s.nAtoms() == 602)

    def testCalc(self):
        s = freesasa.Structure("data/1ubq.pdb")
        r = freesasa.calc(s)
        self.assertTrue(math.fabs(r.totalArea() - 4759.86096) < 1e-5)

if __name__ == '__main__':
    unittest.main()

