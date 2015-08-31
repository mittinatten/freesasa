import freesasa
import unittest
import math
import os
from exceptions import Exception

class DummyClassifier(freesasa.Classifier):
    def classify(self,residueName,atomName):
        return 0

    def radius(self,residueName,atomName):
        return 10

    def class2str(self,classIndex):
        return "bla"

class FreeSASATestCase(unittest.TestCase):
    def testParameters(self):
        d = freesasa.defaultParameters
        p = freesasa.Parameters()
        self.assertTrue(p.algorithm() == freesasa.ShrakeRupley)
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
        self.assertTrue(s.radius(1) == 2.0)
        self.assertTrue(s.chainLabel(1) == 'A')
        self.assertTrue(s.atomName(1) == ' CA ')
        self.assertTrue(s.residueName(1) == 'MET')
        self.assertTrue(s.residueNumber(1) == '   1')

        s2 = freesasa.Structure("data/1ubq.pdb",freesasa.Classifier("data/oons.config"))
        self.assertTrue(s.nAtoms() == 602)
        self.assertTrue(s.radius(1) == 2.0)

        for i in range (0,601):
            self.assertTrue(math.fabs(s.radius(i)- s2.radius(i)) < 1e-5)

        self.assertRaises(Exception,lambda: freesasa.Structure("data/1ubq.pdb","data/err.config"))

        s = freesasa.Structure()
        s.addAtom(' CA ','ALA','   1','A',1,1,1)
        self.assertTrue(s.nAtoms() == 1)
        self.assertTrue(s.atomName(0) == ' CA ');
        self.assertTrue(s.residueName(0) == 'ALA');
        self.assertTrue(s.residueNumber(0) == '   1');
        self.assertTrue(s.chainLabel(0) == 'A');
        self.assertRaises(AssertionError,lambda: s.radius(0))
        self.assertRaises(AssertionError,lambda: s.addAtom('CA','ALA','  12','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom('CA   ','ALA','  12','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ',' ALA','  12','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','AL','  12','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','ALA',' 12','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','ALA',' 12  ','A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','ALA',12345,'A',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','ALA','  12','AB',1,1,1))
        self.assertRaises(AssertionError,lambda: s.addAtom(' CA ','ALA','  12','',1,1,1))
        self.assertTrue(s.nAtoms() == 1)
        s.addAtom(' CB ','ALA',2,'A',2,1,1)
        self.assertTrue(s.nAtoms() == 2)
        self.assertTrue(s.residueNumber(1) == '   2')

        s.setRadiiWithClassifier()
        self.assertTrue(s.radius(0) == 2.0)
        self.assertTrue(s.radius(1) == 2.0)

        s.setRadiiWithClassifier(DummyClassifier())
        self.assertTrue(s.radius(0) == s.radius(1) == 10.0)
        
        s.setRadii([1.0,3.0])
        self.assertTrue(s.radius(0) == 1.0)
        self.assertTrue(s.radius(1) == 3.0)

        self.assertRaises(AssertionError,lambda: s.setRadii([1]))
        self.assertRaises(AssertionError,lambda: s.setRadii([1,2,3]))

        self.assertRaises(AssertionError,lambda: s.atomName(2))
        self.assertRaises(AssertionError,lambda: s.residueName(2))
        self.assertRaises(AssertionError,lambda: s.residueNumber(2))
        self.assertRaises(AssertionError,lambda: s.chainLabel(2))

    def testCalc(self):
        s = freesasa.Structure("data/1ubq.pdb")
        r = freesasa.calc(s)
        self.assertTrue(math.fabs(r.totalArea() - 4759.86096) < 1e-5)
        sasa_classes = freesasa.classifyResults(r,s)
        self.assertTrue(math.fabs(sasa_classes['Polar'] - 2232.23039) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2527.63057) < 1e-5)

        r = freesasa.calc(s,freesasa.Parameters({'algorithm' : freesasa.LeeRichards, 'delta' : 0.25}))
        sasa_classes = freesasa.classifyResults(r,s)
        self.assertTrue(math.fabs(r.totalArea() - 4728.26159) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Polar'] - 2211.41649) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2516.84510) < 1e-5)
        
        sasa_classes = freesasa.classifyResults(r,s,DummyClassifier())
        self.assertTrue(math.fabs(sasa_classes['bla'] - 4728.26159) < 1e-5)
        

if __name__ == '__main__':
    # make sure we're in the right directory (if script is called from
    # outside the directory)
    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    os.chdir(dirname)
    unittest.main()

