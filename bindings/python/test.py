from freesasa import *
import unittest
import math
import os
from exceptions import Exception

# this class tests using derived classes to create custom Classifiers
class DerivedClassifier(Classifier):
    purePython = True

    def classify(self,residueName,atomName):
        return 'bla'

    def radius(self,residueName,atomName):
        return 10

class FreeSASATestCase(unittest.TestCase):
    def testParameters(self):
        d = defaultParameters
        p = Parameters()
        self.assertTrue(p.algorithm() == LeeRichards)
        self.assertTrue(p.algorithm() == d['algorithm'])
        self.assertTrue(p.probeRadius() == d['probe-radius'])
        self.assertTrue(p.nPoints() == d['n-points'])
        self.assertTrue(p.nSlices() == d['n-slices'])
        self.assertTrue(p.nThreads() == d['n-threads'])
        self.assertRaises(AssertionError,lambda: Parameters({'not-an-option' : 1}))
        self.assertRaises(AssertionError,lambda: Parameters({'n-slices' : 50, 'not-an-option' : 1}))
        self.assertRaises(AssertionError,lambda: Parameters({'not-an-option' : 50, 'also-not-an-option' : 1}))

        p.setAlgorithm(ShrakeRupley)
        self.assertTrue(p.algorithm() == ShrakeRupley)
        p.setAlgorithm(LeeRichards)
        self.assertTrue(p.algorithm() == LeeRichards)
        self.assertRaises(AssertionError,lambda: p.setAlgorithm(-10))

        p.setProbeRadius(1.5)
        self.assertTrue(p.probeRadius() == 1.5)
        self.assertRaises(AssertionError,lambda: p.setProbeRadius(-1))

        p.setNPoints(20)
        self.assertTrue(p.nPoints() == 20)
        self.assertRaises(AssertionError,lambda: p.setNPoints(0))

        p.setNSlices(10)
        self.assertTrue(p.nSlices() == 10)
        self.assertRaises(AssertionError,lambda: p.setNSlices(0))

        p.setNThreads(2)
        self.assertTrue(p.nThreads() == 2)
        self.assertRaises(AssertionError, lambda: p.setNThreads(0))
        
    def testResult(self):
        r = Result()
        self.assertRaises(AssertionError,lambda: r.totalArea())
        self.assertRaises(AssertionError,lambda: r.atomArea(0))

    def testClassifier(self):
        c = Classifier()
        self.assertTrue(c._isCClassifier())
        self.assertTrue(c.classify("ALA"," CB ") == apolar)
        self.assertTrue(c.classify("ARG"," NH1") == polar)
        self.assertTrue(c.radius("ALA"," CB ") == 1.88)

        setVerbosity(silent)
        self.assertRaises(Exception,lambda: Classifier("data/err.config"))
        self.assertRaises(IOError,lambda: Classifier(""))
        setVerbosity(normal)

        c = Classifier("data/test.config")
        self.assertTrue(c.classify("AA","aa") == "Polar")
        self.assertTrue(c.classify("BB","bb") == "Apolar")
        self.assertTrue(c.radius("AA","aa") == 1.0)
        self.assertTrue(c.radius("BB","bb") == 2.0)

        c = Classifier("share/oons.config")
        self.assertTrue(c.radius("ALA"," CB ") == 2.00)

        c = DerivedClassifier()
        self.assertTrue(not c._isCClassifier())
        self.assertTrue(c.radius("ALA"," CB ") == 10)
        self.assertTrue(c.radius("ABCDEFG","HIJKLMNO") == 10)
        self.assertTrue(c.classify("ABCDEFG","HIJKLMNO") == "bla")

    def testStructure(self):
        self.assertRaises(IOError,lambda: Structure("xyz#$%"))
        setVerbosity(silent)
        # test any file that's not a PDB file
        self.assertRaises(Exception,lambda: Structure("data/err.config")) 
        self.assertRaises(Exception,lambda: Structure("data/empty.pdb"))
        self.assertRaises(Exception,lambda: Structure("data/empty_model.pdb"))
        setVerbosity(normal)

        s = Structure("data/1ubq.pdb")
        self.assertTrue(s.nAtoms() == 602)
        self.assertTrue(s.radius(1) == 1.88)
        self.assertTrue(s.chainLabel(1) == 'A')
        self.assertTrue(s.atomName(1) == ' CA ')
        self.assertTrue(s.residueName(1) == 'MET')
        self.assertTrue(s.residueNumber(1) == '   1')

        s2 = Structure("data/1ubq.pdb",Classifier("share/oons.config"))
        self.assertTrue(s.nAtoms() == 602)
        self.assertTrue(math.fabs(s2.radius(1) - 2.0) < 1e-5)

        s2 = Structure("data/1ubq.pdb",Classifier("share/protor.config"))
        for i in range (0,601):
            self.assertTrue(math.fabs(s.radius(i)- s2.radius(i)) < 1e-5)

        self.assertRaises(Exception,lambda: Structure("data/1ubq.pdb","data/err.config"))

        s = Structure()
        s.addAtom(' CA ','ALA','   1','A',1,1,1)
        self.assertTrue(s.nAtoms() == 1)
        self.assertTrue(s.atomName(0) == ' CA ')
        self.assertTrue(s.residueName(0) == 'ALA')
        self.assertTrue(s.residueNumber(0) == '   1')
        self.assertTrue(s.chainLabel(0) == 'A')
        self.assertTrue(s.nAtoms() == 1)
        x, y, z = s.coord(0)
        self.assertTrue(x == 1 and y ==1 and z ==1)
        s.addAtom(' CB ','ALA',2,'A',2,1,1)
        self.assertTrue(s.nAtoms() == 2)
        self.assertTrue(s.residueNumber(1) == '2')

        self.assertRaises(AssertionError, lambda: s.atomName(3))
        self.assertRaises(AssertionError, lambda: s.residueName(3))
        self.assertRaises(AssertionError, lambda: s.residueNumber(3))
        self.assertRaises(AssertionError, lambda: s.chainLabel(3))
        self.assertRaises(AssertionError, lambda: s.coord(3))
        self.assertRaises(AssertionError, lambda: s.radius(3))

        s.setRadiiWithClassifier(Classifier())
        self.assertTrue(s.radius(0) == 1.88)
        self.assertTrue(s.radius(1) == 1.88)

        s.setRadiiWithClassifier(DerivedClassifier())
        self.assertTrue(s.radius(0) == s.radius(1) == 10.0)
        
        s.setRadii([1.0,3.0])
        self.assertTrue(s.radius(0) == 1.0)
        self.assertTrue(s.radius(1) == 3.0)

        s.setRadius(0, 10.0)
        self.assertTrue(s.radius(0) == 10.0);

        self.assertRaises(AssertionError,lambda: s.setRadius(2,10));
        self.assertRaises(AssertionError,lambda: s.setRadii([1]))
        self.assertRaises(AssertionError,lambda: s.setRadii([1,2,3]))

        self.assertRaises(AssertionError,lambda: s.atomName(2))
        self.assertRaises(AssertionError,lambda: s.residueName(2))
        self.assertRaises(AssertionError,lambda: s.residueNumber(2))
        self.assertRaises(AssertionError,lambda: s.chainLabel(2))

        setVerbosity(nowarnings)
        s = Structure("data/1d3z.pdb",None,{'hydrogen' : True})
        self.assertTrue(s.nAtoms() == 1231)

        s = Structure("data/1d3z.pdb",None,{'hydrogen' : True, 'join-models' : True})
        self.assertTrue(s.nAtoms() == 12310)

        s = Structure("data/1ubq.pdb",None,{'hetatm' : True})
        self.assertTrue(s.nAtoms() == 660)

        s = Structure("data/1d3z.pdb",None,{'hydrogen' : True, 'skip-unknown' : True})
        self.assertTrue(s.nAtoms() == 602)
        
        setVerbosity(silent)
        self.assertRaises(Exception, lambda : Structure("data/1d3z.pdb",None,{'hydrogen' : True, 'halt-at-unknown' : True}))
        setVerbosity(normal)

    def testStructureArray(self):
        # default separates chains, only uses first model (129 atoms per chain)
        ss = structureArray("data/2jo4.pdb")
        self.assertTrue(len(ss) == 4)
        for s in ss:
            self.assertTrue(s.nAtoms() == 129)

        # include all models, separate chains, and include hydrogen and hetatm (286 atoms per chain)
        setVerbosity(nowarnings)
        ss = structureArray("data/2jo4.pdb",{'separate-models' : True,
                                             'hydrogen' : True,
                                             'hetatm' : True,
                                             'separate-chains' : True})
        self.assertTrue(len(ss) == 4*10)
        for s in ss:
            self.assertTrue(s.nAtoms() == 286)

        # include all models, and include hydrogen and hetatm (286 atoms per chain)
        ss = structureArray("data/2jo4.pdb",{'separate-models' : True,
                                             'hydrogen' : True,
                                             'hetatm' : True})
        self.assertTrue(len(ss) == 10)
        for s in ss:
            self.assertTrue(s.nAtoms() == 286*4)
        setVerbosity(normal)

        # check that the structures initialized this way can be used for calculations
        ss = structureArray("data/1ubq.pdb")
        self.assertTrue(len(ss) == 1)
        self.assertTrue(ss[0].nAtoms() == 602)
        result = calc(ss[0],Parameters({'algorithm' : ShrakeRupley}))
        self.assertTrue(math.fabs(result.totalArea() - 4834.716265) < 1e-5)

        # Test exceptions
        setVerbosity(silent)
        self.assertRaises(AssertionError,lambda: structureArray(None))
        self.assertRaises(IOError,lambda: structureArray(""))
        self.assertRaises(Exception,lambda: structureArray("data/err.config"))
        self.assertRaises(AssertionError,lambda: structureArray("data/2jo4.pdb",{'not-an-option' : True}))
        self.assertRaises(AssertionError,
                          lambda: structureArray("data/2jo4.pdb",
                                                 {'not-an-option' : True, 'hydrogen' : True}))
        self.assertRaises(AssertionError,
                          lambda: structureArray("data/2jo4.pdb",
                                                 {'hydrogen' : True}))
        setVerbosity(normal)

    def testCalc(self):
        # test default settings
        structure = Structure("data/1ubq.pdb")
        result = calc(structure,Parameters({'algorithm' : ShrakeRupley}))
        self.assertTrue(math.fabs(result.totalArea() - 4834.716265) < 1e-5)
        sasa_classes = classifyResults(result,structure)
        self.assertTrue(math.fabs(sasa_classes['Polar'] - 2515.821238) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2318.895027) < 1e-5)

        # test L&R
        result = calc(structure,Parameters({'algorithm' : LeeRichards, 'n-slices' : 20}))
        sasa_classes = classifyResults(result,structure)
        self.assertTrue(math.fabs(result.totalArea() - 4804.055641) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Polar'] - 2504.217302) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2299.838339) < 1e-5)

        # test extending Classifier with derived class
        sasa_classes = classifyResults(result,structure,DerivedClassifier())
        self.assertTrue(math.fabs(sasa_classes['bla'] - 4804.055641) < 1e-5)
        
        ## test calculating with user-defined classifier ##
        classifier = Classifier("share/oons.config")
        # classifier passed to assign user-defined radii, could also have used setRadiiWithClassifier()
        structure = Structure("data/1ubq.pdb",classifier) 
        result = calc(structure,Parameters({'algorithm' : ShrakeRupley}))
        self.assertTrue(math.fabs(result.totalArea() - 4779.5109924) < 1e-5)
        sasa_classes = classifyResults(result,structure,classifier) # classifier passed to get user-classes
        self.assertTrue(math.fabs(sasa_classes['Polar'] - 2236.9298941) < 1e-5)
        self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2542.5810983) < 1e-5)

    def testCalcCoord(self):
        # one unit sphere
        radii = [1]
        coord = [0,0,0]
        parameters = Parameters()
        parameters.setNSlices(5000)
        parameters.setProbeRadius(0)
        result = calcCoord(coord, radii, parameters)
        self.assertTrue(math.fabs(result.totalArea() - 4*math.pi) < 1e-3)

        # two separate unit spheres
        radii = [1,1]
        coord = [0,0,0, 4,4,4]
        result = calcCoord(coord, radii, parameters)
        self.assertTrue(math.fabs(result.totalArea() - 2*4*math.pi) < 1e-3)

        self.assertRaises(AssertionError,
                          lambda: calcCoord(radii, radii))

    def testSelectArea(self):
        structure = Structure("data/1ubq.pdb")
        result = calc(structure,Parameters({'algorithm' : ShrakeRupley}))
        # will only test that this gets through to the C interface,
        # extensive checking of the parser is done in the C unit tests
        selections = selectArea(('s1, resn ala','s2, resi 1'),structure,result)
        self.assertTrue(math.fabs(selections['s1'] - 118.35) < 0.1)
        self.assertTrue(math.fabs(selections['s2'] - 50.77) < 0.1)

    def testBioPDB(self):
        try:
            from Bio.PDB import PDBParser
        except ImportError:
            print "Can't import Bio.PDB, tests skipped"
            pass
        else:
            parser = PDBParser()
            bp_structure = parser.get_structure("Ubiquitin","data/1ubq.pdb")
            s1 = structureFromBioPDB(bp_structure)
            s2 = Structure("data/1ubq.pdb")
            self.assertTrue(s1.nAtoms() == s2.nAtoms())

            for i in range(0, s2.nAtoms()):
                self.assertTrue(s1.radius(i) == s2.radius(i))
                # there can be tiny errors here
                self.assertTrue(math.fabs(s1.coord(i)[0] - s2.coord(i)[0]) < 1e-5)
                self.assertTrue(math.fabs(s1.coord(i)[1] - s2.coord(i)[1]) < 1e-5)
                self.assertTrue(math.fabs(s1.coord(i)[2] - s2.coord(i)[2]) < 1e-5)

            # because Bio.PDB structures will have slightly different
            # coordinates (due to rounding errors) we set the
            # tolerance as high as 1e-3
            result = calc(s1, Parameters({'algorithm' : LeeRichards, 'n-slices' : 20}))
            print result.totalArea()
            self.assertTrue(math.fabs(result.totalArea() - 4804.055641) < 1e-3)
            sasa_classes = classifyResults(result, s1)
            self.assertTrue(math.fabs(sasa_classes['Polar'] - 2504.217302) < 1e-3)
            self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2299.838339) < 1e-3)

            result, sasa_classes = calcBioPDB(bp_structure, Parameters({'algorithm' : ShrakeRupley}))
            self.assertTrue(math.fabs(result.totalArea() - 4834.716265) < 1e-3)
            self.assertTrue(math.fabs(sasa_classes['Polar'] - 2515.821238) < 1e-3)
            self.assertTrue(math.fabs(sasa_classes['Apolar'] - 2318.895027) < 1e-3)
            print result.totalArea()


if __name__ == '__main__':
    # make sure we're in the right directory (if script is called from
    # outside the directory)
    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    os.chdir(dirname)
    unittest.main()

