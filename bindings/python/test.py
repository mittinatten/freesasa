import freesasa
import unittest
from exceptions import Exception

class FreeSASATestCase(unittest.TestCase):
    def testParameters(self):
        d = freesasa.defaultParameters
        p = freesasa.Parameters()
        self.assertTrue(p.algorithm == d['algorithm'])
        self.assertRaises(AssertionError,p.setAlgorithm(-10))
        
if __name__ == '__main__':
    unittest.main()

