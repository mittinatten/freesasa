from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport free 
from libc.string cimport memcpy
from cpython cimport array
from cfreesasa cimport *


polar = FREESASA_POLAR
apolar = FREESASA_APOLAR 
nucleic = FREESASA_NUCLEICACID
unknown = FREESASA_CLASS_UNKNOWN
ShrakeRupley = FREESASA_SHRAKE_RUPLEY
LeeRichards = FREESASA_LEE_RICHARDS
silent = FREESASA_V_SILENT
normal = FREESASA_V_NORMAL

defaultParameters = {'algorithm' : freesasa_default_parameters.alg,
                     'probe-radius' : freesasa_default_parameters.probe_radius,
                     'n-points' :
                       freesasa_default_parameters.shrake_rupley_n_points,
                     'delta' :
                       freesasa_default_parameters.lee_richards_delta,
                     'n-threads' : freesasa_default_parameters.n_threads }

# wrapper for freesasa_parameters
cdef class Parameters:
      cdef freesasa_parameters _c_param
      def __init__(self,param=None):
            self._c_param = freesasa_default_parameters
            if param != None:
                  if 'algorithm' in param:    self.setAlgorithm(param['algorithm'])
                  if 'probe-radius' in param: self.setProbeRadius(param['probe-radius'])
                  if 'n-points' in param:     self.setNPoints(param['n-points'])
                  if 'delta' in param:        self.setDelta(param['delta'])
                  if 'n-threads' in param:    self.setNThreads(param['n-threads'])

      def setAlgorithm(self,alg):
            assert(alg in [ShrakeRupley, LeeRichards])
            self._c_param.alg = alg

      def algorithm(self):
            return self._c_param.alg

      def setProbeRadius(self,r):
            assert(r >= 0)
            self.probe_radius = r

      def probeRadius(self):
            return self._c_param.probe_radius

      def setNPoints(self,n):
            assert(n in [20,50,100,200,500,1000,2000,5000])
            self._c_param.shrake_rupley_n_points = n

      def NPoints(self):
            return self._c_param.shrake_rupley_n_points

      def setDelta(self,delta):
            assert(delta > 0)
            self._c_param.lee_richards_delta = delta

      def delta(self):
            return self._c_param.lee_richars_delta

      def setNThreads(self,n):
            assert(n>0)
            self._c_param.n_threads = n

      def nThreads(self):
            return self.c_param.n_threads

# wrapper for freesasa_result
cdef class Result:
      cdef freesasa_result* _c_result

      def __cinit__ (self):
            self._c_result = NULL

      def __dealloc__(self):
            freesasa_result_free(self._c_result)

      def totalArea(self):
            return self._c_result.total

      def atomArea(self,i):
            assert(i < self._c_result.n_atoms)
            return self._c_result.sasa[i]

# wrapper for freesasa_classifier (should at some point provide facilites to extend this)
cdef class Classifier:
      cdef freesasa_classifier* _c_classifier

      def __cinit__ (self,fileName=None):
            cdef FILE *config
            if fileName is None:
                  self._c_classifier[0] = freesasa_default_classifier
                  self.deallocFlag = False
            else:
                  config = fopen(fileName,'r')
                  if config is NULL:
                        raise IOError("File '%s' could not be opened." % fileName)
                  self._c_classifier = freesasa_classifier_from_file(config)
                  self.deallocFlag = True

      def __dealloc__(self):
            if self.deallocFlag is True:
                  freesasa_classifier_free(self._c_classifier)

# wrapper for freesasa_structure
cdef class Structure:
      cdef freesasa_structure* _c_structure

      def __cinit__(self,fileName,classifier=None,hetatm=0):
            self._c_structure = NULL
            cdef FILE *input
            input = fopen(fileName,'r')
            if input is NULL:
                  raise IOError("File '%s' could not be opened." % fileName)
            self._c_structure = freesasa_structure_from_pdb(input,hetatm)
            fclose(input)
            if self._c_structure is NULL:
                  raise Exception("Error reading '%s' as PDB-file." % fileName)

      def __dealloc__(self):
            freesasa_structure_free(self._c_structure)


# calculate SASA values for structure

def calc(structure,parameters=None,classifier=None):
      cdef freesasa_classifier *c = NULL
      if classifier is not None:
            c = <freesasa_classifier*>classifier._c_classifier
      cdef double* radii = freesasa_structure_radius(<freesasa_structure*>structure._c_structure,NULL)
      if radii is NULL:
            raise Exception("Error calculating radii of atoms using supplied classifier.")

      r = Result()
      r._c_result =  <freesasa_result*>\
          freesasa_calc_structure(<freesasa_structure*>structure._c_structure, 
                                   radii,
                                   <freesasa_parameters*>parameters._c_parameters)
      free(radii)
      return r
            


# classify results of calculation, returns dictionary with class names
# as keys and areas for the corresponding groups of atoms are values
def classifyResults(result,structure,classifier=None):
      if classifier is None:
            classifier = Classifier()
      cdef freesasa_strvp* strvp =  <freesasa_strvp*> \
            freesasa_result_classify(<freesasa_result*> result._c_result,
                                     <freesasa_structure*> structure._c_structure,
                                     <freesasa_classifier*> classifier._c_classifier)
      ret = {}
      for i in range(0,strvp.n):
            ret[strvp.string[i]] = strvp.value[i]
      return ret

# set verbosity for FreeSASA
def setVerbosity(verbosity):
      assert(verbosity in [silent, normal])
      res = freesasa_set_verbosity(verbosity)
      if res is FREESASA_FAIL:
            raise(Warning("Invalid verbosity level %d" % verbosity))

def getVerbosity():
    return freesasa_get_verbosity()

