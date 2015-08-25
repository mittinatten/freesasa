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

defaultParameters = {'algorithm' : FREESASA_SHRAKE_RUPLEY, #freesasa_default_parameters.alg,
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
            self._c_param.probe_radius = r

      def probeRadius(self):
            return self._c_param.probe_radius

      def setNPoints(self,n):
            assert(n in [20,50,100,200,500,1000,2000,5000])
            self._c_param.shrake_rupley_n_points = n

      def nPoints(self):
            return self._c_param.shrake_rupley_n_points

      def setDelta(self,delta):
            assert(delta > 0)
            self._c_param.lee_richards_delta = delta

      def delta(self):
            return self._c_param.lee_richards_delta

      def setNThreads(self,n):
            assert(n>0)
            self._c_param.n_threads = n

      def nThreads(self):
            return self._c_param.n_threads

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_parameters **p = <freesasa_parameters**> ptr2ptr
            p[0] = &self._c_param

# wrapper for freesasa_result
cdef class Result:
      cdef freesasa_result* _c_result

      def __cinit__ (self):
            self._c_result = NULL

      def __dealloc__(self):
            if self._c_result is not NULL:
                  freesasa_result_free(self._c_result)

      def totalArea(self):
            assert(self._c_result is not NULL)
            return self._c_result.total

      def atomArea(self,i):
            assert(self._c_result is not NULL)
            assert(i < self._c_result.n_atoms)
            return self._c_result.sasa[i]

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_result **p = <freesasa_result**> ptr2ptr
            p[0] = self._c_result

# wrapper for freesasa_classifier (should at some point provide facilites to extend this)
cdef class Classifier:
      cdef freesasa_classifier* _c_classifier

      def __cinit__ (self,fileName=None):
            cdef FILE *config
            self._c_classifier = NULL
            if fileName is not None:
                  config = fopen(fileName,'r')
                  if config is NULL:
                        raise IOError("File '%s' could not be opened." % fileName)
                  self._c_classifier = freesasa_classifier_from_file(config)
                  if self._c_classifier is NULL:
                        raise Exception("Error parsing configuration in '%s'." % fileName)

      def classify(self,residueName,atomName):
            if self._c_classifier is not NULL:
                  return self._c_classifier.sasa_class(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.sasa_class(residueName,atomName,&freesasa_default_classifier)

      def radius(self,residueName,atomName):
            if self._c_classifier is not NULL:
                  return self._c_classifier.radius(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.radius(residueName,atomName,&freesasa_default_classifier)
                  

      def class2str(self,classIndex):
            if self._c_classifier is not NULL:
                  return self._c_classifier.class2str(classIndex,self._c_classifier)
            return freesasa_default_classifier.class2str(classIndex,&freesasa_default_classifier)

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_classifier **p = <freesasa_classifier**> ptr2ptr
            p[0] = self._c_classifier

      def __dealloc__(self):
            if self._c_classifier is not NULL:
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

      def nAtoms(self):
            assert(self._c_structure is not NULL)
            return freesasa_structure_n(self._c_structure)

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            p[0] = self._c_structure

      def __dealloc__(self):
            if self._c_structure is not NULL:
                  freesasa_structure_free(self._c_structure)

# calculate SASA values for structure

def calc(structure,parameters=None,classifier=None):
      cdef const freesasa_classifier *c = NULL
      cdef const freesasa_structure *s = NULL
      cdef const freesasa_parameters *p = NULL
      if classifier is not None:
            classifier._assign_ptr(<size_t>&c)
      if parameters is not None:
            parameters._assign_ptr(<size_t>&p)
      structure._assign_ptr(<size_t>&s)
      cdef double* radii = freesasa_structure_radius(s,c)
      if radii is NULL:
            raise Exception("Error calculating radii of atoms using supplied classifier.")

      r = Result()
      r._c_result = <freesasa_result*> freesasa_calc_structure(s, radii,p)
      free(radii)
      return r
            

# classify results of calculation, returns dictionary with class names
# as keys and areas for the corresponding groups of atoms are values
def classifyResults(result,structure,classifier=None):
      cdef const freesasa_classifier *c = NULL
      cdef const freesasa_structure *s = NULL
      cdef const freesasa_result *r = NULL
      if classifier is not None:
            classifier._assign_ptr(<size_t>&c)
      structure._assign_ptr(<size_t>&s)
      result._assign_ptr(<size_t>&r)
      cdef freesasa_strvp* strvp =  <freesasa_strvp*> freesasa_result_classify(r,s,c)
      ret = dict()
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

