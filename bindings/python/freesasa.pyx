cimport cfreesasa
from libc.stdio cimport FILE, fopen, fclose

class FreeSASAFail(Exception):
      pass
      
class FreeSASAWarn(Warning):
      pass

polar = cfreesasa.FREESASA_POLAR
apolar = cfreesasa.FREESASA_APOLAR 
nucleic = cfreesasa.FREESASA_NUCLEICACID
unknown = cfreesasa.FREESASA_CLASS_UNKNOWN
ShrakeRupley = cfreesasa.FREESASA_SHRAKE_RUPLEY
LeeRichards = cfreesasa.FREESASA_LEE_RICHARDS
silent = cfreesasa.FREESASA_V_SILENT
normal = cfreesasa.FREESASA_V_NORMAL

def setVerbosity(verbosity):
     res = cfreesasa.freesasa_set_verbosity(verbosity)
     if res is cfreesasa.FREESASA_FAIL:
        raise(FreeSASAWarn("Invalid verbosity level %d" % verbosity));
def getVerbosity():
    return cfreesasa.freesasa_get_verbosity()

cdef class FreeSASA:
     cdef cfreesasa.freesasa_t* _c_freesasa

     def __cinit__(self):
         self._c_freesasa = cfreesasa.freesasa_init()
         if self._c_freesasa is NULL:
             raise MemoryError()

     def __dealloc__(self):
         if self._c_freesasa is not NULL:
            cfreesasa.freesasa_free(self._c_freesasa)

     def calcPDB(self, fileName):
         cdef FILE *input = fopen(fileName,'r')
         if input is NULL:
            raise IOError("File '%s' could not be opened." % fileName)
         result = cfreesasa.freesasa_calc_pdb(self._c_freesasa,input)
         fclose(input)
         if result is cfreesasa.FREESASA_FAIL:
            raise FreeSASAFail("Not able to calculate SASA, input valid?")
         if result is cfreesasa.FREESASA_WARN:     
            raise FreeSASAWarn("Warning: There may be errors in calculation, " \
                               "probably due to unrecognized input.")

     def totalArea(self):
         area = cfreesasa.freesasa_area_total(self._c_freesasa)
         if area < 0:
            raise FreeSASAWarn("Tried to access total area for a FreeSASA object " \
                               "before valid calculation has been performed")
         return area

     def areaOfClass(self,theClass):
         area = cfreesasa.freesasa_area_class(self._c_freesasa,theClass)
         if area < 0:
             raise FreeSASAWarn("Tried to access area for a FreeSASA object before " \
                               "valid calculation has been performed")
         return area
         