cdef extern from "freesasa.h":
     ctypedef struct freesasa_t:
         pass
     extern freesasa_t* freesasa_init()
