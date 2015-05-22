from libc.stdio cimport FILE


cdef extern from "freesasa.h":
    ctypedef struct freesasa:
        pass
    ctypedef enum freesasa_algorithm:
        FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY
    ctypedef enum freesasa_class:
        FREESASA_POLAR, FREESASA_APOLAR,
        FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
    ctypedef enum freesasa_verbosity:
        FREESASA_V_NORMAL, FREESASA_V_SILENT

    cdef int FREESASA_NAME_LIMIT
    cdef int FREESASA_DEF_PROBE_RADIUS
    cdef int FREESASA_DEF_SR_N
    cdef double FREESASA_DEF_LR_D
    cdef int FREESASA_SUCCESS
    cdef int FREESASA_FAIL
    cdef int FREESASA_WARN
    
    # init
    freesasa* freesasa_new()
    void freesasa_free(freesasa *self)

    # calculate
    void freesasa_copy_param(freesasa *target, const freesasa *source)
    int freesasa_calc_coord(freesasa  *self, const double *coord,
                            const double *r, size_t n)
    int freesasa_calc_pdb(freesasa *self, FILE *pdb_file)
    int freesasa_calc_atoms(freesasa *self, const double *coord,
                            const char **residueNames, 
                            const char **atomNames, size_t n)
    int freesasa_link_coord(freesasa *self, const double *coord,
                            double *radii, size_t n)
    int freesasa_refresh(freesasa *self)

    # protein info
    double freesasa_radius(const char *residueName, const char *atomName)
    size_t freesasa_n_atoms(const freesasa *self)

    # calculation settings
    int freesasa_set_algorithm(freesasa *s, freesasa_algorithm alg)
    freesasa_algorithm freesasa_get_algorithm(const freesasa *s)
    const char* freesasa_algorithm_name(const freesasa *self)
    int freesasa_set_probe_radius(freesasa *self,double r)
    double freesasa_get_probe_radius(const freesasa *self)
    int freesasa_set_sr_points(freesasa *self, int n)
    int freesasa_get_sr_points(const freesasa *self)
    int freesasa_set_lr_delta(freesasa *self, double delta)
    double freesasa_get_lr_delta(const freesasa *self)
    int freesasa_set_nthreads(freesasa *self,int n)
    int freesasa_get_nthreads(const freesasa *self)
    void freesasa_set_proteinname(freesasa *self,const char *name)
    const char* freesasa_get_proteinname(const freesasa *self)

    # access results
    double freesasa_area_total(const freesasa *self)
    double freesasa_area_class(const freesasa *self, freesasa_class c)
    double freesasa_area_residue(const freesasa *self, const char *residueName)
    double freesasa_area_atom(const freesasa *self, int atom)
    const double* freesasa_area_atom_array(const freesasa *self)
    double freesasa_radius_atom(const freesasa *self, int atom)
    const double* freesasa_radius_atom_array(const freesasa *self)

    # write results to output
    int freesasa_write_pdb(const freesasa *self, FILE *output)
    int freesasa_per_residue_type(const freesasa *self, FILE *output)
    int freesasa_per_residue(const freesasa *self, FILE *output)

    int freesasa_set_verbosity(freesasa_verbosity v)
    freesasa_verbosity freesasa_get_verbosity() 
