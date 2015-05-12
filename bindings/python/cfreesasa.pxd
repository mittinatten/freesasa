cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass

cdef extern from "freesasa.h":
    ctypedef struct freesasa_t:
        pass
   
    enum freesasa_algorithm:
        FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY
    enum freesasa_class:
        FREESASA_POLAR=0, FREESASA_APOLAR,
        FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN

    freesasa_t* freesasa_init()
 
    cdef int FREESASA_NAME_LIMIT
    cdef int FREESASA_DEF_PROBE_RADIUS
    cdef int FREESASA_DEF_SR_N
    cdef double FREESAAS_DEF_LR_D
    cdef int FREESASA_SUCCESS
    cdef int FREESASA_FAIL
    cdef int FREESASA_WARN

    void freesasa_free(freesasa_t *s)
    void freesasa_copy_param(freesasa_t *target, const freesasa_t *source)
    int freesasa_calc_coord(freesasa_t  *s, const double *coord,
                            const double *r, size_t n)
    int freesasa_calc_pdb(freesasa_t *s, FILE *pdb_file)
    int freesasa_calc_atoms(freesasa_t *s, const double *coord,
                            const char **resnames, 
                            const char **atomnames, size_t n)
    double freesasa_radius(const char *residue_name, const char *atom_name)
    int freesasa_link_coord(freesasa_t *s, const double *coord,
                            double *r, size_t n)
    int freesasa_refresh(freesasa_t *s)
    size_t freesasa_n_atoms(const freesasa_t *s)
    const char* freesasa_algorithm_name(const freesasa_t *s)
    int freesasa_set_probe_radius(freesasa_t *s,double r)
    double freesasa_get_probe_radius(const freesasa_t *s)
    int freesasa_set_sr_points(freesasa_t *s, int n)
    int freesasa_get_sr_points(const freesasa_t *s)
    int freesasa_set_lr_delta(freesasa_t *s, double d)
    double freesasa_get_lr_delta(const freesasa_t *s)
    int freesasa_set_nthreads(freesasa_t *s,int n)
    int freesasa_get_nthreads(const freesasa_t *s)
    void freesasa_set_proteinname(freesasa_t *s,const char *name)
    const char* freesasa_get_proteinname(const freesasa_t *s)
    double freesasa_area_total(const freesasa_t *s)
    double freesasa_area_class(const freesasa_t *s, freesasa_class c)
    int freesasa_per_residue_type(FILE *output, const freesasa_t *s)
    int freesasa_per_residue(FILE *output, const freesasa_t *s)
    double freesasa_area_residue(const freesasa_t *s, const char *res_name)
    int freesasa_write_pdb(FILE *output, const freesasa_t *s);
    double freesasa_area_atom(const freesasa_t *s, int i)
    const double* freesasa_area_atom_array(const freesasa_t *s)
    double freesasa_radius_atom(const freesasa_t *s, int i)
    const double* freesasa_radius_atom_array(const freesasa_t *s)
    int freesasa_log(FILE *log, const freesasa_t *s)

