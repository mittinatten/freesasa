from libc.stdio cimport FILE

#cdef extern from "stdint.h" nogil:
#    ctypedef size_t uintptr_t:
#        pass


cdef extern from "freesasa.h":
    ctypedef enum freesasa_algorithm:
        FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY
    ctypedef enum freesasa_class:
        FREESASA_POLAR, FREESASA_APOLAR,
        FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
    ctypedef enum freesasa_verbosity:
        FREESASA_V_NORMAL, FREESASA_V_SILENT

    cdef int FREESASA_DEF_PROBE_RADIUS
    cdef int FREESASA_DEF_SR_N
    cdef double FREESASA_DEF_LR_D
    cdef int FREESASA_SUCCESS
    cdef int FREESASA_FAIL
    cdef int FREESASA_WARN

    ctypedef struct freesasa_parameters:
        freesasa_algorithm alg
        double probe_radius
        int shrake_rupley_n_points
        double lee_richards_delta
        int n_threads

    ctypedef struct freesasa_result:
        double total
        double *sasa
        int n_atoms

    ctypedef struct freesasa_strvp:
        double *value
        char **string
        int n

    ctypedef struct freesasa_classifier:
        int n_classes
        void *config
        double (*radius)(const char* res_name,
                         const char* atom_name,
                         const freesasa_classifier *c)
        int (*sasa_class)(const char* res_name,
                          const char* atom_name,
                          const freesasa_classifier *c)
        const char* (*class2str)(int the_class,
                                 const freesasa_classifier *c)
        void (*free_config)(void*)

    ctypedef struct freesasa_structure:
        pass

    cdef extern const freesasa_parameters freesasa_default_parameters
    cdef extern const freesasa_classifier freesasa_default_classifier
    cdef extern const freesasa_classifier freesasa_residue_classifier

    freesasa_result* freesasa_calc_structure(const freesasa_structure *structure,
                                             const double *radii,
                                             const freesasa_parameters *parameters)

    freesasa_result* freesasa_calc_coord(const double *xyz,
                                         const double *radii,
                                         int n,
                                         const freesasa_parameters *parameters)

    void freesasa_result_free(freesasa_result *result)

    freesasa_structure* freesasa_structure_from_pdb(FILE *pdb,
                                                    int include_hetatm)

    int freesasa_structure_n(freesasa_structure *structure)

    void freesasa_structure_free(freesasa_structure* structure)

    double* freesasa_structure_radius(freesasa_structure *structure,
                                      freesasa_classifier *classifier)

    freesasa_classifier* freesasa_classifier_from_file(FILE *file)

    void freesasa_classifier_free(freesasa_classifier *classifier)

    freesasa_strvp* freesasa_result_classify(freesasa_result *result,
                                             const freesasa_structure *structure,
                                             const freesasa_classifier *classifier)

    void freesasa_strvp_free(freesasa_strvp *strvp)

    int freesasa_write_pdb(FILE *output,
                           freesasa_result *result,
                           const freesasa_structure *structure,
                           const double *radii)

    int freesasa_per_residue_type(FILE *output,
                                  freesasa_result *result,
                                  const freesasa_structure *structure)

    int freesasa_per_residue(FILE *output,
                             freesasa_result *result,
                             const freesasa_structure *structure)

    int freesasa_log(FILE *log,
                     freesasa_result *result,
                     const char *name,
                     const freesasa_parameters *parameters,
                     const freesasa_strvp* class_sasa)

    int freesasa_set_verbosity(freesasa_verbosity v)

    freesasa_verbosity freesasa_get_verbosity()

    freesasa_structure* freesasa_structure_new()

    int freesasa_structure_add_atom(freesasa_structure *s,
                                    const char* atom_name,
                                    const char* residue_name,
                                    const char* residue_number,
                                    char chain_label,
                                    double x, double y, double z)

    double* freesasa_structure_radius(const freesasa_structure *structure,
                                      const freesasa_classifier *classifier)

    const char* freesasa_structure_atom_name(const freesasa_structure *s,
                                             int i)

    const char* freesasa_structure_atom_res_name(const freesasa_structure *s,
                                                 int i)

    const char* freesasa_structure_atom_res_number(const freesasa_structure *s,
                                                   int i)

    char freesasa_structure_atom_chain(const freesasa_structure *s, int i)
