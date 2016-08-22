#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <freesasa.h>
#include <freesasa_internal.h>
#include <pdb.h>
#include "tools.h"

#define N 6
const char an[N][PDB_ATOM_NAME_STRL+1] =  {" C  "," CA "," O  "," CB "," SD ", "SE  "};
const char rna[N][PDB_ATOM_RES_NAME_STRL+1] = {"MET", "MET", "MET", "MET", "MET", "SEC",};
const char rnu[N][PDB_ATOM_RES_NUMBER_STRL+1] = {"   1","   1","   1","   1","   1","   2"};
const char symbol[N][PDB_ATOM_SYMBOL_STRL+1] = {" C"," C"," O"," C"," S","SE"};
const char cl[N] = {'A','A','A','A','A','A'};
const double bfactors[N] = {1., 1., 1., 1., 1., 1.};

freesasa_structure *s;

START_TEST (test_structure_api)
{
    s = freesasa_structure_new();
    for (int i = 0; i < N; ++i) {
        ck_assert_int_eq(freesasa_structure_add_atom(s,an[i],rna[i],rnu[i],cl[i],i,i,i),
                         FREESASA_SUCCESS);
    }
    for (int i = 0; i < N; ++i) {
        ck_assert_str_eq(freesasa_structure_atom_name(s,i),an[i]);
        ck_assert_str_eq(freesasa_structure_atom_res_name(s,i),rna[i]);
        ck_assert_str_eq(freesasa_structure_atom_res_number(s,i),rnu[i]);
        ck_assert_str_eq(freesasa_structure_atom_symbol(s,i),symbol[i]);
        ck_assert_int_eq(freesasa_structure_atom_chain(s,i),cl[i]);
    }
    freesasa_structure_atom_set_radius(s, 0, 10.0);
    ck_assert(float_eq(freesasa_structure_atom_radius(s, 0), 10.0, 1e-10));

    const coord_t *c = freesasa_structure_xyz(s);
    ck_assert(freesasa_structure_coord_array(s) == freesasa_coord_all(c));
    for (int i = 0; i < N; ++i) {
        const double *xyz = freesasa_coord_i(c, i);
        ck_assert(fabs(xyz[0]+xyz[1]+xyz[2]-3*i) < 1e-10);
    }
    ck_assert_int_eq(freesasa_structure_n(s), N);
    ck_assert_int_eq(freesasa_structure_n_residues(s), 2);
    ck_assert_int_eq(freesasa_structure_n_chains(s), 1);

    ck_assert_str_eq(freesasa_structure_residue_name(s,0), rna[0]);
    ck_assert_str_eq(freesasa_structure_residue_number(s,0), rnu[0]);
    ck_assert_int_eq(freesasa_structure_residue_chain(s,0), cl[0]);

    ck_assert_int_eq(freesasa_structure_chain_index(s, 'A'), 0);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_structure_chain_index(s, 'B'), FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
        
    int first, last;
    ck_assert(freesasa_structure_residue_atoms(s, 0, &first, &last) == FREESASA_SUCCESS);
    ck_assert(first == 0 && last == N-2);
    ck_assert(freesasa_structure_chain_atoms(s, 'A', &first, &last) == FREESASA_SUCCESS);
    ck_assert(first == 0 && last == N-1);
    ck_assert(freesasa_structure_chain_residues(s, 'A', &first, &last) == FREESASA_SUCCESS);
    ck_assert_int_eq(first, 0);
    ck_assert_int_eq(last, 1);
    freesasa_structure_free(s);
    s = NULL;
}
END_TEST

START_TEST (test_add_atom)
{
    freesasa_structure *s = freesasa_structure_new();
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_structure_add_atom(s,"HABC","ALA","   1",'A',0,0,0), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom(s,"SE  ","SEC","   1",'A',0,0,0), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom(s,"CL  ","ABC","   1",'A',0,0,0), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_n(s), 3);
    ck_assert_str_eq(freesasa_structure_atom_symbol(s,0)," H");
    ck_assert_str_eq(freesasa_structure_atom_symbol(s,1),"SE");
    ck_assert_str_eq(freesasa_structure_atom_symbol(s,2),"CL");

    ck_assert_int_eq(freesasa_structure_add_atom(s,"A","ALA","   1",'A',0,0,0),FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom(s," C  ","AL","   1",'A',0,0,0),FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom(s," C  ","ALA"," 1",'A',0,0,0),FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_n(s), 6);

    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"A","ALA","   1",'A',0,0,0,NULL,FREESASA_SKIP_UNKNOWN),FREESASA_WARN);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"HABC","ALA","   1",'A',0,0,0,NULL,FREESASA_SKIP_UNKNOWN), FREESASA_WARN);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"SE  ","SEC","   1",'A',0,0,0,NULL,FREESASA_SKIP_UNKNOWN), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"CL  ","ABC","   1",'A',0,0,0,NULL,FREESASA_SKIP_UNKNOWN), FREESASA_WARN);
    ck_assert_int_eq(freesasa_structure_n(s), 7);

    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"HABC","ALA","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN), FREESASA_FAIL);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"SE  ","SEC","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"CL  ","ABC","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN), FREESASA_FAIL);
    ck_assert_int_eq(freesasa_structure_n(s), 8);

    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"HABC","ALA","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN | FREESASA_SKIP_UNKNOWN), FREESASA_FAIL);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"SE  ","SEC","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN | FREESASA_SKIP_UNKNOWN), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_structure_add_atom_wopt(s,"CL  ","ABC","   1",'A',0,0,0,NULL,FREESASA_HALT_AT_UNKNOWN | FREESASA_SKIP_UNKNOWN), FREESASA_FAIL);
    ck_assert_int_eq(freesasa_structure_n(s), 9);

    freesasa_set_verbosity(FREESASA_V_NORMAL);

    freesasa_structure_free(s);
}
END_TEST

double a2r(const char *rn, const char *am)
{
    return 1.0;
}

void setup_1ubq(void)
{
    errno = 0;
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    if (s) freesasa_structure_free(s);
    s = freesasa_structure_from_pdb(pdb, NULL, 0);
    fclose(pdb);
}

void teardown_1ubq(void)
{
    freesasa_structure_free(s);
    s = NULL;
}

START_TEST (test_structure_1ubq)
{
    ck_assert(freesasa_structure_n(s) == 602);
    ck_assert(freesasa_structure_n_residues(s) == 76);
    // check at random atom to see that parsing was correct
    ck_assert_str_eq(freesasa_structure_atom_res_name(s,8),"GLN");
    ck_assert_str_eq(freesasa_structure_atom_name(s,8)," N  ");
    ck_assert_str_eq(freesasa_structure_atom_res_number(s,8),"   2");
    ck_assert_int_eq(freesasa_structure_atom_chain(s,8),'A');
    ck_assert_str_eq(freesasa_structure_atom_symbol(s,8)," N");

    // check coordinates of that random atom
    const coord_t *c = freesasa_structure_xyz(s);
    ck_assert(c != NULL);
    const double *x = freesasa_coord_i(c,8);
    ck_assert(x != NULL);
    ck_assert(fabs(x[0]-26.335+x[1]-27.770+x[2]-3.258) < 1e-10);
}
END_TEST

START_TEST (test_pdb)
{
    const char *file_names[] = {DATADIR "alt_model_twochain.pdb",
                                DATADIR "empty.pdb",
                                DATADIR "empty_model.pdb"};
    const int result_null[] = {0,1,1};
    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 0; i < 3; ++i) {
        FILE *pdb = fopen(file_names[i],"r");
        ck_assert(pdb != NULL);
        freesasa_structure *s = freesasa_structure_from_pdb(pdb,NULL,0);
        if (result_null[i]) ck_assert(s == NULL);
        else ck_assert(s != NULL);
        fclose(pdb);
        freesasa_structure_free(s);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_hydrogen) 
{
    FILE *pdb = fopen(DATADIR "1d3z.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_structure* s = freesasa_structure_from_pdb(pdb,NULL,0);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 602);
    freesasa_structure_free(s);
    rewind(pdb);

    freesasa_set_verbosity(FREESASA_V_SILENT);

    s = freesasa_structure_from_pdb(pdb,NULL,FREESASA_INCLUDE_HYDROGEN);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 1231);
    freesasa_structure_free(s);

    rewind(pdb);
    ck_assert_ptr_eq(freesasa_structure_from_pdb(pdb,NULL,FREESASA_INCLUDE_HYDROGEN | FREESASA_HALT_AT_UNKNOWN), NULL);

    rewind(pdb);
    s = freesasa_structure_from_pdb(pdb,NULL,FREESASA_INCLUDE_HYDROGEN | FREESASA_SKIP_UNKNOWN);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 602);
    freesasa_structure_free(s);

    freesasa_set_verbosity(FREESASA_V_NORMAL);
    fclose(pdb);
}
END_TEST

START_TEST (test_hetatm) 
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_set_verbosity(FREESASA_V_SILENT); // unknown atoms warnings suppressed
    freesasa_structure* s = freesasa_structure_from_pdb(pdb,NULL,FREESASA_INCLUDE_HETATM);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 660);
    freesasa_structure_free(s);
    fclose(pdb);
}
END_TEST

START_TEST (test_structure_array)
{
    FILE *pdb;
    int n = 0;
    freesasa_structure **ss;
    
    freesasa_set_verbosity(FREESASA_V_SILENT);
    pdb = fopen(DATADIR "err.config", "r");
    ck_assert_ptr_eq(freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS), NULL);
    fclose(pdb);

    pdb = fopen(DATADIR "1d3z.pdb", "r");
    ck_assert_ptr_eq(freesasa_structure_array(pdb, &n, NULL, 0), NULL);

    rewind(pdb);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS);

    ck_assert(ss != NULL);
    ck_assert(n == 10);

    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 602);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    rewind(pdb);
    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_CHAINS);
    ck_assert(ss != NULL);
    ck_assert(n == 1);
    ck_assert(ss[0] != NULL);
    ck_assert(freesasa_structure_n(ss[0]) == 602);
    freesasa_structure_free(ss[0]);
    free(ss);
    
    rewind(pdb);
    freesasa_set_verbosity(FREESASA_V_SILENT); // Silence Hydrogen warnings
    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_CHAINS | FREESASA_SEPARATE_MODELS |
                                  FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 1231);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    fclose(pdb);
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS |
                                  FREESASA_INCLUDE_HETATM | FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 286*4);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    fclose(pdb);
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS | FREESASA_SEPARATE_CHAINS |
                                  FREESASA_INCLUDE_HETATM | FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10*4);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 286);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    rewind(pdb);
    freesasa_structure *s = freesasa_structure_from_pdb(pdb, NULL, FREESASA_INCLUDE_HETATM | 
                                                        FREESASA_INCLUDE_HYDROGEN |
                                                        FREESASA_JOIN_MODELS);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 286*4*10);
    freesasa_structure_free(s);
    
    fclose(pdb);
    pdb = fopen(DATADIR "1ubq.pdb","r");
    ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_CHAINS);
    
    ck_assert(ss != NULL);
    ck_assert(n == 1);
    ck_assert(freesasa_structure_n(ss[0]) == 602);
    freesasa_structure_free(ss[0]);
    free(ss);
    fclose(pdb);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    pdb = fopen(DATADIR "err.config","r");
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ss = freesasa_structure_array(pdb, &n, NULL,  FREESASA_SEPARATE_CHAINS);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    ck_assert(ss == NULL);
}
END_TEST

START_TEST (test_get_chains) {
    FILE *pdb = fopen(DATADIR "2jo4.pdb","r");
    freesasa_structure *s = freesasa_structure_from_pdb(pdb, NULL, 0);
    int first, last;
    ck_assert_int_eq(freesasa_structure_n(s), 4*129);
    ck_assert_int_eq(freesasa_structure_n_chains(s), 4);
    ck_assert_str_eq(freesasa_structure_chain_labels(s), "ABCD");
    ck_assert_int_eq(freesasa_structure_chain_atoms(s, 'A', &first, &last),FREESASA_SUCCESS);
    ck_assert_int_eq(first, 0);
    ck_assert_int_eq(last, 128);
    ck_assert_int_eq(freesasa_structure_chain_atoms(s, 'B', &first, &last),FREESASA_SUCCESS);
    ck_assert_int_eq(first, 129);
    ck_assert_int_eq(last, 129*2-1);
    ck_assert_int_eq(freesasa_structure_chain_atoms(s, 'C', &first, &last),FREESASA_SUCCESS);
    ck_assert_int_eq(first, 129*2);
    ck_assert_int_eq(last, 129*3-1);
    ck_assert_int_eq(freesasa_structure_chain_atoms(s, 'D', &first, &last),FREESASA_SUCCESS);
    ck_assert_int_eq(first, 129*3);
    ck_assert_int_eq(last, 129*4-1);
    
    freesasa_structure *s2 = freesasa_structure_get_chains(s,"");
    ck_assert(s2 == NULL);
    s2 = freesasa_structure_get_chains(s,"X");
    ck_assert(s2 == NULL);

    s2 = freesasa_structure_get_chains(s,"A");
    ck_assert(freesasa_structure_n(s2) == 129);
    ck_assert(freesasa_structure_atom_chain(s2,0) == 'A');
    ck_assert_str_eq(freesasa_structure_chain_labels(s2),"A");
    freesasa_structure_free(s2);

    s2 = freesasa_structure_get_chains(s,"D");
    ck_assert(freesasa_structure_n(s2) == 129);
    ck_assert(freesasa_structure_atom_chain(s2,0) == 'D');
    ck_assert_str_eq(freesasa_structure_chain_labels(s2),"D");
    freesasa_structure_free(s2);

    s2 = freesasa_structure_get_chains(s,"AC");
    ck_assert(freesasa_structure_n(s2) == 2*129);
    ck_assert(freesasa_structure_atom_chain(s2,0) == 'A');
    ck_assert(freesasa_structure_atom_chain(s2,129) == 'C');
    ck_assert_str_eq(freesasa_structure_chain_labels(s2),"AC");
    freesasa_structure_free(s2);
    freesasa_structure_free(s);
}
END_TEST

START_TEST (test_occupancy)
{
    FILE *pdb = fopen(DATADIR "1ubq.occ.pdb", "r");
    ck_assert_ptr_ne(pdb, NULL);
    freesasa_structure *s = freesasa_structure_from_pdb(pdb, NULL, FREESASA_RADIUS_FROM_OCCUPANCY);
    fclose(pdb);
    ck_assert_ptr_ne(s, NULL);
    const double *r = freesasa_structure_radius(s);
    ck_assert_ptr_ne(r, NULL);
    ck_assert(float_eq(r[0], 3.64, 1e-6));
    ck_assert(float_eq(r[1], 1.88, 1e-6));
    ck_assert(float_eq(r[2], 1.61, 1e-6));
    freesasa_structure_free(s);
    
}
END_TEST


START_TEST (test_structure_node)
{
    s = freesasa_structure_new();
    for (int i = 0; i < N; ++i) {
        freesasa_structure_add_atom(s,an[i],rna[i],rnu[i],cl[i],i,i,i);
    }
    freesasa_result *result = freesasa_calc_structure(s, NULL);
    freesasa_structure_node *tree =
        freesasa_result2tree(result, s, NULL, "test");
    const freesasa_structure_node *next, *structure, *chain, *residue, *atom;
    
    const freesasa_nodearea *area;
    ck_assert_ptr_ne(tree, NULL);
    ck_assert_ptr_eq(freesasa_structure_node_parent(tree), NULL);
    ck_assert_int_eq(freesasa_structure_node_type(tree), FREESASA_NODE_ROOT);
    ck_assert_str_eq(freesasa_structure_node_name(tree), "test");
    ck_assert_ptr_ne((structure = freesasa_structure_node_children(tree)), NULL);
    ck_assert_ptr_ne((chain = freesasa_structure_node_children(structure)), NULL);
    ck_assert_ptr_ne((residue = freesasa_structure_node_children(chain)), NULL);
    ck_assert_ptr_ne((atom = freesasa_structure_node_children(residue)), NULL);

    ck_assert_int_eq(freesasa_structure_node_type(structure), FREESASA_NODE_STRUCTURE);
    ck_assert_int_eq(freesasa_structure_node_type(chain), FREESASA_NODE_CHAIN);
    ck_assert_int_eq(freesasa_structure_node_type(residue), FREESASA_NODE_RESIDUE);
    ck_assert_int_eq(freesasa_structure_node_type(atom), FREESASA_NODE_ATOM);

    ck_assert_str_eq(freesasa_structure_node_name(structure), "A");
    ck_assert_str_eq(freesasa_structure_node_name(chain), "A");
    ck_assert_str_eq(freesasa_structure_node_name(residue), "MET");
    ck_assert_str_eq(freesasa_structure_node_name(atom), " C  ");

    // iterate
    next = freesasa_structure_node_next(tree);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_structure_node_next(structure);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_structure_node_next(chain);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_structure_node_next(residue);
    ck_assert_ptr_ne(next, NULL);
    ck_assert_str_eq(freesasa_structure_node_name(next), "SEC");
    next = freesasa_structure_node_next(next);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_structure_node_next(atom);
    ck_assert_str_eq(freesasa_structure_node_name(next), " CA ");

    ck_assert_ptr_eq(freesasa_structure_node_area(tree), NULL);
    ck_assert_ptr_ne((area = freesasa_structure_node_area(structure)), NULL);
    ck_assert(float_eq(result->total, area->total, 1e-10));
    ck_assert_ptr_ne((area = freesasa_structure_node_area(chain)), NULL);
    ck_assert(float_eq(result->total, area->total, 1e-10));
    next = freesasa_structure_node_next(residue);
    ck_assert_ptr_ne((area = freesasa_structure_node_area(next)), NULL);
    ck_assert(float_eq(result->sasa[5], area->total, 1e-10));
    ck_assert_str_eq(area->name, "SEC");
    ck_assert_ptr_ne((area = freesasa_structure_node_area(atom)), NULL);
    ck_assert(float_eq(result->sasa[0], area->total, 1e-10));
    next = freesasa_structure_node_next(atom);
    ck_assert_ptr_ne((area = freesasa_structure_node_area(next)), NULL);
    ck_assert(float_eq(result->sasa[1], area->total, 1e-10));

    freesasa_structure_node_free(tree);
    freesasa_structure_free(s);
    freesasa_result_free(result);
}
END_TEST

START_TEST (test_memerr)
{
    FILE *file = fopen(DATADIR "1ubq.pdb","r");
    void *ptr;
    int n;
    freesasa_set_verbosity(FREESASA_V_SILENT);
    set_fail_after(1);
    ptr = freesasa_structure_new();
    set_fail_after(0);
    ck_assert_ptr_eq(ptr, NULL);
    for (int i = 1; i < 100; ++i) {
        rewind(file);
        set_fail_after(i);
        ptr = freesasa_structure_from_pdb(file, NULL, 0);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
    }
    freesasa_structure *structure = freesasa_structure_from_pdb(file, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(structure, NULL);
    for (int i = 1; i < 100; ++i) {
        set_fail_after(i);
        ptr = freesasa_result2tree(result, structure, NULL, "test");
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
    }
    freesasa_structure_free(structure);
    freesasa_result_free(result);
    fclose(file);


    file = fopen(DATADIR "2jo4.pdb", "r");
    for (int i = 1; i < 100; ++i) {
        rewind(file);
        set_fail_after(i);
        ptr = freesasa_structure_array(file, &n, NULL, FREESASA_SEPARATE_MODELS);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr,NULL);

        rewind(file);
        set_fail_after(i);
        ptr = freesasa_structure_array(file, &n, NULL, FREESASA_SEPARATE_MODELS | FREESASA_SEPARATE_CHAINS);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr,NULL);
    }
    set_fail_after(0);
    fclose(file);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

Suite* structure_suite() {
    // what goes in what Case is kind of arbitrary
    Suite *s = suite_create("Structure");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_structure_api);
    tcase_add_test(tc_core, test_add_atom);
    tcase_add_test(tc_core, test_memerr);

    TCase *tc_pdb = tcase_create("PDB");
    tcase_add_test(tc_pdb,test_pdb);
    tcase_add_test(tc_pdb,test_hydrogen);
    tcase_add_test(tc_pdb,test_hetatm);
    tcase_add_test(tc_pdb,test_structure_array);
    tcase_add_test(tc_pdb,test_get_chains);
    tcase_add_test(tc_pdb,test_occupancy);

    TCase *tc_1ubq = tcase_create("1UBQ");
    tcase_add_checked_fixture(tc_1ubq,setup_1ubq,teardown_1ubq);
    tcase_add_test(tc_1ubq,test_structure_1ubq);

    TCase *tc_node = tcase_create("Structure-node");
    tcase_add_test(tc_node, test_structure_node);

    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_pdb);
    suite_add_tcase(s, tc_1ubq);
    suite_add_tcase(s, tc_node);
    return s;
}

