#include "freesasa.h"
#include "freesasa_internal.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

struct residue_sasa {
    const char *name;
    double total;
    double main_chain;
    double side_chain;
    double polar;
    double apolar;
};

static const int n_ref = 20;
static struct residue_sasa sasa_ref[20] = {
    {.name = "ALA", .total = 103.10, .main_chain = 46.51, .side_chain = 56.60, .polar = 29.89, .apolar = 73.21},
    {.name = "CYS", .total = 125.02, .main_chain = 45.47, .side_chain = 79.55, .polar = 79.68, .apolar = 45.33},
    {.name = "ASP", .total = 135.76, .main_chain = 44.65, .side_chain = 91.11, .polar = 88.93, .apolar = 46.83},
    {.name = "GLU", .total = 167.95, .main_chain = 45.12, .side_chain = 122.83, .polar = 113.74, .apolar = 54.21},
    {.name = "PHE", .total = 193.68, .main_chain = 43.52, .side_chain = 150.16, .polar = 29.89, .apolar = 163.79},
    {.name = "GLY", .total = 71.84, .main_chain = 71.84, .side_chain = 0.00, .polar = 31.58, .apolar = 40.26},
    {.name = "HIS", .total = 173.43, .main_chain = 44.26, .side_chain = 129.18, .polar = 69.25, .apolar = 104.18},
    {.name = "ILE", .total = 167.30, .main_chain = 39.09, .side_chain = 128.22, .polar = 24.70, .apolar = 142.60},
    {.name = "LYS", .total = 197.47, .main_chain = 45.10, .side_chain = 152.38, .polar = 87.44, .apolar = 110.04},
    {.name = "LEU", .total = 160.87, .main_chain = 44.85, .side_chain = 116.01, .polar = 29.89, .apolar = 130.98},
    {.name = "MET", .total = 185.43, .main_chain = 45.08, .side_chain = 140.35, .polar = 67.61, .apolar = 117.83},
    {.name = "ASN", .total = 138.45, .main_chain = 43.80, .side_chain = 94.65, .polar = 93.13, .apolar = 45.32},
    {.name = "PRO", .total = 132.32, .main_chain = 29.83, .side_chain = 102.49, .polar = 16.16, .apolar = 116.16},
    {.name = "GLN", .total = 172.69, .main_chain = 45.09, .side_chain = 127.60, .polar = 123.13, .apolar = 49.56},
    {.name = "ARG", .total = 231.99, .main_chain = 45.09, .side_chain = 186.90, .polar = 153.92, .apolar = 78.07},
    {.name = "SER", .total = 111.49, .main_chain = 46.10, .side_chain = 65.39, .polar = 58.63, .apolar = 52.86},
    {.name = "THR", .total = 133.09, .main_chain = 40.38, .side_chain = 92.71, .polar = 49.91, .apolar = 83.18},
    {.name = "VAL", .total = 146.72, .main_chain = 44.24, .side_chain = 102.48, .polar = 29.89, .apolar = 116.83},
    {.name = "TRP", .total = 226.55, .main_chain = 40.50, .side_chain = 186.05, .polar = 61.19, .apolar = 165.37},
    {.name = "TYR", .total = 208.08, .main_chain = 43.49, .side_chain = 164.58, .polar = 76.46, .apolar = 131.62},
};

int
freesasa_print_rsa(FILE* output,
                   const freesasa_result *result,
                   const freesasa_structure *structure,
                   const char *name)
{
    assert(output);
    assert(result);
    assert(structure);

    int naa = freesasa_structure_n_residues(structure), first, last, resi, i_ref;
    double total, main_chain, side_chain, polar, apolar;
    char namebuf[FREESASA_MAX_SELECTION_NAME], selbuf[100];
    const char *resi_str, *resn;

    fprintf(output, "REM  Using default relative accessibilites\n");
    fprintf(output, "REM  Absolute and relative SASAs for %s\n", name);
    fprintf(output, "REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar\n");
    fprintf(output, "REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n");

    for (int i = 0; i < naa; ++i) {
        freesasa_structure_residue_atoms(structure, i, &first, &last);
        resi_str = freesasa_structure_atom_res_number(structure, first);
        resn = freesasa_structure_atom_res_name(structure, first);
        i_ref = -1;
        for (int j = 0; j < n_ref; ++j) {
            if (strcmp(sasa_ref[j].name,resn) == 0) { 
                i_ref = j;
                break; 
            }
        }
        resi = atoi(resi_str);
        fprintf(output, "RES %s  %s  ", resn, resi_str);
        total = freesasa_single_residue_sasa(result, structure, i);
        sprintf(selbuf,"bb, resi %d and name c+o+n+ca", resi);
        freesasa_select_area(selbuf, namebuf, &main_chain, structure, result);
        sprintf(selbuf,"sc, resi %d and not name c+o+n+ca", resi);
        freesasa_select_area(selbuf, namebuf, &side_chain, structure, result);
        sprintf(selbuf,"polar, resi %d and not symbol c", resi);
        freesasa_select_area(selbuf, namebuf, &polar, structure, result);
        sprintf(selbuf,"apolar, resi %d and symbol c", resi);
        freesasa_select_area(selbuf, namebuf, &apolar, structure, result);
        if (i_ref >= 0) {
            fprintf(output, "%7.2f%6.1f%7.2f",
                    total, 100*total/sasa_ref[i_ref].total, side_chain);
            if (sasa_ref[i_ref].side_chain > 0)
                fprintf(output, "%6.1f", 100*side_chain/sasa_ref[i_ref].side_chain);
            else fprintf(output, "   N/A");
            fprintf(output, "%7.2f%6.1f%7.2f%6.1f%7.2f%6.1f\n",
                    main_chain, 100*main_chain/sasa_ref[i_ref].main_chain,
                    apolar, 100*apolar/sasa_ref[i_ref].apolar,
                    polar, 100*polar/sasa_ref[i_ref].polar);
        } else {
            fprintf(output,
                    "%7.2f   N/A%7.2f   N/A%7.2f   N/A%7.2f   N/A%7.2f   N/A\n",
                    total, side_chain, main_chain, apolar, polar);
        }
    }
    
    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    
    return FREESASA_SUCCESS;
}
