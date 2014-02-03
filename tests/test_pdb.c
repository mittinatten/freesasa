#include <string.h>
#include <math.h>
#include <stdio.h>
#include <pdb.h>

int test_pdb_empty_lines()
{
    char buf[80];
    double x[3];
    int n_err = 0;

    printf("Testing error reporting for empty PDB-lines.\n");
    if (sasalib_pdb_get_atom_name(buf,"") != SASALIB_FAIL) {
	++n_err;
	printf("Error: Extracting atom-name from empty string did not generate error.\n");
    }
    if (strlen(buf) > 0) {
	++n_err;
	printf("Error: Extracted non-empty atom-name from empty string.\n");
    }
    if (sasalib_pdb_get_res_name(buf,"") != SASALIB_FAIL) {
	++n_err;
	printf("Error: Extracting res-name from empty string did not generate error.\n");
    }
    if (strlen(buf) > 0) {
	++n_err;
	printf("Error: Extracted non-empty res-name from empty string.\n");
    }
    if (sasalib_pdb_get_res_number(buf,"") != SASALIB_FAIL) {
	++n_err;
	printf("Error: Extracting res-number from empty string did not generate error.\n");
    } 
    if (strlen(buf) > 0) {
	++n_err;
	printf("Error: Extracted non-empty res-number from empty string.\n");
    }
    if (sasalib_pdb_get_coord(x,"") != SASALIB_FAIL) {
	++n_err;
	printf("Error: Extracting coordinates from empty string did not generate error.\n");
    }
    if (sasalib_pdb_get_chain_label("") != '\0') {
	++n_err;
	printf("Error: Extracting chain label from empty string did not generate error.\n");
    }
    if (sasalib_pdb_get_alt_coord_label("") != '\0') {
	++n_err;
	printf("Error: Extracting alt coord label from empty string did not generate error.\n");
    }
    if (sasalib_pdb_ishydrogen("") != SASALIB_FAIL) {
	++n_err;
	printf("Error: Checking for hydrogen in empty string did not generate error.\n");
    }
    return n_err;
}

int test_pdb_string_err(const char *result, const char* desired, 
		 const char *input, const char *descriptor)
{
    if (strcmp(result,desired) != 0) {
	printf("Error: failed extracting %s from string:\n  >>%s<<\n"
	       "       Result was '%s', expected value is '%s'.\n",
	       descriptor,input,result,desired);
	return 1;
    }
    return 0;
}

int test_pdb_lines() 
{
    int n_err = 0;
    char buf[80];
    double x[3];
    printf("Testing parsing of single PDB lines.\n");
    const char lines[][80] = { 
	"ATOM    585  C   ARG A  74      41.765  34.829  30.944  0.45 36.22           C",
	"ATOM    573  NH1AARG A  72      34.110  28.437  27.768  1.00 35.02           N",
	"HETATM  610  O   HOH A  83      27.707  15.908   4.653  1.00 20.30           O",
	"ATOM    573  H   ARG A  72      34.110  28.437  27.768  1.00 35.02           H",};
    //Atom-name
    sasalib_pdb_get_atom_name(buf,lines[0]);
    n_err += test_pdb_string_err(buf," C  ",lines[0],"atom-name");
    sasalib_pdb_get_atom_name(buf,lines[1]);
    n_err += test_pdb_string_err(buf," NH1",lines[1],"atom-name");
    sasalib_pdb_get_atom_name(buf,lines[2]);
    n_err += test_pdb_string_err(buf," O  ",lines[2],"atom-name");

    //Res-name
    sasalib_pdb_get_res_name(buf,lines[0]);
    n_err += test_pdb_string_err(buf,"ARG",lines[0],"res-name");
    sasalib_pdb_get_res_name(buf,lines[2]);
    n_err += test_pdb_string_err(buf,"HOH",lines[2],"res-name");

    //Res-number
    sasalib_pdb_get_res_number(buf,lines[0]);
    n_err += test_pdb_string_err(buf,"  74",lines[0],"res-number");
    sasalib_pdb_get_res_number(buf,lines[1]);
    n_err += test_pdb_string_err(buf,"  72",lines[1],"res-number");
    sasalib_pdb_get_res_number(buf,lines[2]);
    n_err += test_pdb_string_err(buf,"  83",lines[2],"res-number");

    //coordinates
    sasalib_pdb_get_coord(x,lines[0]);
    if (fabs(x[0] - 41.765) > 1e-6 || fabs(x[1] - 34.829) > 1e-6 ||
	fabs(x[2] - 30.944) > 1e-6) {
	printf("Error: failed to extract coordinates from line:\n  >>%s<<\n"
	       "       Got (%f,%f,%f)\n",lines[0],x[0],x[1],x[2]);
	++n_err;
    }
    
    //chain label
    char c;
    if ((c = sasalib_pdb_get_chain_label(lines[0])) != 'A') {
	printf("Error: failed to extract chain label from line:\n  >>%s<<\n"
	       "       Got '%c'\n",lines[0],c);
	++n_err;
    }

    // alt coord labels
    if ((c = sasalib_pdb_get_alt_coord_label(lines[0])) != ' ') {
	printf("Error: failed to extract alt coord label from line:\n  >>%s<<\n"
	       "       Got '%c'\n",lines[0],c);
	++n_err;
    }
    if ((c = sasalib_pdb_get_alt_coord_label(lines[1])) != 'A') {
	printf("Error: failed to extract alt coord label from line:\n  >>%s<<\n"
	       "       Got '%c'\n",lines[1],c);
	++n_err;
    }

    // is hydrogen 
    if (sasalib_pdb_ishydrogen(lines[0])) {
	printf("Error: the following was reported as a hydrogen:\n >>%s<<\n",lines[0]);
	++n_err;
    }
    if (!sasalib_pdb_ishydrogen(lines[3])) {
	printf("Error: the following was reported as not hydrogen:\n >>%s<<\n",lines[3]);
	++n_err;
    }

    return n_err;
    
}

int test_pdb() {
    int n_err = 0;
    n_err += test_pdb_empty_lines();
    n_err += test_pdb_lines();

    return n_err;
}
