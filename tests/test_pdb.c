/*
  Copyright Simon Mitternacht 2013-2014.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <pdb.h>
#include <check.h>

START_TEST (test_pdb_empty_lines)
{
    char buf[80];
    double x[3];

    // check string parsing
    ck_assert_int_eq(sasalib_pdb_get_atom_name(buf,""),SASALIB_FAIL);
    ck_assert_str_eq(buf,"");
    ck_assert_int_eq(sasalib_pdb_get_res_name(buf,""),SASALIB_FAIL);
    ck_assert_str_eq(buf,"");
    ck_assert_int_eq(sasalib_pdb_get_res_number(buf,""),SASALIB_FAIL);
    ck_assert_str_eq(buf,"");

    // check coordinate parsing
    ck_assert_int_eq(sasalib_pdb_get_coord(x,""),SASALIB_FAIL);

    // check label parsing
    ck_assert_int_eq(sasalib_pdb_get_chain_label(""),'\0');
    ck_assert_int_eq(sasalib_pdb_get_alt_coord_label(""),'\0');

    // check element parsing
    ck_assert_int_eq(sasalib_pdb_ishydrogen(""),SASALIB_FAIL);
}
END_TEST

START_TEST (test_pdb_lines)
{
    char buf[80];
    double x[3];
    const char lines[][80] = { 
	"ATOM    585  C   ARG A  74      41.765  34.829  30.944  0.45 36.22           C",
	"ATOM    573  NH1AARG A  72      34.110  28.437  27.768  1.00 35.02           N",
	"HETATM  610  O   HOH A  83      27.707  15.908   4.653  1.00 20.30           O",
	"ATOM    573  H   ARG A  72      34.110  28.437  27.768  1.00 35.02           H",};
    //Atom-name
    sasalib_pdb_get_atom_name(buf,lines[0]);
    ck_assert_str_eq(buf," C  ");
    sasalib_pdb_get_atom_name(buf,lines[1]);
    ck_assert_str_eq(buf," NH1");
    sasalib_pdb_get_atom_name(buf,lines[2]);
    ck_assert_str_eq(buf," O  ");

    //Res-name
    sasalib_pdb_get_res_name(buf,lines[0]);
    ck_assert_str_eq(buf,"ARG");
    sasalib_pdb_get_res_name(buf,lines[2]);
    ck_assert_str_eq(buf,"HOH");

    //Res-number
    sasalib_pdb_get_res_number(buf,lines[0]);
    ck_assert_str_eq(buf,"  74");
    sasalib_pdb_get_res_number(buf,lines[1]);
    ck_assert_str_eq(buf,"  72");
    sasalib_pdb_get_res_number(buf,lines[2]);
    ck_assert_str_eq(buf,"  83");
    
    //coordinates
    sasalib_pdb_get_coord(x,lines[0]);
    ck_assert_int_eq((fabs(x[0] - 41.765) > 1e-6 || fabs(x[1] - 34.829) > 1e-6 ||
		      fabs(x[2] - 30.944) > 1e-6),0);
    
    //chain label
    ck_assert_int_eq(sasalib_pdb_get_chain_label(lines[0]), 'A');
    
    // alt coord labels
    ck_assert_int_eq(sasalib_pdb_get_alt_coord_label(lines[0]), ' ');
    ck_assert_int_eq(sasalib_pdb_get_alt_coord_label(lines[1]), 'A');

    // is hydrogen 
    ck_assert(!sasalib_pdb_ishydrogen(lines[0]));
    ck_assert(sasalib_pdb_ishydrogen(lines[3]));

}
END_TEST

Suite *pdb_suite() {
    Suite *s = suite_create("PDB-parser");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_pdb_empty_lines);
    tcase_add_test(tc_core, test_pdb_lines);
    suite_add_tcase(s, tc_core);

    return s;
}
