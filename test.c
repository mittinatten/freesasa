/*
  Copyright Simon Mitternacht 2013.

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

#include <stdio.h>
#include <stdlib.h>
#include "src/protein.h"
#include "src/sasa.h"

/* For now this program just outputs some raw data for sanity
   check. */

int main (int argc, char **argv) {
    protein *p = protein_init_from_pdb(stdin);
    double *sasa_sr = (double*) malloc(sizeof(double)*protein_n(p));
    double *sasa_lr = (double*) malloc(sizeof(double)*protein_n(p));
    double *r = (double*) malloc(sizeof(double)*protein_n(p));
    protein_r_def(p,r);
    sasa_shrake_rupley(sasa_sr,protein_xyz(p),r,protein_n(p),5000);
    sasa_lee_richards(sasa_lr,protein_xyz(p),r,protein_n(p),0.01);
    printf ("%ul\n",protein_n(p));
    for (int i = 0; i < protein_n(p); ++i) {
	printf("%5i\t%4s\t%3s\t%4s\t%c\t%6.2f\t%6.2f\t\n",
	       i,
	       protein_atom_name(p,i),
	       protein_atom_res_name(p,i),
	       protein_atom_res_number(p,i),
	       protein_atom_chain(p,i),
	       //atom_type2str(a->at),
	       //atom_class2str(a->ac),
	       sasa_sr[i], sasa_lr[i]
	    );
	
    }
    protein_free(p);
    free(sasa_sr);
    free(sasa_lr);
    free(r);
}
