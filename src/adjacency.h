/*
  Copyright Simon Mitternacht 2013-2015.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FREESASA_CELL_H
#define FREESASA_CELL_H

#include <stdlib.h>
#include "coord.h"

#ifndef FREESASA_ATOMS_PER_CELL
#define FREESASA_ATOMS_PER_CELL 60
#endif

//! Adjacency list
typedef struct {
    int **nb; //! neighbors to each element
    size_t *nn; //! number of neighbors
    size_t n; //! number of elements
    double **nb_xyd; //! distance between neighbors in xy-plane
    double **nb_xd; //! signed distance between neighbors along x-axis
    double **nb_yd; //! signed distance between neighbors along y-axis
} freesasa_adjacency;

typedef struct freesasa_cell_list freesasa_cell_list;

freesasa_cell_list* freesasa_cell_list_new(double cell_size,
                                           const freesasa_coord *coord);

void freesasa_cell_list_free(freesasa_cell_list *c);

freesasa_adjacency *freesasa_adjacency_new(const freesasa_coord *coord,
                                           const double *radii);

void freesasa_adjacency_free(freesasa_adjacency *adj);

int freesasa_adjacency_contact(const freesasa_adjacency *adj,
                               int i, int j);

#endif
