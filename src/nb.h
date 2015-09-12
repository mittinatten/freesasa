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

#ifndef FREESASA_NB_H
#define FREESASA_NB_H

#include <stdlib.h>
#include "coord.h"
/**
   @file
   @author Simon Mitternacht
   
   Functions to compute neigbor lists. The function
   freesasa_nb_contact() is mainly intended for checking
   consitency, in performance-critical code it is advisible to use the
   struct (as demonstrated by the
 */

//! neigbor list
typedef struct {
    int n; //!< number of elements
    int **nb; //!< neighbors to each element
    int *nn; //!< number of neighbors to each element
    double **nb_xyd; //!< distance between neighbors in xy-plane
    double **nb_xd; //!< signed distance between neighbors along x-axis
    double **nb_yd; //!< signed distance between neighbors along y-axis
    int *capacity; //!< keeps track of memory chunks (don't change this)
} freesasa_nb;

/**
    Creates a neigbor list based on a set of coordinates with
    corresponding sphere radii. 

    Implemented using neigbor lists, giving O(N) performance. Should be
    freed with freesasa_nb_free(). For efficient calculations
    using this list the members of the returned struct should be used
    directly and not freesasa_nb_contact().

    @param coord a set of coordinates
    @param radii radii for the coordinates
    @return a neigbor list.
 */
freesasa_nb *freesasa_nb_new(const freesasa_coord *coord,
                             const double *radii);

/**
    Frees a neigbor list created by freesasa_nb_new().

    @param adj the neigbor list to free
 */
void freesasa_nb_free(freesasa_nb *adj);

/**
    Checks if two atoms are in contact. Only included for reference.

    @param adj the neigbor list
    @param i index of first coordinate
    @param j index of second coordinate
    @return 1 if contact, 0 else.
 */
int freesasa_nb_contact(const freesasa_nb *adj,
                        int i, int j);

#endif /* FREESASA_NB_H*/
