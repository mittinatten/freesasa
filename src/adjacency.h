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

#ifndef FREESASA_ADJACENCY_H
#define FREESASA_ADJACENCY_H

#include <stdlib.h>
#include "coord.h"

typedef struct freesasa_adjacency_element freesasa_adjacency_element;
//! linked list
struct freesasa_adjacency_element {
    int index; //! coordinate index
    double xd; //! x-distance to root 
    double yd; //! y-distance to root
    double xyd; //! distance to root in xy-plane
    freesasa_adjacency_element *next; //! next element, NULL if end of list
};

//! Adjacency list
typedef struct {
    freesasa_adjacency_element **list;
    freesasa_adjacency_element **tail;
    int n; //! number of elements
    int *nn; //! length of each list
} freesasa_adjacency;

/**
    Creates an adjacency list based on a set of coordinates with
    corresponding sphere radii. 

    Implemented using Verlet lists, giving O(N) performance. Should be
    freed with freesasa_adjacency_free(). For efficient calculations
    using this list the members of the returned struct should be used
    directly and not freesasa_adjacency_contact().

    @param coord a set of coordinates
    @param radii radii for the coordinates
    @return an adjacency list.
 */
freesasa_adjacency *freesasa_adjacency_new(const freesasa_coord *coord,
                                           const double *radii);

/**
    Frees an adjacency list created by freesasa_adjacency_new().

    @param adj the adjacency list to free
 */
void freesasa_adjacency_free(freesasa_adjacency *adj);

/**
    Checks if two atoms are in contact. Only included for reference.

    @param adj the adjacency list
    @param i index of first coordinate
    @param j index of second coordinate
    @return 1 if contact, 0 else.
 */
int freesasa_adjacency_contact(const freesasa_adjacency *adj,
                               int i, int j);

#endif
