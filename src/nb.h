#ifndef FREESASA_NB_H
#define FREESASA_NB_H

#include <stdlib.h>

#include "coord.h"

/**
   @file
   @author Simon Mitternacht

   Functions to compute neighbor lists. The function
   freesasa_nb_contact() is mainly intended for checking consistency,
   in performance-critical code it is advisible to use the struct (as
   demonstrated in sasa_lr.c and sasa_sr.c).
 */

/** Neighbor list */
typedef struct {
    int n;         /**< number of elements */
    int **nb;      /**< neighbors to each element */
    int *nn;       /**< number of neighbors to each element */
    double **xyd;  /**< distance between neighbors in xy-plane */
    double **xd;   /**< signed distance between neighbors along x-axis */
    double **yd;   /**< signed distance between neighbors along y-axis */
    int *capacity; /**< keeps track of memory chunks (don't change this) */
} nb_list;

/**
    Creates a neigbor list based on a set of coordinates with
    corresponding sphere radii.

    Implemented using cell lists, giving O(N) performance. Should be
    freed with freesasa_nb_free(). For efficient calculations
    using this list the members of the returned struct should be used
    directly and not freesasa_nb_contact().

    @param coord a set of coordinates
    @param radii radii for the coordinates
    @return a neigbor list. Returns NULL if either argument is null or
      if there were any problems constructing the list (see error
      messages).
 */
nb_list *
freesasa_nb_new(const coord_t *coord,
                const double *radii);

/**
    Frees a neigbor list created by freesasa_nb_new().

    @param nb The neigbor list to free
 */
void freesasa_nb_free(nb_list *nb);

/**
    Checks if two atoms are in contact. Only included for reference.

    @param nb The neigbor list
    @param i Index of first coordinate
    @param j Index of second coordinate
    @return 1 if contact, 0 else.
 */
int freesasa_nb_contact(const nb_list *nb,
                        int i,
                        int j);

#endif /* FREESASA_NB_H */
