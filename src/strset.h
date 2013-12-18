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

#ifndef SASALIB_STRSET_H
#define SASALIB_STRSET_H

/** Data structure to store set of strings. Typically these strings
    will be names of residues or atoms, i.e. there will be 20-30
    different values, has not been tested or optimized for large
    sets.*/

typedef struct sasalib_strset_ sasalib_strset_t;

sasalib_strset_t* sasalib_strset_new();
int sasalib_strset_free(sasalib_strset_t*);

int sasalib_strset_add(sasalib_strset_t*, const char*);
int sasalib_strset_delete(sasalib_strset_t*, const char*);
int sasalib_strset_exists(const sasalib_strset_t*, const char*);
int sasalib_strset_size(const sasalib_strset_t*);

/** Creates array of all registered strings. Array and each string
    will be dynamically allocated, up to the user to free memory. Size
    of array is the same as that returned by
    sasalib_rehash_size(1). */
char** sasalib_strset_array(const sasalib_strset_t*);

#endif
