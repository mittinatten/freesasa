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

#ifndef SASALIB_STRMAP_H
#define SASALIB_STRMAP_H

/** Data structure to store a map of strings to void
    pointers. Typically these strings will be names of residues or
    atoms, i.e. there will be 20-30 different values, has not been
    tested or optimized for large maps. Leading and trailing
    whitespace is ignored in keys. Wrappers are provided for the case
    when values are real numbers. */

typedef struct sasalib_strmap_ sasalib_strmap_t;

/** Create new strmap_t-object. */
sasalib_strmap_t* sasalib_strmap_new();

/** Free strmap_t-object. */
int sasalib_strmap_free(sasalib_strmap_t *map);

/** Set value for specified key. If the key is an empty string (or
    only whitespace), the function returns SASALIB_WARN and the key is
    ignored. */
int sasalib_strmap_set(sasalib_strmap_t *map, 
		       const char* key, void *value);

/** Delete an entry from a map. Returns SASALIB_WARN if the string
    does not exist in the map. */
int sasalib_strmap_delete(sasalib_strmap_t *map, const char *key);

/** Returns 1 if the string exists in the map, 0 if not. */
int sasalib_strmap_exists(const sasalib_strmap_t *map, const char *key);

/** Returns value corresponding to key. Returns NULL if key is
    unknown. Nothing stops keys from having value NULL though. */
void* sasalib_strmap_value(const sasalib_strmap_t *map, 
			   const char *key);

/** Returns the number of elements in the map. */
int sasalib_strmap_size(const sasalib_strmap_t *map);


/** Creates array of all registered key. The array and each string
    will be dynamically allocated, up to the user to free memory. The
    size of the array is the same as that returned by
    sasalib_rehash_size(1). */
char** sasalib_strmap_keys(const sasalib_strmap_t *map);

/** Wrapper for set-function when the value is a real number.*/
int sasalib_strmap_real_set(sasalib_strmap_t *map,
			    const char* key, double value);

/** Second wrapper for set-function when the value is a real
    number. If key already exists the value is added to the present
    value. If the key is an empty string, the function returns
    SASALIB_WARN. */
int sasalib_strmap_real_add(sasalib_strmap_t *map,
			    const char* key, double value);

/** Wrapper for when the value is a real number. */
double sasalib_strmap_real_value(const sasalib_strmap_t *map,
				 const char* key);

#endif
