/*
  Copyright Simon Mitternacht 2013.

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

#ifndef FREESASA_STRMAP_H
#define FREESASA_STRMAP_H

/** Data structure to store a map of strings to void
    pointers. Typically these strings will be names of residues or
    atoms, i.e. there will be 20-30 different values, has not been
    tested or optimized for large maps. Leading and trailing
    whitespace is ignored in keys. Wrappers are provided for the case
    when values are real numbers. */

typedef struct freesasa_strmap_ freesasa_strmap_t;

/** Create new strmap_t-object. */
freesasa_strmap_t* freesasa_strmap_new();

/** Free strmap_t-object. */
int freesasa_strmap_free(freesasa_strmap_t *map);

/** Set value for specified key. If the key is an empty string (or
    only whitespace), the function returns FREESASA_WARN and the key is
    ignored. */
int freesasa_strmap_set(freesasa_strmap_t *map,
                        const char* key, void *value);

/** Delete an entry from a map. Returns FREESASA_WARN if the string
    does not exist in the map. */
int freesasa_strmap_delete(freesasa_strmap_t *map, const char *key);

/** Returns 1 if the string exists in the map, 0 if not. */
int freesasa_strmap_exists(const freesasa_strmap_t *map, const char *key);

/** Returns value corresponding to key. Returns NULL if key is
    unknown. Nothing stops keys from having value NULL though. */
void* freesasa_strmap_value(const freesasa_strmap_t *map,
                            const char *key);

/** Returns the number of elements in the map. */
int freesasa_strmap_size(const freesasa_strmap_t *map);


/** Creates array of all registered key. The array and each string
    will be dynamically allocated, up to the user to free memory. The
    size of the array is the same as that returned by
    freesasa_rehash_size(1). */
char** freesasa_strmap_keys(const freesasa_strmap_t *map);

/** Wrapper for set-function when the value is a real number.*/
int freesasa_strmap_real_set(freesasa_strmap_t *map,
                             const char* key, double value);

/** Second wrapper for set-function when the value is a real
    number. If key already exists the value is added to the present
    value. If the key is an empty string, the function returns
    FREESASA_WARN. */
int freesasa_strmap_real_add(freesasa_strmap_t *map,
                             const char* key, double value);

/** Wrapper for when the value is a real number. */
double freesasa_strmap_real_value(const freesasa_strmap_t *map,
                                  const char* key);

#endif
