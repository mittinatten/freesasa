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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sasalib.h"
#include "pdb.h"
#include "strmap.h"

#define TABLE_SIZE 27

#define T sasalib_strmap_t

typedef struct element_ {
    char *key;
    void *value;
    int owned; /* owned = 1 -> free_element(1) will free value 
		  owned = 0 -> user has to free value, 
		  cannot be set explicitly by user. */
    struct element_ *next;
} element_t;

struct sasalib_strmap_ {
    element_t *table[TABLE_SIZE];
    size_t n;
};

extern size_t sasalib_trim_whitespace(char *target, const char *src, 
				      size_t length);

static int hash(const char* s) 
{
    //gives approx alphabetical tables
    int ret = ((int)(s[0]-'A') % TABLE_SIZE);
    if (ret < 0) ret += TABLE_SIZE;
    return ret;
}

T* sasalib_strmap_new() 
{
    T *h = (T*) malloc(sizeof(T));
    h->n = 0;
    for (int i = 0; i < TABLE_SIZE; ++i) h->table[i] = NULL;
    return h;
}

static void free_element(element_t *e)
{
    assert(e && "attempting to free null-pointer");
    if (e->next != NULL) free_element(e->next);
    free(e->key);
    if (e->owned) free(e->value);
    free(e);
}

int sasalib_strmap_free(T *map) 
{
    assert(map && "attempting to free null-pointer");
    element_t *e;
    for (int i = 0; i < TABLE_SIZE; ++i) {
	if ((e = map->table[i]))
	    free_element(e);
    }
    free(map);
    return SASALIB_SUCCESS;
}

// the_e will point to NULL if key does not exist
static void strmap_find(const T *map, const char *key, 
			element_t **the_e, element_t **e_before)
{
    assert(map);
    char *buf = (char*) malloc(strlen(key) + 1);
    sasalib_trim_whitespace(buf,key,strlen(key));
    
    int i = hash(buf);
    element_t *e1 = map->table[i], *e2 = NULL;

    //find entry 
    while (e1 != NULL && strcmp(buf,e1->key)) { 
	e2 = e1; 
	e1 = e1->next;
    }
    *the_e = e1;
    *e_before = e2;
    free (buf);
}
static int strmap_internal_set(T *map, const char *key, void *value, 
			       int owned)
{
    assert(map);
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (e1) {
	e1->value = value;
	e1->owned = owned;
	return SASALIB_SUCCESS;
    }

    char *buf = (char*) malloc(strlen(key) + 1);
    sasalib_trim_whitespace(buf,key,strlen(key));
    if (strlen(buf) == 0) return SASALIB_WARN;

    //generate new element
    element_t *e = (element_t*) malloc(sizeof(element_t)), *e_prev;
    e->key = buf;
    e->next = NULL;
    e->value = value;
    e->owned = owned;
    
    // find and fill empty slot
    int i = hash(buf);
    if ((e_prev = map->table[i])) {
	while (e_prev->next != NULL) {
	    e_prev = e_prev->next;
	}
	e_prev->next = e;
    } else {
	map->table[i] = e;
    }
    ++map->n;
    return SASALIB_SUCCESS;
}

int sasalib_strmap_set(T *map, const char *key, void *value)
{
    return strmap_internal_set(map,key,value,0);
}

int sasalib_strmap_delete(T *map, const char *key)
{
    //find element
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (e1 == NULL) return SASALIB_WARN;

    //delete element and relink list
    free(e1->key);
    if (e2) e2->next = e1->next;
    else map->table[hash(key)] = NULL;
    if (e1->owned) free(e1->value);
    free(e1);
    --map->n;
    return SASALIB_SUCCESS;
}

void* sasalib_strmap_value(const T *map, const char *key)
{
    assert(map);
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (!e1) return NULL;
    return e1->value;
}


int sasalib_strmap_exists(const T *map, const char *key)
{
    assert(map);
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (e1 == NULL) return 0;
    return 1;
}
int sasalib_strmap_size(const T *map) 
{
    assert(map);
    return map->n;
}

char** sasalib_strmap_keys(const T *map)
{
    element_t *e;
    int pos = 0;
    char **a = (char**) malloc(sizeof(char*)*map->n);
    for (int i = 0; i < TABLE_SIZE; ++i) {
	e = map->table[i];
	if (e) {
	    while (e) {
		assert(pos < map->n);
		a[pos] = (char*)malloc(strlen(e->key)+1);
		strcpy(a[pos],e->key);
		e = e->next;
		++pos;
	    }
	}
    }
    return a;
}
int sasalib_strmap_real_set(sasalib_strmap_t *map,
			    const char* key, double value)
{
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (!e1) {
	double *v = (double*)malloc(sizeof(double));
	*v = value;
	return strmap_internal_set(map,key,v,1);
    }
    *((double*)e1->value) = value;
    return SASALIB_SUCCESS;
}

int sasalib_strmap_real_add(T *map, const char* key, double value)
{
    element_t *e1, *e2;
    strmap_find(map,key,&e1,&e2);
    if (!e1) {
	double *v = (double*)malloc(sizeof(double));
	*v = value;
	return strmap_internal_set(map,key,v,1);
    }
    // the only difference from sasalib_strmap_real_set(3)
    *((double*)e1->value) += value;
    return SASALIB_SUCCESS;
}
double sasalib_strmap_real_value(const T *map, const char* key)
{
    return *(double*)sasalib_strmap_value(map,key);
}

