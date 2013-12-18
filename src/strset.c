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
#include "strset.h"

#define TABLE_SIZE 27

#define T sasalib_strset_t

typedef struct element_ {
    char *string;
    struct element_ *next;
} element_t;

struct sasalib_strset_ {
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

T* sasalib_strset_new() 
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
    free(e->string);
    free(e);
}

int sasalib_strset_free(T *h) 
{
    assert(h && "attempting to free null-pointer");
    element_t *e;
    for (int i = 0; i < TABLE_SIZE; ++i) {
	if ((e = h->table[i]))
	    free_element(e);
    }
    free(h);
    return SASALIB_SUCCESS;
}

int sasalib_strset_add(T *h, const char *key)
{
    assert(h);

    char *buf = (char*) malloc(strlen(key) + 1);
    sasalib_trim_whitespace(buf,key,strlen(key));
    if (strlen(buf) == 0) return SASALIB_WARN;

    //generate new element
    element_t *e = (element_t*) malloc(sizeof(element_t));
    e->string = buf;
    e->next = NULL;
    
    // find and fill empty slot
    int i = hash(buf);
    if (h->table[i]) {
	element_t *e_prev = h->table[i];
	while (e_prev->next != NULL) {
	    e_prev = e_prev->next;
	}
	e_prev->next = e;
    } else {
	h->table[i] = e;
    }
    ++h->n;
    return SASALIB_SUCCESS;
}

// the_e will point to NULL if key does not exist
static void strset_find(const T *h, const char *key, 
			 element_t **the_e, element_t **e_before)
{
    char *buf = (char*) malloc(strlen(key) + 1);
    sasalib_trim_whitespace(buf,key,strlen(key));
    
    int i = hash(buf);
    element_t *e1 = h->table[i], *e2 = NULL;

    //find entry 
    while (e1 != NULL && strcmp(buf,e1->string)) { 
	e2 = e1; 
	e1 = e1->next;
    }
    *the_e = e1;
    *e_before = e2;
}
int sasalib_strset_delete(T *h, const char *key)
{
    //find element
    element_t *e1, *e2;
    strset_find(h,key,&e1,&e2);
    if (e1 == NULL) return SASALIB_WARN;

    //delete element and relink list
    free(e1->string);
    if (e2) e2->next = e1->next;
    else h->table[hash(key)] = NULL;
    free(e1);
    --h->n;
    return SASALIB_SUCCESS;
}

int sasalib_strset_exists(const T *h, const char *key)
{
    element_t *e1, *e2;
    strset_find(h,key,&e1,&e2);
    if (e1 == NULL) return 0;
    return 1;
}
int sasalib_strset_size(const T *h) 
{
    assert(h);
    return h->n;
}

char** sasalib_strset_array(const T *h)
{
    element_t *e;
    int pos = 0;
    char **a = (char**) malloc(sizeof(char*)*h->n);
    for (int i = 0; i < TABLE_SIZE; ++i) {
	e = h->table[i];
	if (e) {
	    while (e) {
		assert(pos < h->n);
		a[pos] = (char*) malloc(strlen(e->string));
		strcpy(a[pos],e->string);
		e = e->next;
		++pos;
	    }
	}
    }
    return a;
}
