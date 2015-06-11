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
  GNU General Public License forr more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "freesasa.h"
#include "adjacency.h"

typedef struct freesasa_cell freesasa_cell;
struct freesasa_cell {
    freesasa_cell *nb[17]; //! includes self, only forward neighbors
    int atom[FREESASA_ATOMS_PER_CELL];
    int n_nb; //! number of neighbors to cell
    int n_atoms; //! number of atoms in cell
};

struct freesasa_cell_list {
    freesasa_cell *cell; //! indexed geometrically
    size_t n; //! number of cells
    size_t nx, ny, nz; //! number of cells along each axis
    double d; //! cell size
    double x_max, x_min;
    double y_max, y_min;
    double z_max, z_min;
};

extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

void get_bounds(freesasa_cell_list *c, const freesasa_coord *coord)
{
    const size_t n = freesasa_coord_n(coord);
    double d = c->d;
    const double *v = freesasa_coord_i(coord,0);
    double x=v[0],X=v[0],y=v[1],Y=v[1],z=v[2],Z=v[2];
    for (size_t i = 1; i < n; ++i) {
        v = freesasa_coord_i(coord,i);
        x = fmin(v[0],x);
        X = fmax(v[0],X);
        y = fmin(v[1],y);
        Y = fmax(v[1],Y);
        z = fmin(v[2],z);
        Z = fmax(v[2],Z);
    }
    c->x_min = x - d/2.;
    c->x_max = X + d/2.;
    c->y_min = y - d/2.;
    c->y_max = Y + d/2.;
    c->z_min = z - d/2.;
    c->z_max = Z + d/2.;
    c->nx = ceil((c->x_max - c->x_min)/d);
    c->ny = ceil((c->y_max - c->y_min)/d);
    c->nz = ceil((c->z_max - c->z_min)/d);
    c->n = c->nx*c->ny*c->nz;
}

static inline int flat_idx(const freesasa_cell_list *c,
                           int ix, int iy, int iz)
{
    assert(ix >= 0 && ix < c->nx);
    assert(iy >= 0 && iy < c->ny);
    return ix + c->nx*(iy + c->ny*iz);
}

/*
static void calc_centers(const freesasa_cell_list *c)
{
    for (int ix = 0; ix < c->nx; ++ix) {
        for (int iy = 0; iy < c->ny; ++iy) {
            for (int iz = 0; iz < c->nz; ++iz) {
                int idx = flat_idx(c,ix,iy,iz);
                freesasa_cell* cell = &c->cell[idx];
                cell->x = (ix*c->d - c->x_min);
                cell->y = (iy*c->d - c->y_min);
                cell->z = (iz*c->d - c->z_min);
            }
        }
    }
}
*/

// Fill the neighbor list for a given cell
static void fill_nb(freesasa_cell_list *c,
                    int ix, int iy, int iz)
{
    freesasa_cell *cell = &c->cell[flat_idx(c,ix,iy,iz)];
    int n = 0;
    int xmin = ix > 0 ? ix - 1 : 0;
    int xmax = ix < c->nx - 1 ? ix + 1 : ix;
    int ymin = iy > 0 ? iy - 1 : 0;
    int ymax = iy < c->ny - 1 ? iy + 1 : iy;
    int zmin = iz > 0 ? iz - 1 : 0;
    int zmax = iz < c->nz - 1 ? iz + 1 : iz;
    for (int i = xmin; i <= xmax; ++i) {
        for (int j = ymin; j <= ymax; ++j) {
            for (int k = zmin; k <= zmax; ++k) {
                /* Scalar product between (i-ix,j-iy,k-iz) and (1,1,1) should
                   be non-negative. Using only forward neighbors means
                   there's no double counting when comparing cells */
                if (i-ix+j-iy+k-iz >= 0) { 
                    cell->nb[n] = &c->cell[flat_idx(c,i,j,k)];
                    ++n;
                }
            }
        }
    }
    cell->n_nb = n;
    assert(n > 0);
}

static void get_nb(freesasa_cell_list *c)
{
    for (size_t ix = 0; ix < c->nx; ++ix) {
        for (size_t iy = 0; iy < c->ny; ++iy) {
            for (size_t iz = 0; iz < c->nz; ++iz) {
                fill_nb(c,ix,iy,iz);
            }
        }
    }
}

static int fill_cells(freesasa_cell_list *c, const freesasa_coord *coord)
{
    const double d = c->d;
    const int nx = c->nx, ny = c->ny, nz = c->nz;
    for (size_t i = 0; i < c->n; ++i) {
        c->cell[i].n_atoms = 0;
    }
    for (size_t i = 0; i < freesasa_coord_n(coord); ++i) {
        const double *v = freesasa_coord_i(coord,i);
        freesasa_cell *cell;
        int ix = (int)((v[0] - c->x_min)/d);
        int iy = (int)((v[1] - c->y_min)/d);
        int iz = (int)((v[2] - c->z_min)/d);
        int idx = flat_idx(c,ix,iy,iz);
        cell = &c->cell[idx];
        cell->atom[cell->n_atoms] = i;
        ++cell->n_atoms;
        if (cell->n_atoms > FREESASA_ATOMS_PER_CELL) {
            return FREESASA_FAIL;
        }
    }
    return FREESASA_SUCCESS;
}

freesasa_cell_list* freesasa_cell_list_new(double cell_size,
                                           const freesasa_coord *coord)
{
    assert(cell_size > 0);
    assert(coord);
    freesasa_cell_list *c = malloc(sizeof(freesasa_cell_list));
    assert(c);
    c->d = cell_size;
    get_bounds(c,coord);
    c->cell = malloc(sizeof(freesasa_cell)*c->n);
    assert(c->cell);
    if (fill_cells(c,coord) == FREESASA_SUCCESS) {
        get_nb(c);
    } else {
        freesasa_cell_list_free(c);
        c = NULL;
        freesasa_fail("%s: Too many atoms per cell, "
                      "try smaller cell size\n",
                      __func__);
    }
    return c;
}
void freesasa_cell_list_free(freesasa_cell_list *c)
{
    if (c) {
        free(c->cell);
        free(c);
    }
}

//assumes max value in a is positive
static double max_array(const double *a,size_t n)
{
    double max = 0;
    for (size_t i = 0; i < n; ++i) {
        max = fmax(a[i],max);
    }
    return max;
}

static freesasa_adjacency *freesasa_adjacency_alloc(size_t n)
{
    freesasa_adjacency *adj = malloc(sizeof(freesasa_adjacency));
    assert(adj);
    adj->n = n;
    adj->nb = malloc(sizeof(int *)*n);
    adj->nn = malloc(sizeof(size_t)*n);
    adj->nb_xyd = malloc(sizeof(double *)*n);
    adj->nb_xd = malloc(sizeof(double *)*n);
    adj->nb_yd = malloc(sizeof(double *)*n);
    
    assert(adj->nb); assert(adj->nn);
    assert(adj->nb_xyd); assert(adj->nb_xd); assert(adj->nb_yd);

    for (size_t i=0; i < n; ++i) {
        adj->nn[i] = 0;
        adj->nb[i] = NULL;
        adj->nb_xyd[i] = NULL;
        adj->nb_xd[i] = NULL;
        adj->nb_yd[i] = NULL;
    }

    return adj;
}
void freesasa_adjacency_free(freesasa_adjacency *adj)
{
    if (adj) {
        for (size_t i = 0; i < adj->n; ++i) {
            free(adj->nb[i]);
        }
        free(adj->nb);
        free(adj->nn);
        free(adj->nb_xyd);
        free(adj->nb_xd);
        free(adj->nb_yd);
        free(adj);
    }
}


static void add_pair(freesasa_adjacency *adj,int i, int j,
                     double dx, double dy)
{
    assert(i != j);
    int **nb = adj->nb;
    size_t *nn = adj->nn;
    double **nb_xyd = adj->nb_xyd;
    double **nb_xd = adj->nb_xd;
    double **nb_yd = adj->nb_yd;
    ++nn[i]; ++nn[j];
    nb[i] = realloc(nb[i],sizeof(int)*nn[i]); 
    nb[j] = realloc(nb[j],sizeof(int)*nn[j]);
    nb_xyd[i] = realloc(nb_xyd[i],sizeof(double)*nn[i]);
    nb_xyd[j] = realloc(nb_xyd[j],sizeof(double)*nn[j]);
    nb_xd[i] = realloc(nb_xd[i],sizeof(double)*nn[i]);
    nb_xd[j] = realloc(nb_xd[j],sizeof(double)*nn[j]);
    nb_yd[i] = realloc(nb_yd[i],sizeof(double)*nn[i]);
    nb_yd[j] = realloc(nb_yd[j],sizeof(double)*nn[j]);
    assert(nb[i] && nb[j]);
    assert(nb_xyd[i] && nb_xyd[j]);
    assert(nb_xd[i] && nb_xd[i] && nb_yd[i] && nb_yd[j]);
    nb[i][nn[i]-1] = j;
    nb[j][nn[j]-1] = i;
    
    double d = sqrt(dx*dx+dy*dy);
    nb_xyd[i][nn[i]-1] = d;
    nb_xyd[j][nn[j]-1] = d;
    
    nb_xd[i][nn[i]-1] = dx;
    nb_xd[j][nn[j]-1] = -dx;
    nb_yd[i][nn[i]-1] = dy;
    nb_yd[j][nn[j]-1] = -dy;
}

// also checks for internal contacts in first cell (ci);
static void calc_cell_adjacency(freesasa_adjacency *adj,
                                const freesasa_coord* coord,
                                const double *radii,
                                const freesasa_cell *ci,
                                const freesasa_cell *cj)
{
    const double *v = freesasa_coord_all(coord);
    double ri, rj, xi, yi, zi, xj, yj, zj,
        dx, dy, dz, cut2;
    int count_neighbors[ci->n_atoms];
    int i,j,ia,ja;
    for (i = 0; i < ci->n_atoms; ++i) {
        ia = ci->atom[i];
        ri = radii[ia];
        xi = v[ia*3]; yi = v[ia*3+1]; zi = v[ia*3+2];
        if (ci == cj) j = i+1;
        else j = 0;
        // the following loop is performance critical
        for (; j < cj->n_atoms; ++j) {
            ja = cj->atom[j];
            assert (ia != ja);
            rj = radii[ja];
            xj = v[ja*3]; yj = v[ja*3+1]; zj = v[ja*3+2];
            cut2 = (ri+rj)*(ri+rj);
            if ((xj-xi)*(xj-xi) > cut2 ||
                (yj-yi)*(yj-yi) > cut2 ||
                (zj-zi)*(zj-zi) > cut2) {
                continue;
            }
            dx = xj-xi; dy = yj-yi; dz = zj-zi;
            if (dx*dx + dy*dy + dz*dz < cut2) {
                add_pair(adj,ia,ja,dx,dy);
            }
        }
    }
}
                             

static void freesasa_adjacency_calc(freesasa_adjacency *adj,
                                    freesasa_cell_list *cl,
                                    const freesasa_coord *coord,
                                    const double *radii)
{
    size_t nc = cl->n;
    for (size_t ic = 0; ic < nc; ++ic) {
        const freesasa_cell *ci = &cl->cell[ic];
        for (int jc = 0; jc < ci->n_nb; ++jc) {
            if (ic == jc) continue;
            const freesasa_cell *cj = ci->nb[jc];
            calc_cell_adjacency(adj,coord,radii,ci,cj);
        }
    }
}

static int check_consistency(const freesasa_coord* coord,
                             const freesasa_cell_list *celllist,
                             const freesasa_adjacency *adj,
                             const double* radii)
{
    size_t na = freesasa_coord_n(coord);
    assert(adj->n == na);
    size_t n_in_cells = 0;
    for (int i = 0; i < celllist->n; ++i) {
        freesasa_cell *ci = &celllist->cell[i];
        n_in_cells += ci->n_atoms;
        for (int j = 0; j < ci->n_nb; ++j) {
            freesasa_cell *cj = ci->nb[j];
            for (int k = 0; k < ci->n_atoms; ++k) {
                int ka = ci->atom[k];
                for (int l = 0; l < cj->n_atoms; ++l) {
                    int la = cj->atom[l];
                    assert(freesasa_coord_dist(coord,ka,la) <= 2*sqrt(3)*celllist->d);
                }
            }
        }
    }
    assert(n_in_cells == na);
    int c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    for (size_t i = 0; i < na; ++i) {
        for (size_t j = 0; j < na; ++j) {
            if (i == j) continue;
            double d = freesasa_coord_dist(coord,i,j);
            if (d < radii[i] + radii[j]) {
                if(freesasa_adjacency_contact(adj,i,j)==0) {
                    ++c2;
                    printf("contact mismatch %d %d %f %f\n",
                           i, j, d, radii[i] + radii[j]);
                } else { ++c4; }
                ++c1;
            }
            if (freesasa_adjacency_contact(adj,i,j) == 1) ++c3;
        }
    }
    printf("actual: %d, c2: %d, c3: %d, c4: %d, c2+c3: %d\n",c1,c2,c3,c4,c2+c3);
    
}

freesasa_adjacency *freesasa_adjacency_new(const freesasa_coord* coord,
                                           const double *radii)
{
    assert(coord);
    assert(radii);
    size_t n = freesasa_coord_n(coord);
    freesasa_adjacency *adj = freesasa_adjacency_alloc(n);
    double cell_size = 2*max_array(radii,n);
    freesasa_cell_list *c = freesasa_cell_list_new(cell_size,coord);
    assert(c);
    freesasa_adjacency_calc(adj,c,coord,radii);
    check_consistency(coord,c,adj,radii);
    freesasa_cell_list_free(c);
    
    return adj;
}

int freesasa_adjacency_contact(const freesasa_adjacency *adj,
                               int i, int j)
{
    assert(adj);
    assert(i < adj->n && i >= 0);
    assert(j < adj->n && j >= 0);
    for (size_t k = 0; k < adj->nn[i]; ++k) {
        if (adj->nb[i][k] == j) return 1;
    }
    return 0;
}

