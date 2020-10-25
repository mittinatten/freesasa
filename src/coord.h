#ifndef FREESASA_COORD_H
#define FREESASA_COORD_H

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/**
    @file
    @author Simon Mitternacht

    This is only for interal use, error handling is done by asserts.
    The header provides functions to store, copy and modify
    coordinates through the type ::coord_t.

    The distance calculation functions (freesasa_dist(),
    freesasa_dist2() and freesasa_dist2_12()) are useful for code that
    is not performance critical. When efficieny is a priority it is
    better to use freesasa_coord_all() to obtain a pointer to the
    array of coordinates and calculate directly using this..
 */

/** Store a set of 3-dimensional coordinates in a contiguous array */
typedef struct coord_t {
    /** number of 3-vectors */
    int n;

    /** If these coordinates are only a link to an externally stored
        array this is 1, else 0. If it it is set, the coordinates can
        not be changed and the array not freed. */
    int is_linked;

    /** array of all coordinates, dimension 3*n,
        x_1,y_1,z_1,...,x_n,y_n,z_n. */
    double *xyz;
} coord_t;

/**
    Initialize new ::coord_t object.

    Return value is dynamically allocated, should be freed with
    freesasa_coord_free().

    @return An empty ::coord_t object. Returns NULL if out of
    memory.
 */
coord_t *
freesasa_coord_new(void);

/**
   Free resources allocated by ::coord_t object.

   Will not free the coordinate array itself if it was initialized by
   freesasa_coord_new_linked().
 */
void freesasa_coord_free(coord_t *coord);

/**
   Clone coordinates.

   Creates a new ::coord_t object that is a copy of the
   argument `src`.

   Return value is dynamically allocated, should be freed with
   freesasa_coord_free().

   @param src Coordinates to be copied.
   @return Copy of coordinates. NULL if out of memory.
 */
coord_t *
freesasa_coord_clone(const coord_t *src);

/**
   Copy coordinates

   Copies value of source coordinates to target.

   @param target Target coordinates
   @param src Source coordinates
   @return ::FREESASA_FAIL if size mismatch between src and target,
     else ::FREESASA_SUCCESS
 */
int freesasa_coord_copy(coord_t *target, const coord_t *src);

/**
    Creates a `const` ::coord_t-object that is linked to an array of
    coordinates owned by callee.

    This allows repeated calculations on a changing object without
    reinitialization. The returned ::coord_t-pointer is not
    explicitly const, to allow it to be freed later, but objects
    initiated through this interface will not change their
    coordinates. It also allows ::coord_t to be used as a wrapper
    for an array controlled by the caller.

    Return value is dynamically allocated, should be freed with
    freesasa_coord_free().

    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return New linked ::coord_t object. NULL if out of memory.
 */
coord_t *
freesasa_coord_new_linked(const double *xyz,
                          int n);

/**
    Append coordinates to ::coord_t object from one array.

    @param coord A ::coord_t object
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int freesasa_coord_append(coord_t *coord,
                          const double *xyz,
                          int n);

/**
    Append coordinates to ::coord_t object from three
    separate arrays.

    @param coord A ::coord_t object
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int freesasa_coord_append_xyz(coord_t *coord,
                              const double *x,
                              const double *y,
                              const double *z,
                              int n);

/**
    Set given coordinate.

    Bounds-checking of i is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::coord_t object
    @param i Index
    @param xyz Array with coordinates x,y,z.
 */
void freesasa_coord_set_i(coord_t *coord,
                          int i,
                          const double *xyz);

/**
    Set given coordinate.

    Bounds-checking of i is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::coord_t object
    @param i Index
    @param x x-coordinate.
    @param y y-coordinate.
    @param z z-coordinate
 */
void freesasa_coord_set_i_xyz(coord_t *coord,
                              int i,
                              double x,
                              double y,
                              double z);

/**
    Reset everything.

    Allocates memory to allow array of size n.

    @param coord A ::coord_t object
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int freesasa_coord_set_all(coord_t *coord,
                           const double *xyz,
                           int n);

/**
    Reset everything.

    Allocates memory to allow array of size n.

    @param coord A ::coord_t object
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int freesasa_coord_set_all_xyz(coord_t *coord,
                               const double *x,
                               const double *y,
                               const double *z,
                               int n);

/**
    Set length of a given coordinate vector. Useful for test-points in
    S&R.

    Bounds-checking of `i` and `length` is handled by asserts for
    efficiency, i.e. only done in debug-mode.

    @param coord A ::coord_t object
    @param i Index
    @param length Desired length (>= 0)
 */
void freesasa_coord_set_length_i(coord_t *coord,
                                 int i,
                                 double length);

/**
    Set length of all coordinate vectors to the same length.

    This means all coordinates are on the same sphere.

    Bounds-checking of length is handled by asserts for
    efficiency, i.e. only done in debug-mode.

    @param coord A ::coord_t object
    @param length Desired length (>= 0)
 */
void freesasa_coord_set_length_all(coord_t *coord,
                                   double length);

/**
    Coordinates for a given index.

    Bounds-checking of `i` is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::coord_t object
    @param i Index
    @return Array with coordinates x,y,z
 */
const double *
freesasa_coord_i(const coord_t *coord,
                 int i);

/**
    Calculate distance between two coordinate vectors. For speed,
    arguments aren't checked.

    @param coord A ::coord_t object
    @param i first index
    @param j second index
    @return Distance between coorsinate i and j.
*/
double freesasa_coord_dist(const coord_t *coord,
                           int i,
                           int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors. For
    speed, arguments aren't checked.

    @param coord A ::coord_t object
    @param i First index
    @param j Second index
    @return Square distance between coordinate i and j
*/
double
freesasa_coord_dist2(const coord_t *coord,
                     int i,
                     int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors in
    separate coordinate sets. For speed, arguments aren't checked.

    @param c1 First set of coordinates
    @param c2 Second set of coordinates
    @param i1 Index in first set
    @param i2 Index in second set
    @return Square distance between coordinates i1 and i2
*/
double freesasa_coord_dist2_12(const coord_t *c1,
                               const coord_t *c2,
                               int i1,
                               int i2)
    __attrib_pure__;

/**
    All coordinates as an array.

    @param coord A ::coord_t object
    @return Array of coordinates x1,y1,z1,x2,y2,z2,...
 */
const double *
freesasa_coord_all(const coord_t *coord) __attrib_pure__;

/**
    Number of coordinates.

    @param coord A ::coord_t object
    @return Number of coordinates
 */
int freesasa_coord_n(const coord_t *coord) __attrib_pure__;

/**
    Translate all coordinates by same vector.

    @param coord A ::coord_t object
    @param xyz Array describing translation vector x,y,z.
 */
void freesasa_coord_translate(coord_t *coord,
                              const double *xyz);

/**
    Translate all coordinates by same vector.

    @param coord A ::coord_t object
    @param x x-coordinate of translation vector
    @param y y-coordinate of translation vector
    @param z z-coordinate of translation vector
 */
void freesasa_coord_translate_xyz(coord_t *coord,
                                  double x,
                                  double y,
                                  double z);

/**
    Scale all coordinates by given factor.

    @param coord A ::coord_t object
    @param a Factor to scale by
 */
void freesasa_coord_scale(coord_t *coord,
                          double a);

#undef __attrib_pure__

#endif
