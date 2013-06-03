#include "vector3.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

const vector3 vector_zero = {0, 0, 0};
const vector3 vector_ex = {1, 0, 0};
const vector3 vector_ey = {0, 1, 0};
const vector3 vector_ez = {0, 0, 1};

void vector3_set(vector3* v, double x, double y, double z) 
{
    v->x = x;
    v->y = y;
    v->z = z;
}

void vector3_set_polar(vector3* v, double polar, double az, double r) 
{
    v->x = r*cos(az)*sin(polar);
    v->y = r*sin(az)*sin(polar);
    v->z = r*cos(polar);
}

double vector3_mag(const vector3* v) 
{
    return sqrt(vector3_dot(v,v));
}

double vector3_mag2(const vector3* v) 
{
    return vector3_dot(v,v);
}

void vector3_setlength(vector3* v, double L) 
{
    double oldL = vector3_mag(v);
    vector3_multiply(v,L/oldL);
}

void vector3_multiply(vector3* v, double m) 
{
    v->x *= m;
    v->y *= m;
    v->z *= m;
}

void vector3_diff(vector3* u, const vector3* v1, const vector3* v2) 
{
    u->x = v2->x - v1->x;
    u->y = v2->y - v1->y;
    u->z = v2->z - v1->z;
}

void vector3_sum(vector3* sum, const vector3* v1, const vector3* v2) 
{
    sum->x = v1->x + v2->x;
    sum->y = v1->y + v2->y;
    sum->z = v1->z + v2->z;
}

void vector3_add(vector3* v1, const vector3* v2) 
{
    v1->x += v2->x;
    v1->y += v2->y;
    v1->z += v2->z;
}

void vector3_sub(vector3* v1, const vector3* v2) 
{
    v1->x -= v2->x;
    v1->y -= v2->y;
    v1->z -= v2->z;
}

double vector3_dist2(const vector3* v1, const vector3* v2) 
{
    double dx = v1->x - v2->x;
    double dy = v1->y - v2->y;
    double dz = v1->z - v2->z;
    return dx*dx + dy*dy + dz*dz;
}


void vector3_copy(vector3* u, const vector3* v) 
{
    u->x = v->x;
    u->y = v->y;
    u->z = v->z;
}

double vector3_perdist2(const vector3* v1, const vector3* v2, const double L) 
{
    const double L2 = L/2.;
    double dx = fabs(v1->x - v2->x);
    double dy = fabs(v1->y - v2->y);
    double dz = fabs(v1->z - v2->z);
    while (dx > L2) { dx -= L; }
    while (dy > L2) { dy -= L; }
    while (dz > L2) { dz -= L; }
    return dx*dx + dy*dy + dz*dz;
}

double vector3_dot(const vector3* v1, const vector3* v2) 
{
    return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

void vector3_cross(vector3* u, const vector3* v1, const vector3* v2) 
{
    u->x = v1->y*v2->z - v1->z*v2->y;
    u->y = v1->z*v2->x - v1->x*v2->z;
    u->z = v1->x*v2->y - v1->y*v2->x;
}

double vector3_theta(const vector3* v1, const vector3* v2) 
{
    if (vector3_mag(v1) == 0. || vector3_mag(v2) == 0) {
	return 0.;
    }
    return acos(vector3_dot(v1,v2)/(vector3_mag(v1)*vector3_mag(v2)));
}

double vector3_torsion(const vector3* v1, const vector3* v2, const vector3* v3) 
{
    //using the atan2-version from wikipedia
    double b2 = vector3_mag(v2);
    vector3 cross12, cross23, b2b1;
    vector3_cross(&cross12, v1, v2);
    vector3_cross(&cross23, v2, v3);
    b2b1.x = v1->x*b2;
    b2b1.y = v1->y*b2;
    b2b1.z = v1->z*b2;
    return atan2(vector3_dot(&b2b1, &cross23), vector3_dot(&cross12, &cross23));
}

void vector3_rotate_block(vector3* v, int i1, int i2, 
			  const vector3* axis, const vector3* origin, double phi) 
{
    //adapted from Anders' code
    double ex, ey, ez, e;
    double Axx, Axy, Axz, Ayy, Ayz, Azz;
    double Bxx, Bxy, Bxz, Byy, Byz, Bzz;
    double Cxx, Cxy, Cxz, Cyx, Cyy, Cyz, Czx, Czy, Czz;
    double cdph = cos(phi), sdph = sin(phi);
    vector3 dv;
    e = vector3_mag(axis);
    ex = axis->x/e;
    ey = axis->y/e;
    ez = axis->z/e;
    Axx = ex*ex; Axy = ex*ey; Axz = ex*ez; 
    Ayy = ey*ey; Ayz = ey*ez; 
    Azz = ez*ez;
    Bxx = Byy = Bzz = cdph;
    Bxy = -ez*sdph; 
    Byz = -ex*sdph; 
    Bxz = ey*sdph;
    Cxx =  Bxx + (1-Bxx)*Axx - Bxy*Axy - Bxz*Axz;
    Cxy =  Bxy + (1-Bxx)*Axy - Bxy*Ayy - Bxz*Ayz;
    Cxz =  Bxz + (1-Bxx)*Axz - Bxy*Ayz - Bxz*Azz;
    Cyx = -Bxy + Bxy*Axx + (1-Byy)*Axy - Byz*Axz;
    Cyy =  Byy + Bxy*Axy + (1-Byy)*Ayy - Byz*Ayz;
    Cyz =  Byz + Bxy*Axz + (1-Byy)*Ayz - Byz*Azz;
    Czx = -Bxz + Bxz*Axx + Byz*Axy + (1-Bzz)*Axz;
    Czy = -Byz + Bxz*Axy + Byz*Ayy + (1-Bzz)*Ayz;
    Czz =  Bzz + Bxz*Axz + Byz*Ayz + (1-Bzz)*Azz;

    for (int j = i1; j <= i2; ++j) {
	vector3_diff(&dv, origin, &v[j]);
	v[j].x = origin->x + Cxx*dv.x + Cxy*dv.y + Cxz*dv.z;
	v[j].y = origin->y + Cyx*dv.x + Cyy*dv.y + Cyz*dv.z;
	v[j].z = origin->z + Czx*dv.x + Czy*dv.y + Czz*dv.z;
    }  
}
void vector3_distance_to_box(vector3* diff, const vector3* v, double L) {
    double dx = 0, dy = 0, dz = 0;
    while (v->x + dx < 0) dx += L;
    while (v->x + dx > L) dx -= L;
    while (v->y + dy < 0) dy += L;
    while (v->y + dy > L) dy -= L;
    while (v->z + dz < 0) dz += L;
    while (v->z + dz > L) dz -= L;
    vector3_set(diff, dx, dy, dz);
}

double vector3_passage_distance(const vector3 *p1, const vector3 *p2, const vector3 *p3)
{
    double t0 = 0;
    vector3 p1p2, p1p3, D;

    // D is distance-vector between the line p1 + (p1-p2)*t and p3 

    vector3_diff(&p1p2, p1, p2);
    vector3_diff(&p1p3, p1, p3);

    // t that minimizes |D| 
    t0 = - vector3_dot(&p1p2, &p1p3) / vector3_mag2(&p1p2);

    // calculate minimal D
    vector3_multiply(&p1p2, t0);
    vector3_sum(&D, &p1p3, &p1p2);
    
    return vector3_mag(&D);
}
