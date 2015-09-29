@page Geometry Geometry of Lee & Richards' algorithm

This page explains the geometry of the calculations in L&R
and can be used to understand the source code. As far as possible the
code uses similar notation to the formulas here.

We will use the following notation: An atom \f$i\f$ has a van der
Waals radius \f$r_i\f$, the rolling sphere (or *probe*) has radius
\f$r_\text{P}\f$ and when these are added we get an extended radius
\f$R_i = r_i + r_\text{P}\f$. The sphere of radius \f$R_i\f$ centered
at the position of atom \f$i\f$ represents the volume not accessible
to the center of the probe. The SASA for a molecule is then obtained
by calculating the non-buried surface area of the extended spheres.

The L&R algorithm calculates the surface area by slicing the
protein, calculating the length of the solvent exposed contours in
each slice and then adding up the length multiplied by the slice
thickness.

![Slice in atom](../fig/lnr_slice.svg)

Divide atom \f$i\f$ into \f$n\f$ slices, orthogonal to an arbitrary
axis, of thickness \f$\delta = 2R_i/n\f$. The position of the middle
of the slice along that axis is denoted \f$z\f$, and the center of
atom \f$i\f$, along the same axis, is at \f$z_i\f$. In each slice, the
atom is thus a circle of radius \f[R_i^\prime =
\sqrt{R_i^2-(z-z_i)^2}\,.\f] These circles are either completely
buried inside neighboring atoms, completely exposed, or partially
exposed.

![Overlap of circles](../fig/lnr_circles.svg)

The exposed arc lengths for each atom can be calculated exactly. For
each pair of atoms \f$i,j\f$, the distance between their centers
projected on the slice is \f$d_{ij}\f$ (independent of \f$z\f$). If
\f$d_{ij} > R_i^\prime + R_j^\prime\f$, there is no overlap. If
\f$d_{ij} < R_j^\prime - R_i^\prime\f$ circle \f$i\f$ is completely
inside \f$j\f$ (and the other way around). If \f$d_{ij}\f$ lies
between these two cases the angle of circle \f$i\f$ that is buried due
to circle \f$j\f$ is

\f[ \alpha = 2\arccos \bigl[({R_i^\prime}^2_{\,}
+ d_{ij}^2 - {R_{j}^\prime}^2_{\,})/(2R_i^\prime d_{ij})\bigr].  \f]

If the middle point of this arc on the circle is at an angle
\f$\beta\f$, the arc spans the interval
\f$[\beta-\alpha/2,\beta+\alpha/2]\f$. By adding up these arcs and
taking into account any overlap between them we get the total buried
angle \f$\gamma\f$ in this slices. The exposed arc angle for this atom
and slice is thus \f$2\pi-\gamma\f$ and the total SASA of that atom

\f[ A_i =R_i \delta \!\! \sum_{s\in\text{slices}} \!\!
\left[2\pi-\gamma_s\right]\,.  \f]

The angle is multiplied by \f$R_i\f$ (not \f$R_i^\prime\f$) to give
the area of a conical frustum circumscribing the sphere at the
slice. Finally, the total area \f$A\f$ is the sum of all \f$A_i\f$.

In FreeSASA, the L\&R SASA calculation begins by finding overlapping
spheres and storing the contacts in an adjacency list. It then
iterates through all the slices of each atom and checks for overlap
with adjacent atoms in each slice, and adds up the exposed arcs to
calculate the atom's contribution to the SASA of the slice. The
calculations for each atom are completely independent and can thus be
parallelized over an arbitrary number of threads, whereas the
calculation of adjacency lists has not been parallelized.