sasalib
=======

C-library for calculating Solvent Accessible Surface Areas.

License: GPLv3 (see file COPYING). Copyright: Simon Mitternacht 2013.

This code has been tested for functionality for a small set of PDBs. 
A standalone version can be compiled by running 'make'. The command 
'./calc_sasa -h' will give information about usage.

Incomplete description of the algorithms and the API can be found in
doc/manual.tex.

It is my intention to keep the interface so general and flexible that
this library can be used by any C/C++ program and to provide a
standalone program with a simple command line interface with some
configuration capabilities.

So far the algorithms by Lee & Richards and Shrake & Rupley have been
implemented. Verification has been done by comparing the results of
the two calculations and by visual inspection of the surfaces found by
them. For high resolution versions of the algorithms, the calculations
give identical results.

The OONS atom-classification and radii are used by default (Ooi et al. 
PNAS 1987). Users should be easily able to add their own schemes by 
defining new classification functions. Are there other standard radii 
that could be added to library?

Has only been tested with GNU compiler, both in Linux and Mac OS X, but 
should work with any C99-compatible compiler. Only requires standard C 
and GNU libraries. 
