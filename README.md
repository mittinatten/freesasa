sasalib
=======

C-library for calculating Solvent Accessible Surface Areas

This code has not been tested extensively yet, and some pieces are
still missing. A rudimentary standalone version can be compiled by
running 'make' (very primitive Makefile so far). The command
'./calc_sasa -h' will give information about usage.

License: GPLv3 (see file COPYING). Copyright: Simon Mitternacht 2013.

It is my intention to keep the interface so general and felxible that
this library can be used by any C/C++ program and to provide a
standalone program with a simple command line interface with some
configuration capabilities.

So far the Shrake & Rupley algorithm has been implemented. Lee &
Richards' algorithm should also be included. 

The OONS atom-classification and radii are used by default. Users
should be easily able to add their own schemes by defining new
classification functions. Are there other standard radii that could be
added to library?
