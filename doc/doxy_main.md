FreeSASA
========

These pages document the 
    
  - @ref CLI
  - @ref API "FreeSASA C API"
  - @ref Python "FreeSASA Python interface"
  - @ref Config-file
  - @ref Geometry

The library is licensed under [GPLv3](GPL.md).

Installation instructions can be found in the [README](README.md) file.

@page CLI Command-line Interface

Building FreeSASA creates the binary `freesasa`, which is installed by
`make install`. Calling

    $ freesasa -h

displays a help message listing all options. The following text
explains how to use most of them.

@section CLI-Default Run using defaults

In the following we will use the PDB structure 1UBQ as
an example. To run a simple SASA calculation using default parameters,
simply type:

    $ freesasa 1ubq.pdb

This generates the following output

    name: 1ubq.pdb
    n_atoms: 602
    algorithm: Shrake & Rupley
    probe-radius: 1.400000 A
    n_thread: 2
    n_testpoint: 100

       Total:   4779.51 A2
       Polar:   2236.93 A2
      Apolar:   2542.58 A2
     Nucleic:      0.00 A2
     Unknown:      0.00 A2
    
The results are all in the unit Ångström-squared. 

@section parameters Changing parameters

If higher precision is needed, the command

    $ freesasa -n 1000 1ubq.pdb

specifies that the calculation should use 1000 test points instead of
the default 100. The command

    $ freesasa --lee-richards -n 200 --probe-radius 1.2 --n-threads 4 1ubq.pdb

instead calculates the SASA using Lee & Richards algorithm with 200
slices per atom, a probe radius of 1.2 Å, using 4 parallel threads to
speed things up.

If the user wants to use their own atomic radii the command 

    $ freesasa --config-file <file> 1ubq.pdb

Reads a configuration from a file and uses it to assign atomic
radii. The program will halt if it encounters atoms in the PDB input
that are not present in the configuration. See @ref Config-file for
instructions how to write a configuration.

@section Output Other output types

To calculate the SASA of each residue in
the sequence, or each residue type, the commands

    $ freesasa --foreach-residue --no-log 1ubq.pdb
    SEQ: A    1 MET   55.99
    SEQ: A    2 GLN   72.16
    SEQ: A    3 ILE    0.00
    ...

and

    $ freesasa --foreach-residue-type --no-log 1ubq.pdb
    RES: ALA    122.87
    RES: ARG    531.77
    RES: ASN    160.57
    ...

to stdout respectively (`--no-log` suppresses the standard log
message). 

The command-line interface can also be used as a PDB filter:

    $ cat 1ubq.pdb | freesasa --no-log --print-as-B-values 
    ATOM      1  N   MET A   1      27.340  24.430   2.614  1.55 15.31
    ATOM      2  CA  MET A   1      26.266  25.413   2.842  2.00 20.34
    ATOM      3  C   MET A   1      26.913  26.639   3.531  1.55  0.00
    ...

The output is PDB-file where the temperature factors have been replaced by
SASA values, and occupancy numbers by the radius of each atom:

Only the atoms and models used in the calculation will be present in
the output (see @ref Input for how to modify this).

To generate all three results at the same time and write them to
separate files, run

    $ freesasa --residue-file=1ubq.seq --residue-type-file=1ubq.res --B-value-file=1ubq.b 1ubq.pdb

@section CLI-select Selecting groups of atoms

The option `--select` can be used to define groups of atoms whose
integrated SASA we are interested in. It uses a subset of the Pymol
`select` command syntax, see @ref Selection for full
documentation. The following example shows how to calculate the sum of
exposed surface areas of all aromatic residues and of ASP and ASN

    $ freesasa --select "aromatic, resn phe+tyr+trp+his+pro" --select "asx, resn asp+asn" resn 1ubq.pdb
    ...
    Selections:
    freesasa: warning: Found no matches to resn 'TRP', typo?
    aromatic:    348.32 A2
    asx:    549.25 A2

This command adds a 'Selection:' section at the end of the
output. This particular protein did not have any TRP residues, hence
the warning (written to stderr). The warning can be supressed with the
flag `-w`.

@section Input PDB input

@subsection Hetatom-hydrogen Including extra atoms

The user can ask to include hydrogen atoms and HETATM entries in the
calculation using the options `--hydrogen` and `--hetatm`. The default
radius of a hydrogen atom is 0, so including hydrogens will only make
sense if the user has also provided a config-file to specify a
hydrogen radius. By default the program tries to guess the radius of a
HETATM entry, but will halt if the element is not recognized. 

@subsection Chains-models Separating and joining chains and models

If a PDB file has several chains and/or models, by default all chains
of the first model are used, and the rest of the file is ignored. This
behavior can be modified using the following options 

  - `--join-models`: Joins all models in the input into one large
    structure. Useful for biological assembly files were different
    locations of the same chain in the oligomer are represented by
    different MODEL entries.

  - `--separate-models`: Calculate SASA separately for each model in
    the input. Useful when the same file contains several
    conformations of the same molecule.

  - `--separate-chains`: Calculate SASA separately for each chain in
    the input. Can be joined with `--separate-models` to calculate
    SASA of each chain in each model.

  - `--chain-groups`: Define groups of chains that should be treated
    as separate entities, can be repeated. If we for example have a
    tetramer with chains ABCD and want to know how much surface was
    buried when the dimers AB and CD were joined, we can use the option
    `--chain-groups=AB+CD`. The output will contain the SASA for the
    full molecule and one entry for each of the pairs of chains AB and
    CD. Can not be combined with `--separate-chains`.

@page API FreeSASA API

@section Basic-API Basics

The API is found in the header [freesasa.h](freesasa_8h.html) and is
the only header installed by `make install`. The other source-files
and headers in the repository are for internal use, but are also
thoroughly documented in the source itself.

To calculate the SASA of a structure, there are two options:

1. Initialize a structure from a PDB-file, calculate a radius for each atom 
    and then run the calculation.

2. Provide an array of cartesian coordinates and an array containing
   the radii of the corresponding atoms to freesasa_calc_coord().

@subsection API-PDB Calculate SASA for a PDB file

The following explains how to use FreeSASA to calculate the SASA of a
fictive PDB file (1abc.pdb). Possible errors are ignored for
brevity. Default parameters are used at every step, the section @ref
Customizing explains how to configure the calculations.

@subsubsection API-Read-PDB Open PDB file

Begin by opening a file and read a structure from it. The second
argument, 0, indicates use default when selecting atoms, i.e. ignore
Hydrogens and HETATMs, and only include the first MODEL if there are
several.

~~~{.c}
    FILE *fp = fopen("1abc.pdb");
    freesasa_structure *structure = freesasa_structure_from_pdb(fp,0);
~~~

@subsubsection API-Radii Calculate atomic radii

Calculate the atomic radii of that structure using the default classifier.
The argument `NULL` here means use the default classifier. A custom classifier
can be passed as a pointer to freesasa_classifier.

~~~{.c}
    double *radii = freesasa_structure_radius(structure,NULL);
~~~

@subsubsection API-Calc Perform calculation

Calculate SASA using the structure and radii and print the total area. The
argument `NULL` means use default parameters. User-defined parameters can be
defined by passing a pointer to freesasa_parameters.

~~~{.c}
    freesasa_result *result = freesasa_calc_structure(structure,radii,NULL);
    printf("Total area: %f A2\n",result->total);
~~~

@subsubsection API-Classes Get area of classes of atoms

Calculate area of classes (Polar/Apolar/..) and print their
values. The argument `NULL` again means use the default classifier. The
type freesasa_strvp stores pairs of strings and values, as
demonstrated below.

~~~{.c}
    freesasa_strvp *class_area = freesasa_result_classify(result,structure,NULL);
    for (int i = 0; i < class_area->n; ++i)
        printf("%s: %f A2\n",class_area->string[i],
               class_area->value[i]);
~~~

@subsubsection API-Select Get area of custom groups of atoms

Groups of atoms can be defined using freesasa_select_area(), which
takes a selection definition uses a subset of the Pymol select syntax
(@ref Selection).

~~~{.c}
    double area;
    char name[FREESASA_MAX_SELECTION_NAME+1];
    freesasa_select_area("aromatic, resn phe+tyr+trp+his+pro",
                         name,&area,structure,result);
~~~

The last statement stores the string "aromatic" in `name` and the SASA
integrated over the atoms defined by the selection.

@subsection Coordinates

If users wish to supply their own coordinates and radii, these are
accepted as arrays of doubles passed to the function
freesasa_calc_coord(). The coordinate-array should have size 3*n with
coordinates in the order `x1,y1,z1,x2,y2,z2,...,xn,yn,zn`.

@subsection Error-reporting 

Errors due to user or system errors, such as malformatted
config-files, I/O errors or memory allocation errors are reported
through return values, either ::FREESASA_FAIL or ::FREESASA_WARN, or
by NULL pointers, depending on the context. See the documentation for
the individual functions.

Errors that are attributable to programmers using the library, such as
passing null pointers are checked by asserts.

@subsection Thread-safety 

The only global state the library stores is the verbosity level (set
by freesasa_set_verbosity()). It should be clear from the
documentation when the other functions have side effects such as
memory allocation and I/O, and thread-safety should generally not be
an issue (to the extent that your C library has threadsafe I/O and
dynamic memory allocation). The SASA calculation itself can be
parallelized by passing a ::freesasa_parameters struct with
::freesasa_parameters.n_threads set to a value > 1 to
freesasa_calc_structure() or freesasa_calc_coord(). This only gives a
significant effect on performance for large proteins, and because not
all steps are parallelized it is not worth it to go beyond 2 threads.

@section Customizing Customizing behavior

The types ::freesasa_parameters and ::freesasa_classifier can be used
to change the parameters of the calculations. Users who wish to use
the defaults can pass NULL wherever pointers to these are requested.

@subsection Parameters Parameters

Changing parameters is done by passing a ::freesasa_parameters object
with the desired values. It can be initialized to default by

~~~{.c}
freesasa_parameters param = freesasa_default_parameters;
~~~

To allow the user to only change the parameters that are non-default.

The following call would run a Lee & Richards calculation with 200
slices per atom

~~~{.c}
param.alg = FREESASA_LEE_RICHARDS;
param.lee_richards_n_slices = 200;
freesasa_result *result = freesasa_calc_structure(structure,radii,param);
~~~

@subsection Classification Specifying atomic radii and classes

The type ::freesasa_classifier has function pointers to functions that
take residue and atom names as argument (pairs such as "ALA"," CA "),
and returns a radius or a class (polar, apolar, etc). The classifier
can be passed to freesasa_structure_radius() to generate an array of
atomic radii, which can then be used to calculate the SASA of the
structure. It can also be used in freesasa_result_classify() to get
the SASA integrated over the different classes of atoms, i.e. the SASA
of all polar atoms, etc.

Users of the API can provide their own classification by writing their
own functions and providing them via a ::freesasa_classifier object. A
classifier-configuration can also be read from a file using
freesasa_classifier_from_file() (see @ref Config-file).

The default classifier is available as a global const variable
::freesasa_default_classifier. This uses the classes and radii,
defined in the paper by Ooi et al.  ([PNAS 1987, 84:
3086](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC304812/)) for the
standard amino acids and also for some capping groups (ACE/NH2) if
HETATM fields are included when the PDB input is read. For other
residues such as Nucleic Acids or nonstandard amino acids, or
unrecognized HETATM entries the VdW radius of the element is
used. Warnings are emitted in the latter case. If the element can't be
determined or is unknown, a negative radius is returned (this will
make freesasa_structure_radius() return a NULL pointer).

@page Config-file Classifier configuration files

The configuration files read by freesasa_classifier_from_file() should
have two sections: `types:` and `atoms:`. A few example configurations 
are available in the directory `share/`.

The types-section defines what types of atoms are available
(aliphatic, aromatic, hydroxyl, ...), what the radius of that type is
and what class a type belongs to (polar, apolar, ...). The types are
just shorthands to associate an atom with a given combination of class
and radius. The user is free to define as many types and classes as
necessary.

The atoms-section consists of triplets of residue-name, atom-name (as
in the corresponding PDB entries) and type. A prototype file would be

~~~
types:
C_ALIPHATIC 2.00 apolar
C_AROMATIC  1.75 apolar
N 1.55 polar

atoms:
ANY N  N             
ANY CA C_ALIPHATIC
ANY CB C_ALIPHATIC

ARG CG C_ALIPHATIC

PRO CB C_AROMATIC  # overrides ANY CB
~~~

The residue type `ANY` can be used for atoms that are the same in all
or most residues (such as backbone atoms). If there is an exception
for a given amino acid this can be overridden as is shown for `PRO CB`
in the example.

@page Selection Selection syntax

FreeSASA uses a subset of the Pymol select commands to give users an
easy way of summing up the SASA of groups of atoms. This is done by
the function freesasa_select_area() in the C API, selectArea() in the
Python interface and the option `--select` for the command line
tool. All commands are case insensitive. A basic selection has a
selection name, a property selector and a list of arguments

    <selection-name>, <selector> <list>

For example

    aromatic, resn phe+tyr+trp+his+pro

Several selectors can be joined using boolean logic and parentheses,

    <selection-name>, (<s1> <l1>) and not (<s2> <l2> or <s3> <l3>)

where s1, s2 and s3 are selectors and l1, l2 and l3 are lists. The
operator `and` has precedence over `or`, so the second parentheses is
necessary but not the first, in the example above.

The following property selectors are supported

- `resn` Residue names like "ala", "arg", "du", etc 
- `resi` Residue index (positive integers)
- `chain` Chain labels (single characters)
- `name` Atom names, such as "ca", "c", "oxt", etc
- `symbol` Element symbols, such as "C", "O", "Se", "Fe", etc.

A list of residues can be selected using

    resn ala+val+leu+ile+met

and similarly for the other four selectors. In addition `resi` and
`chain` support ranges

    resi 1-10
    resi 1-10+20-30+35
    chain A+C-E

Combining ranges with plus signs, as in the two last lines, is not
allowed in Pymol but supported by FreeSASA.

If a selection list contains elements not found in the molecule that
is analyzed, a warning is printed and that part of the list does not
contribute to the selection. Not finding an a list element can be
because it specifies a residue that does not exist in the particular
molecule, or because of typos. The selector does not keep a list of
valid elements, residue names, etc.

@page Python Python bindings

If Python is enabled using 
    
    $ ./configure --enable-python-bindings

Cython is used to build Python bindings for FreeSASA, and `make
install` will install them.

Below follow some illustrations of how to use the package, see the
@ref freesasa "package documentation" for details.

@section Python-basics Basic calculations

Using defaults everywhere a simple calculation can be carried out as
follows (assuming the PDB structure 1UBQ is available)

~~~{.py} 
import freesasa

structure = freesasa.Structure("1ubq.pdb")
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result,structure)

print "Total : %.2f A2" % result.totalArea()
for key in area_classes:
    print key, ": %.2f A2" % area_classes[key]
~~~

Which would give the following output

    Total : 4779.51 A2
    Polar : 2236.93 A2
    Apolar : 2542.58 A2


The following does a high precision L&R calculation

~~~{.py}
result = freesasa.calc(structure,
                       freesasa.Parameters({'algorithm' : freesasa.LeeRichards,
                                            'n-slices' : 100}))
~~~

@section Python-classification Customizing atom classification

This uses the NACCESS parameters (the file 'naccess.config' is
available in the share/ directory of the repository).

~~~{.py}
classifier = freesasa.Classifier("naccess.config")
structure = freesasa.Structure("1ubq.pdb",classifier) 
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result,structure,classifier)
~~~

Classification can be customized also by extending the Classifier
interface. The code below is an illustration of a classifier that
classes Nitrogens separately, and assigns radii based on element only
(and crudely).

~~~{.py}
import re

class DerivedClassifier(Classifier):
    def classify(self,residueName,atomName):
        if re.match('\s*N',atomName):
            return 'Nitrogen'
        return 'Not-nitrogen'

    def radius(self,residueName,atomName):
        if re.match('\s*N',atomName): # Nitrogen 
            return 1.6
        if re.match('\s*C',atomName): # Carbon
            return 1.7
        if re.match('\s*O',atomName): # Oxygen
            return 1.4
        if re.match('\s*S',atomName): # Sulfur
            return 1.8    
    return 0;                     # everything else (Hydrogen, etc)

classifier = DerivedClassifier()

# use the DerivedClassifier to calculate atom radii (will give same result as the default)
structure = freesasa.Structure("1ubq.pdb",classifier)
result = freesasa.calc(structure)

# use the DerivedClassifier to classify atoms
area_classes = freesasa.classifyResults(result,structure,classifier)
~~~

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