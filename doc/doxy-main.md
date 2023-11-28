# FreeSASA

These pages document the

- @ref CLI
- @ref API "FreeSASA C API"
- @ref Config-file
- @ref Selection
- @ref Geometry

The [FreeSASA Python Module](https://github.com/freesasa/freesasa-python)
is documented [elsewhere](http://freesasa.github.io/python/).

The library is released under the [MIT license](license.md).

Installation instructions can be found in the [README](README.md) file.

@page Python Python Module

This page has moved to http://freesasa.github.io/python/.

@page CLI Command-line Interface

Building FreeSASA creates the binary `freesasa`, which is installed by
`make install`. Calling

    $ freesasa -h

displays a help message listing all options. The following text
explains how to use most of them.

@section CLI-Default Run using defaults

In the following we will use the RNA/protein complex PDB structure
3WBM as an example. It has four protein chains A, B, C and D, and two
RNA strands X and Y. To run a simple SASA calculation using default
parameters, simply type:

    $ freesasa 3wbm.pdb

This generates the following output

    ## FreeSASA 2.1.0 ##

    PARAMETERS
    algorithm    : Lee & Richards
    probe-radius : 1.400
    threads      : 2
    slices       : 20

    INPUT
    source  : 3wbm.pdb
    chains  : ABCDXY
    atoms   : 3714

    RESULTS (A^2)
    Total   :   25190.77
    Apolar  :   11552.38
    Polar   :   13638.39
    CHAIN A :    3785.49
    CHAIN B :    4342.33
    CHAIN C :    3961.12
    CHAIN D :    4904.30
    CHAIN X :    4156.46
    CHAIN Y :    4041.08

The results are all in the unit Ångström-squared.

For mmCIF input files, the flag `--cif` can be used

    $ freesasa --cif 3wbm.cif

which will generate the exact same output.

@section parameters Changing parameters

If higher precision is needed, the command

    $ freesasa -n 100 3wbm.pdb

specifies that the calculation should use 100 slices per atom instead of
the default 20. The command

    $ freesasa --shrake-rupley -n 200 --probe-radius 1.2 --n-threads 4 3wbm.pdb

instead calculates the SASA using Shrake & Rupley's algorithm with 200
test points, a probe radius of 1.2 Å, using 4 parallel threads to
speed things up.

If the user wants to use their own atomic radii the command

    $ freesasa --config-file <file> 3wbm.pdb

Reads a configuration from a file and uses it to assign atomic
radii. The program will halt if it encounters atoms in the PDB input
that are not present in the configuration. See @ref Config-file for
instructions how to write a configuration.

To use the atomic radii from NACCESS call

    $ freesasa --radii=naccess 3wbm.pdb

Another way to specify a custom set of atomic radii is to store them as
occupancies in the input PDB file

    $ freesasa --radius-from-occupancy 3wbm.pdb

This option allows the user to first use the option `--format=pdb` (see @ref CLI-PDB) to
write generate a PDB file with the radii used in the calculation,
modify the radii of individual atoms in that file, and then recalculate
the SASA with these modified radii.

@section Output Output formats

In addition to the standard output format above FreeSASA can export
the results as @ref CLI-JSON, @ref CLI-XML, @ref CLI-PDB, @ref CLI-CIF-OUTPUT,
@ref CLI-RSA, @ref CLI-RES and @ref CLI-SEQ using the option
`--format`. The level of detail of JSON and XML output can be
controlled with the option `--output-depth=<depth>` which takes the
values `atom`, `residue`, `chain` and `structure`. If `atom` is
chosen, SASA values are shown for all levels of the structure,
including individual atoms. With `chain`, only structure and chain
SASA values are printed (this is the default).

The output can include relative SASA values for each residues. To
calculate these a reference SASA value is needed, calculated using the
same atomic radii. At the moment such values are only available for
the ProtOr and NACCESS radii (selected using the option `--radii`), if
other radii are used relative SASA will be excluded (in RSA output all
REL columns will have the value 'N/A').

The reference SASA values for residue X are calculated from Ala-X-Ala
peptides in a stretched out configuration. The reference
configurations are supplied for reference in the directory
`rsa`. Since these are not always the most exposed possible
configuration, and because bond lengths and bond angles vary, the
relative SASA values will sometimes be larger than 100 %. At the
moment there is no interface to supply user-defined reference values.

@subsection CLI-JSON JSON

The command

    $ freesasa --format=xml --output-depth=residue 3wbm.pdb

generates the following

```{.json}
{
  "source":"FreeSASA 2.1.0",
  "length-unit":"Ångström",
  "results":[
    {
      "input":"3wbm.pdb",
      "classifier":"ProtOr",
      "parameters":{
        "algorithm":"Lee & Richards",
        "probe-radius":1.3999999999999999,
        "resolution":20
      },
      "structures":[
        {
          "chain-labels":"ABCDXY",
          "area":{
            "total":25190.768387067546,
            "polar":13638.391677017404,
            "apolar":11552.376710050148,
            "main-chain":3337.1622502425053,
            "side-chain":21853.606136825045
          },
          "chains":[
            {
              "label":"A",
              "n-residues":86,
              "area":{
                "total":3785.4864049452635,
                "polar":1733.8560208488598,
                "apolar":2051.6303840964056,
                "main-chain":723.34358684348558,
                "side-chain":3062.1428181017791
              }
              "residues":[
                {
                  "name":"THR",
                  "number":"5",
                  "area":{
                    "total":138.48216994006549,
                    "polar":56.887951514571867,
                    "apolar":81.594218425493622,
                    "main-chain":38.898190013033592,
                    "side-chain":99.583979927031905
                  },
                  "relative-area":{
                    "total":104.05152148175331,
                    "polar":113.98106895325961,
                    "apolar":98.093554250413092,
                    "main-chain":96.330336832673567,
                    "side-chain":107.414496739329
                  },
                  "n-atoms":7
               },

            ...

            },

            ...

          ]
        }
      ]
    }
  ]
}
```

Where ellipsis indicates the remaining residues and chains.

@subsection CLI-XML XML

The command

    $ freesasa --format=xml 3wbm.pdb

Generates the following

```{.xml}
<?xml version="1.0" encoding="UTF-8"?>
<results xmlns="http://freesasa.github.io/" source="FreeSASA 2.1.0" lengthUnit="&#xC5;ngstr&#xF6;m">
  <result classifier="ProtOr" input="3wbm.pdb">
    <parameters algorithm="Lee &amp; Richards" probeRadius="1.400000" resolution="20"/>
    <structure chains="ABCDXY">
      <area total="25190.768" polar="13638.392" apolar="11552.377" mainChain="3337.162" sideChain="21853.606"/>
      <chain label="A" nResidues="86">
        <area total="3785.486" polar="1733.856" apolar="2051.630" mainChain="723.344" sideChain="3062.143"/>
      </chain>
      <chain label="B" nResidues="84">
        <area total="4342.334" polar="1957.114" apolar="2385.220" mainChain="853.707" sideChain="3488.627"/>
      </chain>
      <chain label="C" nResidues="86">
        <area total="3961.119" polar="1838.724" apolar="2122.395" mainChain="782.652" sideChain="3178.468"/>
      </chain>
      <chain label="D" nResidues="89">
        <area total="4904.298" polar="2332.306" apolar="2571.991" mainChain="977.459" sideChain="3926.838"/>
      </chain>
      <chain label="X" nResidues="25">
        <area total="4156.455" polar="2919.576" apolar="1236.879" mainChain="0.000" sideChain="4156.455"/>
      </chain>
      <chain label="Y" nResidues="25">
        <area total="4041.076" polar="2856.815" apolar="1184.261" mainChain="0.000" sideChain="4041.076"/>
      </chain>
    </structure>
  </result>
</results>
```

@subsection CLI-PDB PDB

The command-line interface can also be used as a PDB filter:

    $ cat 3wbm.pdb | freesasa --format=pdb
    REMARK 999 This PDB file was generated by FreeSASA 2.1.0
    REMARK 999 In the ATOM records temperature factors have been
    REMARK 999 replaced by the SASA of the atom, and the occupancy
    REMARK 999 by the radius used in the calculation.
    MODEL        1
    ATOM      1  N   THR A   5     -19.727  29.259  13.573  1.64  9.44
    ATOM      2  CA  THR A   5     -19.209  28.356  14.602  1.88  5.01
    ATOM      3  C   THR A   5     -18.747  26.968  14.116  1.61  0.40
    ...

The output is a PDB-file where the temperature factors have been
replaced by SASA values (last column), and occupancy numbers by the
radius of each atom (second to last column).

Only the atoms and models used in the calculation will be present in
the output (see @ref Input for how to modify this).

@subsection CLI-CIF-OUTPUT CIF

The option `--output=cif` can be used if the input was in the mmCIF format.
It takes the original input and adds two columns to the `atom_site`s

    loop_
    _atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_alt_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_entity_id
    _atom_site.label_seq_id
    _atom_site.pdbx_PDB_ins_code
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.occupancy
    _atom_site.B_iso_or_equiv
    _atom_site.pdbx_formal_charge
    _atom_site.auth_seq_id
    _atom_site.auth_comp_id
    _atom_site.auth_asym_id
    _atom_site.auth_atom_id
    _atom_site.pdbx_PDB_model_num
    _atom_site.FreeSASA_value
    _atom_site.FreeSASA_radius
    ATOM 1 N N . MET A 1 1 ? 27.340 24.430 2.614 1.00 9.67 ? 1 MET A N 1 19.515388 1.640000
    ...

The parameters used

    _freeSASA_parameters.version  2.1.0
    _freeSASA_parameters.algorithm 'Lee & Richards'
    _freeSASA_parameters.probe-radius 1.400000
    _freeSASA_parameters.slices 20

Results at structure and chain level

    loop_
    _freeSASA_results.model
    _freeSASA_results.chains
    _freeSASA_results.atoms
    _freeSASA_results.type
    _freeSASA_results.surface_area
    1 A 602 Total 4804.055641
    1 A 602 Apolar 2299.838339
    1 A 602 Polar 2504.217302
    1 A 602 'CHAIN A' 4804.055641

And residue level results, including relative values, if available

    loop_
    _freeSASA_rsa.asym_id
    _freeSASA_rsa.seq_id
    _freeSASA_rsa.comp_id
    _freeSASA_rsa.abs_total
    _freeSASA_rsa.rel_total
    _freeSASA_rsa.abs_side_chain
    _freeSASA_rsa.rel_side_chain
    _freeSASA_rsa.abs_main_chain
    _freeSASA_rsa.rel_main_chain
    _freeSASA_rsa.abs_apolar
    _freeSASA_rsa.rel_apolar
    _freeSASA_rsa.abs_polar
    _freeSASA_rsa.rel_polar
    A 1 MET 54.393508 28.168570 19.085046 12.630739 35.308462 84.067766 28.483151 24.216248 25.910357 34.327447
    ...

Only the atoms and models used in the calculation will be present in
the output (see @ref Input for how to modify this).

@subsection CLI-RES SASA of each residue type

Calculate the SASA of each residue type:

    $ freesasa --format=res 3wbm.pdb
    # Residue types in 3wbm.pdb
    RES ALA :     251.57
    RES ARG :    2868.98
    RES ASN :    1218.87
    ...
    RES A :    1581.57
    RES C :    2967.12
    RES G :    1955.16
    RES U :    1693.68

@subsection CLI-SEQ SASA of each residue

Calculate the SASA of each residue in the sequence:

    $ freesasa --format=seq 3wbm.pdb
    # Residues in 3wbm.pdb
    SEQ A    5 THR :  138.48
    SEQ A    6 PRO :   25.53
    SEQ A    7 THR :   99.42
    ...

@subsection CLI-RSA RSA

The CLI can also produce output similar to the RSA format from
NACCESS. This format includes both absolute SASA values (ABS) and
relative ones (REL) compared to a precalculated reference max
value. The main difference between FreeSASA's RSA output format and
that of NACCESS, is that FreeSASA will print the value "N/A" where
NACCESS prints "-99.9", when the a reference value for REL is missing.

    $ freesasa --format=rsa 3wbm.pdb
    REM  FreeSASA 2.1.0
    REM  Absolute and relative SASAs for 3wbm.pdb
    REM  Atomic radii and reference values for relative SASA: ProtOr
    REM  Chains: ABCDXY
    REM  Algorithm: Lee & Richards
    REM  Probe-radius: 1.40
    REM  Slices: 20
    REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
    REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
    RES THR A   5   138.48 104.1  99.58 107.4  38.90  96.3  81.59  98.1  56.89 114.0
    RES PRO A   6    25.53  19.3  11.31  11.0  14.23  47.7  21.67  18.7   3.86  23.9
    ...
    RES GLY A  15     0.64   0.9   0.00   N/A   0.64   0.9   0.00   0.0   0.64   2.0
    ...
    RES   U Y  23   165.16   N/A 165.16   N/A   0.00   N/A  52.01   N/A 113.15   N/A
    RES   C Y  24   165.01   N/A 165.01   N/A   0.00   N/A  46.24   N/A 118.77   N/A
    RES   C Y  25   262.46   N/A 262.46   N/A   0.00   N/A  85.93   N/A 176.52   N/A
    END  Absolute sums over single chains surface
    CHAIN  1 A     3785.5       3062.1        723.3       2051.6       1733.9
    CHAIN  2 B     4342.3       3488.6        853.7       2385.2       1957.1
    CHAIN  3 C     3961.1       3178.5        782.7       2122.4       1838.7
    CHAIN  4 D     4904.3       3926.8        977.5       2572.0       2332.3
    CHAIN  5 X     4156.5       4156.5          0.0       1236.9       2919.6
    CHAIN  6 Y     4041.1       4041.1          0.0       1184.3       2856.8
    END  Absolute sums over all chains
    TOTAL         25190.8      21853.6       3337.2      11552.4      13638.4

Note that each `RES` is a single residue, not a residue type as above
(i.e. has the same meaning as `SEQ` above). This unfortunate confusion
of labels is due to RSA support being added much later than the other
options. Fixing it now would break the interface, and will thus
earliest be dealt with in the next major release.

@subsubsection RSA-naccess Using the NACCESS configuration

The reference values for the NACCESS configuration in FreeSASA are not
exactly the same as those that ship with NACCESS, but have been
calculated from using the same Ala-X-Ala tripeptides (kindly donated by
the creators of NACCESS). By default NACCESS defines amino acid CA
atoms as side-chain, whereas FreeSASA defines them as main-chain. The
flag `-b` defines CA as main-chain in NACCESS.

Calling

    $ freesasa 3wbm.pdb -n 20--format=rsa --radii=naccess

will give an RSA file where the ABS columns should be identical to
NACCESS with flag `-b`. All REL values will differ slightly because of
round-off errors in the references. The difference will be large
for side-/main-chain REL values, since the FreeSASA references were
calculated with CA as main-chain. NACCESS also gives different results
for the nucleic acid main-chain and side-chain. FreeSASA defines the
(deoxy)ribose and phosphate groups as main-chain and the base as
side-chain.

Future versions might allow specifying what atoms are main-chain
through the configuration file, and letting the user set their own
reference values.

@section CLI-select Selecting groups of atoms

The option `--select` can be used to define groups of atoms whose
integrated SASA we are interested in. It uses a subset of the Pymol
`select` command syntax, see @ref Selection for full
documentation. The following example shows how to calculate the sum of
exposed surface areas of all aromatic residues and of the four chains
A, B, C and D (just the sum of the areas above).

    $ freesasa --select "aromatic, resn phe+tyr+trp+his+pro" --select "abcd, chain A+B+C+D" 3wbm.pdb
    ...
    SELECTIONS
    freesasa: warning: Found no matches to resn 'TRP', typo?
    freesasa: warning: Found no matches to resn 'HIS', typo?
    aromatic :    1196.45
    abcd :   16993.24

The lines shown above are appended to the regular output. This
particular protein did not have any TRP or HIS residues, hence the
warnings (written to stderr). The warnings can be supressed with the
flag `-w`.

@section Chain-groups Analyzing groups of chains

Calculating the SASA of a given chain or group of chains separately
from the rest of the structure, can be useful for measuring how buried
a chain is in a given structure. The option `--chain-groups` can be
used to do such a separate calculation, calling

    $ freesasa --chain-groups=ABCD+XY 3wbm.pdb

produces the regular output for the structure 3WBM, but in addition it
runs a separate calculation for the chains A, B, C and D as though X
and Y aren't in the structure, and vice versa:

    PARAMETERS
    algorithm    : Lee & Richards
    probe-radius : 1.400
    threads      : 2
    slices       : 20


    ####################

    INPUT
    source  : 3wbm.pdb
    chains  : ABCDXY
    atoms   : 3714

    RESULTS (A^2)
    Total   :   25190.77
    Apolar  :   11552.38
    Polar   :   13638.39
    CHAIN A :    3785.49
    CHAIN B :    4342.33
    CHAIN C :    3961.12
    CHAIN D :    4904.30
    CHAIN X :    4156.46
    CHAIN Y :    4041.08


    ####################

    INPUT
    source  : 3wbm.pdb
    chains  : ABCD
    atoms   : 2664

    RESULTS (A^2)
    Total   :   18202.78
    Apolar  :    9799.46
    Polar   :    8403.32
    CHAIN A :    4243.12
    CHAIN B :    4595.18
    CHAIN C :    4427.11
    CHAIN D :    4937.38


    ####################

    INPUT
    source  : 3wbm.pdb
    chains  : XY
    atoms   : 1050

    RESULTS (A^2)
    Total   :    9396.28
    Apolar  :    2743.09
    Polar   :    6653.19
    CHAIN X :    4714.45
    CHAIN Y :    4681.83

@section Input PDB input

@subsection Hetatom-hydrogen Including extra atoms

The user can ask to include hydrogen atoms and HETATM entries in the
calculation using the options `--hydrogen` and `--hetatm`. In both
cases adding unknown atoms will emit a warning for each atom. This can
either be amended by using the flag `-w` to suppress warnings, or by
using a custom classifier so that they are recognized (see @ref
Config-file).

@subsection Halt-skip Skipping unknown atoms

By default FreeSASA guesses the element of an unknown atom and uses
that elements VdW radius. If this fails the radius is set to 0 (and
hence the atom will not contribute to the calculated area). Users can
request to either skip unknown atoms completely (i.e. no guessing) or
to halt when unknown atoms are found and exit with an error. This is
done with the option `--unknown` which takes one of the three
arguments `skip`, `halt` or `guess` (default). Whenever an unknown
atom is skipped or its radius is guessed a warning is printed to
stderr.

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

- `--chain-groups`: see @ref Chain-groups

@page API FreeSASA API

@section Basic-API Basics

The API is found in the header [freesasa.h](freesasa_8h.html). The
other source-files and headers in the repository are for internal use,
and are not presented here, but are documented in the source
itself. The file [example.c](example_8c_source.html) contains a simple
program that illustrates how to use the API to read a PDB file from
`stdin` and calculate and print the SASA.

To calculate the SASA of a structure, there are two main options:

1. Initialize a structure from a PDB-file, using either the default
   classifier or a custom one to determine the radius of each atom,
   and then run the calculation.

2. Provide an array of cartesian coordinates and an array containing
   the radii of the corresponding atoms to freesasa_calc_coord().

@subsection API-PDB Calculate SASA for a PDB file

The following explains how to use FreeSASA to calculate the SASA of a
fictive PDB file (1abc.pdb). At each step one or more error checks
should have been done, but these are ignored here for brevity. See
the documentation of each function to see what errors can occur.
Default parameters are used at every step, the section @ref
Customizing explains how to configure the calculations.

@subsubsection API-Read-PDB Open PDB file

The function freesasa_structure_from_pdb() reads the atom
coordinates from a PDB file and assigns a radius to each atom. The
third argument can be used to pass options for how to read the PDB
file.

```{.c}
    FILE *fp = fopen("1abc.pdb");
    const freesasa_classifier *classifier = &freesasa_default_classifier;
    freesasa_structure *structure = freesasa_structure_from_pdb(fp, classifier, 0);
```

@subsubsection API-Calc Perform calculation and get total SASA

Next we use freesasa_calc_structure() to calculate SASA using the
structure we just generated, and then print the total area. The argument
`NULL` means use default freesasa_parameters.

```{.c}
    freesasa_result *result = freesasa_calc_structure(structure, NULL);
    printf("Total area: %f A2\n",result->total);
```

@subsubsection API-Classes Get polar and apolar area

We are commonly interested in the polar and apolar areas of a
molecule, this can be calculated by freesasa_result_classes(). To
get other classes of atoms we can either define our own classifier, or
use freesasa_select_area() defined in the next section. The return
type ::freesasa_nodearea is a struct contains the total area and the
area of all apolar and polar atoms, and main-chain and side-chain
atoms.

```{.c}
    freesasa_nodearea area = freesasa_result_classes(structure, result);
    printf("Total      : %f A2\n", area.total);
    printf("Apolar     : %f A2\n", area.apolar);
    printf("Polar      : %f A2\n", area.polar);
    printf("Main-chain : %f A2\n", area.main_chain);
    printf("Side-chain : %f A2\n", area.side_chain);
```

@see @ref Classification

@subsubsection API-Select Get area of custom groups of atoms

Groups of atoms can be defined using freesasa_selection_new(), which
takes a selection definition uses a subset of the Pymol select syntax

```{.c}
    freesasa_selection *selection =
        freesasa_selection_new("aromatic, resn phe+tyr+trp+his+pro",
                               structure, result);
    printf("Area of selection '%s': %f A2\n",
           freesasa_selection_name(selection), freesasa_selection_area(selection);
```

@see @ref Selection

@subsubsection structure-node Navigating the results as a tree

In addition to the flat array of results in ::freesasa_result, and
the global values returned by freesasa_result_classes(), FreeSASA
has an interface for navigating the results as a tree. The leaf nodes
are individual atoms, and there are parent nodes at the residue,
chain, and structure levels. The function freesasa_calc_tree() does
a SASA calculation and returns the root node of such a tree. (If one
already has a ::freesasa_result the function freesasa_tree_init()
can be used instead). Each node stores a ::freesasa_nodearea for the
sum of all atoms belonging to the node. The tree can be traversed with
freesasa_node_children(), freesasa_node_parent() and
freesasa_node_next(), and the area, type and name using
freesasa_node_area(), freesasa_node_type() and
freesasa_node_name(). Additionally there are special properties for
each level of the tree.

@see node

@subsubsection export-tree Exporting to RSA, JSON and XML

The tree structure can also be exported to an RSA, JSON or XML file
using freesasa_tree_export(). The RSA format is fixed, but the user
can select which levels of the tree to include in JSON and XML. The
following illustrates how one would generate a tree and export it to
XML, including nodes for the whole structure, chains and residues (but
excluding individual atoms).

```{.c}
    freesasa_node *tree = freesasa_calc_tree(structure,
                                             &freesasa_default_parameters,
                                             &freesasa_default_classifier);
    FILE *file = fopen("output.xml", "w");
    freesasa_tree_export(file, tree, FREESASA_XML | FREESASA_OUTPUT_RESIDUE);
    fclose(file);
    freesasa_node_free(tree);
```

@subsection Coordinates

If users wish to supply their own coordinates and radii, these are
accepted as arrays of doubles passed to the function
freesasa_calc_coord(). The coordinate-array should have size 3\*n with
coordinates in the order `x1,y1,z1,x2,y2,z2,...,xn,yn,zn`.

```{.c}
    double coord[] = {1.0, /* x */
                      2.0, /* y */
                      3.0  /* z */ };
    double radius[] = {2.0};
    int n_atoms = 1;
    freesasa_result *result = freesasa_calc_coord(coord, radius, n_atoms, NULL);
```

@subsection Error-handling

The principle for error handling is that unpredictable errors should
not cause a crash, but rather allow the user to exit gracefully or
make another attempt. Therefore, errors due to user or system
failures, such as faulty parameters, malformatted config-files, I/O
errors or out of memory errors, are reported through return values,
either ::FREESASA_FAIL or ::FREESASA_WARN, or by `NULL` pointers,
depending on the context (see the documentation for the individual
functions).

Errors that are attributable to programmers using the library, such as
passing null pointers where not allowed, are checked by asserts.

@subsection Thread-safety

The only global state the library stores is the verbosity level (set
by freesasa_set_verbosity()) and the pointer to the error-log
(defaults to `stderr`, can be changed by freesasa_set_err_out()).

It should be clear from the documentation when the other functions
have side effects such as memory allocation and I/O, and thread-safety
should generally not be an issue (to the extent that your C library
has threadsafe I/O and dynamic memory allocation). The SASA
calculation itself can be parallelized by using a
::freesasa_parameters struct with ::freesasa_parameters.n_threads
\> 1 (default is 2) where appropriate. This only gives a significant
effect on performance for large proteins or at high precision, and
because not all steps are parallelized it is usually not worth it to
go beyond 2 threads.

@section Customizing Customizing behavior

The types ::freesasa_parameters and ::freesasa_classifier can be
used to change the parameters of the calculations. Users who wish to
use the defaults can pass `NULL` wherever pointers to these are
requested.

@subsection Parameters Parameters

Calculation parameters can be stored in a ::freesasa_parameters
object. It can be initialized to default by

```{.c}
freesasa_parameters param = freesasa_default_parameters;
```

The following code would run a high precision Shrake & Rupley
calculation with 10000 test points on the provided structure.

```{.c}
freesasa_parameters param = freesasa_default_parameters;
param.alg = FREESASA_SHRAKE_RUPLEY;
param.shrake_rupley_n_points = 10000;
freesasa_result *result = freesasa_calc_structure(structure, param);
```

@subsection Classification Specifying atomic radii and classes

Classifiers are used to determine which atoms are polar or apolar, and
to specify atomic radii. In addition the three standard classifiers
(see below) have reference values for the maximum areas of the 20
standard amino acids which can be used to calculate relative areas of
residues, as in the RSA output.

The default classifier is available through the const variable
::freesasa*default_classifier. This uses the \_ProtOr* radii, defined
in the paper by Tsai et
al. ([JMB 1999, 290: 253](http://www.ncbi.nlm.nih.gov/pubmed/10388571))
for the standard amino acids (20 regular plus SEC, PYL, ASX and GLX),
for some capping groups (ACE/NH2) and the standard nucleic acids. If
the element can't be determined or is unknown, a zero radius is
assigned. It classes all carbons as _apolar_ and all other known atoms
as _polar_.

Early versions of FreeSASA used the atomic radii by Ooi et
al. ([PNAS 1987, 84: 3086-3090](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC304812/)),
this classifier is still available through ::freesasa_oons_classifier.

Users can provide their own classifiers through @ref Config-file. At
the moment these do not allow the user to specify reference values to
calculate relative SASA values for RSA output.

The default behavior of freesasa_structure_from_pdb(),
freesasa_structure_array(), freesasa_structure_add_atom() and
freesasa_structure_add_atom_wopt() is to first try the provided
classifier and then guess the radius if necessary (emitting warnings
if this is done, uses VdW radii defined by [Mantina et al. J Phys Chem
2009, 113:5806](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/)).

See the documentation for these functions for what parameters to use
to change the default behavior.

@page Config-file Classifier configuration files

The configuration files read by freesasa_classifier_from_file() or the
command-line option `-c` should have two sections: `types:` and
`atoms:`, and optionally the section `name:`.

The types-section defines what types of atoms are available
(aliphatic, aromatic, hydroxyl, ...), what the radius of that type is
and what class a type belongs to ('polar' or 'apolar', case
insensitive). The types are just shorthands to associate an atom with
a given combination of class and radius. The user is free to define as
many types and classes as necessary.

The atoms-section consists of triplets of residue-name, atom-name (as
in the corresponding PDB entries) and type. A prototype file would be

```
name: myclassifier  # tag and value must be on the same line (optional)

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
```

The residue type `ANY` can be used for atoms that are the same in all
or most residues (such as backbone atoms). If there is an exception
for a given amino acid this can be overridden as is shown for `PRO CB`
in the example.

A few example configurations are available in the directory
[share/](https://github.com/mittinatten/freesasa/tree/master/share). The
configuration-file
[protor.config](https://github.com/mittinatten/freesasa/tree/master/share/protor.config)
is a copy of the default classifier, and can be used to add extra
atoms that need to be classified, while keeping the defaults for the
standard residues (also see the file
[scripts/chemcomp2config.pl](https://github.com/mittinatten/freesasa/tree/master/scripts/)
for instructions on how to generate configurations for new chemical
components semi-automatically). If something common is missing in the
default classifier, [create an
issue](https://github.com/mittinatten/freesasa/issues) on Github so
that it can be added.

FreeSASA also ships with some configuration-files that mimic other
popular programs, such as
[NACCESS](https://github.com/mittinatten/freesasa/tree/master/share/naccess.config)
and
[DSSP](https://github.com/mittinatten/freesasa/tree/master/share/dssp.config).

The static classifiers in the API were generated using
[scripts/config2c.pl](https://github.com/mittinatten/freesasa/tree/master/scripts/)
to convert the correspoding configurations in `share` to C code.

@page Selection Selection syntax

FreeSASA uses a subset of the Pymol select commands to give users an
easy way of summing up the SASA of groups of atoms. This is done by
the function freesasa_selection_new() in the C API,
freesasa.selectArea() in the Python interface and the option
`--select` for the command line tool. All commands are case
insensitive. A basic selection has a selection name, a property
selector and a list of arguments

    <selection-name>, <selector> <list>

For example

    aromatic, resn phe+tyr+trp+his+pro

Several selectors can be joined using boolean logic and parentheses,

    <selection-name>, (<s1> <l1>) and not (<s2> <l2> or <s3> <l3>)

where s1, s2 and s3 are selectors and l1, l2 and l3 are lists. The
operator `and` has precedence over `or`, so the second parentheses is
necessary but not the first, in the example above. The selection name
can include letters, numbers and underscores. The name can't be longer
than ::FREESASA_MAX_SELECTION_NAME characters.

The following property selectors are supported

- `resn` Residue names like "ala", "arg", "du", etc
- `resi` Residue index (positive or negative integers)
- `chain` Chain labels (single characters)
- `name` Atom names, such as "ca", "c", "oxt", etc
- `symbol` Element symbols, such as "C", "O", "Se", "Fe", etc.

A list of residues can be selected using

    resn ala+val+leu+ile+met

and similarly for the other four selectors. In addition `resi` and
`chain` support ranges

    resi 1-10             (residues 1 to 10)
    resi -10              (residues indices < 10)
    resi 10-              (residues indices > 10)
    resi 1-10+20-30+35-   (residues 1 to 10, 20 to 30 and above 35)
    resi \-20-\-15+\-10-5 (residues -20 to -15 and -10 to 5)
    chain A+C-E           (chains A and C to E, no open intervals allowed here)

Combining ranges with plus signs, as in the three last lines, is not
allowed in Pymol but supported by FreeSASA.

If a selection list contains elements not found in the molecule that
is analyzed, a warning is printed and that part of the list does not
contribute to the selection. Not finding a list element can be because
it specifies a residue that does not exist in the particular molecule,
or because of typos. The selector does not keep a list of valid
elements, residue names, etc.

@page Geometry Geometry of Lee & Richards' algorithm

This page explains the geometry of the calculations in L&R
and can be used to understand the source code. As far as possible the
code uses similar notation to the formulas here.

We will use the following notation: An atom \f$i\f$ has a van der
Waals radius \f$r_i\f$, the rolling sphere (or _probe_) has radius
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

\f[ \alpha = 2\arccos \bigl[({R_i^\prime}^2{\,}

- d_{ij}^2 - {R_{j}^\prime}^2{\,})/(2R_i^\prime d_{ij})\bigr]. \f]

If the middle point of this arc on the circle is at an angle
\f$\beta\f$, the arc spans the interval
\f$[\beta-\alpha/2,\beta+\alpha/2]\f$. By adding up these arcs and
taking into account any overlap between them we get the total buried
angle \f$\gamma\f$ in this slices. The exposed arc angle for this atom
and slice is thus \f$2\pi-\gamma\f$ and the total SASA of that atom

\f[ A_i =R_i \delta \!\! \sum_{s\in\text{slices}} \!\!
\left[2\pi-\gamma_s\right]\,. \f]

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
