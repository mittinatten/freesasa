# Changelog

FreeSASA uses semantic versioning. Changelog added for versions 2.x

## Upcoming release

- Breaking: Support 3-letter chain labels (auth_sym_id) in CIF-files. Default and RSA output formats 
  have changes in whitespace as a consequence, to give space to the larger labels.

## 2.1.2

- Fix error where CLI options weren't parsed properly on ARM processors.
- Fix error in man pages.

## 2.1.1

- Fix error where compiling without multithread support still set the default number of threads to 2.
- Fix issue when stripping quotes from atom names in CIF-files, which sometimes caused memory errors.

## 2.1.0

The new features in this release were created on initiative by and in collaboration
with [Danny Diaz](https://github.com/danny305).

### Added

- Support for mmCIF input with the CLI option `--cif`.
  - The CLI was ported to C++ to allow using Gemmi for CIF import.
  - Gemmi is imported as a git submodule, see README for details.
  - A test runner was added to verify that CIF and PDB input files
    give the same result.
  - The C API does not support CIF for now (this would require conversion to C++).
- Add output option `--format=cif` that can be used when input is mmCIF.
  See documentation for an example.
- Added VdW-radii for all elements.
  - Using values from Gemmi when missing in paper by Mantina et al.
  - Correct Mg radius from 1.74 to 1.73 Å.

### Fixed

- Fix bug in JSON output where relative SASA for amino acids without sidechan
  were written as `NaN`, which is not valid JSON. These values are now
  simply skipped.
- Fix bug where elements with H or D as second letter, such as CD, were classified
  as hydrogens.
- Change some type definitions to be more C++ friendly.

## 2.0.3

This version separates the Python bindings into a separate
[module](https://github.com/freesasa/freesasa-python).
To be able to release Windows binaries of the python module,
the code has been changed to be C89 compatible where necessary,
with a few macros to allow using `restrict` and `inline` when
compiled with a C99 compiler. There are no changes to the behavior
of the CLI or C API in this release.

## 2.0.2

- Relative SASA values are now calculated using the same reference
  configurations as NACCESS. This means relative SASA values will
  differ slightly from those of previous versions.
- Python bindings compatible with Python 3.
- Man page added
- Bugs fixed:
  - CLI option `--separate-models` now outputs all models in PDB output.
  - On some platforms compilation failed when libxml and/or json-c was
    not present. Now fixed.
  - Some memory allocations were not checked for failure in S&R
    calculations. These are now properly checked, and done more
    seldomly, which should improve performance slightly.
- Residue numbers now include the iCode field for insterted residues
- Compatibility with Microsoft C Compiler:
  - No variable length arrays
  - Some macros and keywords are redefined when necessary
  - Lexer does not depend on header `unistd.h`
  - GNU getline not used, replaced by fgets()

## 2.0.1

- Add function `calcCoord()` to Python interface, wraps C function
  `freesasa_calc_coord()`.
- Expand selection syntax:
  - Allow negative residue numbers `resi \-10` (escaped with
    backslash), and expression such as `resi \-10-\-5+\-3-5`
  - Allow open-ended residue ranges `resi -10`, `resi 10-` and
    combined as `resi -10+15-20+30-`
  - Allow primes in atom names (i.e. `O5'`)
  - Allow selection names to be numbers or start with a number

## 2.0

### Changed

#### General

- Classifier configuration files now have new name field

#### C API

- Classifier interface changed.
  - Pointer is now opaque and classification done via
    `freesasa_classifier_class()` instead of
    `classifier->sasa_class()`.
  - Classifiers only classify polar/apolar/unknown. Arbitrary classes
    no longer allowed.
  - The static classifiers in `freesasa.h` have reference values for
    all residue types, allowing calculation of relative SASA for for
    example RSA output. At the moment no support for defining these
    values in classifier configuration file (TODO?).
  - Interface for classifying using string-value pairs has been
    removed (type `freesasa_strvp` and functions
    `freesasa_result_classifiy()` and `freesasa_strvp_free()`). The
    new tree-interface should be used instead.
- Structure interface changed: A structure is initialized with a
  classifier and the relevant classifier values are stored directly in
  the structure. This assures that the structure and results are
  always analyzed with consistent parameters. One effect is that the
  function `freesasa_calc_structure()` no longer takes the atomic
  radii as an explicit parameters.
- The logging functions `freesasa_log()`, `freesasa_write_result()`,
  `freesasa_write_parameters()`, `freesasa_write_rsa()`,
  `freesasa_write_pdb()`, `freesasa_per_residue()`,
  `freesasa_per_residue_type()` and `freesasa_per_chain()` have been
  removed from the interface. The same functionality can be achieved
  using `freesasa_tree_export()`.

#### CLI

- CLI output format is controlled by the option `-f` or
  `--format`. Options specifying files for specific outpt types have
  been removed, and the separate options controlling output format
  have been deprecated (see below). An option `--deprecated` can be
  used to list these. This has also made the option `-l` or `--no-log`
  obsolete.

### Added

#### General

- New output formats
  - XML using libmxl2 (optional)
  - JSON using JSON-C (optional)
  - RSA (was present already in 1.1, but interface now consolidated)
- To build FreeSASA 2.0 the libraries libxml2 and JSON-C need to be
  installed. These can be disabled by `./configure --disable-json --disable-xml` to build a dependency-free version of the library.

#### C API

- New node interface. Results in tree form can be generated using
  `freesasa_calc_tree()` or `freesasa_tree_init()`. The feature was
  added to facilitate generation of XML and JSON output, but can be
  used directly to analyze results as well. The results of a
  calculation are to now best exported using `freesasa_tree_export()`.

### Fixed

- `structure.c` was refactored, hopefully code is slightly more
  transparent now
- general bug cleaning, some minor memory leaks removed.
- Memory error mocking is now more sophisticated: uses dlsym instead
  of macros, allowing more uniform test structure.

### Deprecated

#### CLI

- options `-B`, `-r`, `-R`, `--rsa` and `-l` are deprecated (use
  `-f` or `--format` instead)

#### C API

- Logging using structure and results now deprecated in favor of
  using tree interface
- Function `freesasa_select_area()` is deprecated
