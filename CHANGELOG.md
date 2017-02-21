# Changelog
FreeSASA uses semantic versioning. Changelog added for versions 2.x

## [unreleased]

* Add function `calcCoord()` to Python interface, wraps C function
  `freesasa_calc_coord()`.
* Expand selection syntax:
  * Allow negative residue numbers `resi \-10` (escaped with
    backslash), and expression such as `resi \-10-\-5+\-3-5`
  * Allow open-ended residue ranges `resi -10`, `resi 10-` and
    combined as `resi -10+15-20+30-`
  * Allow primes in atom names (i.e. `O5'`)
  * Allow selection names to be numbers or start with a number

## 2.0

### Changed

#### General
* Classifier configuration files now have new name field

#### C API
* Classifier interface changed.
  * Pointer is now opaque and classification done via
    `freesasa_classifier_class()` instead of
    `classifier->sasa_class()`.
  * Classifiers only classify polar/apolar/unknown. Arbitrary classes
    no longer allowed.
  * The static classifiers in `freesasa.h` have reference values for
    all residue types, allowing calculation of relative SASA for for
    example RSA output. At the moment no support for defining these
    values in classifier configuration file (TODO?).
  * Interface for classifying using string-value pairs has been
    removed (type `freesasa_strvp` and functions
    `freesasa_result_classifiy()` and `freesasa_strvp_free()`). The
    new tree-interface should be used instead.
* Structure interface changed: A structure is initialized with a
  classifier and the relevant classifier values are stored directly in
  the structure. This assures that the structure and results are
  always analyzed with consistent parameters. One effect is that the
  function `freesasa_calc_structure()` no longer takes the atomic
  radii as an explicit parameters.
* The logging functions `freesasa_log()`, `freesasa_write_result()`,
  `freesasa_write_parameters()`, `freesasa_write_rsa()`,
  `freesasa_write_pdb()`, `freesasa_per_residue()`,
  `freesasa_per_residue_type()` and `freesasa_per_chain()` have been
  removed from the interface. The same functionality can be achieved
  using `freesasa_tree_export()`.

#### CLI
* CLI output format is controlled by the option `-f` or
  `--format`. Options specifying files for specific outpt types have
  been removed, and the separate options controlling output format
  have been deprecated (see below). An option `--deprecated` can be
  used to list these. This has also made the option `-l` or `--no-log`
  obsolete.

### Added

#### General
* New output formats
  * XML using libmxl2 (optional)
  * JSON using JSON-C (optional)
  * RSA (was present already in 1.1, but interface now consolidated)
* To build FreeSASA 2.0 the libraries libxml2 and JSON-C need to be
  installed. These can be disabled by `./configure --disable-json
  --disable-xml` to build a dependency-free version of the library.

#### C API
* New node interface. Results in tree form can be generated using
  `freesasa_calc_tree()` or `freesasa_tree_init()`. The feature was
  added to facilitate generation of XML and JSON output, but can be
  used directly to analyze results as well. The results of a
  calculation are to now best exported using `freesasa_tree_export()`.

### Fixed
* `structure.c` was refactored, hopefully code is slightly more
  transparent now
* general bug cleaning, some minor memory leaks removed.
* Memory error mocking is now more sophisticated: uses dlsym instead
  of macros, allowing more uniform test structure.

### Deprecated

#### CLI
* options `-B`, `-r`, `-R`, `--rsa` and `-l` are deprecated (use
  `-f` or `--format` instead)

#### C API
* Logging using structure and results now deprecated in favor of
  using tree interface
* Function `freesasa_select_area()` is deprecated
