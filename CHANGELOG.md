# Changelog
FreeSASA uses semantic versioning. Changelog added for versions 2.x

## [unreleased]
* Classifier interface changed.
  * Pointer is now opaque and classification done via freesasa_classifier_class() instead of classifier->sasa_class().
  * Classifiers only classify polar/apolar/unknown. Arbitrary classes no longer allowed.
  * The static classifiers in freesasa.h have reference values for all residue types, allowing calculation of relative SASA for for example RSA output. At the moment no support for defining these values in classifier configuration file (TODO?).
  * Classifier configuration files now have new name field
* New output formats
  * XML using libmxl2 (optional)
  * JSON using JSON-C (optional)
  * RSA (was present already in 1.1, but interface now consolidated)
* New node interface. Results in tree form can be generated using freesasa_calc_tree() or freesasa_tree_init(). The features was added to facilitate generation of XML and JSON output, but can be used directly to analyze results as well.
* structure.c was refactored, hopefully code is slightly more transparent now
* general bug cleaning, some minor memory leaks removed
* interface for classifying using string-value pairs is deprecated
* logging using structure and results now deprecated in favor of using tree interface
* Travis CI and Coveralls now used to check commits

