/* config.h — compile-time feature flags for FreeSASA.
   Values here are defaults; CMake overrides them via -D flags at build time. */

#ifndef FREESASA_CONFIG_H
#define FREESASA_CONFIG_H

/* Package metadata */
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.1.3"
#endif
#ifndef PACKAGE_STRING
#define PACKAGE_STRING  "FreeSASA 2.1.3"
#endif
#ifndef PACKAGE_NAME
#define PACKAGE_NAME    "freesasa"
#endif

/* Feature flags — set by CMake at build time */
#ifndef USE_OPENMP
#define USE_OPENMP 1
#endif

/* USE_THREADS is kept as an alias for USE_OPENMP for backward-compat with
   existing tests and Python bindings that check #if USE_THREADS */
#ifndef USE_THREADS
#define USE_THREADS 1
#endif

/* Optional output formats (off by default) */
#ifndef USE_JSON
#define USE_JSON 0
#endif

#ifndef USE_XML
#define USE_XML 0
#endif

/* Check unit test framework (enabled when building tests) */
#ifndef USE_CHECK
#define USE_CHECK 0
#endif

/* Memory-error tests require malloc interposition — only safe in test builds */
#ifndef INCLUDE_MEMERR_TESTS
#define INCLUDE_MEMERR_TESTS 0
#endif

/* Data paths for tests — overridden at compile time */
#ifndef DATADIR
#define DATADIR "tests/data/"
#endif

#ifndef SHAREDIR
#define SHAREDIR "share/freesasa/classifications/"
#endif


/* Standard POSIX headers available */
#define HAVE_CONFIG_H 1
#define STDC_HEADERS  1

/* URL strings for --help and --version output */
#define REPORTBUG "Report bugs to <https://github.com/mittinatten/freesasa/issues>"
#define HOMEPAGE  "<http://freesasa.github.io>"

/* Misc defines */
#define FREESASA_XMLNS "http://freesasa.github.io/"

#endif /* FREESASA_CONFIG_H */
