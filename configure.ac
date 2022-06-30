#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([FreeSASA], [2.1.2])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
m4_include([m4/ax_compare_version.m4])
m4_include([m4/ax_cxx_compile_stdcxx.m4])
AC_DEFINE([REPORTBUG], ["Report bugs to <https://github.com/mittinatten/freesasa/issues>"], [Bug report string])
AC_DEFINE([HOMEPAGE], ["<http://freesasa.github.io>"], [Project homepage])
AC_DEFINE([FREESASA_XMLNS], ["http://freesasa.github.io/"], [XML namespace])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([dl], [dlsym])

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([inttypes.h libintl.h malloc.h stddef.h stdlib.h string.h strings.h sys/time.h unistd.h dlfcn.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_RESTRICT
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# use C99
AC_PROG_CC_C99
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset mkdir sqrt strchr strdup strerror strncasecmp getopt_long getline])

# C++ 14
AX_CXX_COMPILE_STDCXX([14])

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile doc/Doxyfile doc/man/freesasa.1 tests/Makefile share/Makefile])
AC_CONFIG_FILES([tests/test-cli], [chmod +x tests/test-cli])

# Optionally disable thread support
AC_ARG_ENABLE([threads],
  AS_HELP_STRING([--disable-threads],
    [Build without support for multiple threads]))

if test "x$enable_threads" != "xno"; then
  AC_CHECK_LIB([pthread], [pthread_create])
  AC_CHECK_HEADER([pthread.h])
fi
if test "$ac_cv_header_pthread_h" != yes || \
   test "$ac_cv_lib_pthread_pthread_create" != yes; then
  AC_DEFINE([USE_THREADS], [0], [])
  AM_CONDITIONAL([USE_THREADS], false)
  echo "Building without thread support"
else
  AC_DEFINE([USE_THREADS], [1], [Define if threads should be used.])
  AM_CONDITIONAL([USE_THREADS], true)
fi

# disable XML
AC_ARG_ENABLE([xml],
  AS_HELP_STRING([--disable-xml],
     [Build without support for XML output]))

AC_DEFINE([USE_XML], [0], [Define if XML should be included.])
AM_CONDITIONAL([USE_XML], false)

if test "x$enable_xml" != "xno" ; then
  AC_CHECK_PROG([XMLLINT],[xmllint],[xmllint],[])
  PKG_CHECK_MODULES([libxml2], [libxml-2.0],
    [AC_DEFINE([USE_XML], [1])
     AC_SUBST([USE_XML], [yes])
     AM_CONDITIONAL([USE_XML], true)],
    [AC_MSG_ERROR([Library libxml2 not found, install library or build without (using --disable-xml).])])
else
    AC_MSG_NOTICE([Building without support for XML output.])
fi

# Disable json-c
AC_ARG_ENABLE([json],
  AS_HELP_STRING([--enable-json],
    [Build without support for JSON output]))

AC_DEFINE([USE_JSON], [0], [Define if JSON should be included])
AM_CONDITIONAL([USE_JSON], false)

if test "x$enable_json" != "xno" ; then
  AC_CHECK_LIB([json-c], [json_object_get],
     [AC_CHECK_HEADER([json-c/json_object.h],
        [AC_DEFINE([USE_JSON], [1])
         AC_SUBST([USE_JSON], [yes])
         AM_CONDITIONAL([USE_JSON], true)])])
  AM_COND_IF([USE_JSON],[],
    [AC_MSG_ERROR([Library JSON-C not found, install library or build without (using --disable-json).])])
  AC_CHECK_PROG([JSONLINT],[jsonlint],[jsonlint],[])
else
  AC_MSG_NOTICE([Building without support for JSON output.])
fi

# Enable parser generation with Flex/Bison
AC_ARG_ENABLE([parser-generator],
  [AS_HELP_STRING([--enable-parser-generator],
    [Regenerate parser/lexer using Flex and Bison (for developers)])])

AM_CONDITIONAL([GENERATE_PARSER],false)

if test "x$enable_parser_generator" = "xyes" ; then
  AC_PROG_LEX
  AC_CHECK_PROGS([YACC], [bison byacc yacc])
  AM_CONDITIONAL([GENERATE_PARSER], [test -n $LEX && test "x$LEX" != "x:" && test -n $YACC])
  AM_COND_IF([GENERATE_PARSER],[],
    [AC_MSG_ERROR([Lexer and/or parser generator missing. Install Flex/lex and/or Bison/Yacc or configure without flag --enable-parser-generator.])])
fi

# Enable Doxygen
AC_ARG_ENABLE([doxygen],
  [AS_HELP_STRING([--enable-doxygen],
    [Build Doxygen html-documentation])])

AM_CONDITIONAL([HAVE_DOXYGEN],false)

if test "x$enable_doxygen" = "xyes" ; then
  AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)
  AM_CONDITIONAL([HAVE_DOXYGEN],[test -n "$DOXYGEN"])
  AM_COND_IF([HAVE_DOXYGEN],[build_doc=yes],
  [AC_MSG_ERROR([Doxygen not available, can not build documentation. Install library or configure without --enable-doxygen.])])

  if test -n "$DOXYGEN"; then
     DOXYGEN_V=`doxygen --version`
     DOXYGEN_REQ=1.8.8
     AX_COMPARE_VERSION([$DOXYGEN_V],[lt],[$DOXYGEN_REQ],
       [echo "Warning: Doxygen should preferrably be of version $DOXYGEN_REQ or higher (found $DOXYGEN_V)"])
  fi
fi

AM_CONDITIONAL([BUILD_DOC],[test -n "$build_doc"])

# Enable testing
AC_ARG_ENABLE([check],
  [AS_HELP_STRING([--enable-check],
    [Enable the Check unit testing framework])])

AM_CONDITIONAL([USE_CHECK],false)
AC_DEFINE([USE_CHECK],[0],[Build with Check unit testing framework])
AM_CONDITIONAL([USE_SUBUNIT], false)

if test "x$enable_check" = "xyes" ; then
  AC_CHECK_LIB(rt, setitimer) # necessary on some platforms
  AC_CHECK_LIB([subunit], [subunit_test_start]) # necessary on some platforms
  AC_CHECK_LIB([check],[suite_create],
    [],
    [AC_MSG_ERROR[Library Check not available, install libary or configure without flag --enable-check]])
  AC_DEFINE([USE_CHECK], [1])
  AM_CONDITIONAL([USE_CHECK], true)
fi

# Enable gcov
AC_ARG_ENABLE([gcov],
  [AS_HELP_STRING([--enable-gcov],
    [Compile with Gcov to measure coverage])])
AM_CONDITIONAL([COND_GCOV],false)

if test "x$enable_gcov" = "xyes" ; then
  AC_CHECK_LIB([gcov], [__gcov_init], [], [AC_MSG_ERROR([Missing gcov library])])
  AM_CONDITIONAL([COND_GCOV],true)
fi

# Disable test for CLI
AC_ARG_ENABLE([CLI-tests],
  AS_HELP_STRING([--disable-CLI-tests],
    [Disable testing of command-line interface. These take time and
     disabling them can help speed up development when only unit
     tests are needed.]))
if test "x$enable_CLI_tests" != "xno" ; then
   AM_CONDITIONAL([RUN_CLI_TESTS],true)
else
   AM_CONDITIONAL([RUN_CLI_TESTS],false)
fi

AC_OUTPUT
