#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([FreeSASA], [0.3.3])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/freesasa.c])
AC_CONFIG_HEADERS([config.h])
m4_include([m4/ax_compare_version.m4])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [pthread_create])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h pthread.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_RESTRICT

# use C99
AC_PROG_CC_C99

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday sqrt strerror])
AC_CHECK_FUNCS([mkdir])

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile doc/Doxyfile
                 doc/fig/Makefile tests/Makefile bindings/Makefile])

# Python-bindings
AC_ARG_WITH([python],
  AS_HELP_STRING([--with-python=<python-binary>], [Specify which python to use.]),
    [python_bin=$withval],[python_bin=python])

AC_ARG_ENABLE([python-binding],
  [AS_HELP_STRING([--enable-python-binding],
    [build Python bindings using Cython])],
  [enable_python_binding=yes], [enable_python_binding=no])

if test x$enable_python_binding = xyes ; then
  AC_CHECK_PROG(PYTHON, python, $python_bin)
  AC_CHECK_PROG(CYTHON, cython, cython) 
  AM_CONDITIONAL([HAVE_CYTHON], [test -n "$CYTHON" && test -n "$PYTHON"])
  AM_COND_IF([HAVE_CYTHON],[],[echo "Cython not available, will not build Python bindings"])
else
  AM_CONDITIONAL([HAVE_CYTHON],FALSE)
fi

# Enable LaTeX
AC_ARG_ENABLE([latex],
 [AS_HELP_STRING([--enable-latex],
   [build documentation from latex sources])],
 [enable_latex=yes], [enable_latex=no])

if test x$enable_latex = xyes ; then
  AC_CHECK_PROG(LATEX, latex, latex)
  AC_CHECK_PROG(DVIPS, dvips, dvips)
  AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
  AC_CHECK_PROG(EPSTOPDF, epstopdf, epstopdf)
  AM_CONDITIONAL([HAVE_LATEX], [test -n "$LATEX" && test -n "$DVIPS" \
                                && test -n "$EPSTOPDF" && test -n "$PDFLATEX"])
  AM_COND_IF([HAVE_LATEX],[],
    [echo "LaTeX not available, option '--enable-latex' ignored"])
  AM_COND_IF([HAVE_LATEX],[build_doc=yes])
else
  AM_CONDITIONAL([HAVE_LATEX],FALSE)
fi

# Enable Doxygen
AC_ARG_ENABLE([doxygen],
  [AS_HELP_STRING([--enable-doxygen],
    [build doxygen html-documentation])],
  [enable_doxygen=yes],
  [enable_doxygen=no])
if test x$enable_doxygen = xyes ; then
  AC_CHECK_PROG(DOXYGEN, doxygen, doxygen)
  AM_CONDITIONAL([HAVE_DOXYGEN],[test -n "$DOXYGEN"])
  AM_COND_IF([HAVE_DOXYGEN],[build_doc=yes])

  if test -n "$DOXYGEN"; then
     DOXYGEN_V=`doxygen --version`
     DOXYGEN_REQ=1.8.8
     AX_COMPARE_VERSION([$DOXYGEN_V],[lt],[$DOXYGEN_REQ],
       [echo "Warning: Doxygen should preferrably be of version $DOXYGEN_REQ or higher (found $DOXYGEN_V)"])
  fi
else
  AM_CONDITIONAL([HAVE_DOXYGEN],FALSE)
fi

AM_CONDITIONAL([BUILD_DOC],[test -n "$build_doc"])

# Enable testing
AC_ARG_ENABLE([check],
  [AS_HELP_STRING([--enable-check],
    [enable the Check unit testing framework])],
    [enable_check=yes],
    [enable_check=no])

if test x$enable_check = xyes ; then
  AC_CHECK_LIB([check],[suite_create],
    [have_check=yes],
    [echo "Check not availabe, option '--enable-check' ignored"])
  AM_CONDITIONAL([USE_CHECK],[test x$have_check = xyes])
else
  AM_CONDITIONAL([USE_CHECK],FALSE)
fi

# Enable gcov
AC_ARG_ENABLE([gcov],
  [AS_HELP_STRING([--enable-gcov],
    [compile with Gcov to measure coverage])],
  [enable_gcov=yes],[enable_gcov=no])
AM_CONDITIONAL([COND_GCOV],[test "x$enable_gcov" = xyes])


AC_OUTPUT
