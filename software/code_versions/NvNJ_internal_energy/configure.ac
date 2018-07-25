AC_PREREQ(2.61)                                                                                                     
AC_INIT(dvm, 0.1.0, peter.clarke29@gmail.com)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([subdir-objects])

# Release versioning
GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=1
GENERIC_MICRO_VERSION=0
GENERIC_LIBRARY_NAME=dvm

#shared library versioning
AC_SUBST(GENERIC_LIBRARY_VERSION)
PACKAGE=$GENERIC_LIBRARY_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])

# Checks for programs.
AC_PROG_CC

AC_PROG_FC
AC_PROG_LIBTOOL
AC_FC_SRCEXT(f90)

#-----------------------
# Check for debug mode
#-----------------------
AC_MSG_CHECKING(whether to enable debugging)
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                     [default=$debug_default]],, enable_debug=$debug_default)

if test "x$enable_debug" = "xyes"; then
   FCFLAGS="$FCFLAGS -g -DDEBUG"
   AM_FCFLAGS="$FCFLAGS -g -DDEBUG"
   AC_MSG_RESULT(yes)
else
   FCFLAGS="$FCFLAGS -O3"
   AC_MSG_RESULT(no)
fi


#----------------
# libGRVY Config
#----------------
AX_PATH_GRVY_NEW([0.29.1],[no])

# Check for GSL
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))

# Query configuration environment
AX_SUMMARIZE_ENV

AC_OUTPUT(    Makefile 
          src/Makefile
          test/Makefile )

# Final summary
AX_SUMMARIZE_CONFIG
