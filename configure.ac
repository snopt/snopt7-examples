#-------------------------------------------------------------------------------
AC_PREREQ(2.69)
AC_INIT([SNOPT-EXAMPLES],[7.7],[optimizers@ccom.ucsd.edu],[snopt])

AC_CONFIG_SRCDIR([fortran/snmaina.f])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT([`pwd`])

user_FCFLAGS=${FCFLAGS}
user_CFLAGS=${CFLAGS}
user_CPPFLAGS=${CPPFLAGS}
user_CXXFLAGS=${CXXFLAGS}
user_LDFLAGS=${LDFLAGS}

# Libtool
LT_INIT
AC_SUBST([LIBTOOL_DEPS])


# Check compilers.
AC_PROG_FC([gfortran nagfor ifort])
AC_PROG_CC
AC_PROG_CXX

AC_FC_LIBRARY_LDFLAGS
AC_SUBST([FCLIBS])

AC_SEARCH_LIBS([pow], [m])
AC_SUBST([LIBS])

AC_FC_MODULE_OUTPUT_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODOUT="-mdir ";;
    esac
  ])
AC_FC_MODULE_FLAG([],
  [ case $FC in
      *nagfor* ) FC_MODINC="-mdir ";;
    esac
  ])
AC_SUBST(FC_MODOUT)
AC_SUBST(FC_MODINC)

AC_PROG_MAKE_SET
AC_SUBST([SET_MAKE])


dnl ----------------------------------------------------------------------------

AC_ARG_ENABLE([openmp],
            [AS_HELP_STRING([--enable-openmp],
              [enable OpenMP support (disabled by default)])],
            [],
            [OPENMP_FC_FLAGS=""; enable_openmp=no;])

if test x"${enable_openmp}" = xyes; then
  AC_LANG_PUSH(Fortran)
  AC_OPENMP
  AC_LANG_POP(Fortran)
fi
AC_SUBST([OPENMP_FCFLAGS])
AC_SUBST([enable_openmp])

dnl ----------------------------------------------------------------------------

# SNOPT library
AC_MSG_CHECKING(for snopt)
AC_ARG_WITH([snopt],
	    [AS_HELP_STRING([--with-snopt="-L/path/to/snopt/lib -lsnopt7"],
	     [location of SNOPT library])],
	    [ snopt_lib="${withval}" ],
	    [ snopt_lib=${libdir} ] )
AC_ARG_WITH([snopt-cpp],
	    [AS_HELP_STRING([--with-snopt-cpp="-L/path/to/snopt/lib -lsnopt7_cpp"],
	     [location of SNOPT C++ library])],
	    [ snopt_lib="${withval}"; with_cpp=yes; ],
	    [] )
AC_MSG_RESULT(${snopt_lib})
AC_SUBST([snopt_lib])
AC_SUBST([with_cpp])

dnl ----------------------------------------------------------------------------

# BLAS library.
AC_ARG_WITH([blas],
	    [AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])],
	    [with_blas=${withval} ], [ blas=no ] )
case $with_blas in
	yes | "" | no) blas=no ;;
	-* | */* | *.a | *.so | *.so.* | *.o)
	  blas=yes; BLAS_LIBS="$with_blas" ;;
        matlab )
	  blas=yes; BLAS_LIBS="matlab" ;;
	*)
	  blas=yes; BLAS_LIBS="-l$with_blas" ;;
esac
AC_SUBST([blas])
AC_SUBST([BLAS_LIBS])

dnl ----------------------------------------------------------------------------

AC_ARG_WITH([64],
            [AS_HELP_STRING( [--with-64],
                              [64-bit integer support (default is 32-bit)]) ],
 	    [ prec=64;
              AC_MSG_NOTICE([64-bit integers enabled]) ], [])

AC_ARG_WITH([128],
            [AS_HELP_STRING( [--with-128],
                              [quad-precision support (64-bit ints, 128-bit reals)]) ],
 	    [ prec=128;
              AC_MSG_NOTICE([quad-precision enabled (64-bit ints, 128-bit reals)]) ], [])

AC_ARG_WITH([debug],
            [AS_HELP_STRING([--with-debug],
                            [Compile with debug flags])],
 	    [DBGFLAG="-g"; OPTFLAG="";],
	    [OPTFLAG="-O"; DBGFLAG="";])

dnl ----------------------------------------------------------------------------

FCFLAGS="${user_FCFLAGS}"
CFLAGS="${user_CFLAGS}"
CXXFLAGS="${user_CXXFLAGS}"
LDFLAGS="${user_LDFLAGS}"

if test x"${prec}" == x; then
   prec=32
fi

if test x"${prec}" = x64; then
   PC_FLAGS="-DD64"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS="-fdefault-integer-8"
   ;;
   *nagfor*)
     PF_FLAGS="-i8"
   ;;
   *ifort*)
     PF_FLAGS="-i8"
   ;;
   esac
elif test x"${prec}" = x128; then
   PC_FLAGS=" -DD128"

   case "${FC}" in
   *gfortran*)
     PF_FLAGS=" -fdefault-integer-8"
   ;;
   *nagfor* )
     PF_FLAGS=" -i8"
   ;;
   *ifort* )
     PF_FLAGS=" -i8"
   ;;
   esac
fi

case "${FC}" in
  *nagfor* )
    PF_FLAGS+=" -w"
  ;;
esac

if test -z "${FCFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    FCFLAGS="${OPTFLAG} ${PF_FLAGS}"
  else
    FCFLAGS="${DBGFLAG} ${PF_FLAGS}"
  fi
else
  FCFLAGS="${FCFLAGS} ${PF_FLAGS}"
fi

if test -z "${CFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CFLAGS="${OPTFLAG} ${PC_FLAGS}"
  else
    CFLAGS="${DBGFLAG} ${PC_FLAGS}"
  fi
else
  CFLAGS="${CFLAGS} ${PC_FLAGS}"
fi

if test -z "${CXXFLAGS}"; then
  if test -z "${DBGFLAG}"; then
    CXXFLAGS="${OPTFLAG}"
  else
    CXXFLAGS="${DBGFLAG}"
  fi
fi

AC_SUBST([prec])
AC_SUBST([FCFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([CXXFLAGS])

dnl ----------------------------------------------------------------------------
AC_SUBST([enable_shared])
AC_SUBST([enable_static])

AC_OUTPUT([GNUmakefile])

echo ""
echo "SNOPT Examples:"
if test x"${with_cpp}" == xyes; then
echo "  SNOPT C++ library               ${snopt_lib}"
else
echo "  SNOPT library                   ${snopt_lib}"
fi
echo ""
echo "To build all examples             'make' or 'make all'"
echo "To build Fortran examples         'make fortran'"
echo "To build C examples               'make c'"
if test x"${with_cpp}" == xyes; then
echo "To build C++ examples             'make cpp'"
fi
echo "To check the build                'make check'"
echo ""
