#!/bin/bash

#
#  Simple build script for bspline-fortran.
#
#  Requires: FoBiS.py and ROBODoc
#

PROJECTNAME='bspline-fortran'   # project name for robodoc
MODCODE='bspline_module.f90'    # module file name
LIBOUT='libbspline.a'           # name of library
DOCDIR='./doc/'                 # build directory for documentation
SRCDIR='./src/'                 # library source directory
TESTSRCDIR='./src/tests/'       # unit test source directory
BINDIR='./bin/'                 # build directory for unit tests
LIBDIR='./lib/'                 # build directory for library

#compiler flags:

FCOMPILER='gnu' #Set default compiler to gfortran
FCOMPILERFLAGS='-c -O2 -std=f2008'

#build the library:

FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./ -t ${MODCODE} -o ${LIBOUT} -mklib static -colors

# build the test program:

FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTSRCDIR} -dmod ./ -dobj ./ -colors -libs ${LIBDIR}${LIBOUT} --include ${LIBDIR}

# build the documentation:

if hash robodoc 2>/dev/null; then
	echo "Building documentation..."
	robodoc --rc ./robodoc.rc --src ${SRCDIR} --doc ${DOCDIR} --documenttitle ${PROJECTNAME}
else
	echo "ROBODoc not found! Cannot build documentation. ROBODoc can be obtained from: http://www.xs4all.nl/~rfsber/Robo/"
fi
