project: bspline-fortran
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/bspline-fortran
summary: BSPLINE-FORTRAN -- Multidimensional B-Spline Interpolation of Data on a Regular Grid
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark: #
exclude: test.f90

Brief description
---------------

The library provides subroutines for 2D-6D interpolation using b-splines. The code is written in modern Fortran (i.e., Fortran 2003+).

License
--------

The bspline-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).

See also
---------------

* This code includes the public domain DBSPLIN and DTENSBS code from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm) (CMLIB).  The 2D and 3D routines are extensively refactored versions of the original routines from DTENSBS.  The 4D-6D routines are new, and are simply extensions of the same algorithm into higher dimensions.