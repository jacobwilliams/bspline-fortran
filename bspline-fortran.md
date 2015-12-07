project: bspline-fortran
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/bspline-fortran
summary: BSPLINE-FORTRAN -- Multidimensional B-Spline Interpolation of Data on a Regular Grid
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
source: true
graph: true
exclude: pyplot_module.f90
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran

Brief description
---------------

The library provides subroutines for 1D-6D interpolation using B-splines. The code is written in modern Fortran (i.e., Fortran 2003+).

License
---------------

The bspline-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).

See also
---------------

* This library includes the public domain DBSPLIN and DTENSBS code from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm) (CMLIB).

Keywords
---------------
* Bspline, spline, interpolation, data fitting, multivariate interpolation, multidimensional interpolation
