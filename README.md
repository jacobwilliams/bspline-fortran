bspline-fortran [![GitHub release](https://img.shields.io/github/release/jacobwilliams/bspline-fortran.svg?style=plastic)](https://github.com/jacobwilliams/bspline-fortran/releases/latest) [![DOI](https://zenodo.org/badge/31299552.svg)](https://zenodo.org/badge/latestdoi/31299552)
============

Multidimensional B-Spline Interpolation of Data on a Regular Grid

# Status

[![Build Status](https://img.shields.io/travis/jacobwilliams/bspline-fortran/master.svg?style=plastic)](https://travis-ci.org/jacobwilliams/bspline-fortran)

# Brief description

The library provides subroutines for 1D-6D interpolation and extrapolation using B-splines. The code is written in modern Fortran (i.e., Fortran 2003+). There are two ways to use the module, via a basic subroutine interface and an object-oriented interface. Both are thread safe.

## Subroutine interface

The core routines for the subroutine interface are:

```Fortran

!f(x)
subroutine db1ink(x,nx,fcn,kx,iknot,tx,bcoef,iflag)
subroutine db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

!f(x,y)
subroutine db2ink(x,nx,y,ny,fcn,kx,ky,iknot,tx,ty,bcoef,iflag)
subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,w1,w0,extrap)

!f(x,y,z)
subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,iknot,tx,ty,tz,bcoef,iflag)
subroutine db3val(xval,yval,zval,idx,idy,idz,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef,f,iflag,inbvx,inbvy,inbvz,iloy,iloz,w2,w1,w0,extrap)

!f(x,y,z,q)
subroutine db4ink(x,nx,y,ny,z,nz,q,nq,fcn,kx,ky,kz,kq,iknot,tx,ty,tz,tq,bcoef,iflag)
subroutine db4val(xval,yval,zval,qval,idx,idy,idz,idq,tx,ty,tz,tq,nx,ny,nz,nq,kx,ky,kz,kq,bcoef,f,iflag,inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq,w3,w2,w1,w0,extrap)

!f(x,y,z,q,r)
subroutine db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn,kx,ky,kz,kq,kr,iknot,tx,ty,tz,tq,tr,bcoef,iflag)
subroutine db5val(xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,tx,ty,tz,tq,tr,nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,bcoef,f,iflag,inbvx,inbvy,inbvz,inbvq,inbvr,iloy,iloz,iloq,ilor,w4,w3,w2,w1,w0,extrap)

!f(x,y,z,q,r,s)
subroutine db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn,kx,ky,kz,kq,kr,ks,iknot,tx,ty,tz,tq,tr,ts,bcoef,iflag)
subroutine db6val(xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,tx,ty,tz,tq,tr,ts,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,bcoef,f,iflag,inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,iloy,iloz,iloq,ilor,ilos,w5,w4,w3,w2,w1,w0,extrap)
```

The ```ink``` routines compute the interpolant coefficients, and the ```val``` routines evalute the interpolant at the specified value of each coordinate. The 2D and 3D routines are extensively refactored versions of the original routines from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).  The others are new, and are simply extensions of the same algorithm into the other dimensions.

## Object-oriented interface

In addition to the main subroutines, an object-oriented interface is also provided. For example, for the 3D case:

```Fortran
type(bspline_3d) :: s
call s%initialize(x,y,z,fcn,kx,ky,kz,iflag,extrap)
call s%evaluate(xval,yval,zval,idx,idy,idz,f,iflag)
call s%destroy()
```
Which uses the default "not-a-knot" end conditions. You can also specify the knot vectors (in this case, `tx`, `ty`, and `tz`) manually during class initialization:

```Fortran
call s%initialize(x,y,z,fcn,kx,ky,kz,tx,ty,tz,iflag,extrap)
```

The various bspline classes can also be initialized using constructors, which have similar interfaces as the `initialize` methods. For example:

```Fortran
type(bspline_3d) :: s
s = bspline_3d(x,y,z,fcn,kx,ky,kz,iflag,extrap)
```

## Extrapolation

The library optionally supports extrapolation for points outside the range of the coefficients. This is disabled by default (in which case an error code is returned for points outside the bounds). To enable extrapolation, use the optional `extrap` input to the various `db*val` subroutines or the `initialize` methods from the object-oriented interface.

## Integration

The library also contains routines for computing definite integrals of bsplines. There are two methods (currently only for 1D):

* Basic version: `db1sqad` (`integral` in the object-oriented interface) -- Computes the integral on `(x1,x2)` of a b-spline by applying a 2, 6, or 10 point Gauss formula on subintervals of `(x1,x2)`. This is only valid for orders <= 20.
* More general version: `db1fqad` (`fintegral` in the object-oriented interface) -- Computes the integral on `(x1,x2)` of a product of a user-defined function `fun(x)` and the ith derivative of a b-spline with an adaptive 8-point Legendre-Gauss algorithm.

Note that extrapolation is not currently supported for these.

# Examples

See the [examples](https://github.com/jacobwilliams/bspline-fortran/tree/master/src/tests) for more details. Note that, to compile and run some of the test programs, the [pyplot_module.f90](https://github.com/jacobwilliams/pyplot-fortran) file (which is used to generate plots) must be copied into the `src/tests` directory.

# Compiling

A simple bash script ```build.sh``` is provided for building bspline-fortran with gfortran using [FoBiS](https://github.com/szaghi/FoBiS). It also builds the API documentation using [FORD](https://github.com/Fortran-FOSS-Programmers/ford). The library can also be compiled with the Intel Fortran Compiler (and presumably any other Fortran compiler that supports modern standards).

A basic CMake configuration file is also included. For example, to build a static library:

```bash
 mkdir build
 cd build
 cmake ..
 make
```

Or, to build a shared library:

```bash
 cmake -DBUILD_SHARED_LIBS=ON ..
```

For a debug build:
```bash
 cmake -DCMAKE_BUILD_TYPE=DEBUG ..
```

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`bspline-fortran.fobis`) is also provided that can also build the library and examples. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f bspline-fortran.fobis -mode tests-gnu`
  * To build all the examples using ifort: `FoBiS.py build -f bspline-fortran.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f bspline-fortran.fobis -mode static-gnu`
  * To build a static library using ifort: `FoBiS.py build -f bspline-fortran.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`

  To generate the documentation using [ford](https://github.com/cmacmackin/ford), run: ```FoBis.py rule --execute makedoc -f bspline-fortran.fobis```

# Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/bspline-fortran/). This was generated from the source code using [FORD](https://github.com/cmacmackin/ford) (note that the build script will also generate these files).

# License

The bspline-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).
