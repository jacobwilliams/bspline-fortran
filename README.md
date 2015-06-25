bspline-fortran [![GitHub release](https://img.shields.io/github/release/jacobwilliams/bspline-fortran.svg?style=plastic)](https://github.com/jacobwilliams/bspline-fortran/releases/latest)
============

Multidimensional B-Spline Interpolation of Data on a Regular Grid

# Status

[![Build Status](https://img.shields.io/travis/jacobwilliams/bspline-fortran/master.svg?style=plastic)](https://travis-ci.org/jacobwilliams/bspline-fortran)

# Brief description

The library provides subroutines for 2D-6D interpolation using b-splines. The code is written in modern Fortran (i.e., Fortran 2003+).

The core routines are:

```Fortran

!f(x,y)
subroutine db2ink(x,nx,y,ny,fcn,kx,ky,tx,ty,bcoef,iflag)
subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag)

!f(x,y,z)
subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,tx,ty,tz,bcoef,iflag)
subroutine db3val(xval,yval,zval,idx,idy,idz,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef,f,iflag)

!f(x,y,z,q)
subroutine db4ink(x,nx,y,ny,z,nz,q,nq,fcn,kx,ky,kz,kq,tx,ty,tz,tq,bcoef,iflag)
subroutine db4val(xval,yval,zval,qval,idx,idy,idz,idq,tx,ty,tz,tq,nx,ny,nz,nq,kx,ky,kz,kq,bcoef,f,iflag)

!f(x,y,z,q,r)
subroutine db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,bcoef,iflag)
subroutine db5val(xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,tx,ty,tz,tq,tr,nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,bcoef,f,iflag)

!f(x,y,z,q,r,s)
subroutine db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,bcoef,iflag)
subroutine db6val(xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,tx,ty,tz,tq,tr,ts,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,bcoef,f,iflag)
```

The ```ink``` routines compute the interpolant coefficients, and the ```val``` routines evalute the interpolant at the specified value of each coordinate.  Eventually, object-oriented wrappers will be created for these routines for ease of use.  The 2D and 3D routines are extensively refactored versions of the original routines from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).  The others are new, and are simply extensions of the same algorithm into higher dimensions.

In addition to the main subroutines, an object-oriented wrapper is also provided. For example:

```Fortran
type(bspline_3d) :: s3
call s3%initialize(x,y,z,fcn,kx,ky,kz,iflag)
call s3%evaluate(xval,yval,zval,idx,idy,idz,f,iflag)
call s3%destroy()
```
See the examples for more details.

# Compiling

A simple bash script ```build.sh``` is provided for building bspline-fortran with gfortran using [FoBiS](https://github.com/szaghi/FoBiS). It also builds the API documentation using [FORD](https://github.com/cmacmackin/ford).

# Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/bspline-fortran/). This was generated using [FORD](https://github.com/cmacmackin/ford).

# License

The bspline-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).
