bspline-fortran
============

Multidimensional B-Spline Interpolation of Data on a Regular Grid

Brief description
---------------

The library provides subroutines for 2D-6D interpolation using b-splines. The code is written in modern Fortran.

The core routines are:

```Fortran

!f(x,y)
subroutine db2ink(x,nx,y,ny,fcn,kx,ky,tx,ty,bcoef,iflag)
function db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef)

!f(x,y,z)
subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,tx,ty,tz,bcoef,iflag)
function db3val(xval,yval,zval,idx,idy,idz,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef)

!f(x,y,z,q)
subroutine db4ink(x,nx,y,ny,z,nz,q,nq,fcn,kx,ky,kz,kq,tx,ty,tz,tq,bcoef,iflag)
function db4val(xval,yval,zval,qval,idx,idy,idz,idq,tx,ty,tz,tq,nx,ny,nz,nq,kx,ky,kz,kq,bcoef)

!f(x,y,z,q,r)
subroutine db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,bcoef,iflag)
function db5val(xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,tx,ty,tz,tq,tr,nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,bcoef)

!f(x,y,z,q,r,s)
subroutine db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,bcoef,iflag)
function db6val(xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,tx,ty,tz,tq,tr,ts,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,bcoef)
```
The ```ink``` routines compute the interpolant coefficients, and the ```val``` routines evalute the interpolant at the specified value of each coordinate.  Eventually, object-oriented wrappers will be created for these routines for ease of use.  The 2D and 3D routines are extensively refactored versions of the original routines from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).  The others are new, and are simply extensions of the same algorithm into higher dimensions.

License
--------

The bspline-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).

See also
---------------

* This code includes the public domain DBSPLIN and DTENSBS code from the [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm) (CMLIB)
