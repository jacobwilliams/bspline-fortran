!> author: Jacob Williams
!  license: BSD
! 
! Description
!---------
!  Multidimensional (2D-6D) B-Spline interpolation of data on a regular grid.
! 
! Notes
!---------
!  This module is based on the bspline and spline routines from [1].
!  The original Fortran 77 routines were converted to free-form source.
!  Some of them are relatively unchanged from the originals, but some have
!  been extensively refactored. In addition, new routines for 
!  4d, 5d, and 6d interpolation were also created (these are simply 
!  extensions of the same algorithm into higher dimensions).
! 
! References
!---------
!  1. DBSPLIN and DTENSBS from the
!     [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).
!     Original code is public domain.
!  2. Carl de Boor, "A Practical Guide to Splines",
!     Springer-Verlag, New York, 1978.
!  3. Carl de Boor, [Efficient Computer Manipulation of Tensor
!     Products](http://dl.acm.org/citation.cfm?id=355831), 
!     ACM Transactions on Mathematical Software,
!     Vol. 5 (1979), p. 173-182.
!  4. D.E. Amos, "Computation with Splines and B-Splines",
!     SAND78-1968, Sandia Laboratories, March, 1979.
!  5. Carl de Boor, 
!     [Package for calculating with B-splines](http://epubs.siam.org/doi/abs/10.1137/0714026),
!     SIAM Journal on Numerical Analysis 14, 3 (June 1977), p. 441-472.

    module bspline_module

    use,intrinsic :: iso_fortran_env, only: real64
    
    implicit none
    
    private

    integer,parameter :: wp = real64  !! Real precision

    !main routines:
    public :: db2ink, db2val
    public :: db3ink, db3val
    public :: db4ink, db4val
    public :: db5ink, db5val
    public :: db6ink, db6val
    
    contains
    
!*****************************************************************************************

!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the two-dimensional gridded data
!  $$ [x(i),y(j),\mathrm{fcn}(i,j)] ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
!  The interpolating function and its derivatives may
!  subsequently be evaluated by the function [[db2val]].
! 
!  The interpolating function is a piecewise polynomial function
!  represented as a tensor product of one-dimensional b-splines. the
!  form of this function is
! 
!  $$ s(x,y) = \sum_{i=1}^{n_x} \sum_{j=1}^{n_y} a_{ij} u_i(x) v_j(y) $$
! 
!  where the functions \(u_i\) and \(v_j\) are one-dimensional b-spline
!  basis functions. the coefficients \( a_{ij} \) are chosen so that
! 
!  $$ s(x(i),y(j)) = \mathrm{fcn}(i,j) ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
! 
!  Note that for each fixed value of y, \( s(x,y) \) is a piecewise
!  polynomial function of x alone, and for each fixed value of x, \( s(x,y) \)
!  is a piecewise polynomial function of y alone. in one dimension
!  a piecewise polynomial may be created by partitioning a given
!  interval into subintervals and defining a distinct polynomial piece
!  on each one. the points where adjacent subintervals meet are called
!  knots. each of the functions \(u_i\) and \(v_j\) above is a piecewise
!  polynomial.
! 
!  Users of db2ink choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the x and
!  y directions (kx and ky). users also may define their own knot
!  sequence in x and y separately (tx and ty). if iflag=0, however,
!  db2ink will choose sequences of knots that result in a piecewise
!  polynomial interpolant with kx-2 continuous partial derivatives in
!  x and ky-2 continuous partial derivatives in y. (kx knots are taken
!  near each endpoint in the x direction, not-a-knot end conditions
!  are used, and the remaining knots are placed at data points if kx
!  is even or at midpoints between data points if kx is odd. the y
!  direction is treated similarly.)
! 
!  After a call to db2ink, all information necessary to define the
!  interpolating function are contained in the parameters nx, ny, kx,
!  ky, tx, ty, and bcoef. These quantities should not be altered until
!  after the last call of the evaluation routine [[db2val]].
! 
! History
!---------
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    subroutine db2ink(x,nx,y,ny,fcn,kx,ky,tx,ty,bcoef,iflag)

    implicit none

    integer,intent(in)                      :: nx     !! Number of x abcissae
    integer,intent(in)                      :: ny     !! Number of y abcissae  
    integer,intent(in)                      :: kx     !! The order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                      :: ky     !! The order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    real(wp),dimension(nx),intent(in)       :: x      !! Array of x abcissae. Must be strictly increasing.
    real(wp),dimension(ny),intent(in)       :: y      !! Array of y abcissae. Must be strictly increasing.   
    real(wp),dimension(nx,ny),intent(in)    :: fcn    !! Array of function values to interpolate. fcn(i,j) should
                                                      !!    contain the function value at the point (x(i),y(j))
    real(wp),dimension(nx+kx),intent(inout) :: tx     !! The knots in the x direction for the spline interpolant.
                                                      !!    If iflag=0 these are chosen by db2ink.
                                                      !!    If iflag=1 these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(ny+ky),intent(inout) :: ty     !! The knots in the y direction for the spline interpolant.
                                                      !!    If iflag=0 these are chosen by db2ink.
                                                      !!    If iflag=1 these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(nx,ny),intent(out)   :: bcoef  !! Array of coefficients of the b-spline interpolant.
    integer,intent(inout)                   :: iflag  !! **on input:**  0 = knot sequence chosen by db2ink.
                                                      !!                1 = knot sequence chosen by user.
                                                      !! **on output:** 1 = successful execution.
                                                      !!                2 = iflag out of range.
                                                      !!                3 = nx out of range.
                                                      !!                4 = kx out of range.
                                                      !!                5 = x not strictly increasing.
                                                      !!                6 = tx not non-decreasing.
                                                      !!                7 = ny out of range.
                                                      !!                8 = ky out of range.
                                                      !!                9 = y not strictly increasing.
                                                      !!               10 = ty not non-decreasing.


    real(wp),dimension(nx*ny) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1))) :: work
    logical :: status_ok
  
    !check validity of inputs
    
    call check_inputs('db2ink',&
                        iflag,&
                        nx=nx,ny=ny,&
                        kx=kx,ky=ky,&
                        x=x,y=y,&
                        tx=tx,ty=ty,&
                        status_ok=status_ok)
                        
    if (status_ok) then

        !choose knots

        if (iflag == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
        end if

        !construct b-spline coefficients

        call dbtpcf(x,nx,fcn, nx,ny,tx,kx,temp, work)
        call dbtpcf(y,ny,temp,ny,nx,ty,ky,bcoef,work)

        iflag = 1
    
    end if

    end subroutine db2ink
!*****************************************************************************************
     
!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db2ink]] or one of its
!  derivatives at the point (xval,yval). 
!
!  To evaluate the interpolant
!  itself, set idx=idy=0, to evaluate the first partial with respect
!  to x, set idx=1,idy=0, and so on.
!
!  db2val returns 0.0 if (xval,yval) is out of range. that is, if
!```fortran
!   xval < tx(1) .or. xval > tx(nx+kx) .or. 
!   yval < ty(1) .or. yval > ty(ny+ky)
!```
!  if the knots tx and ty were chosen by [[db2ink]], then this is equivalent to:
!```fortran
!   xval < x(1) .or. xval > x(nx)+epsx .or. 
!   yval < y(1) .or. yval > y(ny)+epsy
!```
!  where
!```fortran 
!   epsx = 0.1*(x(nx)-x(nx-1))
!   epsy = 0.1*(y(ny)-y(ny-1)) 
!```
!
!  The input quantities tx, ty, nx, ny, kx, ky, and bcoef should be
!  unchanged since the last call of [[db2ink]].
!
! History
!---------
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag)

    implicit none

    integer,intent(in)                   :: idx   !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: idy   !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: nx    !! the number of interpolation points in x. (same as in last call to db2ink)
    integer,intent(in)                   :: ny    !! the number of interpolation points in y. (same as in last call to db2ink)
    integer,intent(in)                   :: kx    !! order of polynomial pieces in x. (same as in last call to db2ink)
    integer,intent(in)                   :: ky    !! order of polynomial pieces in y. (same as in last call to db2ink)
    real(wp),intent(in)                  :: xval  !! x coordinate of evaluation point.
    real(wp),intent(in)                  :: yval  !! y coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx    !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to db2ink)
    real(wp),dimension(ny+ky),intent(in) :: ty    !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to db2ink)
    real(wp),dimension(nx,ny),intent(in) :: bcoef !! the b-spline coefficients computed by db2ink.
    real(wp),intent(out)                 :: f     !! interpolated value
    integer,intent(out)                  :: iflag !! status flag: 0 : no errors, /=0 : error

    integer :: inbv, k, lefty, mflag, kcol
    real(wp),dimension(ky) :: temp
    real(wp),dimension(3*max(kx,ky)) :: work

    integer,save :: iloy = 1
    integer,save :: inbvx = 1
    
    f = 0.0_wp
    iflag = 1

    if (xval<tx(1) .or. xval>tx(nx+kx)) return
    if (yval<ty(1) .or. yval>ty(ny+ky)) return
    
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    
    iflag = 0
    inbv = 1

    kcol = lefty - ky
    do k=1,ky
        kcol = kcol + 1
        temp(k) = dbvalu(tx,bcoef(:,kcol),nx,kx,idx,xval,inbvx,work)
    end do
    
    kcol = lefty - ky + 1
    f = dbvalu(ty(kcol:),temp,ky,ky,idy,yval,inbv,work)
         
    end subroutine db2val
!*****************************************************************************************
      
!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the three-dimensional gridded data 
!  $$ [x(i),y(j),z(k),\mathrm{fcn}(i,j,k)] ~\mathrm{for}~ 
!     i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y, ~\mathrm{and}~ k=1,..,n_z $$
!  The interpolating function and
!  its derivatives may subsequently be evaluated by the function
!  [[db3val]].
! 
!  The interpolating function is a piecewise polynomial function
!  represented as a tensor product of one-dimensional b-splines. the
!  form of this function is
!  $$ s(x,y,z) = \sum_{i=1}^{n_x} \sum_{j=1}^{n_y} \sum_{k=1}^{n_z} 
!                a_{ijk} u_i(x) v_j(y) w_k(z) $$
! 
!  where the functions \(u_i\), \(v_j\), and \(w_k\) are one-dimensional b-
!  spline basis functions. the coefficients \(a_{ijk}\) are chosen so that:
! 
!  $$ s(x(i),y(j),z(k)) = \mathrm{fcn}(i,j,k) 
!     ~\mathrm{for}~ i=1,..,n_x , j=1,..,n_y , k=1,..,n_z $$
! 
!  Note that for fixed values of y and z s(x,y,z) is a piecewise
!  polynomial function of x alone, for fixed values of x and z s(x,y,z) 
!  is a piecewise polynomial function of y alone, and for fixed
!  values of x and y s(x,y,z) is a function of z alone. in one
!  dimension a piecewise polynomial may be created by partitioning a
!  given interval into subintervals and defining a distinct polynomial
!  piece on each one. the points where adjacent subintervals meet are
!  called knots. each of the functions \(u_i\), \(v_j\), and \(w_k\) above is a
!  piecewise polynomial.
! 
!  Users of db3ink choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the x, y,
!  and z directions (kx, ky, and kz). users also may define their own
!  knot sequence in x, y, and z separately (tx, ty, and tz). if iflag=
!  0, however, db3ink will choose sequences of knots that result in a
!  piecewise polynomial interpolant with kx-2 continuous partial
!  derivatives in x, ky-2 continuous partial derivatives in y, and kz-
!  2 continuous partial derivatives in z. (kx knots are taken near
!  each endpoint in x, not-a-knot end conditions are used, and the
!  remaining knots are placed at data points if kx is even or at
!  midpoints between data points if kx is odd. the y and z directions
!  are treated similarly.)
! 
!  After a call to db3ink, all information necessary to define the
!  interpolating function are contained in the parameters nx, ny, nz,
!  kx, ky, kz, tx, ty, tz, and bcoef. these quantities should not be
!  altered until after the last call of the evaluation routine [[db3val]].
! 
! History
!---------
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,tx,ty,tz,bcoef,iflag)

    implicit none

    integer,intent(in)                       :: nx    !! number of x abcissae (>= 3)
    integer,intent(in)                       :: ny    !! number of y abcissae (>= 3)
    integer,intent(in)                       :: nz    !! number of z abcissae (>= 3)
    integer,intent(in)                       :: kx    !! the order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                       :: ky    !! the order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    integer,intent(in)                       :: kz    !! the order of spline pieces in z (>= 2, < nz). (order = polynomial degree + 1)
    real(wp),dimension(nx),intent(in)        :: x     !! array of x abcissae. must be strictly increasing.
    real(wp),dimension(ny),intent(in)        :: y     !! array of y abcissae. must be strictly increasing.
    real(wp),dimension(nz),intent(in)        :: z     !! array of z abcissae. must be strictly increasing.
    real(wp),dimension(nx,ny,nz),intent(in)  :: fcn   !! array of function values to interpolate. fcn(i,j,k) should
                                                      !!   contain the function value at the point (x(i),y(j),z(k))
    real(wp),dimension(nx+kx),intent(inout)  :: tx    !! The knots in the x direction for the spline interpolant.
                                                      !!   If iflag=0 these are chosen by db3ink.
                                                      !!   If iflag=1 these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(ny+ky),intent(inout)  :: ty    !! The knots in the y direction for the spline interpolant.
                                                      !!    If iflag=0 these are chosen by db3ink.
                                                      !!    If iflag=1 these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(nz+kz),intent(inout)  :: tz    !! The knots in the z direction for the spline interpolant.
                                                      !!    If iflag=0 these are chosen by db3ink.
                                                      !!    If iflag=1 these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(nx,ny,nz),intent(out) :: bcoef !! array of coefficients of the b-spline interpolant.
    integer,intent(inout)                    :: iflag !! **on input**    0 = knot sequence chosen by db3ink.
                                                      !!                 1 = knot sequence chosen by user.
                                                      !! **on output**   1 = successful execution.
                                                      !!                 2 = iflag out of range.
                                                      !!                 3 = nx out of range.
                                                      !!                 4 = kx out of range.
                                                      !!                 5 = x not strictly increasing.
                                                      !!                 6 = tx not non-decreasing.
                                                      !!                 7 = ny out of range.
                                                      !!                 8 = ky out of range.
                                                      !!                 9 = y not strictly increasing.
                                                      !!                10 = ty not non-decreasing.
                                                      !!                11 = nz out of range.
                                                      !!                12 = kz out of range.
                                                      !!                13 = z not strictly increasing.
                                                      !!                14 = ty not non-decreasing.

    real(wp),dimension(nx*ny*nz) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))) :: work
    logical :: status_ok
      
    ! check validity of input

    call check_inputs('db3ink',&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,&
                        kx=kx,ky=ky,kz=kz,&
                        x=x,y=y,z=z,&
                        tx=tx,ty=ty,tz=tz,&
                        status_ok=status_ok)
                        
    if (status_ok) then

        ! choose knots

        if (iflag == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
        end if
    
        ! copy fcn to work in packed for dbtpcf
        temp(1:nx*ny*nz) = reshape( fcn, [nx*ny*nz] )

        ! construct b-spline coefficients
        
        call dbtpcf(x,nx,temp, nx,ny*nz,tx,kx,bcoef,work)
        call dbtpcf(y,ny,bcoef,ny,nx*nz,ty,ky,temp, work)
        call dbtpcf(z,nz,temp, nz,nx*ny,tz,kz,bcoef,work)
    
        iflag = 1
    
    end if

    end subroutine db3ink
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db3ink]] or one of its
!  derivatives at the point (xval,yval,zval). 
! 
!  To evaluate the
!  interpolant itself, set idx=idy=idz=0, to evaluate the first
!  partial with respect to x, set idx=1,idy=idz=0, and so on.
!
!  db3val returns 0.0 if (xval,yval,zval) is out of range. that is,
!```fortran
! xval<tx(1) .or. xval>tx(nx+kx) .or.
! yval<ty(1) .or. yval>ty(ny+ky) .or.
! zval<tz(1) .or. zval>tz(nz+kz)
!```
!  if the knots tx, ty, and tz were chosen by [[db3ink]], then this is
!  equivalent to
!```fortran
! xval<x(1) .or. xval>x(nx)+epsx .or.
! yval<y(1) .or. yval>y(ny)+epsy .or.
! zval<z(1) .or. zval>z(nz)+epsz
!```
!  where
!```fortran
! epsx = 0.1*(x(nx)-x(nx-1)) 
! epsy = 0.1*(y(ny)-y(ny-1))
! epsz = 0.1*(z(nz)-z(nz-1))
!```
!
!  The input quantities tx, ty, tz, nx, ny, nz, kx, ky, kz, and bcoef
!  should remain unchanged since the last call of [[db3ink]].
!
! History
!---------
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    subroutine db3val(xval,yval,zval,idx,idy,idz,&
                                     tx,ty,tz,&
                                     nx,ny,nz,kx,ky,kz,bcoef,f,iflag)

    implicit none

    integer,intent(in)                      :: idx    !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idy    !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idz    !! z derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: nx     !! the number of interpolation points in x. (same as in last call to db3ink)
    integer,intent(in)                      :: ny     !! the number of interpolation points in y. (same as in last call to db3ink) 
    integer,intent(in)                      :: nz     !! the number of interpolation points in z. (same as in last call to db3ink)
    integer,intent(in)                      :: kx     !! order of polynomial pieces in x. (same as in last call to db3ink)
    integer,intent(in)                      :: ky     !! order of polynomial pieces in y. (same as in last call to db3ink)
    integer,intent(in)                      :: kz     !! order of polynomial pieces in z. (same as in last call to db3ink)
    real(wp),intent(in)                     :: xval   !! x coordinate of evaluation point.
    real(wp),intent(in)                     :: yval   !! y coordinate of evaluation point.
    real(wp),intent(in)                     :: zval   !! z coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)    :: tx     !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to db3ink)
    real(wp),dimension(ny+ky),intent(in)    :: ty     !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to db3ink) 
    real(wp),dimension(nz+kz),intent(in)    :: tz     !! sequence of knots defining the piecewise polynomial in the z direction. (same as in last call to db3ink) 
    real(wp),dimension(nx,ny,nz),intent(in) :: bcoef  !! the b-spline coefficients computed by db3ink.
    real(wp),intent(out)                    :: f      !! interpolated value
    integer,intent(out)                     :: iflag  !! status flag: 0 : no errors, /=0 : error

    real(wp),dimension(ky,kz)              :: temp1
    real(wp),dimension(kz)                 :: temp2
    real(wp),dimension(3*max(kx,ky,kz))    :: work

    integer :: inbv1, inbv2, lefty, leftz, mflag,&
                kcoly, kcolz, izm1, j, k

    integer,save :: iloy  = 1
    integer,save :: iloz  = 1
    integer,save :: inbvx = 1
    
    f = 0.0_wp
    iflag = 1

    if (xval<tx(1) .or. xval>tx(nx+kx)) return
    if (yval<ty(1) .or. yval>ty(ny+ky)) return
    if (zval<tz(1) .or. zval>tz(nz+kz)) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag); if (mflag /= 0) return

    iflag = 0
    inbv1 = 1
    inbv2 = 1

    kcolz = leftz - kz
    do k=1,kz
        kcolz = kcolz + 1
        kcoly = lefty - ky
        do j=1,ky
            kcoly = kcoly + 1
            temp1(j,k) = dbvalu(tx,bcoef(:,kcoly,kcolz),nx,kx,idx,xval,inbvx,work)
        end do
    end do

    kcoly = lefty - ky + 1
    do k=1,kz
        temp2(k) = dbvalu(ty(kcoly:),temp1(:,k),ky,ky,idy,yval,inbv1,work)
    end do

    kcolz = leftz - kz + 1
    f = dbvalu(tz(kcolz:),temp2,kz,kz,idz,zval,inbv2,work)

    end subroutine db3val
!*****************************************************************************************

!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the four-dimensional gridded data
!  $$ [x(i),y(j),z(k),q(l),\mathrm{fcn}(i,j,k,l)] ~\mathrm{for}~ 
!     i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y, ~\mathrm{and}~ k=1,..,n_z, 
!     ~\mathrm{and}~ l=1,..,n_q $$
!  The interpolating function and its derivatives may 
!  subsequently be evaluated by the function [[db4val]].
!  
!  See [[db3ink]] header for more details.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db4ink(x,nx,y,ny,z,nz,q,nq,&
                        fcn,&
                        kx,ky,kz,kq,&
                        tx,ty,tz,tq,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                          :: nx    !! number of x abcissae (>= 3)
    integer,intent(in)                          :: ny    !! number of y abcissae (>= 3)
    integer,intent(in)                          :: nz    !! number of z abcissae (>= 3)
    integer,intent(in)                          :: nq    !! number of q abcissae (>= 3)
    integer,intent(in)                          :: kx    !! the order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                          :: ky    !! the order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    integer,intent(in)                          :: kz    !! the order of spline pieces in z (>= 2, < nz). (order = polynomial degree + 1)
    integer,intent(in)                          :: kq    !! the order of spline pieces in q (>= 2, < nq). (order = polynomial degree + 1)
    real(wp),dimension(nx),intent(in)           :: x     !! array of x abcissae. must be strictly increasing.
    real(wp),dimension(ny),intent(in)           :: y     !! array of y abcissae. must be strictly increasing.
    real(wp),dimension(nz),intent(in)           :: z     !! array of z abcissae. must be strictly increasing.
    real(wp),dimension(nq),intent(in)           :: q     !! array of q abcissae. must be strictly increasing.
    real(wp),dimension(nx,ny,nz,nq),intent(in)  :: fcn   !! array of function values to interpolate. fcn(i,j,k,q) should
                                                         !!   contain the function value at the point (x(i),y(j),z(k),q(l))
    real(wp),dimension(nx+kx),intent(inout)     :: tx    !! The knots in the x direction for the spline interpolant.
                                                         !!   If iflag=0 these are chosen by db4ink.
                                                         !!   If iflag=1 these are specified by the user.
                                                         !!    Must be non-decreasing.
    real(wp),dimension(ny+ky),intent(inout)     :: ty    !! The knots in the y direction for the spline interpolant.
                                                         !!    If iflag=0 these are chosen by db4ink.
                                                         !!    If iflag=1 these are specified by the user.
                                                         !!    Must be non-decreasing.
    real(wp),dimension(nz+kz),intent(inout)     :: tz    !! The knots in the z direction for the spline interpolant.
                                                         !!    If iflag=0 these are chosen by db4ink.
                                                         !!    If iflag=1 these are specified by the user.
                                                         !!    Must be non-decreasing.
    real(wp),dimension(nq+kq),intent(inout)     :: tq    !! The knots in the q direction for the spline interpolant.
                                                         !!    If iflag=0 these are chosen by db4ink.
                                                         !!    If iflag=1 these are specified by the user.
                                                         !!    Must be non-decreasing.
    real(wp),dimension(nx,ny,nz,nq),intent(out) :: bcoef !! array of coefficients of the b-spline interpolant.
    integer,intent(inout)                       :: iflag !! **on input**    0 = knot sequence chosen by db4ink.
                                                         !!                 1 = knot sequence chosen by user.
                                                         !! **on output**   1 = successful execution.
                                                         !!                 2 = iflag out of range.
                                                         !!                 3 = nx out of range.
                                                         !!                 4 = kx out of range.
                                                         !!                 5 = x not strictly increasing.
                                                         !!                 6 = tx not non-decreasing.
                                                         !!                 7 = ny out of range.
                                                         !!                 8 = ky out of range.
                                                         !!                 9 = y not strictly increasing.
                                                         !!                10 = ty not non-decreasing.
                                                         !!                11 = nz out of range.
                                                         !!                12 = kz out of range.
                                                         !!                13 = z not strictly increasing.
                                                         !!                14 = ty not non-decreasing.

           
    real(wp),dimension(nx*ny*nz*nq) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1),2*kq*(nq+1))) :: work
    logical :: status_ok
      
    ! check validity of input
    
    call check_inputs('db4ink',&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,&
                        kx=kx,ky=ky,kz=kz,kq=kq,&
                        x=x,y=y,z=z,q=q,&
                        tx=tx,ty=ty,tz=tz,tq=tq,&
                        status_ok=status_ok)
                        
    if (status_ok) then
    
        ! choose knots

        if (iflag == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
        end if

        ! construct b-spline coefficients

        call dbtpcf(x,nx,fcn,  nx,ny*nz*nq,tx,kx,temp, work)
        call dbtpcf(y,ny,temp, ny,nx*nz*nq,ty,ky,bcoef,work)
        call dbtpcf(z,nz,bcoef,nz,nx*ny*nq,tz,kz,temp, work)
        call dbtpcf(q,nq,temp, nq,nx*ny*nz,tq,kq,bcoef,work)
      
        iflag = 1
     
     end if

    end subroutine db4ink
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db4ink]] or one of its
!  derivatives at the point (xval,yval,zval,qval). 
! 
!  To evaluate the
!  interpolant itself, set idx=idy=idz=idq=0, to evaluate the first
!  partial with respect to x, set idx=1,idy=idz=idq=0, and so on.
! 
!  See [[db3val]] header for more information.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db4val(xval,yval,zval,qval,&
                                idx,idy,idz,idq,&
                                tx,ty,tz,tq,&
                                nx,ny,nz,nq,&
                                kx,ky,kz,kq,&
                                bcoef,f,iflag)

    implicit none

    integer,intent(in)                         :: idx    !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idy    !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idz    !! z derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idq    !! q derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: nx     !! the number of interpolation points in x. (same as in last call to db4ink)
    integer,intent(in)                         :: ny     !! the number of interpolation points in y. (same as in last call to db4ink) 
    integer,intent(in)                         :: nz     !! the number of interpolation points in z. (same as in last call to db4ink)
    integer,intent(in)                         :: nq     !! the number of interpolation points in q. (same as in last call to db4ink)
    integer,intent(in)                         :: kx     !! order of polynomial pieces in x. (same as in last call to db4ink)
    integer,intent(in)                         :: ky     !! order of polynomial pieces in y. (same as in last call to db4ink)
    integer,intent(in)                         :: kz     !! order of polynomial pieces in z. (same as in last call to db4ink)
    integer,intent(in)                         :: kq     !! order of polynomial pieces in q. (same as in last call to db4ink)
    real(wp),intent(in)                        :: xval   !! x coordinate of evaluation point.
    real(wp),intent(in)                        :: yval   !! y coordinate of evaluation point.
    real(wp),intent(in)                        :: zval   !! z coordinate of evaluation point.
    real(wp),intent(in)                        :: qval   !! q coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)       :: tx     !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to db4ink)
    real(wp),dimension(ny+ky),intent(in)       :: ty     !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to db4ink) 
    real(wp),dimension(nz+kz),intent(in)       :: tz     !! sequence of knots defining the piecewise polynomial in the z direction. (same as in last call to db4ink) 
    real(wp),dimension(nq+kq),intent(in)       :: tq     !! sequence of knots defining the piecewise polynomial in the q direction. (same as in last call to db4ink) 
    real(wp),dimension(nx,ny,nz,nq),intent(in) :: bcoef  !! the b-spline coefficients computed by db4ink.
    real(wp),intent(out)                       :: f      !! interpolated value
    integer,intent(out)                        :: iflag  !! status flag: 0 : no errors, /=0 : error
    
    real(wp),dimension(ky,kz,kq)             :: temp1
    real(wp),dimension(kz,kq)                :: temp2
    real(wp),dimension(kq)                   :: temp3
    real(wp),dimension(3*max(kx,ky,kz,kq))   :: work
    integer :: inbv1, inbv2, inbv3, lefty, leftz, leftq, mflag,&
                kcoly, kcolz, kcolq, i, j, k, q

    integer,save :: iloy  = 1
    integer,save :: iloz  = 1
    integer,save :: iloq  = 1
    integer,save :: inbvx = 1
    
    f = 0.0_wp
    iflag = 1

    if (xval<tx(1) .or. xval>tx(nx+kx) ) return
    if (yval<ty(1) .or. yval>ty(ny+ky) ) return
    if (zval<tz(1) .or. zval>tz(nz+kz) ) return
    if (qval<tq(1) .or. qval>tq(nq+kq) ) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag); if (mflag /= 0) return

    iflag = 0
    inbv1 = 1
    inbv2 = 1
    inbv3 = 1

    ! x -> y, z, q
    kcolq = leftq - kq
    do q=1,kq
        kcolq = kcolq + 1
        kcolz = leftz - kz
        do k=1,kz
            kcolz = kcolz + 1
            kcoly = lefty - ky
            do j=1,ky
                kcoly = kcoly + 1
                temp1(j,k,q) = dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq),&
                                        nx,kx,idx,xval,inbvx,work)
            end do
        end do
    end do

    ! y -> z, q
    kcoly = lefty - ky + 1
    do q=1,kq
        do k=1,kz
            temp2(k,q) = dbvalu(ty(kcoly:),temp1(:,k,q),ky,ky,idy,yval,inbv1,work)    
        end do
    end do

    ! z -> q
    kcolz = leftz - kz + 1
    do q=1,kq
        temp3(q) = dbvalu(tz(kcolz:),temp2(:,q),kz,kz,idz,zval,inbv2,work)    
    end do

    ! q
    kcolq = leftq - kq + 1
    f = dbvalu(tq(kcolq:),temp3,kq,kq,idq,qval,inbv3,work) 

    end subroutine db4val
!*****************************************************************************************

!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the five-dimensional gridded data (x(i),y(j),z(k),q(l),r(m),fcn(i,j,k,l,m)) for
!  i=1,..,nx, j=1,..,ny, k=1,..,nz, l=1,..,nq, and m=1,..,nr. 
!  The interpolating function and its derivatives may subsequently be evaluated 
!  by the function [[db5val]].
! 
!  See [[db3ink]] header for more details.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                        fcn,&
                        kx,ky,kz,kq,kr,&
                        tx,ty,tz,tq,tr,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                             :: nx    !! number of x abcissae (>= 3)
    integer,intent(in)                             :: ny    !! number of y abcissae (>= 3)
    integer,intent(in)                             :: nz    !! number of z abcissae (>= 3)
    integer,intent(in)                             :: nq    !! number of q abcissae (>= 3)
    integer,intent(in)                             :: nr    !! number of r abcissae (>= 3)
    integer,intent(in)                             :: kx    !! the order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                             :: ky    !! the order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    integer,intent(in)                             :: kz    !! the order of spline pieces in z (>= 2, < nz). (order = polynomial degree + 1)
    integer,intent(in)                             :: kq    !! the order of spline pieces in q (>= 2, < nq). (order = polynomial degree + 1)
    integer,intent(in)                             :: kr    !! the order of spline pieces in r (>= 2, < nr). (order = polynomial degree + 1)
    real(wp),dimension(nx),intent(in)              :: x     !! array of x abcissae. must be strictly increasing.
    real(wp),dimension(ny),intent(in)              :: y     !! array of y abcissae. must be strictly increasing.
    real(wp),dimension(nz),intent(in)              :: z     !! array of z abcissae. must be strictly increasing.
    real(wp),dimension(nq),intent(in)              :: q     !! array of q abcissae. must be strictly increasing.
    real(wp),dimension(nr),intent(in)              :: r     !! array of r abcissae. must be strictly increasing.
    real(wp),dimension(nx,ny,nz,nq,nr),intent(in)  :: fcn   !! array of function values to interpolate. fcn(i,j,k,q,r) should
                                                            !!   contain the function value at the point (x(i),y(j),z(k),q(l),r(m))
    real(wp),dimension(nx+kx),intent(inout)        :: tx    !! The knots in the x direction for the spline interpolant.
                                                            !!   If iflag=0 these are chosen by db5ink.
                                                            !!   If iflag=1 these are specified by the user.
                                                            !!    Must be non-decreasing.
    real(wp),dimension(ny+ky),intent(inout)        :: ty    !! The knots in the y direction for the spline interpolant.
                                                            !!    If iflag=0 these are chosen by db5ink.
                                                            !!    If iflag=1 these are specified by the user.
                                                            !!    Must be non-decreasing.
    real(wp),dimension(nz+kz),intent(inout)        :: tz    !! The knots in the z direction for the spline interpolant.
                                                            !!    If iflag=0 these are chosen by db5ink.
                                                            !!    If iflag=1 these are specified by the user.
                                                            !!    Must be non-decreasing.
    real(wp),dimension(nq+kq),intent(inout)        :: tq    !! The knots in the q direction for the spline interpolant.
                                                            !!    If iflag=0 these are chosen by db5ink.
                                                            !!    If iflag=1 these are specified by the user.
                                                            !!    Must be non-decreasing.
    real(wp),dimension(nr+kr),intent(inout)        :: tr    !! The knots in the r direction for the spline interpolant.
                                                            !!    If iflag=0 these are chosen by db5ink.
                                                            !!    If iflag=1 these are specified by the user.
                                                            !!    Must be non-decreasing.
    real(wp),dimension(nx,ny,nz,nq,nr),intent(out) :: bcoef !! array of coefficients of the b-spline interpolant.
    integer,intent(inout)                          :: iflag !! **on input**    0 = knot sequence chosen by db5ink.
                                                            !!                 1 = knot sequence chosen by user.
                                                            !! **on output**   1 = successful execution.
                                                            !!                 2 = iflag out of range.
                                                            !!                 3 = nx out of range.
                                                            !!                 4 = kx out of range.
                                                            !!                 5 = x not strictly increasing.
                                                            !!                 6 = tx not non-decreasing.
                                                            !!                 7 = ny out of range.
                                                            !!                 8 = ky out of range.
                                                            !!                 9 = y not strictly increasing.
                                                            !!                10 = ty not non-decreasing.
                                                            !!                11 = nz out of range.
                                                            !!                12 = kz out of range.
                                                            !!                13 = z not strictly increasing.
                                                            !!                14 = ty not non-decreasing.
           
    real(wp),dimension(nx*ny*nz*nq*nr) :: temp
    real(wp),dimension(max( 2*kx*(nx+1),&
                            2*ky*(ny+1),&
                            2*kz*(nz+1),&
                            2*kq*(nq+1),&
                            2*kr*(nr+1) )) :: work
    logical :: status_ok
      
    !  check validity of input
    
    call check_inputs('db5ink',&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,nr=nr,&
                        kx=kx,ky=ky,kz=kz,kq=kq,kr=kr,&
                        x=x,y=y,z=z,q=q,r=r,&
                        tx=tx,ty=ty,tz=tz,tq=tq,tr=tr,&
                        status_ok=status_ok)
                        
    if (status_ok) then
    
        !  choose knots

        if (iflag == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
        end if

        ! copy fcn to work in packed for dbtpcf
    
        temp(1:nx*ny*nz*nq*nr) = reshape( fcn, [nx*ny*nz*nq*nr] )

        !  construct b-spline coefficients
    
        call dbtpcf(x,nx,temp,  nx,ny*nz*nq*nr,tx,kx,bcoef, work)
        call dbtpcf(y,ny,bcoef, ny,nx*nz*nq*nr,ty,ky,temp,  work)
        call dbtpcf(z,nz,temp,  nz,nx*ny*nq*nr,tz,kz,bcoef, work)
        call dbtpcf(q,nq,bcoef, nq,nx*ny*nz*nr,tq,kq,temp,  work)
        call dbtpcf(r,nr,temp,  nr,nx*ny*nz*nq,tr,kr,bcoef, work)
      
        iflag = 1
     
     end if

    end subroutine db5ink
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db5ink]] or one of its
!  derivatives at the point (xval,yval,zval,qval,rval).
! 
!  To evaluate the
!  interpolant itself, set idx=idy=idz=idq=idr=0, to evaluate the first
!  partial with respect to x, set idx=1,idy=idz=idq=idr=0, and so on.
!  
!  See [[db3val]] header for more information.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db5val(xval,yval,zval,qval,rval,&
                                idx,idy,idz,idq,idr,&
                                tx,ty,tz,tq,tr,&
                                nx,ny,nz,nq,nr,&
                                kx,ky,kz,kq,kr,&
                                bcoef,f,iflag)

    implicit none

    integer,intent(in)                            :: idx    !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idy    !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idz    !! z derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idq    !! q derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idr    !! r derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: nx     !! the number of interpolation points in x. (same as in last call to db5ink)
    integer,intent(in)                            :: ny     !! the number of interpolation points in y. (same as in last call to db5ink) 
    integer,intent(in)                            :: nz     !! the number of interpolation points in z. (same as in last call to db5ink)
    integer,intent(in)                            :: nq     !! the number of interpolation points in q. (same as in last call to db5ink)
    integer,intent(in)                            :: nr     !! the number of interpolation points in r. (same as in last call to db5ink)
    integer,intent(in)                            :: kx     !! order of polynomial pieces in x. (same as in last call to db5ink)
    integer,intent(in)                            :: ky     !! order of polynomial pieces in y. (same as in last call to db5ink)
    integer,intent(in)                            :: kz     !! order of polynomial pieces in z. (same as in last call to db5ink)
    integer,intent(in)                            :: kq     !! order of polynomial pieces in q. (same as in last call to db5ink)
    integer,intent(in)                            :: kr     !! order of polynomial pieces in r. (same as in last call to db5ink)
    real(wp),intent(in)                           :: xval   !! x coordinate of evaluation point.
    real(wp),intent(in)                           :: yval   !! y coordinate of evaluation point.
    real(wp),intent(in)                           :: zval   !! z coordinate of evaluation point.
    real(wp),intent(in)                           :: qval   !! q coordinate of evaluation point.
    real(wp),intent(in)                           :: rval   !! r coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)          :: tx     !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to db5ink)
    real(wp),dimension(ny+ky),intent(in)          :: ty     !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to db5ink) 
    real(wp),dimension(nz+kz),intent(in)          :: tz     !! sequence of knots defining the piecewise polynomial in the z direction. (same as in last call to db5ink) 
    real(wp),dimension(nq+kq),intent(in)          :: tq     !! sequence of knots defining the piecewise polynomial in the q direction. (same as in last call to db5ink) 
    real(wp),dimension(nr+kr),intent(in)          :: tr     !! sequence of knots defining the piecewise polynomial in the r direction. (same as in last call to db5ink) 
    real(wp),dimension(nx,ny,nz,nq,nr),intent(in) :: bcoef  !! the b-spline coefficients computed by db5ink.
    real(wp),intent(out)                          :: f      !! interpolated value
    integer,intent(out)                           :: iflag  !! status flag: 0 : no errors, /=0 : error
    
    real(wp),dimension(ky,kz,kq,kr)           :: temp1
    real(wp),dimension(kz,kq,kr)              :: temp2
    real(wp),dimension(kq,kr)                 :: temp3
    real(wp),dimension(kr)                    :: temp4
    real(wp),dimension(3*max(kx,ky,kz,kq,kr)) :: work
    integer :: inbv1, inbv2, inbv3, inbv4,&
                lefty, leftz, leftq, leftr, mflag,&
                kcoly, kcolz, kcolq, kcolr, i, j, k, q, r

    integer,save :: iloy  = 1
    integer,save :: iloz  = 1
    integer,save :: iloq  = 1
    integer,save :: ilor  = 1
    integer,save :: inbvx = 1
    
    f = 0.0_wp
    iflag = 1
        
    if ( xval<tx(1) .or. xval>tx(nx+kx) ) return
    if ( yval<ty(1) .or. yval>ty(ny+ky) ) return
    if ( zval<tz(1) .or. zval>tz(nz+kz) ) return
    if ( qval<tq(1) .or. qval>tq(nq+kq) ) return
    if ( rval<tr(1) .or. rval>tr(nr+kr) ) return
    
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag); if (mflag /= 0) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,mflag); if (mflag /= 0) return
        
    iflag = 0
    inbv1 = 1
    inbv2 = 1
    inbv3 = 1
    inbv4 = 1

    ! x -> y, z, q, r
    kcolr = leftr - kr
    do r=1,kr
        kcolr = kcolr + 1
        kcolq = leftq - kq
        do q=1,kq
            kcolq = kcolq + 1
            kcolz = leftz - kz
            do k=1,kz
                kcolz = kcolz + 1
                kcoly = lefty - ky
                do j=1,ky
                    kcoly = kcoly + 1
                    temp1(j,k,q,r) = dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr),&
                                            nx,kx,idx,xval,inbvx,work)
                end do
            end do
        end do    
    end do
    
    ! y -> z, q, r
    kcoly = lefty - ky + 1
    do r=1,kr
        do q=1,kq
            do k=1,kz
                temp2(k,q,r) = dbvalu(ty(kcoly:),temp1(:,k,q,r),ky,ky,idy,yval,inbv1,work)    
            end do
        end do
    end do
    
    ! z -> q, r
    kcolz = leftz - kz + 1    
    do r=1,kr
        do q=1,kq
            temp3(q,r) = dbvalu(tz(kcolz:),temp2(:,q,r),kz,kz,idz,zval,inbv2,work)    
        end do    
    end do
    
    ! q -> r
    kcolq = leftq - kq + 1
    do r=1,kr
        temp4(r) = dbvalu(tq(kcolq:),temp3(:,r),kq,kq,idq,qval,inbv3,work)    
    end do
    
    ! r
    kcolr = leftr - kr + 1
    f = dbvalu(tr(kcolr:),temp4,kr,kr,idr,rval,inbv4,work) 
    
    end subroutine db5val
!*****************************************************************************************

!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the six-dimensional gridded data (x(i),y(j),z(k),q(l),r(m),s(n),fcn(i,j,k,l,m,n)) for
!  i=1,..,nx, j=1,..,ny, k=1,..,nz, l=1,..,nq, m=1,..,nr, n=1,..,ns. 
!  the interpolating function and its derivatives may subsequently be evaluated 
!  by the function [[db6val]].
!  
!  See [[db3ink]] header for more details.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                        fcn,&
                        kx,ky,kz,kq,kr,ks,&
                        tx,ty,tz,tq,tr,ts,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                                :: nx    !! number of x abcissae (>= 3)
    integer,intent(in)                                :: ny    !! number of y abcissae (>= 3)
    integer,intent(in)                                :: nz    !! number of z abcissae (>= 3)
    integer,intent(in)                                :: nq    !! number of q abcissae (>= 3)
    integer,intent(in)                                :: nr    !! number of r abcissae (>= 3)
    integer,intent(in)                                :: ns    !! number of s abcissae (>= 3)
    integer,intent(in)                                :: kx    !! the order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                                :: ky    !! the order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    integer,intent(in)                                :: kz    !! the order of spline pieces in z (>= 2, < nz). (order = polynomial degree + 1)
    integer,intent(in)                                :: kq    !! the order of spline pieces in q (>= 2, < nq). (order = polynomial degree + 1)
    integer,intent(in)                                :: kr    !! the order of spline pieces in r (>= 2, < nr). (order = polynomial degree + 1)
    integer,intent(in)                                :: ks    !! the order of spline pieces in s (>= 2, < ns). (order = polynomial degree + 1)
    real(wp),dimension(nx),intent(in)                 :: x     !! array of x abcissae. must be strictly increasing.
    real(wp),dimension(ny),intent(in)                 :: y     !! array of y abcissae. must be strictly increasing.
    real(wp),dimension(nz),intent(in)                 :: z     !! array of z abcissae. must be strictly increasing.
    real(wp),dimension(nq),intent(in)                 :: q     !! array of q abcissae. must be strictly increasing.
    real(wp),dimension(nr),intent(in)                 :: r     !! array of r abcissae. must be strictly increasing.
    real(wp),dimension(ns),intent(in)                 :: s     !! array of s abcissae. must be strictly increasing.
    real(wp),dimension(nx,ny,nz,nq,nr,ns),intent(in)  :: fcn   !! array of function values to interpolate. fcn(i,j,k,q,r,s) should
                                                               !!   contain the function value at the point (x(i),y(j),z(k),q(l),r(m),s(n))
    real(wp),dimension(nx+kx),intent(inout)           :: tx    !! The knots in the x direction for the spline interpolant.
                                                               !!   If iflag=0 these are chosen by db6ink.
                                                               !!   If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(ny+ky),intent(inout)           :: ty    !! The knots in the y direction for the spline interpolant.
                                                               !!    If iflag=0 these are chosen by db6ink.
                                                               !!    If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(nz+kz),intent(inout)           :: tz    !! The knots in the z direction for the spline interpolant.
                                                               !!    If iflag=0 these are chosen by db6ink.
                                                               !!    If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(nq+kq),intent(inout)           :: tq    !! The knots in the q direction for the spline interpolant.
                                                               !!    If iflag=0 these are chosen by db6ink.
                                                               !!    If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(nr+kr),intent(inout)           :: tr    !! The knots in the r direction for the spline interpolant.
                                                               !!    If iflag=0 these are chosen by db6ink.
                                                               !!    If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(ns+ks),intent(inout)           :: ts    !! The knots in the s direction for the spline interpolant.
                                                               !!    If iflag=0 these are chosen by db6ink.
                                                               !!    If iflag=1 these are specified by the user.
                                                               !!    Must be non-decreasing.
    real(wp),dimension(nx,ny,nz,nq,nr,ns),intent(out) :: bcoef !! array of coefficients of the b-spline interpolant.
    integer,intent(inout)                             :: iflag !! **on input**    0 = knot sequence chosen by db6ink.
                                                               !!                 1 = knot sequence chosen by user.
                                                               !! **on output**   1 = successful execution.
                                                               !!                 2 = iflag out of range.
                                                               !!                 3 = nx out of range.
                                                               !!                 4 = kx out of range.
                                                               !!                 5 = x not strictly increasing.
                                                               !!                 6 = tx not non-decreasing.
                                                               !!                 7 = ny out of range.
                                                               !!                 8 = ky out of range.
                                                               !!                 9 = y not strictly increasing.
                                                               !!                10 = ty not non-decreasing.
                                                               !!                11 = nz out of range.
                                                               !!                12 = kz out of range.
                                                               !!                13 = z not strictly increasing.
                                                               !!                14 = ty not non-decreasing.
           
    real(wp),dimension(nx*ny*nz*nq*nr*ns) :: temp
    real(wp),dimension(max( 2*kx*(nx+1),&
                            2*ky*(ny+1),&
                            2*kz*(nz+1),&
                            2*kq*(nq+1),&
                            2*kr*(nr+1),&
                            2*ks*(ns+1))) :: work
    logical :: status_ok
      
    ! check validity of input
    
    call check_inputs('db6ink',&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,nr=nr,ns=ns,&
                        kx=kx,ky=ky,kz=kz,kq=kq,kr=kr,ks=ks,&
                        x=x,y=y,z=z,q=q,r=r,s=s,&
                        tx=tx,ty=ty,tz=tz,tq=tq,tr=tr,ts=ts,&
                        status_ok=status_ok)
                        
    if (status_ok) then
    
        ! choose knots

        if (iflag == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
            call dbknot(s,ns,ks,ts)
        end if

        ! construct b-spline coefficients
    
        call dbtpcf(x,nx,fcn,  nx,ny*nz*nq*nr*ns,tx,kx,temp, work)
        call dbtpcf(y,ny,temp, ny,nx*nz*nq*nr*ns,ty,ky,bcoef,work)
        call dbtpcf(z,nz,bcoef,nz,nx*ny*nq*nr*ns,tz,kz,temp, work)
        call dbtpcf(q,nq,temp, nq,nx*ny*nz*nr*ns,tq,kq,bcoef,work)
        call dbtpcf(r,nr,bcoef,nr,nx*ny*nz*nq*ns,tr,kr,temp, work)
        call dbtpcf(s,ns,temp, ns,nx*ny*nz*nq*nr,ts,ks,bcoef,work)
      
        iflag = 1
     
     end if

    end subroutine db6ink
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db6ink]] or one of its
!  derivatives at the point (xval,yval,zval,qval,rval,sval).
! 
!  To evaluate the
!  interpolant itself, set idx=idy=idz=idq=idr=ids=0, to evaluate the first
!  partial with respect to x, set idx=1,idy=idz=idq=idr=ids=0, and so on.
!  
!  See [[db3val]] header for more information.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine db6val(xval,yval,zval,qval,rval,sval,&
                                idx,idy,idz,idq,idr,ids,&
                                tx,ty,tz,tq,tr,ts,&
                                nx,ny,nz,nq,nr,ns,&
                                kx,ky,kz,kq,kr,ks,&
                                bcoef,f,iflag)

    implicit none

    integer,intent(in)                               :: idx    !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idy    !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idz    !! z derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idq    !! q derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idr    !! r derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: ids    !! s derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: nx     !! the number of interpolation points in x. (same as in last call to db6ink)
    integer,intent(in)                               :: ny     !! the number of interpolation points in y. (same as in last call to db6ink) 
    integer,intent(in)                               :: nz     !! the number of interpolation points in z. (same as in last call to db6ink)
    integer,intent(in)                               :: nq     !! the number of interpolation points in q. (same as in last call to db6ink)
    integer,intent(in)                               :: nr     !! the number of interpolation points in r. (same as in last call to db6ink)
    integer,intent(in)                               :: ns     !! the number of interpolation points in s. (same as in last call to db6ink)
    integer,intent(in)                               :: kx     !! order of polynomial pieces in x. (same as in last call to db6ink)
    integer,intent(in)                               :: ky     !! order of polynomial pieces in y. (same as in last call to db6ink)
    integer,intent(in)                               :: kz     !! order of polynomial pieces in z. (same as in last call to db6ink)
    integer,intent(in)                               :: kq     !! order of polynomial pieces in q. (same as in last call to db6ink)
    integer,intent(in)                               :: kr     !! order of polynomial pieces in r. (same as in last call to db6ink)
    integer,intent(in)                               :: ks     !! order of polynomial pieces in s. (same as in last call to db6ink)
    real(wp),intent(in)                              :: xval   !! x coordinate of evaluation point.
    real(wp),intent(in)                              :: yval   !! y coordinate of evaluation point.
    real(wp),intent(in)                              :: zval   !! z coordinate of evaluation point.
    real(wp),intent(in)                              :: qval   !! q coordinate of evaluation point.
    real(wp),intent(in)                              :: rval   !! r coordinate of evaluation point.
    real(wp),intent(in)                              :: sval   !! s coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)             :: tx     !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to db6ink)
    real(wp),dimension(ny+ky),intent(in)             :: ty     !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to db6ink) 
    real(wp),dimension(nz+kz),intent(in)             :: tz     !! sequence of knots defining the piecewise polynomial in the z direction. (same as in last call to db6ink) 
    real(wp),dimension(nq+kq),intent(in)             :: tq     !! sequence of knots defining the piecewise polynomial in the q direction. (same as in last call to db6ink) 
    real(wp),dimension(nr+kr),intent(in)             :: tr     !! sequence of knots defining the piecewise polynomial in the r direction. (same as in last call to db6ink) 
    real(wp),dimension(ns+ks),intent(in)             :: ts     !! sequence of knots defining the piecewise polynomial in the s direction. (same as in last call to db6ink) 
    real(wp),dimension(nx,ny,nz,nq,nr,ns),intent(in) :: bcoef  !! the b-spline coefficients computed by db6ink.
    real(wp),intent(out)                             :: f      !! interpolated value
    integer,intent(out)                              :: iflag  !! status flag: 0 : no errors, /=0 : error
   
    real(wp),dimension(ky,kz,kq,kr,ks)            :: temp1
    real(wp),dimension(kz,kq,kr,ks)               :: temp2
    real(wp),dimension(kq,kr,ks)                  :: temp3
    real(wp),dimension(kr,ks)                     :: temp4
    real(wp),dimension(ks)                        :: temp5
    real(wp),dimension(3*max(kx,ky,kz,kq,kr,ks))  :: work
    
    integer :: inbv1,inbv2,inbv3,inbv4,inbv5,&
                lefty,leftz,leftq,leftr,lefts,&
                mflag,&
                kcoly,kcolz,kcolq,kcolr,kcols,&
                i,j,k,q,r,s

    integer,save :: iloy  = 1
    integer,save :: iloz  = 1
    integer,save :: iloq  = 1
    integer,save :: ilor  = 1
    integer,save :: ilos  = 1
    integer,save :: inbvx = 1
    
    f = 0.0_wp
    iflag = 1

    if (xval<tx(1) .or. xval>tx(nx+kx) ) return
    if (yval<ty(1) .or. yval>ty(ny+ky) ) return
    if (zval<tz(1) .or. zval>tz(nz+kz) ) return
    if (qval<tq(1) .or. qval>tq(nq+kq) ) return
    if (rval<tr(1) .or. rval>tr(nr+kr) ) return
    if (sval<ts(1) .or. sval>ts(ns+ks) ) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag); if (mflag /= 0) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,mflag); if (mflag /= 0) return
    call dintrv(ts,ns+ks,sval,ilos,lefts,mflag); if (mflag /= 0) return

    iflag = 0
    inbv1 = 1
    inbv2 = 1
    inbv3 = 1
    inbv4 = 1
    inbv5 = 1

    ! x -> y, z, q, r, s    
    kcols = lefts - ks
    do s=1,ks
        kcols = kcols + 1
        kcolr = leftr - kr
        do r=1,kr
            kcolr = kcolr + 1
            kcolq = leftq - kq
            do q=1,kq
                kcolq = kcolq + 1
                kcolz = leftz - kz
                do k=1,kz
                    kcolz = kcolz + 1
                    kcoly = lefty - ky
                    do j=1,ky
                        kcoly = kcoly + 1
                        temp1(j,k,q,r,s) = dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr,kcols),&
                                                    nx,kx,idx,xval,inbvx,work)
                    end do
                end do
            end do    
        end do
    end do
    
    ! y -> z, q, r, s
    kcoly = lefty - ky + 1
    do s=1,ks
        do r=1,kr
            do q=1,kq
                do k=1,kz
                    temp2(k,q,r,s) = dbvalu(ty(kcoly:),temp1(:,k,q,r,s),ky,ky,idy,yval,inbv1,work)    
                end do
            end do
        end do
    end do
    
    ! z -> q, r, s
    kcolz = leftz - kz + 1    
    do s=1,ks
        do r=1,kr
            do q=1,kq
                temp3(q,r,s) = dbvalu(tz(kcolz:),temp2(:,q,r,s),kz,kz,idz,zval,inbv2,work)    
            end do    
        end do
    end do
    
    ! q -> r, s
    kcolq = leftq - kq + 1
    do s=1,ks
        do r=1,kr
            temp4(r,s) = dbvalu(tq(kcolq:),temp3(:,r,s),kq,kq,idq,qval,inbv3,work)    
        end do
    end do
    
    ! r -> s
    kcolr = leftr - kr + 1
    do s=1,ks
        temp5(s) = dbvalu(tr(kcolr:),temp4(:,s),kr,kr,idr,rval,inbv4,work)    
    end do
        
    ! s
    kcols = lefts - ks + 1
    f = dbvalu(ts(kcols:),temp5,ks,ks,ids,sval,inbv5,work) 

    end subroutine db6val
!*****************************************************************************************

!*****************************************************************************************
!> Check the validity of the inputs to the "ink" routines.
!  Prints warning message if there is an error, 
!  and also sets iflag and status_ok.
! 
!  Supports up to 6D: x,y,z,q,r,s
! 
! Notes
!---------
!  The code is new, but the logic is based on the original
!  logic in the CMLIB routines db2ink and db3ink.
! 
! History
!---------
!  * Jacob Williams, 2/24/2015 : Created this routine.

    subroutine check_inputs(routine,&
                            iflag,&
                            nx,ny,nz,nq,nr,ns,&
                            kx,ky,kz,kq,kr,ks,&
                            x,y,z,q,r,s,&
                            tx,ty,tz,tq,tr,ts,&
                            status_ok)

    implicit none
    
    character(len=*),intent(in)                :: routine
    integer,intent(inout)                      :: iflag
    integer,intent(in),optional                :: nx,ny,nz,nq,nr,ns
    integer,intent(in),optional                :: kx,ky,kz,kq,kr,ks
    real(wp),dimension(:),intent(in),optional  :: x,y,z,q,r,s
    real(wp),dimension(:),intent(in),optional  :: tx,ty,tz,tq,tr,ts
    logical,intent(out)                        :: status_ok
    
    logical :: error
    
    status_ok = .false.
    
    if ((iflag < 0) .or. (iflag > 1)) then
    
          write(*,*) trim(routine)//' - iflag is out of range: ',iflag
          iflag = 2
          
    else
    
        call check('x',nx,kx,x,tx,[3,4,5,6],    error); if (error) return
        call check('y',ny,ky,y,ty,[7,8,9,10],   error); if (error) return
        call check('z',nz,kz,z,tz,[11,12,13,14],error); if (error) return
        call check('q',nq,kq,q,tq,[15,16,17,18],error); if (error) return
        call check('r',nr,kr,r,tr,[19,20,21,22],error); if (error) return
        call check('s',ns,ks,s,ts,[23,24,25,26],error); if (error) return

        status_ok = .true.
    
    end if
          
    contains
    
        subroutine check(s,n,k,x,t,ierrs,error)  !check t,x,n,k for validity
        
        implicit none
        
        character(len=1),intent(in),optional       :: s        !! coordinate string: 'x','y','z','q','r','s'
        integer,intent(in),optional                :: n        !! size of x
        integer,intent(in),optional                :: k        !! order
        real(wp),dimension(:),intent(in),optional  :: x        !! abcissae vector
        real(wp),dimension(:),intent(in),optional  :: t        !! knot vector size(n+k)
        integer,dimension(4),intent(in)            :: ierrs    !! int error codes for n,k,x,t checks
        logical,intent(out)                        :: error    !! true if there was an error
        
        if (present(n)) then
            call check_n('n'//s,n,ierrs(1),error); if (error) return
            if (present(k)) then
                call check_k('k'//s,k,n,ierrs(2),error); if (error) return
            end if
            if (present(x)) then
                call check_x(s,n,x,ierrs(3),error); if (error) return
            end if
            if (iflag /= 0) then
                if (present(k) .and. present(t)) then
                    call check_t('t'//s,n,k,t,ierrs(4),error); if (error) return
                end if
            end if
        end if
        
        end subroutine check
        
        subroutine check_n(s,n,ierr,error)
        
        implicit none
        
        character(len=*),intent(in) :: s
        integer,intent(in)          :: n
        integer,intent(in)          :: ierr
        logical,intent(out)         :: error
        
        if (n < 3) then
            write(*,*) trim(routine)//' - '//trim(s)//' is out of range: ',n
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if        
          
          end subroutine check_n

        subroutine check_k(s,k,n,ierr,error)
        
        implicit none
        
        character(len=*),intent(in) :: s
        integer,intent(in)          :: k
        integer,intent(in)          :: n
        integer,intent(in)          :: ierr
        logical,intent(out)         :: error
        
        if ((k < 2) .or. (k >= n)) then              
              write(*,*) trim(routine)//' - '//trim(s)//' is out of range: ',k
              iflag = ierr
              error = .true.
          else
              error = .false.
        end if
                  
          end subroutine check_k
          
          subroutine check_x(s,n,x,ierr,error)
          
        implicit none
        
        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        real(wp),dimension(:),intent(in)  :: x
        integer,intent(in)                :: ierr
        logical,intent(out)               :: error
        
        integer :: i
        
        error = .true.
        do i=2,n
            if (x(i) <= x(i-1)) then
                  iflag = ierr
                write(*,*) trim(routine)//' - '//trim(s)//&
                            ' array must be strictly increasing'
                return
            end if 
        end do          
          error = .false.
                  
          end subroutine check_x
          
          subroutine check_t(s,n,k,t,ierr,error)
          
        implicit none
        
        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        integer,intent(in)                :: k
        real(wp),dimension(:),intent(in)  :: t
        integer,intent(in)                :: ierr
        logical,intent(out)               :: error
        
        integer :: i
        
        error = .true.
        do i=2,n + k
            if (t(i) < t(i-1))  then
                  iflag = ierr
                write(*,*) trim(routine)//' - '//trim(s)//&
                            ' array must be non-decreasing'
                return
            end if 
        end do          
          error = .false.
                  
          end subroutine check_t
          
    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
!> dbknot chooses a knot sequence for interpolation of order k at the
!  data points x(i), i=1,..,n.  the n+k knots are placed in the array
!  t.  k knots are placed at each endpoint and not-a-knot end
!  conditions are used.  the remaining knots are placed at data points
!  if n is even and between data points if n is odd.  the rightmost
!  knot is shifted slightly to the right to insure proper interpolation
!  at x(n) (see page 350 of the reference).
!
! History
!---------
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    subroutine dbknot(x,n,k,t)

    implicit none
    
    integer,intent(in)                 :: n
    integer,intent(in)                 :: k
    real(wp),dimension(n),intent(in)   :: x
    real(wp),dimension(:),intent(out)  :: t

    integer  :: i, j, ipj, npj, ip1, jstrt
    real(wp) :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1_wp*( x(n)-x(n-1) )
    do j=1,k
        t(j)   = x(1)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2) == 1)  then

        !case of odd k --  knots between data points
        
        i = (k-1)/2 - k
        ip1 = i + 1
        jstrt = k + 1
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5_wp*( x(ipj) + x(ipj+1) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2) - k
        jstrt = k+1
        do j=jstrt,n
            ipj = i + j
            t(j) = x(ipj)
        end do

    end if

    end subroutine dbknot
!*****************************************************************************************

!*****************************************************************************************
!> dbtpcf computes b-spline interpolation coefficients for nf sets
!  of data stored in the columns of the array fcn. the b-spline
!  coefficients are stored in the rows of bcoef however.
!  each interpolation is based on the n abcissa stored in the
!  array x, and the n+k knots stored in the array t. the order
!  of each interpolation is k. the work array must be of length
!  at least 2*k*(n+1).
!
! History
!---------
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work)

    integer  :: n, nf, ldf, k
    real(wp) :: x(n), fcn(ldf,nf), t(*), bcoef(nf,n), work(*)

    integer :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0)  then

        ! partition work array
        m1 = k - 1
        m2 = m1 + k
        iq = 1 + n
        iw = iq + m2*n+1

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw))
        do i=1,n
            bcoef(1,i) = work(i)
        end do

        !  all remaining data sets by back-substitution

        if (nf == 1)  return
        do j=2,nf
            do i=1,n
                work(i) = fcn(i,j)
            end do
            call dbnslv(work(iq),m2,n,m1,m1,work)
            do i=1,n
                bcoef(j,i) = work(i)
            end do
        end do
    
    end if

    end subroutine dbtpcf
!*****************************************************************************************

!***************************************************************************************** 
!> dbintk produces the b-spline coefficients, bcoef, of the
!  b-spline of order k with knots t(i), i=1,...,n+k, which
!  takes on the value y(i) at x(i), i=1,...,n.  the spline or
!  any of its derivatives can be evaluated by calls to [[dbvalu]].
!  
!  the i-th equation of the linear system a*bcoef = b for the
!  coefficients of the interpolant enforces interpolation at
!  x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
!  a band matrix with 2k-1 bands if a is invertible.  the matrix
!  a is generated row by row and stored, diagonal by diagonal,
!  in the rows of q, with the main diagonal going into row k.
!  the banded system is then solved by a call to dbnfac (which
!  constructs the triangular factorization for a and stores it
!  again in q), followed by a call to dbnslv (which then
!  obtains the solution bcoef by substitution).  dbnfac does no
!  pivoting, since the total positivity of the matrix a makes
!  this unnecessary.  the linear system to be solved is
!  (theoretically) invertible if and only if
!          t(i) < x(i) < t(i+k),        for all i.
!  equality is permitted on the left for i=1 and on the right
!  for i=n when k knots are used at x(1) or x(n).  otherwise,
!  violation of this condition is certain to lead to an error.
!  
! Error conditions
!---------
!  * improper input is a fatal error
!  * singular system of equations is a fatal error
! 
! History
!---------
!  * splint written by carl de boor [5]
!  * dbintk author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations. (jec)
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    subroutine dbintk(x,y,t,n,k,bcoef,q,work)

    implicit none

    integer,intent(in)                :: n      !!number of data points, n >= k
    real(wp),dimension(n),intent(in)  :: x      !!vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(wp),dimension(n),intent(in)  :: y      !!corresponding vector of length n containing data
                                                !! point ordinates.
    real(wp),dimension(*),intent(in)  :: t      !!knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer,intent(in)                :: k      !!order of the spline, k >= 1
    real(wp),dimension(n),intent(out) :: bcoef  !!a vector of length n containing the b-spline coefficients
    real(wp),dimension(*),intent(out) :: q      !!a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real(wp),dimension(*),intent(out) :: work   !!work vector of length 2*k

    integer :: iflag, iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(wp) :: xi

    if (k<1) then
        call xerror( ' dbintk,  k does not satisfy k>=1', 35, 2, 1)
        return
    end if

    if (n<k) then
        call xerror( ' dbintk,  n does not satisfy n>=k', 35, 2, 1)
        return
    end if

    jj = n - 1
    if (jj/=0) then
        do i=1,jj
            if (x(i)>=x(i+1)) then
                call xerror( ' dbintk,  x(i) does not satisfy x(i)<x(i+1) for some i', 57, 2, 1)
                return
            end if
        end do
    end if

    np1 = n + 1
    km1 = k - 1
    kpkm2 = 2*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1,lenq
        q(i) = 0.0_wp
    end do

    ! loop over i to construct the n interpolation equations
    do i=1,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! *** find  left  in the closed interval (i,i+k-1) such that
        !         t(left) <= x(i) < t(left+1)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            call xerror( ' dbintk,  some abscissa was not in the support of the'//&
                         ' corresponding basis function and the system is singular.',109,2,1)
            return
        end if
        do
            if (xi<t(left+1)) go to 30
            left = left + 1
            if (left>=ilp1mx) exit
        end do
        left = left - 1
        if (xi>t(left+1)) then
            call xerror( ' dbintk,  some abscissa was not in the support of the'//&
                         ' corresponding basis function and the system is singular.',109,2,1)
            return
        end if

        ! *** the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
 30     call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork)

        ! we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
        ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
        ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
        ! as a two-dim. array , with  2*k-1  rows (see comments in
        ! dbnfac). in the present program, we treat  q  as an equivalent
        ! one-dimensional array (because of fortran restrictions on
        ! dimension statements) . we therefore want  bcoef(j) to go into
        ! entry
        !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
        !            = i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
        ! of q.
        jj = i - left + 1 + (left-k)*(k+km1)
        do j=1,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! *** solve  a*bcoef = y  by backsubstitution
        do i=1,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
    else  !failure
        call xerror( ' dbintk,  the system of solver detects a singular system'//&
                     ' although the theoretical conditions for a solution were satisfied.',123,8,1)
    end if

    end subroutine dbintk
!*****************************************************************************************

!*****************************************************************************************
!> Returns in w the LU-factorization (without pivoting) of the banded
!  matrix a of order nrow with (nbandl + 1 + nbandu) bands or diagonals 
!  in the work array w .
!
!  gauss elimination without pivoting is used. the routine is
!  intended for use with matrices a which do not require row inter-
!  changes during factorization, especially for the totally
!  positive matrices which occur in spline calculations.
!  the routine should not be used for an arbitrary banded matrix.
!
! Work array
!---------
!
! **Input**
!
!        w array of size (nroww,nrow) contains the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!
! **Output**
!
!  * if  iflag = 1, then
!        w contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call dbnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!  * if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
! 
! History
!---------
!  * banfac written by carl de boor [5]
!  * dbnfac from CMLIB [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

      subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)
      
      integer,intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
      integer,intent(in) :: nrow    !! matrix order
      integer,intent(in) :: nbandl  !! number of bands of a below the main diagonal
      integer,intent(in) :: nbandu  !! number of bands of a above the main diagonal
      integer,intent(out) :: iflag  !! indicating success(=1) or failure (=2)
      real(wp),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

      integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
      real(wp) :: factor, pivot

      iflag = 1
      middle = nbandu + 1
!                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1

      !if (nrowm1) 120, 110, 10    !JW removed obsolescent arithmetic IF statement
      if (nrowm1 < 0) then;      goto 120
      elseif (nrowm1 == 0) then; goto 110
      else;                      goto 10
      end if

   10 if (nbandl>0) go to 30
!                a is upper triangular. check that diagonal is nonzero .
      do i=1,nrowm1
        if (w(middle,i)==0.0_wp) go to 120
      end do
      go to 110
   30 if (nbandu>0) go to 60
!              a is lower triangular. check that diagonal is nonzero and
!                 divide each column by its diagonal .
      do i=1,nrowm1
        pivot = w(middle,i)
        if (pivot==0.0_wp) go to 120
        jmax = min(nbandl,nrow-i)
        do j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
        end do
      end do
      return
!
!        a  is not just a triangular matrix. construct lu factorization
   60 do i=1,nrowm1
!                                  w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0_wp) go to 120
!                 jmax  is the number of (nonzero) entries in column  i
!                     below the diagonal .
        jmax = min(nbandl,nrow-i)
!              divide each entry in column  i  below diagonal by pivot .
        do j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
        end do
!                 kmax  is the number of (nonzero) entries in row  i  to
!                     the right of the diagonal .
        kmax = min(nbandu,nrow-i)
!                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
!                  (below row  i ) .
        do k=1,kmax
          ipk = i + k
          midmk = middle - k
          factor = w(midmk,ipk)
          do j=1,jmax
            w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
          end do
        end do
      end do
!                                       check the last diagonal entry .
  110 if (w(middle,nrow)/=0.0_wp) return
  120 iflag = 2

      end subroutine dbnfac
      
!*****************************************************************************************
!> Companion routine to [[dbnfac]]. it returns the solution x of the
!  linear system a*x = b in place of b, given the lu-factorization
!  for a in the work array w from dbnfac.
! 
!  (with \( a = l*u \), as stored in w), the unit lower triangular system
!  \( l(u*x) = b \) is solved for \( y = u*x \), and y stored in b. then the
!  upper triangular system \(u*x = y \) is solved for x. the calculations 
!  are so arranged that the innermost loops stay within columns.
! 
! History
!---------
!  * banslv written by carl de boor [5]
!  * dbnslv from SLATEC library [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)
    
    integer,intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    real(wp),dimension(nroww,nrow),intent(in) :: w    !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    real(wp),dimension(nrow),intent(inout) :: b  !! **in**: right side of the system to be solved
                                                 !! **out**: the solution x, of order nrow

    integer :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1
    if (nrow/=1) then

        nrowm1 = nrow - 1
        if (nbandl/=0) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l) 
            !                       from right side (below i-th row).
            do i=1,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0) then
            ! a is lower triangular.
            do i=1,nrow
                b(i) = b(i)/w(1,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1)
            do j=1,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1
            if (i<=1) exit
        end do

    end if

    b(1) = b(1)/w(middle,1)

    end subroutine dbnslv
!*****************************************************************************************

!*****************************************************************************************
!> Calculates the value of all (possibly) nonzero basis
!  functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
!  <= x <= t(n+1) and j=iwork is set inside the routine on
!  the first call when index=1.  ileft is such that t(ileft) <=
!  x < t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
!  produces the proper ileft.  dbspvn calculates using the basic
!  algorithm needed in dbspvd.  if only basis functions are
!  desired, setting jhigh=k and index=1 can be faster than
!  calling dbspvd, but extra coding is required for derivatives
!  (index=2) and dbspvd is set up for this purpose.
!  
!  left limiting values are set up as described in dbspvd.
! 
! Error Conditions
!---------
!  * improper input is a fatal error.
! 
! History
!---------
!  * bsplvn written by carl de boor [5]
!  * dbspvn author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork)

    implicit none

    real(wp),intent(in)  :: t(*)     !!knot vector of length n+k, where
                                     !! n = number of b-spline basis functions
                                     !! n = sum of knot multiplicities-k
                                     !!    dimension t(ileft+jhigh)
    integer,intent(in)   :: jhigh    !!order of b-spline, 1 <= jhigh <= k
    integer,intent(in)   :: k        !!highest possible order
    integer,intent(in)   :: index    !!index = 1 gives basis functions of order jhigh
                                     !!      = 2 denotes previous entry with work, iwork
                                     !!          values saved for subsequent calls to
                                     !!          dbspvn.
    real(wp),intent(in)  :: x        !!argument of basis functions, t(k) <= x <= t(n+1)
    integer,intent(in)   :: ileft    !!largest integer such that t(ileft) <= x < t(ileft+1)
    real(wp),intent(out) :: vnikx(k) !!vector of length k for spline values.
    real(wp),intent(out) :: work(*)  !!a work vector of length 2*k
    integer,intent(out)  :: iwork    !!a work parameter.  both work and iwork contain
                                     !! information necessary to continue for index = 2.
                                     !! when index = 1 exclusively, these are scratch
                                     !! variables and can be used for other purposes.

    integer :: imjp1, ipj, jp1, jp1ml, l
    real(wp) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i), 
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1) then
        call xerror( ' dbspvn,  k does not satisfy k>=1', 35, 2, 1)
        return
    end if
    if (jhigh>k .or. jhigh<1) then
        call xerror( ' dbspvn,  jhigh does not satisfy 1<=jhigh<=k', 48, 2, 1)
        return
    end if
    if (index<1 .or. index>2) then
        call xerror( ' dbspvn,  index is not 1 or 2',29,2,1)
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1)) then
        call xerror( ' dbspvn,  x does not satisfy t(ileft)<=x<=t(ileft+1)', 56, 2, 1)
        return
    end if

    if (index==1) then
        iwork = 1
        vnikx(1) = 1.0_wp
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0_wp
        jp1 = iwork + 1
        do l=1,iwork
            jp1ml = jp1 - l
            vm = vnikx(l)/(work(l)+work(k+jp1ml))
            vnikx(l) = vm*work(l) + vmprev
            vmprev = vm*work(k+jp1ml)
        end do
        vnikx(jp1) = vmprev
        iwork = jp1
        if (iwork>=jhigh) exit
    end do

    end subroutine dbspvn
!*****************************************************************************************
      
!*****************************************************************************************
!> Evaluates the b-representation (t,a,n,k) of a b-spline
!  at x for the function value on ideriv=0 or any of its
!  derivatives on ideriv=1,2,...,k-1.  right limiting values
!  (right derivatives) are returned except at the right end
!  point x=t(n+1) where left limiting values are computed.  the
!  spline is defined on t(k) <= x <= t(n+1).  dbvalu returns
!  a fatal error message when x is outside of this interval.
!  
!  to compute left derivatives or left limiting values at a
!  knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
! 
! Error Conditions
!---------
!  * improper input is a fatal error.
! 
! History
!---------
!  * bvalue written by carl de boor [5]
!  * dbvalu author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    real(wp) function dbvalu(t,a,n,k,ideriv,x,inbv,work)  !!value of the ideriv-th derivative at x

    implicit none

    integer,intent(in)               :: n       !!number of b-spline coefficients.
                                                !!  (sum of knot multiplicities-k)
    real(wp),dimension(:),intent(in) :: t       !!knot vector of length n+k
    real(wp),dimension(n),intent(in) :: a       !!b-spline coefficient vector of length n 
    integer,intent(in)               :: k       !!order of the b-spline, k >= 1
    integer,intent(in)               :: ideriv  !!order of the derivative, 0 <= ideriv <= k-1.
                                                !!  ideriv = 0 returns the b-spline value
    real(wp),intent(in)              :: x       !!argument, t(k) <= x <= t(n+1)
    integer,intent(inout)            :: inbv    !!an initialization parameter which must be set
                                                !!  to 1 the first time dbvalu is called.
                                                !!  inbv contains information for efficient process-
                                                !!  ing after the initial call and inbv must not
                                                !!  be changed by the user.  distinct splines require
                                                !!  distinct inbv parameters.
    real(wp),dimension(:)            :: work    !!work vector of length 3*k

    integer :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(wp) :: fkmj

    dbvalu = 0.0_wp

    if (k<1) then
        call xerror( ' dbvalu,  k does not satisfy k>=1' )
        return
    end if

    if (n<k) then
        call xerror( ' dbvalu,  n does not satisfy n>=k' )
        return
    end if

    if (ideriv<0 .or. ideriv>=k) then
        call xerror( ' dbvalu,  ideriv does not satisfy 0<=ideriv<k' )
        return
    end if

    kmider = k - ideriv

    ! *** find *i* in (k,n) such that t(i) <= x < t(i+1)
    !     (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1
    call dintrv(t, n+1, x, inbv, i, mflag)
    if (x<t(k)) then
        call xerror( ' dbvalu,  x is not greater than or equal to t(k)' )
        return
    end if

    if (mflag/=0) then

        if (x>t(i)) then
            call xerror( ' dbvalu,  x is not less than or equal to t(n+1)' )
            return
        end if

        do
            if (i==k) then
                call xerror( ' dbvalu,  a left limiting value cannot be obtained at t(k)' )
                return
            end if
            i = i - 1
            if (x/=t(i)) exit
        end do

    end if

    ! *** difference the coefficients *ideriv* times
    !     work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0) then
        do j=1,ideriv
            kmj = k - j
            fkmj = real(kmj,wp)
            do jj=1,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! *** compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    !     given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1
        kpk = k + k
        j1 = k + 1
        j2 = kpk + 1
        do j=1,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1
            j2 = j2 + 1
        end do
        iderp1 = ideriv + 1
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1,kmj
                work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    dbvalu = work(1)

    end function dbvalu
!*****************************************************************************************
      
!*****************************************************************************************
!> Computes the largest integer ileft in 1 <= ileft <= lxt 
!  such that xt(ileft) <= x where xt(*) is a subdivision of
!  the x interval.  
!  precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   mflag=0
!         if xt(lxt) <= x           then ileft=lxt, mflag=1
!```
!
!  that is, when multiplicities are present in the break point
!  to the left of x, the largest index is taken for ileft.
! 
! History
!---------
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.

    subroutine dintrv(xt,lxt,x,ilo,ileft,mflag)
    
    implicit none
    
    integer,intent(in)                 :: lxt    !! length of the xt vector
    real(wp),dimension(lxt),intent(in) :: xt     !! a knot or break point vector of length lxt
    real(wp),intent(in)                :: x      !! argument
    integer,intent(inout)              :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array xt is
                                                 !! processed by dintrv. ilo contains information for 
                                                 !! efficient processing after the initial call and ilo
                                                 !! must not be changed by the user.  distinct splines 
                                                 !! require distinct ilo parameters.
    integer,intent(out)                :: ileft  !! largest integer satisfying xt(ileft) <= x
    integer,intent(out)                :: mflag  !! signals when x lies out of bounds
    
    integer :: ihi, istep, middle

      ihi = ilo + 1
      if (ihi<lxt) go to 10
      if (x>=xt(lxt)) go to 110
      if (lxt<=1) go to 90
      ilo = lxt - 1
      ihi = lxt

   10 if (x>=xt(ihi)) go to 40
      if (x>=xt(ilo)) go to 100

! *** now x < xt(ihi) . find lower bound
      istep = 1
   20 ihi = ilo
      ilo = ihi - istep
      if (ilo<=1) go to 30
      if (x>=xt(ilo)) go to 70
      istep = istep*2
      go to 20

   30 ilo = 1
      if (x<xt(1)) go to 90
      go to 70

! *** now x >= xt(ilo) . find upper bound
   40 istep = 1
   50 ilo = ihi
      ihi = ilo + istep
      if (ihi>=lxt) go to 60
      if (x<xt(ihi)) go to 70
      istep = istep*2
      go to 50

   60 if (x>=xt(lxt)) go to 110
      ihi = lxt

! *** now xt(ilo) <= x < xt(ihi) . narrow the interval
   70 middle = (ilo+ihi)/2
      if (middle==ilo) go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x<xt(middle)) go to 80
      ilo = middle
      go to 70

   80 ihi = middle
      go to 70

! *** set output and return
   90 mflag = -1
      ileft = 1
      return

  100 mflag = 0
      ileft = ilo
      return

  110 mflag = 1
      ileft = lxt

    end subroutine dintrv      
!*****************************************************************************************
    
!*****************************************************************************************
!> Just a replacement for the CMLIB XERROR routine.

    subroutine xerror(messg,nmessg,nerr,level)

    use iso_fortran_env, only: error_unit

    implicit none

    character(len=*),intent(in) :: messg
    integer,intent(in),optional :: nmessg,nerr,level
    
    write(error_unit,'(A)') trim(messg)

    end subroutine xerror
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_module
!*****************************************************************************************