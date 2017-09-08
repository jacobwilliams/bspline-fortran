!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!
!  Multidimensional (1D-6D) B-spline interpolation of data on a regular grid.
!  Basic pure subroutine interface.
!
!### Notes
!
!  This module is based on the B-spline and spline routines from [1].
!  The original Fortran 77 routines were converted to free-form source.
!  Some of them are relatively unchanged from the originals, but some have
!  been extensively refactored. In addition, new routines for
!  1d, 4d, 5d, and 6d interpolation were also created (these are simply
!  extensions of the same algorithm into higher dimensions).
!
!### See also
!  * An object-oriented interface can be found in [[bspline_oo_module]].
!
!### References
!
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
!  6. D.E. Amos, "Quadrature subroutines for splines and B-splines",
!     Report SAND79-1825, Sandia Laboratories, December 1979.

    module bspline_sub_module

    use,intrinsic :: iso_fortran_env, only: real64
    use,intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    integer,parameter :: wp = real64  !! Real precision

    abstract interface
        function int_func(x) result(f)
        !! interface for the input function in [[dbfqad]]
        import :: wp
        real(wp),intent(in) :: x
        real(wp)            :: f  !! f(x)
        end function int_func
    end interface

    !Spline function order (order = polynomial degree + 1)
    integer,parameter,public :: bspline_order_quadratic = 3
    integer,parameter,public :: bspline_order_cubic     = 4
    integer,parameter,public :: bspline_order_quartic   = 5
    integer,parameter,public :: bspline_order_quintic   = 6

    !main routines:
    public :: db1ink, db1val, db1qad
    public :: db2ink, db2val
    public :: db3ink, db3val
    public :: db4ink, db4val
    public :: db5ink, db5val
    public :: db6ink, db6val

    public :: get_status_message

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
!  the one-dimensional gridded data
!  $$ [x(i),\mathrm{fcn}(i)] ~\mathrm{for}~ i=1,..,n_x $$
!  The interpolating function and its derivatives may
!  subsequently be evaluated by the function [[db1val]].
!
!### History
!  * Jacob Williams, 10/30/2015 : Created 1D routine.

    pure subroutine db1ink(x,nx,fcn,kx,iknot,tx,bcoef,iflag)

    implicit none

    integer,intent(in)                      :: nx     !! Number of \(x\) abcissae
    integer,intent(in)                      :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)        :: fcn    !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                                      !! contain the function value at the point `x(i)`
    integer,intent(in)                      :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant:
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db1ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(out)       :: bcoef  !! `(nx)` array of coefficients of the b-spline interpolant.
    integer,intent(out)                     :: iflag  !! status flag:
                                                      !!
                                                      !! * 0 = successful execution.
                                                      !! * 2 = `iknot` out of range.
                                                      !! * 3 = `nx` out of range.
                                                      !! * 4 = `kx` out of range.
                                                      !! * 5 = `x` not strictly increasing.
                                                      !! * 6 = `tx` not non-decreasing.
                                                      !! * 700 = `size(x)` \( \ne \) `size(fcn,1)`.
                                                      !! * 706 = `size(x)` \( \ne \) `nx`.
                                                      !! * 712 = `size(tx)` \( \ne \) `nx+kx`.
                                                      !! * 800 = `size(x)` \( \ne \) `size(bcoef,1)`.

    real(wp),dimension(2*kx*(nx+1)) :: work
    logical :: status_ok

    !check validity of inputs

    call check_inputs('db1ink',&
                        iknot,&
                        iflag,&
                        nx=nx,&
                        kx=kx,&
                        x=x,&
                        f1=fcn,&
                        bcoef1=bcoef,&
                        tx=tx,&
                        status_ok=status_ok)

    if (status_ok) then

        !choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
        end if

        !construct b-spline coefficients

        call dbtpcf(x,nx,fcn,nx,1,tx,kx,bcoef,work,iflag)

    end if

    end subroutine db1ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db1ink]] or one of its
!  derivatives at the point `xval`.
!
!  To evaluate the interpolant itself, set `idx=0`,
!  to evaluate the first partial with respect to `x`, set `idx=1`, and so on.
!
!  [[db1val]] returns 0.0 if (`xval`,`yval`) is out of range. that is, if
!```fortran
!   xval < tx(1) .or. xval > tx(nx+kx)
!```
!  if the knots `tx` were chosen by [[db1ink]], then this is equivalent to:
!```fortran
!   xval < x(1) .or. xval > x(nx)+epsx
!```
!  where
!```fortran
!   epsx = 0.1*(x(nx)-x(nx-1))
!```
!
!  The input quantities `tx`, `nx`, `kx`, and `bcoef` should be
!  unchanged since the last call of [[db1ink]].
!
!### History
!  * Jacob Williams, 10/30/2015 : Created 1D routine.

    pure subroutine db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,extrap)

    implicit none

    integer,intent(in)                   :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    integer,intent(in)                   :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction. (same as in last call to [[db1ink]])
    real(wp),dimension(nx),intent(in)    :: bcoef    !! the b-spline coefficients computed by [[db1ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer,intent(out)                  :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer,intent(inout)                :: inbvx    !! initialization parameter which must be set
                                                     !! to 1 the first time this routine is called,
                                                     !! and must not be changed by the user.
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    real(wp),dimension(3*kx) :: work

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return

    call dbvalu(tx,bcoef,nx,kx,idx,xval,inbvx,work,iflag,f,extrap)

    end subroutine db1val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the integral on `(x1,x2)` of a `kx`-th order b-spline.
!  Orders `kx` as high as 20 are permitted by applying a 2, 6, or 10
!  point gauss formula on subintervals of `(x1,x2)` which are
!  formed by included (distinct) knots.
!
!### See also
!  * [[dbsqad]] -- the core routine.

    pure subroutine db1qad(tx,bcoef,nx,kx,x1,x2,f,iflag)

    implicit none

    integer,intent(in)                   :: nx      !! length of coefficient array
    integer,intent(in)                   :: kx      !! order of b-spline, `1 <= k <= 20`
    real(wp),dimension(nx+kx),intent(in) :: tx      !! knot array
    real(wp),dimension(nx),intent(in)    :: bcoef   !! b-spline coefficient array
    real(wp),intent(in)                  :: x1      !! end point of quadrature interval in `t(kx) <= x <= t(nx+1)`
    real(wp),intent(in)                  :: x2      !! end point of quadrature interval in `t(kx) <= x <= t(nx+1)`
    real(wp),intent(out)                 :: f       !! integral of the b-spline over (`x1`,`x2`)
    integer,intent(out)                  :: iflag   !! status flag:
                                                    !!
                                                    !! * \( = 0 \)   : no errors
                                                    !! * \( \ne 0 \) : error

    real(wp),dimension(3*kx) :: work !! work array for [[dbsqad]]

    call dbsqad(tx,bcoef,nx,kx,x1,x2,f,work,iflag)

    end subroutine db1qad
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
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
!  Note that for each fixed value of \(y\), \( s(x,y) \) is a piecewise
!  polynomial function of \(x\) alone, and for each fixed value of \(x\), \( s(x,y) \)
!  is a piecewise polynomial function of \(y\) alone. in one dimension
!  a piecewise polynomial may be created by partitioning a given
!  interval into subintervals and defining a distinct polynomial piece
!  on each one. the points where adjacent subintervals meet are called
!  knots. each of the functions \(u_i\) and \(v_j\) above is a piecewise
!  polynomial.
!
!  Users of [[db2ink]] choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the \(x\) and
!  \(y\) directions (`kx` and `ky`). users also may define their own knot
!  sequence in \(x\) and \(y\) separately (`tx` and `ty`). if `iflag=0`, however,
!  [[db2ink]] will choose sequences of knots that result in a piecewise
!  polynomial interpolant with `kx-2` continuous partial derivatives in
!  \(x\) and `ky-2` continuous partial derivatives in \(y\). (`kx` knots are taken
!  near each endpoint in the \(x\) direction, not-a-knot end conditions
!  are used, and the remaining knots are placed at data points if `kx`
!  is even or at midpoints between data points if `kx` is odd. the \(y\)
!  direction is treated similarly.)
!
!  After a call to [[db2ink]], all information necessary to define the
!  interpolating function are contained in the parameters `nx`, `ny`, `kx`,
!  `ky`, `tx`, `ty`, and `bcoef`. These quantities should not be altered until
!  after the last call of the evaluation routine [[db2val]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2ink(x,nx,y,ny,fcn,kx,ky,iknot,tx,ty,bcoef,iflag)

    implicit none

    integer,intent(in)                      :: nx     !! Number of \(x\) abcissae
    integer,intent(in)                      :: ny     !! Number of \(y\) abcissae
    integer,intent(in)                      :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                      :: ky     !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)        :: y      !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in)      :: fcn    !! `(nx,ny)` matrix of function values to interpolate.
                                                      !! `fcn(i,j)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`)
    integer,intent(in)                      :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)     :: ty     !! The `(ny+ky)` knots in the \(y\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:,:),intent(out)     :: bcoef  !! `(nx,ny)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                     :: iflag  !! *  0 = successful execution.
                                                      !! *  2 = `iknot` out of range.
                                                      !! *  3 = `nx` out of range.
                                                      !! *  4 = `kx` out of range.
                                                      !! *  5 = `x` not strictly increasing.
                                                      !! *  6 = `tx` not non-decreasing.
                                                      !! *  7 = `ny` out of range.
                                                      !! *  8 = `ky` out of range.
                                                      !! *  9 = `y` not strictly increasing.
                                                      !! * 10 = `ty` not non-decreasing.
                                                      !! * 700 = `size(x)`  \( \ne \) `size(fcn,1)`
                                                      !! * 701 = `size(y)`  \( \ne \) `size(fcn,2)`
                                                      !! * 706 = `size(x)`  \( \ne \) `nx`
                                                      !! * 707 = `size(y)`  \( \ne \) `ny`
                                                      !! * 712 = `size(tx)` \( \ne \) `nx+kx`
                                                      !! * 713 = `size(ty)` \( \ne \) `ny+ky`
                                                      !! * 800 = `size(x)`  \( \ne \) `size(bcoef,1)`
                                                      !! * 801 = `size(y)`  \( \ne \) `size(bcoef,2)`

    real(wp),dimension(nx*ny) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1))) :: work
    logical :: status_ok

    !check validity of inputs

    call check_inputs('db2ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,&
                        kx=kx,ky=ky,&
                        x=x,y=y,&
                        tx=tx,ty=ty,&
                        f2=fcn,&
                        bcoef2=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        !choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
        end if

        !construct b-spline coefficients

                      call dbtpcf(x,nx,fcn, nx,ny,tx,kx,temp, work,iflag)
        if (iflag==0) call dbtpcf(y,ny,temp,ny,nx,ty,ky,bcoef,work,iflag)

    end if

    end subroutine db2ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db2ink]] or one of its
!  derivatives at the point (`xval`,`yval`).
!
!  To evaluate the interpolant
!  itself, set `idx=idy=0`, to evaluate the first partial with respect
!  to `x`, set `idx=1,idy=0`, and so on.
!
!  [[db2val]] returns 0.0 if `(xval,yval)` is out of range. that is, if
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
!  The input quantities `tx`, `ty`, `nx`, `ny`, `kx`, `ky`, and `bcoef` should be
!  unchanged since the last call of [[db2ink]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,extrap)

    implicit none

    integer,intent(in)                   :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: ny       !! the number of interpolation points in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: ky       !! order of polynomial pieces in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                  :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(wp),dimension(ny+ky),intent(in) :: ty       !! sequence of knots defining the piecewise
                                                     !! polynomial in the \(y\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(wp),dimension(nx,ny),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db2ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer,intent(out)                  :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer,intent(inout)                :: inbvx    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer,intent(inout)                :: inbvy    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer,intent(inout)                :: iloy     !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    integer :: k, lefty, mflag, kcol
    real(wp),dimension(ky) :: temp
    real(wp),dimension(3*max(kx,ky)) :: work

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return

    kcol = lefty - ky
    do k=1,ky
        kcol = kcol + 1
        call dbvalu(tx,bcoef(:,kcol),nx,kx,idx,xval,inbvx,work,iflag,temp(k),extrap)
        if (iflag/=0) return !error
    end do

    kcol = lefty - ky + 1
    call dbvalu(ty(kcol:),temp,ky,ky,idy,yval,inbvy,work,iflag,f,extrap)

    end subroutine db2val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
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
!  Note that for fixed values of \(y\) and \(z\) \(s(x,y,z)\) is a piecewise
!  polynomial function of \(x\) alone, for fixed values of \(x\) and \(z\) \(s(x,y,z)\)
!  is a piecewise polynomial function of \(y\) alone, and for fixed
!  values of \(x\) and \(y\) \(s(x,y,z)\) is a function of \(z\) alone. in one
!  dimension a piecewise polynomial may be created by partitioning a
!  given interval into subintervals and defining a distinct polynomial
!  piece on each one. the points where adjacent subintervals meet are
!  called knots. each of the functions \(u_i\), \(v_j\), and \(w_k\) above is a
!  piecewise polynomial.
!
!  Users of [[db3ink]] choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the \(x\), \(y\),
!  and \(z\) directions (`kx`, `ky`, and `kz`). users also may define their own
!  knot sequence in \(x\), \(y\), \(z\) separately (`tx`, `ty`, and `tz`). if `iflag=0`,
!  however, [[db3ink]] will choose sequences of knots that result in a
!  piecewise polynomial interpolant with `kx-2` continuous partial
!  derivatives in \(x\), `ky-2` continuous partial derivatives in \(y\), and `kz-2`
!  continuous partial derivatives in \(z\). (`kx` knots are taken near
!  each endpoint in \(x\), not-a-knot end conditions are used, and the
!  remaining knots are placed at data points if `kx` is even or at
!  midpoints between data points if `kx` is odd. the \(y\) and \(z\) directions
!  are treated similarly.)
!
!  After a call to [[db3ink]], all information necessary to define the
!  interpolating function are contained in the parameters `nx`, `ny`, `nz`,
!  `kx`, `ky`, `kz`, `tx`, `ty`, `tz`, and `bcoef`. these quantities should not be
!  altered until after the last call of the evaluation routine [[db3val]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,iknot,tx,ty,tz,bcoef,iflag)

    implicit none

    integer,intent(in)                       :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                       :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                       :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                       :: kx    !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                       :: ky    !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                       :: kz    !! the order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)         :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in)     :: fcn   !! `(nx,ny,nz)` matrix of function values to interpolate. `fcn(i,j,k)` should
                                                      !! contain the function value at the point (`x(i)`,`y(j)`,`z(k)`)
    integer,intent(in)                       :: iknot !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db3ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)      :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:,:,:),intent(out)    :: bcoef !! `(nx,ny,nz)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                      :: iflag !! *  0 = successful execution.
                                                      !! *  2 = `iknot` out of range.
                                                      !! *  3 = `nx` out of range.
                                                      !! *  4 = `kx` out of range.
                                                      !! *  5 = `x` not strictly increasing.
                                                      !! *  6 = `tx` not non-decreasing.
                                                      !! *  7 = `ny` out of range.
                                                      !! *  8 = `ky` out of range.
                                                      !! *  9 = `y` not strictly increasing.
                                                      !! * 10 = `ty` not non-decreasing.
                                                      !! * 11 = `nz` out of range.
                                                      !! * 12 = `kz` out of range.
                                                      !! * 13 = `z` not strictly increasing.
                                                      !! * 14 = `ty` not non-decreasing.
                                                      !! * 700 = `size(x) ` \(\ne\) `size(fcn,1)`
                                                      !! * 701 = `size(y) ` \(\ne\) `size(fcn,2)`
                                                      !! * 702 = `size(z) ` \(\ne\) `size(fcn,3)`
                                                      !! * 706 = `size(x) ` \(\ne\) `nx`
                                                      !! * 707 = `size(y) ` \(\ne\) `ny`
                                                      !! * 708 = `size(z) ` \(\ne\) `nz`
                                                      !! * 712 = `size(tx)` \(\ne\) `nx+kx`
                                                      !! * 713 = `size(ty)` \(\ne\) `ny+ky`
                                                      !! * 714 = `size(tz)` \(\ne\) `nz+kz`
                                                      !! * 800 = `size(x) ` \(\ne\) `size(bcoef,1)`
                                                      !! * 801 = `size(y) ` \(\ne\) `size(bcoef,2)`
                                                      !! * 802 = `size(z) ` \(\ne\) `size(bcoef,3)`

    real(wp),dimension(nx*ny*nz) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))) :: work
    logical :: status_ok

    ! check validity of input

    call check_inputs('db3ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,&
                        kx=kx,ky=ky,kz=kz,&
                        x=x,y=y,z=z,&
                        tx=tx,ty=ty,tz=tz,&
                        f3=fcn,&
                        bcoef3=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        ! choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
        end if

        ! copy fcn to work in packed for dbtpcf
        temp(1:nx*ny*nz) = reshape( fcn, [nx*ny*nz] )

        ! construct b-spline coefficients

                      call dbtpcf(x,nx,temp, nx,ny*nz,tx,kx,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(y,ny,bcoef,ny,nx*nz,ty,ky,temp, work,iflag)
        if (iflag==0) call dbtpcf(z,nz,temp, nz,nx*ny,tz,kz,bcoef,work,iflag)

    end if

    end subroutine db3ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db3ink]] or one of its
!  derivatives at the point (`xval`,`yval`,`zval`).
!
!  To evaluate the
!  interpolant itself, set `idx=idy=idz=0`, to evaluate the first
!  partial with respect to `x`, set `idx=1`,`idy=idz=0`, and so on.
!
!  [[db3val]] returns 0.0 if (`xval`,`yval`,`zval`) is out of range. that is,
!```fortran
! xval<tx(1) .or. xval>tx(nx+kx) .or.
! yval<ty(1) .or. yval>ty(ny+ky) .or.
! zval<tz(1) .or. zval>tz(nz+kz)
!```
!  if the knots `tx`, `ty`, and `tz` were chosen by [[db3ink]], then this is
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
!  The input quantities `tx`, `ty`, `tz`, `nx`, `ny`, `nz`, `kx`, `ky`, `kz`, and `bcoef`
!  should remain unchanged since the last call of [[db3ink]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db3val(xval,yval,zval,idx,idy,idz,&
                                     tx,ty,tz,&
                                     nx,ny,nz,kx,ky,kz,bcoef,f,iflag,&
                                     inbvx,inbvy,inbvz,iloy,iloz,extrap)

    implicit none

    integer,intent(in)                      :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: nx       !! the number of interpolation points in \(x\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ny       !! the number of interpolation points in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: nz       !! the number of interpolation points in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kx       !! order of polynomial pieces in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ky       !! order of polynomial pieces in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kz       !! order of polynomial pieces in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    real(wp),intent(in)                     :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                     :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)                     :: zval     !! \(z\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)    :: tx       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(x\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(ny+ky),intent(in)    :: ty       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(y\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nz+kz),intent(in)    :: tz       !! sequence of knots defining the piecewise polynomial
                                                        !! in the \(z\) direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nx,ny,nz),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db3ink]].
    real(wp),intent(out)                    :: f        !! interpolated value
    integer,intent(out)                     :: iflag    !! status flag:
                                                        !!
                                                        !! * \( = 0 \)   : no errors
                                                        !! * \( \ne 0 \) : error
    integer,intent(inout)                   :: inbvx    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: inbvy    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: inbvz    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: iloy     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer,intent(inout)                   :: iloz     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    logical,intent(in),optional             :: extrap   !! if extrapolation is allowed
                                                        !! (if not present, default is False)

    real(wp),dimension(ky,kz)           :: temp1
    real(wp),dimension(kz)              :: temp2
    real(wp),dimension(3*max(kx,ky,kz)) :: work

    integer :: lefty, leftz, mflag,&
                kcoly, kcolz, j, k

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return
    iflag = check_value(zval,tz,3,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag,extrap); if (mflag /= 0) return

    iflag = 0

    kcolz = leftz - kz
    do k=1,kz
        kcolz = kcolz + 1
        kcoly = lefty - ky
        do j=1,ky
            kcoly = kcoly + 1
            call dbvalu(tx,bcoef(:,kcoly,kcolz),nx,kx,idx,xval,inbvx,work,iflag,temp1(j,k),extrap)
            if (iflag/=0) return
        end do
    end do

    kcoly = lefty - ky + 1
    do k=1,kz
        call dbvalu(ty(kcoly:),temp1(:,k),ky,ky,idy,yval,inbvy,work,iflag,temp2(k),extrap)
        if (iflag/=0) return
    end do

    kcolz = leftz - kz + 1
    call dbvalu(tz(kcolz:),temp2,kz,kz,idz,zval,inbvz,work,iflag,f,extrap)

    end subroutine db3val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
!  the four-dimensional gridded data
!  $$ [x(i),y(j),z(k),q(l),\mathrm{fcn}(i,j,k,l)] ~\mathrm{for}~
!     i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y, ~\mathrm{and}~ k=1,..,n_z,
!     ~\mathrm{and}~ l=1,..,n_q $$
!  The interpolating function and its derivatives may
!  subsequently be evaluated by the function [[db4val]].
!
!  See [[db3ink]] header for more details.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db4ink(x,nx,y,ny,z,nz,q,nq,&
                        fcn,&
                        kx,ky,kz,kq,&
                        iknot,&
                        tx,ty,tz,tq,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                          :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                          :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                          :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                          :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                          :: kx    !! the order of spline pieces in \(x\)
                                                         !! ( \( 2 \le k_x < n_x \) ).
                                                         !! (order = polynomial degree + 1)
    integer,intent(in)                          :: ky    !! the order of spline pieces in \(y\)
                                                         !! ( \( 2 \le k_y < n_y \) ).
                                                         !! (order = polynomial degree + 1)
    integer,intent(in)                          :: kz    !! the order of spline pieces in \(z\)
                                                         !! ( \( 2 \le k_z < n_z \) ).
                                                         !! (order = polynomial degree + 1)
    integer,intent(in)                          :: kq    !! the order of spline pieces in \(q\)
                                                         !! ( \( 2 \le k_q < n_q \) ).
                                                         !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)            :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: q     !! `(nq)` array of \(q\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)      :: fcn   !! `(nx,ny,nz,nq)` matrix of function values to interpolate. `fcn(i,j,k,q)` should
                                                         !! contain the function value at the point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer,intent(in)                          :: iknot !! knot sequence flag:
                                                         !!
                                                         !! * 0 = knot sequence chosen by [[db4ink]].
                                                         !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)         :: tx    !! The `(nx+kx)` knots in the x direction for the spline interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: ty    !! The `(ny+ky)` knots in the y direction for the spline interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: tz    !! The `(nz+kz)` knots in the z direction for the spline interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: tq    !! The `(nq+kq)` knots in the q direction for the spline interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:),intent(out)     :: bcoef !! `(nx,ny,nz,nq)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                         :: iflag !! *  0 = successful execution.
                                                         !! *  2 = `iknot` out of range.
                                                         !! *  3 = `nx` out of range.
                                                         !! *  4 = `kx` out of range.
                                                         !! *  5 = `x` not strictly increasing.
                                                         !! *  6 = `tx` not non-decreasing.
                                                         !! *  7 = `ny` out of range.
                                                         !! *  8 = `ky` out of range.
                                                         !! *  9 = `y` not strictly increasing.
                                                         !! * 10 = `ty` not non-decreasing.
                                                         !! * 11 = `nz` out of range.
                                                         !! * 12 = `kz` out of range.
                                                         !! * 13 = `z` not strictly increasing.
                                                         !! * 14 = `tz` not non-decreasing.
                                                         !! * 15 = `nq` out of range.
                                                         !! * 16 = `kq` out of range.
                                                         !! * 17 = `q` not strictly increasing.
                                                         !! * 18 = `tq` not non-decreasing.
                                                         !! * 700 = `size(x)`  \( \ne \) `size(fcn,1)`
                                                         !! * 701 = `size(y)`  \( \ne \) `size(fcn,2)`
                                                         !! * 702 = `size(z)`  \( \ne \) `size(fcn,3)`
                                                         !! * 703 = `size(q)`  \( \ne \) `size(fcn,4)`
                                                         !! * 706 = `size(x)`  \( \ne \) `nx`
                                                         !! * 707 = `size(y)`  \( \ne \) `ny`
                                                         !! * 708 = `size(z)`  \( \ne \) `nz`
                                                         !! * 709 = `size(q)`  \( \ne \) `nq`
                                                         !! * 712 = `size(tx`) \( \ne \) `nx+kx`
                                                         !! * 713 = `size(ty`) \( \ne \) `ny+ky`
                                                         !! * 714 = `size(tz`) \( \ne \) `nz+kz`
                                                         !! * 715 = `size(tq`) \( \ne \) `nq+kq`
                                                         !! * 800 = `size(x)`  \( \ne \) `size(bcoef,1)`
                                                         !! * 801 = `size(y)`  \( \ne \) `size(bcoef,2)`
                                                         !! * 802 = `size(z)`  \( \ne \) `size(bcoef,3)`
                                                         !! * 803 = `size(q)`  \( \ne \) `size(bcoef,4)`

    real(wp),dimension(nx*ny*nz*nq) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1),2*kq*(nq+1))) :: work
    logical :: status_ok

    ! check validity of input

    call check_inputs('db4ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,&
                        kx=kx,ky=ky,kz=kz,kq=kq,&
                        x=x,y=y,z=z,q=q,&
                        tx=tx,ty=ty,tz=tz,tq=tq,&
                        f4=fcn,&
                        bcoef4=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        ! choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
        end if

        ! construct b-spline coefficients

                      call dbtpcf(x,nx,fcn,  nx,ny*nz*nq,tx,kx,temp, work,iflag)
        if (iflag==0) call dbtpcf(y,ny,temp, ny,nx*nz*nq,ty,ky,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(z,nz,bcoef,nz,nx*ny*nq,tz,kz,temp, work,iflag)
        if (iflag==0) call dbtpcf(q,nq,temp, nq,nx*ny*nz,tq,kq,bcoef,work,iflag)

     end if

    end subroutine db4ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db4ink]] or one of its
!  derivatives at the point (`xval`,`yval`,`zval`,`qval`).
!
!  To evaluate the
!  interpolant itself, set `idx=idy=idz=idq=0`, to evaluate the first
!  partial with respect to `x`, set `idx=1,idy=idz=idq=0`, and so on.
!
!  See [[db3val]] header for more information.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db4val(xval,yval,zval,qval,&
                                idx,idy,idz,idq,&
                                tx,ty,tz,tq,&
                                nx,ny,nz,nq,&
                                kx,ky,kz,kq,&
                                bcoef,f,iflag,&
                                inbvx,inbvy,inbvz,inbvq,&
                                iloy,iloz,iloq,extrap)

    implicit none

    integer,intent(in)                         :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                         :: nx       !! the number of interpolation points in \(x\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: ny       !! the number of interpolation points in \(y\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: nz       !! the number of interpolation points in \(z\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: nq       !! the number of interpolation points in \(q\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: kx       !! order of polynomial pieces in \(x\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: ky       !! order of polynomial pieces in \(y\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: kz       !! order of polynomial pieces in \(z\).
                                                           !! (same as in last call to [[db4ink]])
    integer,intent(in)                         :: kq       !! order of polynomial pieces in \(q\).
                                                           !! (same as in last call to [[db4ink]])
    real(wp),intent(in)                        :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                        :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)                        :: zval     !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)                        :: qval     !! \(q\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)       :: tx       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(x\) direction. (same as in last call to [[db4ink]])
    real(wp),dimension(ny+ky),intent(in)       :: ty       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(y\) direction. (same as in last call to [[db4ink]])
    real(wp),dimension(nz+kz),intent(in)       :: tz       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(z\) direction. (same as in last call to [[db4ink]])
    real(wp),dimension(nq+kq),intent(in)       :: tq       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(q\) direction. (same as in last call to [[db4ink]])
    real(wp),dimension(nx,ny,nz,nq),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db4ink]].
    real(wp),intent(out)                       :: f        !! interpolated value
    integer,intent(out)                        :: iflag    !! status flag:
                                                           !!
                                                           !! * \( = 0 \)   : no errors
                                                           !! * \( \ne 0 \) : error
    integer,intent(inout)                      :: inbvx    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: inbvy    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: inbvz    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: inbvq    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: iloy     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: iloz     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer,intent(inout)                      :: iloq     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    logical,intent(in),optional                :: extrap   !! if extrapolation is allowed
                                                           !! (if not present, default is False)

    real(wp),dimension(ky,kz,kq)             :: temp1
    real(wp),dimension(kz,kq)                :: temp2
    real(wp),dimension(kq)                   :: temp3
    real(wp),dimension(3*max(kx,ky,kz,kq))   :: work
    integer :: lefty, leftz, leftq, mflag,&
                kcoly, kcolz, kcolq, j, k, q

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return
    iflag = check_value(zval,tz,3,extrap); if (iflag/=0) return
    iflag = check_value(qval,tq,4,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag,extrap); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag,extrap); if (mflag /= 0) return

    iflag = 0

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
                call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq),&
                                     nx,kx,idx,xval,inbvx,work,iflag,&
                                     temp1(j,k,q),extrap)
                if (iflag/=0) return
            end do
        end do
    end do

    ! y -> z, q
    kcoly = lefty - ky + 1
    do q=1,kq
        do k=1,kz
            call dbvalu(ty(kcoly:),temp1(:,k,q),&
                        ky,ky,idy,yval,inbvy,work,iflag,&
                        temp2(k,q),extrap)
            if (iflag/=0) return
        end do
    end do

    ! z -> q
    kcolz = leftz - kz + 1
    do q=1,kq
        call dbvalu(tz(kcolz:),temp2(:,q),&
                    kz,kz,idz,zval,inbvz,work,iflag,&
                    temp3(q),extrap)
        if (iflag/=0) return
    end do

    ! q
    kcolq = leftq - kq + 1
    call dbvalu(tq(kcolq:),temp3,kq,kq,idq,qval,inbvq,work,iflag,f,extrap)

    end subroutine db4val
!*****************************************************************************************

!*****************************************************************************************
!>
! Determines the parameters of a function that interpolates
! the five-dimensional gridded data:
!
! $$ [x(i),y(j),z(k),q(l),r(m),\mathrm{fcn}(i,j,k,l,m)] $$
!
! for:
!
! $$ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y, ~\mathrm{and}~ k=1,..,n_z,
!    ~\mathrm{and}~ l=1,..,n_q, ~\mathrm{and}~ m=1,..,n_r $$
!
! The interpolating function and its derivatives may subsequently be evaluated
! by the function [[db5val]].
!
! See [[db3ink]] header for more details.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                        fcn,&
                        kx,ky,kz,kq,kr,&
                        iknot,&
                        tx,ty,tz,tq,tr,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                             :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                             :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                             :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                             :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                             :: nr    !! number of \(r\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                             :: kx    !! the order of spline pieces in \(x\)
                                                            !! ( \( 2 \le k_x < n_x \) ).
                                                            !! (order = polynomial degree + 1)
    integer,intent(in)                             :: ky    !! the order of spline pieces in \(y\)
                                                            !! ( \( 2 \le k_y < n_y \) ).
                                                            !! (order = polynomial degree + 1)
    integer,intent(in)                             :: kz    !! the order of spline pieces in \(z\)
                                                            !! ( \( 2 \le k_z < n_z \) ).
                                                            !! (order = polynomial degree + 1)
    integer,intent(in)                             :: kq    !! the order of spline pieces in \(q\)
                                                            !! ( \( 2 \le k_q < n_q \) ).
                                                            !! (order = polynomial degree + 1)
    integer,intent(in)                             :: kr    !! the order of spline pieces in \(r\)
                                                            !! ( \( 2 \le k_r < n_r \) ).
                                                            !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)               :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: q     !! `(nq)` array of \(q\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: r     !! `(nr)` array of \(r\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)       :: fcn   !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate. `fcn(i,j,k,q,r)` should
                                                            !! contain the function value at the point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer,intent(in)                             :: iknot !! knot sequence flag:
                                                            !!
                                                            !! * 0 = knot sequence chosen by [[db5ink]].
                                                            !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)            :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the spline interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the spline interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the spline interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tq    !! The `(nq+kq)` knots in the \(q\) direction for the spline interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tr    !! The `(nr+kr)` knots in the \(r\) direction for the spline interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:,:),intent(out)      :: bcoef !! `(nx,ny,nz,nq,nr)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                            :: iflag !! *  0 = successful execution.
                                                            !! *  2 = `iknot` out of range.
                                                            !! *  3 = `nx` out of range.
                                                            !! *  4 = `kx` out of range.
                                                            !! *  5 = `x` not strictly increasing.
                                                            !! *  6 = `tx` not non-decreasing.
                                                            !! *  7 = `ny` out of range.
                                                            !! *  8 = `ky` out of range.
                                                            !! *  9 = `y` not strictly increasing.
                                                            !! * 10 = `ty` not non-decreasing.
                                                            !! * 11 = `nz` out of range.
                                                            !! * 12 = `kz` out of range.
                                                            !! * 13 = `z` not strictly increasing.
                                                            !! * 14 = `tz` not non-decreasing.
                                                            !! * 15 = `nq` out of range.
                                                            !! * 16 = `kq` out of range.
                                                            !! * 17 = `q` not strictly increasing.
                                                            !! * 18 = `tq` not non-decreasing.
                                                            !! * 19 = `nr` out of range.
                                                            !! * 20 = `kr` out of range.
                                                            !! * 21 = `r` not strictly increasing.
                                                            !! * 22 = `tr` not non-decreasing.
                                                            !! * 700 = `size(x)`  \( \ne \) `size(fcn,1)`
                                                            !! * 701 = `size(y)`  \( \ne \) `size(fcn,2)`
                                                            !! * 702 = `size(z)`  \( \ne \) `size(fcn,3)`
                                                            !! * 703 = `size(q)`  \( \ne \) `size(fcn,4)`
                                                            !! * 704 = `size(r)`  \( \ne \) `size(fcn,5)`
                                                            !! * 706 = `size(x)`  \( \ne \) `nx`
                                                            !! * 707 = `size(y)`  \( \ne \) `ny`
                                                            !! * 708 = `size(z)`  \( \ne \) `nz`
                                                            !! * 709 = `size(q)`  \( \ne \) `nq`
                                                            !! * 710 = `size(r)`  \( \ne \) `nr`
                                                            !! * 712 = `size(tx)` \( \ne \) `nx+kx`
                                                            !! * 713 = `size(ty)` \( \ne \) `ny+ky`
                                                            !! * 714 = `size(tz)` \( \ne \) `nz+kz`
                                                            !! * 715 = `size(tq)` \( \ne \) `nq+kq`
                                                            !! * 716 = `size(tr)` \( \ne \) `nr+kr`
                                                            !! * 800 = `size(x)`  \( \ne \) `size(bcoef,1)`
                                                            !! * 801 = `size(y)`  \( \ne \) `size(bcoef,2)`
                                                            !! * 802 = `size(z)`  \( \ne \) `size(bcoef,3)`
                                                            !! * 803 = `size(q)`  \( \ne \) `size(bcoef,4)`
                                                            !! * 804 = `size(r)`  \( \ne \) `size(bcoef,5)`

    real(wp),dimension(nx*ny*nz*nq*nr) :: temp
    real(wp),dimension(max( 2*kx*(nx+1),&
                            2*ky*(ny+1),&
                            2*kz*(nz+1),&
                            2*kq*(nq+1),&
                            2*kr*(nr+1) )) :: work
    logical :: status_ok

    !  check validity of input

    call check_inputs('db5ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,nr=nr,&
                        kx=kx,ky=ky,kz=kz,kq=kq,kr=kr,&
                        x=x,y=y,z=z,q=q,r=r,&
                        tx=tx,ty=ty,tz=tz,tq=tq,tr=tr,&
                        f5=fcn,&
                        bcoef5=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        !  choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
        end if

        ! copy fcn to work in packed for dbtpcf

        temp(1:nx*ny*nz*nq*nr) = reshape( fcn, [nx*ny*nz*nq*nr] )

        !  construct b-spline coefficients

                      call dbtpcf(x,nx,temp,  nx,ny*nz*nq*nr,tx,kx,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(y,ny,bcoef, ny,nx*nz*nq*nr,ty,ky,temp, work,iflag)
        if (iflag==0) call dbtpcf(z,nz,temp,  nz,nx*ny*nq*nr,tz,kz,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(q,nq,bcoef, nq,nx*ny*nz*nr,tq,kq,temp, work,iflag)
        if (iflag==0) call dbtpcf(r,nr,temp,  nr,nx*ny*nz*nq,tr,kr,bcoef,work,iflag)

     end if

    end subroutine db5ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db5ink]] or one of its
!  derivatives at the point (`xval`,`yval`,`zval`,`qval`,`rval`).
!
!  To evaluate the
!  interpolant itself, set `idx=idy=idz=idq=idr=0`, to evaluate the first
!  partial with respect to `x`, set `idx=1,idy=idz=idq=idr=0,` and so on.
!
!  See [[db3val]] header for more information.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db5val(xval,yval,zval,qval,rval,&
                                idx,idy,idz,idq,idr,&
                                tx,ty,tz,tq,tr,&
                                nx,ny,nz,nq,nr,&
                                kx,ky,kz,kq,kr,&
                                bcoef,f,iflag,&
                                inbvx,inbvy,inbvz,inbvq,inbvr,&
                                iloy,iloz,iloq,ilor,extrap)

    implicit none

    integer,intent(in)                            :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: idr      !! \(r\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                            :: nx       !! the number of interpolation points in \(x\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: ny       !! the number of interpolation points in \(y\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: nz       !! the number of interpolation points in \(z\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: nq       !! the number of interpolation points in \(q\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: nr       !! the number of interpolation points in \(r\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: kx       !! order of polynomial pieces in \(x\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: ky       !! order of polynomial pieces in \(y\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: kz       !! order of polynomial pieces in \(z\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: kq       !! order of polynomial pieces in \(q\).
                                                              !! (same as in last call to [[db5ink]])
    integer,intent(in)                            :: kr       !! order of polynomial pieces in \(r\).
                                                              !! (same as in last call to [[db5ink]])
    real(wp),intent(in)                           :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                           :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)                           :: zval     !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)                           :: qval     !! \(q\) coordinate of evaluation point.
    real(wp),intent(in)                           :: rval     !! \(r\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)          :: tx       !! sequence of knots defining the piecewise polynomial
                                                              !! in the \(x\) direction.
                                                              !! (same as in last call to [[db5ink]])
    real(wp),dimension(ny+ky),intent(in)          :: ty       !! sequence of knots defining the piecewise polynomial
                                                              !! in the \(y\) direction.
                                                              !! (same as in last call to [[db5ink]])
    real(wp),dimension(nz+kz),intent(in)          :: tz       !! sequence of knots defining the piecewise polynomial
                                                              !! in the \(z\) direction.
                                                              !! (same as in last call to [[db5ink]])
    real(wp),dimension(nq+kq),intent(in)          :: tq       !! sequence of knots defining the piecewise polynomial
                                                              !! in the \(q\) direction.
                                                              !! (same as in last call to [[db5ink]])
    real(wp),dimension(nr+kr),intent(in)          :: tr       !! sequence of knots defining the piecewise polynomial
                                                              !! in the \(r\) direction.
                                                              !! (same as in last call to [[db5ink]])
    real(wp),dimension(nx,ny,nz,nq,nr),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db5ink]].
    real(wp),intent(out)                          :: f        !! interpolated value
    integer,intent(out)                           :: iflag    !! status flag:
                                                              !!
                                                              !! * \( = 0 \)   : no errors
                                                              !! * \( \ne 0 \) : error
    integer,intent(inout)                         :: inbvx    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: inbvy    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: inbvz    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: inbvq    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: inbvr    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: iloy     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: iloz     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: iloq     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer,intent(inout)                         :: ilor     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    logical,intent(in),optional                   :: extrap   !! if extrapolation is allowed
                                                              !! (if not present, default is False)

    real(wp),dimension(ky,kz,kq,kr)           :: temp1
    real(wp),dimension(kz,kq,kr)              :: temp2
    real(wp),dimension(kq,kr)                 :: temp3
    real(wp),dimension(kr)                    :: temp4
    real(wp),dimension(3*max(kx,ky,kz,kq,kr)) :: work
    integer :: lefty, leftz, leftq, leftr, mflag,&
               kcoly, kcolz, kcolq, kcolr, j, k, q, r

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return
    iflag = check_value(zval,tz,3,extrap); if (iflag/=0) return
    iflag = check_value(qval,tq,4,extrap); if (iflag/=0) return
    iflag = check_value(rval,tr,5,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag,extrap); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag,extrap); if (mflag /= 0) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,mflag,extrap); if (mflag /= 0) return

    iflag = 0

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
                    call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr),&
                                nx,kx,idx,xval,inbvx,work,iflag,temp1(j,k,q,r),&
                                extrap)
                    if (iflag/=0) return
                end do
            end do
        end do
    end do

    ! y -> z, q, r
    kcoly = lefty - ky + 1
    do r=1,kr
        do q=1,kq
            do k=1,kz
                call dbvalu(ty(kcoly:),temp1(:,k,q,r),ky,ky,idy,yval,inbvy,&
                            work,iflag,temp2(k,q,r),extrap)
                if (iflag/=0) return
            end do
        end do
    end do

    ! z -> q, r
    kcolz = leftz - kz + 1
    do r=1,kr
        do q=1,kq
            call dbvalu(tz(kcolz:),temp2(:,q,r),kz,kz,idz,zval,inbvz,&
                        work,iflag,temp3(q,r),extrap)
            if (iflag/=0) return
        end do
    end do

    ! q -> r
    kcolq = leftq - kq + 1
    do r=1,kr
        call dbvalu(tq(kcolq:),temp3(:,r),kq,kq,idq,qval,inbvq,&
                    work,iflag,temp4(r),extrap)
        if (iflag/=0) return
    end do

    ! r
    kcolr = leftr - kr + 1
    call dbvalu(tr(kcolr:),temp4,kr,kr,idr,rval,inbvr,work,iflag,f,extrap)

    end subroutine db5val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
!  the six-dimensional gridded data:
!
!  $$ [x(i),y(j),z(k),q(l),r(m),s(n),\mathrm{fcn}(i,j,k,l,m,n)] $$
!
!  for:
!
!  $$ i=1,..,n_x, j=1,..,n_y, k=1,..,n_z, l=1,..,n_q, m=1,..,n_r, n=1,..,n_s $$
!
!  the interpolating function and its derivatives may subsequently be evaluated
!  by the function [[db6val]].
!
!  See [[db3ink]] header for more details.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                        fcn,&
                        kx,ky,kz,kq,kr,ks,&
                        iknot,&
                        tx,ty,tz,tq,tr,ts,&
                        bcoef,iflag)

    implicit none

    integer,intent(in)                                :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: nr    !! number of \(r\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: ns    !! number of \(s\) abcissae ( \( \ge 3 \) )
    integer,intent(in)                                :: kx    !! the order of spline pieces in \(x\)
                                                               !! ( \( 2 \le k_x < n_x \) )
                                                               !! (order = polynomial degree + 1)
    integer,intent(in)                                :: ky    !! the order of spline pieces in \(y\)
                                                               !! ( \( 2 \le k_y < n_y \) )
                                                               !! (order = polynomial degree + 1)
    integer,intent(in)                                :: kz    !! the order of spline pieces in \(z\)
                                                               !! ( \( 2 \le k_z < n_z \) )
                                                               !! (order = polynomial degree + 1)
    integer,intent(in)                                :: kq    !! the order of spline pieces in \(q\)
                                                               !! ( \( 2 \le k_q < n_q \) )
                                                               !! (order = polynomial degree + 1)
    integer,intent(in)                                :: kr    !! the order of spline pieces in \(r\)
                                                               !! ( \( 2 \le k_r < n_r \) )
                                                               !! (order = polynomial degree + 1)
    integer,intent(in)                                :: ks    !! the order of spline pieces in \(s\)
                                                               !! ( \( 2 \le k_s < n_s \) )
                                                               !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)                  :: x     !! `(nx)` array of \(x\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:),intent(in)                  :: y     !! `(ny)` array of \(y\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:),intent(in)                  :: z     !! `(nz)` array of \(z\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:),intent(in)                  :: q     !! `(nq)` array of \(q\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:),intent(in)                  :: r     !! `(nr)` array of \(r\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:),intent(in)                  :: s     !! `(ns)` array of \(s\) abcissae.
                                                               !! must be strictly increasing.
    real(wp),dimension(:,:,:,:,:,:),intent(in)        :: fcn   !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to interpolate.
                                                               !! `fcn(i,j,k,q,r,s)` should contain the function value
                                                               !! at the point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer,intent(in)                                :: iknot !! knot sequence flag:
                                                               !!
                                                               !! * 0 = knot sequence chosen by [[db6ink]].
                                                               !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)               :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the spline interpolant.
                                                               !!
                                                               !! * f `iknot=0` these are chosen by [[db6ink]].
                                                               !! * f `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tq    !! The `(nq+kq)` knots in the \(q\) direction for the spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tr    !! The `(nr+kr)` knots in the \(r\) direction for the spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: ts    !! The `(ns+ks)` knots in the \(s\) direction for the spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:,:,:),intent(out)       :: bcoef !! `(nx,ny,nz,nq,nr,ns)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                               :: iflag !! *  0 = successful execution.
                                                               !! *  2 = `iknot` out of range.
                                                               !! *  3 = `nx` out of range.
                                                               !! *  4 = `kx` out of range.
                                                               !! *  5 = `x` not strictly increasing.
                                                               !! *  6 = `tx` not non-decreasing.
                                                               !! *  7 = `ny` out of range.
                                                               !! *  8 = `ky` out of range.
                                                               !! *  9 = `y` not strictly increasing.
                                                               !! * 10 = `ty` not non-decreasing.
                                                               !! * 11 = `nz` out of range.
                                                               !! * 12 = `kz` out of range.
                                                               !! * 13 = `z` not strictly increasing.
                                                               !! * 14 = `tz` not non-decreasing.
                                                               !! * 15 = `nq` out of range.
                                                               !! * 16 = `kq` out of range.
                                                               !! * 17 = `q` not strictly increasing.
                                                               !! * 18 = `tq` not non-decreasing.
                                                               !! * 19 = `nr` out of range.
                                                               !! * 20 = `kr` out of range.
                                                               !! * 21 = `r` not strictly increasing.
                                                               !! * 22 = `tr` not non-decreasing.
                                                               !! * 23 = `ns` out of range.
                                                               !! * 24 = `ks` out of range.
                                                               !! * 25 = `s` not strictly increasing.
                                                               !! * 26 = `ts` not non-decreasing.
                                                               !! * 700 = `size(x) ` \( \ne \) `size(fcn,1)`
                                                               !! * 701 = `size(y) ` \( \ne \) `size(fcn,2)`
                                                               !! * 702 = `size(z) ` \( \ne \) `size(fcn,3)`
                                                               !! * 703 = `size(q) ` \( \ne \) `size(fcn,4)`
                                                               !! * 704 = `size(r) ` \( \ne \) `size(fcn,5)`
                                                               !! * 705 = `size(s) ` \( \ne \) `size(fcn,6)`
                                                               !! * 706 = `size(x) ` \( \ne \) `nx`
                                                               !! * 707 = `size(y) ` \( \ne \) `ny`
                                                               !! * 708 = `size(z) ` \( \ne \) `nz`
                                                               !! * 709 = `size(q) ` \( \ne \) `nq`
                                                               !! * 710 = `size(r) ` \( \ne \) `nr`
                                                               !! * 711 = `size(s) ` \( \ne \) `ns`
                                                               !! * 712 = `size(tx)` \( \ne \) `nx+kx`
                                                               !! * 713 = `size(ty)` \( \ne \) `ny+ky`
                                                               !! * 714 = `size(tz)` \( \ne \) `nz+kz`
                                                               !! * 715 = `size(tq)` \( \ne \) `nq+kq`
                                                               !! * 716 = `size(tr)` \( \ne \) `nr+kr`
                                                               !! * 717 = `size(ts)` \( \ne \) `ns+ks`
                                                               !! * 800 = `size(x) ` \( \ne \) `size(bcoef,1)`
                                                               !! * 801 = `size(y) ` \( \ne \) `size(bcoef,2)`
                                                               !! * 802 = `size(z) ` \( \ne \) `size(bcoef,3)`
                                                               !! * 803 = `size(q) ` \( \ne \) `size(bcoef,4)`
                                                               !! * 804 = `size(r) ` \( \ne \) `size(bcoef,5)`
                                                               !! * 805 = `size(s) ` \( \ne \) `size(bcoef,6)`

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
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,nq=nq,nr=nr,ns=ns,&
                        kx=kx,ky=ky,kz=kz,kq=kq,kr=kr,ks=ks,&
                        x=x,y=y,z=z,q=q,r=r,s=s,&
                        tx=tx,ty=ty,tz=tz,tq=tq,tr=tr,ts=ts,&
                        f6=fcn,&
                        bcoef6=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        ! choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
            call dbknot(s,ns,ks,ts)
        end if

        ! construct b-spline coefficients

                      call dbtpcf(x,nx,fcn,  nx,ny*nz*nq*nr*ns,tx,kx,temp, work,iflag)
        if (iflag==0) call dbtpcf(y,ny,temp, ny,nx*nz*nq*nr*ns,ty,ky,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(z,nz,bcoef,nz,nx*ny*nq*nr*ns,tz,kz,temp, work,iflag)
        if (iflag==0) call dbtpcf(q,nq,temp, nq,nx*ny*nz*nr*ns,tq,kq,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(r,nr,bcoef,nr,nx*ny*nz*nq*ns,tr,kr,temp, work,iflag)
        if (iflag==0) call dbtpcf(s,ns,temp, ns,nx*ny*nz*nq*nr,ts,ks,bcoef,work,iflag)

     end if

    end subroutine db6ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db6ink]] or one of its
!  derivatives at the point (`xval`,`yval`,`zval`,`qval`,`rval`,`sval`).
!
!  To evaluate the
!  interpolant itself, set `idx=idy=idz=idq=idr=ids=0`, to evaluate the first
!  partial with respect to `x`, set `idx=1,idy=idz=idq=idr=ids=0`, and so on.
!
!  See [[db3val]] header for more information.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine db6val(xval,yval,zval,qval,rval,sval,&
                                idx,idy,idz,idq,idr,ids,&
                                tx,ty,tz,tq,tr,ts,&
                                nx,ny,nz,nq,nr,ns,&
                                kx,ky,kz,kq,kr,ks,&
                                bcoef,f,iflag,&
                                inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,&
                                iloy,iloz,iloq,ilor,ilos,extrap)

    implicit none

    integer,intent(in)                               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: idr      !! \(r\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: ids      !! \(s\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                               :: nx       !! the number of interpolation points in \(x\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: ny       !! the number of interpolation points in \(y\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: nz       !! the number of interpolation points in \(z\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: nq       !! the number of interpolation points in \(q\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: nr       !! the number of interpolation points in \(r\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: ns       !! the number of interpolation points in \(s\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: kx       !! order of polynomial pieces in \(x\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: ky       !! order of polynomial pieces in \(y\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: kz       !! order of polynomial pieces in \(z\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: kq       !! order of polynomial pieces in \(q\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: kr       !! order of polynomial pieces in \(r\).
                                                                 !! (same as in last call to [[db6ink]])
    integer,intent(in)                               :: ks       !! order of polynomial pieces in \(s\).
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),intent(in)                              :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                              :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)                              :: zval     !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)                              :: qval     !! \(q\) coordinate of evaluation point.
    real(wp),intent(in)                              :: rval     !! \(r\) coordinate of evaluation point.
    real(wp),intent(in)                              :: sval     !! \(s\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)             :: tx       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(x\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(ny+ky),intent(in)             :: ty       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(y\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(nz+kz),intent(in)             :: tz       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(z\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(nq+kq),intent(in)             :: tq       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(q\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(nr+kr),intent(in)             :: tr       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(r\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(ns+ks),intent(in)             :: ts       !! sequence of knots defining the piecewise polynomial
                                                                 !! in the \(s\) direction.
                                                                 !! (same as in last call to [[db6ink]])
    real(wp),dimension(nx,ny,nz,nq,nr,ns),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db6ink]].
    real(wp),intent(out)                             :: f        !! interpolated value
    integer,intent(out)                              :: iflag    !! status flag:
                                                                 !!
                                                                 !! * \( = 0 \)   : no errors
                                                                 !! * \( \ne 0 \) : error
    integer,intent(inout)                            :: inbvx    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: inbvy    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: inbvz    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: inbvq    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: inbvr    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: inbvs    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: iloy     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: iloz     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: iloq     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: ilor     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer,intent(inout)                            :: ilos     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    logical,intent(in),optional                      :: extrap   !! if extrapolation is allowed
                                                                 !! (if not present, default is False)

    real(wp),dimension(ky,kz,kq,kr,ks)            :: temp1
    real(wp),dimension(kz,kq,kr,ks)               :: temp2
    real(wp),dimension(kq,kr,ks)                  :: temp3
    real(wp),dimension(kr,ks)                     :: temp4
    real(wp),dimension(ks)                        :: temp5
    real(wp),dimension(3*max(kx,ky,kz,kq,kr,ks))  :: work

    integer :: lefty,leftz,leftq,leftr,lefts,&
               mflag,&
               kcoly,kcolz,kcolq,kcolr,kcols,&
               j,k,q,r,s

    f = 0.0_wp

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return
    iflag = check_value(zval,tz,3,extrap); if (iflag/=0) return
    iflag = check_value(qval,tq,4,extrap); if (iflag/=0) return
    iflag = check_value(rval,tr,5,extrap); if (iflag/=0) return
    iflag = check_value(sval,ts,6,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag,extrap); if (mflag /= 0) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,mflag,extrap); if (mflag /= 0) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,mflag,extrap); if (mflag /= 0) return
    call dintrv(ts,ns+ks,sval,ilos,lefts,mflag,extrap); if (mflag /= 0) return

    iflag = 0

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
                        call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr,kcols),&
                                             nx,kx,idx,xval,inbvx,work,iflag,&
                                             temp1(j,k,q,r,s),extrap)
                        if (iflag/=0) return
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
                    call dbvalu(ty(kcoly:),temp1(:,k,q,r,s),&
                                ky,ky,idy,yval,inbvy,work,iflag,&
                                temp2(k,q,r,s),extrap)
                    if (iflag/=0) return
                end do
            end do
        end do
    end do

    ! z -> q, r, s
    kcolz = leftz - kz + 1
    do s=1,ks
        do r=1,kr
            do q=1,kq
                call dbvalu(tz(kcolz:),temp2(:,q,r,s),&
                            kz,kz,idz,zval,inbvz,work,iflag,&
                            temp3(q,r,s),extrap)
                if (iflag/=0) return
            end do
        end do
    end do

    ! q -> r, s
    kcolq = leftq - kq + 1
    do s=1,ks
        do r=1,kr
            call dbvalu(tq(kcolq:),temp3(:,r,s),&
                        kq,kq,idq,qval,inbvq,work,iflag,&
                        temp4(r,s),extrap)
            if (iflag/=0) return
        end do
    end do

    ! r -> s
    kcolr = leftr - kr + 1
    do s=1,ks
        call dbvalu(tr(kcolr:),temp4(:,s),&
                    kr,kr,idr,rval,inbvr,work,iflag,&
                    temp5(s),extrap)
        if (iflag/=0) return
    end do

    ! s
    kcols = lefts - ks + 1
    call dbvalu(ts(kcols:),temp5,ks,ks,ids,sval,inbvs,work,iflag,f,extrap)

    end subroutine db6val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Checks if the value is withing the range of the knot vectors.
!  This is called by the various `db*val` routines.

    pure function check_value(x,t,i,extrap) result(iflag)

    implicit none

    integer :: iflag  !! returns 0 if value is OK, otherwise returns `600+i`
    real(wp),intent(in) :: x !! the value to check
    integer,intent(in) :: i !! 1=x, 2=y, 3=z, 4=q, 5=r, 6=s
    real(wp),dimension(:),intent(in) :: t  !! the knot vector
    logical,intent(in),optional :: extrap  !! if extrapolation is allowed
                                           !! (if not present, default is False)

    logical :: allow_extrapolation  !! if extrapolation is allowed

    if (present(extrap)) then
        allow_extrapolation = extrap
    else
        allow_extrapolation = .false.
    end if

    if (allow_extrapolation) then
        ! in this case all values are OK
        iflag = 0
    else
        if (x<t(1) .or. x>t(size(t))) then
            iflag = 600 + i  ! value out of bounds (601, 602, etc.)
        else
            iflag = 0
        end if
    end if

    end function check_value
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the `db*ink` routines.
!  Prints warning message if there is an error,
!  and also sets iflag and status_ok.
!
!  Supports up to 6D: `x`,`y`,`z`,`q`,`r`,`s`
!
!### Notes
!
!  The code is new, but the logic is based on the original
!  logic in the CMLIB routines `db2ink` and `db3ink`.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine check_inputs(routine,&
                            iknot,&
                            iflag,&
                            nx,ny,nz,nq,nr,ns,&
                            kx,ky,kz,kq,kr,ks,&
                            x,y,z,q,r,s,&
                            tx,ty,tz,tq,tr,ts,&
                            f1,f2,f3,f4,f5,f6,&
                            bcoef1,bcoef2,bcoef3,bcoef4,bcoef5,bcoef6,&
                            status_ok)

    implicit none

    character(len=*),intent(in)                         :: routine
    integer,intent(in)                                  :: iknot !! = 0 if the `INK` routine is computing the knots.
    integer,intent(out)                                 :: iflag
    integer,intent(in),optional                         :: nx,ny,nz,nq,nr,ns
    integer,intent(in),optional                         :: kx,ky,kz,kq,kr,ks
    real(wp),dimension(:),intent(in),optional           :: x,y,z,q,r,s
    real(wp),dimension(:),intent(in),optional           :: tx,ty,tz,tq,tr,ts
    real(wp),dimension(:),intent(in),optional           :: f1,bcoef1
    real(wp),dimension(:,:),intent(in),optional         :: f2,bcoef2
    real(wp),dimension(:,:,:),intent(in),optional       :: f3,bcoef3
    real(wp),dimension(:,:,:,:),intent(in),optional     :: f4,bcoef4
    real(wp),dimension(:,:,:,:,:),intent(in),optional   :: f5,bcoef5
    real(wp),dimension(:,:,:,:,:,:),intent(in),optional :: f6,bcoef6
    logical,intent(out)                                 :: status_ok

    logical :: error

    status_ok = .false.

    if ((iknot < 0) .or. (iknot > 1)) then

        !write(error_unit,'(A,1X,I5)') &
        !    trim(routine)//' - iknot is out of range: ',iflag
        iflag = 2

    else

        call check('x',nx,kx,x,tx,[3,  4, 5, 6,706,712],iflag,error); if (error) return
        call check('y',ny,ky,y,ty,[7,  8, 9,10,707,713],iflag,error); if (error) return
        call check('z',nz,kz,z,tz,[11,12,13,14,708,714],iflag,error); if (error) return
        call check('q',nq,kq,q,tq,[15,16,17,18,709,715],iflag,error); if (error) return
        call check('r',nr,kr,r,tr,[19,20,21,22,710,716],iflag,error); if (error) return
        call check('s',ns,ks,s,ts,[23,24,25,26,711,717],iflag,error); if (error) return

        if (present(x) .and. present(f1) .and. present(bcoef1)) then
            if (size(x)/=size(f1,1))     then; iflag = 700; return; end if
            if (size(x)/=size(bcoef1,1)) then; iflag = 800; return; end if
        end if
        if (present(x) .and. present(y) .and. present(f2) .and. present(bcoef2)) then
            if (size(x)/=size(f2,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f2,2))     then; iflag = 701; return; end if
            if (size(x)/=size(bcoef2,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef2,2)) then; iflag = 801; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(f3) .and. &
            present(bcoef3)) then
            if (size(x)/=size(f3,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f3,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f3,3))     then; iflag = 702; return; end if
            if (size(x)/=size(bcoef3,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef3,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef3,3)) then; iflag = 802; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(f4) .and. present(bcoef4)) then
            if (size(x)/=size(f4,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f4,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f4,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f4,4))     then; iflag = 703; return; end if
            if (size(x)/=size(bcoef4,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef4,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef4,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef4,4)) then; iflag = 803; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(f5) .and. present(bcoef5)) then
            if (size(x)/=size(f5,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f5,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f5,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f5,4))     then; iflag = 703; return; end if
            if (size(r)/=size(f5,5))     then; iflag = 704; return; end if
            if (size(x)/=size(bcoef5,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef5,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef5,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef5,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef5,5)) then; iflag = 804; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(s) .and. present(f6) .and. present(bcoef6)) then
            if (size(x)/=size(f6,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f6,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f6,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f6,4))     then; iflag = 703; return; end if
            if (size(r)/=size(f6,5))     then; iflag = 704; return; end if
            if (size(s)/=size(f6,6))     then; iflag = 705; return; end if
            if (size(x)/=size(bcoef6,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef6,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef6,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef6,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef6,5)) then; iflag = 804; return; end if
            if (size(s)/=size(bcoef6,6)) then; iflag = 805; return; end if
        end if

        status_ok = .true.
        iflag = 0

    end if

    contains

        pure subroutine check(s,n,k,x,t,ierrs,iflag,error)  !! check `t`,`x`,`n`,`k` for validity

        implicit none

        character(len=1),intent(in)               :: s     !! coordinate string: 'x','y','z','q','r','s'
        integer,intent(in)              ,optional :: n     !! size of `x`
        integer,intent(in)              ,optional :: k     !! order
        real(wp),dimension(:),intent(in),optional :: x     !! abcissae vector
        real(wp),dimension(:),intent(in),optional :: t     !! knot vector `size(n+k)`
        integer,dimension(:),intent(in)           :: ierrs !! int error codes for `n`,`k`,`x`,`t`,
                                                           !! `size(x)`,`size(t)` checks
        integer,intent(out)                       :: iflag !! status return code
        logical,intent(out)                       :: error !! true if there was an error

        if (present(n) .and. present(k) .and. present(x) .and. present(t)) then
            call check_n('n'//s,n,x,[ierrs(1),ierrs(5)],iflag,error); if (error) return
            call check_k('k'//s,k,n,ierrs(2),iflag,error); if (error) return
            call check_x(s,n,x,ierrs(3),iflag,error); if (error) return
            if (iknot /= 0) then
                call check_t('t'//s,n,k,t,[ierrs(4),ierrs(6)],iflag,error); if (error) return
            end if
        end if

        end subroutine check

        pure subroutine check_n(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)     :: s
        integer,intent(in)              :: n
        real(wp),dimension(:),intent(in):: x     !! abcissae vector
        integer,dimension(2),intent(in) :: ierr  !! [n<3 check, size(x)==n check]
        integer,intent(out)             :: iflag !! status return code
        logical,intent(out)             :: error

        if (n < 3) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',n
            iflag = ierr(1)
            error = .true.
        else
            if (size(x)/=n) then
                !write(error_unit,'(A,1X,I5)') &
                !    trim(routine)//' - '//trim(s)//' is not abscissa vector size'
                iflag = ierr(2)
                error = .true.
            else
                error = .false.
            end if
        end if

        end subroutine check_n

        pure subroutine check_k(s,k,n,ierr,iflag,error)

        implicit none

        character(len=*),intent(in) :: s
        integer,intent(in)          :: k
        integer,intent(in)          :: n
        integer,intent(in)          :: ierr
        integer,intent(out)         :: iflag !! status return code
        logical,intent(out)         :: error

        if ((k < 2) .or. (k >= n)) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',k
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if

        end subroutine check_k

        pure subroutine check_x(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        real(wp),dimension(:),intent(in)  :: x
        integer,intent(in)                :: ierr
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.
        do i=2,n
            if (x(i) <= x(i-1)) then
                iflag = ierr
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be strictly increasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_x

        pure subroutine check_t(s,n,k,t,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        integer,intent(in)                :: k
        real(wp),dimension(:),intent(in)  :: t
        integer,dimension(2),intent(in)   :: ierr  !! [non-decreasing check, size check]
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.

        if (size(t)/=(n+k)) then
            !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
            !            ' array is not the correct size'
            iflag = ierr(2)
            return
        end if

        do i=2,n + k
            if (t(i) < t(i-1))  then
                iflag = ierr(1)
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be non-decreasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_t

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbknot chooses a knot sequence for interpolation of order k at the
!  data points x(i), i=1,..,n.  the n+k knots are placed in the array
!  t.  k knots are placed at each endpoint and not-a-knot end
!  conditions are used.  the remaining knots are placed at data points
!  if n is even and between data points if n is odd.  the rightmost
!  knot is shifted slightly to the right to insure proper interpolation
!  at x(n) (see page 350 of the reference).
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbknot(x,n,k,t)

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
!>
!  dbtpcf computes b-spline interpolation coefficients for nf sets
!  of data stored in the columns of the array fcn. the b-spline
!  coefficients are stored in the rows of bcoef however.
!  each interpolation is based on the n abcissa stored in the
!  array x, and the n+k knots stored in the array t. the order
!  of each interpolation is k.
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)

    integer,intent(in)                    :: n
    integer,intent(in)                    :: nf
    integer,intent(in)                    :: ldf
    integer,intent(in)                    :: k
    real(wp),dimension(n),intent(in)      :: x
    real(wp),dimension(ldf,nf),intent(in) :: fcn
    real(wp),dimension(*),intent(in)      :: t
    real(wp),dimension(nf,n),intent(out)  :: bcoef
    real(wp),dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer,intent(out)                   :: iflag  !!   0: no errors
                                                    !! 301: n should be >0

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

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0) then
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

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301
    end if

    end subroutine dbtpcf
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbintk produces the b-spline coefficients, bcoef, of the
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
!# Error conditions
!
!  * improper input
!  * singular system of equations
!
!### History
!  * splint written by carl de boor [5]
!  * dbintk author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations. (jec)
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbintk(x,y,t,n,k,bcoef,q,work,iflag)

    implicit none

    integer,intent(in)                :: n      !! number of data points, n >= k
    real(wp),dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(wp),dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real(wp),dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer,intent(in)                :: k      !! order of the spline, k >= 1
    real(wp),dimension(n),intent(out) :: bcoef  !! a vector of length n containing the b-spline coefficients
    real(wp),dimension(*),intent(out) :: q      !! a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real(wp),dimension(*),intent(out) :: work   !! work vector of length 2*k
    integer,intent(out)               :: iflag  !! *   0: no errors.
                                                !! * 100: k does not satisfy k>=1.
                                                !! * 101: n does not satisfy n>=k.
                                                !! * 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! * 103: some abscissa was not in the support of the.
                                                !! corresponding basis function and the system is singular.
                                                !! * 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(wp) :: xi
    logical :: found

    if (k<1) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101
        return
    end if

    jj = n - 1
    if (jj/=0) then
        do i=1,jj
            if (x(i)>=x(i+1)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102
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
        ! find left in the closed interval (i,i+k-1) such that
        !         t(left) <= x(i) < t(left+1)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1))
            if (found) exit
            left = left + 1
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1
            if (xi>t(left+1)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0) return

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
        ! solve  a*bcoef = y  by backsubstitution
        do i=1,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104
    end if

    end subroutine dbintk
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns in w the LU-factorization (without pivoting) of the banded
!  matrix a of order nrow with (nbandl + 1 + nbandu) bands or diagonals
!  in the work array w .
!
!  gauss elimination without pivoting is used. the routine is
!  intended for use with matrices a which do not require row inter-
!  changes during factorization, especially for the totally
!  positive matrices which occur in spline calculations.
!  the routine should not be used for an arbitrary banded matrix.
!
!### Work array
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
!### History
!  * banfac written by carl de boor [5]
!  * dbnfac from CMLIB [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)

    integer,intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer,intent(in) :: nrow    !! matrix order
    integer,intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer,intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer,intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real(wp),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real(wp) :: factor, pivot

    iflag = 1
    middle = nbandu + 1   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1

    if (nrowm1 < 0) then
        iflag = 2
        return
    else if (nrowm1 == 0) then
        if (w(middle,nrow)==0.0_wp) iflag = 2
        return
    end if

    if (nbandl<=0) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1,nrowm1
            if (w(middle,i)==0.0_wp) then
                iflag = 2
                return
            end if
        end do
        if (w(middle,nrow)==0.0_wp) iflag = 2
        return
    end if

    if (nbandu<=0) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0_wp) then
                iflag = 2
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0_wp) then
            iflag = 2
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0_wp) iflag = 2

    end subroutine dbnfac
!*****************************************************************************************

!*****************************************************************************************
!>
!  Companion routine to [[dbnfac]]. it returns the solution x of the
!  linear system a*x = b in place of b, given the lu-factorization
!  for a in the work array w from dbnfac.
!
!  (with \( a = l*u \), as stored in w), the unit lower triangular system
!  \( l(u*x) = b \) is solved for \( y = u*x \), and y stored in b. then the
!  upper triangular system \(u*x = y \) is solved for x. the calculations
!  are so arranged that the innermost loops stay within columns.
!
!### History
!  * banslv written by carl de boor [5]
!  * dbnslv from SLATEC library [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)

    integer,intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    real(wp),dimension(nroww,nrow),intent(in) :: w    !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    real(wp),dimension(nrow),intent(inout) :: b  !! * **in**: right side of the system to be solved
                                                 !! * **out**: the solution x, of order nrow

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
!>
!  Calculates the value of all (possibly) nonzero basis
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
!### Error Conditions
!
!  * improper input
!
!### History
!  * bsplvn written by carl de boor [5]
!  * dbspvn author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)

    implicit none

    real(wp),dimension(*),intent(in)  :: t        !! knot vector of length n+k, where
                                                  !! n = number of b-spline basis functions
                                                  !! n = sum of knot multiplicities-k
                                                  !! dimension t(ileft+jhigh)
    integer,intent(in)                :: jhigh    !! order of b-spline, 1 <= jhigh <= k
    integer,intent(in)                :: k        !! highest possible order
    integer,intent(in)                :: index    !! index = 1 gives basis functions of order jhigh
                                                  !!       = 2 denotes previous entry with work, iwork
                                                  !!         values saved for subsequent calls to
                                                  !!         dbspvn.
    real(wp),intent(in)               :: x        !! argument of basis functions, t(k) <= x <= t(n+1)
    integer,intent(in)                :: ileft    !! largest integer such that t(ileft) <= x < t(ileft+1)
    real(wp),dimension(k),intent(out) :: vnikx    !! vector of length k for spline values.
    real(wp),dimension(*),intent(out) :: work     !! a work vector of length 2*k
    integer,intent(out)               :: iwork    !! a work parameter.  both work and iwork contain
                                                  !! information necessary to continue for index = 2.
                                                  !! when index = 1 exclusively, these are scratch
                                                  !! variables and can be used for other purposes.
    integer,intent(out)               :: iflag    !!   0: no errors
                                                  !! 201: k does not satisfy k>=1
                                                  !! 202: jhigh does not satisfy 1<=jhigh<=k
                                                  !! 203: index is not 1 or 2
                                                  !! 204: x does not satisfy t(ileft)<=x<=t(ileft+1)

    integer :: imjp1, ipj, jp1, jp1ml, l
    real(wp) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201
        return
    end if
    if (jhigh>k .or. jhigh<1) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202
        return
    end if
    if (index<1 .or. index>2) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204
        return
    end if

    iflag = 0

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
!>
!  Evaluates the b-representation (`t`,`a`,`n`,`k`) of a b-spline
!  at `x` for the function value on `ideriv=0` or any of its
!  derivatives on `ideriv=1,2,...,k-1`.  right limiting values
!  (right derivatives) are returned except at the right end
!  point `x=t(n+1)` where left limiting values are computed.  the
!  spline is defined on `t(k)` \( \le \) `x` \( \le \) `t(n+1)`.
!  dbvalu returns a fatal error message when `x` is outside of this
!  interval.
!
!  To compute left derivatives or left limiting values at a
!  knot `t(i)`, replace `n` by `i-1` and set `x=t(i), i=k+1,n+1`.
!
!### Error Conditions
!
!  * improper input
!
!### History
!  * bvalue written by carl de boor [5]
!  * dbvalu author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val,extrap)

    implicit none

    real(wp),intent(out)             :: val     !! the interpolated value
    integer,intent(in)               :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-`k`)
    real(wp),dimension(:),intent(in) :: t       !! knot vector of length `n+k`
    real(wp),dimension(n),intent(in) :: a       !! b-spline coefficient vector of length `n`
    integer,intent(in)               :: k       !! order of the b-spline, `k >= 1`
    integer,intent(in)               :: ideriv  !! order of the derivative, `0 <= ideriv <= k-1`.
                                                !! `ideriv = 0` returns the b-spline value
    real(wp),intent(in)              :: x       !! argument, `t(k) <= x <= t(n+1)`
    integer,intent(inout)            :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time [[dbvalu]] is called.
                                                !! `inbv` contains information for efficient processing
                                                !! after the initial call and `inbv` must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct `inbv` parameters.
    real(wp),dimension(:),intent(inout) :: work !! work vector of length at least `3*k`
    integer,intent(out)              :: iflag   !! status flag:
                                                !!
                                                !! * 0: no errors
                                                !! * 401: `k` does not satisfy `k` \( \ge \) 1
                                                !! * 402: `n` does not satisfy `n` \( \ge \) `k`
                                                !! * 403: `ideriv` does not satisfy 0 \( \le \) `ideriv` \(<\) `k`
                                                !! * 404: `x` is not greater than or equal to `t(k)`
                                                !! * 405: `x` is not less than or equal to `t(n+1)`
                                                !! * 406: a left limiting value cannot be obtained at `t(k)`
    logical,intent(in),optional :: extrap   !! if extrapolation is allowed
                                            !! (if not present, default is False)

    integer :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(wp) :: fkmj
    real(wp) :: xt
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! make a temp copy of x (for computing the
    ! interval) in case extrapolation is allowed
    if (extrapolation_allowed) then
        if (x<t(1)) then
            xt = t(1)
        else if(x>t(n+k)) then
            xt = t(n+k)
        else
            xt = x
        end if
    else
        xt = x
    end if

    val = 0.0_wp

    if (k<1) then
        iflag = 401  ! dbvalu - k does not satisfy k>=1
        return
    end if

    if (n<k) then
        iflag = 402  ! dbvalu - n does not satisfy n>=k
        return
    end if

    if (ideriv<0 .or. ideriv>=k) then
        iflag = 403  ! dbvalu - ideriv does not satisfy 0<=ideriv<k
        return
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1
    call dintrv(t, n+1, xt, inbv, i, mflag)
    if (xt<t(k)) then
        iflag = 404  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0) then

        if (xt>t(i)) then
            iflag = 405  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (i==k) then
                iflag = 406  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            i = i - 1
            if (xt/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

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

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

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

    iflag = 0
    val = work(1)

    end subroutine dbvalu
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the largest integer `ileft` in 1 \( \le \) `ileft` \( \le \) `lxt`
!  such that `xt(ileft)` \( \le \) `x` where `xt(*)` is a subdivision of
!  the `x` interval.
!  precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   mflag=0
!         if xt(lxt) <= x           then ileft=lxt, mflag=1
!```
!
!  that is, when multiplicities are present in the break point
!  to the left of `x`, the largest index is taken for `ileft`.
!
!### History
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 3/4/2017 : added extrapolation option.

    pure subroutine dintrv(xt,lxt,xx,ilo,ileft,mflag,extrap)

    implicit none

    integer,intent(in)                 :: lxt    !! length of the `xt` vector
    real(wp),dimension(lxt),intent(in) :: xt     !! a knot or break point vector of length `lxt`
    real(wp),intent(in)                :: xx     !! argument
    integer,intent(inout)              :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct `ilo` parameters.
    integer,intent(out)                :: ileft  !! largest integer satisfying `xt(ileft)` \( \le \) `x`
    integer,intent(out)                :: mflag  !! signals when `x` lies out of bounds
    logical,intent(in),optional        :: extrap !! if extrapolation is allowed
                                                 !! (if not present, default is False)

    integer :: ihi, istep, middle
    real(wp) :: x

    x = get_temp_x_for_extrap(xx,xt,extrap)

    ihi = ilo + 1
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = 1
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1
            ileft = 1
            return
        end if
        ilo = lxt - 1
        ihi = lxt
    end if

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = 1
                    ileft = lxt
                    return
                end if
                ihi = lxt
            else if ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            end if
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    return
                end if
            else if ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            end if
            exit
        end do

    end if

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2
        if ( middle==ilo ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! note. it is assumed that middle = ilo in case ihi = ilo+1
        if ( x<xt(middle) ) then
            ihi = middle
        else
            ilo = middle
        end if
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSQAD computes the integral on `(X1,X2)` of a `K`-th order
!  B-spline using the B-representation `(T,BCOEF,N,K)`.  Orders
!  `K` as high as 20 are permitted by applying a 2, 6, or 10
!  point Gauss formula on subintervals of `(X1,X2)` which are
!  formed by included (distinct) knots.
!
!  If orders `K` greater than 20 are needed, use [[dbfqad]] with
!  `F(X) = 1`.
!
!### Note
!  * The maximum number of significant digits obtainable in
!    DBSQAD is the smaller of 18 and the number of digits
!    carried in `real(wp)` arithmetic.
!
!### References
!  * D. E. Amos, "Quadrature subroutines for splines and
!    B-splines", Report SAND79-1825, Sandia Laboratories,
!    December 1979.
!
!### History
!  * Author: Amos, D. E., (SNLA)
!  * 800901  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, 9/6/2017 : refactored to modern Fortran.
!
!@note Extrapolation is not enabled for this routine.

    pure subroutine dbsqad(t,bcoef,n,k,x1,x2,bquad,work,iflag)

    implicit none

    real(wp),dimension(:),intent(in) :: t       !! knot array of length `n+k`
    real(wp),dimension(:),intent(in) :: bcoef   !! b-spline coefficient array of length `n`
    integer,intent(in)  :: n                    !! length of coefficient array
    integer,intent(in)  :: k                    !! order of b-spline, `1 <= k <= 20`
    real(wp),intent(in) :: x1                   !! end point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in) :: x2                   !! end point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(out) :: bquad               !! integral of the b-spline over (`x1`,`x2`)
    real(wp),dimension(:),intent(inout) :: work !! work vector of length `3*k`
    integer,intent(out) :: iflag   !! status flag:
                                   !!
                                   !! * 0: no errors
                                   !! * 901: `k` does not satisfy `1<=k<=20`
                                   !! * 902: `n` does not satisfy `n>=k`
                                   !! * 903: `x1` or `x2` or both do not satisfy `t(k)<=x<=t(n+1)`

    integer :: i,il1,il2,ilo,inbv,jf,left,m,mf,mflag,npk,np1
    real(wp) :: a,aa,b,bb,bma,bpa,c1,gx,q,ta,tb,y1,y2
    real(wp),dimension(5) :: s  !! sum

    real(wp),dimension(9),parameter :: gpts = [ &
        5.77350269189625764e-01_wp , 2.38619186083196909e-01_wp , 6.61209386466264514e-01_wp , &
        9.32469514203152028e-01_wp , 1.48874338981631211e-01_wp , 4.33395394129247191e-01_wp , &
        6.79409568299024406e-01_wp , 8.65063366688984511e-01_wp , 9.73906528517171720e-01_wp ]

    real(wp),dimension(9),parameter :: gwts = [ &
        1.00000000000000000e+00_wp , 4.67913934572691047e-01_wp , 3.60761573048138608e-01_wp , &
        1.71324492379170345e-01_wp , 2.95524224714752870e-01_wp , 2.69266719309996355e-01_wp , &
        2.19086362515982044e-01_wp , 1.49451349150580593e-01_wp , 6.66713443086881376e-02_wp ]

    iflag = 0
    bquad = 0.0_wp

    if ( k<1 .or. k>20 ) then

        iflag = 901 ! error return

    else if ( n<k ) then

        iflag = 902 ! error return

    else

        aa = min(x1,x2)
        bb = max(x1,x2)
        if ( aa>=t(k) ) then
            np1 = n + 1
            if ( bb<=t(np1) ) then
            if ( aa==bb ) return
            npk = n + k
            ! selection of 2, 6, or 10 point gauss formula
            jf = 0
            mf = 1
            if ( k>4 ) then
                jf = 1
                mf = 3
                if ( k>12 ) then
                    jf = 4
                    mf = 5
                end if
            end if
            do i = 1 , mf
                s(i) = 0.0_wp
            end do
            ilo = 1
            inbv = 1
            call dintrv(t,npk,aa,ilo,il1,mflag)
            call dintrv(t,npk,bb,ilo,il2,mflag)
            if ( il2>=np1 ) il2 = n
            do left = il1 , il2
                ta = t(left)
                tb = t(left+1)
                if ( ta/=tb ) then
                    a = max(aa,ta)
                    b = min(bb,tb)
                    bma = 0.5_wp*(b-a)
                    bpa = 0.5_wp*(b+a)
                    do m = 1 , mf
                        c1 = bma*gpts(jf+m)
                        gx = -c1 + bpa
                        call dbvalu(t,bcoef,n,k,0,gx,inbv,work,iflag,y2)
                        if (iflag/=0) return
                        gx = c1 + bpa
                        call dbvalu(t,bcoef,n,k,0,gx,inbv,work,iflag,y1)
                        if (iflag/=0) return
                        s(m) = s(m) + (y1+y2)*bma
                    end do
                end if
            end do
            q = 0.0_wp
            do m = 1 , mf
                q = q + gwts(jf+m)*s(m)
            end do
            if ( x1>x2 ) q = -q
                bquad = q
                return
            end if
        end if

        iflag = 903 ! error return

    end if

    end subroutine dbsqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBFQAD computes the integral on (X1,X2) of a product of a
!  function F and the ID-th derivative of a K-th order B-spline,
!  using the B-representation (T,BCOEF,N,K).  (X1,X2) must be a
!  subinterval of T(K) .LE. X .LE. T(N+1).  An integration rou-
!  tine, DBSGQ8 (a modification of GAUS8), integrates the product
!  on subintervals of (X1,X2) formed by included (distinct) knots
!
!  The maximum number of significant digits obtainable in
!  DBSQAD is the smaller of 18 and the number of digits
!  carried in real(wp) arithmetic.
!
!### Reference
!  * D. E. Amos, "Quadrature subroutines for splines and
!    B-splines", Report SAND79-1825, Sandia Laboratories,
!    December 1979.
!
!### History
!  * 800901  Amos, D. E., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, 9/6/2017 : refactored to modern Fortran. Some changes.

    subroutine dbfqad(f,t,bcoef,n,k,id,x1,x2,tol,quad,iflag,work)

    implicit none

    procedure(int_func)                 :: f      !! external function of one argument for the
                                                  !! integrand `bf(x)=f(x)*dbvalu(t,bcoef,n,k,id,x,inbv,work)`
    integer,intent(in)                  :: n      !! length of coefficient array
    integer,intent(in)                  :: k      !! order of b-spline, `k >= 1`
    real(wp),dimension(n+k),intent(in)  :: t      !! knot array
    real(wp),dimension(n),intent(in)    :: bcoef  !! coefficient array
    integer,intent(in)                  :: id     !! order of the spline derivative, `0 <= id <= k-1`
                                                  !! `id=0` gives the spline function
    real(wp),intent(in)                 :: x1     !! left point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: x2     !! right point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: tol    !! desired accuracy for the quadrature, suggest
                                                  !! `10*dtol < tol <= 0.1` where `dtol` is the maximum
                                                  !! of `1.0e-18` and real(wp) unit roundoff for
                                                  !! the machine
    real(wp),intent(out)                :: quad   !! integral of `bf(x)` on `(x1,x2)`
    real(wp),dimension(:),intent(inout) :: work   !! work vector of length `3*k`
    integer,intent(out)                 :: iflag  !! status flag:
                                                  !!
                                                  !! * 0: no errors
                                                  !! * 1001: `k` does not satisfy `k>=1`
                                                  !! * 1002: `n` does not satisfy `n>=k`
                                                  !! * 1003: `d` does not satisfy `0<=id<k`
                                                  !! * 1004: `x1` or `x2` or both do not
                                                  !!   satisfy `t(k)<=x<=t(n+1)`
                                                  !! * 1005: `tol` is less than machine
                                                  !!   epsilon or greater than 0.1

    integer :: inbv,ilo,il1,il2,left,mflag,npk,np1
    real(wp) :: a,aa,ans,b,bb,q,ta,tb,err

    real(wp),parameter :: eps = epsilon(1.0_wp) !! replaces d1mach(4) in original code.

    iflag = 0
    quad = 0.0_wp
    err = tol
    if ( k<1 ) then
        iflag = 1001     ! error
    elseif ( n<k ) then
        iflag = 1002     ! error
    elseif ( id<0 .or. id>=k ) then
        iflag = 1003     ! error
    else
        if ( tol>=eps .and. tol<=0.1_wp ) then
            aa = min(x1,x2)
            bb = max(x1,x2)
            if ( aa>=t(k) ) then
                np1 = n + 1
                if ( bb<=t(np1) ) then
                    if ( aa==bb ) return
                    npk = n + k
                    ilo = 1
                    call dintrv(t,npk,aa,ilo,il1,mflag)
                    call dintrv(t,npk,bb,ilo,il2,mflag)
                    if ( il2>=np1 ) il2 = n
                    inbv = 1
                    q = 0.0_wp
                    do left = il1 , il2
                        ta = t(left)
                        tb = t(left+1)
                        if ( ta/=tb ) then
                            a = max(aa,ta)
                            b = min(bb,tb)
                            call dbsgq8(f,t,bcoef,n,k,id,a,b,inbv,err,ans,iflag,work)
                            if ( iflag/=0 .and. iflag/=1101 ) return
                            q = q + ans
                        end if
                    end do
                    if ( x1>x2 ) q = -q
                    quad = q
                end if
            else
                iflag = 1004  ! error
            end if
        else
            iflag = 1005  ! error
        end if
    end if

    end subroutine dbfqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSGQ8, a modification of `GAUS8`, integrates the
!  product of FUN(X) by the ID-th derivative of a spline
!  [[DBVALU]]() between limits A and B.
!
!### See also
!  * [[dbfqad]]
!
!### History
!  * 800901  Jones, R. E., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890911  Removed unnecessary intrinsics.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 900328  Added TYPE section.  (WRB)
!  * 910408  Updated the AUTHOR section.  (WRB)

    subroutine dbsgq8(fun,xt,bc,n,kk,id,a,b,inbv,err,ans,iflag,work)

    implicit none

    procedure(int_func)                 :: fun     !! name of external function of one argument which multiplies dbvalu.
    integer,intent(in)                  :: n       !! number of b-coefficients for dbvalu
    integer,intent(in)                  :: kk      !! order of the spline, kk>=1
    real(wp),dimension(:),intent(in)    :: xt      !! knot array for dbvalu
    real(wp),dimension(n),intent(in)    :: bc      !! b-coefficient array for dbvalu
    integer,intent(in)                  :: id      !! Order of the spline derivative, 0<=ID<=KK-1
    real(wp),intent(in)                 :: a       !! lower limit of integral
    real(wp),intent(in)                 :: b       !! upper limit of integral (may be less than a)
    integer,intent(inout)               :: inbv    !! initialization parameter for [[dbvalu]]
    real(wp),intent(inout)              :: err     !! in: is a requested pseudorelative error tolerance.  normally
                                                   !! pick a value of abs(err)<1e-3.  ans will normally
                                                   !! have no more error than abs(err) times the integral of
                                                   !! the absolute value of fun(x)*dbvalu().
                                                   !!
                                                   !! out: will be an estimate of the absolute error in ans if the
                                                   !! input value of err was negative.  (err is unchanged if
                                                   !! the input value of err was nonnegative.)  the estimated
                                                   !! error is solely for information to the user and should
                                                   !! not be used as a correction to the computed integral.
    real(wp),intent(out)                :: ans     !! computed value of integral
    integer,intent(out)                 :: iflag   !! a status code:
                                                   !!
                                                   !! * 0: ans most likely meets requested
                                                   !!   error tolerance, or a=b.
                                                   !! * 1101: a and b are too nearly equal
                                                   !!   to allow normal integration.
                                                   !!   ans is set to zero.
                                                   !! * 1102: ans probably does not meet
                                                   !!   requested error tolerance.
    real(wp),dimension(:),intent(inout) :: work    !! work vector of length `3*k` for [[dbvalu]]

    integer :: k,l,lmn,lmx,mxl,nbits,nib,nlmx
    real(wp) :: ae,anib,area,c,ce,ee,ef,eps,est,gl,glr,tol,vr,x,h

    real(wp),dimension(60) :: aa
    real(wp),dimension(60) :: hh
    integer,dimension(60)  :: lr
    real(wp),dimension(60) :: vl
    real(wp),dimension(60) :: gr

    integer,parameter  :: i1mach14 = digits(1.0_wp)            !! i1mach(14)
    real(wp),parameter :: d1mach5  = log10(real(radix(x),wp))  !! d1mach(5)
    real(wp),parameter :: ln2      = log(2.0_wp)               !! 0.69314718d0
    real(wp),parameter :: sq2      = sqrt(2.0_wp)
    integer,parameter  :: nlmn     = 1
    integer,parameter  :: kmx      = 5000
    integer,parameter  :: kml      = 6
    real(wp),parameter :: x1       = 1.83434642495649805e-01_wp
    real(wp),parameter :: x2       = 5.25532409916328986e-01_wp
    real(wp),parameter :: x3       = 7.96666477413626740e-01_wp
    real(wp),parameter :: x4       = 9.60289856497536232e-01_wp
    real(wp),parameter :: w1       = 3.62683783378361983e-01_wp
    real(wp),parameter :: w2       = 3.13706645877887287e-01_wp
    real(wp),parameter :: w3       = 2.22381034453374471e-01_wp
    real(wp),parameter :: w4       = 1.01228536290376259e-01_wp

    ! initialize
    inbv  = 1
    iflag = 0
    k     = i1mach14
    anib  = d1mach5*k/0.30102000_wp
    nbits = int(anib)
    nlmx  = min((nbits*5)/8,60)
    ans   = 0.0_wp
    ce    = 0.0_wp

    if ( a==b ) then
        if ( err<0.0_wp ) err = ce
    else
        lmx = nlmx
        lmn = nlmn
        if ( b/=0.0_wp ) then
            if ( sign(1.0_wp,b)*a>0.0_wp ) then
                c = abs(1.0_wp-a/b)
                if ( c<=0.1_wp ) then
                    if ( c<=0.0_wp ) then
                        if ( err<0.0_wp ) err = ce
                        return
                    else
                        anib = 0.5_wp - log(c)/ln2
                        nib = int(anib)
                        lmx = min(nlmx,nbits-nib-7)
                        if ( lmx<1 ) then
                            ! a and b are too nearly equal to allow normal integration
                            iflag = 1101
                            if ( err<0.0_wp ) err = ce
                            return
                        else
                            lmn = min(lmn,lmx)
                        end if
                    end if
                end if
            end if
        end if
        tol = max(abs(err),2.0_wp**(5-nbits))/2.0_wp
        if ( err==0.0_wp ) tol = sqrt(epsilon(1.0_wp))  ! was d1mach(4)
        eps = tol
        hh(1) = (b-a)/4.0_wp
        aa(1) = a
        lr(1) = 1
        l = 1
        call g8(aa(l)+2.0_wp*hh(l),2.0_wp*hh(l),est,iflag)    !est = g8(aa(l)+2.0_wp*hh(l),2.0_wp*hh(l))
        if (iflag/=0) return
        k = 8
        area = abs(est)
        ef = 0.5_wp
        mxl = 0
    end if

    do
        ! compute refined estimates, estimate the error, etc.
        call g8(aa(l)+hh(l),hh(l),gl,iflag)               ! gl = g8(aa(l)+hh(l),hh(l))
        if (iflag/=0) return
        call g8(aa(l)+3.0_wp*hh(l),hh(l),gr(l),iflag)     ! gr(l) = g8(aa(l)+3.0_wp*hh(l),hh(l))
        if (iflag/=0) return
        k = k + 16
        area = area + (abs(gl)+abs(gr(l))-abs(est))
        glr = gl + gr(l)
        ee = abs(est-glr)*ef
        ae = max(eps*area,tol*abs(glr))
        if ( ee>ae ) then
            ! consider the left half of this level
            if ( k>kmx ) lmx = kml
            if ( l>=lmx ) then
                mxl = 1
            else
                l = l + 1
                eps = eps*0.5_wp
                ef = ef/sq2
                hh(l) = hh(l-1)*0.5_wp
                lr(l) = -1
                aa(l) = aa(l-1)
                est = gl
                cycle
            end if
        end if
        ce = ce + (est-glr)
        if ( lr(l)<=0 ) then
            ! proceed to right half at this level
            vl(l) = glr
        else
            ! return one level
            vr = glr
            do
                if ( l<=1 ) then
                    ! exit
                    ans = vr
                    if ( (mxl/=0) .and. (abs(ce)>2.0_wp*tol*area) ) then
                        iflag = 1102
                    end if
                    if ( err<0.0_wp ) err = ce
                    return
                else
                    l = l - 1
                    eps = eps*2.0_wp
                    ef = ef*sq2
                    if ( lr(l)<=0 ) then
                        vl(l) = vl(l+1) + vr
                        exit
                    else
                        vr = vl(l+1) + vr
                    end if
                end if
            end do
        end if
        est = gr(l-1)
        lr(l) = 1
        aa(l) = aa(l) + 4.0_wp*hh(l)
    end do

    contains

        subroutine g8(x,h,res,iflag)

        implicit none

        real(wp),intent(in)  :: x
        real(wp),intent(in)  :: h
        real(wp),intent(out) :: res
        integer,intent(out)  :: iflag

        integer,dimension(8) :: iflags
        real(wp),dimension(8) :: f

        res = 0.0_wp

        call dbvalu(xt,bc,n,kk,id,x-x1*h,inbv,work,iflag,f(1)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x+x1*h,inbv,work,iflag,f(2)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x-x2*h,inbv,work,iflag,f(3)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x+x2*h,inbv,work,iflag,f(4)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x-x3*h,inbv,work,iflag,f(5)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x+x3*h,inbv,work,iflag,f(6)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x-x4*h,inbv,work,iflag,f(7)); if (iflag/=0) return
        call dbvalu(xt,bc,n,kk,id,x+x4*h,inbv,work,iflag,f(8)); if (iflag/=0) return

        res = h*((w1*(fun(x-x1*h)*f(1)+fun(x+x1*h)*f(2)) + &
                 w2*(fun(x-x2*h)*f(3)+fun(x+x2*h)*f(4))) + &
                (w3*(fun(x-x3*h)*f(5)+fun(x+x3*h)*f(6))  + &
                 w4*(fun(x-x4*h)*f(7)+fun(x+x4*h)*f(8))))

        end subroutine g8

    end subroutine dbsgq8
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the value of `x` to use for computing the interval
!  in `t`, depending on if extrapolation is allowed or not.
!
!  If extrapolation is allowed and x is > t(1) or < t(n), then either
!  t(1) or t(n) is returned. Otherwise, `x` is returned.

    pure function get_temp_x_for_extrap(x,t,extrap) result(xt)

    implicit none

    real(wp),intent(in)              :: x       !! variable value
    real(wp),dimension(:),intent(in) :: t       !! knot vector for b-splines
    real(wp)                         :: xt      !! The value returned (it will either
                                                !! be `t(1)`, `x`, or `t(n)`)
    logical,intent(in),optional      :: extrap  !! if extrapolation is allowed
                                                !! (if not present, default is False)

    integer :: n  !! size of `t`
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    n = size(t)

    if (extrapolation_allowed) then
        if (x<t(1)) then
            xt = t(1)
        else if (x>t(n)) then
            xt = t(n)
        else
            xt = x
        end if
    else
        xt = x
    end if

    end function get_temp_x_for_extrap
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a message string associated with the status code.

    pure function get_status_message(iflag) result(msg)

    implicit none

    integer,intent(in)           :: iflag  !! return code from one of the routines
    character(len=:),allocatable :: msg    !! status message associated with the flag

    character(len=10) :: istr   !! for integer to string conversion
    integer           :: istat  !! for write statement

    select case (iflag)

    case(  0); msg='Successful execution'

    case(  1); msg='Error in evaluate_*d: class is not initialized'

    case(  2); msg='Error in db*ink: iknot out of range'
    case(  3); msg='Error in db*ink: nx out of range'
    case(  4); msg='Error in db*ink: kx out of range'
    case(  5); msg='Error in db*ink: x not strictly increasing'
    case(  6); msg='Error in db*ink: tx not non-decreasing'
    case(  7); msg='Error in db*ink: ny out of range'
    case(  8); msg='Error in db*ink: ky out of range'
    case(  9); msg='Error in db*ink: y not strictly increasing'
    case( 10); msg='Error in db*ink: ty not non-decreasing'
    case( 11); msg='Error in db*ink: nz out of range'
    case( 12); msg='Error in db*ink: kz out of range'
    case( 13); msg='Error in db*ink: z not strictly increasing'
    case( 14); msg='Error in db*ink: tz not non-decreasing'
    case( 15); msg='Error in db*ink: nq out of range'
    case( 16); msg='Error in db*ink: kq out of range'
    case( 17); msg='Error in db*ink: q not strictly increasing'
    case( 18); msg='Error in db*ink: tq not non-decreasing'
    case( 19); msg='Error in db*ink: nr out of range'
    case( 20); msg='Error in db*ink: kr out of range'
    case( 21); msg='Error in db*ink: r not strictly increasing'
    case( 22); msg='Error in db*ink: tr not non-decreasing'
    case( 23); msg='Error in db*ink: ns out of range'
    case( 24); msg='Error in db*ink: ks out of range'
    case( 25); msg='Error in db*ink: s not strictly increasing'
    case( 26); msg='Error in db*ink: ts not non-decreasing'
    case(700); msg='Error in db*ink: size(x) /= size(fcn,1)'
    case(701); msg='Error in db*ink: size(y) /= size(fcn,2)'
    case(702); msg='Error in db*ink: size(z) /= size(fcn,3)'
    case(703); msg='Error in db*ink: size(q) /= size(fcn,4)'
    case(704); msg='Error in db*ink: size(r) /= size(fcn,5)'
    case(705); msg='Error in db*ink: size(s) /= size(fcn,6)'
    case(706); msg='Error in db*ink: size(x) /= nx'
    case(707); msg='Error in db*ink: size(y) /= ny'
    case(708); msg='Error in db*ink: size(z) /= nz'
    case(709); msg='Error in db*ink: size(q) /= nq'
    case(710); msg='Error in db*ink: size(r) /= nr'
    case(711); msg='Error in db*ink: size(s) /= ns'
    case(712); msg='Error in db*ink: size(tx) /= nx+kx'
    case(713); msg='Error in db*ink: size(ty) /= ny+ky'
    case(714); msg='Error in db*ink: size(tz) /= nz+kz'
    case(715); msg='Error in db*ink: size(tq) /= nq+kq'
    case(716); msg='Error in db*ink: size(tr) /= nr+kr'
    case(717); msg='Error in db*ink: size(ts) /= ns+ks'
    case(800); msg='Error in db*ink: size(x) /= size(bcoef,1)'
    case(801); msg='Error in db*ink: size(y) /= size(bcoef,2)'
    case(802); msg='Error in db*ink: size(z) /= size(bcoef,3)'
    case(803); msg='Error in db*ink: size(q) /= size(bcoef,4)'
    case(804); msg='Error in db*ink: size(r) /= size(bcoef,5)'
    case(805); msg='Error in db*ink: size(s) /= size(bcoef,6)'

    case(100); msg='Error in dbintk: k does not satisfy k>=1'
    case(101); msg='Error in dbintk: n does not satisfy n>=k'
    case(102); msg='Error in dbintk: x(i) does not satisfy x(i)<x(i+1) for some i'
    case(103); msg='Error in dbintk: some abscissa was not in the support of the '//&
                    'corresponding basis function and the system is singular'
    case(104); msg='Error in dbintk: the system of solver detects a singular system '//&
                   'although the theoretical conditions for a solution were satisfied'

    case(201); msg='Error in dbspvn: k does not satisfy k>=1'
    case(202); msg='Error in dbspvn: jhigh does not satisfy 1<=jhigh<=k'
    case(203); msg='Error in dbspvn: index is not 1 or 2'
    case(204); msg='Error in dbspvn: x does not satisfy t(ileft)<=x<=t(ileft+1)'

    case(301); msg='Error in dbtpcf: n should be > 0'

    case(401); msg='Error in dbvalu: k does not satisfy k>=1'
    case(402); msg='Error in dbvalu: n does not satisfy n>=k'
    case(403); msg='Error in dbvalu: ideriv does not satisfy 0<=ideriv<k'
    case(404); msg='Error in dbvalu: x is not greater than or equal to t(k)'
    case(405); msg='Error in dbvalu: x is not less than or equal to t(n+1)'
    case(406); msg='Error in dbvalu: a left limiting value cannot be obtained at t(k)'

    case(501); msg='Error in initialize_*d_specify_knots: tx is not the correct size (kx+nx)'
    case(502); msg='Error in initialize_*d_specify_knots: ty is not the correct size (ky+ny)'
    case(503); msg='Error in initialize_*d_specify_knots: tz is not the correct size (kz+nz)'
    case(504); msg='Error in initialize_*d_specify_knots: tq is not the correct size (kq+nq)'
    case(505); msg='Error in initialize_*d_specify_knots: tr is not the correct size (kr+nr)'
    case(506); msg='Error in initialize_*d_specify_knots: ts is not the correct size (ks+ns)'

    case(601); msg='Error in db*val: x value out of bounds'
    case(602); msg='Error in db*val: y value out of bounds'
    case(603); msg='Error in db*val: z value out of bounds'
    case(604); msg='Error in db*val: q value out of bounds'
    case(605); msg='Error in db*val: r value out of bounds'
    case(606); msg='Error in db*val: s value out of bounds'

    case(901); msg='Error in dbsqad: k does not satisfy 1<=k<=20'
    case(902); msg='Error in dbsqad: n does not satisfy n>=k'
    case(903); msg='Error in dbsqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'

    case(1001); msg='Error in dbfqad: k does not satisfy k>=1'
    case(1002); msg='Error in dbfqad: n does not satisfy n>=k'
    case(1003); msg='Error in dbfqad: d does not satisfy 0<=id<k'
    case(1004); msg='Error in dbfqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'
    case(1005); msg='Error in dbfqad: is less than machine epsilon or greater than 0.1'

    case(1101); msg='Warning in dbsgq8: a and b are too nearly equal to allow normal integration.'
    case(1102); msg='Error in dbsgq8: ans is probably insufficiently accurate.'

    case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select

    end function get_status_message
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_sub_module
!*****************************************************************************************
