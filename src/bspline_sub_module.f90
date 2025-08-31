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

    use bspline_kinds_module, only: wp, ip
    use,intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    abstract interface
        function b1fqad_func(x) result(f)
        !! interface for the input function in [[dbfqad]]
        import :: wp
        implicit none
        real(wp),intent(in) :: x
        real(wp)            :: f  !! f(x)
        end function b1fqad_func
    end interface
    public :: b1fqad_func

    integer(ip),parameter,public :: bspline_order_linear    = 2_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_quadratic = 3_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_cubic     = 4_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_quartic   = 5_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_quintic   = 6_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_hexic     = 7_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_heptic    = 8_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]
    integer(ip),parameter,public :: bspline_order_octic     = 9_ip !! spline order `k` parameter
                                                                   !! (for input to the `db*ink` routines)
                                                                   !! [order = polynomial degree + 1]

    interface db1ink
        !! 1D initialization routines.
        module procedure :: db1ink_default, db1ink_alt, db1ink_alt_2
    end interface
    interface db1val
        !! 1D evaluation routines.
        module procedure :: db1val_default, db1val_alt
    end interface

    !main routines:
    public :: db1ink, db1val, db1sqad, db1fqad
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

    pure subroutine db1ink_default(x,nx,fcn,kx,iknot,tx,bcoef,iflag)

    implicit none

    integer(ip),intent(in)                  :: nx     !! Number of \(x\) abcissae
    integer(ip),intent(in)                  :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)        :: fcn    !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                                      !! contain the function value at the point `x(i)`
    integer(ip),intent(in)                  :: iknot  !! knot sequence flag:
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
    integer(ip),intent(out)                 :: iflag  !! status flag:
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: work   !! work array of dimension `2*kx*(nx+1)`

    !check validity of inputs

    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
        end if

        allocate(work(2_ip*kx*(nx+1_ip)))

        !construct b-spline coefficients
        call dbtpcf(x,nx,fcn,nx,1_ip,tx,kx,bcoef,work,iflag)

        deallocate(work)

    end if

    end subroutine db1ink_default
!*****************************************************************************************

!*****************************************************************************************
!>
!  Alternate version of [[db1ink_default]], where the boundary conditions can be specified.
!
!### History
!  * Jacob Williams, 9/4/2018 : created this routine.
!
!### See also
!  * [[dbint4]] -- the main routine that is called here.
!
!@note Currently, this only works for 3rd order (k=4).

    pure subroutine db1ink_alt(x,nx,fcn,kx,ibcl,ibcr,fbcl,fbcr,kntopt,tx,bcoef,iflag)

    implicit none

    real(wp),dimension(:),intent(in)   :: x       !! \(x\) vector of abscissae of length `nx`, distinct
                                                  !! and in increasing order
    integer(ip),intent(in)             :: nx      !! number of data points, \( n_x \ge 2 \)
    real(wp),dimension(:),intent(in)   :: fcn     !! \(y\) vector of ordinates of length `nx`
    integer(ip),intent(in)             :: kx      !! spline order (Currently, this must be `4`)
    integer(ip),intent(in)             :: ibcl    !! selection parameter for left boundary condition:
                                                  !!
                                                  !! * `ibcl = 1` constrain the first derivative at `x(1)` to `fbcl`
                                                  !! * `ibcl = 2` constrain the second derivative at `x(1)` to `fbcl`
    integer(ip),intent(in)             :: ibcr    !! selection parameter for right boundary condition:
                                                  !!
                                                  !! * `ibcr = 1` constrain first derivative at `x(nx)` to `fbcr`
                                                  !! * `ibcr = 2` constrain second derivative at `x(nx)` to `fbcr`
    real(wp),intent(in)                :: fbcl    !! left boundary values governed by `ibcl`
    real(wp),intent(in)                :: fbcr    !! right boundary values governed by `ibcr`
    integer(ip),intent(in)             :: kntopt  !! knot selection parameter:
                                                  !!
                                                  !! * `kntopt = 1` sets knot multiplicity at `t(4)` and
                                                  !!   `t(nx+3)` to 4
                                                  !! * `kntopt = 2` sets a symmetric placement of knots
                                                  !!   about `t(4)` and `t(nx+3)`
    real(wp),dimension(:),intent(out)  :: tx      !! knot array of length `nx+6`
    real(wp),dimension(:),intent(out)  :: bcoef   !! b spline coefficient array of length `nx+2`
    integer(ip),intent(out)            :: iflag   !! status flag:
                                                  !!
                                                  !! * 0: no errors
                                                  !! * 806: [[dbint4]] can only be used when `k=4`

    real(wp),dimension(:,:),allocatable :: w         !! work array of dimension `5,nx+2`
    integer(ip)                         :: n         !! number of coefficients (n=nx+2)
    integer(ip)                         :: k         !! order of spline (k=4)
    logical                             :: status_ok !! status flag for error checking

    real(wp),dimension(3),parameter :: tleft = 0.0_wp   !! not used for this case (see [[dbint4]])
    real(wp),dimension(3),parameter :: tright = 0.0_wp  !! not used for this case (see [[dbint4]])


    if (kx /= 4_ip) then
        iflag = 806_ip
    else

        call check_inputs(  1_ip,& ! so it will check size of t
                            iflag,&
                            nx=nx,&
                            kx=kx,&
                            x=x,&
                            f1=fcn,&
                            bcoef1=bcoef,&
                            tx=tx,&
                            status_ok=status_ok,&
                            alt=.true.)

        if (status_ok) then
            allocate(w(5_ip,nx+2_ip))
            call dbint4(x,fcn,nx,ibcl,ibcr,fbcl,fbcr,kntopt,tleft,tright,tx,bcoef,n,k,w,iflag)
            deallocate(w)
        end if

    end if

    end subroutine db1ink_alt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Alternate version of [[db1ink_alt]], where the first and
!  last 3 knots are specified by the user.
!
!### History
!  * Jacob Williams, 9/4/2018 : created this routine.
!
!### See also
!  * [[dbint4]] -- the main routine that is called here.
!
!@note Currently, this only works for 3rd order (k=4).

    pure subroutine db1ink_alt_2(x,nx,fcn,kx,ibcl,ibcr,fbcl,fbcr,tleft,tright,tx,bcoef,iflag)

    implicit none

    real(wp),dimension(:),intent(in)   :: x       !! \(x\) vector of abscissae of length `nx`, distinct
                                                  !! and in increasing order
    integer(ip),intent(in)             :: nx      !! number of data points, \( n_x \ge 2 \)
    real(wp),dimension(:),intent(in)   :: fcn     !! \(y\) vector of ordinates of length `nx`
    integer(ip),intent(in)             :: kx      !! spline order (Currently, this must be `4`)
    integer(ip),intent(in)             :: ibcl    !! selection parameter for left boundary condition:
                                                  !!
                                                  !! * `ibcl = 1` constrain the first derivative at `x(1)` to `fbcl`
                                                  !! * `ibcl = 2` constrain the second derivative at `x(1)` to `fbcl`
    integer(ip),intent(in)             :: ibcr    !! selection parameter for right boundary condition:
                                                  !!
                                                  !! * `ibcr = 1` constrain first derivative at `x(nx)` to `fbcr`
                                                  !! * `ibcr = 2` constrain second derivative at `x(nx)` to `fbcr`
    real(wp),intent(in)                :: fbcl    !! left boundary values governed by `ibcl`
    real(wp),intent(in)                :: fbcr    !! right boundary values governed by `ibcr`
    real(wp),dimension(3),intent(in)   :: tleft   !! `t(1:3)` in increasing order supplied by the user.
    real(wp),dimension(3),intent(in)   :: tright  !! `t(nx+4:nx+6)` in increasing order supplied by the user.
    real(wp),dimension(:),intent(out)  :: tx      !! knot array of length `nx+6`
    real(wp),dimension(:),intent(out)  :: bcoef   !! b spline coefficient array of length `nx+2`
    integer(ip),intent(out)            :: iflag   !! status flag:
                                                  !!
                                                  !! * 0: no errors
                                                  !! * 806: [[dbint4]] can only be used when k=4

    real(wp),dimension(:,:),allocatable :: w         !! work array of dimension `5,nx+2`
    integer(ip)                         :: n         !! number of coefficients (`n=nx+2`)
    integer(ip)                         :: k         !! order of spline (`k=4`)
    logical                             :: status_ok !! status flag for error checking

    integer(ip),parameter :: kntopt = 3 !! use `tleft` and `tright` in [[dbint4]]

    if (kx /= 4_ip) then
        iflag = 806_ip
    else

        call check_inputs(  1_ip,& ! so it will check size of t
                            iflag,&
                            nx=nx,&
                            kx=kx,&
                            x=x,&
                            f1=fcn,&
                            bcoef1=bcoef,&
                            tx=tx,&
                            status_ok=status_ok,&
                            alt=.true.)

        if (status_ok) then
            allocate(w(5,nx+2))
            call dbint4(x,fcn,nx,ibcl,ibcr,fbcl,fbcr,kntopt,tleft,tright,tx,bcoef,n,k,w,iflag)
            deallocate(w)
        end if

    end if

    end subroutine db1ink_alt_2
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

    pure subroutine db1val_default(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

    implicit none

    integer(ip),intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    integer(ip),intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction. (same as in last call to [[db1ink]])
    real(wp),dimension(nx),intent(in)    :: bcoef    !! the b-spline coefficients computed by [[db1ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer(ip),intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer(ip),intent(inout)            :: inbvx    !! initialization parameter which must be set
                                                     !! to 1 the first time this routine is called,
                                                     !! and must not be changed by the user.
    real(wp),dimension(3_ip*kx),intent(inout) :: w0  !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return

    call dbvalu(tx,bcoef,nx,kx,idx,xval,inbvx,w0,iflag,f,extrap)

    end subroutine db1val_default
!*****************************************************************************************

!*****************************************************************************************
!>
!  Alternate version of [[db1val_default]] for use with [[db1ink_alt]] and [[db1ink_alt_2]].

    pure subroutine db1val_alt(xval,idx,tx,nx,n,kx,bcoef,f,iflag,inbvx,w0,extrap)

    implicit none

    real(wp),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    integer(ip),intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: nx       !! the number of interpolation points in \(x\).
    integer(ip),intent(in)               :: n        !! length of `bcoef`: `nx+2`
    integer(ip),intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db1ink]])
    real(wp),dimension(n+kx),intent(in)  :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction.
    real(wp),dimension(n),intent(in)     :: bcoef    !! the b-spline coefficients computed by [[db1ink]].
    real(wp),intent(out)                 :: f        !! interpolated value
    integer(ip),intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer(ip),intent(inout)            :: inbvx    !! initialization parameter which must be set
                                                     !! to 1 the first time this routine is called,
                                                     !! and must not be changed by the user.
    real(wp),dimension(3_ip*kx),intent(inout) :: w0  !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return

    call dbvalu(tx,bcoef,n,kx,idx,xval,inbvx,w0,iflag,f,extrap)

    end subroutine db1val_alt
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

    pure subroutine db1sqad(tx,bcoef,nx,kx,x1,x2,f,iflag,w0)

    implicit none

    integer(ip),intent(in)               :: nx      !! length of coefficient array
    integer(ip),intent(in)               :: kx      !! order of b-spline, `1 <= k <= 20`
    real(wp),dimension(nx+kx),intent(in) :: tx      !! knot array
    real(wp),dimension(nx),intent(in)    :: bcoef   !! b-spline coefficient array
    real(wp),intent(in)                  :: x1      !! left point of quadrature interval in `t(kx) <= x <= t(nx+1)`
    real(wp),intent(in)                  :: x2      !! right point of quadrature interval in `t(kx) <= x <= t(nx+1)`
    real(wp),intent(out)                 :: f       !! integral of the b-spline over (`x1`,`x2`)
    integer(ip),intent(out)              :: iflag   !! status flag:
                                                    !!
                                                    !! * \( = 0 \)   : no errors
                                                    !! * \( \ne 0 \) : error
    real(wp),dimension(3*kx),intent(inout) :: w0    !! work array for [[dbsqad]]

    call dbsqad(tx,bcoef,nx,kx,x1,x2,f,w0,iflag)

    end subroutine db1sqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the integral on `(x1,x2)` of a product of a
!  function `fun` and the `idx`-th derivative of a `kx`-th order b-spline,
!  using the b-representation `(tx,bcoef,nx,kx)`, with an adaptive
!  8-point Legendre-Gauss algorithm.
!  `(x1,x2)` must be a subinterval of `t(kx) <= x <= t(nx+1)`.
!
!### See also
!  * [[dbfqad]] -- the core routine.
!
!@note This one is not pure, because we are not enforcing
!      that the user function `fun` be pure.

    subroutine db1fqad(fun,tx,bcoef,nx,kx,idx,x1,x2,tol,f,iflag,w0)

    implicit none

    procedure(b1fqad_func)              :: fun    !! external function of one argument for the
                                                  !! integrand `bf(x)=fun(x)*dbvalu(tx,bcoef,nx,kx,id,x,inbv,work)`
    integer(ip),intent(in)              :: nx     !! length of coefficient array
    integer(ip),intent(in)              :: kx     !! order of b-spline, `kx >= 1`
    real(wp),dimension(nx+kx),intent(in):: tx     !! knot array
    real(wp),dimension(nx),intent(in)   :: bcoef  !! b-spline coefficient array
    integer(ip),intent(in)              :: idx    !! order of the spline derivative, `0 <= idx <= k-1`
                                                  !! `idx=0` gives the spline function
    real(wp),intent(in)                 :: x1     !! left point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: x2     !! right point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: tol    !! desired accuracy for the quadrature, suggest
                                                  !! `10*dtol < tol <= 0.1` where `dtol` is the maximum
                                                  !! of `1.0e-300` and real(wp) unit roundoff for
                                                  !! the machine
    real(wp),intent(out)                :: f      !! integral of `bf(x)` on `(x1,x2)`
    integer(ip),intent(out)             :: iflag  !! status flag:
                                                  !!
                                                  !! * \( = 0 \)   : no errors
                                                  !! * \( \ne 0 \) : error
    real(wp),dimension(3_ip*kx),intent(inout) :: w0  !! work array for [[dbfqad]]

    call dbfqad(fun,tx,bcoef,nx,kx,idx,x1,x2,tol,f,iflag,w0)

    end subroutine db1fqad
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

    integer(ip),intent(in)                  :: nx     !! Number of \(x\) abcissae
    integer(ip),intent(in)                  :: ny     !! Number of \(y\) abcissae
    integer(ip),intent(in)                  :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                  :: ky     !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)        :: y      !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in)      :: fcn    !! `(nx,ny)` matrix of function values to interpolate.
                                                      !! `fcn(i,j)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)                  :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)     :: ty     !! The `(ny+ky)` knots in the \(y\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:,:),intent(out)     :: bcoef  !! `(nx,ny)` matrix of coefficients of the b-spline interpolant.
    integer(ip),intent(out)                 :: iflag  !! *  0 = successful execution.
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of length `nx*ny`
    real(wp),dimension(:),allocatable :: work !! work array of length `max(2*kx*(nx+1),2*ky*(ny+1))`

    !check validity of inputs

    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
        end if

        allocate(temp(nx*ny))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip))))

        !construct b-spline coefficients
                         call dbtpcf(x,nx,fcn, nx,ny,tx,kx,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,temp,ny,nx,ty,ky,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

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

    pure subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)               :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)               :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: ny       !! the number of interpolation points in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer(ip),intent(in)               :: ky       !! order of polynomial pieces in \(y\).
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
    integer(ip),intent(out)              :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer(ip),intent(inout)            :: inbvx    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer(ip),intent(inout)            :: inbvy    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer(ip),intent(inout)            :: iloy     !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    real(wp),dimension(ky),intent(inout)              :: w1 !! work array
    real(wp),dimension(3_ip*max(kx,ky)),intent(inout) :: w0 !! work array
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    integer(ip) :: k, lefty, kcol

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return

    kcol = lefty - ky
    do k=1_ip,ky
        kcol = kcol + 1_ip
        call dbvalu(tx,bcoef(:,kcol),nx,kx,idx,xval,inbvx,w0,iflag,w1(k),extrap)
        if (iflag/=0_ip) return !error
    end do

    kcol = lefty - ky + 1_ip
    call dbvalu(ty(kcol:),w1,ky,ky,idy,yval,inbvy,w0,iflag,f,extrap)

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

    integer(ip),intent(in)                   :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                   :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                   :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                   :: kx    !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                   :: ky    !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                   :: kz    !! the order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)         :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in)     :: fcn   !! `(nx,ny,nz)` matrix of function values to interpolate. `fcn(i,j,k)` should
                                                      !! contain the function value at the point (`x(i)`,`y(j)`,`z(k)`)
    integer(ip),intent(in)                   :: iknot !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db3ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)      :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the spline
                                                      !! interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db3ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(wp),dimension(:,:,:),intent(out)    :: bcoef !! `(nx,ny,nz)` matrix of coefficients of the b-spline interpolant.
    integer(ip),intent(out)                  :: iflag !! *  0 = successful execution.
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of length `nx*ny*nz`
    real(wp),dimension(:),allocatable :: work !! work array of length `max(2*kx*(nx+1),
                                              !! 2*ky*(ny+1),2*kz*(nz+1))`
    integer(ip) :: i, j, k, ii !! counter

    ! check validity of input

    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
        end if

        allocate(temp(nx*ny*nz))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip),2_ip*kz*(nz+1_ip))))

        ! copy fcn to work in packed for dbtpcf
        !temp = reshape( fcn, [nx*ny*nz] )
        ! replaced with loops to avoid stack
        ! overflow for large data set:
        ii = 0_ip
        do k = 1_ip, nz
            do j = 1_ip, ny
                do i = 1_ip, nx
                    ii = ii + 1_ip
                    temp(ii) = fcn(i,j,k)
                end do
            end do
        end do

        ! construct b-spline coefficients
                         call dbtpcf(x,nx,temp, nx,ny*nz,tx,kx,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,bcoef,ny,nx*nz,ty,ky,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(z,nz,temp, nz,nx*ny,tz,kz,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

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
                           inbvx,inbvy,inbvz,iloy,iloz,w2,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)                  :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                  :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                  :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                  :: nx       !! the number of interpolation points in \(x\).
                                                        !! (same as in last call to [[db3ink]])
    integer(ip),intent(in)                  :: ny       !! the number of interpolation points in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer(ip),intent(in)                  :: nz       !! the number of interpolation points in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer(ip),intent(in)                  :: kx       !! order of polynomial pieces in \(z\).
                                                        !! (same as in last call to [[db3ink]])
    integer(ip),intent(in)                  :: ky       !! order of polynomial pieces in \(y\).
                                                        !! (same as in last call to [[db3ink]])
    integer(ip),intent(in)                  :: kz       !! order of polynomial pieces in \(z\).
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
    integer(ip),intent(out)                 :: iflag    !! status flag:
                                                        !!
                                                        !! * \( = 0 \)   : no errors
                                                        !! * \( \ne 0 \) : error
    integer(ip),intent(inout)               :: inbvx    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer(ip),intent(inout)               :: inbvy    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer(ip),intent(inout)               :: inbvz    !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer(ip),intent(inout)               :: iloy     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    integer(ip),intent(inout)               :: iloz     !! initialization parameter which must be
                                                        !! set to 1 the first time this routine is called,
                                                        !! and must not be changed by the user.
    real(wp),dimension(ky,kz),intent(inout)              :: w2  !! work array
    real(wp),dimension(kz),intent(inout)                 :: w1  !! work array
    real(wp),dimension(3_ip*max(kx,ky,kz)),intent(inout) :: w0  !! work array
    logical,intent(in),optional             :: extrap   !! if extrapolation is allowed
                                                        !! (if not present, default is False)

    integer(ip) :: lefty, leftz, kcoly, kcolz, j, k

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(zval,tz,3_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,iflag,extrap); if (iflag/=0_ip) return

    iflag = 0_ip

    kcolz = leftz - kz
    do k=1_ip,kz
        kcolz = kcolz + 1_ip
        kcoly = lefty - ky
        do j=1_ip,ky
            kcoly = kcoly + 1_ip
            call dbvalu(tx,bcoef(:,kcoly,kcolz),nx,kx,idx,xval,inbvx,w0,iflag,w2(j,k),extrap)
            if (iflag/=0_ip) return
        end do
    end do

    kcoly = lefty - ky + 1_ip
    do k=1_ip,kz
        call dbvalu(ty(kcoly:),w2(:,k),ky,ky,idy,yval,inbvy,w0,iflag,w1(k),extrap)
        if (iflag/=0_ip) return
    end do

    kcolz = leftz - kz + 1_ip
    call dbvalu(tz(kcolz:),w1,kz,kz,idz,zval,inbvz,w0,iflag,f,extrap)

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

    integer(ip),intent(in)                      :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                      :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                      :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                      :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                      :: kx    !! the order of spline pieces in \(x\)
                                                         !! ( \( 2 \le k_x < n_x \) ).
                                                         !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: ky    !! the order of spline pieces in \(y\)
                                                         !! ( \( 2 \le k_y < n_y \) ).
                                                         !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: kz    !! the order of spline pieces in \(z\)
                                                         !! ( \( 2 \le k_z < n_z \) ).
                                                         !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: kq    !! the order of spline pieces in \(q\)
                                                         !! ( \( 2 \le k_q < n_q \) ).
                                                         !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)            :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)            :: q     !! `(nq)` array of \(q\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)      :: fcn   !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                         !! `fcn(i,j,k,q)` should contain the function value at the
                                                         !!  point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer(ip),intent(in)                      :: iknot !! knot sequence flag:
                                                         !!
                                                         !! * 0 = knot sequence chosen by [[db4ink]].
                                                         !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)         :: tx    !! The `(nx+kx)` knots in the x direction for the spline
                                                         !! interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: ty    !! The `(ny+ky)` knots in the y direction for the spline
                                                         !! interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: tz    !! The `(nz+kz)` knots in the z direction for the spline
                                                         !! interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)         :: tq    !! The `(nq+kq)` knots in the q direction for the spline
                                                         !! interpolant.
                                                         !!
                                                         !! * If `iknot=0` these are chosen by [[db4ink]].
                                                         !! * If `iknot=1` these are specified by the user.
                                                         !!
                                                         !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:),intent(out)     :: bcoef !! `(nx,ny,nz,nq)` matrix of coefficients of the b-spline
                                                         !! interpolant.
    integer(ip),intent(out)                     :: iflag !! *  0 = successful execution.
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of dimension `nx*ny*nz*nq`
    real(wp),dimension(:),allocatable :: work !! work array of dimension `max(2*kx*(nx+1),
                                              !! 2*ky*(ny+1),2*kz*(nz+1),2*kq*(nq+1))`

    ! check validity of input

    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
        end if

        allocate(temp(nx*ny*nz*nq))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip),2_ip*kz*(nz+1_ip),2_ip*kq*(nq+1_ip))))

        ! construct b-spline coefficients
                         call dbtpcf(x,nx,fcn,  nx,ny*nz*nq,tx,kx,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,temp, ny,nx*nz*nq,ty,ky,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(z,nz,bcoef,nz,nx*ny*nq,tz,kz,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(q,nq,temp, nq,nx*ny*nz,tq,kq,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

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
                           iloy,iloz,iloq,w3,w2,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)                     :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                     :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                     :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                     :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                     :: nx       !! the number of interpolation points in \(x\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: ny       !! the number of interpolation points in \(y\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: nz       !! the number of interpolation points in \(z\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: nq       !! the number of interpolation points in \(q\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: kx       !! order of polynomial pieces in \(x\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: ky       !! order of polynomial pieces in \(y\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: kz       !! order of polynomial pieces in \(z\).
                                                           !! (same as in last call to [[db4ink]])
    integer(ip),intent(in)                     :: kq       !! order of polynomial pieces in \(q\).
                                                           !! (same as in last call to [[db4ink]])
    real(wp),intent(in)                        :: xval     !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)                        :: yval     !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)                        :: zval     !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)                        :: qval     !! \(q\) coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)       :: tx       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(x\) direction. (same as in last call to
                                                           !! [[db4ink]])
    real(wp),dimension(ny+ky),intent(in)       :: ty       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(y\) direction. (same as in last call to
                                                           !! [[db4ink]])
    real(wp),dimension(nz+kz),intent(in)       :: tz       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(z\) direction. (same as in last call to
                                                           !! [[db4ink]])
    real(wp),dimension(nq+kq),intent(in)       :: tq       !! sequence of knots defining the piecewise polynomial
                                                           !! in the \(q\) direction. (same as in last call to
                                                           !! [[db4ink]])
    real(wp),dimension(nx,ny,nz,nq),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db4ink]].
    real(wp),intent(out)                       :: f        !! interpolated value
    integer(ip),intent(out)                    :: iflag    !! status flag:
                                                           !!
                                                           !! * \( = 0 \)   : no errors
                                                           !! * \( \ne 0 \) : error
    integer(ip),intent(inout)                  :: inbvx    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: inbvy    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: inbvz    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: inbvq    !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: iloy     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: iloz     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    integer(ip),intent(inout)                  :: iloq     !! initialization parameter which must be set
                                                           !! to 1 the first time this routine is called,
                                                           !! and must not be changed by the user.
    real(wp),dimension(ky,kz,kq),intent(inout)           :: w3 !! work array
    real(wp),dimension(kz,kq),intent(inout)              :: w2 !! work array
    real(wp),dimension(kq),intent(inout)                 :: w1 !! work array
    real(wp),dimension(3_ip*max(kx,ky,kz,kq)),intent(inout) :: w0 !! work array
    logical,intent(in),optional                :: extrap   !! if extrapolation is allowed
                                                           !! (if not present, default is False)

    integer(ip) :: lefty, leftz, leftq, &
               kcoly, kcolz, kcolq, j, k, q

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(zval,tz,3_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(qval,tq,4_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,iflag,extrap); if (iflag/=0_ip) return

    iflag = 0_ip

    ! x -> y, z, q
    kcolq = leftq - kq
    do q=1_ip,kq
        kcolq = kcolq + 1_ip
        kcolz = leftz - kz
        do k=1_ip,kz
            kcolz = kcolz + 1_ip
            kcoly = lefty - ky
            do j=1_ip,ky
                kcoly = kcoly + 1_ip
                call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq),&
                                     nx,kx,idx,xval,inbvx,w0,iflag,&
                                     w3(j,k,q),extrap)
                if (iflag/=0_ip) return
            end do
        end do
    end do

    ! y -> z, q
    kcoly = lefty - ky + 1_ip
    do q=1_ip,kq
        do k=1_ip,kz
            call dbvalu(ty(kcoly:),w3(:,k,q),&
                        ky,ky,idy,yval,inbvy,w0,iflag,&
                        w2(k,q),extrap)
            if (iflag/=0_ip) return
        end do
    end do

    ! z -> q
    kcolz = leftz - kz + 1_ip
    do q=1_ip,kq
        call dbvalu(tz(kcolz:),w2(:,q),&
                    kz,kz,idz,zval,inbvz,w0,iflag,&
                    w1(q),extrap)
        if (iflag/=0_ip) return
    end do

    ! q
    kcolq = leftq - kq + 1_ip
    call dbvalu(tq(kcolq:),w1,kq,kq,idq,qval,inbvq,w0,iflag,f,extrap)

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

    integer(ip),intent(in)                         :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                         :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                         :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                         :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                         :: nr    !! number of \(r\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                         :: kx    !! the order of spline pieces in \(x\)
                                                            !! ( \( 2 \le k_x < n_x \) ).
                                                            !! (order = polynomial degree + 1)
    integer(ip),intent(in)                         :: ky    !! the order of spline pieces in \(y\)
                                                            !! ( \( 2 \le k_y < n_y \) ).
                                                            !! (order = polynomial degree + 1)
    integer(ip),intent(in)                         :: kz    !! the order of spline pieces in \(z\)
                                                            !! ( \( 2 \le k_z < n_z \) ).
                                                            !! (order = polynomial degree + 1)
    integer(ip),intent(in)                         :: kq    !! the order of spline pieces in \(q\)
                                                            !! ( \( 2 \le k_q < n_q \) ).
                                                            !! (order = polynomial degree + 1)
    integer(ip),intent(in)                         :: kr    !! the order of spline pieces in \(r\)
                                                            !! ( \( 2 \le k_r < n_r \) ).
                                                            !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)               :: x     !! `(nx)` array of \(x\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: y     !! `(ny)` array of \(y\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: z     !! `(nz)` array of \(z\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: q     !! `(nq)` array of \(q\) abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)               :: r     !! `(nr)` array of \(r\) abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)       :: fcn   !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate.
                                                            !! `fcn(i,j,k,q,r)` should contain the function value at the
                                                            !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer(ip),intent(in)                         :: iknot !! knot sequence flag:
                                                            !!
                                                            !! * 0 = knot sequence chosen by [[db5ink]].
                                                            !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)            :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the spline
                                                            !! interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the spline
                                                            !! interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the spline
                                                            !! interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tq    !! The `(nq+kq)` knots in the \(q\) direction for the spline
                                                            !! interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)            :: tr    !! The `(nr+kr)` knots in the \(r\) direction for the spline
                                                            !! interpolant.
                                                            !!
                                                            !! * If `iknot=0` these are chosen by [[db5ink]].
                                                            !! * If `iknot=1` these are specified by the user.
                                                            !!
                                                            !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:,:),intent(out)      :: bcoef !! `(nx,ny,nz,nq,nr)` matrix of coefficients of the b-spline
                                                            !! interpolant.
    integer(ip),intent(out)                        :: iflag !! *  0 = successful execution.
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of length `nx*ny*nz*nq*nr`
    real(wp),dimension(:),allocatable :: work !! work array of length `max(2*kx*(nx+1),
                                              !! 2*ky*(ny+1),2*kz*(nz+1),2*kq*(nq+1),2*kr*(nr+1))`
    integer(ip) :: i, j, k, l, m, ii !! counter

    !  check validity of input
    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
        end if

        allocate(temp(nx*ny*nz*nq*nr))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip),2_ip*kz*(nz+1_ip),2_ip*kq*(nq+1_ip),2_ip*kr*(nr+1_ip))))

        ! copy fcn to work in packed for dbtpcf
        !temp(1:nx*ny*nz*nq*nr) = reshape( fcn, [nx*ny*nz*nq*nr] )
        ! replaced with loops to avoid stack
        ! overflow for large data set:
        ii = 0_ip
        do m = 1_ip, nr
            do l = 1_ip, nq
                do k = 1_ip, nz
                    do j = 1_ip, ny
                        do i = 1_ip, nx
                            ii = ii + 1_ip
                            temp(ii) = fcn(i,j,k,l,m)
                        end do
                    end do
                end do
            end do
        end do

        !  construct b-spline coefficients
                         call dbtpcf(x,nx,temp,  nx,ny*nz*nq*nr,tx,kx,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,bcoef, ny,nx*nz*nq*nr,ty,ky,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(z,nz,temp,  nz,nx*ny*nq*nr,tz,kz,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(q,nq,bcoef, nq,nx*ny*nz*nr,tq,kq,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(r,nr,temp,  nr,nx*ny*nz*nq,tr,kr,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

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
                           iloy,iloz,iloq,ilor,&
                           w4,w3,w2,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)                        :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                        :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                        :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                        :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                        :: idr      !! \(r\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                        :: nx       !! the number of interpolation points in \(x\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: ny       !! the number of interpolation points in \(y\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: nz       !! the number of interpolation points in \(z\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: nq       !! the number of interpolation points in \(q\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: nr       !! the number of interpolation points in \(r\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: kx       !! order of polynomial pieces in \(x\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: ky       !! order of polynomial pieces in \(y\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: kz       !! order of polynomial pieces in \(z\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: kq       !! order of polynomial pieces in \(q\).
                                                              !! (same as in last call to [[db5ink]])
    integer(ip),intent(in)                        :: kr       !! order of polynomial pieces in \(r\).
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
    integer(ip),intent(out)                       :: iflag    !! status flag:
                                                              !!
                                                              !! * \( = 0 \)   : no errors
                                                              !! * \( \ne 0 \) : error
    integer(ip),intent(inout)                     :: inbvx    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: inbvy    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: inbvz    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: inbvq    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: inbvr    !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: iloy     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: iloz     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: iloq     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    integer(ip),intent(inout)                     :: ilor     !! initialization parameter which must be set
                                                              !! to 1 the first time this routine is called,
                                                              !! and must not be changed by the user.
    real(wp),dimension(ky,kz,kq,kr),intent(inout)           :: w4  !! work array
    real(wp),dimension(kz,kq,kr),intent(inout)              :: w3  !! work array
    real(wp),dimension(kq,kr),intent(inout)                 :: w2  !! work array
    real(wp),dimension(kr),intent(inout)                    :: w1  !! work array
    real(wp),dimension(3_ip*max(kx,ky,kz,kq,kr)),intent(inout) :: w0  !! work array
    logical,intent(in),optional                   :: extrap   !! if extrapolation is allowed
                                                              !! (if not present, default is False)

    integer(ip) :: lefty, leftz, leftq, leftr, &
               kcoly, kcolz, kcolq, kcolr, j, k, q, r

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(zval,tz,3_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(qval,tq,4_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(rval,tr,5_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,iflag,extrap); if (iflag/=0_ip) return

    iflag = 0_ip

    ! x -> y, z, q, r
    kcolr = leftr - kr
    do r=1_ip,kr
        kcolr = kcolr + 1_ip
        kcolq = leftq - kq
        do q=1_ip,kq
            kcolq = kcolq + 1_ip
            kcolz = leftz - kz
            do k=1_ip,kz
                kcolz = kcolz + 1_ip
                kcoly = lefty - ky
                do j=1_ip,ky
                    kcoly = kcoly + 1_ip
                    call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr),&
                                nx,kx,idx,xval,inbvx,w0,iflag,w4(j,k,q,r),&
                                extrap)
                    if (iflag/=0_ip) return
                end do
            end do
        end do
    end do

    ! y -> z, q, r
    kcoly = lefty - ky + 1_ip
    do r=1_ip,kr
        do q=1_ip,kq
            do k=1_ip,kz
                call dbvalu(ty(kcoly:),w4(:,k,q,r),ky,ky,idy,yval,inbvy,&
                            w0,iflag,w3(k,q,r),extrap)
                if (iflag/=0_ip) return
            end do
        end do
    end do

    ! z -> q, r
    kcolz = leftz - kz + 1_ip
    do r=1_ip,kr
        do q=1_ip,kq
            call dbvalu(tz(kcolz:),w3(:,q,r),kz,kz,idz,zval,inbvz,&
                        w0,iflag,w2(q,r),extrap)
            if (iflag/=0_ip) return
        end do
    end do

    ! q -> r
    kcolq = leftq - kq + 1_ip
    do r=1_ip,kr
        call dbvalu(tq(kcolq:),w2(:,r),kq,kq,idq,qval,inbvq,&
                    w0,iflag,w1(r),extrap)
        if (iflag/=0_ip) return
    end do

    ! r
    kcolr = leftr - kr + 1_ip
    call dbvalu(tr(kcolr:),w1,kr,kr,idr,rval,inbvr,w0,iflag,f,extrap)

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

    integer(ip),intent(in)                            :: nx    !! number of \(x\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: ny    !! number of \(y\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: nz    !! number of \(z\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: nq    !! number of \(q\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: nr    !! number of \(r\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: ns    !! number of \(s\) abcissae ( \( \ge 3 \) )
    integer(ip),intent(in)                            :: kx    !! the order of spline pieces in \(x\)
                                                               !! ( \( 2 \le k_x < n_x \) )
                                                               !! (order = polynomial degree + 1)
    integer(ip),intent(in)                            :: ky    !! the order of spline pieces in \(y\)
                                                               !! ( \( 2 \le k_y < n_y \) )
                                                               !! (order = polynomial degree + 1)
    integer(ip),intent(in)                            :: kz    !! the order of spline pieces in \(z\)
                                                               !! ( \( 2 \le k_z < n_z \) )
                                                               !! (order = polynomial degree + 1)
    integer(ip),intent(in)                            :: kq    !! the order of spline pieces in \(q\)
                                                               !! ( \( 2 \le k_q < n_q \) )
                                                               !! (order = polynomial degree + 1)
    integer(ip),intent(in)                            :: kr    !! the order of spline pieces in \(r\)
                                                               !! ( \( 2 \le k_r < n_r \) )
                                                               !! (order = polynomial degree + 1)
    integer(ip),intent(in)                            :: ks    !! the order of spline pieces in \(s\)
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
    real(wp),dimension(:,:,:,:,:,:),intent(in)        :: fcn   !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to
                                                               !! interpolate. `fcn(i,j,k,q,r,s)` should contain the
                                                               !! function value at the point
                                                               !! (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer(ip),intent(in)                            :: iknot !! knot sequence flag:
                                                               !!
                                                               !! * 0 = knot sequence chosen by [[db6ink]].
                                                               !! * 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)               :: tx    !! The `(nx+kx)` knots in the \(x\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * f `iknot=0` these are chosen by [[db6ink]].
                                                               !! * f `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: ty    !! The `(ny+ky)` knots in the \(y\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tz    !! The `(nz+kz)` knots in the \(z\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tq    !! The `(nq+kq)` knots in the \(q\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: tr    !! The `(nr+kr)` knots in the \(r\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:),intent(inout)               :: ts    !! The `(ns+ks)` knots in the \(s\) direction for the
                                                               !! spline interpolant.
                                                               !!
                                                               !! * If `iknot=0` these are chosen by [[db6ink]].
                                                               !! * If `iknot=1` these are specified by the user.
                                                               !!
                                                               !! Must be non-decreasing.
    real(wp),dimension(:,:,:,:,:,:),intent(out)       :: bcoef !! `(nx,ny,nz,nq,nr,ns)` matrix of coefficients of the
                                                               !! b-spline interpolant.
    integer(ip),intent(out)                           :: iflag !! *  0 = successful execution.
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

    logical :: status_ok
    real(wp),dimension(:),allocatable :: temp !! work array of size `nx*ny*nz*nq*nr*ns`
    real(wp),dimension(:),allocatable :: work !! work array of size `max(2*kx*(nx+1),
                                              !! 2*ky*(ny+1),2*kz*(nz+1),2*kq*(nq+1),
                                              !! 2*kr*(nr+1),2*ks*(ns+1))`

    ! check validity of input
    call check_inputs(  iknot,&
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
        if (iknot == 0_ip) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
            call dbknot(q,nq,kq,tq)
            call dbknot(r,nr,kr,tr)
            call dbknot(s,ns,ks,ts)
        end if

        allocate(temp(nx*ny*nz*nq*nr*ns))
        allocate(work(max(2_ip*kx*(nx+1_ip),2_ip*ky*(ny+1_ip),&
                          2_ip*kz*(nz+1_ip),2_ip*kq*(nq+1_ip),&
                          2_ip*kr*(nr+1_ip),2_ip*ks*(ns+1_ip))))

        ! construct b-spline coefficients
                         call dbtpcf(x,nx,fcn,  nx,ny*nz*nq*nr*ns,tx,kx,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(y,ny,temp, ny,nx*nz*nq*nr*ns,ty,ky,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(z,nz,bcoef,nz,nx*ny*nq*nr*ns,tz,kz,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(q,nq,temp, nq,nx*ny*nz*nr*ns,tq,kq,bcoef,work,iflag)
        if (iflag==0_ip) call dbtpcf(r,nr,bcoef,nr,nx*ny*nz*nq*ns,tr,kr,temp, work,iflag)
        if (iflag==0_ip) call dbtpcf(s,ns,temp, ns,nx*ny*nz*nq*nr,ts,ks,bcoef,work,iflag)

        deallocate(temp)
        deallocate(work)

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
                           iloy,iloz,iloq,ilor,ilos,&
                           w5,w4,w3,w2,w1,w0,extrap)

    implicit none

    integer(ip),intent(in)                           :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: idz      !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: idq      !! \(q\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: idr      !! \(r\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: ids      !! \(s\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)                           :: nx       !! the number of interpolation points in \(x\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: ny       !! the number of interpolation points in \(y\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: nz       !! the number of interpolation points in \(z\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: nq       !! the number of interpolation points in \(q\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: nr       !! the number of interpolation points in \(r\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: ns       !! the number of interpolation points in \(s\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: kx       !! order of polynomial pieces in \(x\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: ky       !! order of polynomial pieces in \(y\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: kz       !! order of polynomial pieces in \(z\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: kq       !! order of polynomial pieces in \(q\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: kr       !! order of polynomial pieces in \(r\).
                                                                 !! (same as in last call to [[db6ink]])
    integer(ip),intent(in)                           :: ks       !! order of polynomial pieces in \(s\).
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
    integer(ip),intent(out)                          :: iflag    !! status flag:
                                                                 !!
                                                                 !! * \( = 0 \)   : no errors
                                                                 !! * \( \ne 0 \) : error
    integer(ip),intent(inout)                        :: inbvx    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: inbvy    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: inbvz    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: inbvq    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: inbvr    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: inbvs    !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: iloy     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: iloz     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: iloq     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: ilor     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    integer(ip),intent(inout)                        :: ilos     !! initialization parameter which must be set
                                                                 !! to 1 the first time this routine is called,
                                                                 !! and must not be changed by the user.
    real(wp),dimension(ky,kz,kq,kr,ks),intent(inout)              :: w5 !! work array
    real(wp),dimension(kz,kq,kr,ks),intent(inout)                 :: w4 !! work array
    real(wp),dimension(kq,kr,ks),intent(inout)                    :: w3 !! work array
    real(wp),dimension(kr,ks),intent(inout)                       :: w2 !! work array
    real(wp),dimension(ks),intent(inout)                          :: w1 !! work array
    real(wp),dimension(3_ip*max(kx,ky,kz,kq,kr,ks)),intent(inout) :: w0 !! work array
    logical,intent(in),optional                      :: extrap   !! if extrapolation is allowed
                                                                 !! (if not present, default is False)

    integer(ip) :: lefty,leftz,leftq,leftr,lefts,&
                   kcoly,kcolz,kcolq,kcolr,kcols,&
                   j,k,q,r,s

    f = 0.0_wp

    iflag = check_value(xval,tx,1_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(yval,ty,2_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(zval,tz,3_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(qval,tq,4_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(rval,tr,5_ip,extrap); if (iflag/=0_ip) return
    iflag = check_value(sval,ts,6_ip,extrap); if (iflag/=0_ip) return

    call dintrv(ty,ny+ky,yval,iloy,lefty,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tq,nq+kq,qval,iloq,leftq,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(tr,nr+kr,rval,ilor,leftr,iflag,extrap); if (iflag/=0_ip) return
    call dintrv(ts,ns+ks,sval,ilos,lefts,iflag,extrap); if (iflag/=0_ip) return

    iflag = 0_ip

    ! x -> y, z, q, r, s
    kcols = lefts - ks
    do s=1_ip,ks
        kcols = kcols + 1_ip
        kcolr = leftr - kr
        do r=1_ip,kr
            kcolr = kcolr + 1_ip
            kcolq = leftq - kq
            do q=1_ip,kq
                kcolq = kcolq + 1_ip
                kcolz = leftz - kz
                do k=1_ip,kz
                    kcolz = kcolz + 1_ip
                    kcoly = lefty - ky
                    do j=1_ip,ky
                        kcoly = kcoly + 1_ip
                        call dbvalu(tx,bcoef(:,kcoly,kcolz,kcolq,kcolr,kcols),&
                                             nx,kx,idx,xval,inbvx,w0,iflag,&
                                             w5(j,k,q,r,s),extrap)
                        if (iflag/=0_ip) return
                    end do
                end do
            end do
        end do
    end do

    ! y -> z, q, r, s
    kcoly = lefty - ky + 1_ip
    do s=1_ip,ks
        do r=1_ip,kr
            do q=1_ip,kq
                do k=1_ip,kz
                    call dbvalu(ty(kcoly:),w5(:,k,q,r,s),&
                                ky,ky,idy,yval,inbvy,w0,iflag,&
                                w4(k,q,r,s),extrap)
                    if (iflag/=0_ip) return
                end do
            end do
        end do
    end do

    ! z -> q, r, s
    kcolz = leftz - kz + 1_ip
    do s=1_ip,ks
        do r=1_ip,kr
            do q=1_ip,kq
                call dbvalu(tz(kcolz:),w4(:,q,r,s),&
                            kz,kz,idz,zval,inbvz,w0,iflag,&
                            w3(q,r,s),extrap)
                if (iflag/=0_ip) return
            end do
        end do
    end do

    ! q -> r, s
    kcolq = leftq - kq + 1_ip
    do s=1_ip,ks
        do r=1_ip,kr
            call dbvalu(tq(kcolq:),w3(:,r,s),&
                        kq,kq,idq,qval,inbvq,w0,iflag,&
                        w2(r,s),extrap)
            if (iflag/=0_ip) return
        end do
    end do

    ! r -> s
    kcolr = leftr - kr + 1_ip
    do s=1_ip,ks
        call dbvalu(tr(kcolr:),w2(:,s),&
                    kr,kr,idr,rval,inbvr,w0,iflag,&
                    w1(s),extrap)
        if (iflag/=0_ip) return
    end do

    ! s
    kcols = lefts - ks + 1_ip
    call dbvalu(ts(kcols:),w1,ks,ks,ids,sval,inbvs,w0,iflag,f,extrap)

    end subroutine db6val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Checks if the value is withing the range of the knot vectors.
!  This is called by the various `db*val` routines.

    pure function check_value(x,t,i,extrap) result(iflag)

    implicit none

    integer(ip)                      :: iflag   !! returns 0 if value is OK, otherwise returns `600+i`
    real(wp),intent(in)              :: x       !! the value to check
    integer(ip),intent(in)           :: i       !! 1=x, 2=y, 3=z, 4=q, 5=r, 6=s
    real(wp),dimension(:),intent(in) :: t       !! the knot vector
    logical,intent(in),optional      :: extrap  !! if extrapolation is allowed
                                                !! (if not present, default is False)

    logical :: allow_extrapolation  !! if extrapolation is allowed

    if (present(extrap)) then
        allow_extrapolation = extrap
    else
        allow_extrapolation = .false.
    end if

    if (allow_extrapolation) then
        ! in this case all values are OK
        iflag = 0_ip
    else
        if (x<t(1_ip) .or. x>t(size(t,kind=ip))) then
            iflag = 600_ip + i  ! value out of bounds (601, 602, etc.)
        else
            iflag = 0_ip
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

    pure subroutine check_inputs(iknot,&
                                 iflag,&
                                 nx,ny,nz,nq,nr,ns,&
                                 kx,ky,kz,kq,kr,ks,&
                                 x,y,z,q,r,s,&
                                 tx,ty,tz,tq,tr,ts,&
                                 f1,f2,f3,f4,f5,f6,&
                                 bcoef1,bcoef2,bcoef3,bcoef4,bcoef5,bcoef6,&
                                 alt,&
                                 status_ok)

    implicit none

    integer(ip),intent(in)                              :: iknot !! = 0 if the `INK` routine is computing the knots.
    integer(ip),intent(out)                             :: iflag
    integer(ip),intent(in),optional                     :: nx,ny,nz,nq,nr,ns
    integer(ip),intent(in),optional                     :: kx,ky,kz,kq,kr,ks
    real(wp),dimension(:),intent(in),optional           :: x,y,z,q,r,s
    real(wp),dimension(:),intent(in),optional           :: tx,ty,tz,tq,tr,ts
    real(wp),dimension(:),intent(in),optional           :: f1,bcoef1
    real(wp),dimension(:,:),intent(in),optional         :: f2,bcoef2
    real(wp),dimension(:,:,:),intent(in),optional       :: f3,bcoef3
    real(wp),dimension(:,:,:,:),intent(in),optional     :: f4,bcoef4
    real(wp),dimension(:,:,:,:,:),intent(in),optional   :: f5,bcoef5
    real(wp),dimension(:,:,:,:,:,:),intent(in),optional :: f6,bcoef6
    logical,intent(in),optional                         :: alt !! using the alt routine where 1st or
                                                               !! 2nd deriv is fixed at endpoints
                                                               !! [default is False]
    logical,intent(out)                                 :: status_ok

    logical :: error
    integer :: iex  !! extra points for the alt case (in `t` and `bcoef`)
                    !! [currently, only allowed for the 1D case & `k=4`]

    status_ok = .false.

    iex = 0_ip ! default
    if (present(alt)) then
        if (alt) iex = 2_ip  ! for "alt" mode
    end if

    if ((iknot < 0_ip) .or. (iknot > 1_ip)) then

        iflag = 2_ip ! iknot is out of range

    else

        call check('x',nx,kx,x,tx,[3_ip,  4_ip, 5_ip, 6_ip,706_ip,712_ip],iflag,error,iex); if (error) return
        call check('y',ny,ky,y,ty,[7_ip,  8_ip, 9_ip,10_ip,707_ip,713_ip],iflag,error,iex); if (error) return
        call check('z',nz,kz,z,tz,[11_ip,12_ip,13_ip,14_ip,708_ip,714_ip],iflag,error,iex); if (error) return
        call check('q',nq,kq,q,tq,[15_ip,16_ip,17_ip,18_ip,709_ip,715_ip],iflag,error,iex); if (error) return
        call check('r',nr,kr,r,tr,[19_ip,20_ip,21_ip,22_ip,710_ip,716_ip],iflag,error,iex); if (error) return
        call check('s',ns,ks,s,ts,[23_ip,24_ip,25_ip,26_ip,711_ip,717_ip],iflag,error,iex); if (error) return

        if (present(x) .and. present(f1) .and. present(bcoef1)) then
            if (size(x,kind=ip)/=size(f1,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef1,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(f2) .and. present(bcoef2)) then
            if (size(x,kind=ip)/=size(f2,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f2,2_ip,kind=ip))         then; iflag = 701_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef2,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)+iex/=size(bcoef2,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(f3) .and. &
            present(bcoef3)) then
            if (size(x,kind=ip)/=size(f3,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f3,2_ip,kind=ip))         then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f3,3_ip,kind=ip))         then; iflag = 702_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef3,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)+iex/=size(bcoef3,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)+iex/=size(bcoef3,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(f4) .and. present(bcoef4)) then
            if (size(x,kind=ip)/=size(f4,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f4,2_ip,kind=ip))         then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f4,3_ip,kind=ip))         then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f4,4_ip,kind=ip))         then; iflag = 703_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef4,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)+iex/=size(bcoef4,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)+iex/=size(bcoef4,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)+iex/=size(bcoef4,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(f5) .and. present(bcoef5)) then
            if (size(x,kind=ip)/=size(f5,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f5,2_ip,kind=ip))         then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f5,3_ip,kind=ip))         then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f5,4_ip,kind=ip))         then; iflag = 703_ip; return; end if
            if (size(r,kind=ip)/=size(f5,5_ip,kind=ip))         then; iflag = 704_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef5,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)+iex/=size(bcoef5,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)+iex/=size(bcoef5,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)+iex/=size(bcoef5,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
            if (size(r,kind=ip)+iex/=size(bcoef5,5_ip,kind=ip)) then; iflag = 804_ip; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(s) .and. present(f6) .and. present(bcoef6)) then
            if (size(x,kind=ip)/=size(f6,1_ip,kind=ip))         then; iflag = 700_ip; return; end if
            if (size(y,kind=ip)/=size(f6,2_ip,kind=ip))         then; iflag = 701_ip; return; end if
            if (size(z,kind=ip)/=size(f6,3_ip,kind=ip))         then; iflag = 702_ip; return; end if
            if (size(q,kind=ip)/=size(f6,4_ip,kind=ip))         then; iflag = 703_ip; return; end if
            if (size(r,kind=ip)/=size(f6,5_ip,kind=ip))         then; iflag = 704_ip; return; end if
            if (size(s,kind=ip)/=size(f6,6_ip,kind=ip))         then; iflag = 705_ip; return; end if
            if (size(x,kind=ip)+iex/=size(bcoef6,1_ip,kind=ip)) then; iflag = 800_ip; return; end if
            if (size(y,kind=ip)+iex/=size(bcoef6,2_ip,kind=ip)) then; iflag = 801_ip; return; end if
            if (size(z,kind=ip)+iex/=size(bcoef6,3_ip,kind=ip)) then; iflag = 802_ip; return; end if
            if (size(q,kind=ip)+iex/=size(bcoef6,4_ip,kind=ip)) then; iflag = 803_ip; return; end if
            if (size(r,kind=ip)+iex/=size(bcoef6,5_ip,kind=ip)) then; iflag = 804_ip; return; end if
            if (size(s,kind=ip)+iex/=size(bcoef6,6_ip,kind=ip)) then; iflag = 805_ip; return; end if

        end if

        status_ok = .true.
        iflag = 0_ip

    end if

    contains

        pure subroutine check(s,n,k,x,t,ierrs,iflag,error,ik)  !! check `t`,`x`,`n`,`k` for validity

        implicit none

        character(len=1),intent(in)               :: s     !! coordinate string: 'x','y','z','q','r','s'
        integer(ip),intent(in),optional           :: n     !! size of `x`
        integer(ip),intent(in),optional           :: k     !! order
        real(wp),dimension(:),intent(in),optional :: x     !! abcissae vector
        real(wp),dimension(:),intent(in),optional :: t     !! knot vector `size(n+k)`
        integer(ip),dimension(:),intent(in)       :: ierrs !! int error codes for `n`,`k`,`x`,`t`,
                                                           !! `size(x)`,`size(t)` checks
        integer(ip),intent(out)                   :: iflag !! status return code
        logical,intent(out)                       :: error !! true if there was an error
        integer,intent(in)                        :: ik    !! add this value to k

        integer(ip),dimension(2) :: itmp !! temp integer array

        if (present(n) .and. present(k) .and. present(x) .and. present(t)) then
            itmp = [ierrs(1_ip),ierrs(5)]
            call check_n('n'//s,n,x,itmp,iflag,error);     if (error) return
            call check_k('k'//s,k+ik,n,ierrs(2),iflag,error); if (error) return
            call check_x(s,n,x,ierrs(3),iflag,error);      if (error) return
            if (iknot /= 0_ip) then
                itmp = [ierrs(4),ierrs(6)]
                call check_t('t'//s,n,k+ik,t,itmp,iflag,error); if (error) return
            end if
        end if

        end subroutine check

        pure subroutine check_n(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)         :: s
        integer(ip),intent(in)              :: n
        real(wp),dimension(:),intent(in)    :: x     !! abcissae vector
        integer(ip),dimension(2),intent(in) :: ierr  !! [n<3 check, size(x)==n check]
        integer(ip),intent(out)             :: iflag !! status return code
        logical,intent(out)                 :: error

        if (n < 3_ip) then
            iflag = ierr(1_ip)
            error = .true.
        else
            if (size(x)/=n) then
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
        integer(ip),intent(in)      :: k
        integer(ip),intent(in)      :: n
        integer(ip),intent(in)      :: ierr
        integer(ip),intent(out)     :: iflag !! status return code
        logical,intent(out)         :: error

        if ((k < 2_ip) .or. (k >= n)) then
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if

        end subroutine check_k

        pure subroutine check_x(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer(ip),intent(in)            :: n
        real(wp),dimension(:),intent(in)  :: x
        integer(ip),intent(in)            :: ierr
        integer(ip),intent(out)           :: iflag !! status return code
        logical,intent(out)               :: error

        integer(ip) :: i

        error = .true.
        do i=2_ip,n
            if (x(i) <= x(i-1_ip)) then
                iflag = ierr
                return
            end if
        end do
        error = .false.

        end subroutine check_x

        pure subroutine check_t(s,n,k,t,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)         :: s
        integer(ip),intent(in)              :: n
        integer(ip),intent(in)              :: k
        real(wp),dimension(:),intent(in)    :: t
        integer(ip),dimension(2),intent(in) :: ierr  !! [non-decreasing check, size check]
        integer(ip),intent(out)             :: iflag !! status return code
        logical,intent(out)                 :: error

        integer(ip) :: i

        error = .true.

        if (size(t)/=(n+k)) then
            iflag = ierr(2)
            return
        end if

        if (iex==0_ip) then ! don't do this for "alt" mode since they haven't been computed yet
            do i=2_ip,n + k
                if (t(i) < t(i-1_ip))  then
                    iflag = ierr(1_ip)
                    return
                end if
            end do
        end if

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

    integer(ip),intent(in)             :: n  !! dimension of `x`
    integer(ip),intent(in)             :: k
    real(wp),dimension(:),intent(in)   :: x
    real(wp),dimension(:),intent(out)  :: t

    integer(ip) :: i, j, ipj, npj, ip1, jstrt
    real(wp) :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1_wp*( x(n)-x(n-1_ip) )
    do j=1_ip,k
        t(j)   = x(1_ip)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2_ip) == 1_ip)  then

        !case of odd k --  knots between data points

        i = (k-1_ip)/2_ip - k
        ip1 = i + 1_ip
        jstrt = k + 1_ip
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5_wp*( x(ipj) + x(ipj+1_ip) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2_ip) - k
        jstrt = k+1_ip
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

    integer(ip),intent(in)                :: n  !! dimension of `x`
    integer(ip),intent(in)                :: nf
    integer(ip),intent(in)                :: ldf
    integer(ip),intent(in)                :: k
    real(wp),dimension(:),intent(in)      :: x
    real(wp),dimension(ldf,nf),intent(in) :: fcn
    real(wp),dimension(:),intent(in)      :: t
    real(wp),dimension(nf,n),intent(out)  :: bcoef
    real(wp),dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer(ip),intent(out)               :: iflag  !! status flag:
                                                    !!
                                                    !! * 0: no errors
                                                    !! * 301: n should be >0

    integer(ip) :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0_ip)  then

        ! partition work array
        m1 = k - 1_ip
        m2 = m1 + k
        iq = 1_ip + n
        iw = iq + m2*n+1_ip

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0_ip) then
            do i=1_ip,n
                bcoef(1_ip,i) = work(i)
            end do

            !  all remaining data sets by back-substitution

            if (nf == 1_ip)  return
            do j=2_ip,nf
                do i=1_ip,n
                    work(i) = fcn(i,j)
                end do
                call dbnslv(work(iq),m2,n,m1,m1,work)
                do i=1_ip,n
                    bcoef(j,i) = work(i)
                end do
            end do
        end if

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301_ip
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
!### Error conditions
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

    integer(ip),intent(in)            :: n      !! number of data points, n >= k
    real(wp),dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(wp),dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real(wp),dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer(ip),intent(in)            :: k      !! order of the spline, k >= 1
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
    integer(ip),intent(out)           :: iflag  !! *   0: no errors.
                                                !! * 100: k does not satisfy k>=1.
                                                !! * 101: n does not satisfy n>=k.
                                                !! * 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! * 103: some abscissa was not in the support of the
                                                !! corresponding basis function and the system is singular.
                                                !! * 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer(ip) :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(wp) :: xi
    logical :: found

    if (k<1_ip) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100_ip
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101_ip
        return
    end if

    jj = n - 1_ip
    if (jj/=0_ip) then
        do i=1_ip,jj
            if (x(i)>=x(i+1_ip)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102_ip
                return
            end if
        end do
    end if

    np1 = n + 1_ip
    km1 = k - 1_ip
    kpkm2 = 2_ip*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1_ip,lenq
        q(i) = 0.0_wp
    end do

    ! loop over i to construct the n interpolation equations
    do i=1_ip,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! find left in the closed interval (i,i+k-1_ip) such that
        !         t(left) <= x(i) < t(left+1_ip)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103_ip
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1_ip))
            if (found) exit
            left = left + 1_ip
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1_ip
            if (xi>t(left+1_ip)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103_ip
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1_ip, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0_ip) return

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
        jj = i - left + 1_ip + (left-k)*(k+km1)
        do j=1_ip,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! solve  a*bcoef = y  by backsubstitution
        do i=1_ip,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0_ip
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104_ip
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

    integer(ip),intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer(ip),intent(in) :: nrow    !! matrix order
    integer(ip),intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer(ip),intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer(ip),intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real(wp),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer(ip) :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real(wp) :: factor, pivot

    iflag = 1_ip
    middle = nbandu + 1_ip   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1_ip

    if (nrowm1 < 0_ip) then
        iflag = 2_ip
        return
    else if (nrowm1 == 0_ip) then
        if (w(middle,nrow)==0.0_wp) iflag = 2_ip
        return
    end if

    if (nbandl<=0_ip) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1_ip,nrowm1
            if (w(middle,i)==0.0_wp) then
                iflag = 2_ip
                return
            end if
        end do
        if (w(middle,nrow)==0.0_wp) iflag = 2_ip
        return
    end if

    if (nbandu<=0_ip) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1_ip,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0_wp) then
                iflag = 2_ip
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1_ip,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1_ip,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0_wp) then
            iflag = 2_ip
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1_ip,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1_ip,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1_ip,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0_wp) iflag = 2_ip

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

    integer(ip),intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    integer(ip),intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order `nrow`
                                      !! as constructed in [[dbnfac]].
    real(wp),dimension(nroww,nrow),intent(in) :: w !! describes the lu-factorization of a banded matrix a of
                                                   !! order `nrow` as constructed in [[dbnfac]].
    real(wp),dimension(nrow),intent(inout) :: b  !! * **in**: right side of the system to be solved
                                                 !! * **out**: the solution x, of order nrow

    integer(ip) :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1_ip
    if (nrow/=1_ip) then

        nrowm1 = nrow - 1_ip
        if (nbandl/=0_ip) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l)
            !                       from right side (below i-th row).
            do i=1_ip,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1_ip,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0_ip) then
            ! a is lower triangular.
            do i=1_ip,nrow
                b(i) = b(i)/w(1_ip,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1_ip)
            do j=1_ip,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1_ip
            if (i<=1_ip) exit
        end do

    end if

    b(1_ip) = b(1_ip)/w(middle,1_ip)

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

    real(wp),dimension(*),intent(in)    :: t        !! knot vector of length `n+k`, where
                                                    !! `n` = number of b-spline basis functions
                                                    !! `n` = sum of knot multiplicities-`k`
                                                    !! dimension `t(ileft+jhigh)`
    integer(ip),intent(in)              :: jhigh    !! order of b-spline, `1 <= jhigh <= k`
    integer(ip),intent(in)              :: k        !! highest possible order
    integer(ip),intent(in)              :: index    !! index = 1 gives basis functions of order `jhigh`
                                                    !!       = 2 denotes previous entry with `work`, `iwork`
                                                    !!         values saved for subsequent calls to
                                                    !!         dbspvn.
    real(wp),intent(in)                 :: x        !! argument of basis functions, `t(k) <= x <= t(n+1)`
    integer(ip),intent(in)              :: ileft    !! largest integer such that `t(ileft) <= x < t(ileft+1)`
    real(wp),dimension(k),intent(out)   :: vnikx    !! vector of length `k` for spline values.
    real(wp),dimension(*),intent(inout) :: work     !! a work vector of length `2*k`
    integer(ip),intent(inout)           :: iwork    !! a work parameter.  both `work` and `iwork` contain
                                                    !! information necessary to continue for `index = 2`.
                                                    !! when `index = 1` exclusively, these are scratch
                                                    !! variables and can be used for other purposes.
    integer(ip),intent(out)             :: iflag    !! *   0: no errors
                                                    !! * 201: `k` does not satisfy `k>=1`
                                                    !! * 202: `jhigh` does not satisfy `1<=jhigh<=k`
                                                    !! * 203: `index` is not 1 or 2
                                                    !! * 204: `x` does not satisfy `t(ileft)<=x<=t(ileft+1)`

    integer(ip) :: imjp1, ipj, jp1, jp1ml, l
    real(wp) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1_ip) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201_ip
        return
    end if
    if (jhigh>k .or. jhigh<1_ip) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202_ip
        return
    end if
    if (index<1_ip .or. index>2_ip) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203_ip
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1_ip)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204_ip
        return
    end if

    iflag = 0_ip

    if (index==1_ip) then
        iwork = 1_ip
        vnikx(1_ip) = 1.0_wp
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1_ip
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0_wp
        jp1 = iwork + 1_ip
        do l=1_ip,iwork
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
    integer(ip),intent(in)           :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-`k`)
    real(wp),dimension(:),intent(in) :: t       !! knot vector of length `n+k`
    real(wp),dimension(n),intent(in) :: a       !! b-spline coefficient vector of length `n`
    integer(ip),intent(in)           :: k       !! order of the b-spline, `k >= 1`
    integer(ip),intent(in)           :: ideriv  !! order of the derivative, `0 <= ideriv <= k-1`.
                                                !! `ideriv = 0` returns the b-spline value
    real(wp),intent(in)              :: x       !! argument, `t(k) <= x <= t(n+1)`
    integer(ip),intent(inout)        :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time [[dbvalu]] is called.
                                                !! `inbv` contains information for efficient processing
                                                !! after the initial call and `inbv` must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct `inbv` parameters.
    real(wp),dimension(:),intent(inout) :: work !! work vector of length at least `3*k`
    integer(ip),intent(out)          :: iflag   !! status flag:
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

    integer(ip) :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(wp) :: fkmj
    real(wp) :: xt
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    val = 0.0_wp

    if (k<1_ip) then
        iflag = 401_ip  ! dbvalu - k does not satisfy k>=1
        return
    end if

    if (n<k) then
        iflag = 402_ip  ! dbvalu - n does not satisfy n>=k
        return
    end if

    if (ideriv<0_ip .or. ideriv>=k) then
        iflag = 403_ip  ! dbvalu - ideriv does not satisfy 0<=ideriv<k
        return
    end if

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! make a temp copy of x (for computing the
    ! interval) in case extrapolation is allowed
    if (extrapolation_allowed) then
        if (x<t(k)) then
            xt = t(k)
        else if (x>t(n+1_ip)) then
            xt = t(n+1_ip)
        else
            xt = x
        end if
    else
        xt = x
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1_ip
    call dintrv(t, n+1, xt, inbv, i, mflag)
    if (xt<t(k)) then
        iflag = 404_ip  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0_ip) then

        if (xt>t(i)) then
            iflag = 405_ip  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (i==k) then
                iflag = 406_ip  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            i = i - 1_ip
            if (xt/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1_ip,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0_ip) then
        do j=1_ip,ideriv
            kmj = k - j
            fkmj = real(kmj,wp)
            do jj=1_ip,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1_ip)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1_ip
        kpk = k + k
        j1 = k + 1_ip
        j2 = kpk + 1_ip
        do j=1_ip,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1_ip
            j2 = j2 + 1_ip
        end do
        iderp1 = ideriv + 1_ip
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1_ip,kmj
                work(jj) = (work(jj+1_ip)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    iflag = 0_ip
    val = work(1_ip)

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
!         if xt(lxt) <= x           then ileft=lxt, mflag=-2
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

    integer(ip),intent(in)             :: lxt    !! length of the `xt` vector
    real(wp),dimension(:),intent(in)   :: xt     !! a knot or break point vector of length `lxt`
    real(wp),intent(in)                :: xx     !! argument
    integer(ip),intent(inout)          :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct `ilo` parameters.
    integer(ip),intent(out)            :: ileft  !! largest integer satisfying `xt(ileft)` \( \le \) `x`
    integer(ip),intent(out)            :: mflag  !! signals when `x` lies out of bounds
    logical,intent(in),optional        :: extrap !! if extrapolation is allowed
                                                 !! (if not present, default is False)

    integer(ip) :: ihi, istep, middle
    real(wp) :: x

    x = get_temp_x_for_extrap(xx,xt(1_ip),xt(lxt),extrap)

    ihi = ilo + 1_ip
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = -2_ip
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1_ip
            ileft = 1_ip
            return
        end if
        ilo = lxt - 1_ip
        ihi = lxt
    end if

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1_ip
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = -2_ip
                    ileft = lxt
                    return
                end if
                ihi = lxt
            else if ( x>=xt(ihi) ) then
                istep = istep*2_ip
                cycle
            end if
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0_ip
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1_ip
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1_ip ) then
                ilo = 1_ip
                if ( x<xt(1_ip) ) then
                    mflag = -1_ip
                    ileft = 1_ip
                    return
                end if
            else if ( x<xt(ilo) ) then
                istep = istep*2_ip
                cycle
            end if
            exit
        end do

    end if

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2_ip
        if ( middle==ilo ) then
            mflag = 0_ip
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
!  DBINT4 computes the B representation (`t`,`bcoef`,`n`,`k`) of a
!  cubic spline (`k=4`) which interpolates data (`x(i)`,`y(i)`),`i=1,ndata`.
!
!  Parameters `ibcl`, `ibcr`, `fbcl`, `fbcr` allow the specification of the spline
!  first or second derivative at both `x(1)` and `x(ndata)`.  When this data is not specified
!  by the problem, it is common practice to use a natural spline by setting second
!  derivatives at `x(1)` and `x(ndata)` to zero (`ibcl=ibcr=2`,`fbcl=fbcr=0.0`).
!
!  The spline is defined on `t(4) <= x <= t(n+1)` with (ordered) interior knots at
!  `x(i)` values where n=ndata+2.  The knots `t(1)`,`t(2)`,`t(3)` lie to the left of
!  `t(4)=x(1)` and the knots `t(n+2)`, `t(n+3)`, `t(n+4)` lie to the right of `t(n+1)=x(ndata)`
!  in increasing order.
!
!  * If no extrapolation outside (`x(1)`,`x(ndata)`) is anticipated, the
!    knots `t(1)=t(2)=t(3)=t(4)=x(1)` and `t(n+2)=t(n+3)=t(n+4)=t(n+1)=x(ndata)`
!    can be specified by `kntopt=1`.
!  * `kntopt=2` selects a knot placement for `t(1)`, `t(2)`, `t(3)` to make the
!    first 7 knots symmetric about `t(4)=x(1)` and similarly for
!    `t(n+2)`, `t(n+3)`, `t(n+4)` about `t(n+1)=x(ndata)`.
!  * `kntopt=3` allows the user to make his own selection, in increasing order,
!    for `t(1)`, `t(2)`, `t(3)` to the left of `x(1)` and `t(n+2)`, `t(n+3)`, `t(n+4)` to
!    the right of x(ndata).
!
!  In any case, the interpolation on `t(4) <= x <= t(n+1)`
!  by using function [[dbvalu]] is unique for given boundary
!  conditions.
!
!### Error conditions
!  * improper input
!  * singular system of equations
!
!### See also
!  * [[dbintk]]
!
!### History
!  * Written by D. E. Amos (SNLA), August, 1979.
!  * date written 800901
!  * revision date 820801
!  * 000330  Modified array declarations.  (JEC)
!  * Jacob Williams, 8/30/2018 : refactored to modern Fortran.

    pure subroutine dbint4(x,y,ndata,ibcl,ibcr,fbcl,fbcr,kntopt,tleft,tright,t,bcoef,n,k,w,iflag)

    implicit none

    real(wp),dimension(:),intent(in)   :: x       !! x vector of abscissae of length `ndata`, distinct
                                                  !! and in increasing order
    real(wp),dimension(:),intent(in)   :: y       !! y vector of ordinates of length ndata
    integer(ip),intent(in)             :: ndata   !! number of data points, `ndata >= 2`
    integer(ip),intent(in)             :: ibcl    !! selection parameter for left boundary condition:
                                                  !!
                                                  !! * `ibcl = 1` constrain the first derivative at `x(1)` to `fbcl`
                                                  !! * `ibcl = 2` constrain the second derivative at `x(1)` to `fbcl`
    integer(ip),intent(in)             :: ibcr    !! selection parameter for right boundary condition:
                                                  !!
                                                  !! * `ibcr = 1` constrain first derivative at `x(ndata)` to `fbcr`
                                                  !! * `ibcr = 2` constrain second derivative at `x(ndata)` to `fbcr`
    real(wp),intent(in)                :: fbcl    !! left boundary values governed by `ibcl`
    real(wp),intent(in)                :: fbcr    !! right boundary values governed by `ibcr`
    integer(ip),intent(in)             :: kntopt  !! knot selection parameter:
                                                  !!
                                                  !! * `kntopt = 1` sets knot multiplicity at `t(4)` and
                                                  !!   `t(n+1)` to 4
                                                  !! * `kntopt = 2` sets a symmetric placement of knots
                                                  !!   about `t(4)` and `t(n+1)`
                                                  !! * `kntopt = 3` sets `t(i)=tleft(i)` and
                                                  !!   `t(n+1+i)=tright(i)`,`i=1,3`
    real(wp),dimension(3),intent(in)  :: tleft    !! when `kntopt = 3`: `t(1:3)` in increasing
                                                  !! order to be supplied by the user.
    real(wp),dimension(3),intent(in)  :: tright   !! when `kntopt = 3`: `t(n+2:n+4)` in increasing
                                                  !! order to be supplied by the user.
    real(wp),dimension(:),intent(out)  :: t       !! knot array of length `n+4`
    real(wp),dimension(:),intent(out)  :: bcoef   !! b spline coefficient array of length `n`
    integer(ip),intent(out)            :: n       !! number of coefficients, `n=ndata+2`
    integer(ip),intent(out)            :: k       !! order of spline, `k=4`
    real(wp),dimension(5,ndata+2),intent(inout) :: w  !! work array
    integer(ip),intent(out)            :: iflag   !! status flag:
                                                  !!
                                                  !! * 0: no errors
                                                  !! * 2001: `ndata` is less than 2
                                                  !! * 2002: `x` values are not distinct or not ordered
                                                  !! * 2003: `ibcl` is not 1 or 2
                                                  !! * 2004: `ibcr` is not 1 or 2
                                                  !! * 2005: `kntopt` is not 1, 2, or 3
                                                  !! * 2006: knot input through `tleft`, `tright` is
                                                  !!   not ordered properly
                                                  !! * 2007: the system of equations is singular

    integer(ip)  :: i, ilb, ileft, it, iub, iw, iwp, j, jw, ndm, np, nwrow
    real(wp) :: txn, tx1, xl
    real(wp),dimension(4,4) :: vnikx
    real(wp),dimension(15) :: work  !! work array for [[dbspvd]] -- length `(k+1)*(k+2)/2`

    real(wp),parameter :: wdtol = epsilon(1.0_wp) !! d1mach(4)
    real(wp),parameter :: tol = sqrt(wdtol)

    if (ndata<2_ip) then
        iflag = 2001_ip ! ndata is less than 2
        return
    end if

    ndm = ndata - 1_ip
    do i=1_ip,ndm
        if (x(i)>=x(i+1_ip)) then
            iflag = 2002_ip ! x values are not distinct or not ordered
            return
        end if
    end do

    if (ibcl<1_ip .or. ibcl>2_ip) then
        iflag = 2003_ip ! ibcl is not 1 or 2
        return
    end if

    if (ibcr<1_ip .or. ibcr>2_ip) then
        iflag = 2004_ip ! ibcr is not 1 or 2
        return
    end if

    if (kntopt<1_ip .or. kntopt>3_ip) then
        iflag = 2005_ip ! kntopt is not 1, 2, or 3
        return
    end if

    iflag = 0_ip
    k = 4_ip
    n = ndata + 2_ip
    np = n + 1_ip
    do i=1_ip,ndata
        t(i+3) = x(i)
    end do

    select case (kntopt)
    case(1_ip)
        ! set up knot array with multiplicity 4 at x(1) and x(ndata)
        do i=1,3_ip
            t(4-i) = x(1)
            t(np+i) = x(ndata)
        end do
    case(2_ip)
        !set up knot array with symmetric placement about end points
        if (ndata>3) then
            tx1 = x(1) + x(1)
            txn = x(ndata) + x(ndata)
            do i=1,3
                t(4-i) = tx1 - x(i+1)
                t(np+i) = txn - x(ndata-i)
            end do
        else
            xl = (x(ndata)-x(1))/3.0_wp
            do i=1,3
                t(4-i) = t(5-i) - xl
                t(np+i) = t(np+i-1) + xl
            end do
        end if
    case(3)
        ! set up knot array less than x(1) and greater than x(ndata) to be
        ! supplied by user in tleft & tright when kntopt=3
        t(1:3)             = tleft
        t(ndata+4:ndata+6) = tright
        do i=1,3
            if ((t(4-i)>t(5-i)) .or. (t(np+i)<t(np+i-1))) then
                iflag = 2006_ip ! knot input through tleft, tright is not ordered properly
                return
            end if
        end do
    end select

    w = 0.0_wp

    ! set up left interpolation point and left boundary condition for
    ! right limits
    it = ibcl + 1
    call dbspvd(t, k, it, x(1), k, 4_ip, vnikx, work, iflag)
    if (iflag/=0_ip) return ! error check
    iw = 0_ip
    if (abs(vnikx(3,1))<tol) iw = 1_ip
    do j=1,3
        w(j+1,4-j) = vnikx(4-j,it)
        w(j,4-j) = vnikx(4-j,1)
    end do
    bcoef(1) = y(1)
    bcoef(2) = fbcl
    ! set up interpolation equations for points i=2 to i=ndata-1
    ileft = 4_ip
    if (ndm>=2) then
        do i=2,ndm
            ileft = ileft + 1_ip
            call dbspvd(t, k, 1_ip, x(i), ileft, 4_ip, vnikx, work, iflag)
            if (iflag/=0_ip) return ! error check
            do j=1,3
                w(j+1,3+i-j) = vnikx(4-j,1)
            end do
            bcoef(i+1) = y(i)
        end do
    end if

    ! set up right interpolation point and right boundary condition for
    ! left limits(ileft is associated with t(n)=x(ndata-1))
    it = ibcr + 1_ip
    call dbspvd(t, k, it, x(ndata), ileft, 4_ip, vnikx, work, iflag)
    if (iflag/=0_ip) return ! error check
    jw = 0_ip
    if (abs(vnikx(2,1))<tol) jw = 1_ip
    do j=1,3
        w(j+1,3+ndata-j) = vnikx(5-j,it)
        w(j+2,3+ndata-j) = vnikx(5-j,1)
    end do
    bcoef(n-1) = fbcr
    bcoef(n) = y(ndata)
    ! solve system of equations
    ilb = 2_ip - jw
    iub = 2_ip - iw
    nwrow = 5_ip
    iwp = iw + 1_ip
    call dbnfac(w(iwp,1), nwrow, n, ilb, iub, iflag)
    if (iflag==2_ip) then
        iflag = 2007_ip  ! the system of equations is singular
    else
        iflag = 0_ip  ! success
        call dbnslv(w(iwp,1), nwrow, n, ilb, iub, bcoef)
    end if

    end subroutine dbint4
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSPVD calculates the value and all derivatives of order
!  less than `nderiv` of all basis functions which do not
!  (possibly) vanish at `x`.  `ileft` is input such that
!  `t(ileft) <= x < t(ileft+1)`.  A call to [[dintrv]](`t`,`n+1`,`x`,
!  `ilo`,`ileft`,`mflag`) will produce the proper `ileft`.  The output of
!  dbspvd is a matrix `vnikx(i,j)` of dimension at least `(k,nderiv)`
!  whose columns contain the `k` nonzero basis functions and
!  their `nderiv-1` right derivatives at `x`, `i=1,k, j=1,nderiv`.
!  These basis functions have indices `ileft-k+i`, `i=1,k,
!  k <= ileft <= n`.  The nonzero part of the `i`-th basis
!  function lies in `(t(i),t(i+k)), i=1,n)`.
!
!  If `x=t(ileft+1)` then `vnikx` contains left limiting values
!  (left derivatives) at `t(ileft+1)`.  In particular, `ileft = n`
!  produces left limiting values at the right end point
!  `x=t(n+1)`.  To obtain left limiting values at `t(i)`, `i=k+1,n+1`,
!  set `x` = next lower distinct knot, call [[dintrv]] to get `ileft`,
!  set `x=t(i)`, and then call dbspvd.
!
!### History
!  * Written by Carl de Boor and modified by D. E. Amos
!  * date written 800901
!  * revision date 820801
!  * 000330  Modified array declarations.  (JEC)
!  * Jacob Williams, 8/30/2018 : refactored to modern Fortran.
!
!@note `DBSPVD` is the `BSPLVD` routine of the reference.

    pure subroutine dbspvd(t,k,nderiv,x,ileft,ldvnik,vnikx,work,iflag)

    implicit none

    real(wp),dimension(:),intent(in)              :: t       !! knot vector of length `n+k`, where
                                                             !! `n` = number of b-spline basis functions
                                                             !! `n` = sum of knot multiplicities-k
    integer(ip),intent(in)                        :: k       !! order of the b-spline, `k >= 1`
    integer(ip),intent(in)                        :: nderiv  !! number of derivatives = `nderiv-1`,
                                                             !! `1 <= nderiv <= k`
    real(wp),intent(in)                           :: x       !! argument of basis functions,
                                                             !! `t(k) <= x <= t(n+1)`
    integer(ip),intent(in)                        :: ileft   !! largest integer such that
                                                             !! `t(ileft) <= x < t(ileft+1)`
    integer(ip),intent(in)                        :: ldvnik  !! leading dimension of matrix `vnikx`
    real(wp),dimension(ldvnik,nderiv),intent(out) :: vnikx   !! matrix of dimension at least `(k,nderiv)`
                                                             !! containing the nonzero basis functions
                                                             !! at `x` and their derivatives columnwise.
    real(wp),dimension(*),intent(out)             :: work    !! a work vector of length `(k+1)*(k+2)/2`
    integer(ip),intent(out)                       :: iflag   !! status flag:
                                                             !!
                                                             !! * 0: no errors
                                                             !! * 3001: `k` does not satisfy `k>=1`
                                                             !! * 3002: `nderiv` does not satisfy `1<=nderiv<=k`
                                                             !! * 3003: `ldvnik` does not satisfy `ldvnik>=k`

    integer(ip) :: i,ideriv,ipkmd,j,jj,jlow,jm,jp1mid,kmd,kp1,l,ldummy,m,mhigh,iwork
    real(wp) :: factor, fkmd, v

    ! dimension t(ileft+k), work((k+1)*(k+2)/2)
    ! a(i,j) = work(i+j*(j+1)/2),  i=1,j+1  j=1,k-1
    ! a(i,k) = work(i+k*(k-1)/2)  i=1.k
    ! work(1) and work((k+1)*(k+2)/2) are not used.

    if (k<1) then
        iflag = 3001_ip ! k does not satisfy k>=1
        return
    end if

    if (nderiv<1 .or. nderiv>k) then
        iflag = 3002_ip ! nderiv does not satisfy 1<=nderiv<=k
        return
    end if

    if (ldvnik<k) then
        iflag = 3003_ip ! ldvnik does not satisfy ldvnik>=k
        return
    end if

    iflag = 0_ip

    ideriv = nderiv
    kp1 = k + 1
    jj = kp1 - ideriv
    call dbspvn(t, jj, k, 1_ip, x, ileft, vnikx, work, iwork, iflag)
    if (iflag/=0 .or. ideriv==1) return
    mhigh = ideriv
    do m=2,mhigh
        jp1mid = 1
        do j=ideriv,k
            vnikx(j,ideriv) = vnikx(jp1mid,1)
            jp1mid = jp1mid + 1
        end do
        ideriv = ideriv - 1
        jj = kp1 - ideriv
        call dbspvn(t, jj, k, 2_ip, x, ileft, vnikx, work, iwork, iflag)
        if (iflag/=0) return
    end do

    jm = kp1*(kp1+1)/2
    do l = 1,jm
        work(l) = 0.0_wp
    end do
    ! a(i,i) = work(i*(i+3)/2) = 1.0       i = 1,k
    l = 2
    j = 0
    do i = 1,k
        j = j + l
        work(j) = 1.0_wp
        l = l + 1
    end do
    kmd = k
    do m=2,mhigh
        kmd = kmd - 1
        fkmd = real(kmd,wp)
        i = ileft
        j = k
        jj = j*(j+1)/2
        jm = jj - j
        do ldummy=1,kmd
            ipkmd = i + kmd
            factor = fkmd/(t(ipkmd)-t(i))
            do l=1,j
                work(l+jj) = (work(l+jj)-work(l+jm))*factor
            end do
            i = i - 1
            j = j - 1
            jj = jm
            jm = jm - j
        end do

        do i=1,k
            v = 0.0_wp
            jlow = max(i,m)
            jj = jlow*(jlow+1)/2
            do j=jlow,k
                v = work(i+jj)*vnikx(j,m) + v
                jj = jj + j + 1
            end do
            vnikx(i,m) = v
        end do
    end do

    end subroutine dbspvd
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSQAD computes the integral on `(x1,x2)` of a `k`-th order
!  b-spline using the b-representation `(t,bcoef,n,k)`.  orders
!  `k` as high as 20 are permitted by applying a 2, 6, or 10
!  point gauss formula on subintervals of `(x1,x2)` which are
!  formed by included (distinct) knots.
!
!  If orders `k` greater than 20 are needed, use [[dbfqad]] with
!  `f(x) = 1`.
!
!### Note
!  * The maximum number of significant digits obtainable in
!    DBSQAD is the smaller of ~300 and the number of digits
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
!    Added higher precision coefficients.
!
!@note Extrapolation is not enabled for this routine.

    pure subroutine dbsqad(t,bcoef,n,k,x1,x2,bquad,work,iflag)

    implicit none

    real(wp),dimension(:),intent(in)    :: t       !! knot array of length `n+k`
    real(wp),dimension(:),intent(in)    :: bcoef   !! b-spline coefficient array of length `n`
    integer(ip),intent(in)              :: n       !! length of coefficient array
    integer(ip),intent(in)              :: k       !! order of b-spline, `1 <= k <= 20`
    real(wp),intent(in)                 :: x1      !! end point of quadrature interval
                                                   !! in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: x2      !! end point of quadrature interval
                                                   !! in `t(k) <= x <= t(n+1)`
    real(wp),intent(out)                :: bquad   !! integral of the b-spline over (`x1`,`x2`)
    real(wp),dimension(:),intent(inout) :: work    !! work vector of length `3*k`
    integer(ip),intent(out)             :: iflag   !! status flag:
                                                   !!
                                                   !! * 0: no errors
                                                   !! * 901: `k` does not satisfy `1<=k<=20`
                                                   !! * 902: `n` does not satisfy `n>=k`
                                                   !! * 903: `x1` or `x2` or both do
                                                   !!   not satisfy `t(k)<=x<=t(n+1)`

    integer(ip) :: i,il1,il2,ilo,inbv,jf,left,m,mf,mflag,npk,np1
    real(wp) :: a,aa,b,bb,bma,bpa,c1,gx,q,ta,tb,y1,y2
    real(wp),dimension(5) :: s  !! sum

    real(wp),dimension(9),parameter :: gpts = [ &
        &0.577350269189625764509148780501957455647601751270126876018602326483977&
        &67230293334569371539558574952522520871380513556767665664836499965082627&
        &05518373647912161760310773007685273559916067003615583077550051041144223&
        &01107628883557418222973945990409015710553455953862673016662179126619796&
        &4892168_wp,&
        &0.238619186083196908630501721680711935418610630140021350181395164574274&
        &93427563984224922442725734913160907222309701068720295545303507720513526&
        &28872175189982985139866216812636229030578298770859440976999298617585739&
        &46921613621659222233462641640013936777894532787145324672151888999339900&
        &0945406150514997832_wp,&
        &0.661209386466264513661399595019905347006448564395170070814526705852183&
        &49660714310094428640374646145642988837163927514667955734677222538043817&
        &23198010093367423918538864300079016299442625145884902455718821970386303&
        &22362011735232135702218793618906974301231555871064213101639896769013566&
        &1651261150514997832_wp,&
        &0.932469514203152027812301554493994609134765737712289824872549616526613&
        &50084420019627628873992192598504786367972657283410658797137951163840419&
        &21786180750210169211578452038930846310372961174632524612619760497437974&
        &07422632089671621172178385230505104744277222209386367655366917903888025&
        &2326771150514997832_wp,&
        &0.148874338981631210884826001129719984617564859420691695707989253515903&
        &61735566852137117762979946369123003116080525533882610289018186437654023&
        &16761969968090913050737827720371059070942475859422743249837177174247346&
        &21691485290294292900319346665908243383809435507599683357023000500383728&
        &0634351_wp,&
        &0.433395394129247190799265943165784162200071837656246496502701513143766&
        &98907770350122510275795011772122368293504099893794727422475772324920512&
        &67741032822086200952319270933462032011328320387691584063411149801129823&
        &14148878744320432476641442157678880770848387945248811854979703928792696&
        &4254222_wp,&
        &0.679409568299024406234327365114873575769294711834809467664817188952558&
        &57539507492461507857357048037949983390204739931506083674084257663009076&
        &82741718202923543197852846977409718369143712013552962837733153108679126&
        &93254495485472934132472721168027426848661712101171203022718105101071880&
        &4444161_wp,&
        &0.865063366688984510732096688423493048527543014965330452521959731845374&
        &75513805556135679072894604577069440463108641176516867830016149345356373&
        &92729396890950011571349689893051612072435760480900979725923317923795535&
        &73929059587977695683242770223694276591148364371481692378170157259728913&
        &9322313_wp,&
        &0.973906528517171720077964012084452053428269946692382119231212066696595&
        &20323463615962572356495626855625823304251877421121502216860143447777992&
        &05409587259942436704413695764881258799146633143510758737119877875210567&
        &06745243536871368303386090938831164665358170712568697066873725922944928&
        &4383797_wp]

    real(wp),dimension(9),parameter :: gwts = [ &
        &1.0_wp,&
        &0.467913934572691047389870343989550994811655605769210535311625319963914&
        &20162039812703111009258479198230476626878975479710092836255417350295459&
        &35635592733866593364825926382559018030281273563502536241704619318259000&
        &99756987095900533474080074634376824431808173206369174103416261765346292&
        &7888917150514997832_wp,&
        &0.360761573048138607569833513837716111661521892746745482289739240237140&
        &03783726171832096220198881934794311720914037079858987989027836432107077&
        &67872114085818922114502722525757771126000732368828591631602895111800517&
        &40813685547074482472486101183259931449817216402425586777526768199930950&
        &3106873150514997832_wp,&
        &0.171324492379170345040296142172732893526822501484043982398635439798945&
        &76054234015464792770542638866975211652206987440430919174716746217597462&
        &96492293180314484520671351091683210843717994067668872126692485569940481&
        &59429327357024984053433824182363244118374610391205239119044219703570297&
        &7497812150514997832_wp,&
        &0.295524224714752870173892994651338329421046717026853601354308029755995&
        &93821715232927035659579375421672271716440125255838681849078955200582600&
        &19363424941869666095627186488841680432313050615358674090830512706638652&
        &87483901746874726597515954450775158914556548308329986393605934912382356&
        &670244_wp,&
        &0.269266719309996355091226921569469352859759938460883795800563276242153&
        &43231917927676422663670925276075559581145036869830869292346938114524155&
        &64658846634423711656014432259960141729044528030344411297902977067142537&
        &53480628460839927657500691168674984281408628886853320804215041950888191&
        &6391898_wp,&
        &0.219086362515982043995534934228163192458771870522677089880956543635199&
        &91065295128124268399317720219278659121687281288763476662690806694756883&
        &09211843316656677105269915322077536772652826671027878246851010208832173&
        &32006427348325475625066841588534942071161341022729156547776892831330068&
        &8702802_wp,&
        &0.149451349150580593145776339657697332402556639669427367835477268753238&
        &65472663001094594726463473195191400575256104543633823445170674549760147&
        &13716011937109528798134828865118770953566439639333773939909201690204649&
        &08381561877915752257830034342778536175692764212879241228297015017259084&
        &2897331_wp,&
        &0.066671344308688137593568809893331792857864834320158145128694881613412&
        &06408408710177678550968505887782109005471452041933148750712625440376213&
        &93049873169940416344953637064001870112423155043935262424506298327181987&
        &18647480566044117862086478449236378557180717569208295026105115288152794&
        &421677_wp]

    iflag = 0_ip
    bquad = 0.0_wp

    if ( k<1_ip .or. k>20_ip ) then

        iflag = 901_ip ! error return

    else if ( n<k ) then

        iflag = 902_ip ! error return

    else

        aa = min(x1,x2)
        bb = max(x1,x2)
        if ( aa>=t(k) ) then
            np1 = n + 1_ip
            if ( bb<=t(np1) ) then
            if ( aa==bb ) return
            npk = n + k
            ! selection of 2, 6, or 10 point gauss formula
            jf = 0_ip
            mf = 1_ip
            if ( k>4_ip ) then
                jf = 1_ip
                mf = 3_ip
                if ( k>12_ip ) then
                    jf = 4_ip
                    mf = 5_ip
                end if
            end if
            do i = 1_ip , mf
                s(i) = 0.0_wp
            end do
            ilo = 1_ip
            inbv = 1_ip
            call dintrv(t,npk,aa,ilo,il1,mflag)
            call dintrv(t,npk,bb,ilo,il2,mflag)
            if ( il2>=np1 ) il2 = n
            do left = il1 , il2
                ta = t(left)
                tb = t(left+1_ip)
                if ( ta/=tb ) then
                    a = max(aa,ta)
                    b = min(bb,tb)
                    bma = 0.5_wp*(b-a)
                    bpa = 0.5_wp*(b+a)
                    do m = 1_ip , mf
                        c1 = bma*gpts(jf+m)
                        gx = -c1 + bpa
                        call dbvalu(t,bcoef,n,k,0_ip,gx,inbv,work,iflag,y2)
                        if (iflag/=0_ip) return
                        gx = c1 + bpa
                        call dbvalu(t,bcoef,n,k,0_ip,gx,inbv,work,iflag,y1)
                        if (iflag/=0_ip) return
                        s(m) = s(m) + (y1+y2)*bma
                    end do
                end if
            end do
            q = 0.0_wp
            do m = 1_ip , mf
                q = q + gwts(jf+m)*s(m)
            end do
            if ( x1>x2 ) q = -q
                bquad = q
                return
            end if
        end if

        iflag = 903_ip ! error return

    end if

    end subroutine dbsqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbfqad computes the integral on `(x1,x2)` of a product of a
!  function `f` and the `id`-th derivative of a `k`-th order b-spline,
!  using the b-representation `(t,bcoef,n,k)`.  `(x1,x2)` must be a
!  subinterval of `t(k) <= x <= t(n+1)`.  an integration routine,
!  [[dbsgq8]] (a modification of `gaus8`), integrates the product
!  on subintervals of `(x1,x2)` formed by included (distinct) knots
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
!
!@note the maximum number of significant digits obtainable in
!      [[dbsqad]] is the smaller of ~300 and the number of digits
!      carried in `real(wp)` arithmetic.
!
!@note Extrapolation is not enabled for this routine.

    subroutine dbfqad(f,t,bcoef,n,k,id,x1,x2,tol,quad,iflag,work)

    implicit none

    procedure(b1fqad_func)              :: f      !! external function of one argument for the
                                                  !! integrand `bf(x)=f(x)*dbvalu(t,bcoef,n,k,id,x,inbv,work)`
    integer(ip),intent(in)              :: n      !! length of coefficient array
    integer(ip),intent(in)              :: k      !! order of b-spline, `k >= 1`
    real(wp),dimension(n+k),intent(in)  :: t      !! knot array
    real(wp),dimension(n),intent(in)    :: bcoef  !! coefficient array
    integer(ip),intent(in)              :: id     !! order of the spline derivative, `0 <= id <= k-1`
                                                  !! `id=0` gives the spline function
    real(wp),intent(in)                 :: x1     !! left point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: x2     !! right point of quadrature interval in `t(k) <= x <= t(n+1)`
    real(wp),intent(in)                 :: tol    !! desired accuracy for the quadrature, suggest
                                                  !! `10*dtol < tol <= 0.1` where `dtol` is the maximum
                                                  !! of `1.0e-300` and real(wp) unit roundoff for
                                                  !! the machine
    real(wp),intent(out)                :: quad   !! integral of `bf(x)` on `(x1,x2)`
    real(wp),dimension(:),intent(inout) :: work   !! work vector of length `3*k`
    integer(ip),intent(out)             :: iflag  !! status flag:
                                                  !!
                                                  !! * 0: no errors
                                                  !! * 1001: `k` does not satisfy `k>=1`
                                                  !! * 1002: `n` does not satisfy `n>=k`
                                                  !! * 1003: `d` does not satisfy `0<=id<k`
                                                  !! * 1004: `x1` or `x2` or both do not
                                                  !!   satisfy `t(k)<=x<=t(n+1)`
                                                  !! * 1005: `tol` is less than `dtol`
                                                  !!   or greater than 0.1

    integer(ip) :: inbv,ilo,il1,il2,left,mflag,npk,np1
    real(wp) :: a,aa,ans,b,bb,q,ta,tb,err

    real(wp),parameter :: min_tol = max(epsilon(1.0_wp),1.0e-300_wp) !! minimum allowed `tol`

    iflag = 0_ip
    quad = 0.0_wp
    err = tol
    if ( k<1_ip ) then
        iflag = 1001_ip     ! error
    elseif ( n<k ) then
        iflag = 1002_ip     ! error
    elseif ( id<0_ip .or. id>=k ) then
        iflag = 1003_ip     ! error
    else
        if ( tol>=min_tol .and. tol<=0.1_wp ) then
            aa = min(x1,x2)
            bb = max(x1,x2)
            if ( aa>=t(k) ) then
                np1 = n + 1_ip
                if ( bb<=t(np1) ) then
                    if ( aa==bb ) return
                    npk = n + k
                    ilo = 1_ip
                    call dintrv(t,npk,aa,ilo,il1,mflag)
                    call dintrv(t,npk,bb,ilo,il2,mflag)
                    if ( il2>=np1 ) il2 = n
                    inbv = 1_ip
                    q = 0.0_wp
                    do left = il1 , il2
                        ta = t(left)
                        tb = t(left+1_ip)
                        if ( ta/=tb ) then
                            a = max(aa,ta)
                            b = min(bb,tb)
                            call dbsgq8(f,t,bcoef,n,k,id,a,b,inbv,err,ans,iflag,work)
                            if ( iflag/=0_ip .and. iflag/=1101_ip ) return
                            q = q + ans
                        end if
                    end do
                    if ( x1>x2 ) q = -q
                    quad = q
                end if
            else
                iflag = 1004_ip  ! error
            end if
        else
            iflag = 1005_ip  ! error
        end if
    end if

    end subroutine dbfqad
!*****************************************************************************************

!*****************************************************************************************
!>
!  DBSGQ8, a modification of [gaus8](http://netlib.sandia.gov/slatec/src/gaus8.f),
!  integrates the product of `fun(x)` by the `id`-th derivative of a spline
!  [[dbvalu]] between limits `a` and `b` using an adaptive 8-point Legendre-Gauss
!  algorithm.
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
!  * Jacob Williams, 9/6/2017 : refactored to modern Fortran. Some changes.
!    Added higher precision coefficients.

    subroutine dbsgq8(fun,xt,bc,n,kk,id,a,b,inbv,err,ans,iflag,work)

    implicit none

    procedure(b1fqad_func)              :: fun     !! name of external function of one
                                                   !! argument which multiplies [[dbvalu]].
    integer(ip),intent(in)              :: n       !! number of b-coefficients for [[dbvalu]]
    integer(ip),intent(in)              :: kk      !! order of the spline, `kk>=1`
    real(wp),dimension(:),intent(in)    :: xt      !! knot array for [[dbvalu]]
    real(wp),dimension(n),intent(in)    :: bc      !! b-coefficient array for [[dbvalu]]
    integer(ip),intent(in)              :: id      !! Order of the spline derivative, `0<=id<=kk-1`
    real(wp),intent(in)                 :: a       !! lower limit of integral
    real(wp),intent(in)                 :: b       !! upper limit of integral (may be less than `a`)
    integer(ip),intent(inout)           :: inbv    !! initialization parameter for [[dbvalu]]
    real(wp),intent(inout)              :: err     !! **IN:** is a requested pseudorelative error
                                                   !! tolerance.  normally pick a value of
                                                   !! `abs(err)<1e-3`.  `ans` will normally
                                                   !! have no more error than `abs(err)` times
                                                   !! the integral of the absolute value of
                                                   !! `fun(x)*[[dbvalu]]()`.
                                                   !!
                                                   !! **OUT:** will be an estimate of the absolute
                                                   !! error in ans if the input value of `err`
                                                   !! was negative.  (`err` is unchanged if
                                                   !! the input value of `err` was nonnegative.)
                                                   !! the estimated error is solely for information
                                                   !! to the user and should not be used as a
                                                   !! correction to the computed integral.
    real(wp),intent(out)                :: ans     !! computed value of integral
    integer(ip),intent(out)             :: iflag   !! a status code:
                                                   !!
                                                   !! * 0: `ans` most likely meets requested
                                                   !!   error tolerance, or `a=b`.
                                                   !! * 1101: `a` and `b` are too nearly equal
                                                   !!   to allow normal integration.
                                                   !!   `ans` is set to zero.
                                                   !! * 1102: `ans` probably does not meet
                                                   !!   requested error tolerance.
    real(wp),dimension(:),intent(inout) :: work    !! work vector of length `3*k` for [[dbvalu]]

    integer(ip) :: k,l,lmn,lmx,mxl,nbits,nib,nlmx
    real(wp) :: ae,anib,area,c,ce,ee,ef,eps,est,gl,glr,tol,vr,x
    integer(ip),dimension(60)  :: lr
    real(wp),dimension(60) :: aa,hh,vl,gr

    integer(ip),parameter  :: i1mach14 = digits(1.0_wp)            !! i1mach(14)
    real(wp),parameter     :: d1mach5  = log10(real(radix(x),wp))  !! d1mach(5)
    real(wp),parameter     :: ln2      = log(2.0_wp)               !! 0.69314718d0
    real(wp),parameter     :: sq2      = sqrt(2.0_wp)
    integer(ip),parameter  :: nlmn     = 1
    integer(ip),parameter  :: kmx      = 5000
    integer(ip),parameter  :: kml      = 6

    ! initialize
    inbv  = 1_ip
    iflag = 0_ip
    k     = i1mach14
    anib  = d1mach5*k/0.30102000_wp
    nbits = int(anib,ip)
    nlmx  = min((nbits*5_ip)/8_ip,60_ip)
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
                        nib = int(anib,ip)
                        lmx = min(nlmx,nbits-nib-7_ip)
                        if ( lmx<1_ip ) then
                            ! a and b are too nearly equal
                            ! to allow normal integration
                            iflag = 1101_ip
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
        if ( err==0.0_wp ) tol = sqrt(epsilon(1.0_wp))
        eps = tol
        hh(1_ip) = (b-a)/4.0_wp
        aa(1_ip) = a
        lr(1_ip) = 1_ip
        l = 1_ip
        call g8(aa(l)+2.0_wp*hh(l),2.0_wp*hh(l),est,iflag)
        if (iflag/=0_ip) return
        k = 8_ip
        area = abs(est)
        ef = 0.5_wp
        mxl = 0_ip
    end if

    do
        ! compute refined estimates, estimate the error, etc.
        call g8(aa(l)+hh(l),hh(l),gl,iflag)
        if (iflag/=0_ip) return
        call g8(aa(l)+3.0_wp*hh(l),hh(l),gr(l),iflag)
        if (iflag/=0_ip) return
        k = k + 16_ip
        area = area + (abs(gl)+abs(gr(l))-abs(est))
        glr = gl + gr(l)
        ee = abs(est-glr)*ef
        ae = max(eps*area,tol*abs(glr))
        if ( ee>ae ) then
            ! consider the left half of this level
            if ( k>kmx ) lmx = kml
            if ( l>=lmx ) then
                mxl = 1_ip
            else
                l = l + 1_ip
                eps = eps*0.5_wp
                ef = ef/sq2
                hh(l) = hh(l-1)*0.5_wp
                lr(l) = -1_ip
                aa(l) = aa(l-1_ip)
                est = gl
                cycle
            end if
        end if
        ce = ce + (est-glr)
        if ( lr(l)<=0_ip ) then
            ! proceed to right half at this level
            vl(l) = glr
        else
            ! return one level
            vr = glr
            do
                if ( l<=1_ip ) then
                    ! exit
                    ans = vr
                    if ( (mxl/=0_ip) .and. (abs(ce)>2.0_wp*tol*area) ) then
                        iflag = 1102_ip
                    end if
                    if ( err<0.0_wp ) err = ce
                    return
                else
                    l = l - 1_ip
                    eps = eps*2.0_wp
                    ef = ef*sq2
                    if ( lr(l)<=0 ) then
                        vl(l) = vl(l+1_ip) + vr
                        exit
                    else
                        vr = vl(l+1_ip) + vr
                    end if
                end if
            end do
        end if
        est = gr(l-1_ip)
        lr(l) = 1_ip
        aa(l) = aa(l) + 4.0_wp*hh(l)
    end do

    contains

        subroutine g8(x,h,res,iflag)

        !! 8-point formula.
        !!
        !!@note Replaced the original double precision abscissa and weight
        !!      coefficients with the higher precision versions from here:
        !!      http://pomax.github.io/bezierinfo/legendre-gauss.html
        !!      So, if `wp` is changed to say, `real128`, more precision
        !!      can be obtained. These coefficients have about 300 digits.

        implicit none

        real(wp),intent(in)     :: x
        real(wp),intent(in)     :: h
        real(wp),intent(out)    :: res
        integer(ip),intent(out) :: iflag

        real(wp),dimension(8) :: f
        real(wp),dimension(8) :: v

        ! abscissa and weight coefficients:
        real(wp),parameter :: x1 = &
        &0.1834346424956498049394761423601839806667578129129737823171884736992044&
        &742215421141160682237111233537452676587642867666089196012523876865683788&
        &569995160663568104475551617138501966385810764205532370882654749492812314&
        &961247764619363562770645716456613159405134052985058171969174306064445289&
        &638150514997832_wp
        real(wp),parameter :: x2 = &
        &0.5255324099163289858177390491892463490419642431203928577508570992724548&
        &207685612725239614001936319820619096829248252608507108793766638779939805&
        &395303668253631119018273032402360060717470006127901479587576756241288895&
        &336619643528330825624263470540184224603688817537938539658502113876953598&
        &879150514997832_wp
        real(wp),parameter :: x3 = &
        &0.7966664774136267395915539364758304368371717316159648320701702950392173&
        &056764730921471519272957259390191974534530973092653656494917010859602772&
        &562074621689676153935016290342325645582634205301545856060095727342603557&
        &415761265140428851957341933710803722783136113628137267630651413319993338&
        &002150514997832_wp
        real(wp),parameter :: x4 = &
        &0.9602898564975362316835608685694729904282352343014520382716397773724248&
        &977434192844394389592633122683104243928172941762102389581552171285479373&
        &642204909699700433982618326637346808781263553346927867359663480870597542&
        &547603929318533866568132868842613474896289232087639988952409772489387324&
        &25615051499783203_wp
        real(wp),parameter :: w1 = &
        &0.3626837833783619829651504492771956121941460398943305405248230675666867&
        &347239066773243660420848285095502587699262967065529258215569895173844995&
        &576007862076842778350382862546305771007553373269714714894268328780431822&
        &779077846722965535548199601402487767505928976560993309027632737537826127&
        &502150514997832_wp
        real(wp),parameter :: w2 = &
        &0.3137066458778872873379622019866013132603289990027349376902639450749562&
        &719421734969616980762339285560494275746410778086162472468322655616056890&
        &624276469758994622503118776562559463287222021520431626467794721603822601&
        &295276898652509723185157998353156062419751736972560423953923732838789657&
        &919150514997832_wp
        real(wp),parameter :: w3 = &
        &0.2223810344533744705443559944262408844301308700512495647259092892936168&
        &145704490408536531423771979278421592661012122181231114375798525722419381&
        &826674532090577908613289536840402789398648876004385697202157482063253247&
        &195590228631570651319965589733545440605952819880671616779621183704306688&
        &233150514997832_wp
        real(wp),parameter :: w4 = &
        &0.1012285362903762591525313543099621901153940910516849570590036980647401&
        &787634707848602827393040450065581543893314132667077154940308923487678731&
        &973041136073584690533208824050731976306575729205467961435779467552492328&
        &730055025992954089946676810510810729468366466585774650346143712142008566&
        &866150514997832_wp

        res = 0.0_wp

        v(1_ip) = x-x1*h
        v(2_ip) = x+x1*h
        v(3_ip) = x-x2*h
        v(4_ip) = x+x2*h
        v(5_ip) = x-x3*h
        v(6_ip) = x+x3*h
        v(7_ip) = x-x4*h
        v(8_ip) = x+x4*h

        call dbvalu(xt,bc,n,kk,id,v(1_ip),inbv,work,iflag,f(1_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(2_ip),inbv,work,iflag,f(2_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(3_ip),inbv,work,iflag,f(3_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(4_ip),inbv,work,iflag,f(4_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(5_ip),inbv,work,iflag,f(5_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(6_ip),inbv,work,iflag,f(6_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(7_ip),inbv,work,iflag,f(7_ip)); if (iflag/=0_ip) return
        call dbvalu(xt,bc,n,kk,id,v(8_ip),inbv,work,iflag,f(8_ip)); if (iflag/=0_ip) return

        res = h*((w1*(fun(v(1_ip))*f(1_ip) + fun(v(2_ip))*f(2_ip))  + &
                  w2*(fun(v(3_ip))*f(3_ip) + fun(v(4_ip))*f(4_ip))) + &
                 (w3*(fun(v(5_ip))*f(5_ip) + fun(v(6_ip))*f(6_ip))  + &
                  w4*(fun(v(7_ip))*f(7_ip) + fun(v(8_ip))*f(8_ip))))

        end subroutine g8

    end subroutine dbsgq8
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the value of `x` to use for computing the interval
!  in `t`, depending on if extrapolation is allowed or not.
!
!  If extrapolation is allowed and x is < tmin or > tmax, then either
!  `tmin` or `tmax - 2.0_wp*spacing(tmax)` is returned.
!  Otherwise, `x` is returned.

    pure function get_temp_x_for_extrap(x,tmin,tmax,extrap) result(xt)

    implicit none

    real(wp),intent(in) :: x    !! variable value
    real(wp),intent(in) :: tmin !! first knot vector element for b-splines
    real(wp),intent(in) :: tmax !! last knot vector element for b-splines
    real(wp)            :: xt   !! The value returned (it will either
                                !! be `tmin`, `x`, or `tmax`)
    logical,intent(in),optional :: extrap  !! if extrapolation is allowed
                                           !! (if not present, default is False)

    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    if (extrapolation_allowed) then
        if (x<tmin) then
            xt = tmin
        else if (x>tmax) then
            ! Put it just inside the upper bound.
            ! This is sort of a hack to get
            ! extrapolation to work.
            xt = tmax - 2.0_wp*spacing(tmax)
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

    integer(ip),intent(in)       :: iflag  !! return code from one of the routines
    character(len=:),allocatable :: msg    !! status message associated with the flag

    character(len=10) :: istr   !! for integer to string conversion
    integer(ip)       :: istat  !! for write statement

    select case (iflag)

    case(  0_ip); msg='Successful execution'

    case( -1_ip); msg='Error in dintrv: x < xt(1_ip)'
    case( -2_ip); msg='Error in dintrv: x >= xt(lxt)'

    case(  1_ip); msg='Error in evaluate_*d: class is not initialized'

    case(  2_ip); msg='Error in db*ink: iknot out of range'
    case(  3_ip); msg='Error in db*ink: nx out of range'
    case(  4_ip); msg='Error in db*ink: kx out of range'
    case(  5_ip); msg='Error in db*ink: x not strictly increasing'
    case(  6_ip); msg='Error in db*ink: tx not non-decreasing'
    case(  7_ip); msg='Error in db*ink: ny out of range'
    case(  8_ip); msg='Error in db*ink: ky out of range'
    case(  9_ip); msg='Error in db*ink: y not strictly increasing'
    case( 10_ip); msg='Error in db*ink: ty not non-decreasing'
    case( 11_ip); msg='Error in db*ink: nz out of range'
    case( 12_ip); msg='Error in db*ink: kz out of range'
    case( 13_ip); msg='Error in db*ink: z not strictly increasing'
    case( 14_ip); msg='Error in db*ink: tz not non-decreasing'
    case( 15_ip); msg='Error in db*ink: nq out of range'
    case( 16_ip); msg='Error in db*ink: kq out of range'
    case( 17_ip); msg='Error in db*ink: q not strictly increasing'
    case( 18_ip); msg='Error in db*ink: tq not non-decreasing'
    case( 19_ip); msg='Error in db*ink: nr out of range'
    case( 20_ip); msg='Error in db*ink: kr out of range'
    case( 21_ip); msg='Error in db*ink: r not strictly increasing'
    case( 22_ip); msg='Error in db*ink: tr not non-decreasing'
    case( 23_ip); msg='Error in db*ink: ns out of range'
    case( 24_ip); msg='Error in db*ink: ks out of range'
    case( 25_ip); msg='Error in db*ink: s not strictly increasing'
    case( 26_ip); msg='Error in db*ink: ts not non-decreasing'
    case(700_ip); msg='Error in db*ink: size(x) /= size(fcn,1)'
    case(701_ip); msg='Error in db*ink: size(y) /= size(fcn,2)'
    case(702_ip); msg='Error in db*ink: size(z) /= size(fcn,3)'
    case(703_ip); msg='Error in db*ink: size(q) /= size(fcn,4)'
    case(704_ip); msg='Error in db*ink: size(r) /= size(fcn,5)'
    case(705_ip); msg='Error in db*ink: size(s) /= size(fcn,6)'
    case(706_ip); msg='Error in db*ink: size(x) /= nx'
    case(707_ip); msg='Error in db*ink: size(y) /= ny'
    case(708_ip); msg='Error in db*ink: size(z) /= nz'
    case(709_ip); msg='Error in db*ink: size(q) /= nq'
    case(710_ip); msg='Error in db*ink: size(r) /= nr'
    case(711_ip); msg='Error in db*ink: size(s) /= ns'
    case(712_ip); msg='Error in db*ink: size(tx) /= nx+kx'
    case(713_ip); msg='Error in db*ink: size(ty) /= ny+ky'
    case(714_ip); msg='Error in db*ink: size(tz) /= nz+kz'
    case(715_ip); msg='Error in db*ink: size(tq) /= nq+kq'
    case(716_ip); msg='Error in db*ink: size(tr) /= nr+kr'
    case(717_ip); msg='Error in db*ink: size(ts) /= ns+ks'
    case(800_ip); msg='Error in db*ink: size(x) /= size(bcoef,1)'
    case(801_ip); msg='Error in db*ink: size(y) /= size(bcoef,2)'
    case(802_ip); msg='Error in db*ink: size(z) /= size(bcoef,3)'
    case(803_ip); msg='Error in db*ink: size(q) /= size(bcoef,4)'
    case(804_ip); msg='Error in db*ink: size(r) /= size(bcoef,5)'
    case(805_ip); msg='Error in db*ink: size(s) /= size(bcoef,6)'

    case(806_ip); msg='Error in dbint4: currently, only k=4 can be used'

    case(100_ip); msg='Error in dbintk: k does not satisfy k>=1'
    case(101_ip); msg='Error in dbintk: n does not satisfy n>=k'
    case(102_ip); msg='Error in dbintk: x(i) does not satisfy x(i)<x(i+1) for some i'
    case(103_ip); msg='Error in dbintk: some abscissa was not in the support of the '//&
                      'corresponding basis function and the system is singular'
    case(104_ip); msg='Error in dbintk: the system of solver detects a singular system '//&
                      'although the theoretical conditions for a solution were satisfied'

    case(201_ip); msg='Error in dbspvn: k does not satisfy k>=1'
    case(202_ip); msg='Error in dbspvn: jhigh does not satisfy 1<=jhigh<=k'
    case(203_ip); msg='Error in dbspvn: index is not 1 or 2'
    case(204_ip); msg='Error in dbspvn: x does not satisfy t(ileft)<=x<=t(ileft+1)'

    case(301_ip); msg='Error in dbtpcf: n should be > 0'

    case(401_ip); msg='Error in dbvalu: k does not satisfy k>=1'
    case(402_ip); msg='Error in dbvalu: n does not satisfy n>=k'
    case(403_ip); msg='Error in dbvalu: ideriv does not satisfy 0<=ideriv<k'
    case(404_ip); msg='Error in dbvalu: x is not greater than or equal to t(k)'
    case(405_ip); msg='Error in dbvalu: x is not less than or equal to t(n+1)'
    case(406_ip); msg='Error in dbvalu: a left limiting value cannot be obtained at t(k)'

    case(501_ip); msg='Error in initialize_*d_specify_knots: tx is not the correct size (kx+nx)'
    case(502_ip); msg='Error in initialize_*d_specify_knots: ty is not the correct size (ky+ny)'
    case(503_ip); msg='Error in initialize_*d_specify_knots: tz is not the correct size (kz+nz)'
    case(504_ip); msg='Error in initialize_*d_specify_knots: tq is not the correct size (kq+nq)'
    case(505_ip); msg='Error in initialize_*d_specify_knots: tr is not the correct size (kr+nr)'
    case(506_ip); msg='Error in initialize_*d_specify_knots: ts is not the correct size (ks+ns)'

    case(601_ip); msg='Error in db*val: x value out of bounds'
    case(602_ip); msg='Error in db*val: y value out of bounds'
    case(603_ip); msg='Error in db*val: z value out of bounds'
    case(604_ip); msg='Error in db*val: q value out of bounds'
    case(605_ip); msg='Error in db*val: r value out of bounds'
    case(606_ip); msg='Error in db*val: s value out of bounds'

    case(901_ip); msg='Error in dbsqad: k does not satisfy 1<=k<=20'
    case(902_ip); msg='Error in dbsqad: n does not satisfy n>=k'
    case(903_ip); msg='Error in dbsqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'

    case(1001_ip); msg='Error in dbfqad: k does not satisfy k>=1'
    case(1002_ip); msg='Error in dbfqad: n does not satisfy n>=k'
    case(1003_ip); msg='Error in dbfqad: d does not satisfy 0<=id<k'
    case(1004_ip); msg='Error in dbfqad: x1 or x2 or both do not satisfy t(k)<=x<=t(n+1)'
    case(1005_ip); msg='Error in dbfqad: tol is less than dtol or greater than 0.1'

    case(1101_ip); msg='Warning in dbsgq8: a and b are too nearly equal to allow normal integration.'
    case(1102_ip); msg='Error in dbsgq8: ans is probably insufficiently accurate.'

    case(2001_ip); msg='Error in dbint4: ndata is less than 2'
    case(2002_ip); msg='Error in dbint4: x values are not distinct or not ordered'
    case(2003_ip); msg='Error in dbint4: ibcl is not 1 or 2'
    case(2004_ip); msg='Error in dbint4: ibcr is not 1 or 2'
    case(2005_ip); msg='Error in dbint4: kntopt is not 1, 2, or 3'
    case(2006_ip); msg='Error in dbint4: knot input through tleft, tright is not ordered properly'
    case(2007_ip); msg='Error in dbint4: the system of equations is singular'

    case(3001_ip); msg='Error in dbspvd: k does not satisfy k>=1'
    case(3002_ip); msg='Error in dbspvd: nderiv does not satisfy 1<=nderiv<=k'
    case(3003_ip); msg='Error in dbspvd: ldvnik does not satisfy ldvnik>=k'

    case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select

    end function get_status_message
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_sub_module
!*****************************************************************************************
