!*****************************************************************************************
!>
!  Regression tests for specific bug fixes. Each test reproduces a
!  case that used to fail (or crash with bounds checking enabled).
!
!  * integral with upper limit exactly on an interior knot ([[dbsgq8]])
!  * integral with upper limit past the end of the spline ([[dbfqad]])
!  * integral with upper limit just above a knot ([[dbfqad]] flag 1101)
!  * integral over an interval only a few ulps wide ([[dbfqad]] flag 1101)
!  * 2D-6D evaluation at the last knot in every direction
!  * 2D-6D evaluation with knots that are not clamped (out of bounds access)
!  * [[db1ink]] alt mode with a small number of points
!  * 4D-6D destructors resetting the class state
!  * integral over whole knot intervals that are narrow relative to `|x|`
!    ([[dbfqad]], which used the midpoint rule for these)
!  * [[db1ink]] alt mode with very uneven data spacing near an end
!    (reported as singular, 2007)

    program bspline_regression_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: n = 6  !! number of points in each dimension
    integer(ip),parameter :: k = 4  !! spline order in each dimension
    real(wp),parameter :: tol = 1.0e-10_wp  !! tolerance for value checks

    real(wp),dimension(n) :: x  !! data points (same in each dimension)
    real(wp),dimension(n+k) :: t_clamped  !! clamped knot vector (ends at `x(n)`)
    real(wp),dimension(n+k) :: t_open  !! knot vector that is not clamped on the right
    integer :: nfail  !! number of failed checks
    integer(ip) :: i  !! counter

    nfail = 0
    x = [(real(i,wp), i=1,n)]
    t_clamped = [1,1,1,1,3,4,6,6,6,6]*1.0_wp
    t_open    = [1,1,1,1,3,4,6,7,8,9]*1.0_wp

    call test_integral_to_interior_knot()
    call test_integral_past_end()
    call test_integral_just_above_knot()
    call test_integral_tiny_interval()
    call test_eval_at_last_knot()
    call test_eval_open_knots()
    call test_db1ink_alt_small_n()
    call test_destroy()
    call test_integral_narrow_knot_intervals()
    call test_db1ink_alt_uneven_spacing()

    write(*,*) ''
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'regression test failed'
    else
        write(*,*) 'all regression tests PASSED'
    end if

    contains
!*****************************************************************************************

    subroutine check(name, ok)

    !! record a pass/fail result

    character(len=*),intent(in) :: name
    logical,intent(in) :: ok

    if (ok) then
        write(*,'(A)') 'PASS: '//name
    else
        write(*,'(A)') 'FAIL: '//name
        nfail = nfail + 1
    end if

    end subroutine check

    pure function linear(p) result(f)

    !! linear test function, which every spline of order >= 2 reproduces exactly

    real(wp),dimension(:),intent(in) :: p
    real(wp) :: f

    integer(ip) :: j

    f = 1.0_wp
    do j = 1, size(p)
        f = f + real(j,wp)*p(j)
    end do

    end function linear

    function one(xx) result(f)

    !! integrand multiplier for `fintegral`

    real(wp),intent(in) :: xx
    real(wp) :: f

    f = 1.0_wp + 0.0_wp*xx

    end function one

!*****************************************************************************************
!>
!  `fintegral` with the upper limit exactly on an interior knot.
!  [[dbsgq8]] used to fall through its `a==b` branch with uninitialized
!  variables (crash, or garbage result).

    subroutine test_integral_to_interior_knot()

    type(bspline_1d) :: s
    integer(ip) :: iflag
    real(wp) :: f
    real(wp),parameter :: exact = 26.0_wp/3.0_wp  !! integral of x**2 from 1 to 3

    call s%initialize(x,x**2,k,t_clamped,iflag)   ! t(5)=3 is an interior knot
    call check('integral to interior knot: initialize', iflag==0)

    call s%fintegral(one,0_ip,1.0_wp,3.0_wp,1.0e-12_wp,f,iflag)
    call check('integral to interior knot: fintegral iflag', iflag==0)
    call check('integral to interior knot: fintegral value', abs(f-exact)<tol)

    call s%integral(1.0_wp,3.0_wp,f,iflag)
    call check('integral to interior knot: integral value', iflag==0 .and. abs(f-exact)<tol)

    end subroutine test_integral_to_interior_knot

!*****************************************************************************************
!>
!  `fintegral` with the upper limit past the end of the spline
!  used to return 0 with no error.

    subroutine test_integral_past_end()

    type(bspline_1d) :: s
    integer(ip) :: iflag
    real(wp) :: f

    call s%initialize(x,x**2,k,iflag)
    call check('integral past end: initialize', iflag==0)

    call s%fintegral(one,0_ip,1.0_wp,10.0_wp,1.0e-12_wp,f,iflag)
    call check('integral past end: fintegral returns 1004', iflag==1004)

    call s%fintegral(one,0_ip,10.0_wp,1.0_wp,1.0e-12_wp,f,iflag)
    call check('integral past end: fintegral reversed returns 1004', iflag==1004)

    call s%integral(1.0_wp,10.0_wp,f,iflag)
    call check('integral past end: integral returns error', iflag/=0)

    end subroutine test_integral_past_end

!*****************************************************************************************
!>
!  `fintegral` with the upper limit a few ulps above a knot.
!  The last subinterval is too small for [[dbsgq8]] (flag 1101), which
!  [[dbfqad]] now integrates with the midpoint rule, but the flag used to
!  be returned anyway.

    subroutine test_integral_just_above_knot()

    type(bspline_1d) :: s
    integer(ip) :: iflag
    real(wp) :: f, b
    real(wp),parameter :: exact = 26.0_wp/3.0_wp  !! integral of x**2 from 1 to 3

    call s%initialize(x,x**2,k,t_clamped,iflag)
    b = 3.0_wp + 4.0_wp*spacing(3.0_wp)
    call s%fintegral(one,0_ip,1.0_wp,b,1.0e-12_wp,f,iflag)
    call check('integral just above knot: iflag', iflag==0)
    call check('integral just above knot: value', abs(f-exact)<tol)

    end subroutine test_integral_just_above_knot

!*****************************************************************************************
!>
!  `fintegral` over an interval only a few ulps wide, where every
!  subinterval is too small for [[dbsgq8]] (flag 1101). The result used
!  to be zero, even though the integral is significant for a large integrand.

    subroutine test_integral_tiny_interval()

    type(bspline_1d) :: s
    integer(ip) :: iflag
    real(wp) :: f, a, b, exact
    real(wp),parameter :: c = 1.0e15_wp  !! constant spline value

    call s%initialize(x,[(c, i=1,n)],k,t_clamped,iflag)
    call check('integral tiny interval: initialize', iflag==0)

    ! whole interval inside one knot interval
    a = 1.0_wp
    b = 1.0_wp + 4.0_wp*spacing(1.0_wp)
    exact = c*(b-a)
    call s%fintegral(one,0_ip,a,b,1.0e-12_wp,f,iflag)
    call check('integral tiny interval: iflag', iflag==0)
    call check('integral tiny interval: value', abs(f-exact)<=1.0e-12_wp*exact)
    call s%fintegral(one,0_ip,b,a,1.0e-12_wp,f,iflag)
    call check('integral tiny interval: reversed value', iflag==0 .and. abs(f+exact)<=1.0e-12_wp*exact)

    ! interval straddling the interior knot at 3
    a = 3.0_wp - 2.0_wp*spacing(3.0_wp)
    b = 3.0_wp + 2.0_wp*spacing(3.0_wp)
    exact = c*(b-a)
    call s%fintegral(one,0_ip,a,b,1.0e-12_wp,f,iflag)
    call check('integral tiny interval across knot: iflag', iflag==0)
    call check('integral tiny interval across knot: value', abs(f-exact)<=1.0e-12_wp*exact)

    end subroutine test_integral_tiny_interval

!*****************************************************************************************
!>
!  2D-6D evaluation at the last knot (here also the last data point).
!  This used to fail with `iflag=-2` in every direction except `x`.

    subroutine test_eval_at_last_knot()

    real(wp),dimension(6) :: p
    integer :: idim, j
    logical :: extrap
    character(len=100) :: name

    do j = 0, 1
        extrap = j==1
        do idim = 1, 6
            ! last knot in direction idim, interior elsewhere
            p = 2.5_wp
            p(idim) = x(n)
            write(name,'(A,I0,A,L1)') 'eval at last knot: direction ', idim, ', extrap=', extrap
            call eval_all_dims(trim(name), t_clamped, p, extrap, expect_ok=.true.)
        end do
        p = x(n)
        write(name,'(A,L1)') 'eval at last knot: all directions, extrap=', extrap
        call eval_all_dims(trim(name), t_clamped, p, extrap, expect_ok=.true.)
    end do

    end subroutine test_eval_at_last_knot

!*****************************************************************************************
!>
!  2D-6D evaluation with knots that are not clamped on the right.
!  Points in `(t(n+1), t(n+k)]` used to pass the range check and then
!  index `bcoef` and the knot vector out of bounds.

    subroutine test_eval_open_knots()

    real(wp),dimension(6) :: p
    integer :: idim
    character(len=100) :: name

    do idim = 1, 6
        p = 2.5_wp
        p(idim) = 8.5_wp   ! inside [t(1),t(n+k)] but past t(n+1)
        write(name,'(A,I0)') 'open knots, past t(n+1): direction ', idim
        call eval_all_dims(trim(name), t_open, p, .false., expect_ok=.false.)
        write(name,'(A,I0)') 'open knots, past t(n+1), extrapolated: direction ', idim
        call eval_all_dims(trim(name), t_open, p, .true., expect_ok=.true.)
        p(idim) = 5.5_wp   ! normal interior point
        write(name,'(A,I0)') 'open knots, interior: direction ', idim
        call eval_all_dims(trim(name), t_open, p, .false., expect_ok=.true.)
    end do

    end subroutine test_eval_open_knots

!*****************************************************************************************
!>
!  Initialize 2D-6D splines of the linear function with knots `t` in each
!  direction, evaluate them at the first `d` components of `p`, and check either
!  that the value is exact (`expect_ok`) or that an error is returned.

    subroutine eval_all_dims(name, t, p, extrap, expect_ok)

    character(len=*),intent(in) :: name
    real(wp),dimension(n+k),intent(in) :: t
    real(wp),dimension(6),intent(in) :: p
    logical,intent(in) :: extrap
    logical,intent(in) :: expect_ok

    type(bspline_2d) :: s2
    type(bspline_3d) :: s3
    type(bspline_4d) :: s4
    type(bspline_5d) :: s5
    type(bspline_6d) :: s6
    real(wp),dimension(n,n) :: f2
    real(wp),dimension(n,n,n) :: f3
    real(wp),dimension(n,n,n,n) :: f4
    real(wp),dimension(n,n,n,n,n) :: f5
    real(wp),dimension(n,n,n,n,n,n) :: f6
    integer(ip) :: i1,i2,i3,i4,i5,i6,iflag,iflag_init
    real(wp) :: f

    do i6=1,n; do i5=1,n; do i4=1,n; do i3=1,n; do i2=1,n; do i1=1,n
        f6(i1,i2,i3,i4,i5,i6) = linear([x(i1),x(i2),x(i3),x(i4),x(i5),x(i6)])
    end do; end do; end do; end do; end do; end do
    f5 = f6(:,:,:,:,:,1) - 6.0_wp*x(1)
    f4 = f5(:,:,:,:,1)   - 5.0_wp*x(1)
    f3 = f4(:,:,:,1)     - 4.0_wp*x(1)
    f2 = f3(:,:,1)       - 3.0_wp*x(1)

    call s2%initialize(x,x,f2,k,k,t,t,iflag_init,extrap)
    call s2%evaluate(p(1),p(2),0_ip,0_ip,f,iflag)
    call report(name, 2, p, t, expect_ok, iflag_init, iflag, f)

    call s3%initialize(x,x,x,f3,k,k,k,t,t,t,iflag_init,extrap)
    call s3%evaluate(p(1),p(2),p(3),0_ip,0_ip,0_ip,f,iflag)
    call report(name, 3, p, t, expect_ok, iflag_init, iflag, f)

    call s4%initialize(x,x,x,x,f4,k,k,k,k,t,t,t,t,iflag_init,extrap)
    call s4%evaluate(p(1),p(2),p(3),p(4),0_ip,0_ip,0_ip,0_ip,f,iflag)
    call report(name, 4, p, t, expect_ok, iflag_init, iflag, f)

    call s5%initialize(x,x,x,x,x,f5,k,k,k,k,k,t,t,t,t,t,iflag_init,extrap)
    call s5%evaluate(p(1),p(2),p(3),p(4),p(5),0_ip,0_ip,0_ip,0_ip,0_ip,f,iflag)
    call report(name, 5, p, t, expect_ok, iflag_init, iflag, f)

    call s6%initialize(x,x,x,x,x,x,f6,k,k,k,k,k,k,t,t,t,t,t,t,iflag_init,extrap)
    call s6%evaluate(p(1),p(2),p(3),p(4),p(5),p(6),0_ip,0_ip,0_ip,0_ip,0_ip,0_ip,f,iflag)
    call report(name, 6, p, t, expect_ok, iflag_init, iflag, f)

    end subroutine eval_all_dims

!*****************************************************************************************
!>
!  Check the result of one of the [[eval_all_dims]] cases for a `d`-dimensional spline.

    subroutine report(name, d, p, t, expect_ok, iflag_init, iflag, f)

    character(len=*),intent(in) :: name
    integer,intent(in) :: d
    real(wp),dimension(6),intent(in) :: p
    real(wp),dimension(n+k),intent(in) :: t
    logical,intent(in) :: expect_ok
    integer(ip),intent(in) :: iflag_init, iflag
    real(wp),intent(in) :: f

    character(len=:),allocatable :: label

    label = name//' ['//achar(iachar('0')+d)//'D]'

    if (iflag_init/=0) then
        call check(label//' initialize', .false.)
    else if (expect_ok .or. all(p(1:d)<=t(n+1))) then
        ! (only the directions present in this spline can be out of range)
        call check(label, iflag==0 .and. abs(f-linear(p(1:d)))<tol)
    else
        call check(label//' returns error', iflag/=0)
    end if

    end subroutine report

!*****************************************************************************************
!>
!  [[db1ink]] alt mode (derivative boundary conditions) with few points.
!  The order check used to reject any `nx <= 6`.

    subroutine test_db1ink_alt_small_n()

    integer(ip) :: nx, iflag, inbvx, j
    real(wp),dimension(:),allocatable :: xx, yy, tx, bcoef
    real(wp),dimension(3*k) :: w0
    real(wp) :: f, df
    character(len=100) :: name

    do nx = 2, 6
        xx = [(real(j,wp), j=1,nx)]
        yy = xx**3
        allocate(tx(nx+6), bcoef(nx+2))
        ! first derivative constrained at both ends (exact for x**3)
        call db1ink(xx,nx,yy,k,1_ip,1_ip,3.0_wp*xx(1)**2,3.0_wp*xx(nx)**2,1_ip,tx,bcoef,iflag)
        write(name,'(A,I0)') 'db1ink alt mode: nx=', nx
        call check(trim(name)//' initialize', iflag==0)
        if (iflag==0) then
            ! the cubic is reproduced exactly
            inbvx = 1
            call db1val(2.5_wp-1.0_wp/real(nx,wp),0_ip,tx,nx,nx+2,k,bcoef,f,iflag,inbvx,w0)
            call check(trim(name)//' value', &
                       iflag==0 .and. abs(f-(2.5_wp-1.0_wp/real(nx,wp))**3)<tol)
            call db1val(xx(nx),1_ip,tx,nx,nx+2,k,bcoef,df,iflag,inbvx,w0)
            call check(trim(name)//' end derivative', iflag==0 .and. abs(df-3.0_wp*xx(nx)**2)<tol)
        end if
        deallocate(tx, bcoef)
    end do

    end subroutine test_db1ink_alt_small_n

!*****************************************************************************************
!>
!  4D-6D destructors used to leave the class marked as initialized,
!  so evaluating after `destroy` would use deallocated arrays.

    subroutine test_destroy()

    type(bspline_4d) :: s4
    type(bspline_5d) :: s5
    type(bspline_6d) :: s6
    real(wp),dimension(n,n,n,n) :: f4
    real(wp),dimension(n,n,n,n,n) :: f5
    real(wp),dimension(n,n,n,n,n,n) :: f6
    integer(ip) :: iflag
    real(wp) :: f

    f4 = 1.0_wp; f5 = 1.0_wp; f6 = 1.0_wp

    call s4%initialize(x,x,x,x,f4,k,k,k,k,iflag)
    call s4%destroy()
    call check('destroy 4D: status_ok is false', .not. s4%status_ok())
    call s4%evaluate(2.0_wp,2.0_wp,2.0_wp,2.0_wp,0_ip,0_ip,0_ip,0_ip,f,iflag)
    call check('destroy 4D: evaluate returns error', iflag/=0)

    call s5%initialize(x,x,x,x,x,f5,k,k,k,k,k,iflag)
    call s5%destroy()
    call check('destroy 5D: status_ok is false', .not. s5%status_ok())
    call s5%evaluate(2.0_wp,2.0_wp,2.0_wp,2.0_wp,2.0_wp,0_ip,0_ip,0_ip,0_ip,0_ip,f,iflag)
    call check('destroy 5D: evaluate returns error', iflag/=0)

    call s6%initialize(x,x,x,x,x,x,f6,k,k,k,k,k,k,iflag)
    call s6%destroy()
    call check('destroy 6D: status_ok is false', .not. s6%status_ok())
    call s6%evaluate(2.0_wp,2.0_wp,2.0_wp,2.0_wp,2.0_wp,2.0_wp,0_ip,0_ip,0_ip,0_ip,0_ip,0_ip,f,iflag)
    call check('destroy 6D: evaluate returns error', iflag/=0)

    end subroutine test_destroy

!*****************************************************************************************
!>
!  Integrals over whole knot intervals that are narrow relative to `|x|`
!  (too narrow for [[dbsgq8]]): the spline is not constant on them, so the
!  midpoint rule (used in an earlier version) gave the wrong result with
!  no error flag.

    subroutine test_integral_narrow_knot_intervals()

    real(wp) :: h         !! width of the knot interval
    real(wp) :: f         !! integral
    real(wp) :: x0        !! offset of the knots
    real(wp) :: exact     !! exact integral
    real(wp),dimension(3*4) :: w0  !! work array
    real(wp),dimension(15) :: tx      !! knots for the offset cubic
    real(wp),dimension(11) :: bcoef   !! coefficients for the offset cubic
    integer(ip) :: iflag  !! status flag
    integer(ip) :: i      !! counter

    ! s(x) = ((x-1)/h)**2 on [1,1+h] (k=3): the integral is h/3
    h = 64.0_wp*spacing(1.0_wp)
    call db1fqad(one,[1.0_wp,1.0_wp,1.0_wp,1.0_wp+h,1.0_wp+h,1.0_wp+h],[0.0_wp,0.0_wp,1.0_wp], &
                 3_ip,3_ip,0_ip,1.0_wp,1.0_wp+h,1.0e-12_wp,f,iflag,w0)
    call check('narrow quadratic interval: iflag', iflag==0)
    call check('narrow quadratic interval: value', abs(f-h/3.0_wp)<=1.0e-12_wp*h)

    ! cubic spline with knots near 1e9, spaced 2**-12 apart (exactly
    ! representable, and every knot interval is too narrow for dbsgq8).
    ! The integral over the whole domain is sum(c(i)*(t(i+k)-t(i))/k).
    x0 = 1.0e9_wp
    tx(1:4) = x0
    do i = 5, 11
        tx(i) = x0 + real(i-4,wp)*2.0_wp**(-12)
    end do
    tx(12:15) = x0 + 8.0_wp*2.0_wp**(-12)
    bcoef = [(sin(real(i,wp)), i=1,11)]
    exact = 0.0_wp
    do i = 1, 11
        exact = exact + bcoef(i)*(tx(i+4)-tx(i))/4.0_wp
    end do
    call db1fqad(one,tx,bcoef,11_ip,4_ip,0_ip,tx(4),tx(12),1.0e-12_wp,f,iflag,w0)
    call check('narrow cubic intervals: iflag', iflag==0)
    call check('narrow cubic intervals: value', abs(f-exact)<=1.0e-12_wp*abs(exact))

    end subroutine test_integral_narrow_knot_intervals

!*****************************************************************************************

!*****************************************************************************************
!>
!  Interpolation with end derivative conditions when the data spacing
!  near an end is very uneven: an outer basis function is small at all
!  of the data points, which made the system look singular (2007).

    subroutine test_db1ink_alt_uneven_spacing()

    real(wp),dimension(3) :: xd = [0.0_wp, 1.0_wp, 1.0_wp+1.0e-7_wp]  !! data abscissae
    real(wp),dimension(3) :: yd = [1.0_wp, 2.0_wp, 0.5_wp]  !! data
    real(wp),dimension(3+6) :: tx    !! knots
    real(wp),dimension(3+2) :: bcoef !! coefficients
    real(wp),dimension(3*4) :: w0    !! work array
    real(wp) :: f        !! spline value
    real(wp) :: err      !! maximum interpolation error
    integer(ip) :: iflag !! status flag
    integer(ip) :: inbvx !! initialization parameter for [[db1val]]
    integer(ip) :: i     !! counter

    call db1ink(xd,3_ip,yd,4_ip,1_ip,1_ip,0.3_wp,-0.7_wp,2_ip,tx,bcoef,iflag)
    call check('alt mode, uneven spacing: iflag', iflag==0)
    if (iflag==0) then
        err = 0.0_wp
        do i = 1, 3
            inbvx = 1
            call db1val(xd(i),0_ip,tx,5_ip,4_ip,bcoef,f,iflag,inbvx,w0)
            err = max(err, abs(f-yd(i)))
        end do
        call check('alt mode, uneven spacing: interpolates', err<1.0e-8_wp)
        inbvx = 1
        call db1val(xd(1),1_ip,tx,5_ip,4_ip,bcoef,f,iflag,inbvx,w0)
        call check('alt mode, uneven spacing: left derivative', abs(f-0.3_wp)<1.0e-6_wp)
    end if

    end subroutine test_db1ink_alt_uneven_spacing
!*****************************************************************************************

    end program bspline_regression_test
!*****************************************************************************************
