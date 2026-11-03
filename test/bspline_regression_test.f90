!*****************************************************************************************
!>
!  Regression tests for specific bug fixes. Each test reproduces a
!  case that used to fail (or crash with bounds checking enabled).
!
!  * integral with upper limit exactly on an interior knot ([[dbsgq8]])
!  * integral with upper limit past the end of the spline ([[dbfqad]])
!  * integral with upper limit just above a knot ([[dbfqad]] flag 1101)
!  * 2D-6D evaluation at the last knot in every direction
!  * 2D-6D evaluation with knots that are not clamped (out of bounds access)
!  * [[db1ink]] alt mode with a small number of points
!  * 4D-6D destructors resetting the class state

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
    call test_eval_at_last_knot()
    call test_eval_open_knots()
    call test_db1ink_alt_small_n()
    call test_destroy()

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
!  [[dbfqad]] treats as negligible, but the flag used to be returned anyway.

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

    end program bspline_regression_test
!*****************************************************************************************
