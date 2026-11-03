!*****************************************************************************************
!>
!  Tests for [[db1root]], [[db1min]], [[db1max]], and the
!  `roots`, `inverse`, `minimum`, and `maximum` methods of [[bspline_1d]].
!
!  * polynomials with known roots and extrema (reproduced exactly by the spline)
!  * roots of `s(x)=a` for `a /= 0`
!  * double and triple roots
!  * roots at knots and at the ends of the domain
!  * splines that are constant over part of the domain
!  * random splines, compared with a dense search for sign changes
!  * error codes and the object-oriented interface

    program bspline_roots_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    real(wp),parameter :: pi = acos(-1.0_wp)

    integer :: nfail   !! number of failed checks
    integer :: ncheck  !! number of checks

    nfail = 0
    ncheck = 0

    call test_cubic_roots()
    call test_level()
    call test_sin()
    call test_linear()
    call test_multiple_roots()
    call test_knot_and_end_roots()
    call test_flat()
    call test_cubic_extrema()
    call test_sin_extrema()
    call test_random()
    call test_errors()
    call test_oo()
    call test_review_cases()
    call test_review_cases_2()

    write(*,*) ''
    write(*,*) ncheck, ' checks'
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'roots test failed'
    else
        write(*,*) 'all roots tests PASSED'
    end if

    contains
!*****************************************************************************************

    subroutine check(name, ok)

    !! record a pass/fail result

    character(len=*),intent(in) :: name
    logical,intent(in) :: ok

    ncheck = ncheck + 1
    if (.not. ok) then
        write(*,'(A)') 'FAIL: '//name
        nfail = nfail + 1
    end if

    end subroutine check

    subroutine check_roots(name, roots, expected, tol)

    !! check that the roots match the expected values

    character(len=*),intent(in) :: name
    real(wp),dimension(:),intent(in) :: roots
    real(wp),dimension(:),intent(in) :: expected
    real(wp),intent(in) :: tol

    if (size(roots)/=size(expected)) then
        call check(name//': number of roots', .false.)
        write(*,*) '   found:   ', roots
        write(*,*) '   expected:', expected
    else
        call check(name//': root values', all(abs(roots-expected)<=tol))
        if (any(abs(roots-expected)>tol)) write(*,*) '   errors: ', roots-expected
    end if

    end subroutine check_roots

    function istr(i) result(s)

    !! integer to string

    integer(ip),intent(in) :: i
    character(len=:),allocatable :: s

    character(len=20) :: tmp

    write(tmp,'(I0)') i
    s = trim(tmp)

    end function istr

    subroutine make_spline(x, y, k, tx, bcoef)

    !! interpolate `y(x)` with the default (not-a-knot) end conditions

    real(wp),dimension(:),intent(in) :: x, y
    integer(ip),intent(in) :: k
    real(wp),dimension(:),allocatable,intent(out) :: tx, bcoef

    integer(ip) :: n, iflag

    n = size(x)
    allocate(tx(n+k), bcoef(n))
    call db1ink(x,n,y,k,0_ip,tx,bcoef,iflag)
    if (iflag/=0) error stop 'db1ink failed'

    end subroutine make_spline

    function linspace(a, b, n) result(x)

    real(wp),intent(in) :: a, b
    integer(ip),intent(in) :: n
    real(wp),dimension(n) :: x

    integer(ip) :: i

    x = [(a + (b-a)*real(i-1,wp)/real(n-1,wp), i=1,n)]

    end function linspace

    function spline_value(tx, bcoef, k, x) result(f)

    real(wp),dimension(:),intent(in) :: tx, bcoef
    integer(ip),intent(in) :: k
    real(wp),intent(in) :: x
    real(wp) :: f

    integer(ip) :: inbvx, iflag
    real(wp),dimension(3*k) :: w0

    inbvx = 1
    call db1val(x,0_ip,tx,size(bcoef,kind=ip),k,bcoef,f,iflag,inbvx,w0)

    end function spline_value

!*****************************************************************************************
!>
!  A cubic polynomial with three known roots, for all orders that reproduce it.

    subroutine test_cubic_roots()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(15) :: x
    integer(ip) :: k, iflag

    x = linspace(-1.0_wp, 2.5_wp, 15_ip)
    do k = 4, 8
        call make_spline(x, cubic(x), k, tx, bcoef)
        call db1root(tx,15_ip,k,bcoef,0.0_wp,roots,iflag)
        call check('cubic roots k='//istr(k)//': iflag', iflag==0)
        call check_roots('cubic roots k='//istr(k), roots, [-0.7_wp, 0.3_wp, 1.9_wp], 1.0e-12_wp)
    end do

    end subroutine test_cubic_roots

    elemental function cubic(x) result(f)

    real(wp),intent(in) :: x
    real(wp) :: f

    f = (x+0.7_wp)*(x-0.3_wp)*(x-1.9_wp)

    end function cubic

!*****************************************************************************************
!>
!  Roots of `s(x) = a` for nonzero `a`: `x**2 = a` on `[-2,2]`.

    subroutine test_level()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(9) :: x
    integer(ip) :: k, iflag

    x = linspace(-2.0_wp, 2.0_wp, 9_ip)
    do k = 3, 6
        call make_spline(x, x**2, k, tx, bcoef)
        call db1root(tx,9_ip,k,bcoef,2.0_wp,roots,iflag)
        call check_roots('x**2=2 k='//istr(k), roots, [-sqrt(2.0_wp), sqrt(2.0_wp)], 1.0e-12_wp)
        call db1root(tx,9_ip,k,bcoef,-1.0_wp,roots,iflag)
        call check('x**2=-1 k='//istr(k)//': no roots', iflag==0 .and. size(roots)==0)
        call db1root(tx,9_ip,k,bcoef,4.0_wp,roots,iflag)
        call check_roots('x**2=4 k='//istr(k)//' (roots at the ends)', roots, [-2.0_wp, 2.0_wp], 1.0e-12_wp)
    end do

    end subroutine test_level

!*****************************************************************************************
!>
!  `sin(x)` on `[0.5, 4*pi-0.5]`: roots near `pi`, `2*pi`, `3*pi`.
!  The roots must be close to the true roots (interpolation error)
!  and must be roots of the spline itself (to rounding).

    subroutine test_sin()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(60) :: x
    integer(ip) :: k, iflag, i

    x = linspace(0.5_wp, 4.0_wp*pi-0.5_wp, 60_ip)
    do k = 3, 8
        call make_spline(x, sin(x), k, tx, bcoef)
        call db1root(tx,60_ip,k,bcoef,0.0_wp,roots,iflag)
        call check_roots('sin k='//istr(k), roots, [pi, 2.0_wp*pi, 3.0_wp*pi], 1.0e-4_wp)
        do i = 1, size(roots)
            call check('sin k='//istr(k)//': s(root)=0', abs(spline_value(tx,bcoef,k,roots(i)))<1.0e-14_wp)
        end do
    end do

    end subroutine test_sin

!*****************************************************************************************
!>
!  Linear spline (`k=2`): the roots are where the polygon crosses zero.

    subroutine test_linear()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp) :: xmin, fmin, xmax, fmax
    integer(ip) :: iflag

    call make_spline([0.0_wp,1.0_wp,2.0_wp,3.0_wp,4.0_wp], [1.0_wp,-1.0_wp,-3.0_wp,1.0_wp,2.0_wp], 2_ip, tx, bcoef)
    call db1root(tx,5_ip,2_ip,bcoef,0.0_wp,roots,iflag)
    call check_roots('linear', roots, [0.5_wp, 2.75_wp], 1.0e-14_wp)

    call db1min(tx,5_ip,2_ip,bcoef,0.0_wp,4.0_wp,xmin,fmin,iflag)
    call check('linear min', iflag==0 .and. xmin==2.0_wp .and. abs(fmin+3.0_wp)<1.0e-14_wp)
    call db1max(tx,5_ip,2_ip,bcoef,0.0_wp,4.0_wp,xmax,fmax,iflag)
    call check('linear max', iflag==0 .and. xmax==4.0_wp .and. abs(fmax-2.0_wp)<1.0e-14_wp)

    ! the polygon touches the level at a corner (a knot)
    call make_spline([0.0_wp,1.0_wp,2.0_wp,3.0_wp], [0.0_wp,0.3_wp,0.1_wp,0.2_wp], 2_ip, tx, bcoef)
    call db1root(tx,4_ip,2_ip,bcoef,0.3_wp,roots,iflag)
    call check_roots('linear corner', roots, [1.0_wp], 1.0e-14_wp)

    end subroutine test_linear

!*****************************************************************************************
!>
!  Double root (the spline touches zero) and triple root (flat crossing).

    subroutine test_multiple_roots()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(11) :: x
    integer(ip) :: k, iflag, i
    real(wp) :: c

    real(wp),dimension(*),parameter :: centers = [0.37_wp, 1.3_wp, 2.11_wp, 3.0_wp, 3.6789_wp]

    x = linspace(0.0_wp, 4.0_wp, 11_ip)
    do k = 4, 6
        ! double root at 1.3, simple root at 3
        call make_spline(x, (x-1.3_wp)**2*(x-3.0_wp), k, tx, bcoef)
        call db1root(tx,11_ip,k,bcoef,0.0_wp,roots,iflag)
        call check_roots('double root k='//istr(k), roots, [1.3_wp, 3.0_wp], 1.0e-12_wp)

        ! triple root at 1.3 (conditioning limits the accuracy to about eps**(1/3))
        call make_spline(x, (x-1.3_wp)**3, k, tx, bcoef)
        call db1root(tx,11_ip,k,bcoef,0.0_wp,roots,iflag)
        call check('triple root k='//istr(k)//': found', size(roots)>=1)
        if (size(roots)>=1) call check('triple root k='//istr(k)//': values', all(abs(roots-1.3_wp)<1.0e-4_wp))
    end do

    ! a touching root at a level a /= 0 must be accurate to about eps
    ! (it is the root of s'), not just sqrt(eps)
    do k = 4, 8
        do i = 1, size(centers)
            c = centers(i)
            call make_spline(x, (x-c)**2*(1.0_wp+0.1_wp*x) + 0.25_wp, k, tx, bcoef)
            call db1root(tx,11_ip,k,bcoef,0.25_wp,roots,iflag)
            call check_roots('touching root k='//istr(k)//' c='//istr(i), roots, [c], 1.0e-12_wp)
        end do
    end do

    end subroutine test_multiple_roots

!*****************************************************************************************
!>
!  Roots exactly at a knot (an interior data point) are returned once,
!  and roots at the ends of the domain are found.

    subroutine test_knot_and_end_roots()

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(9) :: x
    integer(ip) :: k, iflag

    x = linspace(0.0_wp, 8.0_wp, 9_ip)
    do k = 2, 6
        call make_spline(x, x-4.0_wp, k, tx, bcoef)
        call db1root(tx,9_ip,k,bcoef,0.0_wp,roots,iflag)
        call check_roots('root at knot k='//istr(k), roots, [4.0_wp], 1.0e-13_wp)

        call make_spline(x, x*(x-8.0_wp)*(x-3.0_wp)/10.0_wp, k, tx, bcoef)
        call db1root(tx,9_ip,k,bcoef,0.0_wp,roots,iflag)
        if (k>=4) then
            call check_roots('roots at ends k='//istr(k), roots, [0.0_wp, 3.0_wp, 8.0_wp], 1.0e-12_wp)
        else
            ! not reproduced exactly, but the end roots are at data points
            call check('roots at ends k='//istr(k), size(roots)==3)
            if (size(roots)==3) call check('roots at ends k='//istr(k)//': values', &
                abs(roots(1))<1.0e-13_wp .and. abs(roots(3)-8.0_wp)<1.0e-13_wp)
        end if
    end do

    end subroutine test_knot_and_end_roots

!*****************************************************************************************
!>
!  A linear spline equal to 1 on `[1,3]`.

    subroutine test_flat()

    real(wp),dimension(:),allocatable :: roots
    real(wp),dimension(7),parameter :: tx = [0.0_wp,0.0_wp,1.0_wp,2.0_wp,3.0_wp,4.0_wp,4.0_wp]
    real(wp),dimension(5),parameter :: bcoef = [0.0_wp,1.0_wp,1.0_wp,1.0_wp,2.0_wp]
    real(wp),dimension(:),allocatable :: tc, bc
    integer(ip) :: iflag

    call db1root(tx,5_ip,2_ip,bcoef,1.0_wp,roots,iflag)
    call check('flat: iflag=4004', iflag==4004)
    call check('flat: only the ends of the flat segment', all(roots==1.0_wp .or. roots==3.0_wp))
    call db1root(tx,5_ip,2_ip,bcoef,0.5_wp,roots,iflag)
    call check_roots('flat: below', roots, [0.5_wp], 1.0e-14_wp)
    call check('flat: below iflag', iflag==0)
    call db1root(tx,5_ip,2_ip,bcoef,1.5_wp,roots,iflag)
    call check_roots('flat: above', roots, [3.5_wp], 1.0e-14_wp)

    ! constant cubic spline
    call make_spline(linspace(0.0_wp,1.0_wp,6_ip), [(2.0_wp, iflag=1,6)], 4_ip, tc, bc)
    call db1root(tc,6_ip,4_ip,bc,2.0_wp,roots,iflag)
    call check('constant: iflag=4004', iflag==4004)
    call check('constant: only the ends of the domain', all(roots==tc(4) .or. roots==tc(7)))
    call db1root(tc,6_ip,4_ip,bc,1.0_wp,roots,iflag)
    call check('constant: no roots', iflag==0 .and. size(roots)==0)

    end subroutine test_flat

!*****************************************************************************************
!>
!  Extrema of the cubic, which are known in closed form.

    subroutine test_cubic_extrema()

    real(wp),dimension(:),allocatable :: tx, bcoef
    real(wp),dimension(15) :: x
    real(wp) :: xc1, xc2, s1, s2, disc, xmin, fmin, xmax, fmax
    integer(ip) :: k, iflag

    ! critical points: p'(x) = 3x**2 - 2*s1*x + s2 = 0
    s1 = -0.7_wp + 0.3_wp + 1.9_wp
    s2 = -0.7_wp*0.3_wp - 0.7_wp*1.9_wp + 0.3_wp*1.9_wp
    disc = sqrt(4.0_wp*s1**2 - 12.0_wp*s2)
    xc1 = (2.0_wp*s1 - disc) / 6.0_wp  ! local max
    xc2 = (2.0_wp*s1 + disc) / 6.0_wp  ! local min

    x = linspace(-1.0_wp, 2.5_wp, 15_ip)
    do k = 4, 8
        call make_spline(x, cubic(x), k, tx, bcoef)

        ! on [-0.6,1.8] the extrema are the critical points
        call db1min(tx,15_ip,k,bcoef,-0.6_wp,1.8_wp,xmin,fmin,iflag)
        call check('cubic min k='//istr(k), iflag==0 .and. abs(xmin-xc2)<1.0e-10_wp .and. &
                                             abs(fmin-cubic(xc2))<1.0e-13_wp)
        call db1max(tx,15_ip,k,bcoef,-0.6_wp,1.8_wp,xmax,fmax,iflag)
        call check('cubic max k='//istr(k), iflag==0 .and. abs(xmax-xc1)<1.0e-10_wp .and. &
                                             abs(fmax-cubic(xc1))<1.0e-13_wp)

        ! on the whole domain the minimum is still interior, but the maximum is at the end
        call db1min(tx,15_ip,k,bcoef,-1.0_wp,2.5_wp,xmin,fmin,iflag)
        call check('cubic min (domain) k='//istr(k), iflag==0 .and. abs(xmin-xc2)<1.0e-10_wp)
        call db1min(tx,15_ip,k,bcoef,-1.0_wp,0.0_wp,xmin,fmin,iflag)
        call check('cubic min (end) k='//istr(k), iflag==0 .and. xmin==-1.0_wp)
        call db1max(tx,15_ip,k,bcoef,2.5_wp,-1.0_wp,xmax,fmax,iflag)  ! reversed limits
        call check('cubic max (end) k='//istr(k), iflag==0 .and. xmax==2.5_wp)
    end do

    end subroutine test_cubic_extrema

!*****************************************************************************************
!>
!  Extrema of an interpolated `sin(x)`, compared with a dense sampling
!  of the spline itself (the result must be at least as good as every sample).

    subroutine test_sin_extrema()

    real(wp),dimension(:),allocatable :: tx, bcoef
    real(wp),dimension(25) :: x
    real(wp),dimension(20001) :: xs
    real(wp) :: xmin, fmin, xmax, fmax, smin, smax, f
    integer(ip) :: k, iflag, i

    x = linspace(0.0_wp, 2.0_wp*pi, 25_ip)
    xs = linspace(0.0_wp, 2.0_wp*pi, 20001_ip)
    xs(size(xs)) = x(size(x))  ! avoid rounding past the end
    do k = 3, 8
        call make_spline(x, sin(x), k, tx, bcoef)
        call db1min(tx,25_ip,k,bcoef,0.0_wp,2.0_wp*pi,xmin,fmin,iflag)
        call db1max(tx,25_ip,k,bcoef,0.0_wp,2.0_wp*pi,xmax,fmax,iflag)
        smin = huge(1.0_wp)
        smax = -huge(1.0_wp)
        do i = 1, size(xs)
            f = spline_value(tx,bcoef,k,xs(i))
            smin = min(smin,f)
            smax = max(smax,f)
        end do
        call check('sin min k='//istr(k), fmin<=smin+1.0e-15_wp .and. abs(xmin-1.5_wp*pi)<1.0e-3_wp)
        call check('sin max k='//istr(k), fmax>=smax-1.0e-15_wp .and. abs(xmax-0.5_wp*pi)<1.0e-3_wp)
        call check('sin min value k='//istr(k), abs(fmin-spline_value(tx,bcoef,k,xmin))==0.0_wp)
    end do

    end subroutine test_sin_extrema

!*****************************************************************************************
!>
!  Random splines (random knots and coefficients). Every sign change
!  found by a dense search must be found, every root returned must be
!  a root of the spline, and the minimum must not be above any sample.
!  Half of the cases use the spline value at a knot as the level, so
!  there is a root exactly at a knot.

    subroutine test_random()

    integer(ip),parameter :: n = 20
    integer(ip),parameter :: nsamp = 4001

    real(wp),dimension(:),allocatable :: tx, bcoef, roots
    real(wp),dimension(nsamp) :: xs, fs
    real(wp) :: a, xr, xmin, fmin, xmax, fmax, scale, u
    integer(ip) :: k, itrial, iflag, i, l, nmissed, nbad
    integer,dimension(:),allocatable :: seed
    integer :: nseed

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 12345
    call random_seed(put=seed)

    do k = 2, 8
        nmissed = 0
        nbad = 0
        do itrial = 1, 200
            allocate(tx(n+k), bcoef(n))
            ! clamped knots with random interior spacing
            tx(1:k) = 0.0_wp
            call random_number(tx(k+1:n))
            do i = k+1, n
                tx(i) = tx(i-1) + 0.1_wp + tx(i)
            end do
            tx(n+1:n+k) = tx(n) + 0.1_wp + 0.5_wp
            call random_number(bcoef)
            bcoef = 2.0_wp*bcoef - 1.0_wp
            call random_number(u)
            if (mod(itrial,2)==0) then
                l = min(n, k + 1 + int(u*real(n-k,wp)))
                a = spline_value(tx,bcoef,k,tx(l))
            else
                a = 0.5_wp*u - 0.25_wp
            end if
            scale = maxval(abs(bcoef))

            call db1root(tx,n,k,bcoef,a,roots,iflag)
            if (iflag/=0) nbad = nbad + 1
            if (mod(itrial,2)==0) then
                if (.not. any(abs(roots-tx(l))<1.0e-12_wp*tx(l))) nmissed = nmissed + 1
            end if

            ! every root returned is a root
            do i = 1, size(roots)
                if (abs(spline_value(tx,bcoef,k,roots(i))-a)>1.0e-12_wp*scale) nbad = nbad + 1
            end do

            ! every sign change in a dense sampling is found
            xs = linspace(tx(k), tx(n+1), nsamp)
            xs(nsamp) = tx(n+1)  ! avoid rounding past the end
            do i = 1, nsamp
                fs(i) = spline_value(tx,bcoef,k,xs(i)) - a
            end do
            do i = 1, nsamp-1
                if (fs(i)*fs(i+1)<0.0_wp) then
                    if (.not. any(roots>=xs(i) .and. roots<=xs(i+1))) nmissed = nmissed + 1
                end if
            end do

            ! the extrema are at least as good as every sample
            call db1min(tx,n,k,bcoef,tx(k),tx(n+1),xmin,fmin,iflag)
            if (iflag/=0 .or. fmin>minval(fs)+a+1.0e-14_wp) nbad = nbad + 1
            call db1max(tx,n,k,bcoef,tx(k),tx(n+1),xmax,fmax,iflag)
            if (iflag/=0 .or. fmax<maxval(fs)+a-1.0e-14_wp) nbad = nbad + 1

            ! a random interior point is a root of s(x)=s(xr)
            call random_number(xr)
            xr = tx(k) + xr*(tx(n+1)-tx(k))
            call db1root(tx,n,k,bcoef,spline_value(tx,bcoef,k,xr),roots,iflag)
            if (.not. any(abs(roots-xr)<1.0e-9_wp)) nbad = nbad + 1

            deallocate(tx, bcoef)
        end do
        call check('random k='//istr(k)//': no missed sign changes', nmissed==0)
        call check('random k='//istr(k)//': roots and extrema valid', nbad==0)
        if (nmissed>0 .or. nbad>0) write(*,*) '   missed: ', nmissed, ' bad: ', nbad
    end do

    end subroutine test_random

!*****************************************************************************************
!>
!  Invalid inputs return the expected error codes.

    subroutine test_errors()

    real(wp),dimension(:),allocatable :: roots
    real(wp),dimension(7),parameter :: tx = [0.0_wp,0.0_wp,1.0_wp,2.0_wp,3.0_wp,4.0_wp,4.0_wp]
    real(wp),dimension(5),parameter :: bcoef = [0.0_wp,1.0_wp,1.0_wp,1.0_wp,2.0_wp]
    real(wp) :: xmin, fmin
    integer(ip) :: iflag

    call db1root(tx(2:6),5_ip,0_ip,bcoef,0.0_wp,roots,iflag)
    call check('k=0 returns 4001', iflag==4001 .and. size(roots)==0)
    call db1root(tx(2:6),4_ip,1_ip,bcoef(1:4),0.0_wp,roots,iflag)
    call check('k=1 returns 4001', iflag==4001)
    call db1root(tx(1:5),2_ip,3_ip,bcoef(1:2),0.0_wp,roots,iflag)
    call check('n<k returns 4002', iflag==4002)
    call db1min(tx,5_ip,2_ip,bcoef,-1.0_wp,2.0_wp,xmin,fmin,iflag)
    call check('min below domain returns 4003', iflag==4003)
    call db1max(tx,5_ip,2_ip,bcoef,1.0_wp,4.5_wp,xmin,fmin,iflag)
    call check('max above domain returns 4003', iflag==4003)
    call check('status messages', get_status_message(4004_ip)/='Unknown status flag: 4004')

    end subroutine test_errors

!*****************************************************************************************
!>
!  Object-oriented interface.

    subroutine test_oo()

    type(bspline_1d) :: s, s_uninit
    real(wp),dimension(:),allocatable :: roots
    real(wp),dimension(30) :: x
    real(wp) :: xr, xmin, fmin, xmax, fmax, f
    integer(ip) :: iflag

    ! monotonic: exp(x) on [0,2]
    x = linspace(0.0_wp, 2.0_wp, 30_ip)
    call s%initialize(x, exp(x), 4_ip, iflag)

    call s%inverse(exp(1.0_wp), xr, iflag)
    call check('oo inverse: iflag', iflag==0 .and. s%status_ok())
    call check('oo inverse: value', abs(xr-1.0_wp)<1.0e-6_wp)
    call s%evaluate(xr, 0_ip, f, iflag)
    call check('oo inverse: s(x)=y', abs(f-exp(1.0_wp))<1.0e-14_wp)
    call s%inverse(100.0_wp, xr, iflag)
    call check('oo inverse: no solution returns 4005', iflag==4005 .and. .not. s%status_ok())

    call s%minimum(xmin, fmin, iflag)
    call check('oo minimum (domain)', iflag==0 .and. xmin==0.0_wp .and. abs(fmin-1.0_wp)<1.0e-14_wp)
    call s%maximum(xmax, fmax, iflag, xa=0.5_wp, xb=1.5_wp)
    call check('oo maximum (interval)', iflag==0 .and. xmax==1.5_wp)

    ! not monotonic: cos(x) on [0,4*pi]
    x = linspace(0.0_wp, 4.0_wp*pi, 30_ip)
    call s%initialize(x, cos(x), 5_ip, iflag)
    call s%roots(0.0_wp, roots, iflag)
    call check_roots('oo roots', roots, [0.5_wp,1.5_wp,2.5_wp,3.5_wp]*pi, 1.0e-3_wp)
    call s%inverse(0.0_wp, xr, iflag)
    call check('oo inverse: multiple solutions returns 4006', iflag==4006 .and. abs(xr-roots(1))==0.0_wp)
    call s%minimum(xmin, fmin, iflag, xa=2.0_wp, xb=4.0_wp)
    call check('oo minimum (interval)', iflag==0 .and. abs(xmin-pi)<1.0e-2_wp)

    ! not initialized
    call s_uninit%roots(0.0_wp, roots, iflag)
    call check('oo roots: not initialized', iflag==1 .and. size(roots)==0)
    call s_uninit%minimum(xmin, fmin, iflag)
    call check('oo minimum: not initialized', iflag==1)

    end subroutine test_oo

!*****************************************************************************************
!>
!  Cases from the code review of PR #107.

    subroutine test_review_cases()

    real(wp),dimension(:),allocatable :: roots
    real(wp) :: z, xmin, fmin, xmax, fmax
    integer(ip) :: iflag

    ! the left end of a flat region is a root, even if the spline reaches
    ! it with a root of high multiplicity: s(x) = (1-x)**4 on [0,1], 0 on [1,2]
    call db1root([0,0,0,0,0,1,2,2,2,2,2]*1.0_wp, 6_ip, 5_ip, [1,0,0,0,0,0]*1.0_wp, 0.0_wp, roots, iflag)
    call check('flat region after high multiplicity root: iflag', iflag==4004)
    call check_roots('flat region after high multiplicity root', roots, [1.0_wp, 2.0_wp], 0.0_wp)
    ! same, with the flat region in the middle of the domain
    call db1root([0,0,0,0,0,1,2,3,3,3,3,3]*1.0_wp, 7_ip, 5_ip, [1,0,0,0,0,0,1]*1.0_wp, 0.0_wp, roots, iflag)
    call check_roots('flat region in the middle', roots, [1.0_wp, 2.0_wp], 0.0_wp)

    ! discontinuous linear spline: -1 -> 0 on [0,1), 1 -> 2 on [1,2].
    ! s(1) = 1 (the value from the right), so there is no root
    call db1root([0,0,1,1,2,2]*1.0_wp, 4_ip, 2_ip, [-1,0,1,2]*1.0_wp, 0.0_wp, roots, iflag)
    call check('discontinuous: left limit is not a root', iflag==0 .and. size(roots)==0)
    ! but a root of the right piece at the knot is a root
    call db1root([0,0,1,1,2,2]*1.0_wp, 4_ip, 2_ip, [-1,-2,0,1]*1.0_wp, 0.0_wp, roots, iflag)
    call check_roots('discontinuous: right value is a root', roots, [1.0_wp], 0.0_wp)
    ! jump across the level: no root
    call db1root([0,0,1,1,2,2]*1.0_wp, 4_ip, 2_ip, [-2,-1,1,2]*1.0_wp, 0.0_wp, roots, iflag)
    call check('discontinuous: jump across the level is not a root', size(roots)==0)
    ! quadratic with a discontinuity (triple knot): left piece ends at zero
    call db1root([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [-1,-1,0,1,1,2]*1.0_wp, 0.0_wp, roots, iflag)
    call check('discontinuous quadratic: left limit is not a root', size(roots)==0)

    ! narrow domain far from zero: the roots z+2 and z+6 are distinct
    z = 1.0e16_wp
    call db1root([z,z,z+4,z+8,z+8], 3_ip, 2_ip, [1,-1,1]*1.0_wp, 0.0_wp, roots, iflag)
    call check_roots('narrow domain far from zero', roots, [z+2.0_wp, z+6.0_wp], 0.0_wp)
    z = 1.0e6_wp
    call db1root([z,z,z+1.0e-9_wp,z+2.0e-9_wp,z+2.0e-9_wp], 3_ip, 2_ip, [1,-1,1]*1.0_wp, 0.0_wp, roots, iflag)
    call check('narrow domain far from zero (2): two roots', size(roots)==2)

    ! tiny coefficients: 1e-200*(x-0.25)*(x-0.75) on [0,1] (Bernstein form)
    call db1root([0,0,0,1,1,1]*1.0_wp, 3_ip, 3_ip, 1.0e-200_wp*[0.1875_wp,-0.3125_wp,0.1875_wp], &
                 0.0_wp, roots, iflag)
    call check_roots('tiny quadratic', roots, [0.25_wp, 0.75_wp], 1.0e-15_wp)
    ! huge coefficients
    call db1root([0,0,0,1,1,1]*1.0_wp, 3_ip, 3_ip, 1.0e200_wp*[0.1875_wp,-0.3125_wp,0.1875_wp], &
                 0.0_wp, roots, iflag)
    call check_roots('huge quadratic', roots, [0.25_wp, 0.75_wp], 1.0e-15_wp)
    ! extrema of a tiny cubic (quadratic derivative): Bernstein [0,1,-1,0]*1e-200,
    ! minimum -sqrt(3)/6*1e-200 at (3+sqrt(3))/6, maximum at (3-sqrt(3))/6
    call db1min([0,0,0,0,1,1,1,1]*1.0_wp, 4_ip, 4_ip, 1.0e-200_wp*[0,1,-1,0], 0.0_wp, 1.0_wp, xmin, fmin, iflag)
    call check('tiny cubic minimum', iflag==0 .and. abs(xmin-(3.0_wp+sqrt(3.0_wp))/6.0_wp)<1.0e-14_wp .and. &
               abs(fmin/1.0e-200_wp+sqrt(3.0_wp)/6.0_wp)<1.0e-14_wp)
    call db1max([0,0,0,0,1,1,1,1]*1.0_wp, 4_ip, 4_ip, 1.0e-200_wp*[0,1,-1,0], 0.0_wp, 1.0_wp, xmax, fmax, iflag)
    call check('tiny cubic maximum', iflag==0 .and. abs(xmax-(3.0_wp-sqrt(3.0_wp))/6.0_wp)<1.0e-14_wp .and. &
               abs(fmax/1.0e-200_wp-sqrt(3.0_wp)/6.0_wp)<1.0e-14_wp)

    ! double root at 1/3 and simple root at 0.75: no extra root near 1/3
    ! (Bernstein coefficients [-9,13,-16,12] are 108*(x-1/3)**2*(x-3/4) on [0,1])
    call db1root([0,0,0,0,1,1,1,1]*1.0_wp, 4_ip, 4_ip, [-9,13,-16,12]*1.0_wp, 0.0_wp, roots, iflag)
    call check_roots('double and simple root', roots, [1.0_wp/3.0_wp, 0.75_wp], 1.0e-14_wp)

    end subroutine test_review_cases

!*****************************************************************************************
!>
!  Cases from the review of the changes since 7.5.1.

    subroutine test_review_cases_2()

    real(wp),dimension(:),allocatable :: roots
    real(wp) :: xmin, fmin, xmax, fmax, xref, fref
    integer(ip) :: iflag

    ! discontinuous quadratic: (1-x)**2 on [0,1), 2 on [1,2]. The minimum is
    ! the limit from the left at x=1 (0, not attained)
    call db1min([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [1,0,0,2,2,2]*1.0_wp, 0.0_wp, 2.0_wp, xmin, fmin, iflag)
    call check('discontinuous minimum: left limit', iflag==4007 .and. xmin==1.0_wp .and. fmin==0.0_wp)
    call db1max([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [1,0,0,2,2,2]*1.0_wp, 0.0_wp, 2.0_wp, xmax, fmax, iflag)
    call check('discontinuous maximum: attained', iflag==0 .and. fmax==2.0_wp)
    ! the maximum is the limit from the left: 3*x**2 on [0,1), 0 on [1,2]
    call db1max([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [0,0,3,0,0,0]*1.0_wp, 0.0_wp, 2.0_wp, xmax, fmax, iflag)
    call check('discontinuous maximum: left limit', iflag==4007 .and. xmax==1.0_wp .and. fmax==3.0_wp)
    ! the jump is at the right end of [xa,xb]: the limit from the left counts
    call db1min([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [1,0,0,2,2,2]*1.0_wp, 0.5_wp, 1.0_wp, xmin, fmin, iflag)
    call check('discontinuous minimum: jump at xb', iflag==4007 .and. fmin==0.0_wp)
    ! but not at the left end (only the value from the right is in [xa,xb])
    call db1min([0,0,0,1,1,1,2,2,2]*1.0_wp, 6_ip, 3_ip, [1,0,0,2,2,2]*1.0_wp, 1.0_wp, 2.0_wp, xmin, fmin, iflag)
    call check('discontinuous minimum: jump at xa', iflag==0 .and. fmin==2.0_wp)
    call check('4007 message', index(get_status_message(4007_ip),'limit from the left')>0)

    ! large coefficients on a narrow interval: the derivative coefficients
    ! overflowed. The minimum is the same as for the unscaled spline.
    call db1min([0,0,0,0,1,1,1,1]*1.0_wp, 4_ip, 4_ip, [1.0_wp,-0.2_wp,-0.6_wp,1.0_wp], &
                0.0_wp, 1.0_wp, xref, fref, iflag)
    call db1min([0,0,0,0,1,1,1,1]*1.0e-10_wp, 4_ip, 4_ip, [1.0_wp,-0.2_wp,-0.6_wp,1.0_wp]*1.0e300_wp, &
                0.0_wp, 1.0e-10_wp, xmin, fmin, iflag)
    call check('scaled minimum: iflag', iflag==0)
    call check('scaled minimum: location', abs(xmin-xref*1.0e-10_wp)<=1.0e-12_wp*1.0e-10_wp)
    call check('scaled minimum: value', abs(fmin-fref*1.0e300_wp)<=1.0e-12_wp*abs(fref)*1.0e300_wp)

    ! a domain wider than huge(): the duplicate tolerance overflowed
    call db1root([-1,-1,0,1,1]*1.0e308_wp, 3_ip, 2_ip, [1,-1,1]*1.0_wp, 0.0_wp, roots, iflag)
    call check_roots('domain wider than huge()', roots, [-0.5e308_wp, 0.5e308_wp], 1.0e294_wp)

    end subroutine test_review_cases_2

    end program bspline_roots_test
!*****************************************************************************************
