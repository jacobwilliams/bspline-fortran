!*****************************************************************************************
!>
!  Tests for [[greville_abscissae]] and the `get_knots`, `get_coefficients`,
!  `initialize_from_coefficients`, and `control_points` methods of the
!  object-oriented interface (issue #86).
!
!  * Greville abscissae for known knot vectors, and the linear precision
!    property (the spline with the Greville abscissae as coefficients is
!    \( s(x) = x \))
!  * the knots and coefficients returned match the subroutine interface
!  * a spline rebuilt from its knots and coefficients gives the same
!    values (and derivatives) as the original, for 1D to 6D
!  * the control points of a 1D spline, and the convex hull property
!  * error codes for invalid knots or coefficients, and uninitialized classes

    program bspline_control_points_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer :: nfail   !! number of failed checks
    integer :: ncheck  !! number of checks

    nfail = 0
    ncheck = 0

    call test_greville()
    call test_1d()
    call test_2d()
    call test_3d()
    call test_4d()
    call test_5d()
    call test_6d()
    call test_errors()

    write(*,*) ''
    write(*,*) ncheck, ' checks'
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'control points test failed'
    else
        write(*,*) 'all control points tests PASSED'
    end if

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Record a pass/fail result.

    subroutine check(name, ok)

    character(len=*),intent(in) :: name  !! name of the check
    logical,intent(in)          :: ok    !! true if the check passed

    ncheck = ncheck + 1
    if (.not. ok) then
        write(*,'(A)') 'FAIL: '//name
        nfail = nfail + 1
    end if

    end subroutine check
!*****************************************************************************************

!*****************************************************************************************
!>
!  `n` equally spaced points from `a` to `b`.

    pure function linspace(a, b, n) result(x)

    real(wp),intent(in)    :: a  !! first point
    real(wp),intent(in)    :: b  !! last point
    integer(ip),intent(in) :: n  !! number of points
    real(wp),dimension(n)  :: x  !! the points

    integer(ip) :: i  !! counter

    do i = 1, n
        x(i) = a + (b-a)*real(i-1,wp)/real(n-1,wp)
    end do
    x(n) = b

    end function linspace
!*****************************************************************************************

!*****************************************************************************************
!>
!  Greville abscissae.

    subroutine test_greville()

    real(wp),dimension(6) :: xi        !! Greville abscissae (cubic example)
    real(wp),dimension(:),allocatable :: tx     !! knots
    real(wp),dimension(:),allocatable :: bcoef  !! coefficients
    real(wp),dimension(:),allocatable :: g      !! Greville abscissae
    real(wp),dimension(3*10) :: w0     !! work array for [[db1val]]
    real(wp) :: x       !! evaluation point
    real(wp) :: f       !! spline value
    real(wp) :: err     !! maximum error
    integer(ip) :: k     !! spline order
    integer(ip) :: i     !! counter
    integer(ip) :: iflag !! status flag
    integer(ip) :: inbvx !! initialization parameter for [[db1val]]
    integer(ip),parameter :: n = 12  !! number of coefficients

    ! clamped cubic: knots [0,0,0,0,1,2,3,3,3,3]
    call greville_abscissae([0,0,0,0,1,2,3,3,3,3]*1.0_wp, 6_ip, 4_ip, xi)
    call check('greville: clamped cubic', &
               all(abs(xi-[0.0_wp, 1.0_wp/3.0_wp, 1.0_wp, 2.0_wp, 8.0_wp/3.0_wp, 3.0_wp])<1.0e-15_wp))
    ! linear: the interior knots
    call greville_abscissae([0,0,1,2,3,3]*1.0_wp, 4_ip, 2_ip, xi(1:4))
    call check('greville: linear', all(xi(1:4)==[0.0_wp,1.0_wp,2.0_wp,3.0_wp]))
    ! order 1: midpoints of the supports
    call greville_abscissae([0,1,3,6]*1.0_wp, 3_ip, 1_ip, xi(1:3))
    call check('greville: order 1', all(xi(1:3)==[0.5_wp,2.0_wp,4.5_wp]))

    ! linear precision: the spline with the Greville abscissae as its
    ! coefficients is s(x) = x (random interior knots, orders 2 to 10)
    do k = 2, 10
        allocate(tx(n+k), bcoef(n), g(n))
        tx(1:k) = 0.0_wp
        do i = k+1, n
            tx(i) = tx(i-1) + 0.1_wp + 0.37_wp*mod(real(i*k,wp),1.7_wp)
        end do
        tx(n+1:n+k) = tx(n) + 0.5_wp
        call greville_abscissae(tx, n, k, g)
        call check('greville: in the domain k='//istr(k), all(g>=tx(k) .and. g<=tx(n+1)))
        call check('greville: increasing k='//istr(k), all(g(2:)>g(:n-1)))
        bcoef = g
        err = 0.0_wp
        inbvx = 1
        do i = 0, 200
            x = tx(k) + (tx(n+1)-tx(k))*real(i,wp)/200.0_wp
            call db1val(x,0_ip,tx,n,k,bcoef,f,iflag,inbvx,w0)
            err = max(err, abs(f-x))
        end do
        call check('greville: linear precision k='//istr(k), err<1.0e-13_wp)
        deallocate(tx, bcoef, g)
    end do

    end subroutine test_greville
!*****************************************************************************************

!*****************************************************************************************
!>
!  1D: knots and coefficients, rebuilding, and control points.

    subroutine test_1d()

    integer(ip),parameter :: nx = 15  !! number of data points
    integer(ip),parameter :: kx = 4   !! spline order

    type(bspline_1d) :: s1   !! spline from the data
    type(bspline_1d) :: s2   !! spline rebuilt from the knots and coefficients
    real(wp),dimension(nx) :: x       !! data abscissae
    real(wp),dimension(nx) :: fcn     !! data
    real(wp),dimension(nx+kx) :: tx_ref     !! knots from [[db1ink]]
    real(wp),dimension(nx) :: bcoef_ref     !! coefficients from [[db1ink]]
    real(wp),dimension(:),allocatable :: tx     !! knots from `get_knots`
    real(wp),dimension(:),allocatable :: bcoef  !! coefficients from `get_coefficients`
    real(wp),dimension(:),allocatable :: xc     !! control point locations
    real(wp),dimension(:),allocatable :: yc     !! control point values
    real(wp),dimension(:),allocatable :: g      !! Greville abscissae
    real(wp) :: xv    !! evaluation point
    real(wp) :: f1    !! value from `s1`
    real(wp) :: f2    !! value from `s2`
    real(wp) :: err   !! maximum difference
    logical :: hull   !! convex hull property satisfied
    integer(ip) :: i      !! counter
    integer(ip) :: idx    !! derivative
    integer(ip) :: iflag  !! status flag

    x = linspace(0.0_wp, 3.0_wp, nx)
    fcn = sin(2.0_wp*x) + 0.1_wp*x**2

    call db1ink(x,nx,fcn,kx,0_ip,tx_ref,bcoef_ref,iflag)
    call s1%initialize(x,fcn,kx,iflag)
    call check('1d: initialize', iflag==0)

    call s1%get_knots(tx,iflag)
    call check('1d: get_knots', iflag==0 .and. size(tx)==nx+kx)
    if (size(tx)==nx+kx) call check('1d: knots match db1ink', all(tx==tx_ref))
    call s1%get_coefficients(bcoef,iflag)
    call check('1d: get_coefficients', iflag==0 .and. size(bcoef)==nx)
    if (size(bcoef)==nx) call check('1d: coefficients match db1ink', all(bcoef==bcoef_ref))

    ! rebuild from the knots and coefficients
    call s2%initialize_from_coefficients(tx,bcoef,kx,iflag)
    call check('1d: initialize_from_coefficients', iflag==0 .and. s2%status_ok())
    err = 0.0_wp
    do idx = 0, kx-1
        do i = 0, 300
            xv = 3.0_wp*real(i,wp)/300.0_wp
            call s1%evaluate(xv,idx,f1,iflag)
            call s2%evaluate(xv,idx,f2,iflag)
            err = max(err, abs(f1-f2))
        end do
    end do
    call check('1d: rebuilt spline is the same', err==0.0_wp)
    call s1%integral(0.5_wp,2.5_wp,f1,iflag)
    call s2%integral(0.5_wp,2.5_wp,f2,iflag)
    call check('1d: rebuilt spline integral', f1==f2)

    ! extrapolation flag
    call s2%initialize_from_coefficients(tx,bcoef,kx,iflag,extrap=.true.)
    call s2%evaluate(3.5_wp,0_ip,f2,iflag)
    call check('1d: rebuilt spline with extrapolation', iflag==0)
    call s2%initialize_from_coefficients(tx,bcoef,kx,iflag)
    call s2%evaluate(3.5_wp,0_ip,f2,iflag)
    call check('1d: rebuilt spline without extrapolation', iflag/=0)

    ! control points
    call s1%control_points(xc,yc,iflag)
    call check('1d: control_points', iflag==0 .and. size(xc)==nx .and. size(yc)==nx)
    allocate(g(nx))
    call greville_abscissae(tx_ref,nx,kx,g)
    call check('1d: control point locations', all(xc==g))
    call check('1d: control point values', all(yc==bcoef_ref))
    ! convex hull: on each knot interval, s(x) is between the min and
    ! max of the kx control point values that are nonzero there
    hull = .true.
    do i = 0, 300
        xv = 3.0_wp*real(i,wp)/300.0_wp
        call s1%evaluate(xv,0_ip,f1,iflag)
        if (f1<minval(yc)-1.0e-14_wp .or. f1>maxval(yc)+1.0e-14_wp) hull = .false.
    end do
    call check('1d: convex hull', hull)

    end subroutine test_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D: rebuilding from the knots and coefficients.

    subroutine test_2d()

    integer(ip),parameter :: nx = 8, ny = 7       !! number of data points
    integer(ip),parameter :: kx = 4, ky = 3       !! spline orders

    type(bspline_2d) :: s1   !! spline from the data
    type(bspline_2d) :: s2   !! spline rebuilt from the knots and coefficients
    real(wp),dimension(nx) :: x   !! data abscissae
    real(wp),dimension(ny) :: y   !! data abscissae
    real(wp),dimension(nx,ny) :: fcn  !! data
    real(wp),dimension(:),allocatable :: tx, ty      !! knots
    real(wp),dimension(:,:),allocatable :: bcoef     !! coefficients
    real(wp) :: f1, f2, err  !! values and maximum difference
    real(wp),dimension(2) :: v  !! evaluation point
    integer(ip) :: i, j, iflag  !! counters and status flag

    x = linspace(0.0_wp,1.0_wp,nx); y = linspace(-1.0_wp,1.0_wp,ny)
    do j = 1, ny
        do i = 1, nx
            fcn(i,j) = sin(x(i)+2.0_wp*y(j))
        end do
    end do
    call s1%initialize(x,y,fcn,kx,ky,iflag)
    call s1%get_knots(tx,ty,iflag)
    call check('2d: get_knots', iflag==0 .and. size(tx)==nx+kx .and. size(ty)==ny+ky)
    call s1%get_coefficients(bcoef,iflag)
    call check('2d: get_coefficients', iflag==0 .and. all(shape(bcoef)==[nx,ny]))
    call s2%initialize_from_coefficients(tx,ty,bcoef,kx,ky,iflag)
    call check('2d: initialize_from_coefficients', iflag==0)
    err = 0.0_wp
    do i = 1, 200
        v = [real(mod(i*37,101),wp)/100.0_wp, -1.0_wp+2.0_wp*real(mod(i*53,97),wp)/96.0_wp]
        call s1%evaluate(v(1),v(2),0_ip,1_ip,f1,iflag)
        call s2%evaluate(v(1),v(2),0_ip,1_ip,f2,iflag)
        err = max(err, abs(f1-f2))
    end do
    call check('2d: rebuilt spline is the same', err==0.0_wp)

    end subroutine test_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  3D: rebuilding from the knots and coefficients.

    subroutine test_3d()

    integer(ip),parameter :: n = 6, k = 4  !! number of data points and order

    type(bspline_3d) :: s1, s2  !! original and rebuilt splines
    real(wp),dimension(n) :: x  !! data abscissae (all directions)
    real(wp),dimension(n,n,n) :: fcn  !! data
    real(wp),dimension(:),allocatable :: tx, ty, tz  !! knots
    real(wp),dimension(:,:,:),allocatable :: bcoef   !! coefficients
    real(wp),dimension(3) :: v  !! evaluation point
    real(wp) :: f1, f2, err     !! values and maximum difference
    integer(ip) :: i, j, l, iflag  !! counters and status flag

    x = linspace(0.0_wp,1.0_wp,n)
    do l = 1, n; do j = 1, n; do i = 1, n
        fcn(i,j,l) = exp(x(i)) * cos(x(j)) + x(l)**3
    end do; end do; end do
    call s1%initialize(x,x,x,fcn,k,k,k,iflag)
    call s1%get_knots(tx,ty,tz,iflag)
    call check('3d: get_knots', iflag==0 .and. size(tz)==n+k)
    call s1%get_coefficients(bcoef,iflag)
    call check('3d: get_coefficients', iflag==0 .and. all(shape(bcoef)==[n,n,n]))
    call s2%initialize_from_coefficients(tx,ty,tz,bcoef,k,k,k,iflag)
    call check('3d: initialize_from_coefficients', iflag==0)
    err = 0.0_wp
    do i = 1, 200
        v = [real(mod(i*37,101),wp), real(mod(i*53,97),wp), real(mod(i*29,89),wp)] / [100.0_wp, 96.0_wp, 88.0_wp]
        call s1%evaluate(v(1),v(2),v(3),1_ip,0_ip,1_ip,f1,iflag)
        call s2%evaluate(v(1),v(2),v(3),1_ip,0_ip,1_ip,f2,iflag)
        err = max(err, abs(f1-f2))
    end do
    call check('3d: rebuilt spline is the same', err==0.0_wp)

    end subroutine test_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  4D: rebuilding from the knots and coefficients.

    subroutine test_4d()

    integer(ip),parameter :: n = 5, k = 3  !! number of data points and order

    type(bspline_4d) :: s1, s2  !! original and rebuilt splines
    real(wp),dimension(n) :: x  !! data abscissae (all directions)
    real(wp),dimension(n,n,n,n) :: fcn  !! data
    real(wp),dimension(:),allocatable :: tx, ty, tz, tq  !! knots
    real(wp),dimension(:,:,:,:),allocatable :: bcoef     !! coefficients
    real(wp),dimension(4) :: v  !! evaluation point
    real(wp) :: f1, f2, err     !! values and maximum difference
    integer(ip) :: i, j, l, m, iflag  !! counters and status flag

    x = linspace(0.0_wp,1.0_wp,n)
    do m = 1, n; do l = 1, n; do j = 1, n; do i = 1, n
        fcn(i,j,l,m) = sin(x(i)+x(j)) + x(l)*x(m)
    end do; end do; end do; end do
    call s1%initialize(x,x,x,x,fcn,k,k,k,k,iflag)
    call s1%get_knots(tx,ty,tz,tq,iflag)
    call check('4d: get_knots', iflag==0 .and. size(tq)==n+k)
    call s1%get_coefficients(bcoef,iflag)
    call check('4d: get_coefficients', iflag==0 .and. all(shape(bcoef)==[n,n,n,n]))
    call s2%initialize_from_coefficients(tx,ty,tz,tq,bcoef,k,k,k,k,iflag)
    call check('4d: initialize_from_coefficients', iflag==0)
    err = 0.0_wp
    do i = 1, 200
        v = real([mod(i*37,101), mod(i*53,97), mod(i*29,89), mod(i*17,83)],wp) / &
            [100.0_wp, 96.0_wp, 88.0_wp, 82.0_wp]
        call s1%evaluate(v(1),v(2),v(3),v(4),0_ip,1_ip,0_ip,0_ip,f1,iflag)
        call s2%evaluate(v(1),v(2),v(3),v(4),0_ip,1_ip,0_ip,0_ip,f2,iflag)
        err = max(err, abs(f1-f2))
    end do
    call check('4d: rebuilt spline is the same', err==0.0_wp)

    end subroutine test_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  5D: rebuilding from the knots and coefficients.

    subroutine test_5d()

    integer(ip),parameter :: n = 4, k = 3  !! number of data points and order

    type(bspline_5d) :: s1, s2  !! original and rebuilt splines
    real(wp),dimension(n) :: x  !! data abscissae (all directions)
    real(wp),dimension(n,n,n,n,n) :: fcn  !! data
    real(wp),dimension(:),allocatable :: tx, ty, tz, tq, tr  !! knots
    real(wp),dimension(:,:,:,:,:),allocatable :: bcoef       !! coefficients
    real(wp),dimension(5) :: v  !! evaluation point
    real(wp) :: f1, f2, err     !! values and maximum difference
    integer(ip) :: i, j, l, m, o, iflag  !! counters and status flag

    x = linspace(0.0_wp,1.0_wp,n)
    do o = 1, n; do m = 1, n; do l = 1, n; do j = 1, n; do i = 1, n
        fcn(i,j,l,m,o) = cos(x(i)) + x(j)*x(l) - x(m)**2 + x(o)
    end do; end do; end do; end do; end do
    call s1%initialize(x,x,x,x,x,fcn,k,k,k,k,k,iflag)
    call s1%get_knots(tx,ty,tz,tq,tr,iflag)
    call check('5d: get_knots', iflag==0 .and. size(tr)==n+k)
    call s1%get_coefficients(bcoef,iflag)
    call check('5d: get_coefficients', iflag==0 .and. all(shape(bcoef)==[n,n,n,n,n]))
    call s2%initialize_from_coefficients(tx,ty,tz,tq,tr,bcoef,k,k,k,k,k,iflag)
    call check('5d: initialize_from_coefficients', iflag==0)
    err = 0.0_wp
    do i = 1, 200
        v = real([mod(i*37,101), mod(i*53,97), mod(i*29,89), mod(i*17,83), mod(i*13,79)],wp) / &
            [100.0_wp, 96.0_wp, 88.0_wp, 82.0_wp, 78.0_wp]
        call s1%evaluate(v(1),v(2),v(3),v(4),v(5),0_ip,0_ip,1_ip,0_ip,0_ip,f1,iflag)
        call s2%evaluate(v(1),v(2),v(3),v(4),v(5),0_ip,0_ip,1_ip,0_ip,0_ip,f2,iflag)
        err = max(err, abs(f1-f2))
    end do
    call check('5d: rebuilt spline is the same', err==0.0_wp)

    end subroutine test_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  6D: rebuilding from the knots and coefficients.

    subroutine test_6d()

    integer(ip),parameter :: n = 4, k = 2  !! number of data points and order

    type(bspline_6d) :: s1, s2  !! original and rebuilt splines
    real(wp),dimension(n) :: x  !! data abscissae (all directions)
    real(wp),dimension(n,n,n,n,n,n) :: fcn  !! data
    real(wp),dimension(:),allocatable :: tx, ty, tz, tq, tr, ts  !! knots
    real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef         !! coefficients
    real(wp),dimension(6) :: v  !! evaluation point
    real(wp) :: f1, f2, err     !! values and maximum difference
    integer(ip) :: i, j, l, m, o, p, iflag  !! counters and status flag

    x = linspace(0.0_wp,1.0_wp,n)
    do p = 1, n; do o = 1, n; do m = 1, n; do l = 1, n; do j = 1, n; do i = 1, n
        fcn(i,j,l,m,o,p) = x(i) + 2.0_wp*x(j)*x(l) + sin(x(m)+x(o)) - x(p)
    end do; end do; end do; end do; end do; end do
    call s1%initialize(x,x,x,x,x,x,fcn,k,k,k,k,k,k,iflag)
    call s1%get_knots(tx,ty,tz,tq,tr,ts,iflag)
    call check('6d: get_knots', iflag==0 .and. size(ts)==n+k)
    call s1%get_coefficients(bcoef,iflag)
    call check('6d: get_coefficients', iflag==0 .and. all(shape(bcoef)==[n,n,n,n,n,n]))
    call s2%initialize_from_coefficients(tx,ty,tz,tq,tr,ts,bcoef,k,k,k,k,k,k,iflag)
    call check('6d: initialize_from_coefficients', iflag==0)
    err = 0.0_wp
    do i = 1, 200
        v = real([mod(i*37,101), mod(i*53,97), mod(i*29,89), mod(i*17,83), mod(i*13,79), mod(i*11,73)],wp) / &
            [100.0_wp, 96.0_wp, 88.0_wp, 82.0_wp, 78.0_wp, 72.0_wp]
        call s1%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),0_ip,0_ip,0_ip,0_ip,0_ip,0_ip,f1,iflag)
        call s2%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),0_ip,0_ip,0_ip,0_ip,0_ip,0_ip,f2,iflag)
        err = max(err, abs(f1-f2))
    end do
    call check('6d: rebuilt spline is the same', err==0.0_wp)

    end subroutine test_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Error codes.

    subroutine test_errors()

    type(bspline_1d) :: s1  !! 1D spline
    type(bspline_2d) :: s2  !! 2D spline
    type(bspline_6d) :: s6  !! 6D spline (not initialized)
    real(wp),dimension(:),allocatable :: tx, ty, tz, tq, tr, ts  !! knots
    real(wp),dimension(:),allocatable :: bcoef     !! coefficients
    real(wp),dimension(:),allocatable :: xc        !! control point locations
    real(wp),dimension(:,:),allocatable :: bcoef2  !! 2D coefficients
    real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef6  !! 6D coefficients
    real(wp) :: f  !! spline value
    integer(ip) :: iflag  !! status flag

    ! not initialized
    call s1%get_knots(tx,iflag)
    call check('error: get_knots not initialized', iflag==1 .and. size(tx)==0)
    call s1%get_coefficients(bcoef,iflag)
    call check('error: get_coefficients not initialized', iflag==1 .and. size(bcoef)==0)
    call s1%control_points(xc,bcoef,iflag)
    call check('error: control_points not initialized', iflag==1 .and. size(xc)==0)
    call s6%get_knots(tx,ty,tz,tq,tr,ts,iflag)
    call check('error: 6d get_knots not initialized', iflag==1 .and. size(ts)==0)
    call s6%get_coefficients(bcoef6,iflag)
    call check('error: 6d get_coefficients not initialized', iflag==1 .and. size(bcoef6)==0)

    ! invalid knots or coefficients (the class is left uninitialized)
    call s1%initialize_from_coefficients([0,0,1,1]*1.0_wp,[1,2]*1.0_wp,0_ip,iflag)
    call check('error: k < 1', iflag==5001 .and. .not. s1%status_ok())
    call s1%initialize_from_coefficients([0,0,0,1,1,1]*1.0_wp,[1,2]*1.0_wp,4_ip,iflag)
    call check('error: n < k', iflag==5002)
    call s1%initialize_from_coefficients([0,0,1,1,1]*1.0_wp,[1,2]*1.0_wp,2_ip,iflag)
    call check('error: wrong number of knots', iflag==5003)
    call s1%initialize_from_coefficients([0,0,2,1,3,3]*1.0_wp,[1,2,3,4]*1.0_wp,2_ip,iflag)
    call check('error: decreasing knots', iflag==5004)
    call s1%initialize_from_coefficients([1,1,1,1]*1.0_wp,[1,2]*1.0_wp,2_ip,iflag)
    call check('error: empty domain', iflag==5005)
    call check('error: message', index(get_status_message(5004_ip),'non-decreasing')>0)
    call s1%evaluate(0.5_wp,0_ip,f,iflag)
    call check('error: evaluate after failed initialization', iflag/=0)

    ! a valid spline, then an invalid one: the object is destroyed
    call s1%initialize_from_coefficients([0,0,1,1]*1.0_wp,[1,2]*1.0_wp,2_ip,iflag)
    call s1%evaluate(0.5_wp,0_ip,f,iflag)
    call check('valid linear spline', iflag==0 .and. abs(f-1.5_wp)<1.0e-15_wp)
    call s1%initialize_from_coefficients([0,0,1,1,1]*1.0_wp,[1,2]*1.0_wp,2_ip,iflag)
    call s1%get_knots(tx,iflag)
    call check('error: invalid data destroys the previous spline', iflag==1)

    ! the error is found in the second direction
    allocate(bcoef2(2,2))
    bcoef2 = 1.0_wp
    call s2%initialize_from_coefficients([0,0,1,1]*1.0_wp,[0,0,1]*1.0_wp,bcoef2,2_ip,2_ip,iflag)
    call check('error: 2d second direction', iflag==5003)

    end subroutine test_errors
!*****************************************************************************************

!*****************************************************************************************
!>
!  Integer to string.

    pure function istr(i) result(s)

    integer(ip),intent(in)       :: i  !! the integer
    character(len=:),allocatable :: s  !! the string

    character(len=20) :: tmp  !! work string

    write(tmp,'(I0)') i
    s = trim(tmp)

    end function istr
!*****************************************************************************************

    end program bspline_control_points_test
!*****************************************************************************************
