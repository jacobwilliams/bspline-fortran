!*****************************************************************************************
!>
!  Checks [[db1root]], [[db1min]], and [[db1max]] against an independent
!  polynomial root finder: the Jenkins-Traub `rpoly` routine from
!  [polyroots-fortran](https://github.com/jacobwilliams/polyroots-fortran).
!
!  For random splines, each polynomial piece is expanded about the
!  midpoint of its knot interval (from the derivatives computed by
!  [[db1val]]), and all of its complex roots are computed with `rpoly`
!  (and polished with `newton_root_polish`).
!  Then:
!
!  * every simple, well separated real root of a piece inside its
!    interval must be found by [[db1root]];
!  * every root found by [[db1root]] must be close to a real (or nearly
!    real) root of a piece, must satisfy \( s(x) \approx a \) (at a
!    discontinuity, \( s(x) \) is the value from the right), and must
!    not be returned twice;
!  * when the level `a` is the value of the spline at a local extremum
!    (a root of even multiplicity, which is not a sign change), the
!    extremum must be found, and no other root may be returned near it;
!  * the minimum and maximum from [[db1min]] and [[db1max]] must match the
!    best of the end points, the interior knots, and the real roots of the
!    derivative pieces computed by `rpoly`.
!
!  Some interior knots are repeated, up to multiplicity `k` (where the
!  spline is discontinuous). The level `a` is a random value, the value
!  at a random interior knot, or the value at a random local extremum.

    program bspline_roots_polyroots_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip
    use polyroots_module, only: rpoly, newton_root_polish, polyroots_module_rk

    implicit none

    integer(ip),parameter :: ntrial = 300  !! number of random splines for each order
    integer(ip),parameter :: kmax   = 9    !! maximum spline order to test

    real(wp),parameter :: sep_tol  = 1.0e-3_wp   !! a reference root closer than this (in the scaled
                                                 !! variable \( w \in [-1,1] \)) to another root of the
                                                 !! same piece is not required to be found
    real(wp),parameter :: end_tol  = 1.0e-6_wp   !! a reference root closer than this (in \( w \)) to the
                                                 !! end of its interval is not required to be found
    real(wp),parameter :: real_tol = 1.0e-10_wp  !! a reference root with \( |\mathrm{Im}(w)| \) below this
                                                 !! is real (and must be found)
    real(wp),parameter :: near_tol = 1.0e-6_wp   !! a reference root with \( |\mathrm{Im}(w)| \) below this
                                                 !! may be real (and may be returned)
    real(wp),parameter :: x_tol    = 1.0e-10_wp  !! relative tolerance (scaled by the domain width)
                                                 !! for matching a root
    real(wp),parameter :: f_tol    = 1.0e-13_wp  !! tolerance for matching an extremum value
    real(wp),parameter :: res_tol  = 1.0e-10_wp  !! tolerance on \( |s(x)-a| \) for a root found
                                                 !! (the coefficients are in \( [-1,1] \))

    integer :: nfail   !! number of failed checks
    integer :: ncheck  !! number of checks
    integer,dimension(:),allocatable :: seed  !! random number seed
    integer :: nseed   !! size of `seed`
    integer(ip) :: k   !! spline order

    if (polyroots_module_rk/=wp) error stop 'polyroots and bspline real kinds do not match'

    nfail = 0
    ncheck = 0

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 24680
    call random_seed(put=seed)

    do k = 2_ip, kmax
        call test_order(k)
    end do

    write(*,*) ''
    write(*,*) ncheck, ' checks'
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'polyroots comparison test failed'
    else
        write(*,*) 'all polyroots comparison tests PASSED'
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
!  Run the random spline comparisons for one spline order.

    subroutine test_order(k)

    integer(ip),intent(in) :: k  !! spline order

    real(wp),dimension(:),allocatable :: tx      !! knot array
    real(wp),dimension(:),allocatable :: bcoef   !! b-spline coefficients
    real(wp),dimension(:),allocatable :: roots   !! roots from [[db1root]]
    real(wp),dimension(:),allocatable :: xreal   !! real roots of the pieces that must be found
    real(wp),dimension(:),allocatable :: xnear   !! real or nearly real roots of the pieces
    real(wp),dimension(:),allocatable :: xcrit   !! real roots of the derivative pieces
    complex(wp),dimension(:),allocatable :: z    !! roots of a piece (in \( w \))
    real(wp) :: a          !! level for \( s(x) = a \)
    real(wp) :: u          !! random number
    real(wp) :: span       !! width of the domain
    real(wp) :: xt         !! location of a local extremum used as the level (or a knot)
    real(wp) :: xmin       !! location of the minimum from [[db1min]]
    real(wp) :: fmin       !! minimum from [[db1min]]
    real(wp) :: xmax       !! location of the maximum from [[db1max]]
    real(wp) :: fmax       !! maximum from [[db1max]]
    real(wp) :: fref_min   !! reference minimum
    real(wp) :: fref_max   !! reference maximum
    real(wp) :: f          !! spline value
    integer(ip) :: n       !! number of b-spline coefficients
    integer(ip) :: itrial  !! trial counter
    integer(ip) :: iflag   !! status flag
    integer(ip) :: i       !! counter
    integer(ip) :: l       !! knot interval index
    integer(ip) :: mode    !! how the level is chosen: 0 = random, 1 = value at a knot,
                           !! 2 = value at a local extremum
    integer(ip) :: nmissed !! number of reference roots not found
    integer(ip) :: nextra  !! number of roots found with no reference root nearby
    integer(ip) :: nknot   !! number of knot roots not found
    integer(ip) :: ntouch  !! number of local extremum (even multiplicity) roots not found
    integer(ip) :: nntouch !! number of local extremum cases tested
    integer(ip) :: nref    !! number of reference roots compared
    integer(ip) :: next    !! number of extrema that do not match the reference
    integer(ip) :: nres    !! number of roots found with \( |s(x)-a| \) above `res_tol`
    integer(ip) :: ndup    !! number of pairs of roots found that are duplicates

    nmissed = 0; nextra = 0; nknot = 0; ntouch = 0; nntouch = 0; nref = 0; next = 0
    nres = 0; ndup = 0

    do itrial = 1_ip, ntrial

        ! random clamped knots (some interior knots repeated, up to
        ! multiplicity k, where the spline is discontinuous), random coefficients
        call random_number(u)
        n = k + int(u*16.0_wp, ip)
        allocate(tx(n+k), bcoef(n))
        tx(1:k) = 0.0_wp
        do i = k+1_ip, n
            call random_number(u)
            if (i>k+1_ip .and. u<0.12_wp .and. count(tx(k+1_ip:i-1_ip)==tx(i-1_ip))<k) then
                tx(i) = tx(i-1_ip)  ! repeated knot (lower continuity)
            else
                tx(i) = tx(i-1_ip) + 0.05_wp + u
            end if
        end do
        tx(n+1_ip:n+k) = tx(n) + 0.05_wp + 0.5_wp
        call random_number(bcoef)
        bcoef = 2.0_wp*bcoef - 1.0_wp
        span = tx(n+1_ip) - tx(k)

        ! the level
        mode = mod(itrial, 3_ip)
        if (mode==1_ip .and. n==k) mode = 0_ip  ! no interior knots
        if (mode==2_ip .and. k==2_ip) mode = 0_ip  ! no smooth extrema
        select case (mode)
        case(0_ip)
            call random_number(u)
            a = u - 0.5_wp
        case(1_ip)
            call random_number(u)
            l = min(n, k + 1_ip + int(u*real(n-k,wp), ip))
            xt = tx(l)
            a = spline_value(tx,bcoef,k,xt)
        case(2_ip)
            call derivative_roots(tx,bcoef,n,k,xcrit,interior=.true.)
            if (size(xcrit)==0) then
                mode = 0_ip
                a = 0.0_wp
            else
                call random_number(u)
                xt = xcrit(1_ip + min(size(xcrit,kind=ip)-1_ip, int(u*size(xcrit), ip)))
                a = spline_value(tx,bcoef,k,xt)
            end if
        end select

        call db1root(tx,n,k,bcoef,a,roots,iflag)
        call check('k='//istr(k)//' db1root iflag', iflag==0_ip)

        ! reference roots of each piece
        allocate(xreal(0), xnear(0))
        do l = k, n
            if (tx(l)>=tx(l+1_ip)) cycle
            call piece_roots(tx,bcoef,n,k,l,0_ip,a,z)
            do i = 1_ip, size(z,kind=ip)
                if (abs(aimag(z(i)))<=near_tol .and. abs(real(z(i)))<=1.0_wp+near_tol) &
                    xnear = [xnear, to_x(tx,l,real(z(i)))]
                if (abs(aimag(z(i)))<=real_tol .and. abs(real(z(i)))<=1.0_wp-end_tol) then
                    if (count(abs(z-z(i))<sep_tol)==1) xreal = [xreal, to_x(tx,l,real(z(i)))]
                end if
            end do
        end do

        ! every real reference root is found
        do i = 1_ip, size(xreal,kind=ip)
            nref = nref + 1_ip
            if (.not. any(abs(roots-xreal(i))<=x_tol*span)) then
                nmissed = nmissed + 1_ip
                write(*,*) '   k=',k,' trial=',itrial,' missed reference root: ', xreal(i)
            end if
        end do

        ! every root found is near a reference root, and is a root
        ! (at a discontinuity, s(x) is the value from the right)
        do i = 1_ip, size(roots,kind=ip)
            if (.not. any(abs(xnear-roots(i))<=sqrt(x_tol)*span)) then
                nextra = nextra + 1_ip
                write(*,*) '   k=',k,' trial=',itrial,' root with no reference: ', roots(i)
            end if
            if (abs(spline_value(tx,bcoef,k,roots(i))-a)>res_tol) then
                nres = nres + 1_ip
                write(*,*) '   k=',k,' trial=',itrial,' not a root: ', roots(i), &
                           ' s(x)-a = ', spline_value(tx,bcoef,k,roots(i))-a
            end if
        end do

        ! no root is returned twice
        do i = 1_ip, size(roots,kind=ip)-1_ip
            if (roots(i+1_ip)-roots(i)<=1.0e-9_wp*span) then
                ndup = ndup + 1_ip
                write(*,*) '   k=',k,' trial=',itrial,' duplicate roots: ', roots(i), roots(i+1_ip)
            end if
        end do

        ! the knot or local extremum used to set the level is found
        select case (mode)
        case(1_ip)
            if (.not. any(abs(roots-xt)<=x_tol*span)) nknot = nknot + 1_ip
        case(2_ip)
            nntouch = nntouch + 1_ip
            ! (exactly once: no spurious roots in the noise band around it)
            if (.not. any(abs(roots-xt)<=x_tol*span) .or. count(abs(roots-xt)<=1.0e-5_wp*span)/=1) then
                ntouch = ntouch + 1_ip
                write(*,*) '   k=',k,' trial=',itrial,' missed extremum root: ', xt, &
                           ' roots nearby: ', pack(roots, abs(roots-xt)<=1.0e-5_wp*span)
            end if
        end select

        ! extrema: the best of the ends, the interior knots, and the critical points
        call db1min(tx,n,k,bcoef,tx(k),tx(n+1_ip),xmin,fmin,iflag)
        call check('k='//istr(k)//' db1min iflag', iflag==0_ip)
        call db1max(tx,n,k,bcoef,tx(k),tx(n+1_ip),xmax,fmax,iflag)
        call check('k='//istr(k)//' db1max iflag', iflag==0_ip)
        call derivative_roots(tx,bcoef,n,k,xcrit,interior=.false.)
        xcrit = [tx(k:n+1_ip), xcrit]
        fref_min = huge(1.0_wp)
        fref_max = -huge(1.0_wp)
        do i = 1_ip, size(xcrit,kind=ip)
            f = spline_value(tx,bcoef,k,xcrit(i))
            fref_min = min(fref_min, f)
            fref_max = max(fref_max, f)
        end do
        if (abs(fmin-fref_min)>f_tol .or. abs(fmax-fref_max)>f_tol .or. &
            abs(spline_value(tx,bcoef,k,xmin)-fmin)>f_tol .or. &
            abs(spline_value(tx,bcoef,k,xmax)-fmax)>f_tol) then
            next = next + 1_ip
            write(*,*) '   k=',k,' trial=',itrial,' extrema: ', fmin-fref_min, fmax-fref_max
        end if

        deallocate(tx, bcoef, xreal, xnear)
    end do

    call check('k='//istr(k)//': all reference roots found', nmissed==0_ip)
    call check('k='//istr(k)//': no roots without a reference root', nextra==0_ip)
    call check('k='//istr(k)//': all roots found are roots', nres==0_ip)
    call check('k='//istr(k)//': no duplicate roots', ndup==0_ip)
    call check('k='//istr(k)//': knot roots found', nknot==0_ip)
    call check('k='//istr(k)//': local extremum roots found', ntouch==0_ip)
    call check('k='//istr(k)//': extrema match reference', next==0_ip)
    write(*,'(A,I2,A,I6,A,I5,A)') ' k=', k, ': ', nref, ' reference roots, ', &
                                  nntouch, ' local extremum levels'

    end subroutine test_order
!*****************************************************************************************

!*****************************************************************************************
!>
!  All the complex roots of \( s^{(d)}(x) = a \) on knot interval `l`,
!  in the scaled variable \( w \in [-1,1] \), where
!  \( x = x_m + \frac{h}{2} w \), \( x_m \) is the midpoint of the
!  interval, and \( h \) is its width. Computed with `rpoly`, and
!  polished with `newton_root_polish`.

    subroutine piece_roots(tx,bcoef,n,k,l,d,a,z)

    integer(ip),intent(in)                 :: n      !! number of b-spline coefficients
    integer(ip),intent(in)                 :: k      !! spline order
    real(wp),dimension(n+k),intent(in)     :: tx     !! knot array
    real(wp),dimension(n),intent(in)       :: bcoef  !! b-spline coefficients
    integer(ip),intent(in)                 :: l      !! knot interval index: `tx(l) <= x < tx(l+1)`
    integer(ip),intent(in)                 :: d      !! derivative of the spline to use
    real(wp),intent(in)                    :: a      !! level
    complex(wp),dimension(:),allocatable,intent(out) :: z  !! the roots (in \( w \))

    real(wp),dimension(0:k-1) :: q    !! power basis coefficients in \( w \)
    real(wp),dimension(:),allocatable :: op  !! coefficients in order of decreasing powers
                                             !! (as required by `rpoly`)
    real(wp),dimension(:),allocatable :: zr  !! real parts of the roots
    real(wp),dimension(:),allocatable :: zi  !! imaginary parts of the roots
    real(wp),dimension(3*k) :: work   !! work array for [[db1val]]
    real(wp) :: xm      !! midpoint of the interval
    real(wp) :: hh      !! half width of the interval
    real(wp) :: fac     !! \( (h/2)^j / j! \)
    real(wp) :: val     !! derivative of the spline
    integer(ip) :: j    !! counter
    integer(ip) :: m    !! degree of the piece
    integer(ip) :: inbv !! initialization parameter for [[db1val]]
    integer(ip) :: iflag  !! status flag from [[db1val]]
    integer :: istat    !! status flag from `rpoly` and `newton_root_polish`
                        !! (not used for `newton_root_polish`, which returns
                        !! the best point found)

    xm = 0.5_wp*(tx(l)+tx(l+1_ip))
    hh = 0.5_wp*(tx(l+1_ip)-tx(l))
    m = k - 1_ip - d
    inbv = 1_ip
    fac = 1.0_wp
    do j = 0_ip, m
        call db1val(xm,d+j,tx,n,k,bcoef,val,iflag,inbv,work)
        if (iflag/=0_ip) error stop 'db1val error'
        if (j>0_ip) fac = fac*hh/real(j,wp)
        q(j) = val*fac
    end do
    q(0) = q(0) - a

    do while (m>0_ip)
        if (q(m)/=0.0_wp) exit  ! leading coefficient must be nonzero for rpoly
        m = m - 1_ip
    end do
    if (m==0_ip) then
        allocate(z(0))  ! constant piece
        return
    end if

    allocate(zr(m), zi(m))
    op = q(m:0:-1)
    call rpoly(op, int(m), zr, zi, istat)
    if (istat/=0) error stop 'rpoly error'
    ! rpoly roots can lose some accuracy to deflation:
    ! polish each one with Newton's method on the original polynomial
    do j = 1_ip, m
        call newton_root_polish(int(m), op, zr(j), zi(j), &
                                0.0_wp, epsilon(1.0_wp), 20, istat)
    end do
    z = cmplx(zr, zi, kind=wp)

    end subroutine piece_roots
!*****************************************************************************************

!*****************************************************************************************
!>
!  The real roots of \( s'(x) \) inside the knot intervals, from `rpoly`.

    subroutine derivative_roots(tx,bcoef,n,k,xcrit,interior)

    integer(ip),intent(in)             :: n         !! number of b-spline coefficients
    integer(ip),intent(in)             :: k         !! spline order
    real(wp),dimension(n+k),intent(in) :: tx        !! knot array
    real(wp),dimension(n),intent(in)   :: bcoef     !! b-spline coefficients
    real(wp),dimension(:),allocatable,intent(out) :: xcrit  !! the critical points
    logical,intent(in)                 :: interior  !! only return well separated simple roots
                                                    !! away from the ends of their intervals
                                                    !! (local extrema, not inflection points)

    complex(wp),dimension(:),allocatable :: z  !! roots of a derivative piece (in \( w \))
    integer(ip) :: l  !! knot interval index
    integer(ip) :: i  !! counter

    allocate(xcrit(0))
    if (k<3_ip) return
    do l = k, n
        if (tx(l)>=tx(l+1_ip)) cycle
        call piece_roots(tx,bcoef,n,k,l,1_ip,0.0_wp,z)
        do i = 1_ip, size(z,kind=ip)
            if (interior) then
                if (abs(aimag(z(i)))<=real_tol .and. abs(real(z(i)))<=0.9_wp .and. &
                    count(abs(z-z(i))<0.1_wp)==1) xcrit = [xcrit, to_x(tx,l,real(z(i)))]
            else
                if (abs(aimag(z(i)))<=near_tol .and. abs(real(z(i)))<=1.0_wp) &
                    xcrit = [xcrit, to_x(tx,l,real(z(i)))]
            end if
        end do
    end do

    end subroutine derivative_roots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert the scaled variable \( w \in [-1,1] \) on knot interval `l` to \( x \).

    pure function to_x(tx,l,w) result(x)

    real(wp),dimension(:),intent(in) :: tx  !! knot array
    integer(ip),intent(in)           :: l   !! knot interval index
    real(wp),intent(in)              :: w   !! scaled variable
    real(wp)                         :: x   !! \( x = x_m + \frac{h}{2} w \)

    x = 0.5_wp*(tx(l)+tx(l+1_ip)) + 0.5_wp*(tx(l+1_ip)-tx(l))*w

    end function to_x
!*****************************************************************************************

!*****************************************************************************************
!>
!  Value of the spline at `x`.

    function spline_value(tx,bcoef,k,x) result(f)

    real(wp),dimension(:),intent(in) :: tx     !! knot array
    real(wp),dimension(:),intent(in) :: bcoef  !! b-spline coefficients
    integer(ip),intent(in)           :: k      !! spline order
    real(wp),intent(in)              :: x      !! the point
    real(wp)                         :: f      !! \( s(x) \)

    integer(ip) :: inbvx  !! initialization parameter for [[db1val]]
    integer(ip) :: iflag  !! status flag from [[db1val]]
    real(wp),dimension(3*k) :: w0  !! work array for [[db1val]]

    inbvx = 1_ip
    call db1val(x,0_ip,tx,size(bcoef,kind=ip),k,bcoef,f,iflag,inbvx,w0)
    if (iflag/=0_ip) error stop 'db1val error'

    end function spline_value
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

    end program bspline_roots_polyroots_test
!*****************************************************************************************
