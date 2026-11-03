!*****************************************************************************************
!>
!  Test of [[db1ink]] with specified end derivatives ([[dbint4]]) for
!  spline orders other than 4.
!
!  A spline of order `k` can represent any polynomial of degree `k-1`
!  exactly, so if the data and the end derivatives are taken from such a
!  polynomial, the interpolant must reproduce it (and its derivatives).
!  This is checked for all orders `3 <= k <= 8`, all allowed end
!  derivative orders, and all knot options.

    program dbint4_order_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    real(wp),parameter :: tol = 1.0e-8_wp  !! relative tolerance (all derivatives up to `k-1`)

    integer :: nfail  !! number of failed checks
    integer :: ncase  !! number of cases checked
    integer(ip) :: k, nx, ibcl, ibcr, kntopt, iflag, i, j
    real(wp),dimension(:),allocatable :: c  !! polynomial coefficients
    real(wp),dimension(:),allocatable :: x, y, tx, bcoef
    real(wp) :: u, err
    integer,dimension(:),allocatable :: seed
    integer :: nseed
    logical :: singular

    nfail = 0
    ncase = 0
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 42
    call random_seed(put=seed)

    do k = 3, 8
        allocate(c(k))
        call random_number(c)
        c = 2.0_wp*c - 1.0_wp
        do nx = max(2_ip,k-2_ip), 12
            do kntopt = 1, 3
                do ibcl = 1, k-1
                    do ibcr = 1, k-1

                        allocate(x(nx), y(nx), tx(nx+k+2), bcoef(nx+2))
                        x(1) = -1.0_wp
                        do i = 2, nx
                            call random_number(u)
                            x(i) = x(i-1) + 0.3_wp + u
                        end do
                        do i = 1, nx
                            y(i) = p(x(i),0_ip)
                        end do

                        if (kntopt==3) then
                            call db1ink(x,nx,y,k,ibcl,ibcr,p(x(1),ibcl),p(x(nx),ibcr),&
                                        [(x(1)-0.5_wp*real(k-j,wp), j=1,k-1)],&
                                        [(x(nx)+0.4_wp*real(j,wp), j=1,k-1)],&
                                        tx,bcoef,iflag)
                        else
                            call db1ink(x,nx,y,k,ibcl,ibcr,p(x(1),ibcl),p(x(nx),ibcr),&
                                        kntopt,tx,bcoef,iflag)
                        end if
                        ncase = ncase + 1

                        ! with no interior knots the spline is a single polynomial, so
                        ! constraining its (constant) k-1 derivative at both ends is singular
                        singular = nx==k-2 .and. ibcl==k-1 .and. ibcr==k-1

                        if (singular) then
                            if (iflag/=2007) call fail('singular system not detected')
                        else if (iflag/=0) then
                            call fail('initialization error: '//get_status_message(iflag))
                        else
                            err = max_error()
                            if (err>tol) call fail('polynomial not reproduced')
                        end if

                        deallocate(x, y, tx, bcoef)

                    end do
                end do
            end do
        end do
        deallocate(c)
    end do

    call test_k4_knots()
    call test_errors()

    write(*,*) ''
    write(*,*) ncase, ' polynomial cases checked'
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'dbint4 order test failed'
    else
        write(*,*) 'all dbint4 order tests PASSED'
    end if

    contains
!*****************************************************************************************

    subroutine fail(msg)

    !! report a failure for the current case

    character(len=*),intent(in) :: msg

    nfail = nfail + 1
    write(*,'(A,5(A,I0))') 'FAIL: '//msg, ': k=',k,' nx=',nx,' kntopt=',kntopt,&
                           ' ibcl=',ibcl,' ibcr=',ibcr

    end subroutine fail

    pure function p(xx,d) result(f)

    !! `d`-th derivative of the test polynomial `sum(c(j)*x**(j-1))`

    real(wp),intent(in)    :: xx
    integer(ip),intent(in) :: d
    real(wp)               :: f

    integer(ip) :: jj, m
    real(wp) :: fac

    f = 0.0_wp
    do jj = d+1, size(c)
        fac = 1.0_wp
        do m = 0, d-1
            fac = fac*real(jj-1-m,wp)
        end do
        f = f + c(jj)*fac*xx**(jj-1-d)
    end do

    end function p

    function max_error() result(errmax)

    !! max error of the spline (and all its derivatives) on `[x(1),x(nx)]`,
    !! relative to the largest magnitude of that derivative on the interval

    real(wp) :: errmax

    integer(ip),parameter :: npts = 41  !! number of points to check

    integer(ip) :: ii, id, inbvx, iflag2
    real(wp),dimension(npts) :: xv
    real(wp),dimension(npts,0:k-1) :: f, fexact
    real(wp),dimension(3*k) :: w0

    errmax = huge(1.0_wp)
    do ii = 1, npts
        xv(ii) = x(1) + (x(nx)-x(1))*real(ii-1,wp)/real(npts-1,wp)
    end do
    xv(npts) = x(nx)

    inbvx = 1
    do ii = 1, npts
        do id = 0, k-1
            call db1val(xv(ii),id,tx,nx,nx+2,k,bcoef,f(ii,id),iflag2,inbvx,w0)
            if (iflag2/=0) return
            fexact(ii,id) = p(xv(ii),id)
        end do
    end do

    errmax = 0.0_wp
    do id = 0, k-1
        errmax = max(errmax, maxval(abs(f(:,id)-fexact(:,id))) / max(1.0_wp,maxval(abs(fexact(:,id)))))
    end do

    end function max_error

!*****************************************************************************************
!>
!  For `k=4` the interior knots must be the interior data points
!  (as in the original cubic-only version of [[dbint4]]).

    subroutine test_k4_knots()

    real(wp),dimension(6) :: xx, yy
    real(wp),dimension(6+4+2) :: t
    real(wp),dimension(6+2) :: b
    integer(ip) :: iflag2

    xx = [1.0_wp, 2.0_wp, 3.5_wp, 4.0_wp, 6.0_wp, 7.0_wp]
    yy = sin(xx)
    call db1ink(xx,6_ip,yy,4_ip,2_ip,2_ip,0.0_wp,0.0_wp,1_ip,t,b,iflag2)
    ncase = ncase + 1
    if (iflag2/=0 .or. any(t(5:8)/=xx(2:5)) .or. any(t(1:4)/=xx(1)) .or. any(t(9:12)/=xx(6))) then
        nfail = nfail + 1
        write(*,'(A)') 'FAIL: k=4 knots are not at the data points'
    end if

    end subroutine test_k4_knots

!*****************************************************************************************
!>
!  Invalid inputs return the expected error codes.

    subroutine test_errors()

    real(wp),dimension(6) :: xx, yy
    real(wp),dimension(6+5+2) :: t
    real(wp),dimension(6+2) :: b
    integer(ip) :: iflag2

    xx = [1.0_wp, 2.0_wp, 3.5_wp, 4.0_wp, 6.0_wp, 7.0_wp]
    yy = sin(xx)

    ! k=2 is not supported
    call db1ink(xx,6_ip,yy,2_ip,1_ip,1_ip,0.0_wp,0.0_wp,1_ip,t(1:10),b,iflag2)
    call check_flag('k=2', iflag2, 806_ip)

    ! derivative order must be < k
    call db1ink(xx,6_ip,yy,5_ip,5_ip,1_ip,0.0_wp,0.0_wp,1_ip,t,b,iflag2)
    call check_flag('ibcl=k', iflag2, 2003_ip)
    call db1ink(xx,6_ip,yy,5_ip,1_ip,0_ip,0.0_wp,0.0_wp,1_ip,t,b,iflag2)
    call check_flag('ibcr=0', iflag2, 2004_ip)

    ! tleft and tright must have size k-1
    call db1ink(xx,6_ip,yy,5_ip,1_ip,1_ip,0.0_wp,0.0_wp,[0.0_wp,0.0_wp,0.0_wp],&
                [8.0_wp,8.0_wp,8.0_wp],t,b,iflag2)
    call check_flag('tleft size', iflag2, 2008_ip)

    end subroutine test_errors

    subroutine check_flag(name, iflag2, expected)

    character(len=*),intent(in) :: name
    integer(ip),intent(in) :: iflag2, expected

    ncase = ncase + 1
    if (iflag2/=expected) then
        nfail = nfail + 1
        write(*,'(A,2(A,I0))') 'FAIL: '//name, ': iflag=', iflag2, ' expected ', expected
    end if

    end subroutine check_flag

    end program dbint4_order_test
!*****************************************************************************************
