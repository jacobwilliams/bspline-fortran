!*****************************************************************************************
!>
!  Linear interpolation test

    program bspline_linear_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip
    use pyplot_module
    use linear_interpolation_module

    implicit none

    integer(ip),parameter :: nx    = 10    !! number of points in x
    integer(ip),parameter :: nxv   = 111 !100   !! number of points to evaluate interpolant

    integer(ip),parameter :: kx    = 2    !! order in x [linear]
    !integer(ip),parameter :: kx    = bspline_order_cubic    !! order in x [cubic]
    !integer(ip),parameter :: kx    = bspline_order_quartic
    integer(ip),parameter :: iknot = 0    !! automatically select the knots

    real(wp) :: x(nx)
    real(wp) :: xval(nxv),fval(nxv),fval_linear(nxv)
    real(wp) :: tx(nx+kx)
    real(wp) :: fcn_1d(nx)
    real(wp) :: val,tru,err,errmax
    integer(ip) :: i,idx,iflag,inbvx,iloy
    type(pyplot) :: plt
    integer :: istat  !! pyplot-fortran status flag
    real(wp),dimension(3*kx) :: w1_1d !! work array
    type(linear_interp_1d) :: s1

    idx = 0

    x = huge(1.0_wp)
    xval = huge(1.0_wp)

    do i=1,nx
        x(i) = real(i-1,wp) * 10.0_wp
        fcn_1d(i) = f1(x(i))
    end do
    do i=1,nxv
        xval(i) = real(i-1,wp) - 10
    end do

    !have to set these before the first evaluate call:
    inbvx = 1
    iloy  = 1

    ! initialize
    call db1ink(x,nx,fcn_1d,kx,iknot,tx,fcn_1d,iflag)
    if (iflag/=0) then
        write(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
    end if
    call s1%initialize(x,fcn_1d,iflag)
    if (iflag/=0) then
        write(*,*) 'Error initializing 1D linear interpolator: ', iflag
    end if

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x (deg)',ylabel='f(x)',&
                        title='Linear Test',legend=.true.)
    call plt%add_plot(x,fcn_1d,&
                        label='Function $f(x) = \sin(x)$',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)

    errmax = 0.0_wp
    do i=1,nxv
        call db1val(xval(i),idx,tx,nx,kx,fcn_1d,val,iflag,inbvx,w1_1d,extrap=.true.)
        fval(i) = val  ! save it for plot
        if (iflag/=0) error stop 'error calling db1val: '//get_status_message(iflag)

        tru    = f1(xval(i))
        err    = abs(tru-val)
        errmax = max(err,errmax)
        !write(*,*) xval(i), val, tru, err, iflag

        ! also compute linear interpolation:
        call s1%evaluate(xval(i),val)
        fval_linear(i) = val ! linear vector for plot
        write(*,*) "error : ", xval(i), fval(i) - fval_linear(i)

    end do

    ! check max error against tolerance
    write(*,*) ''
    write(*,*) '1D: max error:', errmax
    write(*,*) ''

    call plt%add_plot(xval,fval,&
            label='k=2 Spline Interpolated',&
            linestyle='g.-',linewidth=1,istat=istat)
    call plt%add_plot(xval,fval_linear,&
            label='Linear Interpolated',&
            linestyle='r.-',linewidth=1,istat=istat)
    call plt%savefig('bspline_linear_test.png',istat=istat)

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        real(wp),parameter :: a = acos(-1.0_wp)/18.0_wp
        f1 = sin(a*x)
        end function f1

    end program bspline_linear_test
