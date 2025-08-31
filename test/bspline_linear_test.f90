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
    integer(ip),parameter :: nxv   = 111   !! number of points to evaluate interpolant
    integer(ip),parameter :: iknot = 0    !! automatically select the knots
    integer,parameter :: idx = 0 !! compute value for the spline interpolation

    real(wp) :: x(nx)
    real(wp) :: xval(nxv),fval(nxv),fval_linear(nxv),fval3(nxv),fval4(nxv)
    real(wp) :: fcn_1d(nx)
    real(wp) :: val,err,errmax
    integer(ip) :: i,iflag,inbvx,iloy
    integer :: istat  !! pyplot-fortran status flag
    type(pyplot) :: plt
    type(bspline_1d) :: b2, b3, b4
    type(linear_interp_1d) :: s1

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
    call b2%initialize(x,fcn_1d,bspline_order_linear,iflag,extrap=.true.) ! linear
    if (iflag/=0) error stop 'Error initializing 1D linear spline: '//get_status_message(iflag)
    call b3%initialize(x,fcn_1d,bspline_order_quadratic,iflag,extrap=.true.) ! quadratic
    if (iflag/=0) error stop 'Error initializing 1D quadratic spline: '//get_status_message(iflag)
    call b4%initialize(x,fcn_1d,bspline_order_cubic,iflag,extrap=.true.) ! cubic
    if (iflag/=0) error stop 'Error initializing 1D cubic spline: '//get_status_message(iflag)
    call s1%initialize(x,fcn_1d,iflag)
    if (iflag/=0) error stop 'Error initializing 1D linear interpolator'

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x (deg)',ylabel='f(x)',&
                        title='Linear Test',legend=.true.,figsize=[10,5])
    call plt%add_plot(x,fcn_1d,&
                        label='Function $f(x) = \\sin(x)$',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)

    errmax = 0.0_wp
    do i=1,nxv
        call b2%evaluate(xval(i),idx,val,iflag)
        fval(i) = val  ! save it for plot
        if (iflag/=0) error stop 'error evaluating linear spline: '//get_status_message(iflag)

        ! also compute linear interpolation:
        call s1%evaluate(xval(i),val)
        fval_linear(i) = val ! linear vector for plot

        err    = abs(fval(i) - fval_linear(i))
        errmax = max(err,errmax)

        ! also others:
        call b3%evaluate(xval(i),idx,fval3(i),iflag)
        if (iflag/=0) error stop 'error evaluating quadratic spline: '//get_status_message(iflag)
        call b4%evaluate(xval(i),idx,fval4(i),iflag)
        if (iflag/=0) error stop 'error evaluating cubic spline: '//get_status_message(iflag)

    end do

    write(*,*) ''
    write(*,*) 'Max difference (spline - linear):', errmax
    write(*,*) ''

    call plt%add_plot(xval,fval_linear,&
            label='Linear Interpolated',&
            linestyle='r.',linewidth=1,istat=istat)
    call plt%add_plot(xval,fval,&
            label='Linear ($k=2$) Spline Interpolated',&
            linestyle='r-',linewidth=1,istat=istat)
    call plt%add_plot(xval,fval3,&
            label='Quadratic ($k=3$) Spline Interpolated',&
            linestyle='k.-',linewidth=1,istat=istat)
    call plt%add_plot(xval,fval4,&
            label='Cubic ($k=4$) Spline Interpolated',&
            linestyle='c.-',linewidth=1,istat=istat)
    call plt%savefig('bspline_linear_test.png',istat=istat)

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp),intent(in) :: x
        real(wp),parameter :: a = acos(-1.0_wp)/18.0_wp
        f1 = sin(a*x)
        end function f1

    end program bspline_linear_test
