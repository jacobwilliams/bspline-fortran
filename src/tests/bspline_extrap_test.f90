!*****************************************************************************************
!>
!  Extrapolation test
!
!### Results
!  ![Plot of results](https://raw.githubusercontent.com/jacobwilliams/bspline-fortran/master/src/tests/bspline_extrap_test.png)

    program bspline_extrap_test

    use bspline_module
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module

    implicit none

    integer,parameter :: nx    = 400   !! number of points in x
    integer,parameter :: nxv   = 800   !! number of points to evaluate interpolant

    integer,parameter :: kx    = 4    !! order in x
    integer,parameter :: iknot = 0    !! automatically select the knots

    real(wp) :: x(nx)
    real(wp) :: xval(nxv),fval(nxv)
    real(wp) :: tx(nx+kx)
    real(wp) :: fcn_1d(nx)
    real(wp) :: val,tru,err,errmax
    integer  :: i,j,idx,iflag,inbvx,iloy
    logical :: extrap
    type(pyplot) :: plt

    real(wp),parameter :: rad2deg = 180.0_wp / acos(-1.0_wp)  !! deg. to radians conversion factor

    idx = 0

    do i=1,nx
        x(i) = real(i-1,wp)/100.0_wp + 0.0001_wp
        fcn_1d(i) = f1(x(i))
    end do
    do i=1,nxv
        xval(i) = real(i-200,wp)/100.0_wp
    end do

    !have to set these before the first evaluate call:
    inbvx = 1
    iloy  = 1

    ! initialize
    call db1ink(x,nx,fcn_1d,kx,iknot,tx,fcn_1d,iflag)

    if (iflag/=0) then
        write(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
    end if

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x (deg)',ylabel='f(x)',&
                        title='Extrapolation Test',legend=.true.)
    call plt%add_plot(x*rad2deg,fcn_1d,label='Function $f(x) = \sin(x)$',&
                        linestyle='ko',markersize=5,linewidth=2)

    do j=1,2

        ! run once without extrapolation, and once with extrapolation
        extrap = j==2

        errmax = 0.0_wp
        do i=1,nxv
            call db1val(xval(i),idx,tx,nx,kx,fcn_1d,val,iflag,inbvx,extrap=extrap)
            fval(i) = val  ! save it for plot
            tru    = f1(xval(i))
            err    = abs(tru-val)
            errmax = max(err,errmax)
            !write(*,*) xval(i), val, tru, err, iflag
        end do

        ! check max error against tolerance
        write(*,*) ''
        write(*,*) '1D: max error:', errmax
        write(*,*) ''

        if (extrap) then
            call plt%add_plot(xval*rad2deg,fval,&
                    label='Interpolated',&
                    linestyle='g.-',linewidth=1)
            call plt%savefig('bspline_extrap_test.png')
        end if

    end do

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        f1 = sin(x)
        end function f1

    end program bspline_extrap_test
