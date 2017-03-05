!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/6/2015
!
!  Tests for different knot vectors.
!
!  This requires [pyplot_module](https://github.com/jacobwilliams/pyplot-fortran)
!
!### Results
!  ![Plot of results](https://raw.githubusercontent.com/jacobwilliams/bspline-fortran/master/src/tests/knot_tests.png)

    program knot_tests

    use bspline_module
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module

    implicit none

    integer :: i    !! counter

    integer,parameter :: kx  = 4    !! x bspline order
    integer,parameter :: nx  = 6    !! number of points in x dimension
    real(wp),dimension(nx),parameter :: x = [(real(i*10,wp), i=0,11,(10+2)/nx)]    !! [0,20,40,60,80,100]

    real(wp),dimension(nx)    :: fcn
    real(wp),dimension(nx+kx) :: tx
    type(bspline_1d)          :: s_default, s1, s2
    real(wp),dimension(0:100) :: x_new,f_new_default,f1,f2 !,f_actual
    real(wp)                  :: xval
    type(pyplot)              :: plt
    integer                   :: iflag

    !function evaluations for original grid:
    do i=1,nx
        fcn(i) = test_func(x(i))
    end do

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='f(x)',&
                        title='Knot Test',legend=.true.)
    call plt%add_plot(x,fcn,label='Function $f(x) = \sin(x \cdot \pi/18)$ : $x=[0,20,40,60,80,100]$',&
                        linestyle='ko--',markersize=5,linewidth=2)

    !initialize three with different knot sequences:

    ! x = [ 0,      20,40,60,80,100                ] (x points)
    ! t = [ 0,0,0,0,   40,60,      101,101,101,101 ] (tx for not a knot conditions)

    call s_default%initialize(x,fcn,kx,iflag)  !default (not-a-knot)
    if (iflag/=0) error stop 'error initializing s_default'

    !user-specified knots:

    tx = real([0,0,0,0,20,40,101,101,101,101], wp)
    call s1%initialize(x,fcn,kx,tx,iflag)
    if (iflag/=0) error stop 'error initializing s1'

    tx = real([0,0,0,0,60,80,101,101,101,101], wp)
    call s2%initialize(x,fcn,kx,tx,iflag)
    if (iflag/=0) error stop 'error initializing s2'

    do i = 0,100

        xval     = real(i,wp)
        x_new(i) = xval

        !f_actual(i) = test_func(xval)

        call s_default%evaluate(xval,0,f_new_default(i),iflag)
        if (iflag/=0) error stop 'error evaluating s_default'

        call s1%evaluate(xval,0,f1(i),iflag)
        if (iflag/=0) error stop 'error evaluating s1'

        call s2%evaluate(xval,0,f2(i),iflag)
        if (iflag/=0) error stop 'error evaluating s2'

    end do

    !call plt%add_plot(x_new,f_actual,label='Actual function',linestyle='k--',linewidth=2)

    call plt%add_plot(x_new,f_new_default,&
            label='Interpolated : $t_x=[0,0,0,0,40,60,102,102,102,102]$ (Default)',&
            linestyle='b-',linewidth=1)
    call plt%add_plot(x_new,f1,&
            label='Interpolated : $t_x=[0,0,0,0,20,40,101,101,101,101]$',&
            linestyle='r-',linewidth=1)
    call plt%add_plot(x_new,f2,&
            label='Interpolated : $t_x=[0,0,0,0,60,80,101,101,101,101]$',&
            linestyle='g-',linewidth=1)

    !plot the results:
    call plt%savefig('knot_tests.png')

    contains

        pure function test_func(x) result(f)
        !! 1d test function

        implicit none

        real(wp) :: f
        real(wp),intent(in) :: x

        real(wp),parameter :: a = acos(-1.0_wp)/18.0_wp

        f = sin(a*x)

        end function test_func

    end program knot_tests
