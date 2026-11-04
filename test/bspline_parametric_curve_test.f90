!*****************************************************************************************
!>
!  Example: evaluate a parametric B-spline curve from a CAD program
!  (FreeCAD) as a function z(x) (issue #96).
!
!  A CAD B-spline curve is parametric: each coordinate is a B-spline in
!  the parameter `u`, with the same knots, and the coordinates of the
!  poles (control points) are the coefficients. So `x(u)` and `z(u)` are
!  two 1D splines. To get `z` at a given `x`, solve `x(u) = x` for `u`
!  (with `inverse`), and then evaluate `z(u)`.
!
!  The curve here is from the issue: degree 3 (order 4), 7 poles, knots
!  `[0, 229.98, 456.56, 682.58, 912.14]` with multiplicities `[4,1,1,1,4]`.
!  The `y` coordinates of the poles are about zero, so they are not used.

    program bspline_parametric_curve_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: k = 4  !! spline order (degree + 1)

    real(wp),dimension(11),parameter :: t = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                             229.98020932693515_wp, 456.5562217064585_wp, &
                                             682.584333954564_wp, &
                                             912.1417114001321_wp, 912.1417114001321_wp, &
                                             912.1417114001321_wp, 912.1417114001321_wp]
                                             !! knot vector (the knots repeated by their multiplicities)
    real(wp),dimension(7),parameter :: px = [347.119461269_wp, 421.35826684612715_wp, &
                                             568.0349611242125_wp, 787.6981915383186_wp, &
                                             1017.878120964171_wp, 1166.9796636987746_wp, &
                                             1239.265344_wp]  !! x coordinates of the poles
    real(wp),dimension(7),parameter :: pz = [251.048911919_wp, 234.39242900777492_wp, &
                                             190.78154071596157_wp, 130.3077433064477_wp, &
                                             122.11454085418023_wp, 154.29592542750504_wp, &
                                             182.878520542_wp]  !! z coordinates of the poles

    type(bspline_1d) :: sx  !! the spline x(u)
    type(bspline_1d) :: sz  !! the spline z(u)
    real(wp) :: x       !! the abscissa where z is wanted
    real(wp) :: u       !! the parameter where x(u) = x
    real(wp) :: z       !! z(u)
    real(wp) :: xu      !! x(u), to check the solution
    integer(ip) :: iflag  !! status flag
    integer(ip) :: i      !! counter
    integer :: nfail      !! number of failed checks

    nfail = 0

    call sx%initialize_from_coefficients(t, px, k, iflag)
    call check('initialize x(u)', iflag==0)
    call sz%initialize_from_coefficients(t, pz, k, iflag)
    call check('initialize z(u)', iflag==0)

    write(*,'(A)') '         x              u              z'
    do i = 0, 8
        x = px(1) + (px(7)-px(1))*real(i,wp)/8.0_wp
        call sx%inverse(x, u, iflag)        ! solve x(u) = x for u
        call check('inverse', iflag==0)
        call sz%evaluate(u, 0_ip, z, iflag) ! then z = z(u)
        call check('evaluate', iflag==0)
        write(*,'(3F15.6)') x, u, z
        call sx%evaluate(u, 0_ip, xu, iflag)
        call check('x(u) = x', abs(xu-x)<=1.0e-10_wp*abs(x))
        ! the curve lies within the range of its control points
        call check('z within the poles', z>=minval(pz) .and. z<=maxval(pz))
    end do

    ! clamped knots: the curve starts and ends at the first and last poles
    call sx%inverse(px(1), u, iflag)
    call sz%evaluate(u, 0_ip, z, iflag)
    call check('first pole', u==t(1) .and. abs(z-pz(1))<=1.0e-12_wp*abs(pz(1)))
    call sx%inverse(px(7), u, iflag)
    call sz%evaluate(u, 0_ip, z, iflag)
    call check('last pole', abs(u-t(11))<=1.0e-12_wp*t(11) .and. abs(z-pz(7))<=1.0e-12_wp*abs(pz(7)))

    ! outside the curve: no solution
    call sx%inverse(px(1)-1.0_wp, u, iflag)
    call check('no solution outside the curve', iflag==4005)

    write(*,*) ''
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'parametric curve test failed'
    else
        write(*,*) 'all parametric curve tests PASSED'
    end if

    contains

    subroutine check(name, ok)

    !! Record a pass/fail result.

    character(len=*),intent(in) :: name  !! name of the check
    logical,intent(in)          :: ok    !! true if the check passed

    if (.not. ok) then
        write(*,'(A)') 'FAIL: '//name
        nfail = nfail + 1
    end if

    end subroutine check

    end program bspline_parametric_curve_test
!*****************************************************************************************
