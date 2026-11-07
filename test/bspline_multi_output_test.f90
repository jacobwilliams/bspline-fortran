!*****************************************************************************************
!>
!  Tests for splines with multiple outputs (several functions of the same
!  inputs, on the same grid): the multiple-output versions of `db1ink`..`db6ink`
!  and `db1val`..`db6val`, and of the `initialize`, `evaluate`,
!  `initialize_from_coefficients`, and `get_coefficients` methods (issue #97).
!
!  * The values (and derivatives) must be the same as for separate
!    single-output splines: bit for bit in 2D-6D (the same arithmetic is
!    used), and to rounding in 1D (where the single-output version uses
!    a different algorithm).
!  * The coefficients must be the same as for separate splines.
!  * Error codes 6001-6005.
!  * The timing of the evaluation is compared with separate splines
!    (printed only).

    program bspline_multi_output_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer :: nfail   !! number of failed checks
    integer :: ncheck  !! number of checks
    integer :: nseed   !! size of the random seed
    integer,dimension(:),allocatable :: seed  !! random seed

    nfail = 0
    ncheck = 0
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 97
    call random_seed(put=seed)

    call test_1d()
    call test_2d()
    call test_3d()
    call test_4d()
    call test_5d()
    call test_6d()
    call test_1d_methods()
    call timing()

    write(*,*) ''
    write(*,*) ncheck, ' checks'
    if (nfail > 0) then
        write(*,*) nfail, ' checks FAILED'
        error stop 'multiple output test failed'
    else
        write(*,*) 'all multiple output tests PASSED'
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
!  A strictly increasing grid of `n` points with random spacing.

    function grid(n) result(x)

    integer(ip),intent(in) :: n  !! number of points
    real(wp),dimension(n)  :: x  !! the grid

    integer(ip) :: i  !! counter
    real(wp) :: u     !! random number

    x(1) = 0.0_wp
    do i = 2, n
        call random_number(u)
        x(i) = x(i-1) + 0.2_wp + u
    end do

    end function grid
!*****************************************************************************************

!*****************************************************************************************
!>
!  1D: multiple outputs compared with separate splines.

    subroutine test_1d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_1d) :: sm      !! spline with multiple outputs
    type(bspline_1d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_1d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(1) :: v  !! evaluation point
    integer(ip),dimension(1) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx  !! status flag, counters, and search state
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 1
            n(i) = 12 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 5), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        allocate(fcn(n(1),nf))
        call random_number(fcn)
        allocate(tx(n(1)+k(1)))

        call sm%initialize(x,fcn,k(1),iflag)
        call check('1d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,fcn(:,j),k(1),iflag)
        end do
        call s1(1)%get_knots(tx,iflag)
        call sm2%initialize(x,fcn,k(1),tx,iflag)
        call check('1d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),nf))
        call db1ink(x,n(1),fcn,k(1),1_ip,tx,bcoef2,iflag)
        call check('1d db1ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1
        do i = 1, 50
            do j = 1, 1
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            call sm%evaluate(v(1),id(1),fm,iflag)
            if (iflag/=0) call check('1d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),id(1),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db1val(v(1),id(1),tx,n(1),k(1),bcoef2,fm2,iflag,inbvx)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),id(1),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('1d multi = separate splines (to rounding)', err<=1.0e-13_wp)
        call check('1d multi: specified knots give the same values', same)
        call check('1d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('1d multi get_coefficients', iflag==0 .and. size(bcoef,2)==nf)
        call check('1d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('1d multi coefficients = separate splines', all(bcoef(:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,bcoef,k(1),iflag)
        call check('1d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),id(1),fm2,iflag)
        call sm%evaluate(v(1),id(1),fm,iflag)
        call check('1d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),id(1),f1,iflag)
        call check('1d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('1d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),id(1),fbad,iflag)
        call check('1d wrong size of f returns 6002', iflag==6002)
        call db1ink(x,n(1),fcn,k(1),1_ip,tx,bcoef2(:,1:nf-1),iflag)
        call check('1d db1ink multi with mismatched outputs returns 6001', iflag==6001)
        call db1val(v(1),id(1),tx,n(1),k(1),bcoef2(1:n(1)-1,:),fm2,iflag,inbvx)
        call check('1d db1val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,bcoef(:,1:0),k(1),iflag)
        call check('1d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('1d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx)
    end do

    end subroutine test_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  2D: multiple outputs compared with separate splines.

    subroutine test_2d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_2d) :: sm      !! spline with multiple outputs
    type(bspline_2d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_2d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:),allocatable :: y, ty  !! data abscissae and knots in \(y\)
    real(wp),dimension(:,:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:,:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(2) :: v  !! evaluation point
    integer(ip),dimension(2) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx, iloy  !! status flag, counters, and search state
    real(wp),dimension(:),allocatable :: w1  !! work array for [[db2val]]
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 2
            n(i) = 8 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 4), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        y = grid(n(2))
        allocate(fcn(n(1),n(2),nf))
        call random_number(fcn)
        allocate(w1(k(2)))
        allocate(tx(n(1)+k(1)))
        allocate(ty(n(2)+k(2)))

        call sm%initialize(x,y,fcn,k(1),k(2),iflag)
        call check('2d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,y,fcn(:,:,j),k(1),k(2),iflag)
        end do
        call s1(1)%get_knots(tx,ty,iflag)
        call sm2%initialize(x,y,fcn,k(1),k(2),tx,ty,iflag)
        call check('2d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),n(2),nf))
        call db2ink(x,n(1),y,n(2),fcn,k(1),k(2),1_ip,tx,ty,bcoef2,iflag)
        call check('2d db2ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1; iloy = 1
        do i = 1, 50
            do j = 1, 2
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            v(2) = y(1) + v(2)*(y(n(2))-y(1))
            call sm%evaluate(v(1),v(2),id(1),id(2),fm,iflag)
            if (iflag/=0) call check('2d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),v(2),id(1),id(2),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db2val(v(1),v(2),id(1),id(2),tx,ty,n(1),n(2),k(1),k(2),bcoef2,fm2,iflag,inbvx,iloy,w1)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),v(2),id(1),id(2),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('2d multi = separate splines (bit for bit)', err==0.0_wp)
        call check('2d multi: specified knots give the same values', same)
        call check('2d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('2d multi get_coefficients', iflag==0 .and. size(bcoef,3)==nf)
        call check('2d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('2d multi coefficients = separate splines', all(bcoef(:,:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,ty,bcoef,k(1),k(2),iflag)
        call check('2d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),v(2),id(1),id(2),fm2,iflag)
        call sm%evaluate(v(1),v(2),id(1),id(2),fm,iflag)
        call check('2d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),v(2),id(1),id(2),f1,iflag)
        call check('2d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('2d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),v(2),id(1),id(2),fbad,iflag)
        call check('2d wrong size of f returns 6002', iflag==6002)
        call db2ink(x,n(1),y,n(2),fcn,k(1),k(2),1_ip,tx,ty,bcoef2(:,:,1:nf-1),iflag)
        call check('2d db2ink multi with mismatched outputs returns 6001', iflag==6001)
        call db2val(v(1),v(2),id(1),id(2),tx,ty,n(1),n(2),k(1),k(2),bcoef2(1:n(1)-1,:,:),fm2,iflag,inbvx,iloy,w1)
        call check('2d db2val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,ty,bcoef(:,:,1:0),k(1),k(2),iflag)
        call check('2d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('2d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx, ty, w1)
    end do

    end subroutine test_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  3D: multiple outputs compared with separate splines.

    subroutine test_3d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_3d) :: sm      !! spline with multiple outputs
    type(bspline_3d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_3d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:),allocatable :: y, ty  !! data abscissae and knots in \(y\)
    real(wp),dimension(:),allocatable :: z, tz  !! data abscissae and knots in \(z\)
    real(wp),dimension(:,:,:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:,:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:,:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:,:,:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(3) :: v  !! evaluation point
    integer(ip),dimension(3) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx, iloy, iloz  !! status flag, counters, and search state
    real(wp),dimension(:,:),allocatable :: w2  !! work array for [[db3val]]
    real(wp),dimension(:),allocatable :: w1  !! work array for [[db3val]]
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 3
            n(i) = 6 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 3), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        y = grid(n(2))
        z = grid(n(3))
        allocate(fcn(n(1),n(2),n(3),nf))
        call random_number(fcn)
        allocate(w2(k(2),k(3)))
        allocate(w1(k(3)))
        allocate(tx(n(1)+k(1)))
        allocate(ty(n(2)+k(2)))
        allocate(tz(n(3)+k(3)))

        call sm%initialize(x,y,z,fcn,k(1),k(2),k(3),iflag)
        call check('3d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,y,z,fcn(:,:,:,j),k(1),k(2),k(3),iflag)
        end do
        call s1(1)%get_knots(tx,ty,tz,iflag)
        call sm2%initialize(x,y,z,fcn,k(1),k(2),k(3),tx,ty,tz,iflag)
        call check('3d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),n(2),n(3),nf))
        call db3ink(x,n(1),y,n(2),z,n(3),fcn,k(1),k(2),k(3),1_ip,tx,ty,tz,bcoef2,iflag)
        call check('3d db3ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1; iloy = 1; iloz = 1
        do i = 1, 50
            do j = 1, 3
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            v(2) = y(1) + v(2)*(y(n(2))-y(1))
            v(3) = z(1) + v(3)*(z(n(3))-z(1))
            call sm%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),fm,iflag)
            if (iflag/=0) call check('3d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db3val(v(1),v(2),v(3),id(1),id(2),id(3),tx,ty,tz,n(1),n(2),n(3),k(1),k(2),k(3),bcoef2,fm2,iflag,inbvx,iloy,iloz,w2,w1)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('3d multi = separate splines (bit for bit)', err==0.0_wp)
        call check('3d multi: specified knots give the same values', same)
        call check('3d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('3d multi get_coefficients', iflag==0 .and. size(bcoef,4)==nf)
        call check('3d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('3d multi coefficients = separate splines', all(bcoef(:,:,:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,ty,tz,bcoef,k(1),k(2),k(3),iflag)
        call check('3d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),fm2,iflag)
        call sm%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),fm,iflag)
        call check('3d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),f1,iflag)
        call check('3d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('3d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),v(2),v(3),id(1),id(2),id(3),fbad,iflag)
        call check('3d wrong size of f returns 6002', iflag==6002)
        call db3ink(x,n(1),y,n(2),z,n(3),fcn,k(1),k(2),k(3),1_ip,tx,ty,tz,bcoef2(:,:,:,1:nf-1),iflag)
        call check('3d db3ink multi with mismatched outputs returns 6001', iflag==6001)
        call db3val(v(1),v(2),v(3),id(1),id(2),id(3),tx,ty,tz,n(1),n(2),n(3),k(1),k(2),k(3),bcoef2(1:n(1)-1,:,:,:),fm2,iflag,inbvx,iloy,iloz,w2,w1)
        call check('3d db3val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,ty,tz,bcoef(:,:,:,1:0),k(1),k(2),k(3),iflag)
        call check('3d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('3d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx, ty, tz, w2, w1)
    end do

    end subroutine test_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  4D: multiple outputs compared with separate splines.

    subroutine test_4d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_4d) :: sm      !! spline with multiple outputs
    type(bspline_4d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_4d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:),allocatable :: y, ty  !! data abscissae and knots in \(y\)
    real(wp),dimension(:),allocatable :: z, tz  !! data abscissae and knots in \(z\)
    real(wp),dimension(:),allocatable :: q, tq  !! data abscissae and knots in \(q\)
    real(wp),dimension(:,:,:,:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:,:,:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:,:,:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:,:,:,:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(4) :: v  !! evaluation point
    integer(ip),dimension(4) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx, iloy, iloz, iloq  !! status flag, counters, and search state
    real(wp),dimension(:,:,:),allocatable :: w3  !! work array for [[db4val]]
    real(wp),dimension(:,:),allocatable :: w2  !! work array for [[db4val]]
    real(wp),dimension(:),allocatable :: w1  !! work array for [[db4val]]
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 4
            n(i) = 5 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 3), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        y = grid(n(2))
        z = grid(n(3))
        q = grid(n(4))
        allocate(fcn(n(1),n(2),n(3),n(4),nf))
        call random_number(fcn)
        allocate(w3(k(2),k(3),k(4)))
        allocate(w2(k(3),k(4)))
        allocate(w1(k(4)))
        allocate(tx(n(1)+k(1)))
        allocate(ty(n(2)+k(2)))
        allocate(tz(n(3)+k(3)))
        allocate(tq(n(4)+k(4)))

        call sm%initialize(x,y,z,q,fcn,k(1),k(2),k(3),k(4),iflag)
        call check('4d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,y,z,q,fcn(:,:,:,:,j),k(1),k(2),k(3),k(4),iflag)
        end do
        call s1(1)%get_knots(tx,ty,tz,tq,iflag)
        call sm2%initialize(x,y,z,q,fcn,k(1),k(2),k(3),k(4),tx,ty,tz,tq,iflag)
        call check('4d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),n(2),n(3),n(4),nf))
        call db4ink(x,n(1),y,n(2),z,n(3),q,n(4),fcn,k(1),k(2),k(3),k(4),1_ip,tx,ty,tz,tq,bcoef2,iflag)
        call check('4d db4ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1; iloy = 1; iloz = 1; iloq = 1
        do i = 1, 50
            do j = 1, 4
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            v(2) = y(1) + v(2)*(y(n(2))-y(1))
            v(3) = z(1) + v(3)*(z(n(3))-z(1))
            v(4) = q(1) + v(4)*(q(n(4))-q(1))
            call sm%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),fm,iflag)
            if (iflag/=0) call check('4d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db4val(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),tx,ty,tz,tq,n(1),n(2),n(3),n(4),k(1),k(2),k(3),k(4),bcoef2,fm2,iflag,inbvx,iloy,iloz,iloq,w3,w2,w1)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('4d multi = separate splines (bit for bit)', err==0.0_wp)
        call check('4d multi: specified knots give the same values', same)
        call check('4d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('4d multi get_coefficients', iflag==0 .and. size(bcoef,5)==nf)
        call check('4d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('4d multi coefficients = separate splines', all(bcoef(:,:,:,:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,bcoef,k(1),k(2),k(3),k(4),iflag)
        call check('4d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),fm2,iflag)
        call sm%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),fm,iflag)
        call check('4d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),f1,iflag)
        call check('4d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('4d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),fbad,iflag)
        call check('4d wrong size of f returns 6002', iflag==6002)
        call db4ink(x,n(1),y,n(2),z,n(3),q,n(4),fcn,k(1),k(2),k(3),k(4),1_ip,tx,ty,tz,tq,bcoef2(:,:,:,:,1:nf-1),iflag)
        call check('4d db4ink multi with mismatched outputs returns 6001', iflag==6001)
        call db4val(v(1),v(2),v(3),v(4),id(1),id(2),id(3),id(4),tx,ty,tz,tq,n(1),n(2),n(3),n(4),k(1),k(2),k(3),k(4),bcoef2(1:n(1)-1,:,:,:,:),fm2,iflag,inbvx,iloy,iloz,iloq,w3,w2,w1)
        call check('4d db4val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,bcoef(:,:,:,:,1:0),k(1),k(2),k(3),k(4),iflag)
        call check('4d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('4d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx, ty, tz, tq, w3, w2, w1)
    end do

    end subroutine test_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  5D: multiple outputs compared with separate splines.

    subroutine test_5d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_5d) :: sm      !! spline with multiple outputs
    type(bspline_5d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_5d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:),allocatable :: y, ty  !! data abscissae and knots in \(y\)
    real(wp),dimension(:),allocatable :: z, tz  !! data abscissae and knots in \(z\)
    real(wp),dimension(:),allocatable :: q, tq  !! data abscissae and knots in \(q\)
    real(wp),dimension(:),allocatable :: r, tr  !! data abscissae and knots in \(r\)
    real(wp),dimension(:,:,:,:,:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:,:,:,:,:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(5) :: v  !! evaluation point
    integer(ip),dimension(5) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx, iloy, iloz, iloq, ilor  !! status flag, counters, and search state
    real(wp),dimension(:,:,:,:),allocatable :: w4  !! work array for [[db5val]]
    real(wp),dimension(:,:,:),allocatable :: w3  !! work array for [[db5val]]
    real(wp),dimension(:,:),allocatable :: w2  !! work array for [[db5val]]
    real(wp),dimension(:),allocatable :: w1  !! work array for [[db5val]]
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 5
            n(i) = 4 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 2), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        y = grid(n(2))
        z = grid(n(3))
        q = grid(n(4))
        r = grid(n(5))
        allocate(fcn(n(1),n(2),n(3),n(4),n(5),nf))
        call random_number(fcn)
        allocate(w4(k(2),k(3),k(4),k(5)))
        allocate(w3(k(3),k(4),k(5)))
        allocate(w2(k(4),k(5)))
        allocate(w1(k(5)))
        allocate(tx(n(1)+k(1)))
        allocate(ty(n(2)+k(2)))
        allocate(tz(n(3)+k(3)))
        allocate(tq(n(4)+k(4)))
        allocate(tr(n(5)+k(5)))

        call sm%initialize(x,y,z,q,r,fcn,k(1),k(2),k(3),k(4),k(5),iflag)
        call check('5d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,y,z,q,r,fcn(:,:,:,:,:,j),k(1),k(2),k(3),k(4),k(5),iflag)
        end do
        call s1(1)%get_knots(tx,ty,tz,tq,tr,iflag)
        call sm2%initialize(x,y,z,q,r,fcn,k(1),k(2),k(3),k(4),k(5),tx,ty,tz,tq,tr,iflag)
        call check('5d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),n(2),n(3),n(4),n(5),nf))
        call db5ink(x,n(1),y,n(2),z,n(3),q,n(4),r,n(5),fcn,k(1),k(2),k(3),k(4),k(5),1_ip,tx,ty,tz,tq,tr,bcoef2,iflag)
        call check('5d db5ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1; iloy = 1; iloz = 1; iloq = 1; ilor = 1
        do i = 1, 50
            do j = 1, 5
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            v(2) = y(1) + v(2)*(y(n(2))-y(1))
            v(3) = z(1) + v(3)*(z(n(3))-z(1))
            v(4) = q(1) + v(4)*(q(n(4))-q(1))
            v(5) = r(1) + v(5)*(r(n(5))-r(1))
            call sm%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),fm,iflag)
            if (iflag/=0) call check('5d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db5val(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),tx,ty,tz,tq,tr,n(1),n(2),n(3),n(4),n(5),k(1),k(2),k(3),k(4),k(5),bcoef2,fm2,iflag,inbvx,iloy,iloz,iloq,ilor,w4,w3,w2,w1)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('5d multi = separate splines (bit for bit)', err==0.0_wp)
        call check('5d multi: specified knots give the same values', same)
        call check('5d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('5d multi get_coefficients', iflag==0 .and. size(bcoef,6)==nf)
        call check('5d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('5d multi coefficients = separate splines', all(bcoef(:,:,:,:,:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,tr,bcoef,k(1),k(2),k(3),k(4),k(5),iflag)
        call check('5d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),fm2,iflag)
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),fm,iflag)
        call check('5d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),f1,iflag)
        call check('5d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('5d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),fbad,iflag)
        call check('5d wrong size of f returns 6002', iflag==6002)
        call db5ink(x,n(1),y,n(2),z,n(3),q,n(4),r,n(5),fcn,k(1),k(2),k(3),k(4),k(5),1_ip,tx,ty,tz,tq,tr,bcoef2(:,:,:,:,:,1:nf-1),iflag)
        call check('5d db5ink multi with mismatched outputs returns 6001', iflag==6001)
        call db5val(v(1),v(2),v(3),v(4),v(5),id(1),id(2),id(3),id(4),id(5),tx,ty,tz,tq,tr,n(1),n(2),n(3),n(4),n(5),k(1),k(2),k(3),k(4),k(5),bcoef2(1:n(1)-1,:,:,:,:,:),fm2,iflag,inbvx,iloy,iloz,iloq,ilor,w4,w3,w2,w1)
        call check('5d db5val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,tr,bcoef(:,:,:,:,:,1:0),k(1),k(2),k(3),k(4),k(5),iflag)
        call check('5d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('5d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx, ty, tz, tq, tr, w4, w3, w2, w1)
    end do

    end subroutine test_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  6D: multiple outputs compared with separate splines.

    subroutine test_6d()

    integer(ip),parameter :: nf = 3  !! number of outputs
    type(bspline_6d) :: sm      !! spline with multiple outputs
    type(bspline_6d) :: sm2     !! spline with multiple outputs (specified knots, or from coefficients)
    type(bspline_6d),dimension(nf) :: s1  !! separate splines
    real(wp),dimension(:),allocatable :: x, tx  !! data abscissae and knots in \(x\)
    real(wp),dimension(:),allocatable :: y, ty  !! data abscissae and knots in \(y\)
    real(wp),dimension(:),allocatable :: z, tz  !! data abscissae and knots in \(z\)
    real(wp),dimension(:),allocatable :: q, tq  !! data abscissae and knots in \(q\)
    real(wp),dimension(:),allocatable :: r, tr  !! data abscissae and knots in \(r\)
    real(wp),dimension(:),allocatable :: s, ts  !! data abscissae and knots in \(s\)
    real(wp),dimension(:,:,:,:,:,:,:),allocatable :: fcn    !! data for all the outputs
    real(wp),dimension(:,:,:,:,:,:,:),allocatable :: bcoef  !! coefficients for all the outputs
    real(wp),dimension(:,:,:,:,:,:,:),allocatable :: bcoef2 !! coefficients from the subroutine interface
    real(wp),dimension(:,:,:,:,:,:),allocatable :: b1     !! coefficients of one output
    real(wp),dimension(nf) :: fm   !! values of all the outputs
    real(wp),dimension(nf) :: fm2  !! values of all the outputs (subroutine interface)
    real(wp),dimension(nf+1) :: fbad  !! wrong size
    real(wp) :: f1     !! value of one output
    real(wp) :: err    !! maximum difference
    real(wp) :: err2   !! maximum difference (subroutine interface)
    real(wp),dimension(6) :: v  !! evaluation point
    integer(ip),dimension(6) :: n, k, id  !! number of points, orders, derivatives
    integer(ip) :: iflag, j, i, itrial, inbvx, iloy, iloz, iloq, ilor, ilos  !! status flag, counters, and search state
    real(wp),dimension(:,:,:,:,:),allocatable :: w5  !! work array for [[db6val]]
    real(wp),dimension(:,:,:,:),allocatable :: w4  !! work array for [[db6val]]
    real(wp),dimension(:,:,:),allocatable :: w3  !! work array for [[db6val]]
    real(wp),dimension(:,:),allocatable :: w2  !! work array for [[db6val]]
    real(wp),dimension(:),allocatable :: w1  !! work array for [[db6val]]
    logical :: same  !! all values are the same

    do itrial = 1, 5
        do i = 1, 6
            n(i) = 4 - mod(itrial+i, 2)
            k(i) = min(2 + mod(itrial*i, 2), n(i)-1)  ! order < number of points
        end do
        x = grid(n(1))
        y = grid(n(2))
        z = grid(n(3))
        q = grid(n(4))
        r = grid(n(5))
        s = grid(n(6))
        allocate(fcn(n(1),n(2),n(3),n(4),n(5),n(6),nf))
        call random_number(fcn)
        allocate(w5(k(2),k(3),k(4),k(5),k(6)))
        allocate(w4(k(3),k(4),k(5),k(6)))
        allocate(w3(k(4),k(5),k(6)))
        allocate(w2(k(5),k(6)))
        allocate(w1(k(6)))
        allocate(tx(n(1)+k(1)))
        allocate(ty(n(2)+k(2)))
        allocate(tz(n(3)+k(3)))
        allocate(tq(n(4)+k(4)))
        allocate(tr(n(5)+k(5)))
        allocate(ts(n(6)+k(6)))

        call sm%initialize(x,y,z,q,r,s,fcn,k(1),k(2),k(3),k(4),k(5),k(6),iflag)
        call check('6d multi initialize', iflag==0 .and. sm%number_of_outputs()==nf)
        do j = 1, nf
            call s1(j)%initialize(x,y,z,q,r,s,fcn(:,:,:,:,:,:,j),k(1),k(2),k(3),k(4),k(5),k(6),iflag)
        end do
        call s1(1)%get_knots(tx,ty,tz,tq,tr,ts,iflag)
        call sm2%initialize(x,y,z,q,r,s,fcn,k(1),k(2),k(3),k(4),k(5),k(6),tx,ty,tz,tq,tr,ts,iflag)
        call check('6d multi initialize (specified knots)', iflag==0)
        allocate(bcoef2(n(1),n(2),n(3),n(4),n(5),n(6),nf))
        call db6ink(x,n(1),y,n(2),z,n(3),q,n(4),r,n(5),s,n(6),fcn,k(1),k(2),k(3),k(4),k(5),k(6),1_ip,tx,ty,tz,tq,tr,ts,bcoef2,iflag)
        call check('6d db6ink multi', iflag==0)
        err = 0.0_wp; err2 = 0.0_wp; same = .true.
        inbvx = 1; iloy = 1; iloz = 1; iloq = 1; ilor = 1; ilos = 1
        do i = 1, 50
            do j = 1, 6
                call random_number(v(j))
                id(j) = mod(int(v(j)*1000.0_wp,ip), k(j))  ! random derivative order
            end do
            v(1) = x(1) + v(1)*(x(n(1))-x(1))
            v(2) = y(1) + v(2)*(y(n(2))-y(1))
            v(3) = z(1) + v(3)*(z(n(3))-z(1))
            v(4) = q(1) + v(4)*(q(n(4))-q(1))
            v(5) = r(1) + v(5)*(r(n(5))-r(1))
            v(6) = s(1) + v(6)*(s(n(6))-s(1))
            call sm%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),fm,iflag)
            if (iflag/=0) call check('6d multi evaluate iflag', .false.)
            call sm2%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),fm2,iflag)
            same = same .and. all(fm2==fm)
            call db6val(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),tx,ty,tz,tq,tr,ts,n(1),n(2),n(3),n(4),n(5),n(6),k(1),k(2),k(3),k(4),k(5),k(6),bcoef2,fm2,iflag,inbvx,iloy,iloz,iloq,ilor,ilos,w5,w4,w3,w2,w1)
            err2 = max(err2, maxval(abs(fm2-fm)))
            do j = 1, nf
                call s1(j)%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),f1,iflag)
                err = max(err, abs(f1-fm(j))/max(1.0_wp,abs(f1)))
            end do
        end do
        call check('6d multi = separate splines (bit for bit)', err==0.0_wp)
        call check('6d multi: specified knots give the same values', same)
        call check('6d multi: subroutine interface gives the same values', err2==0.0_wp)
        call sm%get_coefficients(bcoef,iflag)
        call check('6d multi get_coefficients', iflag==0 .and. size(bcoef,7)==nf)
        call check('6d multi coefficients = subroutine interface', all(bcoef==bcoef2))
        do j = 1, nf
            call s1(j)%get_coefficients(b1,iflag)
            call check('6d multi coefficients = separate splines', all(bcoef(:,:,:,:,:,:,j)==b1))
        end do
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,tr,ts,bcoef,k(1),k(2),k(3),k(4),k(5),k(6),iflag)
        call check('6d multi initialize_from_coefficients', iflag==0 .and. sm2%number_of_outputs()==nf)
        call sm2%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),fm2,iflag)
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),fm,iflag)
        call check('6d multi from coefficients gives the same values', all(fm2==fm))
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),f1,iflag)
        call check('6d single-output evaluate on multi returns 6003', iflag==6003)
        call sm%get_coefficients(b1,iflag)
        call check('6d single-output get_coefficients on multi returns 6003', iflag==6003)
        call sm%evaluate(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),fbad,iflag)
        call check('6d wrong size of f returns 6002', iflag==6002)
        call db6ink(x,n(1),y,n(2),z,n(3),q,n(4),r,n(5),s,n(6),fcn,k(1),k(2),k(3),k(4),k(5),k(6),1_ip,tx,ty,tz,tq,tr,ts,bcoef2(:,:,:,:,:,:,1:nf-1),iflag)
        call check('6d db6ink multi with mismatched outputs returns 6001', iflag==6001)
        call db6val(v(1),v(2),v(3),v(4),v(5),v(6),id(1),id(2),id(3),id(4),id(5),id(6),tx,ty,tz,tq,tr,ts,n(1),n(2),n(3),n(4),n(5),n(6),k(1),k(2),k(3),k(4),k(5),k(6),bcoef2(1:n(1)-1,:,:,:,:,:,:),fm2,iflag,inbvx,iloy,iloz,iloq,ilor,ilos,w5,w4,w3,w2,w1)
        call check('6d db6val multi with wrong bcoef shape returns 6004', iflag==6004)
        call sm2%initialize_from_coefficients(tx,ty,tz,tq,tr,ts,bcoef(:,:,:,:,:,:,1:0),k(1),k(2),k(3),k(4),k(5),k(6),iflag)
        call check('6d initialize_from_coefficients with no outputs returns 6005', &
                   iflag==6005 .and. sm2%number_of_outputs()==0)
        call sm%destroy()
        call check('6d destroy resets the number of outputs', sm%number_of_outputs()==0)
        deallocate(fcn, bcoef2, tx, ty, tz, tq, tr, ts, w5, w4, w3, w2, w1)
    end do

    end subroutine test_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  The 1D-only methods (integral, roots, etc.) require a single output.

    subroutine test_1d_methods()

    type(bspline_1d) :: sm  !! spline with multiple outputs
    type(bspline_1d) :: s1  !! spline with multiple outputs, but nf=1
    real(wp),dimension(10) :: x  !! data abscissae
    real(wp),dimension(10,2) :: fcn  !! data
    real(wp),dimension(:),allocatable :: r  !! roots
    real(wp) :: f, xm  !! results
    integer(ip) :: iflag  !! status flag

    x = [(real(iflag,wp), iflag=1,10)]
    fcn(:,1) = sin(x)
    fcn(:,2) = cos(x)
    call sm%initialize(x,fcn,4_ip,iflag)
    call sm%integral(1.0_wp,5.0_wp,f,iflag)
    call check('1d integral on multi returns 6003', iflag==6003)
    call sm%roots(0.0_wp,r,iflag)
    call check('1d roots on multi returns 6003', iflag==6003 .and. size(r)==0)
    call sm%minimum(xm,f,iflag)
    call check('1d minimum on multi returns 6003', iflag==6003)
    call sm%inverse(0.5_wp,xm,iflag)
    call check('1d inverse on multi returns 6003', iflag==6003)
    call check('6003 message', index(get_status_message(6003_ip),'single output')>0)

    ! a multiple-output spline with one output can use all the methods
    call s1%initialize(x,fcn(:,1:1),4_ip,iflag)
    call check('1d multi with one output: number_of_outputs', s1%number_of_outputs()==1)
    call s1%roots(0.0_wp,r,iflag)
    call check('1d multi with one output: roots', iflag==0 .and. size(r)==3)
    call s1%evaluate(2.5_wp,0_ip,f,iflag)
    call check('1d multi with one output: evaluate', iflag==0 .and. abs(f-sin(2.5_wp))<1.0e-2_wp)

    end subroutine test_1d_methods
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compare the time to evaluate `nf` outputs with one multiple-output
!  spline and with `nf` separate splines (printed only).

    subroutine timing()

    integer(ip),parameter :: nf = 10      !! number of outputs
    integer(ip),parameter :: neval = 20000 !! number of evaluations
    type(bspline_1d) :: m1          !! 1D multiple-output spline
    type(bspline_1d),dimension(nf) :: s1  !! 1D separate splines
    type(bspline_3d) :: m3          !! 3D multiple-output spline
    type(bspline_3d),dimension(nf) :: s3  !! 3D separate splines
    real(wp),dimension(20) :: x1     !! 1D data abscissae
    real(wp),dimension(20,nf) :: f1  !! 1D data
    real(wp),dimension(10) :: x3     !! 3D data abscissae
    real(wp),dimension(10,10,10,nf) :: f3  !! 3D data
    real(wp),dimension(nf) :: f      !! values
    real(wp),dimension(neval) :: u   !! evaluation points
    real(wp) :: t0, t1, t2, sum1, sum2  !! times and checksums
    integer(ip) :: i, j, iflag  !! counters and status flag

    x1 = grid(20_ip); call random_number(f1)
    x3 = grid(10_ip); call random_number(f3)
    call random_number(u)
    call m1%initialize(x1,f1,4_ip,iflag)
    call m3%initialize(x3,x3,x3,f3,4_ip,4_ip,4_ip,iflag)
    do j = 1, nf
        call s1(j)%initialize(x1,f1(:,j),4_ip,iflag)
        call s3(j)%initialize(x3,x3,x3,f3(:,:,:,j),4_ip,4_ip,4_ip,iflag)
    end do

    write(*,'(/A)') ' evaluation time for 10 outputs (multiple-output vs separate splines):'
    call cpu_time(t0)
    sum1 = 0.0_wp
    do i = 1, neval
        call m1%evaluate(x1(20)*u(i),0_ip,f,iflag); sum1 = sum1 + sum(f)
    end do
    call cpu_time(t1)
    sum2 = 0.0_wp
    do i = 1, neval
        do j = 1, nf
            call s1(j)%evaluate(x1(20)*u(i),0_ip,f(j),iflag); sum2 = sum2 + f(j)
        end do
    end do
    call cpu_time(t2)
    write(*,'(A,F8.4,A,F8.4,A,F6.2,A)') '   1D: ', t1-t0, ' s vs ', t2-t1, ' s (', (t2-t1)/max(t1-t0,1.0e-9_wp), 'x)'
    call check('timing 1d checksum', abs(sum1-sum2)<=1.0e-10_wp*abs(sum2))

    call cpu_time(t0)
    sum1 = 0.0_wp
    do i = 1, neval
        call m3%evaluate(x3(10)*u(i),x3(10)*u(neval+1-i),x3(10)*u(1+mod(i,neval)),0_ip,0_ip,0_ip,f,iflag)
        sum1 = sum1 + sum(f)
    end do
    call cpu_time(t1)
    sum2 = 0.0_wp
    do i = 1, neval
        do j = 1, nf
            call s3(j)%evaluate(x3(10)*u(i),x3(10)*u(neval+1-i),x3(10)*u(1+mod(i,neval)),0_ip,0_ip,0_ip,f(j),iflag)
            sum2 = sum2 + f(j)
        end do
    end do
    call cpu_time(t2)
    write(*,'(A,F8.4,A,F8.4,A,F6.2,A)') '   3D: ', t1-t0, ' s vs ', t2-t1, ' s (', (t2-t1)/max(t1-t0,1.0e-9_wp), 'x)'
    call check('timing 3d checksum', abs(sum1-sum2)<=1.0e-10_wp*abs(sum2))  ! (summed in a different order)

    end subroutine timing
!*****************************************************************************************

    end program bspline_multi_output_test
!*****************************************************************************************
