!*****************************************************************************************
!> author: Jacob Williams
!  date: 10/5/2015
!
!  2D data regridding using the bspline module.

    program bspline_regridding_test

    use bspline_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    integer,parameter :: kx = 4    !! x bspline order
    integer,parameter :: ky = 4    !! y bspline order
    integer,parameter :: idx = 0   !! [[db2val]] input
    integer,parameter :: idy = 0   !! [[db2val]] input
    integer,parameter :: nx = 6    !! number of points in x dimension in original grid
    integer,parameter :: ny = 5    !! number of points in y dimension in original grid
    real(wp),dimension(nx),parameter :: x = [0.0_wp,2.0_wp,4.0_wp,6.0_wp,8.0_wp,10.0_wp]  !! x points in original grid
    real(wp),dimension(ny),parameter :: y = [0.0_wp,2.0_wp,4.0_wp,6.0_wp,8.0_wp]          !! y points in original grid
    integer,parameter :: nx_new = 11  !! number of points in x dimension for new grid
    integer,parameter :: ny_new = 9   !! number of points in y dimension for new grid

    real(wp),dimension(nx_new)    :: x_new   !! new grid x points
    real(wp),dimension(ny_new)    :: y_new   !! new grid y points
    real(wp),dimension(nx_new,ny_new) :: fcn_new  !! new grid function evaluations
    real(wp),dimension(nx+kx) :: tx  !! x knots
    real(wp),dimension(ny+ky) :: ty  !! y knots
    real(wp),dimension(nx,ny) :: fcn_2d  !! original grid function evaluations
    real(wp) :: val,tru,err,errmax
    integer :: i,j
    integer :: iflag  !! status flag
    integer :: inbvx,inbvy,iloy

    !function evaluations for original grid:
    do i=1,nx
       do j=1,ny
           fcn_2d(i,j) = test_func(x(i),y(j))
       end do
    end do

    !display original data:
    write(*,*) '-----------------'
    write(*,*) '  INITIAL DATA:'
    write(*,*) '-----------------'
    write(*,'(A/,*(F12.6,1X))') 'x:', x
    write(*,*) ''
    write(*,'(A/,*(F12.6,1X))') 'y:', y
    write(*,*) ''
    write(*,'(A)') 'fcn(x,y):'
    do i=1,nx
        write(*,'(5F12.6)') fcn_2d(i,:)
    end do
    write(*,*) ''

    !regrid:

    inbvx = 1
    inbvy = 1
    iloy  = 1

    iflag = 0
    call db2ink(x,nx,y,ny,fcn_2d,kx,ky,tx,ty,fcn_2d,iflag)
    if (iflag/=1) error stop 'error calling db2ink'
    errmax = 0.0_wp
    do i=1,nx_new
        x_new(i) = real(i-1,wp)
        do j=1,ny_new
            y_new(j) = real(j-1,wp)
            call db2val(x_new(i),y_new(j),idx,idy,tx,ty,nx,ny,kx,ky,fcn_2d,val,iflag,&
                        inbvx,inbvy,iloy)
            if (iflag/=0) error stop 'error calling db2val'
            fcn_new(i,j) = val
            tru    = test_func(x_new(i),y_new(j))  !truth value
            err    = abs(tru-val)
            errmax = max(err,errmax)
        end do
    end do

    !display new grid:
    write(*,*) '-----------------'
    write(*,*) '  NEW GRID:'
    write(*,*) '-----------------'
    write(*,'(A/,*(F12.6,1X))') 'x:', x_new
    write(*,*) ''
    write(*,'(A/,*(F12.6,1X))') 'y:', y_new
    write(*,*) ''
    write(*,'(A)') 'fcn(x,y):'
    do i=1,nx_new
        write(*,'(11F12.6)') fcn_new(i,:)
    end do
    write(*,*) ''
    write(*,*) ' max error:', errmax
    write(*,*) ''

    contains

        function test_func(x,y) result(f)
        !! 2d test function

        implicit none

        real(wp) :: f
        real(wp),intent(in) :: x,y

        real(wp),parameter :: deg2rad = acos(-1.0_wp)/180.0_wp  !! degrees to radians conversion factor

        f = sin(deg2rad*(x+y))

        end function test_func

    end program bspline_regridding_test
