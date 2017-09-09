!*****************************************************************************************
!> author: Jacob Williams
!  date: 9/6/2017
!
!  1D definite integral test using the bspline module.
!
!  This test case evaluates the integral:
!
!  $$ \int_{0}^{\pi} \sin (x) dx = 2 $$
!
!  using a B-Spline of the points:
!
!  $$ \left[ \sin(0^{\circ}), \sin(1^{\circ}), ..., \sin(180^{\circ}) \right] $$

    program bspline_integrate_test

    use bspline_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    real(wp),parameter :: pi      = acos(-1.0_wp)  !! \( \pi \)
    real(wp),parameter :: deg2rad = pi/180.0_wp    !! degrees to radians
    integer,parameter  :: iknot   = 0              !! automatically select the knots
    real(wp),parameter :: x1      = 0.0_wp         !! left endpoint
    real(wp),parameter :: x2      = pi             !! right endpoint
    integer,parameter  :: nx      = 181            !! number of points in x dimension
                                                   !! in original grid
    real(wp),parameter :: tol     = 1.0e-12_wp     !! tolerance for [[db1fqad]]

    real(wp),dimension(:),allocatable :: tx     !! x knots
    integer                           :: kx     !! x bspline order
    real(wp),dimension(nx)            :: x      !! new grid x points
    real(wp),dimension(nx)            :: fcn    !! original grid
                                                !! function evaluations
    integer                           :: i      !! counter
    integer                           :: iflag  !! status flag
    real(wp)                          :: f      !! the evaluated integral
                                                !! (should be close to 2)

    write(*,*) ''
    write(*,'(A8,1X,A5,1X,A30,1X,A30)') 'Method','Order','Integral','Error'

    do kx = 2, 8

        if (allocated(tx)) deallocate(tx)
        allocate(tx(nx+kx))

        !function evaluations for original grid:
        do i=1,nx
            x(i)   = deg2rad*real(i-1,wp)
            fcn(i) = sin(x(i))
        end do

        ! initialize:
        call db1ink(x,nx,fcn,kx,iknot,tx,fcn,iflag)
        if (iflag/=0) error stop 'error calling db1ink'

        ! now integrate:
        call db1sqad(tx,fcn,nx,kx,x1,x2,f,iflag)

        ! display results:
        if (iflag/=0) then
            write(*,*) ''
            write(*,*) 'iflag: ',iflag
            error stop 'error calling db1qad'
        else
            write(*,'(A8,1X,I5,1X,E30.16,1X,E30.16)') 'db1qad',kx,f,f-2.0_wp
        end if

        ! integrate using adaptive version:
        call db1fqad(test_function,tx,fcn,nx,kx,0,x1,x2,tol,f,iflag)

        ! display results:
        if (iflag/=0) then
            write(*,*) ''
            write(*,*) 'iflag: ',iflag
            error stop 'error calling db1fqad'
        else
            write(*,'(A8,1X,I5,1X,E30.16,1X,E30.16)') 'db1fqad',kx,f,f-4.0_wp
        end if


    end do
    write(*,*) ''

    contains

        function test_function(x) result(f)

        implicit none

        real(wp),intent(in) :: x
        real(wp)            :: f  !! f(x)

        f = 2.0_wp

        end function test_function

    end program bspline_integrate_test
!*****************************************************************************************
