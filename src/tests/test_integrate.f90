!*****************************************************************************************
!> author: Jacob Williams
!  date: 9/6/2017
!
!  1D definite integral test using the bspline module.
!
!  This test case evaluates the integral:
!
!  $$ \int_{0}^{\pi} x^2 \sin (x) dx = 2 $$
!
!  using a B-Spline of the points:
!
!  $$ \left[ \sin(0^{\circ}), \sin(1^{\circ}), ..., \sin(180^{\circ}) \right] $$

    program bspline_integrate_test

    use bspline_module
    use bspline_kinds_module, only: wp

    implicit none

    real(wp),parameter :: pi      = acos(-1.0_wp)  !! \( \pi \)
    real(wp),parameter :: deg2rad = pi/180.0_wp    !! degrees to radians
    integer,parameter  :: iknot   = 0              !! automatically select the knots
    real(wp),parameter :: x1      = 0.0_wp         !! left endpoint
    real(wp),parameter :: x2      = pi             !! right endpoint
    integer,parameter  :: nx      = 181            !! number of points in x dimension
                                                   !! in original grid
    real(wp),parameter :: tol = 10.0_wp*epsilon(1.0_wp)  !! tolerance for [[db1fqad]]

    real(wp),dimension(:),allocatable :: tx     !! x knots
    integer                           :: kx     !! x bspline order
    real(wp),dimension(nx)            :: x      !! new grid x points
    real(wp),dimension(nx)            :: fcn    !! original grid
                                                !! function evaluations
    integer                           :: i      !! counter
    integer                           :: iflag  !! status flag
    real(wp)                          :: f      !! the evaluated integral
    integer                           :: imeth  !! method counter
    character(len=:),allocatable      :: meth   !! method string
    real(wp)                          :: f_true !! the true integral of
                                                !! the analytic function

    do imeth = 1,2  ! the two methods

        write(*,*) ''
        write(*,'(A8,1X,A5,1X,A30,1X,A30)') 'Method','Order','Integral','Error'

        do kx = 2, 15  ! spline orders

            if (allocated(tx)) deallocate(tx)
            allocate(tx(nx+kx))

            ! x^2 * sin(x) function evaluations for original grid:
            do i=1,nx
                x(i) = deg2rad*real(i-1,wp)
                if (imeth==1) then
                    fcn(i) = x(i)**2 * sin(x(i))
                else
                    ! for this one the x^2 will be handled by func
                    fcn(i) = sin(x(i))
                end if
            end do

            f_true = pi**2 - 4.0_wp

            ! initialize:
            ! [note we are overwriting fcn here with the b coeffs]
            call db1ink(x,nx,fcn,kx,iknot,tx,fcn,iflag)
            if (iflag/=0) error stop 'error calling db1ink'

            ! integrate:
            if (imeth==1) then
                if (kx>20) cycle
                meth = 'db1sqad'
                call db1sqad(tx,fcn,nx,kx,x1,x2,f,iflag)
            else
                meth = 'db1fqad'
                call db1fqad(test_function,tx,fcn,nx,kx,0,x1,x2,tol,f,iflag)
            end if

            ! display results:
            if (iflag/=0) then
                write(*,*) ''
                write(*,*) 'iflag: ',iflag
                write(*,*) 'error calling '//meth
                error stop
            else
                write(*,'(A8,1X,I5,1X,E30.16,1X,E30.16)') meth,kx,f,f-f_true
            end if

        end do
    end do
    write(*,*) ''

    contains

        function test_function(x) result(f)

        !! the function \( f(x) = x^2 \)
        !! to use for [[db1fqad]] test.

        implicit none

        real(wp),intent(in) :: x
        real(wp)            :: f  !! f(x)

        f = x*x

        end function test_function

    end program bspline_integrate_test
!*****************************************************************************************
