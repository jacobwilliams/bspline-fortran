!*****************************************************************************************
!>
!  A test of very large data sets for the 3D case.

    program bspline_stack_size_oo_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer :: i !! counter

    ! number of data points to test (n x n x n):
    integer,parameter :: n_start  = 90     !! initial size
    integer,parameter :: n_stop   = 1000   !! final size
    integer,parameter :: n_step   = 20     !! step size

    write(*,*) ''
    write(*,'(A10,1X,A20,1X,A30)') 'num points','size (bytes)', 'error'
    do i = n_start, n_stop, n_step
        call run_test(i)
    end do
    write(*,*) ''

    contains

    subroutine run_test(n)  !! run test for n x n x n data set

    implicit none

    integer,intent(in) :: n

    integer(ip),parameter :: kx = 4    !order
    integer(ip),parameter :: ky = 4
    integer(ip),parameter :: kz = 4
    real(wp),parameter :: tol = 1.0e-14_wp

    type(bspline_3d) :: s3
    integer(ip) :: nx,ny,nz
    real(wp),allocatable :: x(:),y(:),z(:)
    real(wp),allocatable  :: fcn_3d(:,:,:)
    real(wp) :: val,tru,err,errmax
    logical  :: fail
    integer(ip) :: i,j,k,idx,idy,idz
    integer(ip) :: iflag

    nx = n; ny = n; nz = n
    idx = 0; idy = 0; idz = 0
    fail = .false.

    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))
    allocate(fcn_3d(nx,ny,nz))

     do concurrent (i=1:nx)
        x(i) = real(i-1,wp)/real(nx-1,wp)
     end do
     do concurrent (j=1:ny)
        y(j) = real(j-1,wp)/real(ny-1,wp)
     end do
     do concurrent (k=1:nz)
        z(k) = real(k-1,wp)/real(nz-1,wp)
     end do

     do i=1,nx
        do j=1,ny
           do k=1,nz
              fcn_3d(i,j,k) = f3(x(i),y(j),z(k))
           end do
        end do
     end do

     !initialize:
     call s3%initialize(x,y,z,fcn_3d,kx,ky,kz,iflag)

     ! free up memory (don't need anymore):
     deallocate(fcn_3d)

     if (iflag/=0) then
         write(*,*) 'Error initializing spline: '//get_status_message(iflag)
         error stop
     end if

    ! compute max error at interpolation points
     errmax = 0.0_wp
     do i=1,nx
        do j=1,ny
           do k=1,nz
                call s3%evaluate(x(i),y(j),z(k),idx,idy,idz,val,iflag)
                tru    = f3(x(i),y(j),z(k))
                err    = abs(tru-val)
                errmax = max(err,errmax)
           end do
        end do
     end do

    ! check max error against tolerance
    write(*,'(I10,1X,I20,1X,E30.16)') n, s3%size_of()*8, errmax    ! n, size, error
    if (errmax >= tol) error stop  ' ** test failed ** '

    end subroutine run_test

    pure real(wp) function f3 (x,y,z) !! 3d test function
    implicit none
    real(wp),intent(in) :: x,y,z
    real(wp) :: piov2
    piov2 = 2.0_wp*atan(1.0_wp)
    f3 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) )
    end function f3

    end program bspline_stack_size_oo_test
