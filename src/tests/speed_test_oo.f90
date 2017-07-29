!*****************************************************************************************
!>
!  Speed test for 1d-6d tensor product b-spline interpolation (object-oriented version).
!
!### Results
!  ![Plot of results](https://raw.githubusercontent.com/jacobwilliams/bspline-fortran/master/src/tests/speed_test_oo.png)

    program bspline_speed_test_oo

    use bspline_module
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use pyplot_module

    implicit none

    integer,parameter :: nx = 8   !number of points
    integer,parameter :: ny = 8
    integer,parameter :: nz = 8
    integer,parameter :: nq = 8
    integer,parameter :: nr = 8
    integer,parameter :: ns = 8

    integer,parameter :: kx = 4    !order
    integer,parameter :: ky = 4
    integer,parameter :: kz = 4
    integer,parameter :: kq = 4
    integer,parameter :: kr = 4
    integer,parameter :: ks = 4

    integer,parameter :: n_cases = nx*ny*nz*nq*nr*ns

    real(wp) :: x(nx),y(ny),z(nz),q(nq),r(nr),s(ns)
    real(wp) :: fcn_1d(nx)
    real(wp) :: fcn_2d(nx,ny)
    real(wp) :: fcn_3d(nx,ny,nz)
    real(wp) :: fcn_4d(nx,ny,nz,nq)
    real(wp) :: fcn_5d(nx,ny,nz,nq,nr)
    real(wp) :: fcn_6d(nx,ny,nz,nq,nr,ns)

    type(bspline_1d) :: s1
    type(bspline_2d) :: s2
    type(bspline_3d) :: s3
    type(bspline_4d) :: s4
    type(bspline_5d) :: s5
    type(bspline_6d) :: s6

    real(wp) :: val,sumval
    integer  :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids,iflag
    real :: tstart, tend
    type(pyplot) :: plt
    real(wp),dimension(6) :: cases_per_sec
    integer :: istat  !! pyplot-fortran status flag

    idx = 0
    idy = 0
    idz = 0
    idq = 0
    idr = 0
    ids = 0

    x = [ (dble(i), i=1,nx) ]
    y = [ (dble(i), i=1,ny) ]
    z = [ (dble(i), i=1,nz) ]
    q = [ (dble(i), i=1,nq) ]
    r = [ (dble(i), i=1,nr) ]
    s = [ (dble(i), i=1,ns) ]

    !evaluate the functions:

    do i=1,nx
        fcn_1d(i) = f1(x(i))
        do j=1,ny
            fcn_2d(i,j) = f2(x(i),y(j))
            do k=1,nz
                fcn_3d(i,j,k) = f3(x(i),y(j),z(k))
                do l=1,nq
                    fcn_4d(i,j,k,l) = f4(x(i),y(j),z(k),q(l))
                    do m=1,nr
                        fcn_5d(i,j,k,l,m) = f5(x(i),y(j),z(k),q(l),r(m))
                        do n=1,ns
                            fcn_6d(i,j,k,l,m,n) = f6(x(i),y(j),z(k),q(l),r(m),s(n))
                        end do
                    end do
                end do
            end do
        end do
    end do

    !initialize using the constructors:

    s1 = bspline_1d(x,fcn_1d,kx)                         ; if (.not. s1%status_ok()) error stop 'error initializing s1'
    s2 = bspline_2d(x,y,fcn_2d,kx,ky)                    ; if (.not. s1%status_ok()) error stop 'error initializing s2'
    s3 = bspline_3d(x,y,z,fcn_3d,kx,ky,kz)               ; if (.not. s1%status_ok()) error stop 'error initializing s3'
    s4 = bspline_4d(x,y,z,q,fcn_4d,kx,ky,kz,kq)          ; if (.not. s1%status_ok()) error stop 'error initializing s4'
    s5 = bspline_5d(x,y,z,q,r,fcn_5d,kx,ky,kz,kq,kr)     ; if (.not. s1%status_ok()) error stop 'error initializing s5'
    s6 = bspline_6d(x,y,z,q,r,s,fcn_6d,kx,ky,kz,kq,kr,ks); if (.not. s1%status_ok()) error stop 'error initializing s6'

    ! evaluate the interpolants:
    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call s1%evaluate(x(i),idx,val,iflag)
                            sumval = sumval + val
                        end do
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '1D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(1) = n_cases/(tend-tstart)

    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                           call s2%evaluate(x(i),y(j),idx,idy,val,iflag)
                           sumval = sumval + val
                       end do
                   end do
               end do
           end do
       end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '2D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(2) = n_cases/(tend-tstart)

    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call s3%evaluate(x(i),y(j),z(k),idx,idy,idz,val,iflag)
                            sumval = sumval + val
                        end do
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '3D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(3) = n_cases/(tend-tstart)

    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call s4%evaluate(x(i),y(j),z(k),q(l),idx,idy,idz,idq,val,iflag)
                            sumval = sumval + val
                        end do
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '4D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(4) = n_cases/(tend-tstart)

    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call s5%evaluate(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,val,iflag)
                            sumval = sumval + val
                        end do
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '5D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(5) = n_cases/(tend-tstart)

    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call s6%evaluate(x(i),y(j),z(k),q(l),r(m),s(n),idx,idy,idz,idq,idr,ids,val,iflag)
                            sumval = sumval + val
                        end do
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(tend)
    write(*,*) ''
    write(*,*) '6D'
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
    cases_per_sec(6) = n_cases/(tend-tstart)

    !plot results in bar chart:
    call plt%initialize(grid=.false.,xlabel='Dimension',ylabel='Cases Per Second',&
                        title='Speed Test (Object-Oriented Interface)',legend=.false.,&
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20)
    call plt%add_bar(left=real([1,2,3,4,5,6],wp),height=cases_per_sec,label='Speed test runs',&
                        yscale='log',align='center',color='r',istat=istat)
    call plt%savefig('speed_test_oo.png',istat=istat)

    contains

        function f1(x) result(f) !! 1d test function
        implicit none
        real(wp),intent(in) :: x
        real(wp) :: f
        f = x**1.1_wp + x**1.2_wp + x**1.3_wp + x**1.4_wp + x**1.5_wp + x**1.6_wp
        end function f1

        function f2(x,y) result(f) !! 2d test function
        implicit none
        real(wp),intent(in) :: x,y
        real(wp) :: f
        f = x**1.1_wp + y**1.2_wp + x**1.3_wp + y**1.4_wp + x**1.5_wp + y**1.6_wp
        end function f2

        function f3 (x,y,z) result(f) !! 3d test function
        implicit none
        real(wp),intent(in) :: x,y,z
        real(wp) :: f
        f = x**1.1_wp + y**1.2_wp + z**1.3_wp + x**1.4_wp + y**1.5_wp + z**1.6_wp
        end function f3

        function f4 (x,y,z,q) result(f) !! 4d test function
        implicit none
        real(wp),intent(in) :: x,y,z,q
        real(wp) :: f
        f = x**1.1_wp + y**1.2_wp + z**1.3_wp + q**1.4_wp + x**1.5_wp + y**1.6_wp
        end function f4

        function f5 (x,y,z,q,r) result(f) !! 5d test function
        implicit none
        real(wp),intent(in) :: x,y,z,q,r
        real(wp) :: f
        f = x**1.1_wp + y**1.2_wp + z**1.3_wp + q**1.4_wp + r**1.5_wp + x**1.6_wp
        end function f5

        function f6 (x,y,z,q,r,s) result(f) !! 6d test function
        implicit none
        real(wp),intent(in) :: x,y,z,q,r,s
        real(wp) :: f
        f = x**1.1_wp + y**1.2_wp + z**1.3_wp + q**1.4_wp + r**1.5_wp + s**1.6_wp
        end function f6

    end program bspline_speed_test_oo
