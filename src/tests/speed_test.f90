!*****************************************************************************************
!>
!  Speed test for 1d-6d tensor product b-spline interpolation (subroutine version).
!
!### Results
!  ![Plot of results](https://raw.githubusercontent.com/jacobwilliams/bspline-fortran/master/src/tests/speed_test.png)

    program bspline_speed_test

    use bspline_module
    use bspline_kinds_module, only: wp
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
    real(wp),dimension(nx) :: fcn_1d, bcoef_1d
    real(wp),dimension(nx,ny) :: fcn_2d, bcoef_2d
    real(wp),dimension(nx,ny,nz) :: fcn_3d, bcoef_3d
    real(wp),dimension(nx,ny,nz,nq) :: fcn_4d, bcoef_4d
    real(wp),dimension(nx,ny,nz,nq,nr) :: fcn_5d, bcoef_5d
    real(wp),dimension(nx,ny,nz,nq,nr,ns) :: fcn_6d, bcoef_6d
    real(wp) :: tx(nx+kx),ty(ny+ky),tz(nz+kz),tq(nq+kq),tr(nr+kr),ts(ns+ks)

    real(wp) :: val,sumval
    integer  :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids,iknot,iflag
    real :: tstart, tend
    type(pyplot) :: plt
    real(wp),dimension(6) :: cases_per_sec
    integer :: inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,iloy,iloz,iloq,ilor,ilos
    integer :: istat  !! pyplot-fortran status flag

    inbvx = 1
    inbvy = 1
    inbvz = 1
    inbvq = 1
    inbvr = 1
    inbvs = 1
    iloy  = 1
    iloz  = 1
    iloq  = 1
    ilor  = 1
    ilos  = 1

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

    iknot = 0 !auto-compute the knots

    !initialize:

    call db1ink(x,nx,fcn_1d,kx,iknot,tx,bcoef_1d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 1d spline.'
        end if
    call db2ink(x,nx,y,ny,fcn_2d,kx,ky,iknot,tx,ty,bcoef_2d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 2d spline.'
        end if
    call db3ink(x,nx,y,ny,z,nz,fcn_3d,kx,ky,kz,iknot,tx,ty,tz,bcoef_3d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 3d spline.'
        end if
    call db4ink(x,nx,y,ny,z,nz,q,nq,fcn_4d,kx,ky,kz,kq,iknot,tx,ty,tz,tq,bcoef_4d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 4d spline.'
        end if
    call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn_5d,kx,ky,kz,kq,kr,iknot,tx,ty,tz,tq,tr,bcoef_5d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 5d spline.'
        end if
    call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn_6d,kx,ky,kz,kq,kr,ks,iknot,tx,ty,tz,tq,tr,ts,bcoef_6d,iflag)
        if (iflag/=0) then
            write(*,*) 'iflag=',iflag
            error stop 'error initializing 6d spline.'
        end if

    ! evaluate the interpolants:
    sumval = 0.0_wp
    call cpu_time(tstart)
    do i=1,nx
        do j=1,ny
           do k=1,nz
                do l=1,nq
                    do m=1,nr
                        do n=1,ns
                            call db1val(x(i),idx,tx,nx,kx,bcoef_1d,val,iflag,inbvx)
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
                            call db2val(x(i),y(j),idx,idy,tx,ty,&
                                        nx,ny,kx,ky,bcoef_2d,val,iflag,&
                                        inbvx,inbvy,iloy)
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
                            call db3val(x(i),y(j),z(k),idx,idy,idz,tx,ty,tz,&
                                        nx,ny,nz,kx,ky,kz,bcoef_3d,val,iflag,&
                                        inbvx,inbvy,inbvz,iloy,iloz)
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
                            call db4val(x(i),y(j),z(k),q(l),idx,idy,idz,idq,tx,ty,tz,tq,&
                                        nx,ny,nz,nq,kx,ky,kz,kq,bcoef_4d,val,iflag,&
                                        inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq)
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
                            call db5val(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,tx,ty,tz,tq,tr,&
                                        nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,bcoef_5d,val,iflag,&
                                        inbvx,inbvy,inbvz,inbvq,inbvr,iloy,iloz,iloq,ilor)
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
                            call db6val(x(i),y(j),z(k),q(l),r(m),s(n),idx,idy,idz,idq,idr,ids,&
                                        tx,ty,tz,tq,tr,ts,&
                                        nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,bcoef_6d,val,iflag,&
                                        inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,iloy,iloz,iloq,ilor,ilos)
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
                        title='Speed Test (Subroutine Interface)',legend=.false.,&
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20)
    call plt%add_bar(x=real([1,2,3,4,5,6],wp),height=cases_per_sec,label='Speed test runs',&
                        yscale='log',align='center',color='r',istat=istat)
    call plt%savefig('speed_test.png',istat=istat)

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

    end program bspline_speed_test
