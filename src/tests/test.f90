!*****************************************************************************************
!>
!  Units test for 1d-6d tensor product b-spline interpolation.

    program bspline_test

    use bspline_module
    use bspline_kinds_module, only: wp

    implicit none

    integer,parameter :: nx = 6     !! number of points in x
    integer,parameter :: ny = 6     !! number of points in y
    integer,parameter :: nz = 6     !! number of points in z
    integer,parameter :: nq = 6     !! number of points in q
    integer,parameter :: nr = 6     !! number of points in r
    integer,parameter :: ns = 6     !! number of points in s

    integer,parameter :: kx = 4     !! order in x
    integer,parameter :: ky = 4     !! order in y
    integer,parameter :: kz = 4     !! order in z
    integer,parameter :: kq = 4     !! order in q
    integer,parameter :: kr = 4     !! order in r
    integer,parameter :: ks = 4     !! order in s

    integer,parameter :: iknot = 0  !! automatically select the knots

    real(wp) :: x(nx),y(ny),z(nz),q(nq),r(nr),s(ns)
    real(wp) :: tx(nx+kx),ty(ny+ky),tz(nz+kz),tq(nq+kq),tr(nr+kr),ts(ns+ks)
    real(wp) :: fcn_1d(nx)
    real(wp) :: fcn_2d(nx,ny)
    real(wp) :: fcn_3d(nx,ny,nz)
    real(wp) :: fcn_4d(nx,ny,nz,nq)
    real(wp) :: fcn_5d(nx,ny,nz,nq,nr)
    real(wp) :: fcn_6d(nx,ny,nz,nq,nr,ns)

    real(wp) :: tol
    real(wp),dimension(6) :: val,tru,err,errmax
    logical :: fail
    integer :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids
    integer,dimension(6) :: iflag
    integer :: inbvx,inbvy,inbvz,inbvq,inbvr,inbvs
    integer :: iloy,iloz,iloq,ilor,ilos

    fail = .false.
    tol = 1.0e-14_wp
    idx = 0
    idy = 0
    idz = 0
    idq = 0
    idr = 0
    ids = 0

     do i=1,nx
        x(i) = dble(i-1)/dble(nx-1)
     end do
     do j=1,ny
        y(j) = dble(j-1)/dble(ny-1)
     end do
     do k=1,nz
        z(k) = dble(k-1)/dble(nz-1)
     end do
     do l=1,nq
        q(l) = dble(l-1)/dble(nq-1)
     end do
     do m=1,nr
        r(m) = dble(m-1)/dble(nr-1)
     end do
     do n=1,ns
        s(n) = dble(n-1)/dble(ns-1)
     end do
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

    !have to set these before the first evaluate call:
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

    ! initialize
    call db1ink(x,nx,fcn_1d,kx,iknot,tx,fcn_1d,iflag(1))
    call db2ink(x,nx,y,ny,fcn_2d,kx,ky,iknot,tx,ty,fcn_2d,iflag(2))
    call db3ink(x,nx,y,ny,z,nz,fcn_3d,kx,ky,kz,iknot,tx,ty,tz,fcn_3d,iflag(3))
    call db4ink(x,nx,y,ny,z,nz,q,nq,fcn_4d,kx,ky,kz,kq,iknot,tx,ty,tz,tq,fcn_4d,iflag(4))
    call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn_5d,kx,ky,kz,kq,kr,iknot,tx,ty,tz,tq,tr,fcn_5d,iflag(5))
    call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn_6d,kx,ky,kz,kq,kr,ks,iknot,tx,ty,tz,tq,tr,ts,fcn_6d,iflag(6))

    if (any(iflag/=0)) then
        do i=1,6
            if (iflag(i)/=0) then
                write(*,*) 'Error initializing ',i,'D spline: '//get_status_message(iflag(i))
            end if
        end do
    end if

    ! compute max error at interpolation points

     errmax = 0.0_wp
     do i=1,nx
                        call db1val(x(i),idx,&
                                            tx,nx,kx,fcn_1d,val(1),iflag(1),inbvx)
                        tru(1)    = f1(x(i))
                        err(1)    = abs(tru(1)-val(1))
                        errmax(1) = max(err(1),errmax(1))
        do j=1,ny
                        call db2val(x(i),y(j),idx,idy,&
                                            tx,ty,nx,ny,kx,ky,fcn_2d,val(2),iflag(2),&
                                            inbvx,inbvy,iloy)
                        tru(2)    = f2(x(i),y(j))
                        err(2)    = abs(tru(2)-val(2))
                        errmax(2) = max(err(2),errmax(2))
           do k=1,nz
                        call db3val(x(i),y(j),z(k),idx,idy,idz,&
                                            tx,ty,tz,nx,ny,nz,kx,ky,kz,fcn_3d,val(3),iflag(3),&
                                            inbvx,inbvy,inbvz,iloy,iloz)
                        tru(3)    = f3(x(i),y(j),z(k))
                        err(3)    = abs(tru(3)-val(3))
                        errmax(3) = max(err(3),errmax(3))
              do l=1,nq
                        call db4val(x(i),y(j),z(k),q(l),idx,idy,idz,idq,&
                                            tx,ty,tz,tq,nx,ny,nz,nq,kx,ky,kz,kq,fcn_4d,val(4),iflag(4),&
                                            inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq)
                        tru(4)    = f4(x(i),y(j),z(k),q(l))
                        err(4)    = abs(tru(4)-val(4))
                        errmax(4) = max(err(4),errmax(4))
                do m=1,nr
                        call db5val(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,&
                                            tx,ty,tz,tq,tr,nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,fcn_5d,val(5),iflag(5),&
                                            inbvx,inbvy,inbvz,inbvq,inbvr,iloy,iloz,iloq,ilor)
                        tru(5)    = f5(x(i),y(j),z(k),q(l),r(m))
                        err(5)    = abs(tru(5)-val(5))
                        errmax(5) = max(err(5),errmax(5))
                    do n=1,ns
                        call db6val(x(i),y(j),z(k),q(l),r(m),s(n),idx,idy,idz,idq,idr,ids,&
                                            tx,ty,tz,tq,tr,ts,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,fcn_6d,val(6),iflag(6),&
                                            inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,iloy,iloz,iloq,ilor,ilos)
                        tru(6)    = f6(x(i),y(j),z(k),q(l),r(m),s(n))
                        err(6)    = abs(tru(6)-val(6))
                        errmax(6) = max(err(6),errmax(6))
                    end do
                end do
              end do
           end do
        end do
     end do

    ! check max error against tolerance
    do i=1,6
        write(*,*) i,'D: max error:', errmax(i)
        if (errmax(i) >= tol) then
            write(*,*)  ' ** test failed ** '
        else
            write(*,*)  ' ** test passed ** '
        end if
        write(*,*) ''
    end do

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        f1 = 0.5_wp * (x*exp(-x) + sin(x) )
        end function f1

        real(wp) function f2(x,y) !! 2d test function
        implicit none
        real(wp) x,y,piov2
        piov2 = 2.0_wp * atan(1.0_wp)
        f2 = 0.5_wp * (y*exp(-x) + sin(piov2*y) )
        end function f2

        real(wp) function f3 (x,y,z) !! 3d test function
        implicit none
        real(wp) x,y,z,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f3 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) )
        end function f3

        real(wp) function f4 (x,y,z,q) !! 4d test function
        implicit none
        real(wp) x,y,z,q,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f4 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q )
        end function f4

        real(wp) function f5 (x,y,z,q,r) !! 5d test function
        implicit none
        real(wp) x,y,z,q,r,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f5 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q*r )
        end function f5

        real(wp) function f6 (x,y,z,q,r,s) !! 6d test function
        implicit none
        real(wp) x,y,z,q,r,s,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f6 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q*r + 2.0_wp*s )
        end function f6

    end program bspline_test
