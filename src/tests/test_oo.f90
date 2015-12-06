!*****************************************************************************************
!
!> Units test for 2d-6d tensor product b-spline interpolation (object-oriented version).

    program bspline_oo_test

    use bspline_oo_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    integer,parameter :: nx = 6    !number of points
    integer,parameter :: ny = 6
    integer,parameter :: nz = 6
    integer,parameter :: nq = 6
    integer,parameter :: nr = 6
    integer,parameter :: ns = 6

    integer,parameter :: kx = 4    !order
    integer,parameter :: ky = 4
    integer,parameter :: kz = 4
    integer,parameter :: kq = 4
    integer,parameter :: kr = 4
    integer,parameter :: ks = 4

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

    real(wp) :: tol
    real(wp),dimension(6) :: val,tru,err,errmax
    logical  :: fail
    integer  :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids,iflag

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

     !initialize:

     call s1%initialize(x,fcn_1d,kx,iflag)
     call s2%initialize(x,y,fcn_2d,kx,ky,iflag)
     call s3%initialize(x,y,z,fcn_3d,kx,ky,kz,iflag)
     call s4%initialize(x,y,z,q,fcn_4d,kx,ky,kz,kq,iflag)
     call s5%initialize(x,y,z,q,r,fcn_5d,kx,ky,kz,kq,kr,iflag)
     call s6%initialize(x,y,z,q,r,s,fcn_6d,kx,ky,kz,kq,kr,ks,iflag)

    ! compute max error at interpolation points

     errmax = 0.0_wp
     do i=1,nx
                        call s1%evaluate(x(i),idx,val(1),iflag)
                        tru(1)    = f1(x(i))
                        err(1)    = abs(tru(1)-val(1))
                        errmax(1) = max(err(1),errmax(1))
        do j=1,ny
                        call s2%evaluate(x(i),y(j),idx,idy,val(2),iflag)
                        tru(2)    = f2(x(i),y(j))
                        err(2)    = abs(tru(2)-val(2))
                        errmax(2) = max(err(2),errmax(2))
           do k=1,nz
                        call s3%evaluate(x(i),y(j),z(k),idx,idy,idz,val(3),iflag)
                        tru(3)    = f3(x(i),y(j),z(k))
                        err(3)    = abs(tru(3)-val(3))
                        errmax(3) = max(err(3),errmax(3))
              do l=1,nq
                        call s4%evaluate(x(i),y(j),z(k),q(l),idx,idy,idz,idq,val(4),iflag)
                        tru(4)    = f4(x(i),y(j),z(k),q(l))
                        err(4)    = abs(tru(4)-val(4))
                        errmax(4) = max(err(4),errmax(4))
                do m=1,nr
                        call s5%evaluate(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,val(5),iflag)
                        tru(5)    = f5(x(i),y(j),z(k),q(l),r(m))
                        err(5)    = abs(tru(5)-val(5))
                        errmax(5) = max(err(5),errmax(5))
                    do n=1,ns
                        call s6%evaluate(x(i),y(j),z(k),q(l),r(m),s(n),idx,idy,idz,idq,idr,ids,val(6),iflag)
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

    end program bspline_oo_test
