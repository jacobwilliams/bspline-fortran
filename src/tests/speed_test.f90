!*****************************************************************************************
!
!> Speed test for 1d-6d tensor product b-spline interpolation (object-oriented version).

    program bspline_speed_test
    
    use bspline_oo_module
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none
    
    integer,parameter :: nx = 6    !number of points
    integer,parameter :: ny = 8
    integer,parameter :: nz = 10
    integer,parameter :: nq = 9
    integer,parameter :: nr = 7
    integer,parameter :: ns = 8
    
    integer,parameter :: kx = 2    !order
    integer,parameter :: ky = 3
    integer,parameter :: kz = 4
    integer,parameter :: kq = 3
    integer,parameter :: kr = 2
    integer,parameter :: ks = 3
            
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
    integer  :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids,iflag,n_cases
    real :: tstart, tend
        
    idx = 0
    idy = 0
    idz = 0
    idq = 0
    idr = 0
    ids = 0

    x = [(dble(i), i=1,nx)]
    y = [(dble(i), i=1,ny)]
    z = [(dble(i), i=1,nz)]
    q = [(dble(i), i=1,nq)]
    r = [(dble(i), i=1,nr)]
    s = [(dble(i), i=1,ns)]

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
    sumval = 0.0_wp
    n_cases = nx*ny*nz*nq*nr*ns
    
    call cpu_time(tstart)
    do i=1,nx
        call s1%evaluate(x(i),idx,val,iflag)
        sumval = sumval + val
        do j=1,ny
           call s2%evaluate(x(i),y(j),idx,idy,val,iflag)
           sumval = sumval + val
           do k=1,nz
                call s3%evaluate(x(i),y(j),z(k),idx,idy,idz,val,iflag)
                sumval = sumval + val
                  do l=1,nq
                    call s4%evaluate(x(i),y(j),z(k),q(l),idx,idy,idz,idq,val,iflag)
                    sumval = sumval + val
                    do m=1,nr
                        call s5%evaluate(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,val,iflag)
                        sumval = sumval + val
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
    
    write(*,*) 'result         :', sumval
    write(*,*) 'number of cases:', n_cases
    write(*,*) 'cases/sec      :', n_cases/(tend-tstart)
 
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
                  
    end program bspline_speed_test