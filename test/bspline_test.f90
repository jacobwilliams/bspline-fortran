!*****************************************************************************************
!>
!  Units test for 1d-6d tensor product b-spline interpolation.

    program bspline_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: nx = 6     !! number of points in x
    integer(ip),parameter :: ny = 6     !! number of points in y
    integer(ip),parameter :: nz = 6     !! number of points in z
    integer(ip),parameter :: nq = 6     !! number of points in q
    integer(ip),parameter :: nr = 6     !! number of points in r
    integer(ip),parameter :: ns = 6     !! number of points in s

    integer(ip),parameter :: kx = 4     !! order in x
    integer(ip),parameter :: ky = 4     !! order in y
    integer(ip),parameter :: kz = 4     !! order in z
    integer(ip),parameter :: kq = 4     !! order in q
    integer(ip),parameter :: kr = 4     !! order in r
    integer(ip),parameter :: ks = 4     !! order in s

    integer(ip),parameter :: iknot = 0  !! automatically select the knots

    real(wp) :: x(nx),y(ny),z(nz),q(nq),r(nr),s(ns)
    real(wp) :: tx(nx+kx),ty(ny+ky),tz(nz+kz),tq(nq+kq),tr(nr+kr),ts(ns+ks)
    real(wp) :: fcn_1d(nx)               , bcoef_1d(nx)
    real(wp) :: fcn_2d(nx,ny)            , bcoef_2d(nx,ny)
    real(wp) :: fcn_3d(nx,ny,nz)         , bcoef_3d(nx,ny,nz)
    real(wp) :: fcn_4d(nx,ny,nz,nq)      , bcoef_4d(nx,ny,nz,nq)
    real(wp) :: fcn_5d(nx,ny,nz,nq,nr)   , bcoef_5d(nx,ny,nz,nq,nr)
    real(wp) :: fcn_6d(nx,ny,nz,nq,nr,ns), bcoef_6d(nx,ny,nz,nq,nr,ns)

    real(wp),dimension(3*kx)                     :: w1_1d
    real(wp),dimension(ky)                       :: w1_2d
    real(wp),dimension(3*max(kx,ky))             :: w2_2d
    real(wp),dimension(ky,kz)                    :: w1_3d
    real(wp),dimension(kz)                       :: w2_3d
    real(wp),dimension(3*max(kx,ky,kz))          :: w3_3d
    real(wp),dimension(ky,kz,kq)                 :: w1_4d
    real(wp),dimension(kz,kq)                    :: w2_4d
    real(wp),dimension(kq)                       :: w3_4d
    real(wp),dimension(3*max(kx,ky,kz,kq))       :: w4_4d
    real(wp),dimension(ky,kz,kq,kr)              :: w1_5d
    real(wp),dimension(kz,kq,kr)                 :: w2_5d
    real(wp),dimension(kq,kr)                    :: w3_5d
    real(wp),dimension(kr)                       :: w4_5d
    real(wp),dimension(3*max(kx,ky,kz,kq,kr))    :: w5_5d
    real(wp),dimension(ky,kz,kq,kr,ks)           :: w1_6d
    real(wp),dimension(kz,kq,kr,ks)              :: w2_6d
    real(wp),dimension(kq,kr,ks)                 :: w3_6d
    real(wp),dimension(kr,ks)                    :: w4_6d
    real(wp),dimension(ks)                       :: w5_6d
    real(wp),dimension(3*max(kx,ky,kz,kq,kr,ks)) :: w6_6d

    real(wp) :: tol
    real(wp),dimension(6) :: val,tru,err,errmax
    logical :: fail
    integer(ip) :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids
    integer(ip),dimension(6) :: iflag
    integer(ip) :: inbvx,inbvy,inbvz,inbvq,inbvr,inbvs
    integer(ip) :: iloy,iloz,iloq,ilor,ilos
    character(len=:),allocatable :: msg    !! status message associated with `flag`

    integer,dimension(*),parameter :: iflags = [   -1_ip, &
                                                   -2_ip, &
                                                    0_ip, &
                                                    1_ip, &
                                                    2_ip, &
                                                    3_ip, &
                                                    4_ip, &
                                                    5_ip, &
                                                    6_ip, &
                                                    7_ip, &
                                                    8_ip, &
                                                    9_ip, &
                                                    10_ip, &
                                                    11_ip, &
                                                    12_ip, &
                                                    13_ip, &
                                                    14_ip, &
                                                    15_ip, &
                                                    16_ip, &
                                                    17_ip, &
                                                    18_ip, &
                                                    19_ip, &
                                                    20_ip, &
                                                    21_ip, &
                                                    22_ip, &
                                                    23_ip, &
                                                    24_ip, &
                                                    25_ip, &
                                                    26_ip, &
                                                    700_ip, &
                                                    701_ip, &
                                                    702_ip, &
                                                    703_ip, &
                                                    704_ip, &
                                                    705_ip, &
                                                    706_ip, &
                                                    707_ip, &
                                                    708_ip, &
                                                    709_ip, &
                                                    710_ip, &
                                                    711_ip, &
                                                    712_ip, &
                                                    713_ip, &
                                                    714_ip, &
                                                    715_ip, &
                                                    716_ip, &
                                                    717_ip, &
                                                    800_ip, &
                                                    801_ip, &
                                                    802_ip, &
                                                    803_ip, &
                                                    804_ip, &
                                                    805_ip, &
                                                    806_ip, &
                                                    100_ip, &
                                                    101_ip, &
                                                    102_ip, &
                                                    103_ip, &
                                                    104_ip, &
                                                    201_ip, &
                                                    202_ip, &
                                                    203_ip, &
                                                    204_ip, &
                                                    301_ip, &
                                                    401_ip, &
                                                    402_ip, &
                                                    403_ip, &
                                                    404_ip, &
                                                    405_ip, &
                                                    406_ip, &
                                                    501_ip, &
                                                    502_ip, &
                                                    503_ip, &
                                                    504_ip, &
                                                    505_ip, &
                                                    506_ip, &
                                                    601_ip, &
                                                    602_ip, &
                                                    603_ip, &
                                                    604_ip, &
                                                    605_ip, &
                                                    606_ip, &
                                                    901_ip, &
                                                    902_ip, &
                                                    903_ip, &
                                                    1001_ip, &
                                                    1002_ip, &
                                                    1003_ip, &
                                                    1004_ip, &
                                                    1005_ip, &
                                                    1101_ip, &
                                                    1102_ip, &
                                                    2001_ip, &
                                                    2002_ip, &
                                                    2003_ip, &
                                                    2004_ip, &
                                                    2005_ip, &
                                                    2006_ip, &
                                                    2007_ip, &
                                                    3001_ip, &
                                                    3002_ip, &
                                                    3003_ip, &
                                                    -999999_ip ] ! an unknown code

    fail = .false.
    tol = 100 * epsilon(1.0_wp)
    idx = 0
    idy = 0
    idz = 0
    idq = 0
    idr = 0
    ids = 0

     do i=1,nx
        x(i) = real(i-1,wp)/real(nx-1,wp)
     end do
     do j=1,ny
        y(j) = real(j-1,wp)/real(ny-1,wp)
     end do
     do k=1,nz
        z(k) = real(k-1,wp)/real(nz-1,wp)
     end do
     do l=1,nq
        q(l) = real(l-1,wp)/real(nq-1,wp)
     end do
     do m=1,nr
        r(m) = real(m-1,wp)/real(nr-1,wp)
     end do
     do n=1,ns
        s(n) = real(n-1,wp)/real(ns-1,wp)
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
    call db1ink(x,nx,fcn_1d,kx,iknot,tx,bcoef_1d,iflag(1))
    call db2ink(x,nx,y,ny,fcn_2d,kx,ky,iknot,tx,ty,bcoef_2d,iflag(2))
    call db3ink(x,nx,y,ny,z,nz,fcn_3d,kx,ky,kz,iknot,tx,ty,tz,bcoef_3d,iflag(3))
    call db4ink(x,nx,y,ny,z,nz,q,nq,fcn_4d,kx,ky,kz,kq,iknot,tx,ty,tz,tq,bcoef_4d,iflag(4))
    call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,fcn_5d,kx,ky,kz,kq,kr,iknot,tx,ty,tz,tq,tr,bcoef_5d,iflag(5))
    call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,fcn_6d,kx,ky,kz,kq,kr,ks,iknot,tx,ty,tz,tq,tr,ts,bcoef_6d,iflag(6))

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
                                            tx,nx,kx,bcoef_1d,val(1),iflag(1),inbvx,&
                                            w1_1d)
                        tru(1)    = f1(x(i))
                        err(1)    = abs(tru(1)-val(1))
                        errmax(1) = max(err(1),errmax(1))
        do j=1,ny
                        call db2val(x(i),y(j),idx,idy,&
                                            tx,ty,nx,ny,kx,ky,bcoef_2d,val(2),iflag(2),&
                                            inbvx,inbvy,iloy,&
                                            w1_2d,w2_2d)
                        tru(2)    = f2(x(i),y(j))
                        err(2)    = abs(tru(2)-val(2))
                        errmax(2) = max(err(2),errmax(2))
           do k=1,nz
                        call db3val(x(i),y(j),z(k),idx,idy,idz,&
                                            tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef_3d,val(3),iflag(3),&
                                            inbvx,inbvy,inbvz,iloy,iloz,&
                                            w1_3d,w2_3d,w3_3d)
                        tru(3)    = f3(x(i),y(j),z(k))
                        err(3)    = abs(tru(3)-val(3))
                        errmax(3) = max(err(3),errmax(3))
              do l=1,nq
                        call db4val(x(i),y(j),z(k),q(l),idx,idy,idz,idq,&
                                            tx,ty,tz,tq,nx,ny,nz,nq,kx,ky,kz,kq,bcoef_4d,val(4),iflag(4),&
                                            inbvx,inbvy,inbvz,inbvq,iloy,iloz,iloq,&
                                            w1_4d,w2_4d,w3_4d,w4_4d)
                        tru(4)    = f4(x(i),y(j),z(k),q(l))
                        err(4)    = abs(tru(4)-val(4))
                        errmax(4) = max(err(4),errmax(4))
                do m=1,nr
                        call db5val(x(i),y(j),z(k),q(l),r(m),idx,idy,idz,idq,idr,&
                                            tx,ty,tz,tq,tr,nx,ny,nz,nq,nr,kx,ky,kz,kq,kr,bcoef_5d,val(5),iflag(5),&
                                            inbvx,inbvy,inbvz,inbvq,inbvr,iloy,iloz,iloq,ilor,&
                                            w1_5d,w2_5d,w3_5d,w4_5d,w5_5d)
                        tru(5)    = f5(x(i),y(j),z(k),q(l),r(m))
                        err(5)    = abs(tru(5)-val(5))
                        errmax(5) = max(err(5),errmax(5))
                    do n=1,ns
                        call db6val(x(i),y(j),z(k),q(l),r(m),s(n),idx,idy,idz,idq,idr,ids,&
                                            tx,ty,tz,tq,tr,ts,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,bcoef_6d,val(6),iflag(6),&
                                            inbvx,inbvy,inbvz,inbvq,inbvr,inbvs,iloy,iloz,iloq,ilor,ilos,&
                                            w1_6d,w2_6d,w3_6d,w4_6d,w5_6d,w6_6d)
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

    ! a test just to get all the status messages
    ! [see get_status_message]
    do i = 1, size(iflags)
        msg = get_status_message(iflags(i))
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
