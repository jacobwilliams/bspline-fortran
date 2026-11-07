!*****************************************************************************************
!>
!  Test [[defc]], [[dfc]], and [[dcv]] with a B-spline order greater than 20
!  (there used to be a hard-coded limit of `NORD <= 20`).
!
!  A quadratic polynomial lies in the spline space, so the least squares
!  fit should reproduce it (to within roundoff).
!
!  The data includes the right endpoint `x = BKPT(NBKPT-NORD+1)` with repeated
!  end knots, which exercises the `ILEFT <= N` interval search in [[defc]].

    program bspline_defc_high_order_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: ndata = 200 !! number of data points
    integer(ip),parameter :: nord  = 25  !! B-spline order (> old limit of 20)
    integer(ip),parameter :: nbkpt = 2*nord + 2 !! a few interior knots
    integer(ip),parameter :: n     = nbkpt - nord
    integer(ip),parameter :: nconst = 2
    real(wp),parameter :: tol = 1.0e-7_wp

    real(wp),dimension(ndata) :: xdata, ydata, sddata
    real(wp),dimension(nbkpt) :: bkpt
    real(wp),dimension(n) :: coeff
    real(wp),dimension(:),allocatable :: w
    integer(ip),dimension(:),allocatable :: iw
    real(wp),dimension(3_ip*nord) :: w0
    real(wp),dimension(nconst) :: xconst, yconst
    integer(ip),dimension(nconst) :: nderiv
    real(wp) :: f, err, var, dx
    integer(ip) :: i, lw, mdeout, iflag, inbvx, mode, l, iw1, neqcon, nincon

    do i = 1, ndata
        xdata(i) = real(i-1,wp) / real(ndata-1,wp)
        ydata(i) = func(xdata(i))
    end do
    sddata = 0.0_wp

    ! knots: nord-fold end knots, uniform interior knots
    bkpt(1:nord) = 0.0_wp
    bkpt(nbkpt-nord+1:nbkpt) = 1.0_wp
    dx = 1.0_wp / real(nbkpt-2*nord+1,wp)
    do i = nord+1, nbkpt-nord
        bkpt(i) = (i-nord)*dx
    end do

    !--------------------------------------------------------
    ! defc
    lw = (nbkpt-nord+3)*(nord+1) + (nbkpt+1)*(nord+1) + &
         2*max(ndata,nbkpt) + nbkpt + nord**2
    allocate(w(lw)); w = 0.0_wp
    call defc(ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, 1_ip, &
              mdeout, coeff, lw, w)
    if (mdeout /= 1) error stop 'error calling defc'
    err = max_error()
    write(*,*) 'defc max error: ', err
    if (.not. (err <= tol)) error stop 'defc error too large'

    !--------------------------------------------------------
    ! dfc with derivative constraints (exercises dfspvd), and covariance for dcv
    xconst = [0.0_wp, 1.0_wp]
    yconst = [dfunc(0.0_wp), dfunc(1.0_wp)]
    nderiv = 2 + 4*1  ! 1st derivative == y
    neqcon = 2
    nincon = 0
    mode   = 2        ! new problem, with covariance
    l   = nbkpt-nord+1
    iw1 = nincon+2*l
    lw  = (nbkpt-nord+3)*(nord+1)+2*max(ndata,nbkpt)+nbkpt+nord**2 + &
          (l+nconst)*l+2*(neqcon+l)+(nincon+l)+(nincon+2)*(l+6)
    deallocate(w); allocate(w(lw)); w = 0.0_wp
    allocate(iw(iw1)); iw = 0
    iw(1) = lw
    iw(2) = iw1
    call dfc(ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, &
             nconst, xconst, yconst, nderiv, mode, coeff, w, iw)
    if (mode /= 0) then
        write(*,*) 'dfc mode = ', mode
        error stop 'error calling dfc'
    end if
    err = max_error()
    write(*,*) 'dfc  max error: ', err
    if (.not. (err <= tol)) error stop 'dfc error too large'

    !--------------------------------------------------------
    ! dcv: variance should be finite and nonnegative
    do i = 1, ndata, 20
        var = dcv(xdata(i), ndata, nconst, nord, nbkpt, bkpt, w)
        if (var < 0.0_wp .or. var /= var .or. var > huge(1.0_wp)) error stop 'invalid dcv result'
    end do
    write(*,*) 'dcv ok'

    write(*,*) 'PASSED'

    contains

        pure real(wp) function func(x)
        real(wp),intent(in) :: x
        func = 1.0_wp + 2.0_wp*x - 3.0_wp*x**2
        end function func

        pure real(wp) function dfunc(x)
        real(wp),intent(in) :: x
        dfunc = 2.0_wp - 6.0_wp*x
        end function dfunc

        real(wp) function max_error()
        max_error = 0.0_wp
        inbvx = 1
        w0 = 0.0_wp
        do i = 1, ndata
            call db1val(xdata(i), 0_ip, bkpt, n, nord, coeff, f, iflag, inbvx, w0)
            if (iflag /= 0) error stop 'error calling db1val'
            if (f /= f) then
                max_error = huge(1.0_wp) ! NaN
                return
            end if
            max_error = max(max_error, abs(f - ydata(i)))
        end do
        end function max_error

    end program bspline_defc_high_order_test
!*****************************************************************************************
