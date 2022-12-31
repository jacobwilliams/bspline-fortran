!*****************************************************************************************
!>
!  [[defc]] test

    program bspline_defc_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip
    use pyplot_module

    implicit none

    integer(ip),parameter :: ndata = 100 !! number of data points
    integer(ip),parameter :: nord  = 4     !! cubic spline
    integer(ip),parameter :: nbkpt = 2*(nord-1) + (ndata/5.0_wp) !! 200 knots plus endpoint knots
    integer(ip),parameter :: n     = nbkpt - nord

    real(wp),dimension(ndata) :: xdata  !! x data
    real(wp),dimension(ndata) :: ydata  !! y data
    real(wp),dimension(ndata) :: sddata !! Y value standard deviation or uncertainty.
                                        !! A zero value for any entry of
                                        !! SDDATA(*) will weight that data point as 1.
                                        !! Otherwise the weight of that data point is
                                        !! the reciprocal of this entry.
    real(wp),dimension(nbkpt) :: bkpt   !! knots of the B-spline
    real(wp),dimension(n) :: coeff  !! computed B-spline coefficients
    real(wp),dimension(:),allocatable :: w !! workspace array
    real(wp),dimension(ndata) :: ydata_est  !! y data estimated

    real(wp) :: r !! random number
    real(wp) :: dx !! step size for knots
    integer(ip) :: i,j !! counter
    integer(ip) :: lw
    integer(ip) :: mdeout
    integer(ip) :: istat    !! pyplot-fortran status flag
    integer(ip) :: iflag    !! status flag
    integer(ip) :: inbvx
    integer(ip) :: mdein
    integer(ip) :: isize !! for `random_seed`
    integer(ip),dimension(:),allocatable :: iseed
    real(wp),dimension(3_ip*NORD) :: w0  !! work array
    type(pyplot) :: plt
    real(wp),dimension(ndata+nord) :: tx
    real(wp),dimension(ndata) :: bcoef
    real(wp),dimension(ndata*10) :: x_int, y_int

    call random_seed(size=isize)
    allocate(iseed(isize)); iseed = 42_ip
    call random_seed(put=iseed)

    lw = (nbkpt-nord+3)*(nord+1)+ &
         (nbkpt+1)*(nord+1)+ &
         2*max(ndata,nbkpt)+nbkpt+nord**2
    allocate(w(lw)); w = 0.0_wp

    ! generate some slightly randomized data:
    sddata = 0.0_wp
    do i = 1, ndata
        call random_number(r)
        r = (r - 0.5_wp)/ 5.0_wp   ! some random noise
        xdata(i) = i / (ndata/10.0_wp)
        ydata(i) = sin(xdata(i)) + r
        sddata(i) = abs(r)
    end do

    ! first compute the interpolating spline
    call db1ink(xdata,ndata,ydata,&
                kx=nord,iknot=0,tx=tx,bcoef=bcoef,iflag=iflag)
    if (iflag /= 0) then
        write(*,*) 'db1ink iflag = ', iflag
        error stop 'error calling db1ink'
    end if
    inbvx = 1
    w0 = 0.0_wp
    do i = 1, ndata*10
        x_int(i) = max(xdata(1), i / (ndata/10.0_wp) / 10.0_wp)
        call db1val(xval  = x_int(i),&
                    idx   = 0_ip,&
                    tx    = tx,&
                    nx    = ndata,&
                    kx    = NORD,&
                    bcoef = bcoef,&
                    f     = y_int(i),&
                    iflag = iflag,&
                    inbvx = inbvx,&
                    w0    = w0)
    end do
   ! write(*,*) 'x_int=',x_int

    ! now, do the least squares splines:

    !    end knots : bkpt(i),i=1,...,nord-1
    ! data interval: bkpt(nord) .... bkpt(nbkpt-nord+1)
    !   end knots  : bkpt(i),i=nbkpt-nord+2,...,nbkpt

    ! knot arrays [maybe try it with different knots...]
    bkpt(1:nord-1) = xdata(1) ! endpoints
    dx = (xdata(ndata) - xdata(1)) / (nbkpt-2*nord+1)
    j = 0
    do i = nord, nbkpt-nord+1
        bkpt(i) = xdata(1) + j*dx
        j = j + 1
    end do
    bkpt(nbkpt-nord+2:nbkpt) = xdata(ndata) + 0.1_wp ! endpoints

    mdein = 1
    call defc(Ndata, Xdata, Ydata, Sddata, Nord, Nbkpt, Bkpt, Mdein, &
              Mdeout, Coeff, Lw, w)
    if (mdeout /= 1) error stop 'ERROR'

    ! use the splines to interpolate the data:
    inbvx = 1
    w0 = 0.0_wp
    do i = 1, ndata
        call db1val(xval  = xdata(i),&
                    idx   = 0_ip,&
                    tx    = bkpt,&
                    nx    = NBKPT-NORD,&
                    kx    = NORD,&
                    bcoef = COEFF,&
                    f     = ydata_est(i),&
                    iflag = iflag,&
                    inbvx = inbvx,&
                    w0    = w0)
        if (iflag /= 0) error stop 'error calling db1val'
    end do

    ! make a plot:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                        figsize=[20,10],font_size=20,axes_labelsize=20,&
                        xtick_labelsize=20, ytick_labelsize=20,&
                        legend_fontsize=20,&
                        title='bspline_defc_test',legend=.true.)
    call plt%add_plot(xdata,ydata,&
                        label='Original points',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x_int,y_int,&
                        label='Interpolating bspline',&
                        linestyle='b-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(xdata,ydata_est,&
                        label='Least squares bspline',&
                        linestyle='r-',markersize=2,linewidth=2,istat=istat)
    call plt%savefig(pyfile='bspline_defc_test.py', figfile='bspline_defc_test.png',istat=istat)

    end program bspline_defc_test
!*****************************************************************************************