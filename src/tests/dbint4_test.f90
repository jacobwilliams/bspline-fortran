!*****************************************************************************************
!>
!  Test of [[dbint4]] routine

    program dbint4_test

    use bspline_module
    use bspline_kinds_module, only: wp
    use pyplot_module

    implicit none

    integer,parameter :: nx = 7     !! number of points in x
    integer,parameter :: kx = 4     !! order in x

    logical,parameter :: extrap = .true.
    real(wp),parameter :: rad2deg = 180.0_wp / acos(-1.0_wp)  !! deg. to radians conversion factor

    real(wp) :: x(nx)
    real(wp),dimension(:),allocatable :: xval,fval  !(nx+2)
    real(wp) :: fcn_1d(nx)
    real(wp),dimension(:),allocatable :: tx      !(nx+2+kx)
    real(wp),dimension(:),allocatable :: bcoef   !(nx+2)
    real(wp) :: tol,val,tru,err,errmax,fbcr,fbcl
    logical  :: fail
    integer  :: n,i,idx,iflag,inbvx,ibcl,ibcr,kntopt
    real(wp),dimension(5,nx+2) :: w
    type(pyplot) :: plt
    integer :: istat
    integer :: icase
    real(wp),dimension(3) :: tleft, tright
    character(len=*),dimension(7),parameter :: labels = ['not-a-knot [db1ink]          ', &
                                                         '2nd der=0, kntopt=1 [dbint4] ', &
                                                         '1st der=0, kntopt=1 [dbint4] ', &
                                                         '2nd der=0, kntopt=2 [dbint4] ', &
                                                         '1st der=0, kntopt=2 [dbint4] ', &
                                                         '2nd der=0, kntopt=3 [dbint4] ', &
                                                         '1st der=0, kntopt=3 [dbint4] ']
    character(len=*),dimension(7),parameter :: linestyles = ['r-', &
                                                             'g-', &
                                                             'b-', &
                                                             'c-', &
                                                             'm-', &
                                                             'c:', &
                                                             'm:' ]

    fail   = .false.
    tol    = 1.0e-14_wp
    idx    = 0
    x      = real([1,2,3,4,5,6,7], wp)  ! nx points
    fcn_1d = f1(x)

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x (deg)',ylabel='f(x)',&
                        title='B-Spline End Conditions',legend=.true.)
    call plt%add_plot(x*rad2deg,fcn_1d,label='Function $f(x) = \sin(x)$',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)

    if (extrap) then
        ! points to evaluate [with extrapolation]:
        xval = real([(real(i)/100.0_wp, i=50, 750)], wp)
    else
        ! points to evaluate [no extrapolation]:
        xval = real([(real(i)/100.0_wp, i=100, 700, 10)], wp)
    end if
    allocate(fval(size(xval)))

    do icase = 1, 7

        write(*,*) ''
        write(*,*) '==============================================='
        write(*,'(I3,1X,A)') icase, ' ... '//trim(labels(icase))//' ... '

        ! 1 - use db1ink
        ! 2 - use dbint4 - constrain 2nd derivative, kntopt = 1
        ! 3 - use dvint4 - constrain 1st derivative, kntopt = 1
        ! 4 - use dbint4 - constrain 2nd derivative, kntopt = 2
        ! 5 - use dvint4 - constrain 1st derivative, kntopt = 2
        ! 6 - use dbint4 - constrain 2nd derivative, kntopt = 3
        ! 7 - use dvint4 - constrain 1st derivative, kntopt = 3

        if (allocated(bcoef)) deallocate(bcoef)
        if (allocated(tx))    deallocate(tx)

        ! initialize
        ! write(*,*) '==============================================='
        ! write(*,*) 'initialize'
        ! write(*,*) '==============================================='
        if (icase==1) then
            ! use the original init routine
            allocate(bcoef(nx))
            allocate(tx(nx+kx))
            call db1ink(x,nx,fcn_1d,kx,0,tx,bcoef,iflag)
        else
            ! use the dbint4 routine and specify the endpoint derivatives
            allocate(bcoef(nx+2))
            allocate(tx(nx+2+kx))
            if (icase==2 .or. icase==4 .or. icase==6) then
                ! constrain 2nd derivative
                ibcl = 2
                ibcr = 2
            else
                ! constrain 1st derivative
                ibcl = 1
                ibcr = 1
            end if
            fbcl = 0.0_wp
            fbcr = 0.0_wp

            write(*,*) 'kntopt:  ', kntopt

            ! WARNING: kntopt seems to have no effect on result for this example.

            select case (icase)
            case(2,3)
                kntopt = 1
                call db1ink(x,nx,fcn_1d,kx,ibcl,ibcr,fbcl,fbcr,kntopt,tx,bcoef,iflag)
            case(4,5)
                kntopt = 2
                call db1ink(x,nx,fcn_1d,kx,ibcl,ibcr,fbcl,fbcr,kntopt,tx,bcoef,iflag)
            case(6,7)
                kntopt = 3
                w = 0.0_wp
                ! WARNING: the knot values seem to make no difference in the result
                tleft  = [-999.0_wp,-999.0_wp,-999.0_wp]
                tright = [999.0_wp, 999.0_wp, 999.0_wp]
                call db1ink(x,nx,fcn_1d,kx,ibcl,ibcr,fbcl,fbcr,tleft,tright,tx,bcoef,iflag)
            end select

            write(*,*) ''
            write(*,*) 'x:  ', x
            write(*,*) ''
            write(*,*) 'tx: ', tx
            write(*,*) ''

        end if
        if (iflag/=0) then
            write(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
            stop
        end if

        ! compute max error at interpolation points
        ! write(*,*) ''
        ! write(*,*) '==============================================='
        ! write(*,*) 'db1val'
        ! write(*,*) '==============================================='

        inbvx = 1    !have to set this before the first evaluate call
        val = 0.0_wp
        errmax = 0.0_wp
        err = -99999.9_wp
        fval = 0.0_wp
        do i=1,size(xval)
            if (icase==1) then
                call db1val(xval(i),idx,tx,nx,kx,bcoef,val,iflag,inbvx,extrap=extrap)
            else
                call db1val(xval(i),idx,tx,nx+2,kx,bcoef,val,iflag,inbvx,extrap=extrap)
            end if
            fval(i) = val  ! save it for plot
            tru     = f1(xval(i))
            err     = abs(tru-val)
            errmax  = max(err,errmax)
            !write(*,*) '1D: xval(i),err = ', xval(i),err
            if (iflag/=0) then
                write(*,*) 'Error evaluating 1D spline: '//get_status_message(iflag)
                write(*,*) 'x = ', xval(i)
                stop
            end if
        end do
        !write(*,*) '==============================================='

        ! check max error against tolerance
        write(*,*) '1D: max error:', errmax
        if (errmax >= tol) then
            write(*,*)  ' ** test failed ** '
        else
            write(*,*)  ' ** test passed ** '
        end if
        !write(*,*) ''
        !write(*,*) '==============================================='

        call plt%add_plot(xval*rad2deg,fval,&
                label=trim(labels(icase)),&
                linestyle=linestyles(icase),linewidth=2,istat=istat)

    end do

    call plt%savefig('dbint4_test.png',istat=istat)

    contains

        pure elemental real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp),intent(in) :: x
        f1 = sin(x)
        !f1 = 0.5_wp * (x*exp(-x) + sin(x) )
        end function f1

    end program dbint4_test
