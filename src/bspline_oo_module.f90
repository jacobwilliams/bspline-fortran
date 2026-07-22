!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 12/6/2015
!
!  Object-oriented style wrappers to [[bspline_sub_module]].
!  This module provides classes ([[bspline_1d(type)]], [[bspline_2d(type)]],
!  [[bspline_3d(type)]], [[bspline_4d(type)]], [[bspline_5d(type)]], and [[bspline_6d(type)]])
!  which can be used instead of the main subroutine interface.

    module bspline_oo_module

    use bspline_kinds_module, only: wp, ip
    use,intrinsic :: iso_fortran_env, only: error_unit
    use bspline_sub_module

    implicit none

    private

    integer(ip),parameter :: int_size     = storage_size(1_ip,kind=ip)   !! size of a default integer [bits]
    integer(ip),parameter :: logical_size = storage_size(.true.,kind=ip) !! size of a default logical [bits]
    integer(ip),parameter :: real_size    = storage_size(1.0_wp,kind=ip) !! size of a `real(wp)` [bits]

    type,public,abstract :: bspline_class
        !! Base class for the b-spline types
        private
        integer(ip) :: inbvx = 1_ip  !! internal variable used by [[dbvalu]] for efficient processing
        integer(ip) :: iflag = 1_ip  !! saved `iflag` from the list routine call.
        logical :: initialized = .false. !! true if the class is initialized and ready to use
        logical :: extrap = .false. !! if true, then extrapolation is allowed during evaluation
    contains
        private
        procedure,non_overridable :: destroy_base  !! destructor for the abstract type
        procedure,non_overridable :: set_extrap_flag !! internal routine to set the `extrap` flag
        procedure(destroy_func),deferred,public :: destroy  !! destructor
        procedure(size_func),deferred,public :: size_of !! size of the structure in bits
        procedure,public,non_overridable :: status_ok  !! returns true if the last `iflag` status code was `=0`.
        procedure,public,non_overridable :: status_message => get_bspline_status_message  !! retrieve the last
                                                                                          !! status message
        procedure,public,non_overridable :: clear_flag => clear_bspline_flag  !! to reset the `iflag` saved in the class.
    end type bspline_class

    abstract interface

        pure subroutine destroy_func(me)
        !! interface for bspline destructor routines
        import :: bspline_class
        implicit none
        class(bspline_class),intent(inout) :: me
        end subroutine destroy_func

        pure function size_func(me) result(s)
        !! interface for size routines
        import :: bspline_class,ip
        implicit none
        class(bspline_class),intent(in) :: me
        integer(ip) :: s !! size of the structure in bits
        end function size_func

    end interface

    type,extends(bspline_class),public :: bspline_1d
        !! Class for 1d b-spline interpolation.
        !!
        !!@note The 1D class also contains two methods
        !!      for computing definite integrals.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        real(wp),dimension(:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: work_val_1   !! [[db1val] work array of dimension `3*kx`
        contains
        private
        generic,public :: initialize => initialize_1d_auto_knots,initialize_1d_specify_knots
        procedure :: initialize_1d_auto_knots
        procedure :: initialize_1d_specify_knots
        procedure,public :: evaluate => evaluate_1d
        procedure,public :: destroy => destroy_1d
        procedure,public :: size_of => size_1d
        procedure,public :: integral => integral_1d
        procedure,public :: fintegral => fintegral_1d
        final :: finalize_1d
    end type bspline_1d

    type,extends(bspline_class),public :: bspline_2d
        !! Class for 2d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        real(wp),dimension(:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:),allocatable :: work_val_1  !! [[db2val] work array of dimension `ky`
        real(wp),dimension(:),allocatable :: work_val_2  !! [[db2val] work array of dimension `3_ip*max(kx,ky)`
        contains
        private
        generic,public :: initialize => initialize_2d_auto_knots,initialize_2d_specify_knots
        procedure :: initialize_2d_auto_knots
        procedure :: initialize_2d_specify_knots
        procedure,public :: evaluate => evaluate_2d
        procedure,public :: destroy => destroy_2d
        procedure,public :: size_of => size_2d
        final :: finalize_2d
    end type bspline_2d

    type,extends(bspline_class),public :: bspline_3d
        !! Class for 3d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: nz = 0_ip  !! Number of \(z\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        integer(ip) :: kz = 0_ip  !! The order of spline pieces in \(z\)
        real(wp),dimension(:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvz = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloz = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:,:),allocatable :: work_val_1  !! [[db3val] work array of dimension `ky,kz`
        real(wp),dimension(:),allocatable   :: work_val_2  !! [[db3val] work array of dimension `kz`
        real(wp),dimension(:),allocatable   :: work_val_3  !! [[db3val] work array of dimension `3_ip*max(kx,ky,kz)`
        contains
        private
        generic,public :: initialize => initialize_3d_auto_knots,initialize_3d_specify_knots
        procedure :: initialize_3d_auto_knots
        procedure :: initialize_3d_specify_knots
        procedure,public :: evaluate => evaluate_3d
        procedure,public :: destroy => destroy_3d
        procedure,public :: size_of => size_3d
        final :: finalize_3d
    end type bspline_3d

    type,extends(bspline_class),public :: bspline_4d
        !! Class for 4d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: nz = 0_ip  !! Number of \(z\) abcissae
        integer(ip) :: nq = 0_ip  !! Number of \(q\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        integer(ip) :: kz = 0_ip  !! The order of spline pieces in \(z\)
        integer(ip) :: kq = 0_ip  !! The order of spline pieces in \(q\)
        real(wp),dimension(:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvz = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvq = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloz  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloq  = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:,:,:),allocatable :: work_val_1  !! [[db4val]] work array of dimension `ky,kz,kq`
        real(wp),dimension(:,:),allocatable   :: work_val_2  !! [[db4val]] work array of dimension `kz,kq`
        real(wp),dimension(:),allocatable     :: work_val_3  !! [[db4val]] work array of dimension `kq`
        real(wp),dimension(:),allocatable     :: work_val_4  !! [[db4val]] work array of dimension `3_ip*max(kx,ky,kz,kq)`
        contains
        private
        generic,public :: initialize => initialize_4d_auto_knots,initialize_4d_specify_knots
        procedure :: initialize_4d_auto_knots
        procedure :: initialize_4d_specify_knots
        procedure,public :: evaluate => evaluate_4d
        procedure,public :: destroy => destroy_4d
        procedure,public :: size_of => size_4d
        final :: finalize_4d
    end type bspline_4d

    type,extends(bspline_class),public :: bspline_5d
        !! Class for 5d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: nz = 0_ip  !! Number of \(z\) abcissae
        integer(ip) :: nq = 0_ip  !! Number of \(q\) abcissae
        integer(ip) :: nr = 0_ip  !! Number of \(r\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        integer(ip) :: kz = 0_ip  !! The order of spline pieces in \(z\)
        integer(ip) :: kq = 0_ip  !! The order of spline pieces in \(q\)
        integer(ip) :: kr = 0_ip  !! The order of spline pieces in \(r\)
        real(wp),dimension(:,:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tr  !! The knots in the \(r\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvz = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvq = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvr = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloz  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloq  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: ilor  = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:,:,:,:),allocatable :: work_val_1  !! [[db5val]] work array of dimension `ky,kz,kq,kr`
        real(wp),dimension(:,:,:),allocatable   :: work_val_2  !! [[db5val]] work array of dimension `kz,kq,kr`
        real(wp),dimension(:,:),allocatable     :: work_val_3  !! [[db5val]] work array of dimension `kq,kr`
        real(wp),dimension(:),allocatable       :: work_val_4  !! [[db5val]] work array of dimension `kr`
        real(wp),dimension(:),allocatable       :: work_val_5  !! [[db5val]] work array of dimension `3_ip*max(kx,ky,kz,kq,kr)`
        contains
        private
        generic,public :: initialize => initialize_5d_auto_knots,initialize_5d_specify_knots
        procedure :: initialize_5d_auto_knots
        procedure :: initialize_5d_specify_knots
        procedure,public :: evaluate => evaluate_5d
        procedure,public :: destroy => destroy_5d
        procedure,public :: size_of => size_5d
        final :: finalize_5d
    end type bspline_5d

    type,extends(bspline_class),public :: bspline_6d
        !! Class for 6d b-spline interpolation.
        private
        integer(ip) :: nx = 0_ip  !! Number of \(x\) abcissae
        integer(ip) :: ny = 0_ip  !! Number of \(y\) abcissae
        integer(ip) :: nz = 0_ip  !! Number of \(z\) abcissae
        integer(ip) :: nq = 0_ip  !! Number of \(q\) abcissae
        integer(ip) :: nr = 0_ip  !! Number of \(r\) abcissae
        integer(ip) :: ns = 0_ip  !! Number of \(s\) abcissae
        integer(ip) :: kx = 0_ip  !! The order of spline pieces in \(x\)
        integer(ip) :: ky = 0_ip  !! The order of spline pieces in \(y\)
        integer(ip) :: kz = 0_ip  !! The order of spline pieces in \(z\)
        integer(ip) :: kq = 0_ip  !! The order of spline pieces in \(q\)
        integer(ip) :: kr = 0_ip  !! The order of spline pieces in \(r\)
        integer(ip) :: ks = 0_ip  !! The order of spline pieces in \(s\)
        real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tr  !! The knots in the \(r\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ts  !! The knots in the \(s\) direction for the spline interpolant
        integer(ip) :: inbvy = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvz = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvq = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvr = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: inbvs = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloy  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloz  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: iloq  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: ilor  = 1_ip  !! internal variable used for efficient processing
        integer(ip) :: ilos  = 1_ip  !! internal variable used for efficient processing
        real(wp),dimension(:,:,:,:,:),allocatable :: work_val_1  !! [[db6val]] work array of dimension `ky,kz,kq,kr,ks`
        real(wp),dimension(:,:,:,:),allocatable   :: work_val_2  !! [[db6val]] work array of dimension `kz,kq,kr,ks`
        real(wp),dimension(:,:,:),allocatable     :: work_val_3  !! [[db6val]] work array of dimension `kq,kr,ks`
        real(wp),dimension(:,:),allocatable       :: work_val_4  !! [[db6val]] work array of dimension `kr,ks`
        real(wp),dimension(:),allocatable         :: work_val_5  !! [[db6val]] work array of dimension `ks`
        real(wp),dimension(:),allocatable         :: work_val_6  !! [[db6val]] work array of dimension `3_ip*max(kx,ky,kz,kq,kr,ks)`
        contains
        private
        generic,public :: initialize => initialize_6d_auto_knots,initialize_6d_specify_knots
        procedure :: initialize_6d_auto_knots
        procedure :: initialize_6d_specify_knots
        procedure,public :: evaluate => evaluate_6d
        procedure,public :: destroy => destroy_6d
        procedure,public :: size_of => size_6d
        final :: finalize_6d
    end type bspline_6d

    interface bspline_1d
        !! Constructor for [[bspline_1d(type)]]
        procedure :: bspline_1d_constructor_empty,&
                     bspline_1d_constructor_auto_knots,&
                     bspline_1d_constructor_specify_knots
    end interface
    interface bspline_2d
        !! Constructor for [[bspline_2d(type)]]
        procedure :: bspline_2d_constructor_empty,&
                     bspline_2d_constructor_auto_knots,&
                     bspline_2d_constructor_specify_knots
    end interface
    interface bspline_3d
        !! Constructor for [[bspline_3d(type)]]
        procedure :: bspline_3d_constructor_empty,&
                     bspline_3d_constructor_auto_knots,&
                     bspline_3d_constructor_specify_knots
    end interface
    interface bspline_4d
        !! Constructor for [[bspline_4d(type)]]
        procedure :: bspline_4d_constructor_empty,&
                     bspline_4d_constructor_auto_knots,&
                     bspline_4d_constructor_specify_knots
    end interface
    interface bspline_5d
        !! Constructor for [[bspline_5d(type)]]
        procedure :: bspline_5d_constructor_empty,&
                     bspline_5d_constructor_auto_knots,&
                     bspline_5d_constructor_specify_knots
    end interface
    interface bspline_6d
        !! Constructor for [[bspline_6d(type)]]
        procedure :: bspline_6d_constructor_empty,&
                     bspline_6d_constructor_auto_knots,&
                     bspline_6d_constructor_specify_knots
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  This routines returns true if the `iflag` code from the last
!  routine called was `=0`. Maybe of the routines have output `iflag`
!  variables, so they can be checked explicitly, or this routine
!  can be used.
!
!  If the class is initialized using a function constructor, then
!  this is the only way to know if it was properly initialized,
!  since those are pure functions with not output `iflag` arguments.
!
!  If `status_ok=.false.`, then the error message can be
!  obtained from the [[get_bspline_status_message]] routine.
!
!  Note: after an error condition, the [[clear_bspline_flag]] routine
!  can be called to reset the `iflag` to 0.

    elemental function status_ok(me) result(ok)

    implicit none

    class(bspline_class),intent(in) :: me
    logical                         :: ok

    ok = ( me%iflag == 0_ip )

    end function status_ok
!*****************************************************************************************

!*****************************************************************************************
!>
!  This sets the `iflag` variable in the class to `0`
!  (which indicates that everything is OK). It can be used
!  after an error is encountered.

    elemental subroutine clear_bspline_flag(me)

    implicit none

    class(bspline_class),intent(inout) :: me

    me%iflag = 0_ip

    end subroutine clear_bspline_flag
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get the status message from a [[bspline_class]] routine call.
!
!  If `iflag` is not included, then the one in the class is used (which
!  corresponds to the last routine called.)
!  Otherwise, it will convert the
!  input `iflag` argument into the appropriate message.
!
!  This is a wrapper for [[get_status_message]].

    pure function get_bspline_status_message(me,iflag) result(msg)

    implicit none

    class(bspline_class),intent(in) :: me
    character(len=:),allocatable    :: msg    !! status message associated with the flag
    integer(ip),intent(in),optional :: iflag  !! the corresponding status code

    if (present(iflag)) then
        msg = get_status_message(iflag)
    else
        msg = get_status_message(me%iflag)
    end if

    end function get_bspline_status_message
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_1d]] structure in bits.

    pure function size_1d(me) result(s)

    implicit none

    class(bspline_1d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 2_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,kind=ip)

    end function size_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_2d]] structure in bits.

    pure function size_2d(me) result(s)

    implicit none

    class(bspline_2d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 6_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,kind=ip)

    end function size_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_3d]] structure in bits.

    pure function size_3d(me) result(s)

    implicit none

    class(bspline_3d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 10_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)*&
                                                    size(me%bcoef,3_ip,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%tz))         s = s + real_size*size(me%tz,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,1_ip,kind=ip)*&
                                                    size(me%work_val_1,2_ip,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,kind=ip)
    if (allocated(me%work_val_3)) s = s + real_size*size(me%work_val_3,kind=ip)

    end function size_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_4d]] structure in bits.

    pure function size_4d(me) result(s)

    implicit none

    class(bspline_4d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 14_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)*&
                                                    size(me%bcoef,3_ip,kind=ip)*&
                                                    size(me%bcoef,4_ip,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%tz))         s = s + real_size*size(me%tz,kind=ip)
    if (allocated(me%tq))         s = s + real_size*size(me%tq,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,1_ip,kind=ip)*&
                                                    size(me%work_val_1,2_ip,kind=ip)*&
                                                    size(me%work_val_1,3_ip,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,1_ip,kind=ip)*&
                                                    size(me%work_val_2,2_ip,kind=ip)
    if (allocated(me%work_val_3)) s = s + real_size*size(me%work_val_3,kind=ip)
    if (allocated(me%work_val_4)) s = s + real_size*size(me%work_val_4,kind=ip)

    end function size_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_5d]] structure in bits.

    pure function size_5d(me) result(s)

    implicit none

    class(bspline_5d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 18_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)*&
                                                    size(me%bcoef,3_ip,kind=ip)*&
                                                    size(me%bcoef,4_ip,kind=ip)*&
                                                    size(me%bcoef,5_ip,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%tz))         s = s + real_size*size(me%tz,kind=ip)
    if (allocated(me%tq))         s = s + real_size*size(me%tq,kind=ip)
    if (allocated(me%tr))         s = s + real_size*size(me%tr,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,1_ip,kind=ip)*&
                                                    size(me%work_val_1,2_ip,kind=ip)*&
                                                    size(me%work_val_1,3_ip,kind=ip)*&
                                                    size(me%work_val_1,4_ip,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,1_ip,kind=ip)*&
                                                    size(me%work_val_2,2_ip,kind=ip)*&
                                                    size(me%work_val_2,3_ip,kind=ip)
    if (allocated(me%work_val_3)) s = s + real_size*size(me%work_val_3,1_ip,kind=ip)*&
                                                    size(me%work_val_3,2_ip,kind=ip)
    if (allocated(me%work_val_4)) s = s + real_size*size(me%work_val_4,kind=ip)
    if (allocated(me%work_val_5)) s = s + real_size*size(me%work_val_5,kind=ip)

    end function size_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_6d]] structure in bits.

    pure function size_6d(me) result(s)

    implicit none

    class(bspline_6d),intent(in) :: me
    integer(ip) :: s !! size of the structure in bits

    s = 2_ip*int_size + logical_size + 22_ip*int_size

    if (allocated(me%bcoef))      s = s + real_size*size(me%bcoef,1_ip,kind=ip)*&
                                                    size(me%bcoef,2_ip,kind=ip)*&
                                                    size(me%bcoef,3_ip,kind=ip)*&
                                                    size(me%bcoef,4_ip,kind=ip)*&
                                                    size(me%bcoef,5_ip,kind=ip)*&
                                                    size(me%bcoef,6,kind=ip)
    if (allocated(me%tx))         s = s + real_size*size(me%tx,kind=ip)
    if (allocated(me%ty))         s = s + real_size*size(me%ty,kind=ip)
    if (allocated(me%tz))         s = s + real_size*size(me%tz,kind=ip)
    if (allocated(me%tq))         s = s + real_size*size(me%tq,kind=ip)
    if (allocated(me%tr))         s = s + real_size*size(me%tr,kind=ip)
    if (allocated(me%ts))         s = s + real_size*size(me%ts,kind=ip)
    if (allocated(me%work_val_1)) s = s + real_size*size(me%work_val_1,1_ip,kind=ip)*&
                                                    size(me%work_val_1,2_ip,kind=ip)*&
                                                    size(me%work_val_1,3_ip,kind=ip)*&
                                                    size(me%work_val_1,4_ip,kind=ip)*&
                                                    size(me%work_val_1,5_ip,kind=ip)
    if (allocated(me%work_val_2)) s = s + real_size*size(me%work_val_2,1_ip,kind=ip)*&
                                                    size(me%work_val_2,2_ip,kind=ip)*&
                                                    size(me%work_val_2,3_ip,kind=ip)*&
                                                    size(me%work_val_2,4_ip,kind=ip)
    if (allocated(me%work_val_3)) s = s + real_size*size(me%work_val_3,1_ip,kind=ip)*&
                                                    size(me%work_val_3,2_ip,kind=ip)*&
                                                    size(me%work_val_3,3_ip,kind=ip)
    if (allocated(me%work_val_4)) s = s + real_size*size(me%work_val_4,1_ip,kind=ip)*&
                                                    size(me%work_val_4,2_ip,kind=ip)
    if (allocated(me%work_val_5)) s = s + real_size*size(me%work_val_5,kind=ip)
    if (allocated(me%work_val_6)) s = s + real_size*size(me%work_val_6,kind=ip)

    end function size_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for contents of the base [[bspline_class]] class.
!  (this routine is called by the extended classes).

    pure subroutine destroy_base(me)

    implicit none

    class(bspline_class),intent(inout) :: me

    me%inbvx = 1_ip
    me%iflag = 1_ip
    me%initialized = .false.
    me%extrap = .false.

    end subroutine destroy_base
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_1d]] class.

    pure subroutine destroy_1d(me)

    implicit none

    class(bspline_1d),intent(inout) :: me

    call me%destroy_base()

    me%nx = 0_ip
    me%kx = 0_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)

    end subroutine destroy_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_2d]] class.

    pure subroutine destroy_2d(me)

    implicit none

    class(bspline_2d),intent(inout) :: me

    call me%destroy_base()

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%inbvy = 1_ip
    me%iloy  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_3d]] class.

    pure subroutine destroy_3d(me)

    implicit none

    class(bspline_3d),intent(inout) :: me

    call me%destroy_base()

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%nz    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%kz    = 0_ip
    me%inbvy = 1_ip
    me%inbvz = 1_ip
    me%iloy  = 1_ip
    me%iloz  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%tz))         deallocate(me%tz)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)
    if (allocated(me%work_val_3)) deallocate(me%work_val_3)

    end subroutine destroy_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_4d]] class.

    pure subroutine destroy_4d(me)

    implicit none

    class(bspline_4d),intent(inout) :: me

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%nz    = 0_ip
    me%nq    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%kz    = 0_ip
    me%kq    = 0_ip
    me%inbvy = 1_ip
    me%inbvz = 1_ip
    me%inbvq = 1_ip
    me%iloy  = 1_ip
    me%iloz  = 1_ip
    me%iloq  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%tz))         deallocate(me%tz)
    if (allocated(me%tq))         deallocate(me%tq)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)
    if (allocated(me%work_val_3)) deallocate(me%work_val_3)
    if (allocated(me%work_val_4)) deallocate(me%work_val_4)

    end subroutine destroy_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_5d]] class.

    pure subroutine destroy_5d(me)

    implicit none

    class(bspline_5d),intent(inout) :: me

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%nz    = 0_ip
    me%nq    = 0_ip
    me%nr    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%kz    = 0_ip
    me%kq    = 0_ip
    me%kr    = 0_ip
    me%inbvy = 1_ip
    me%inbvz = 1_ip
    me%inbvq = 1_ip
    me%inbvr = 1_ip
    me%iloy  = 1_ip
    me%iloz  = 1_ip
    me%iloq  = 1_ip
    me%ilor  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%tz))         deallocate(me%tz)
    if (allocated(me%tq))         deallocate(me%tq)
    if (allocated(me%tr))         deallocate(me%tr)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)
    if (allocated(me%work_val_3)) deallocate(me%work_val_3)
    if (allocated(me%work_val_4)) deallocate(me%work_val_4)
    if (allocated(me%work_val_5)) deallocate(me%work_val_5)

    end subroutine destroy_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_6d]] class.

    pure subroutine destroy_6d(me)

    implicit none

    class(bspline_6d),intent(inout) :: me

    me%nx    = 0_ip
    me%ny    = 0_ip
    me%nz    = 0_ip
    me%nq    = 0_ip
    me%nr    = 0_ip
    me%ns    = 0_ip
    me%kx    = 0_ip
    me%ky    = 0_ip
    me%kz    = 0_ip
    me%kq    = 0_ip
    me%kr    = 0_ip
    me%ks    = 0_ip
    me%inbvy = 1_ip
    me%inbvz = 1_ip
    me%inbvq = 1_ip
    me%inbvr = 1_ip
    me%inbvs = 1_ip
    me%iloy  = 1_ip
    me%iloz  = 1_ip
    me%iloq  = 1_ip
    me%ilor  = 1_ip
    me%ilos  = 1_ip
    if (allocated(me%bcoef))      deallocate(me%bcoef)
    if (allocated(me%tx))         deallocate(me%tx)
    if (allocated(me%ty))         deallocate(me%ty)
    if (allocated(me%tz))         deallocate(me%tz)
    if (allocated(me%tq))         deallocate(me%tq)
    if (allocated(me%tr))         deallocate(me%tr)
    if (allocated(me%ts))         deallocate(me%ts)
    if (allocated(me%work_val_1)) deallocate(me%work_val_1)
    if (allocated(me%work_val_2)) deallocate(me%work_val_2)
    if (allocated(me%work_val_3)) deallocate(me%work_val_3)
    if (allocated(me%work_val_4)) deallocate(me%work_val_4)
    if (allocated(me%work_val_5)) deallocate(me%work_val_5)
    if (allocated(me%work_val_6)) deallocate(me%work_val_6)

    end subroutine destroy_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finalizer for [[bspline_1d]] class. Just a wrapper for [[destroy_1d]].
    pure elemental subroutine finalize_1d(me)
        type(bspline_1d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_1d
!*****************************************************************************************
!*****************************************************************************************
!>
!  Finalizer for [[bspline_2d]] class. Just a wrapper for [[destroy_2d]].
    pure elemental subroutine finalize_2d(me)
        type(bspline_2d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_2d
!*****************************************************************************************
!*****************************************************************************************
!>
!  Finalizer for [[bspline_3d]] class. Just a wrapper for [[destroy_3d]].
    pure elemental subroutine finalize_3d(me)
        type(bspline_3d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_3d
!*****************************************************************************************
!*****************************************************************************************
!>
!  Finalizer for [[bspline_4d]] class. Just a wrapper for [[destroy_4d]].
    pure elemental subroutine finalize_4d(me)
        type(bspline_4d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_4d
!*****************************************************************************************
!*****************************************************************************************
!>
!  Finalizer for [[bspline_5d]] class. Just a wrapper for [[destroy_5d]].
    pure elemental subroutine finalize_5d(me)
        type(bspline_5d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_5d
!*****************************************************************************************
!*****************************************************************************************
!>
!  Finalizer for [[bspline_6d]] class. Just a wrapper for [[destroy_6d]].
    pure elemental subroutine finalize_6d(me)
        type(bspline_6d),intent(inout) :: me; call me%destroy()
    end subroutine finalize_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sets the `extrap` flag in the class.

    pure subroutine set_extrap_flag(me,extrap)

    implicit none

    class(bspline_class),intent(inout) :: me
    logical,intent(in),optional :: extrap  !! if not present, then False is used

    if (present(extrap)) then
        me%extrap = extrap
    else
        me%extrap = .false.
    end if

    end subroutine set_extrap_flag
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_1d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    pure elemental function bspline_1d_constructor_empty() result(me)

    implicit none

    type(bspline_1d) :: me

    end function bspline_1d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_1d]] type (auto knots).
!  This is a wrapper for [[initialize_1d_auto_knots]].

    pure function bspline_1d_constructor_auto_knots(x,fcn,kx,extrap) result(me)

    implicit none

    type(bspline_1d)                 :: me
    real(wp),dimension(:),intent(in) :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn    !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                               !! contain the function value at the point `x(i)`
    integer(ip),intent(in)           :: kx     !! The order of spline pieces in \(x\)
                                               !! ( \( 2 \le k_x < n_x \) )
                                               !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap !! if true, then extrapolation is allowed
                                               !! (default is false)

    call initialize_1d_auto_knots(me,x,fcn,kx,me%iflag,extrap)

    end function bspline_1d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_1d]] type (user-specified knots).
!  This is a wrapper for [[initialize_1d_specify_knots]].

    pure function bspline_1d_constructor_specify_knots(x,fcn,kx,tx,extrap) result(me)

    implicit none

    type(bspline_1d)                 :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer(ip),intent(in)           :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in) :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                              !! for the spline interpolant.
                                              !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap !! if true, then extrapolation is allowed
                                               !! (default is false)

    call initialize_1d_specify_knots(me,x,fcn,kx,tx,me%iflag,extrap)

    end function bspline_1d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_1d]] type (with automatically-computed knots).
!  This is a wrapper for [[db1ink]].

    pure subroutine initialize_1d_auto_knots(me,x,fcn,kx,iflag,extrap)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer(ip),intent(in)           :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    integer(ip),intent(out)          :: iflag !! status flag (see [[db1ink]])
    logical,intent(in),optional      :: extrap !! if true, then extrapolation is allowed
                                               !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx

    call me%destroy()

    nx = size(x,kind=ip)

    me%nx = nx
    me%kx = kx

    allocate(me%tx(nx+kx))
    allocate(me%bcoef(nx))
    allocate(me%work_val_1(3_ip*kx))

    iknot = 0_ip         !knot sequence chosen by db1ink

    call db1ink(x,nx,fcn,kx,iknot,me%tx,me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_1d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_1d]] type (with user-specified knots).
!  This is a wrapper for [[db1ink]].

    pure subroutine initialize_1d_specify_knots(me,x,fcn,kx,tx,iflag,extrap)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer(ip),intent(in)           :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in) :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                              !! for the spline interpolant.
                                              !! Must be non-decreasing.
    integer(ip),intent(out)          :: iflag !! status flag (see [[db1ink]])
    logical,intent(in),optional      :: extrap !! if true, then extrapolation is allowed
                                               !! (default is false)

    integer(ip) :: nx

    call me%destroy()

    nx = size(x,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%kx = kx

        allocate(me%tx(nx+kx))
        allocate(me%bcoef(nx))
        allocate(me%work_val_1(3_ip*kx))

        me%tx = tx

        call db1ink(x,nx,fcn,kx,1_ip,me%tx,me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_1d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_1d]] interpolate.  This is a wrapper for [[db1val]].

    pure subroutine evaluate_1d(me,xval,idx,f,iflag)

    implicit none

    class(bspline_1d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    integer(ip),intent(in)          :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db1val]])

    if (me%initialized) then
        call db1val(xval,idx,me%tx,me%nx,me%kx,me%bcoef,f,iflag,&
                    me%inbvx,me%work_val_1,extrap=me%extrap)
    else
        iflag = 1_ip
    end if
    me%iflag = iflag

    end subroutine evaluate_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_1d]] definite integral.  This is a wrapper for [[db1sqad]].

    pure subroutine integral_1d(me,x1,x2,f,iflag)

    implicit none

    class(bspline_1d),intent(inout) :: me
    real(wp),intent(in)             :: x1    !! left point of interval
    real(wp),intent(in)             :: x2    !! right point of interval
    real(wp),intent(out)            :: f     !! integral of the b-spline over \( [x_1, x_2] \)
    integer(ip),intent(out)         :: iflag !! status flag (see [[db1sqad]])

    if (me%initialized) then
        call db1sqad(me%tx,me%bcoef,me%nx,me%kx,x1,x2,f,iflag,me%work_val_1)
    else
        iflag = 1_ip
    end if
    me%iflag = iflag

    end subroutine integral_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_1d]] definite integral.  This is a wrapper for [[db1fqad]].

    subroutine fintegral_1d(me,fun,idx,x1,x2,tol,f,iflag)

    implicit none

    class(bspline_1d),intent(inout) :: me
    procedure(b1fqad_func)          :: fun   !! external function of one argument for the
                                             !! integrand `bf(x)=fun(x)*dbvalu(tx,bcoef,nx,kx,idx,x,inbv)`
    integer(ip),intent(in)          :: idx   !! order of the spline derivative, `0 <= idx <= k-1`
                                             !! `idx=0` gives the spline function
    real(wp),intent(in)             :: x1    !! left point of interval
    real(wp),intent(in)             :: x2    !! right point of interval
    real(wp),intent(in)             :: tol   !! desired accuracy for the quadrature
    real(wp),intent(out)            :: f     !! integral of `bf(x)` over \( [x_1, x_2] \)
    integer(ip),intent(out)         :: iflag !! status flag (see [[db1sqad]])

    if (me%initialized) then
        call db1fqad(fun,me%tx,me%bcoef,me%nx,me%kx,idx,x1,x2,tol,f,iflag,me%work_val_1)
    else
        iflag = 1_ip
    end if
    me%iflag = iflag

    end subroutine fintegral_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_2d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_2d_constructor_empty() result(me)

    implicit none

    type(bspline_2d) :: me

    end function bspline_2d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_2d]] type (auto knots).
!  This is a wrapper for [[initialize_2d_auto_knots]].

    pure function bspline_2d_constructor_auto_knots(x,y,fcn,kx,ky,extrap) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_2d_auto_knots(me,x,y,fcn,kx,ky,me%iflag,extrap)

    end function bspline_2d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_2d]] type (user-specified knots).
!  This is a wrapper for [[initialize_2d_specify_knots]].

    pure function bspline_2d_constructor_specify_knots(x,y,fcn,kx,ky,tx,ty,extrap) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,me%iflag,extrap)

    end function bspline_2d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with automatically-computed knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_auto_knots(me,x,y,fcn,kx,ky,iflag,extrap)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(out)            :: iflag !! status flag (see [[db2ink]])
    logical,intent(in),optional        :: extrap !! if true, then extrapolation is allowed
                                                 !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)

    me%nx = nx
    me%ny = ny

    me%kx = kx
    me%ky = ky

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%bcoef(nx,ny))
    allocate(me%work_val_1(ky))
    allocate(me%work_val_2(3_ip*max(kx,ky)))

    iknot = 0_ip         !knot sequence chosen by db2ink

    call db2ink(x,nx,y,ny,fcn,kx,ky,iknot,me%tx,me%ty,me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_2d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with user-specified knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,iflag,extrap)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer(ip),intent(in)             :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)             :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    integer(ip),intent(out)            :: iflag !! status flag (see [[db2ink]])
    logical,intent(in),optional      :: extrap  !! if true, then extrapolation is allowed
                                                !! (default is false)

    integer(ip) :: nx,ny

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny

        me%kx = kx
        me%ky = ky

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%bcoef(nx,ny))
        allocate(me%work_val_1(ky))
        allocate(me%work_val_2(3_ip*max(kx,ky)))

        me%tx = tx
        me%ty = ty

        call db2ink(x,nx,y,ny,fcn,kx,ky,1_ip,me%tx,me%ty,me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_2d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_2d]] interpolate.  This is a wrapper for [[db2val]].

    pure subroutine evaluate_2d(me,xval,yval,idx,idy,f,iflag)

    implicit none

    class(bspline_2d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    integer(ip),intent(in)           :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)           :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db2val]])

    if (me%initialized) then
        call db2val(xval,yval,&
                    idx,idy,&
                    me%tx,me%ty,&
                    me%nx,me%ny,&
                    me%kx,me%ky,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%iloy,&
                    me%work_val_1,me%work_val_2,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_3d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_3d_constructor_empty() result(me)

    implicit none

    type(bspline_3d) :: me

    end function bspline_3d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_3d]] type (auto knots).
!  This is a wrapper for [[initialize_3d_auto_knots]].

    pure function bspline_3d_constructor_auto_knots(x,y,z,fcn,kx,ky,kz,extrap) result(me)

    implicit none

    type(bspline_3d)                     :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer(ip),intent(in)               :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)
    logical,intent(in),optional       :: extrap !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_3d_auto_knots(me,x,y,z,fcn,kx,ky,kz,me%iflag,extrap)

    end function bspline_3d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_3d]] type (user-specified knots).
!  This is a wrapper for [[initialize_3d_specify_knots]].

    pure function bspline_3d_constructor_specify_knots(x,y,z,fcn,kx,ky,kz,tx,ty,tz,extrap) result(me)

    implicit none

    type(bspline_3d)                     :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer(ip),intent(in)               :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)     :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)     :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)     :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    logical,intent(in),optional :: extrap       !! if true, then extrapolation is allowed
                                                !! (default is false)

    call initialize_3d_specify_knots(me,x,y,z,fcn,kx,ky,kz,tx,ty,tz,me%iflag,extrap)

    end function bspline_3d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_3d]] type (with automatically-computed knots).
!  This is a wrapper for [[db3ink]].

    pure subroutine initialize_3d_auto_knots(me,x,y,z,fcn,kx,ky,kz,iflag,extrap)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer(ip),intent(in)               :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(out) :: iflag            !! status flag (see [[db3ink]])
    logical,intent(in),optional :: extrap       !! if true, then extrapolation is allowed
                                                !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny,nz

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)

    me%nx = nx
    me%ny = ny
    me%nz = nz

    me%kx = kx
    me%ky = ky
    me%kz = kz

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%tz(nz+kz))
    allocate(me%bcoef(nx,ny,nz))
    allocate(me%work_val_1(ky,kz))
    allocate(me%work_val_2(kz))
    allocate(me%work_val_3(3_ip*max(kx,ky,kz)))

    iknot = 0_ip         !knot sequence chosen by db3ink

    call db3ink(x,nx,y,ny,z,nz,&
                fcn,&
                kx,ky,kz,&
                iknot,&
                me%tx,me%ty,me%tz,&
                me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_3d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_3d]] type (with user-specified knots).
!  This is a wrapper for [[db3ink]].

    pure subroutine initialize_3d_specify_knots(me,x,y,z,fcn,kx,ky,kz,tx,ty,tz,iflag,extrap)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer(ip),intent(in)               :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer(ip),intent(in)               :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)     :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)     :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)     :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    integer(ip),intent(out)  :: iflag               !! status flag (see [[db3ink]])
    logical,intent(in),optional  :: extrap      !! if true, then extrapolation is allowed
                                                !! (default is false)

    integer(ip) :: nx,ny,nz

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  nz=nz,kz=kz,tz=tz,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny
        me%nz = nz

        me%kx = kx
        me%ky = ky
        me%kz = kz

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%tz(nz+kz))
        allocate(me%bcoef(nx,ny,nz))
        allocate(me%work_val_1(ky,kz))
        allocate(me%work_val_2(kz))
        allocate(me%work_val_3(3_ip*max(kx,ky,kz)))

        me%tx = tx
        me%ty = ty
        me%tz = tz

        call db3ink(x,nx,y,ny,z,nz,&
                    fcn,&
                    kx,ky,kz,&
                    1_ip,&
                    me%tx,me%ty,me%tz,&
                    me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_3d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_3d]] interpolate.  This is a wrapper for [[db3val]].

    pure subroutine evaluate_3d(me,xval,yval,zval,idx,idy,idz,f,iflag)

    implicit none

    class(bspline_3d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)             :: zval  !! \(z\) coordinate of evaluation point.
    integer(ip),intent(in)          :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db3val]])

    if (me%initialized) then
        call db3val(xval,yval,zval,&
                    idx,idy,idz,&
                    me%tx,me%ty,me%tz,&
                    me%nx,me%ny,me%nz,&
                    me%kx,me%ky,me%kz,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,&
                    me%iloy,me%iloz,&
                    me%work_val_1,me%work_val_2,me%work_val_3,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_4d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_4d_constructor_empty() result(me)

    implicit none

    type(bspline_4d) :: me

    end function bspline_4d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_4d]] type (auto knots).
!  This is a wrapper for [[initialize_4d_auto_knots]].

    pure function bspline_4d_constructor_auto_knots(x,y,z,q,fcn,kx,ky,kz,kq,extrap) result(me)

    implicit none

    type(bspline_4d)                       :: me
    real(wp),dimension(:),intent(in)       :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                  !! `fcn(i,j,k,l)` should contain the function value at the
                                                  !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer(ip),intent(in)                 :: kx  !! The order of spline pieces in \(x\)
                                                  !! ( \( 2 \le k_x < n_x \) )
                                                  !! (order = polynomial degree + 1)
    integer(ip),intent(in)                 :: ky  !! The order of spline pieces in \(y\)
                                                  !! ( \( 2 \le k_y < n_y \) )
                                                  !! (order = polynomial degree + 1)
    integer(ip),intent(in)                 :: kz  !! The order of spline pieces in \(z\)
                                                  !! ( \( 2 \le k_z < n_z \) )
                                                  !! (order = polynomial degree + 1)
    integer(ip),intent(in)                 :: kq  !! The order of spline pieces in \(q\)
                                                  !! ( \( 2 \le k_q < n_q \) )
                                                  !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap    !! if true, then extrapolation is allowed
                                                  !! (default is false)

    call initialize_4d_auto_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,me%iflag,extrap)

    end function bspline_4d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_4d]] type (user-specified knots).
!  This is a wrapper for [[initialize_4d_specify_knots]].

    pure function bspline_4d_constructor_specify_knots(x,y,z,q,fcn,kx,ky,kz,kq,&
                                                        tx,ty,tz,tq,extrap) result(me)

    implicit none

    type(bspline_4d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap        !! if true, then extrapolation is allowed
                                                      !! (default is false)

    call initialize_4d_specify_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,tx,ty,tz,tq,me%iflag,extrap)

    end function bspline_4d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_4d]] type (with automatically-computed knots).
!  This is a wrapper for [[db4ink]].

    pure subroutine initialize_4d_auto_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,iflag,extrap)

    implicit none

    class(bspline_4d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db4ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny,nz,nq

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)

    me%nx = nx
    me%ny = ny
    me%nz = nz
    me%nq = nq

    me%kx = kx
    me%ky = ky
    me%kz = kz
    me%kq = kq

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%tz(nz+kz))
    allocate(me%tq(nq+kq))
    allocate(me%bcoef(nx,ny,nz,nq))
    allocate(me%work_val_1(ky,kz,kq))
    allocate(me%work_val_2(kz,kq))
    allocate(me%work_val_3(kq))
    allocate(me%work_val_4(3_ip*max(kx,ky,kz,kq)))

    iknot = 0_ip         !knot sequence chosen by db4ink

    call db4ink(x,nx,y,ny,z,nz,q,nq,&
                fcn,&
                kx,ky,kz,kq,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,&
                me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_4d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_4d]] type (with user-specified knots).
!  This is a wrapper for [[db4ink]].

    pure subroutine initialize_4d_specify_knots(me,x,y,z,q,fcn,&
                                                kx,ky,kz,kq,tx,ty,tz,tq,iflag,extrap)

    implicit none

    class(bspline_4d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db4ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: nx,ny,nz,nq

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  nz=nz,kz=kz,tz=tz,&
                                  nq=nq,kq=kq,tq=tq,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny
        me%nz = nz
        me%nq = nq

        me%kx = kx
        me%ky = ky
        me%kz = kz
        me%kq = kq

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%tz(nz+kz))
        allocate(me%tq(nq+kq))
        allocate(me%bcoef(nx,ny,nz,nq))
        allocate(me%work_val_1(ky,kz,kq))
        allocate(me%work_val_2(kz,kq))
        allocate(me%work_val_3(kq))
        allocate(me%work_val_4(3_ip*max(kx,ky,kz,kq)))

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq

        call db4ink(x,nx,y,ny,z,nz,q,nq,&
                    fcn,&
                    kx,ky,kz,kq,&
                    1_ip,&
                    me%tx,me%ty,me%tz,me%tq,&
                    me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_4d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_4d]] interpolate.  This is a wrapper for [[db4val]].

    pure subroutine evaluate_4d(me,xval,yval,zval,qval,idx,idy,idz,idq,f,iflag)

    implicit none

    class(bspline_4d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)             :: zval  !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)             :: qval  !! \(q\) coordinate of evaluation point.
    integer(ip),intent(in)          :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db4val]])

    if (me%initialized) then
        call db4val(xval,yval,zval,qval,&
                    idx,idy,idz,idq,&
                    me%tx,me%ty,me%tz,me%tq,&
                    me%nx,me%ny,me%nz,me%nq,&
                    me%kx,me%ky,me%kz,me%kq,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,&
                    me%iloy,me%iloz,me%iloq,&
                    me%work_val_1,me%work_val_2,me%work_val_3,me%work_val_4,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_5d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_5d_constructor_empty() result(me)

    implicit none

    type(bspline_5d) :: me

    end function bspline_5d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_5d]] type (auto knots).
!  This is a wrapper for [[initialize_5d_auto_knots]].

    pure function bspline_5d_constructor_auto_knots(x,y,z,q,r,fcn,kx,ky,kz,kq,kr,extrap) result(me)

    implicit none

    type(bspline_5d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)   :: fcn !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap        !! if true, then extrapolation is allowed
                                                      !! (default is false)

    call initialize_5d_auto_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,me%iflag,extrap)

    end function bspline_5d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_5d]] type (user-specified knots).
!  This is a wrapper for [[initialize_5d_specify_knots]].

    pure function bspline_5d_constructor_specify_knots(x,y,z,q,r,fcn,&
                                                        kx,ky,kz,kq,kr,&
                                                        tx,ty,tz,tq,tr,extrap) result(me)

    implicit none

    type(bspline_5d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)   :: fcn !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tr  !! The `(nr+kr)` knots in the \(r\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap        !! if true, then extrapolation is allowed
                                                      !! (default is false)

    call initialize_5d_specify_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,me%iflag,extrap)

    end function bspline_5d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_5d]] type (with automatically-computed knots).
!  This is a wrapper for [[db5ink]].

    pure subroutine initialize_5d_auto_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,iflag,extrap)

    implicit none

    class(bspline_5d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)   :: fcn !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db5ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny,nz,nq,nr

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)
    nr = size(r,kind=ip)

    me%nx = nx
    me%ny = ny
    me%nz = nz
    me%nq = nq
    me%nr = nr

    me%kx = kx
    me%ky = ky
    me%kz = kz
    me%kq = kq
    me%kr = kr

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%tz(nz+kz))
    allocate(me%tq(nq+kq))
    allocate(me%tr(nr+kr))
    allocate(me%bcoef(nx,ny,nz,nq,nr))
    allocate(me%work_val_1(ky,kz,kq,kr))
    allocate(me%work_val_2(kz,kq,kr))
    allocate(me%work_val_3(kq,kr))
    allocate(me%work_val_4(kr))
    allocate(me%work_val_5(3_ip*max(kx,ky,kz,kq,kr)))

    iknot = 0_ip         !knot sequence chosen by db5ink

    call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                fcn,&
                kx,ky,kz,kq,kr,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,me%tr,&
                me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_5d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_5d]] type (with user-specified knots).
!  This is a wrapper for [[db5ink]].

    pure subroutine initialize_5d_specify_knots(me,x,y,z,q,r,fcn,&
                                                kx,ky,kz,kq,kr,&
                                                tx,ty,tz,tq,tr,iflag,extrap)

    implicit none

    class(bspline_5d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:),intent(in)   :: fcn !! `(nx,ny,nz,nq,nr)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tr  !! The `(nr+kr)` knots in the \(r\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db5ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: nx,ny,nz,nq,nr

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)
    nr = size(r,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  nz=nz,kz=kz,tz=tz,&
                                  nq=nq,kq=kq,tq=tq,&
                                  nr=nr,kr=kr,tr=tr,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny
        me%nz = nz
        me%nq = nq
        me%nr = nr

        me%kx = kx
        me%ky = ky
        me%kz = kz
        me%kq = kq
        me%kr = kr

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%tz(nz+kz))
        allocate(me%tq(nq+kq))
        allocate(me%tr(nr+kr))
        allocate(me%bcoef(nx,ny,nz,nq,nr))
        allocate(me%work_val_1(ky,kz,kq,kr))
        allocate(me%work_val_2(kz,kq,kr))
        allocate(me%work_val_3(kq,kr))
        allocate(me%work_val_4(kr))
        allocate(me%work_val_5(3_ip*max(kx,ky,kz,kq,kr)))

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq
        me%tr = tr

        call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                    fcn,&
                    kx,ky,kz,kq,kr,&
                    1_ip,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,&
                    me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_5d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_5d]] interpolate.  This is a wrapper for [[db5val]].

    pure subroutine evaluate_5d(me,xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,f,iflag)

    implicit none

    class(bspline_5d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)             :: zval  !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)             :: qval  !! \(q\) coordinate of evaluation point.
    real(wp),intent(in)             :: rval  !! \(r\) coordinate of evaluation point.
    integer(ip),intent(in)          :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idr   !! \(r\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db5val]])

    if (me%initialized) then
        call db5val(xval,yval,zval,qval,rval,&
                    idx,idy,idz,idq,idr,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,&
                    me%nx,me%ny,me%nz,me%nq,me%nr,&
                    me%kx,me%ky,me%kz,me%kq,me%kr,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,&
                    me%iloy,me%iloz,me%iloq,me%ilor,&
                    me%work_val_1,me%work_val_2,me%work_val_3,me%work_val_4,me%work_val_5,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  It returns an empty [[bspline_6d]] type. Note that INITIALIZE still
!  needs to be called before it can be used.
!  Not really that useful except perhaps in some OpenMP applications.

    elemental function bspline_6d_constructor_empty() result(me)

    implicit none

    type(bspline_6d) :: me

    end function bspline_6d_constructor_empty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_6d]] type (auto knots).
!  This is a wrapper for [[initialize_6d_auto_knots]].

    pure function bspline_6d_constructor_auto_knots(x,y,z,q,r,s,fcn,&
                                                    kx,ky,kz,kq,kr,ks,extrap) result(me)

    implicit none

    type(bspline_6d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: s   !! `(ns)` array of \(s\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m,n)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer(ip),intent(in)                      :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                      :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    logical,intent(in),optional      :: extrap        !! if true, then extrapolation is allowed
                                                      !! (default is false)

    call initialize_6d_auto_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,me%iflag,extrap)

    end function bspline_6d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_6d]] type (user-specified knots).
!  This is a wrapper for [[initialize_6d_specify_knots]].

    pure function bspline_6d_constructor_specify_knots(x,y,z,q,r,s,fcn,&
                                                        kx,ky,kz,kq,kr,ks,&
                                                        tx,ty,tz,tq,tr,ts,extrap) result(me)

    implicit none

    type(bspline_6d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: s   !! `(ns)` array of \(s\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m,n)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tr  !! The `(nr+kr)` knots in the \(r\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ts  !! The `(ns+ks)` knots in the \(s\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    logical,intent(in),optional      :: extrap        !! if true, then extrapolation is allowed
                                                      !! (default is false)

    call initialize_6d_specify_knots(me,x,y,z,q,r,s,fcn,&
                                        kx,ky,kz,kq,kr,ks,&
                                        tx,ty,tz,tq,tr,ts,me%iflag,extrap)

    end function bspline_6d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_6d]] type (with automatically-computed knots).
!  This is a wrapper for [[db6ink]].

    pure subroutine initialize_6d_auto_knots(me,x,y,z,q,r,s,fcn,&
                                                kx,ky,kz,kq,kr,ks,iflag,extrap)

    implicit none

    class(bspline_6d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: s   !! `(ns)` array of \(s\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m,n)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db6ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: iknot
    integer(ip) :: nx,ny,nz,nq,nr,ns

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)
    nr = size(r,kind=ip)
    ns = size(s,kind=ip)

    me%nx = nx
    me%ny = ny
    me%nz = nz
    me%nq = nq
    me%nr = nr
    me%ns = ns

    me%kx = kx
    me%ky = ky
    me%kz = kz
    me%kq = kq
    me%kr = kr
    me%ks = ks

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%tz(nz+kz))
    allocate(me%tq(nq+kq))
    allocate(me%tr(nr+kr))
    allocate(me%ts(ns+ks))
    allocate(me%bcoef(nx,ny,nz,nq,nr,ns))
    allocate(me%work_val_1(ky,kz,kq,kr,ks))
    allocate(me%work_val_2(kz,kq,kr,ks))
    allocate(me%work_val_3(kq,kr,ks))
    allocate(me%work_val_4(kr,ks))
    allocate(me%work_val_5(ks))
    allocate(me%work_val_6(3_ip*max(kx,ky,kz,kq,kr,ks)))

    iknot = 0_ip         !knot sequence chosen by db6ink

    call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                fcn,&
                kx,ky,kz,kq,kr,ks,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                me%bcoef,iflag)

    if (iflag==0_ip) then
        call me%set_extrap_flag(extrap)
    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_6d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_6d]] type (with user-specified knots).
!  This is a wrapper for [[db6ink]].

    pure subroutine initialize_6d_specify_knots(me,x,y,z,q,r,s,fcn,&
                                                kx,ky,kz,kq,kr,ks,&
                                                tx,ty,tz,tq,tr,ts,iflag,extrap)

    implicit none

    class(bspline_6d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: r   !! `(nr)` array of \(r\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: s   !! `(ns)` array of \(s\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq,nr,ns)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l,m,n)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`,`r(m)`,`s(n)`)
    integer(ip),intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer(ip),intent(in)                     :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)           :: tx  !! The `(nx+kx)` knots in the \(x\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ty  !! The `(ny+ky)` knots in the \(y\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tz  !! The `(nz+kz)` knots in the \(z\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tq  !! The `(nq+kq)` knots in the \(q\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: tr  !! The `(nr+kr)` knots in the \(r\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)           :: ts  !! The `(ns+ks)` knots in the \(s\) direction
                                                      !! for the spline interpolant.
                                                      !! Must be non-decreasing.
    integer(ip),intent(out)  :: iflag                     !! status flag (see [[db6ink]])
    logical,intent(in),optional  :: extrap            !! if true, then extrapolation is allowed
                                                      !! (default is false)

    integer(ip) :: nx,ny,nz,nq,nr,ns

    call me%destroy()

    nx = size(x,kind=ip)
    ny = size(y,kind=ip)
    nz = size(z,kind=ip)
    nq = size(q,kind=ip)
    nr = size(r,kind=ip)
    ns = size(s,kind=ip)

    call check_knot_vectors_sizes(nx=nx,kx=kx,tx=tx,&
                                  ny=ny,ky=ky,ty=ty,&
                                  nz=nz,kz=kz,tz=tz,&
                                  nq=nq,kq=kq,tq=tq,&
                                  nr=nr,kr=kr,tr=tr,&
                                  ns=ns,ks=ks,ts=ts,&
                                  iflag=iflag)

    if (iflag == 0_ip) then

        me%nx = nx
        me%ny = ny
        me%nz = nz
        me%nq = nq
        me%nr = nr
        me%ns = ns

        me%kx = kx
        me%ky = ky
        me%kz = kz
        me%kq = kq
        me%kr = kr
        me%ks = ks

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%tz(nz+kz))
        allocate(me%tq(nq+kq))
        allocate(me%tr(nr+kr))
        allocate(me%ts(ns+ks))
        allocate(me%bcoef(nx,ny,nz,nq,nr,ns))
        allocate(me%work_val_1(ky,kz,kq,kr,ks))
        allocate(me%work_val_2(kz,kq,kr,ks))
        allocate(me%work_val_3(kq,kr,ks))
        allocate(me%work_val_4(kr,ks))
        allocate(me%work_val_5(ks))
        allocate(me%work_val_6(3_ip*max(kx,ky,kz,kq,kr,ks)))

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq
        me%tr = tr
        me%ts = ts

        call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                    fcn,&
                    kx,ky,kz,kq,kr,ks,&
                    1_ip,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                    me%bcoef,iflag)

        call me%set_extrap_flag(extrap)

    end if

    me%initialized = iflag==0_ip
    me%iflag = iflag

    end subroutine initialize_6d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a [[bspline_6d]] interpolate.  This is a wrapper for [[db6val]].

    pure subroutine evaluate_6d(me,xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,f,iflag)

    implicit none

    class(bspline_6d),intent(inout) :: me
    real(wp),intent(in)             :: xval  !! \(x\) coordinate of evaluation point.
    real(wp),intent(in)             :: yval  !! \(y\) coordinate of evaluation point.
    real(wp),intent(in)             :: zval  !! \(z\) coordinate of evaluation point.
    real(wp),intent(in)             :: qval  !! \(q\) coordinate of evaluation point.
    real(wp),intent(in)             :: rval  !! \(r\) coordinate of evaluation point.
    real(wp),intent(in)             :: sval  !! \(s\) coordinate of evaluation point.
    integer(ip),intent(in)          :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: idr   !! \(r\) derivative of piecewise polynomial to evaluate.
    integer(ip),intent(in)          :: ids   !! \(s\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer(ip),intent(out)         :: iflag !! status flag (see [[db6val]])

    if (me%initialized) then
        call db6val(xval,yval,zval,qval,rval,sval,&
                    idx,idy,idz,idq,idr,ids,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                    me%nx,me%ny,me%nz,me%nq,me%nr,me%ns,&
                    me%kx,me%ky,me%kz,me%kq,me%kr,me%ks,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,me%inbvs,&
                    me%iloy,me%iloz,me%iloq,me%ilor,me%ilos,&
                    me%work_val_1,me%work_val_2,me%work_val_3,me%work_val_4,me%work_val_5,me%work_val_6,&
                    extrap=me%extrap)
    else
        iflag = 1_ip
    end if

    me%iflag = iflag

    end subroutine evaluate_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Error checks for the user-specified knot vector sizes.
!
!@note If more than one is the wrong size, then the `iflag` error code will
!      correspond to the one with the highest rank.

    pure subroutine check_knot_vectors_sizes(nx,ny,nz,nq,nr,ns,&
                                             kx,ky,kz,kq,kr,ks,&
                                             tx,ty,tz,tq,tr,ts,iflag)

    implicit none

    integer(ip),intent(in),optional           :: nx
    integer(ip),intent(in),optional           :: ny
    integer(ip),intent(in),optional           :: nz
    integer(ip),intent(in),optional           :: nq
    integer(ip),intent(in),optional           :: nr
    integer(ip),intent(in),optional           :: ns
    integer(ip),intent(in),optional           :: kx
    integer(ip),intent(in),optional           :: ky
    integer(ip),intent(in),optional           :: kz
    integer(ip),intent(in),optional           :: kq
    integer(ip),intent(in),optional           :: kr
    integer(ip),intent(in),optional           :: ks
    real(wp),dimension(:),intent(in),optional :: tx
    real(wp),dimension(:),intent(in),optional :: ty
    real(wp),dimension(:),intent(in),optional :: tz
    real(wp),dimension(:),intent(in),optional :: tq
    real(wp),dimension(:),intent(in),optional :: tr
    real(wp),dimension(:),intent(in),optional :: ts
    integer(ip),intent(out)                   :: iflag  !! 0 if everything is OK

    iflag = 0_ip

    if (present(nx) .and. present(kx) .and. present(tx)) then
        if (size(tx,kind=ip)/=(nx+kx)) then
            iflag = 501_ip  ! tx is not the correct size (nx+kx)
        end if
    end if

    if (present(ny) .and. present(ky) .and. present(ty)) then
        if (size(ty,kind=ip)/=(ny+ky)) then
            iflag = 502_ip  ! ty is not the correct size (ny+ky)
        end if
    end if

    if (present(nz) .and. present(kz) .and. present(tz)) then
        if (size(tz,kind=ip)/=(nz+kz)) then
            iflag = 503_ip  ! tz is not the correct size (nz+kz)
        end if
    end if

    if (present(nq) .and. present(kq) .and. present(tq)) then
        if (size(tq,kind=ip)/=(nq+kq)) then
            iflag = 504_ip  ! tq is not the correct size (nq+kq)
        end if
    end if

    if (present(nr) .and. present(kr) .and. present(tr)) then
        if (size(tr,kind=ip)/=(nr+kr)) then
            iflag = 505_ip  ! tr is not the correct size (nr+kr)
        end if
    end if

    if (present(ns) .and. present(ks) .and. present(ts)) then
        if (size(ts,kind=ip)/=(ns+ks)) then
            iflag = 506_ip  ! ts is not the correct size (ns+ks)
        end if
    end if

    end subroutine check_knot_vectors_sizes
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_oo_module
!*****************************************************************************************
