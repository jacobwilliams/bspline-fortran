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

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use,intrinsic :: iso_fortran_env, only: error_unit
    use bspline_sub_module

    implicit none

    private

    integer,parameter :: int_size     = storage_size(1)       !! size of a default integer [bits]
    integer,parameter :: logical_size = storage_size(.true.)  !! size of a default logical [bits]
    integer,parameter :: real_size    = storage_size(1.0_wp)  !! size of a `real(wp)` [bits]

    type,public,abstract :: bspline_class
        !! Base class for the b-spline types
        private
        integer :: inbvx = 1  !! internal variable used by [[dbvalu]] for efficient processing
        integer :: iflag = 1  !! saved `iflag` from the list routine call.
        logical :: initialized = .false. !! true if the class is initialized and ready to use
    contains
        private
        procedure,non_overridable :: destroy_base  !! destructor for the abstract type
        procedure(destroy_func),deferred,public :: destroy  !! destructor
        procedure(size_func),deferred,public :: size_of !! size of the structure in bits
        procedure,public,non_overridable :: status_ok  !! returns true if the last `iflag` status code was `=0`.
        procedure,public,non_overridable :: status_message => get_bspline_status_message  !! retrieve the last status message
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
        import :: bspline_class
        implicit none
        class(bspline_class),intent(in) :: me
        integer :: s !! size of the structure in bits
        end function size_func

    end interface

    type,extends(bspline_class),public :: bspline_1d
        !! Class for 1d b-spline interpolation.
        private
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        real(wp),dimension(:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        contains
        private
        generic,public :: initialize => initialize_1d_auto_knots,initialize_1d_specify_knots
        procedure :: initialize_1d_auto_knots
        procedure :: initialize_1d_specify_knots
        procedure,public :: evaluate => evaluate_1d
        procedure,public :: destroy => destroy_1d
        procedure,public :: size_of => size_1d
        final :: finalize_1d
    end type bspline_1d

    type,extends(bspline_class),public :: bspline_2d
        !! Class for 2d b-spline interpolation.
        private
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: ny  = 0  !! Number of \(y\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        integer :: ky  = 0  !! The order of spline pieces in \(y\)
        real(wp),dimension(:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        integer :: inbvy = 1  !! internal variable used for efficient processing
        integer :: iloy = 1  !! internal variable used for efficient processing
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
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: ny  = 0  !! Number of \(y\) abcissae
        integer :: nz  = 0  !! Number of \(z\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        integer :: ky  = 0  !! The order of spline pieces in \(y\)
        integer :: kz  = 0  !! The order of spline pieces in \(z\)
        real(wp),dimension(:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        integer :: inbvy = 1  !! internal variable used for efficient processing
        integer :: inbvz = 1  !! internal variable used for efficient processing
        integer :: iloy = 1  !! internal variable used for efficient processing
        integer :: iloz = 1  !! internal variable used for efficient processing
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
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: ny  = 0  !! Number of \(y\) abcissae
        integer :: nz  = 0  !! Number of \(z\) abcissae
        integer :: nq  = 0  !! Number of \(q\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        integer :: ky  = 0  !! The order of spline pieces in \(y\)
        integer :: kz  = 0  !! The order of spline pieces in \(z\)
        integer :: kq  = 0  !! The order of spline pieces in \(q\)
        real(wp),dimension(:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        integer :: inbvy = 1  !! internal variable used for efficient processing
        integer :: inbvz = 1  !! internal variable used for efficient processing
        integer :: inbvq = 1  !! internal variable used for efficient processing
        integer :: iloy  = 1  !! internal variable used for efficient processing
        integer :: iloz  = 1  !! internal variable used for efficient processing
        integer :: iloq  = 1  !! internal variable used for efficient processing
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
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: ny  = 0  !! Number of \(y\) abcissae
        integer :: nz  = 0  !! Number of \(z\) abcissae
        integer :: nq  = 0  !! Number of \(q\) abcissae
        integer :: nr  = 0  !! Number of \(r\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        integer :: ky  = 0  !! The order of spline pieces in \(y\)
        integer :: kz  = 0  !! The order of spline pieces in \(z\)
        integer :: kq  = 0  !! The order of spline pieces in \(q\)
        integer :: kr  = 0  !! The order of spline pieces in \(r\)
        real(wp),dimension(:,:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tr  !! The knots in the \(r\) direction for the spline interpolant
        integer :: inbvy = 1  !! internal variable used for efficient processing
        integer :: inbvz = 1  !! internal variable used for efficient processing
        integer :: inbvq = 1  !! internal variable used for efficient processing
        integer :: inbvr = 1  !! internal variable used for efficient processing
        integer :: iloy  = 1  !! internal variable used for efficient processing
        integer :: iloz  = 1  !! internal variable used for efficient processing
        integer :: iloq  = 1  !! internal variable used for efficient processing
        integer :: ilor  = 1  !! internal variable used for efficient processing
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
        integer :: nx  = 0  !! Number of \(x\) abcissae
        integer :: ny  = 0  !! Number of \(y\) abcissae
        integer :: nz  = 0  !! Number of \(z\) abcissae
        integer :: nq  = 0  !! Number of \(q\) abcissae
        integer :: nr  = 0  !! Number of \(r\) abcissae
        integer :: ns  = 0  !! Number of \(s\) abcissae
        integer :: kx  = 0  !! The order of spline pieces in \(x\)
        integer :: ky  = 0  !! The order of spline pieces in \(y\)
        integer :: kz  = 0  !! The order of spline pieces in \(z\)
        integer :: kq  = 0  !! The order of spline pieces in \(q\)
        integer :: kr  = 0  !! The order of spline pieces in \(r\)
        integer :: ks  = 0  !! The order of spline pieces in \(s\)
        real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef  !! array of coefficients of the b-spline interpolant
        real(wp),dimension(:),allocatable :: tx  !! The knots in the \(x\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ty  !! The knots in the \(y\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tz  !! The knots in the \(z\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tq  !! The knots in the \(q\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: tr  !! The knots in the \(r\) direction for the spline interpolant
        real(wp),dimension(:),allocatable :: ts  !! The knots in the \(s\) direction for the spline interpolant
        integer :: inbvy = 1  !! internal variable used for efficient processing
        integer :: inbvz = 1  !! internal variable used for efficient processing
        integer :: inbvq = 1  !! internal variable used for efficient processing
        integer :: inbvr = 1  !! internal variable used for efficient processing
        integer :: inbvs = 1  !! internal variable used for efficient processing
        integer :: iloy  = 1  !! internal variable used for efficient processing
        integer :: iloz  = 1  !! internal variable used for efficient processing
        integer :: iloq  = 1  !! internal variable used for efficient processing
        integer :: ilor  = 1  !! internal variable used for efficient processing
        integer :: ilos  = 1  !! internal variable used for efficient processing
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

    ok = ( me%iflag == 0 )

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

    me%iflag = 0

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
    integer,intent(in),optional     :: iflag  !! the corresponding status code

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
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 2*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef)
    if (allocated(me%tx))    s = s + real_size*size(me%tx)

    end function size_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_2d]] structure in bits.

    pure function size_2d(me) result(s)

    implicit none

    class(bspline_2d),intent(in) :: me
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 6*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef,1)*&
                                               size(me%bcoef,2)
    if (allocated(me%tx)) s = s + real_size*size(me%tx)
    if (allocated(me%ty)) s = s + real_size*size(me%ty)

    end function size_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_3d]] structure in bits.

    pure function size_3d(me) result(s)

    implicit none

    class(bspline_3d),intent(in) :: me
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 10*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef,1)*&
                                               size(me%bcoef,2)*&
                                               size(me%bcoef,3)
    if (allocated(me%tx)) s = s + real_size*size(me%tx)
    if (allocated(me%ty)) s = s + real_size*size(me%ty)
    if (allocated(me%tz)) s = s + real_size*size(me%tz)

    end function size_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_4d]] structure in bits.

    pure function size_4d(me) result(s)

    implicit none

    class(bspline_4d),intent(in) :: me
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 14*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef,1)*&
                                               size(me%bcoef,2)*&
                                               size(me%bcoef,3)*&
                                               size(me%bcoef,4)
    if (allocated(me%tx)) s = s + real_size*size(me%tx)
    if (allocated(me%ty)) s = s + real_size*size(me%ty)
    if (allocated(me%tz)) s = s + real_size*size(me%tz)
    if (allocated(me%tq)) s = s + real_size*size(me%tq)

    end function size_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_5d]] structure in bits.

    pure function size_5d(me) result(s)

    implicit none

    class(bspline_5d),intent(in) :: me
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 18*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef,1)*&
                                               size(me%bcoef,2)*&
                                               size(me%bcoef,3)*&
                                               size(me%bcoef,4)*&
                                               size(me%bcoef,5)
    if (allocated(me%tx)) s = s + real_size*size(me%tx)
    if (allocated(me%ty)) s = s + real_size*size(me%ty)
    if (allocated(me%tz)) s = s + real_size*size(me%tz)
    if (allocated(me%tq)) s = s + real_size*size(me%tq)
    if (allocated(me%tr)) s = s + real_size*size(me%tr)

    end function size_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Actual size of a [[bspline_6d]] structure in bits.

    pure function size_6d(me) result(s)

    implicit none

    class(bspline_6d),intent(in) :: me
    integer :: s !! size of the structure in bits

    s = 2*int_size + logical_size + 22*int_size

    if (allocated(me%bcoef)) s = s + real_size*size(me%bcoef,1)*&
                                               size(me%bcoef,2)*&
                                               size(me%bcoef,3)*&
                                               size(me%bcoef,4)*&
                                               size(me%bcoef,5)*&
                                               size(me%bcoef,6)
    if (allocated(me%tx)) s = s + real_size*size(me%tx)
    if (allocated(me%ty)) s = s + real_size*size(me%ty)
    if (allocated(me%tz)) s = s + real_size*size(me%tz)
    if (allocated(me%tq)) s = s + real_size*size(me%tq)
    if (allocated(me%tr)) s = s + real_size*size(me%tr)
    if (allocated(me%ts)) s = s + real_size*size(me%ts)

    end function size_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for contents of the base [[bspline_class]] class.
!  (this routine is called by the extended classes).

    pure subroutine destroy_base(me)

    implicit none

    class(bspline_class),intent(inout) :: me

    me%inbvx = 1
    me%iflag = 1
    me%initialized = .false.

    end subroutine destroy_base
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_1d]] class.

    pure subroutine destroy_1d(me)

    implicit none

    class(bspline_1d),intent(inout) :: me

    call me%destroy_base()

    me%nx  = 0
    me%kx  = 0
    if (allocated(me%bcoef))    deallocate(me%bcoef)
    if (allocated(me%tx))       deallocate(me%tx)

    end subroutine destroy_1d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_2d]] class.

    pure subroutine destroy_2d(me)

    implicit none

    class(bspline_2d),intent(inout) :: me

    call me%destroy_base()

    me%nx    = 0
    me%ny    = 0
    me%kx    = 0
    me%ky    = 0
    me%inbvy = 1
    me%iloy  = 1
    if (allocated(me%bcoef))    deallocate(me%bcoef)
    if (allocated(me%tx))       deallocate(me%tx)
    if (allocated(me%ty))       deallocate(me%ty)

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_3d]] class.

    pure subroutine destroy_3d(me)

    implicit none

    class(bspline_3d),intent(inout) :: me

    call me%destroy_base()

    me%nx    = 0
    me%ny    = 0
    me%nz    = 0
    me%kx    = 0
    me%ky    = 0
    me%kz    = 0
    me%inbvy = 1
    me%inbvz = 1
    me%iloy  = 1
    me%iloz  = 1
    if (allocated(me%bcoef))    deallocate(me%bcoef)
    if (allocated(me%tx))       deallocate(me%tx)
    if (allocated(me%ty))       deallocate(me%ty)
    if (allocated(me%tz))       deallocate(me%tz)

    end subroutine destroy_3d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_4d]] class.

    pure subroutine destroy_4d(me)

    implicit none

    class(bspline_4d),intent(inout) :: me

    me%nx    = 0
    me%ny    = 0
    me%nz    = 0
    me%nq    = 0
    me%kx    = 0
    me%ky    = 0
    me%kz    = 0
    me%kq    = 0
    me%inbvy = 1
    me%inbvz = 1
    me%inbvq = 1
    me%iloy  = 1
    me%iloz  = 1
    me%iloq  = 1
    if (allocated(me%bcoef))   deallocate(me%bcoef)
    if (allocated(me%tx))      deallocate(me%tx)
    if (allocated(me%ty))      deallocate(me%ty)
    if (allocated(me%tz))      deallocate(me%tz)
    if (allocated(me%tq))      deallocate(me%tq)

    end subroutine destroy_4d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_5d]] class.

    pure subroutine destroy_5d(me)

    implicit none

    class(bspline_5d),intent(inout) :: me

    me%nx    = 0
    me%ny    = 0
    me%nz    = 0
    me%nq    = 0
    me%nr    = 0
    me%kx    = 0
    me%ky    = 0
    me%kz    = 0
    me%kq    = 0
    me%kr    = 0
    me%inbvy = 1
    me%inbvz = 1
    me%inbvq = 1
    me%inbvr = 1
    me%iloy  = 1
    me%iloz  = 1
    me%iloq  = 1
    me%ilor  = 1
    if (allocated(me%bcoef))    deallocate(me%bcoef)
    if (allocated(me%tx))       deallocate(me%tx)
    if (allocated(me%ty))       deallocate(me%ty)
    if (allocated(me%tz))       deallocate(me%tz)
    if (allocated(me%tq))       deallocate(me%tq)
    if (allocated(me%tr))       deallocate(me%tr)

    end subroutine destroy_5d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[bspline_6d]] class.

    pure subroutine destroy_6d(me)

    implicit none

    class(bspline_6d),intent(inout) :: me

    me%nx    = 0
    me%ny    = 0
    me%nz    = 0
    me%nq    = 0
    me%nr    = 0
    me%ns    = 0
    me%kx    = 0
    me%ky    = 0
    me%kz    = 0
    me%kq    = 0
    me%kr    = 0
    me%ks    = 0
    me%inbvy = 1
    me%inbvz = 1
    me%inbvq = 1
    me%inbvr = 1
    me%inbvs = 1
    me%iloy  = 1
    me%iloz  = 1
    me%iloq  = 1
    me%ilor  = 1
    me%ilos  = 1
    if (allocated(me%bcoef))    deallocate(me%bcoef)
    if (allocated(me%tx))       deallocate(me%tx)
    if (allocated(me%ty))       deallocate(me%ty)
    if (allocated(me%tz))       deallocate(me%tz)
    if (allocated(me%tq))       deallocate(me%tq)
    if (allocated(me%tr))       deallocate(me%tr)
    if (allocated(me%ts))       deallocate(me%ts)

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

    pure function bspline_1d_constructor_auto_knots(x,fcn,kx) result(me)

    implicit none

    type(bspline_1d)                 :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer,intent(in)               :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)

    call initialize_1d_auto_knots(me,x,fcn,kx,me%iflag)

    end function bspline_1d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_1d]] type (user-specified knots).
!  This is a wrapper for [[initialize_1d_specify_knots]].

    pure function bspline_1d_constructor_specify_knots(x,fcn,kx,tx) result(me)

    implicit none

    type(bspline_1d)                 :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer,intent(in)               :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in) :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                              !! for the spline interpolant.
                                              !! Must be non-decreasing.

    call initialize_1d_specify_knots(me,x,fcn,kx,tx,me%iflag)

    end function bspline_1d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_1d]] type (with automatically-computed knots).
!  This is a wrapper for [[db1ink]].

    pure subroutine initialize_1d_auto_knots(me,x,fcn,kx,iflag)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer,intent(in)               :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    integer,intent(out)              :: iflag !! status flag (see [[db1ink]])

    integer :: iknot
    integer :: nx

    call me%destroy()

    nx = size(x)

    me%nx = nx
    me%kx = kx

    allocate(me%tx(nx+kx))
    allocate(me%bcoef(nx))

    iknot = 0         !knot sequence chosen by db1ink

    call db1ink(x,nx,fcn,kx,iknot,me%tx,me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_1d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_1d]] type (with user-specified knots).
!  This is a wrapper for [[db1ink]].

    pure subroutine initialize_1d_specify_knots(me,x,fcn,kx,tx,iflag)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in) :: fcn   !! `(nx)` array of function values to interpolate. `fcn(i)` should
                                              !! contain the function value at the point `x(i)`
    integer,intent(in)               :: kx    !! The order of spline pieces in \(x\)
                                              !! ( \( 2 \le k_x < n_x \) )
                                              !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in) :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                              !! for the spline interpolant.
                                              !! Must be non-decreasing.
    integer,intent(out)              :: iflag !! status flag (see [[db1ink]])

    integer :: nx

    call me%destroy()

    nx = size(x)

    call check_knot_vectors_sizes('initialize_1d_specify_knots',nx=nx,kx=kx,tx=tx,iflag=iflag)

    if (iflag == 0) then

        me%nx = nx
        me%kx = kx

        allocate(me%tx(nx+kx))
        allocate(me%bcoef(nx))

        me%tx = tx

        call db1ink(x,nx,fcn,kx,1,me%tx,me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db1val]])

    if (me%initialized) then
        call db1val(xval,idx,me%tx,me%nx,me%kx,me%bcoef,f,iflag,me%inbvx)
    else
        iflag = 1
    end if
    me%iflag = iflag

    end subroutine evaluate_1d
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

    pure function bspline_2d_constructor_auto_knots(x,y,fcn,kx,ky) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer,intent(in)                 :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                 :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)

    call initialize_2d_auto_knots(me,x,y,fcn,kx,ky,me%iflag)

    end function bspline_2d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_2d]] type (user-specified knots).
!  This is a wrapper for [[initialize_2d_specify_knots]].

    pure function bspline_2d_constructor_specify_knots(x,y,fcn,kx,ky,tx,ty) result(me)

    implicit none

    type(bspline_2d)                   :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer,intent(in)                 :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                 :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.

    call initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,me%iflag)

    end function bspline_2d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with automatically-computed knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_auto_knots(me,x,y,fcn,kx,ky,iflag)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer,intent(in)                 :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                 :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(out)                :: iflag !! status flag (see [[db2ink]])

    integer :: iknot
    integer :: nx,ny

    call me%destroy()

    nx = size(x)
    ny = size(y)

    me%nx = nx
    me%ny = ny

    me%kx = kx
    me%ky = ky

    allocate(me%tx(nx+kx))
    allocate(me%ty(ny+ky))
    allocate(me%bcoef(nx,ny))

    iknot = 0         !knot sequence chosen by db2ink

    call db2ink(x,nx,y,ny,fcn,kx,ky,iknot,me%tx,me%ty,me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_2d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_2d]] type (with user-specified knots).
!  This is a wrapper for [[db2ink]].

    pure subroutine initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,iflag)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x     !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)   :: y     !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:),intent(in) :: fcn   !! `(nx,ny)` matrix of function values to interpolate.
                                                !! `fcn(i,j)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`)
    integer,intent(in)                 :: kx    !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                 :: ky    !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)   :: tx    !! The `(nx+kx)` knots in the \(x\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    real(wp),dimension(:),intent(in)   :: ty    !! The `(ny+ky)` knots in the \(y\) direction
                                                !! for the spline interpolant.
                                                !! Must be non-decreasing.
    integer,intent(out)                :: iflag !! status flag (see [[db2ink]])

    integer :: nx,ny

    call me%destroy()

    nx = size(x)
    ny = size(y)

    call check_knot_vectors_sizes('initialize_2d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                iflag=iflag)

    if (iflag == 0) then

        me%nx = nx
        me%ny = ny

        me%kx = kx
        me%ky = ky

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%bcoef(nx,ny))

        me%tx = tx
        me%ty = ty

        call db2ink(x,nx,y,ny,fcn,kx,ky,1,me%tx,me%ty,me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db2val]])

    if (me%initialized) then
        call db2val(xval,yval,&
                    idx,idy,&
                    me%tx,me%ty,&
                    me%nx,me%ny,&
                    me%kx,me%ky,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%iloy)
    else
        iflag = 1
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

    pure function bspline_3d_constructor_auto_knots(x,y,z,fcn,kx,ky,kz) result(me)

    implicit none

    type(bspline_3d)                     :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer,intent(in)                   :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)

    call initialize_3d_auto_knots(me,x,y,z,fcn,kx,ky,kz,me%iflag)

    end function bspline_3d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_3d]] type (user-specified knots).
!  This is a wrapper for [[initialize_3d_specify_knots]].

    pure function bspline_3d_constructor_specify_knots(x,y,z,fcn,kx,ky,kz,tx,ty,tz) result(me)

    implicit none

    type(bspline_3d)                     :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer,intent(in)                   :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: kz  !! The order of spline pieces in \(z\)
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

    call initialize_3d_specify_knots(me,x,y,z,fcn,kx,ky,kz,tx,ty,tz,me%iflag)

    end function bspline_3d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_3d]] type (with automatically-computed knots).
!  This is a wrapper for [[db3ink]].

    pure subroutine initialize_3d_auto_knots(me,x,y,z,fcn,kx,ky,kz,iflag)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer,intent(in)                   :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: kz  !! The order of spline pieces in \(z\)
                                                !! ( \( 2 \le k_z < n_z \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(out)                  :: iflag  !! status flag (see [[db3ink]])

    integer :: iknot
    integer :: nx,ny,nz

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)

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

    iknot = 0         !knot sequence chosen by db3ink

    call db3ink(x,nx,y,ny,z,nz,&
                fcn,&
                kx,ky,kz,&
                iknot,&
                me%tx,me%ty,me%tz,&
                me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_3d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_3d]] type (with user-specified knots).
!  This is a wrapper for [[db3ink]].

    pure subroutine initialize_3d_specify_knots(me,x,y,z,fcn,kx,ky,kz,tx,ty,tz,iflag)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)     :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in) :: fcn !! `(nx,ny,nz)` matrix of function values to interpolate.
                                                !! `fcn(i,j,k)` should contain the function value at the
                                                !! point (`x(i)`,`y(j)`,`z(k)`)
    integer,intent(in)                   :: kx  !! The order of spline pieces in \(x\)
                                                !! ( \( 2 \le k_x < n_x \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: ky  !! The order of spline pieces in \(y\)
                                                !! ( \( 2 \le k_y < n_y \) )
                                                !! (order = polynomial degree + 1)
    integer,intent(in)                   :: kz  !! The order of spline pieces in \(z\)
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
    integer,intent(out)                  :: iflag  !! status flag (see [[db3ink]])

    integer :: nx,ny,nz

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)

    call check_knot_vectors_sizes('initialize_3d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                nz=nz,kz=kz,tz=tz,&
                                                                iflag=iflag)

    if (iflag == 0) then

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

        me%tx = tx
        me%ty = ty
        me%tz = tz

        call db3ink(x,nx,y,ny,z,nz,&
                    fcn,&
                    kx,ky,kz,&
                    1,&
                    me%tx,me%ty,me%tz,&
                    me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db3val]])

    if (me%initialized) then
        call db3val(xval,yval,zval,&
                    idx,idy,idz,&
                    me%tx,me%ty,me%tz,&
                    me%nx,me%ny,me%nz,&
                    me%kx,me%ky,me%kz,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,&
                    me%iloy,me%iloz)
    else
        iflag = 1
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

    pure function bspline_4d_constructor_auto_knots(x,y,z,q,fcn,kx,ky,kz,kq) result(me)

    implicit none

    type(bspline_4d)                       :: me
    real(wp),dimension(:),intent(in)       :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)       :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in) :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                  !! `fcn(i,j,k,l)` should contain the function value at the
                                                  !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer,intent(in)                     :: kx  !! The order of spline pieces in \(x\)
                                                  !! ( \( 2 \le k_x < n_x \) )
                                                  !! (order = polynomial degree + 1)
    integer,intent(in)                     :: ky  !! The order of spline pieces in \(y\)
                                                  !! ( \( 2 \le k_y < n_y \) )
                                                  !! (order = polynomial degree + 1)
    integer,intent(in)                     :: kz  !! The order of spline pieces in \(z\)
                                                  !! ( \( 2 \le k_z < n_z \) )
                                                  !! (order = polynomial degree + 1)
    integer,intent(in)                     :: kq  !! The order of spline pieces in \(q\)
                                                  !! ( \( 2 \le k_q < n_q \) )
                                                  !! (order = polynomial degree + 1)

    call initialize_4d_auto_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,me%iflag)

    end function bspline_4d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_4d]] type (user-specified knots).
!  This is a wrapper for [[initialize_4d_specify_knots]].

    pure function bspline_4d_constructor_specify_knots(x,y,z,q,fcn,kx,ky,kz,kq,tx,ty,tz,tq) result(me)

    implicit none

    type(bspline_4d)                           :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
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

    call initialize_4d_specify_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,tx,ty,tz,tq,me%iflag)

    end function bspline_4d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_4d]] type (with automatically-computed knots).
!  This is a wrapper for [[db4ink]].

    pure subroutine initialize_4d_auto_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,iflag)

    implicit none

    class(bspline_4d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(out)                        :: iflag  !! status flag (see [[db4ink]])

    integer :: iknot
    integer :: nx,ny,nz,nq

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)

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

    iknot = 0         !knot sequence chosen by db4ink

    call db4ink(x,nx,y,ny,z,nz,q,nq,&
                fcn,&
                kx,ky,kz,kq,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,&
                me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_4d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_4d]] type (with user-specified knots).
!  This is a wrapper for [[db4ink]].

    pure subroutine initialize_4d_specify_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,tx,ty,tz,tq,iflag)

    implicit none

    class(bspline_4d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x   !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: y   !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: z   !! `(nz)` array of \(z\) abcissae. Must be strictly increasing.
    real(wp),dimension(:),intent(in)           :: q   !! `(nq)` array of \(q\) abcissae. Must be strictly increasing.
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn !! `(nx,ny,nz,nq)` matrix of function values to interpolate.
                                                      !! `fcn(i,j,k,l)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`,`z(k)`,`q(l)`)
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
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
    integer,intent(out)                    :: iflag   !! status flag (see [[db4ink]])

    integer :: nx,ny,nz,nq

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)

    call check_knot_vectors_sizes('initialize_4d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                nz=nz,kz=kz,tz=tz,&
                                                                nq=nq,kq=kq,tq=tq,&
                                                                iflag=iflag)

    if (iflag == 0) then

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

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq

        call db4ink(x,nx,y,ny,z,nz,q,nq,&
                    fcn,&
                    kx,ky,kz,kq,&
                    1,&
                    me%tx,me%ty,me%tz,me%tq,&
                    me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db4val]])

    if (me%initialized) then
        call db4val(xval,yval,zval,qval,&
                    idx,idy,idz,idq,&
                    me%tx,me%ty,me%tz,me%tq,&
                    me%nx,me%ny,me%nz,me%nq,&
                    me%kx,me%ky,me%kz,me%kq,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,&
                    me%iloy,me%iloz,me%iloq)
    else
        iflag = 1
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

    pure function bspline_5d_constructor_auto_knots(x,y,z,q,r,fcn,kx,ky,kz,kq,kr) result(me)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)

    call initialize_5d_auto_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,me%iflag)

    end function bspline_5d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_5d]] type (user-specified knots).
!  This is a wrapper for [[initialize_5d_specify_knots]].

    pure function bspline_5d_constructor_specify_knots(x,y,z,q,r,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr) result(me)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
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

    call initialize_5d_specify_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,me%iflag)

    end function bspline_5d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_5d]] type (with automatically-computed knots).
!  This is a wrapper for [[db5ink]].

    pure subroutine initialize_5d_auto_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,iflag)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(out)                        :: iflag !! status flag (see [[db5ink]])

    integer :: iknot
    integer :: nx,ny,nz,nq,nr

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)
    nr = size(r)

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

    iknot = 0         !knot sequence chosen by db5ink

    call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                fcn,&
                kx,ky,kz,kq,kr,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,me%tr,&
                me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_5d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_5d]] type (with user-specified knots).
!  This is a wrapper for [[db5ink]].

    pure subroutine initialize_5d_specify_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,iflag)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
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
    integer,intent(out)                      :: iflag !! status flag (see [[db5ink]])

    integer :: nx,ny,nz,nq,nr

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)
    nr = size(r)

    call check_knot_vectors_sizes('initialize_5d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                nz=nz,kz=kz,tz=tz,&
                                                                nq=nq,kq=kq,tq=tq,&
                                                                nr=nr,kr=kr,tr=tr,&
                                                                iflag=iflag)

    if (iflag == 0) then

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

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq
        me%tr = tr

        call db5ink(x,nx,y,ny,z,nz,q,nq,r,nr,&
                    fcn,&
                    kx,ky,kz,kq,kr,&
                    1,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,&
                    me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idr   !! \(r\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db5val]])

    if (me%initialized) then
        call db5val(xval,yval,zval,qval,rval,&
                    idx,idy,idz,idq,idr,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,&
                    me%nx,me%ny,me%nz,me%nq,me%nr,&
                    me%kx,me%ky,me%kz,me%kq,me%kr,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,&
                    me%iloy,me%iloz,me%iloq,me%ilor)
    else
        iflag = 1
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

    pure function bspline_6d_constructor_auto_knots(x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks) result(me)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)

    call initialize_6d_auto_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,me%iflag)

    end function bspline_6d_constructor_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[bspline_6d]] type (user-specified knots).
!  This is a wrapper for [[initialize_6d_specify_knots]].

    pure function bspline_6d_constructor_specify_knots(x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts) result(me)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ks  !! The order of spline pieces in \(z\)
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

    call initialize_6d_specify_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,me%iflag)

    end function bspline_6d_constructor_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_6d]] type (with automatically-computed knots).
!  This is a wrapper for [[db6ink]].

    pure subroutine initialize_6d_auto_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,iflag)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ks  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(out)                        :: iflag !! status flag (see [[db6ink]])

    integer :: iknot
    integer :: nx,ny,nz,nq,nr,ns

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)
    nr = size(r)
    ns = size(s)

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

    iknot = 0         !knot sequence chosen by db6ink

    call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                fcn,&
                kx,ky,kz,kq,kr,ks,&
                iknot,&
                me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                me%bcoef,iflag)

    me%initialized = iflag==0
    me%iflag = iflag

    end subroutine initialize_6d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a [[bspline_6d]] type (with user-specified knots).
!  This is a wrapper for [[db6ink]].

    pure subroutine initialize_6d_specify_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,iflag)

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
    integer,intent(in)                         :: kx  !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ky  !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kz  !! The order of spline pieces in \(z\)
                                                      !! ( \( 2 \le k_z < n_z \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kq  !! The order of spline pieces in \(q\)
                                                      !! ( \( 2 \le k_q < n_q \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: kr  !! The order of spline pieces in \(r\)
                                                      !! ( \( 2 \le k_r < n_r \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                         :: ks  !! The order of spline pieces in \(z\)
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
    integer,intent(out)                        :: iflag !! status flag (see [[db6ink]])

    integer :: nx,ny,nz,nq,nr,ns

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)
    nq = size(q)
    nr = size(r)
    ns = size(s)

    call check_knot_vectors_sizes('initialize_6d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                nz=nz,kz=kz,tz=tz,&
                                                                nq=nq,kq=kq,tq=tq,&
                                                                nr=nr,kr=kr,tr=tr,&
                                                                ns=ns,ks=ks,ts=ts,&
                                                                iflag=iflag)

    if (iflag == 0) then

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

        me%tx = tx
        me%ty = ty
        me%tz = tz
        me%tq = tq
        me%tr = tr
        me%ts = ts

        call db6ink(x,nx,y,ny,z,nz,q,nq,r,nr,s,ns,&
                    fcn,&
                    kx,ky,kz,kq,kr,ks,&
                    1,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                    me%bcoef,iflag)

    end if

    me%initialized = iflag==0
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
    integer,intent(in)              :: idx   !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idy   !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idz   !! \(z\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idq   !! \(q\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: idr   !! \(r\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)              :: ids   !! \(s\) derivative of piecewise polynomial to evaluate.
    real(wp),intent(out)            :: f     !! interpolated value
    integer,intent(out)             :: iflag !! status flag (see [[db6val]])

    if (me%initialized) then
        call db6val(xval,yval,zval,qval,rval,sval,&
                    idx,idy,idz,idq,idr,ids,&
                    me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                    me%nx,me%ny,me%nz,me%nq,me%nr,me%ns,&
                    me%kx,me%ky,me%kz,me%kq,me%kr,me%ks,&
                    me%bcoef,f,iflag,&
                    me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,me%inbvs,&
                    me%iloy,me%iloz,me%iloq,me%ilor,me%ilos)
    else
        iflag = 1
    end if

    me%iflag = iflag

    end subroutine evaluate_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Error checks for the user-specified knot vector sizes.
!  Note that if more than one is the wrong size, then the iflag error code will
!  correspond to the one with for the highest rank.

    pure subroutine check_knot_vectors_sizes(routine,nx,ny,nz,nq,nr,ns,&
                                             kx,ky,kz,kq,kr,ks,&
                                             tx,ty,tz,tq,tr,ts,iflag)

    implicit none

    character(len=*),intent(in)               :: routine
    integer,intent(in),optional               :: nx
    integer,intent(in),optional               :: ny
    integer,intent(in),optional               :: nz
    integer,intent(in),optional               :: nq
    integer,intent(in),optional               :: nr
    integer,intent(in),optional               :: ns
    integer,intent(in),optional               :: kx
    integer,intent(in),optional               :: ky
    integer,intent(in),optional               :: kz
    integer,intent(in),optional               :: kq
    integer,intent(in),optional               :: kr
    integer,intent(in),optional               :: ks
    real(wp),dimension(:),intent(in),optional :: tx
    real(wp),dimension(:),intent(in),optional :: ty
    real(wp),dimension(:),intent(in),optional :: tz
    real(wp),dimension(:),intent(in),optional :: tq
    real(wp),dimension(:),intent(in),optional :: tr
    real(wp),dimension(:),intent(in),optional :: ts
    integer,intent(out)                       :: iflag  !! 0 if everything is OK

    iflag = 0

    if (present(nx) .and. present(kx) .and. present(tx)) then
        if (size(tx)/=(nx+kx)) then
            !write(error_unit,'(A)') trim(routine)//' - tx is not the correct size (nx+kx)'
            iflag = 501
        end if
    end if

    if (present(ny) .and. present(ky) .and. present(ty)) then
        if (size(ty)/=(ny+ky)) then
            !write(error_unit,'(A)') trim(routine)//' - ty is not the correct size (ny+ky)'
            iflag = 502
        end if
    end if

    if (present(nz) .and. present(kz) .and. present(tz)) then
        if (size(tz)/=(nz+kz)) then
            !write(error_unit,'(A)') trim(routine)//' - tz is not the correct size (nz+kz)'
            iflag = 503
        end if
    end if

    if (present(nq) .and. present(kq) .and. present(tq)) then
        if (size(tq)/=(nq+kq)) then
            !write(error_unit,'(A)') trim(routine)//' - tq is not the correct size (nq+kq)'
            iflag = 504
        end if
    end if

    if (present(nr) .and. present(kr) .and. present(tr)) then
        if (size(tr)/=(nr+kr)) then
            !write(error_unit,'(A)') trim(routine)//' - tr is not the correct size (nr+kr)'
            iflag = 505
        end if
    end if

    if (present(ns) .and. present(ks) .and. present(ts)) then
        if (size(ts)/=(ns+ks)) then
            !write(error_unit,'(A)') trim(routine)//' - ts is not the correct size (ns+ks)'
            iflag = 506
        end if
    end if

    end subroutine check_knot_vectors_sizes
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_oo_module
!*****************************************************************************************
