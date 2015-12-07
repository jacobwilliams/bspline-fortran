!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!  date: 12/6/2015
!
!# Description
!
!  Object-oriented style wrappers to [[bspline_sub_module]].
!  This module provides classes ([[bspline_1d]], [[bspline_2d]],
!  [[bspline_3d]], [[bspline_4d]], [[bspline_5d]], and [[bspline_6d]])
!  which can be used instead of the main subroutine interface.

    module bspline_oo_module

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use,intrinsic :: iso_fortran_env, only: error_unit
    use bspline_sub_module

    implicit none

    private

    type,public,abstract :: bspline_class
        !! Base class for the b-spline types
        private
        integer :: inbvx = 1  !! internal variable used by dbvalu for efficient processing
    contains
        private
        procedure(destroy_func),deferred,public :: destroy  !! destructor
    end type bspline_class

    abstract interface
        subroutine destroy_func(me)  !! interface for bspline destructor routines
        import :: bspline_class
        implicit none
        class(bspline_class),intent(out) :: me
        end subroutine destroy_func
    end interface

    type,extends(bspline_class),public :: bspline_1d
        !! Class for 1d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: kx  = 0
        real(wp),dimension(:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        contains
        private
        generic,public :: initialize => initialize_1d_auto_knots,initialize_1d_specify_knots
        procedure :: initialize_1d_auto_knots,initialize_1d_specify_knots
        procedure,public :: evaluate => evaluate_1d
        procedure,public :: destroy => destroy_1d
    end type bspline_1d

    type,extends(bspline_class),public :: bspline_2d
        !! Class for 2d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        real(wp),dimension(:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        real(wp),dimension(:),allocatable :: ty
        integer :: inbvy = 1
        integer :: iloy = 1
        contains
        private
        generic,public :: initialize => initialize_2d_auto_knots,initialize_2d_specify_knots
        procedure :: initialize_2d_auto_knots,initialize_2d_specify_knots
        procedure,public :: evaluate => evaluate_2d
        procedure,public :: destroy => destroy_2d
    end type bspline_2d

    type,extends(bspline_class),public :: bspline_3d
        !! Class for 3d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: nz  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        integer :: kz  = 0
        real(wp),dimension(:,:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        real(wp),dimension(:),allocatable :: ty
        real(wp),dimension(:),allocatable :: tz
        integer :: inbvy = 1
        integer :: inbvz = 1
        integer :: iloy = 1
        integer :: iloz = 1
        contains
        private
        generic,public :: initialize => initialize_3d_auto_knots,initialize_3d_specify_knots
        procedure :: initialize_3d_auto_knots,initialize_3d_specify_knots
        procedure,public :: evaluate => evaluate_3d
        procedure,public :: destroy => destroy_3d
    end type bspline_3d

    type,extends(bspline_class),public :: bspline_4d
        !! Class for 4d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: nz  = 0
        integer :: nq  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        integer :: kz  = 0
        integer :: kq  = 0
        integer,dimension(4) :: k  = 0
        real(wp),dimension(:,:,:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        real(wp),dimension(:),allocatable :: ty
        real(wp),dimension(:),allocatable :: tz
        real(wp),dimension(:),allocatable :: tq
        integer :: inbvy = 1
        integer :: inbvz = 1
        integer :: inbvq = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        contains
        private
        generic,public :: initialize => initialize_4d_auto_knots,initialize_4d_specify_knots
        procedure :: initialize_4d_auto_knots,initialize_4d_specify_knots
        procedure,public :: evaluate => evaluate_4d
        procedure,public :: destroy => destroy_4d
    end type bspline_4d

    type,extends(bspline_class),public :: bspline_5d
        !! Class for 5d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: nz  = 0
        integer :: nq  = 0
        integer :: nr  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        integer :: kz  = 0
        integer :: kq  = 0
        integer :: kr  = 0
        real(wp),dimension(:,:,:,:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        real(wp),dimension(:),allocatable :: ty
        real(wp),dimension(:),allocatable :: tz
        real(wp),dimension(:),allocatable :: tq
        real(wp),dimension(:),allocatable :: tr
        integer :: inbvy = 1
        integer :: inbvz = 1
        integer :: inbvq = 1
        integer :: inbvr = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        contains
        private
        generic,public :: initialize => initialize_5d_auto_knots,initialize_5d_specify_knots
        procedure :: initialize_5d_auto_knots,initialize_5d_specify_knots
        procedure,public :: evaluate => evaluate_5d
        procedure,public :: destroy => destroy_5d
    end type bspline_5d

    type,extends(bspline_class),public :: bspline_6d
        !! Class for 6d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: nz  = 0
        integer :: nq  = 0
        integer :: nr  = 0
        integer :: ns  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        integer :: kz  = 0
        integer :: kq  = 0
        integer :: kr  = 0
        integer :: ks  = 0
        real(wp),dimension(:,:,:,:,:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx
        real(wp),dimension(:),allocatable :: ty
        real(wp),dimension(:),allocatable :: tz
        real(wp),dimension(:),allocatable :: tq
        real(wp),dimension(:),allocatable :: tr
        real(wp),dimension(:),allocatable :: ts
        integer :: inbvy = 1
        integer :: inbvz = 1
        integer :: inbvq = 1
        integer :: inbvr = 1
        integer :: inbvs = 1
        integer :: iloy  = 1
        integer :: iloz  = 1
        integer :: iloq  = 1
        integer :: ilor  = 1
        integer :: ilos  = 1
        contains
        private
        generic,public :: initialize => initialize_6d_auto_knots,initialize_6d_specify_knots
        procedure :: initialize_6d_auto_knots,initialize_6d_specify_knots
        procedure,public :: evaluate => evaluate_6d
        procedure,public :: destroy => destroy_6d
    end type bspline_6d

    contains
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_1d]] type.

    subroutine destroy_1d(me)

    implicit none

    class(bspline_1d),intent(out) :: me

    end subroutine destroy_1d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_1d]] type (with automatically-computed knots).
!  This is a wrapper for [[db1ink]].

    subroutine initialize_1d_auto_knots(me,x,fcn,kx,iflag)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(in) :: fcn
    integer,intent(in)               :: kx
    integer,intent(out)              :: iflag

    integer :: iknot
    integer :: nx

    call me%destroy()

    nx = size(x)

    me%nx = nx
    me%kx = kx

    allocate(me%tx(nx+kx))
    allocate(me%bcoef(nx))

    iknot = 0         !knot sequence chosen by db2ink

    call db1ink(x,nx,fcn,kx,me%tx,me%bcoef,iknot)

    iflag = iknot     !status flag

end subroutine initialize_1d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_1d]] type (with user-specified knots).
!  This is a wrapper for [[db1ink]].

    subroutine initialize_1d_specify_knots(me,x,fcn,kx,tx,iflag)

    implicit none

    class(bspline_1d),intent(inout)  :: me
    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(in) :: fcn
    integer,intent(in)               :: kx
    real(wp),dimension(:),intent(in) :: tx
    integer,intent(out)              :: iflag

    integer :: nx

    call me%destroy()

    nx = size(x)

    call check_knot_vectors_sizes('initialize_1d_specify_knots',nx=nx,kx=kx,tx=tx,iflag=iflag)

    if (iflag == 1) then

        me%nx = nx
        me%kx = kx

        allocate(me%tx(nx+kx))
        allocate(me%bcoef(nx))

        me%tx = tx

        call db1ink(x,nx,fcn,kx,me%tx,me%bcoef,iflag)

    end if

    end subroutine initialize_1d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_1d]] interpolate.  This is a wrapper for [[db1val]].

    subroutine evaluate_1d(me,xval,idx,f,iflag)

    implicit none

    class(bspline_1d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    integer,intent(in)              :: idx
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db1val(xval,idx,me%tx,me%nx,me%kx,me%bcoef,f,iflag,me%inbvx)

    end subroutine evaluate_1d
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_2d]] type.

    subroutine destroy_2d(me)

    implicit none

    class(bspline_2d),intent(out) :: me

    end subroutine destroy_2d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_2d]] type (with automatically-computed knots).
!  This is a wrapper for [[db2ink]].

    subroutine initialize_2d_auto_knots(me,x,y,fcn,kx,ky,iflag)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x
    real(wp),dimension(:),intent(in)   :: y
    real(wp),dimension(:,:),intent(in) :: fcn
    integer,intent(in)                 :: kx
    integer,intent(in)                 :: ky
    integer,intent(out)                :: iflag

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

    call db2ink(x,nx,y,ny,fcn,kx,ky,me%tx,me%ty,me%bcoef,iknot)

    iflag = iknot     !status flag

    end subroutine initialize_2d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_2d]] type (with user-specified knots).
!  This is a wrapper for [[db2ink]].

    subroutine initialize_2d_specify_knots(me,x,y,fcn,kx,ky,tx,ty,iflag)

    implicit none

    class(bspline_2d),intent(inout)    :: me
    real(wp),dimension(:),intent(in)   :: x
    real(wp),dimension(:),intent(in)   :: y
    real(wp),dimension(:,:),intent(in) :: fcn
    integer,intent(in)                 :: kx
    integer,intent(in)                 :: ky
    real(wp),dimension(:),intent(in)   :: tx
    real(wp),dimension(:),intent(in)   :: ty
    integer,intent(out)                :: iflag

    integer :: nx,ny

    call me%destroy()

    nx = size(x)
    ny = size(y)

    call check_knot_vectors_sizes('initialize_2d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                iflag=iflag)

    if (iflag == 1) then

        me%nx = nx
        me%ny = ny

        me%kx = kx
        me%ky = ky

        allocate(me%tx(nx+kx))
        allocate(me%ty(ny+ky))
        allocate(me%bcoef(nx,ny))

        me%tx = tx
        me%ty = ty

        call db2ink(x,nx,y,ny,fcn,kx,ky,me%tx,me%ty,me%bcoef,iflag)

    end if

    end subroutine initialize_2d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_2d]] interpolate.  This is a wrapper for [[db2val]].

    subroutine evaluate_2d(me,xval,yval,idx,idy,f,iflag)

    implicit none

    class(bspline_2d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    real(wp),intent(in)             :: yval
    integer,intent(in)              :: idx
    integer,intent(in)              :: idy
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db2val(xval,yval,&
                idx,idy,&
                me%tx,me%ty,&
                me%nx,me%ny,&
                me%kx,me%ky,&
                me%bcoef,f,iflag,&
                me%inbvx,me%inbvy,me%iloy)

    end subroutine evaluate_2d
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_3d]] type.

    subroutine destroy_3d(me)

    implicit none

    class(bspline_3d),intent(out) :: me

    end subroutine destroy_3d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_3d]] type (with automatically-computed knots).
!  This is a wrapper for [[db3ink]].

    subroutine initialize_3d_auto_knots(me,x,y,z,fcn,kx,ky,kz,iflag)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(in)     :: y
    real(wp),dimension(:),intent(in)     :: z
    real(wp),dimension(:,:,:),intent(in) :: fcn
    integer,intent(in)                   :: kx
    integer,intent(in)                   :: ky
    integer,intent(in)                   :: kz
    integer,intent(out)                  :: iflag

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
                me%tx,me%ty,me%tz,&
                me%bcoef,iknot)

    iflag = iknot     !status flag

    end subroutine initialize_3d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_3d]] type (with user-specified knots).
!  This is a wrapper for [[db3ink]].

    subroutine initialize_3d_specify_knots(me,x,y,z,fcn,kx,ky,kz,tx,ty,tz,iflag)

    implicit none

    class(bspline_3d),intent(inout)      :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(in)     :: y
    real(wp),dimension(:),intent(in)     :: z
    real(wp),dimension(:,:,:),intent(in) :: fcn
    integer,intent(in)                   :: kx
    integer,intent(in)                   :: ky
    integer,intent(in)                   :: kz
    real(wp),dimension(:),intent(in)     :: tx
    real(wp),dimension(:),intent(in)     :: ty
    real(wp),dimension(:),intent(in)     :: tz
    integer,intent(out)                  :: iflag

    integer :: nx,ny,nz

    call me%destroy()

    nx = size(x)
    ny = size(y)
    nz = size(z)

    call check_knot_vectors_sizes('initialize_3d_specify_knots',nx=nx,kx=kx,tx=tx,&
                                                                ny=ny,ky=ky,ty=ty,&
                                                                nz=nz,kz=kz,tz=tz,&
                                                                iflag=iflag)

    if (iflag == 1) then

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
                    me%tx,me%ty,me%tz,&
                    me%bcoef,iflag)

    end if

    end subroutine initialize_3d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_3d]] interpolate.  This is a wrapper for [[db3val]].

    subroutine evaluate_3d(me,xval,yval,zval,idx,idy,idz,f,iflag)

    implicit none

    class(bspline_3d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    real(wp),intent(in)             :: yval
    real(wp),intent(in)             :: zval
    integer,intent(in)              :: idx
    integer,intent(in)              :: idy
    integer,intent(in)              :: idz
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db3val(xval,yval,zval,&
                idx,idy,idz,&
                me%tx,me%ty,me%tz,&
                me%nx,me%ny,me%nz,&
                me%kx,me%ky,me%kz,&
                me%bcoef,f,iflag,&
                me%inbvx,me%inbvy,me%inbvz,&
                me%iloy,me%iloz)

    end subroutine evaluate_3d
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_4d]] type.

    subroutine destroy_4d(me)

    implicit none

    class(bspline_4d),intent(out) :: me

    end subroutine destroy_4d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_4d]] type (with automatically-computed knots).
!  This is a wrapper for [[db4ink]].

    subroutine initialize_4d_auto_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,iflag)

    implicit none

    class(bspline_4d),intent(inout)        :: me
    real(wp),dimension(:),intent(in)       :: x
    real(wp),dimension(:),intent(in)       :: y
    real(wp),dimension(:),intent(in)       :: z
    real(wp),dimension(:),intent(in)       :: q
    real(wp),dimension(:,:,:,:),intent(in) :: fcn
    integer,intent(in)                     :: kx
    integer,intent(in)                     :: ky
    integer,intent(in)                     :: kz
    integer,intent(in)                     :: kq
    integer,intent(out)                    :: iflag

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
                me%tx,me%ty,me%tz,me%tq,&
                me%bcoef,iknot)

    iflag = iknot     !status flag

    end subroutine initialize_4d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_4d]] type (with user-specified knots).
!  This is a wrapper for [[db4ink]].

    subroutine initialize_4d_specify_knots(me,x,y,z,q,fcn,kx,ky,kz,kq,tx,ty,tz,tq,iflag)

    implicit none

    class(bspline_4d),intent(inout)        :: me
    real(wp),dimension(:),intent(in)       :: x
    real(wp),dimension(:),intent(in)       :: y
    real(wp),dimension(:),intent(in)       :: z
    real(wp),dimension(:),intent(in)       :: q
    real(wp),dimension(:,:,:,:),intent(in) :: fcn
    integer,intent(in)                     :: kx
    integer,intent(in)                     :: ky
    integer,intent(in)                     :: kz
    integer,intent(in)                     :: kq
    real(wp),dimension(:),intent(in)       :: tx
    real(wp),dimension(:),intent(in)       :: ty
    real(wp),dimension(:),intent(in)       :: tz
    real(wp),dimension(:),intent(in)       :: tq
    integer,intent(out)                    :: iflag

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

    if (iflag == 1) then

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
                    me%tx,me%ty,me%tz,me%tq,&
                    me%bcoef,iflag)

    end if

    end subroutine initialize_4d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_4d]] interpolate.  This is a wrapper for [[db4val]].

    subroutine evaluate_4d(me,xval,yval,zval,qval,idx,idy,idz,idq,f,iflag)

    implicit none

    class(bspline_4d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    real(wp),intent(in)             :: yval
    real(wp),intent(in)             :: zval
    real(wp),intent(in)             :: qval
    integer,intent(in)              :: idx
    integer,intent(in)              :: idy
    integer,intent(in)              :: idz
    integer,intent(in)              :: idq
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db4val(xval,yval,zval,qval,&
                idx,idy,idz,idq,&
                me%tx,me%ty,me%tz,me%tq,&
                me%nx,me%ny,me%nz,me%nq,&
                me%kx,me%ky,me%kz,me%kq,&
                me%bcoef,f,iflag,&
                me%inbvx,me%inbvy,me%inbvz,me%inbvq,&
                me%iloy,me%iloz,me%iloq)

    end subroutine evaluate_4d
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_5d]] type.

    subroutine destroy_5d(me)

    implicit none

    class(bspline_5d),intent(out) :: me

    end subroutine destroy_5d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_5d]] type (with automatically-computed knots).
!  This is a wrapper for [[db5ink]].

    subroutine initialize_5d_auto_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,iflag)

    implicit none

    class(bspline_5d),intent(inout)          :: me
    real(wp),dimension(:),intent(in)         :: x
    real(wp),dimension(:),intent(in)         :: y
    real(wp),dimension(:),intent(in)         :: z
    real(wp),dimension(:),intent(in)         :: q
    real(wp),dimension(:),intent(in)         :: r
    real(wp),dimension(:,:,:,:,:),intent(in) :: fcn
    integer,intent(in)                       :: kx
    integer,intent(in)                       :: ky
    integer,intent(in)                       :: kz
    integer,intent(in)                       :: kq
    integer,intent(in)                       :: kr
    integer,intent(out)                      :: iflag

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
                me%tx,me%ty,me%tz,me%tq,me%tr,&
                me%bcoef,iknot)

    iflag = iknot     !status flag

    end subroutine initialize_5d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_5d]] type (with user-specified knots).
!  This is a wrapper for [[db5ink]].

    subroutine initialize_5d_specify_knots(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,tx,ty,tz,tq,tr,iflag)

    implicit none

    class(bspline_5d),intent(inout)          :: me
    real(wp),dimension(:),intent(in)         :: x
    real(wp),dimension(:),intent(in)         :: y
    real(wp),dimension(:),intent(in)         :: z
    real(wp),dimension(:),intent(in)         :: q
    real(wp),dimension(:),intent(in)         :: r
    real(wp),dimension(:,:,:,:,:),intent(in) :: fcn
    integer,intent(in)                       :: kx
    integer,intent(in)                       :: ky
    integer,intent(in)                       :: kz
    integer,intent(in)                       :: kq
    integer,intent(in)                       :: kr
    real(wp),dimension(:),intent(in)         :: tx
    real(wp),dimension(:),intent(in)         :: ty
    real(wp),dimension(:),intent(in)         :: tz
    real(wp),dimension(:),intent(in)         :: tq
    real(wp),dimension(:),intent(in)         :: tr
    integer,intent(out)                      :: iflag

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

    if (iflag == 1) then

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
                    me%tx,me%ty,me%tz,me%tq,me%tr,&
                    me%bcoef,iflag)

    end if

    end subroutine initialize_5d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_5d]] interpolate.  This is a wrapper for [[db5val]].

    subroutine evaluate_5d(me,xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,f,iflag)

    implicit none

    class(bspline_5d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    real(wp),intent(in)             :: yval
    real(wp),intent(in)             :: zval
    real(wp),intent(in)             :: qval
    real(wp),intent(in)             :: rval
    integer,intent(in)              :: idx
    integer,intent(in)              :: idy
    integer,intent(in)              :: idz
    integer,intent(in)              :: idq
    integer,intent(in)              :: idr
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db5val(xval,yval,zval,qval,rval,&
                idx,idy,idz,idq,idr,&
                me%tx,me%ty,me%tz,me%tq,me%tr,&
                me%nx,me%ny,me%nz,me%nq,me%nr,&
                me%kx,me%ky,me%kz,me%kq,me%kr,&
                me%bcoef,f,iflag,&
                me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,&
                me%iloy,me%iloz,me%iloq,me%ilor)

    end subroutine evaluate_5d
!*****************************************************************************************

!*****************************************************************************************
!> Destructor for [[bspline_6d]] type.

    subroutine destroy_6d(me)

    implicit none

    class(bspline_6d),intent(out) :: me

    end subroutine destroy_6d
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_6d]] type (with automatically-computed knots).
!  This is a wrapper for [[db6ink]].

    subroutine initialize_6d_auto_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,iflag)

    implicit none

    class(bspline_6d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x
    real(wp),dimension(:),intent(in)           :: y
    real(wp),dimension(:),intent(in)           :: z
    real(wp),dimension(:),intent(in)           :: q
    real(wp),dimension(:),intent(in)           :: r
    real(wp),dimension(:),intent(in)           :: s
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn
    integer,intent(in)                         :: kx
    integer,intent(in)                         :: ky
    integer,intent(in)                         :: kz
    integer,intent(in)                         :: kq
    integer,intent(in)                         :: kr
    integer,intent(in)                         :: ks
    integer,intent(out)                        :: iflag

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
                me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                me%bcoef,iknot)

    iflag = iknot     !status flag

    end subroutine initialize_6d_auto_knots
!*****************************************************************************************

!*****************************************************************************************
!> Initialize a [[bspline_6d]] type (with user-specified knots).
!  This is a wrapper for [[db6ink]].

    subroutine initialize_6d_specify_knots(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,iflag)

    implicit none

    class(bspline_6d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x
    real(wp),dimension(:),intent(in)           :: y
    real(wp),dimension(:),intent(in)           :: z
    real(wp),dimension(:),intent(in)           :: q
    real(wp),dimension(:),intent(in)           :: r
    real(wp),dimension(:),intent(in)           :: s
    real(wp),dimension(:,:,:,:,:,:),intent(in) :: fcn
    integer,intent(in)                         :: kx
    integer,intent(in)                         :: ky
    integer,intent(in)                         :: kz
    integer,intent(in)                         :: kq
    integer,intent(in)                         :: kr
    integer,intent(in)                         :: ks
    real(wp),dimension(:),intent(in)           :: tx
    real(wp),dimension(:),intent(in)           :: ty
    real(wp),dimension(:),intent(in)           :: tz
    real(wp),dimension(:),intent(in)           :: tq
    real(wp),dimension(:),intent(in)           :: tr
    real(wp),dimension(:),intent(in)           :: ts
    integer,intent(out)                        :: iflag

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

    if (iflag == 1) then

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
                    me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                    me%bcoef,iflag)

    end if

    end subroutine initialize_6d_specify_knots
!*****************************************************************************************

!*****************************************************************************************
!> Evaluate a [[bspline_6d]] interpolate.  This is a wrapper for [[db6val]].

    subroutine evaluate_6d(me,xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,f,iflag)

    implicit none

    class(bspline_6d),intent(inout) :: me
    real(wp),intent(in)             :: xval
    real(wp),intent(in)             :: yval
    real(wp),intent(in)             :: zval
    real(wp),intent(in)             :: qval
    real(wp),intent(in)             :: rval
    real(wp),intent(in)             :: sval
    integer,intent(in)              :: idx
    integer,intent(in)              :: idy
    integer,intent(in)              :: idz
    integer,intent(in)              :: idq
    integer,intent(in)              :: idr
    integer,intent(in)              :: ids
    real(wp),intent(out)            :: f
    integer,intent(out)             :: iflag

    call db6val(xval,yval,zval,qval,rval,sval,&
                idx,idy,idz,idq,idr,ids,&
                me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                me%nx,me%ny,me%nz,me%nq,me%nr,me%ns,&
                me%kx,me%ky,me%kz,me%kq,me%kr,me%ks,&
                me%bcoef,f,iflag,&
                me%inbvx,me%inbvy,me%inbvz,me%inbvq,me%inbvr,me%inbvs,&
                me%iloy,me%iloz,me%iloq,me%ilor,me%ilos)

    end subroutine evaluate_6d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Error checks for the user-specified knot vector sizes.
!  Note that if more than one is the wrong size, then the iflag error code will
!  correspond to the one with for the highest rank.

    subroutine check_knot_vectors_sizes(routine,nx,ny,nz,nq,nr,ns,kx,ky,kz,kq,kr,ks,tx,ty,tz,tq,tr,ts,iflag)

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
    integer,intent(out)                       :: iflag  !! 1 if everything is OK

    iflag = 1

    if (present(nx) .and. present(kx) .and. present(tx)) then
        if (size(tx)/=(nx+kx)) then
            write(error_unit,'(A)') trim(routine)//' - tx is not the correct size (nx+kx)'
            iflag = 501
        end if
    end if

    if (present(ny) .and. present(ky) .and. present(ty)) then
        if (size(ty)/=(ny+ky)) then
            write(error_unit,'(A)') trim(routine)//' - ty is not the correct size (ny+ky)'
            iflag = 502
        end if
    end if

    if (present(nz) .and. present(kz) .and. present(tz)) then
        if (size(tz)/=(nz+kz)) then
            write(error_unit,'(A)') trim(routine)//' - tz is not the correct size (nz+kz)'
            iflag = 503
        end if
    end if

    if (present(nq) .and. present(kq) .and. present(tq)) then
        if (size(tq)/=(nq+kq)) then
            write(error_unit,'(A)') trim(routine)//' - tq is not the correct size (nq+kq)'
            iflag = 504
        end if
    end if

    if (present(nr) .and. present(kr) .and. present(tr)) then
        if (size(tr)/=(nr+kr)) then
            write(error_unit,'(A)') trim(routine)//' - tr is not the correct size (nr+kr)'
            iflag = 505
        end if
    end if

    if (present(ns) .and. present(ks) .and. present(ts)) then
        if (size(ts)/=(ns+ks)) then
            write(error_unit,'(A)') trim(routine)//' - ts is not the correct size (ns+ks)'
            iflag = 506
        end if
    end if

    end subroutine check_knot_vectors_sizes
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_oo_module
!*****************************************************************************************
