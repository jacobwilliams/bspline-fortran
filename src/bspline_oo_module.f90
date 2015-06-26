!*****************************************************************************************
!
!> author: Jacob Williams
!  license: BSD
! 
!# Description
!
!  Object-oriented style wrappers to [[bspline_sub_module]].
!  This module provides classes ([[bspline_2d]], [[bspline_3d]], 
!  [[bspline_4d]], [[bspline_5d]], and [[bspline_6d]]) which can
!  be used instead of the main subroutine interface.
!
!# History
!
!  * Jacob Williams, created 6/20/2015.

    module bspline_oo_module
    
    use,intrinsic :: iso_fortran_env, only: wp => real64
    use bspline_sub_module
    
    implicit none
    
    private
    
    type,public :: bspline_2d
        !! Class for 2d b-spline interpolation.
        private
        integer :: nx  = 0
        integer :: ny  = 0
        integer :: kx  = 0
        integer :: ky  = 0
        real(wp),dimension(:,:),allocatable :: bcoef
        real(wp),dimension(:),allocatable :: tx    
        real(wp),dimension(:),allocatable :: ty
        contains
        procedure,public :: initialize => initialize_2d
        procedure,public :: evaluate => evaluate_2d
        procedure,public :: destroy => destroy_2d
    end type bspline_2d

    type,public :: bspline_3d
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
        contains
        procedure,public :: initialize => initialize_3d
        procedure,public :: evaluate => evaluate_3d
        procedure,public :: destroy => destroy_3d
    end type bspline_3d

    type,public :: bspline_4d
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
        contains
        procedure,public :: initialize => initialize_4d
        procedure,public :: evaluate => evaluate_4d
        procedure,public :: destroy => destroy_4d
    end type bspline_4d

    type,public :: bspline_5d
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
        contains
        procedure,public :: initialize => initialize_5d
        procedure,public :: evaluate => evaluate_5d
        procedure,public :: destroy => destroy_5d
    end type bspline_5d

    type,public :: bspline_6d
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
        contains
        procedure,public :: initialize => initialize_6d
        procedure,public :: evaluate => evaluate_6d
        procedure,public :: destroy => destroy_6d
    end type bspline_6d

    contains
!*****************************************************************************************
    
!*****************************************************************************************
!> Destructor for [[bspline_2d]] type.

    subroutine destroy_2d(me)
    
    implicit none
    
    class(bspline_2d),intent(out) :: me
    
    end subroutine destroy_2d
!*****************************************************************************************
    
!*****************************************************************************************
!> Initialize a [[bspline_2d]] type.  This is a wrapper for [[db2ink]].

    subroutine initialize_2d(me,x,y,fcn,kx,ky,iflag)

    implicit none
    
    class(bspline_2d),intent(inout)         :: me
    real(wp),dimension(:),intent(in)        :: x
    real(wp),dimension(:),intent(in)        :: y
    real(wp),dimension(:,:),intent(in)      :: fcn
    integer,intent(in)                      :: kx
    integer,intent(in)                      :: ky
    integer,intent(out)                     :: iflag
    
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
    
    end subroutine initialize_2d
!*****************************************************************************************
    
!*****************************************************************************************
!> Evaluate a [[bspline_2d]] interpolate.  This is a wrapper for [[db2val]].

    subroutine evaluate_2d(me,xval,yval,idx,idy,f,iflag)
    
    implicit none
    
    class(bspline_2d),intent(in) :: me
    real(wp),intent(in)         :: xval
    real(wp),intent(in)         :: yval
    integer,intent(in)          :: idx
    integer,intent(in)          :: idy
    real(wp),intent(out)        :: f
    integer,intent(out)         :: iflag
    
    call db2val(xval,yval,idx,idy,me%tx,me%ty,me%nx,me%ny,me%kx,me%ky,me%bcoef,f,iflag)
    
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
!> Initialize a [[bspline_3d]] type.  This is a wrapper for [[db3ink]].

    subroutine initialize_3d(me,x,y,z,fcn,kx,ky,kz,iflag)

    implicit none
    
    class(bspline_3d),intent(inout)         :: me
    real(wp),dimension(:),intent(in)        :: x
    real(wp),dimension(:),intent(in)        :: y
    real(wp),dimension(:),intent(in)        :: z
    real(wp),dimension(:,:,:),intent(in)    :: fcn
    integer,intent(in)                      :: kx
    integer,intent(in)                      :: ky
    integer,intent(in)                      :: kz
    integer,intent(out)                     :: iflag
    
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
    
    end subroutine initialize_3d
!*****************************************************************************************
    
!*****************************************************************************************
!> Evaluate a [[bspline_3d]] interpolate.  This is a wrapper for [[db3val]].

    subroutine evaluate_3d(me,xval,yval,zval,idx,idy,idz,f,iflag)
    
    implicit none
    
    class(bspline_3d),intent(in) :: me
    real(wp),intent(in)         :: xval
    real(wp),intent(in)         :: yval
    real(wp),intent(in)         :: zval
    integer,intent(in)          :: idx
    integer,intent(in)          :: idy
    integer,intent(in)          :: idz
    real(wp),intent(out)        :: f
    integer,intent(out)         :: iflag
    
    call db3val(xval,yval,zval,&
                idx,idy,idz,&
                me%tx,me%ty,me%tz,&
                me%nx,me%ny,me%nz,&
                me%kx,me%ky,me%kz,&
                me%bcoef,f,iflag)
    
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
!> Initialize a [[bspline_4d]] type.  This is a wrapper for [[db4ink]].

    subroutine initialize_4d(me,x,y,z,q,fcn,kx,ky,kz,kq,iflag)

    implicit none
    
    class(bspline_4d),intent(inout)            :: me
    real(wp),dimension(:),intent(in)           :: x
    real(wp),dimension(:),intent(in)           :: y
    real(wp),dimension(:),intent(in)           :: z
    real(wp),dimension(:),intent(in)           :: q
    real(wp),dimension(:,:,:,:),intent(in)     :: fcn
    integer,intent(in)                         :: kx
    integer,intent(in)                         :: ky
    integer,intent(in)                         :: kz
    integer,intent(in)                         :: kq
    integer,intent(out)                        :: iflag
    
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
    
    end subroutine initialize_4d
!*****************************************************************************************
    
!*****************************************************************************************
!> Evaluate a [[bspline_4d]] interpolate.  This is a wrapper for [[db4val]].

    subroutine evaluate_4d(me,xval,yval,zval,qval,idx,idy,idz,idq,f,iflag)
    
    implicit none
    
    class(bspline_4d),intent(in) :: me
    real(wp),intent(in)         :: xval
    real(wp),intent(in)         :: yval
    real(wp),intent(in)         :: zval
    real(wp),intent(in)         :: qval
    integer,intent(in)          :: idx
    integer,intent(in)          :: idy
    integer,intent(in)          :: idz
    integer,intent(in)          :: idq
    real(wp),intent(out)        :: f
    integer,intent(out)         :: iflag
    
    call db4val(xval,yval,zval,qval,&
                idx,idy,idz,idq,&
                me%tx,me%ty,me%tz,me%tq,&
                me%nx,me%ny,me%nz,me%nq,&
                me%kx,me%ky,me%kz,me%kq,&
                me%bcoef,f,iflag)
    
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
!> Initialize a [[bspline_5d]] type.  This is a wrapper for [[db5ink]].

    subroutine initialize_5d(me,x,y,z,q,r,fcn,kx,ky,kz,kq,kr,iflag)

    implicit none
    
    class(bspline_5d),intent(inout)               :: me
    real(wp),dimension(:),intent(in)              :: x
    real(wp),dimension(:),intent(in)              :: y
    real(wp),dimension(:),intent(in)              :: z
    real(wp),dimension(:),intent(in)              :: q
    real(wp),dimension(:),intent(in)              :: r
    real(wp),dimension(:,:,:,:,:),intent(in)      :: fcn
    integer,intent(in)                            :: kx
    integer,intent(in)                            :: ky
    integer,intent(in)                            :: kz
    integer,intent(in)                            :: kq
    integer,intent(in)                            :: kr
    integer,intent(out)                           :: iflag
    
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
    
    end subroutine initialize_5d
!*****************************************************************************************
    
!*****************************************************************************************
!> Evaluate a [[bspline_5d]] interpolate.  This is a wrapper for [[db5val]].

    subroutine evaluate_5d(me,xval,yval,zval,qval,rval,idx,idy,idz,idq,idr,f,iflag)
    
    implicit none
    
    class(bspline_5d),intent(in) :: me
    real(wp),intent(in)         :: xval
    real(wp),intent(in)         :: yval
    real(wp),intent(in)         :: zval
    real(wp),intent(in)         :: qval
    real(wp),intent(in)         :: rval
    integer,intent(in)          :: idx
    integer,intent(in)          :: idy
    integer,intent(in)          :: idz
    integer,intent(in)          :: idq
    integer,intent(in)          :: idr
    real(wp),intent(out)        :: f
    integer,intent(out)         :: iflag
    
    call db5val(xval,yval,zval,qval,rval,&
                idx,idy,idz,idq,idr,&
                me%tx,me%ty,me%tz,me%tq,me%tr,&
                me%nx,me%ny,me%nz,me%nq,me%nr,&
                me%kx,me%ky,me%kz,me%kq,me%kr,&
                me%bcoef,f,iflag)
    
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
!> Initialize a [[bspline_6d]] type.  This is a wrapper for [[db6ink]].

    subroutine initialize_6d(me,x,y,z,q,r,s,fcn,kx,ky,kz,kq,kr,ks,iflag)

    implicit none
    
    class(bspline_6d),intent(inout)                  :: me
    real(wp),dimension(:),intent(in)                 :: x
    real(wp),dimension(:),intent(in)                 :: y
    real(wp),dimension(:),intent(in)                 :: z
    real(wp),dimension(:),intent(in)                 :: q
    real(wp),dimension(:),intent(in)                 :: r
    real(wp),dimension(:),intent(in)                 :: s
    real(wp),dimension(:,:,:,:,:,:),intent(in)       :: fcn
    integer,intent(in)                               :: kx
    integer,intent(in)                               :: ky
    integer,intent(in)                               :: kz
    integer,intent(in)                               :: kq
    integer,intent(in)                               :: kr
    integer,intent(in)                               :: ks
    integer,intent(out)                              :: iflag
    
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
    
    end subroutine initialize_6d
!*****************************************************************************************
    
!*****************************************************************************************
!> Evaluate a [[bspline_6d]] interpolate.  This is a wrapper for [[db6val]].

    subroutine evaluate_6d(me,xval,yval,zval,qval,rval,sval,idx,idy,idz,idq,idr,ids,f,iflag)
    
    implicit none
    
    class(bspline_6d),intent(in) :: me
    real(wp),intent(in)         :: xval
    real(wp),intent(in)         :: yval
    real(wp),intent(in)         :: zval
    real(wp),intent(in)         :: qval
    real(wp),intent(in)         :: rval
    real(wp),intent(in)         :: sval
    integer,intent(in)          :: idx
    integer,intent(in)          :: idy
    integer,intent(in)          :: idz
    integer,intent(in)          :: idq
    integer,intent(in)          :: idr
    integer,intent(in)          :: ids
    real(wp),intent(out)        :: f
    integer,intent(out)         :: iflag
    
    call db6val(xval,yval,zval,qval,rval,sval,&
                idx,idy,idz,idq,idr,ids,&
                me%tx,me%ty,me%tz,me%tq,me%tr,me%ts,&
                me%nx,me%ny,me%nz,me%nq,me%nr,me%ns,&
                me%kx,me%ky,me%kz,me%kq,me%kr,me%ks,&
                me%bcoef,f,iflag)
    
    end subroutine evaluate_6d
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_oo_module
!*****************************************************************************************