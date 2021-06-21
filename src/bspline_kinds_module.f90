!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!  Numeric kind definitions for BSpline-Fortran.

    module bspline_kinds_module

    use,intrinsic :: iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter,public :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter,public :: wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter,public :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

#ifdef INT8
    integer,parameter,public :: ip = int8     !! Integer working precision [1 byte]
#elif INT16
    integer,parameter,public :: ip = int16    !! Integer working precision [2 bytes]
#elif INT32
    integer,parameter,public :: ip = int32    !! Integer working precision [4 bytes]
#elif INT64
    integer,parameter,public :: ip = int64    !! Integer working precision [8 bytes]
#else
    integer,parameter,public :: ip = int32    !! Integer working precision if not specified [4 bytes]
#endif

!*****************************************************************************************
    end module bspline_kinds_module
!*****************************************************************************************
