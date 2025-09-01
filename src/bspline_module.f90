!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!
!  Multidimensional (1D-6D) B-Spline interpolation of data on a regular grid.
!  This module uses both the subroutine and object-oriented modules.

    module bspline_module

    use bspline_kinds_module, only: bspline_wp => wp
    use bspline_oo_module
    use bspline_sub_module
    use bspline_defc_module

    implicit none

    public

!*****************************************************************************************
    end module bspline_module
!*****************************************************************************************
