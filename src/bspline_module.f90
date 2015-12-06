!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!# Description
!
!  Multidimensional (2D-6D) B-Spline interpolation of data on a regular grid.
!  This module uses both the subroutine and object-oriented modules.
!
!# iflag Codes
!
! The following are the `iflag` codes that can be issued by the various
! routines. They are divided into the initialize routines
! (e.g., [[db1ink]], [[db2ink]], etc.)
! and the evaluate routines (e.g., [[db1val]], [[db2val]], etc.).
!
!## Initialize Routines
!```
!db*ink    1: successful execution.
!          2: iflag out of range
!          3: nx out of range
!          4: kx out of range
!          5: x not strictly increasing
!          6: tx not non-decreasing
!          7: ny out of range
!          8: ky out of range
!          9: y not strictly increasing
!         10: ty not non-decreasing
!         11: nz out of range
!         12: kz out of range
!         13: z not strictly increasing
!         14: tz not non-decreasing
!         15: nq out of range
!         16: kq out of range
!         17: q not strictly increasing
!         18: tq not non-decreasing
!         19: nr out of range
!         20: kr out of range
!         21: r not strictly increasing
!         22: tr not non-decreasing
!
!dbintk  100: k does not satisfy k>=1
!        101: n does not satisfy n>=k
!        102: x(i) does not satisfy x(i)<x(i+1) for some i
!        103: some abscissa was not in the support of the
!             corresponding basis function and the system is singular
!        104: the system of solver detects a singular system
!             although the theoretical conditions for a solution were satisfied
!
!dbspvn  201: k does not satisfy k>=1
!        202: jhigh does not satisfy 1<=jhigh<=k
!        203: index is not 1 or 2
!        204: x does not satisfy t(ileft)<=x<=t(ileft+1)
!
!dbtpcf  301: n should be >0
!```
!
!## Evaluate Routines
!```
!db*val   0: successful execution
!         1: x value out of bounds
!         2: y value out of bounds
!         3: z value out of bounds
!         4: q value out of bounds
!         5: r value out of bounds
!         6: s value out of bounds
!
!dbvalu 401: k does not satisfy k>=1
!       402: n does not satisfy n>=k
!       403: ideriv does not satisfy 0<=ideriv<k
!       404: x is not greater than or equal to t(k)
!       405: x is not less than or equal to t(n+1)
!       406: a left limiting value cannot be obtained at t(k)
!```

    module bspline_module

    use bspline_oo_module
    use bspline_sub_module

    implicit none

    public

    end module bspline_module
!*****************************************************************************************
