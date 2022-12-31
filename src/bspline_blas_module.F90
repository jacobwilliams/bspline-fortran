!*****************************************************************************************
!>
!  BLAS procedures, which can be use used if not linking with a BLAS library,
!  if one is not available, or if a real kind /= `real64` is required.
!
!  The original code has been slightly modernized.
!
!### Notes
!```
!  reference blas level1 routines
!  reference blas is a software package provided by univ. of tennessee,
!  univ. of california berkeley, univ. of colorado denver and nag ltd..
!```
!
!### See also
!  * [BLAS Sourcecode](https://github.com/Reference-LAPACK/lapack/tree/master/BLAS/SRC)

module bspline_blas_module
#ifndef HAS_BLAS

   use bspline_kinds_module, only: wp, ip

   implicit none

   private

   public :: daxpy,dcopy,dscal,dswap,ddot

contains

   subroutine daxpy(n, da, dx, incx, dy, incy)
         !! DAXPY constant times a vector plus a vector.
         !! uses unrolled loops for increments equal to one.

      real(wp) :: da
      integer(ip) :: incx, incy, n
      real(wp) :: dx(*), dy(*)

      integer(ip) :: i, ix, iy, m, mp1

      if (n <= 0_ip) return
      if (da == 0.0_wp) return
      if (incx == 1_ip .and. incy == 1_ip) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n, 4_ip)
         if (m /= 0_ip) then
            do i = 1_ip, m
               dy(i) = dy(i) + da*dx(i)
            end do
         end if
         if (n < 4_ip) return
         mp1 = m + 1_ip
         do i = mp1, n, 4_ip
            dy(i) = dy(i) + da*dx(i)
            dy(i + 1_ip) = dy(i + 1_ip) + da*dx(i + 1_ip)
            dy(i + 2_ip) = dy(i + 2_ip) + da*dx(i + 2_ip)
            dy(i + 3_ip) = dy(i + 3_ip) + da*dx(i + 3_ip)
         end do
      else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1_ip
         iy = 1_ip
         if (incx < 0_ip) ix = (-n + 1_ip)*incx + 1_ip
         if (incy < 0_ip) iy = (-n + 1_ip)*incy + 1_ip
         do i = 1_ip, n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do
      end if

   end subroutine daxpy

   subroutine dcopy(n, dx, incx, dy, incy)
         !! DCOPY copies a vector, x, to a vector, y.
         !! uses unrolled loops for increments equal to 1.

      integer(ip) :: incx, incy, n
      real(wp) :: dx(*), dy(*)

      integer(ip) :: i, ix, iy, m, mp1

      if (n <= 0_ip) return
      if (incx == 1_ip .and. incy == 1_ip) then
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n, 7_ip)
         if (m /= 0_ip) then
            do i = 1_ip, m
               dy(i) = dx(i)
            end do
            if (n < 7_ip) return
         end if
         mp1 = m + 1_ip
         do i = mp1, n, 7_ip
            dy(i) = dx(i)
            dy(i + 1_ip) = dx(i + 1_ip)
            dy(i + 2_ip) = dx(i + 2_ip)
            dy(i + 3_ip) = dx(i + 3_ip)
            dy(i + 4_ip) = dx(i + 4_ip)
            dy(i + 5_ip) = dx(i + 5_ip)
            dy(i + 6_ip) = dx(i + 6_ip)
         end do
      else
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1_ip
         iy = 1_ip
         if (incx < 0_ip) ix = (-n + 1_ip)*incx + 1_ip
         if (incy < 0_ip) iy = (-n + 1_ip)*incy + 1_ip
         do i = 1_ip, n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do
      end if

   end subroutine dcopy

   subroutine dscal(n, da, dx, incx)
         !! DSCAL scales a vector by a constant.
         !! uses unrolled loops for increment equal to 1.

      real(wp) :: da
      integer(ip) :: incx, n
      real(wp) :: dx(*)

      integer i, m, mp1, nincx

      if (n <= 0_ip .or. incx <= 0_ip) return
      if (incx == 1_ip) then
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = mod(n, 5_ip)
         if (m /= 0_ip) then
            do i = 1_ip, m
               dx(i) = da*dx(i)
            end do
            if (n < 5_ip) return
         end if
         mp1 = m + 1_ip
         do i = mp1, n, 5_ip
            dx(i) = da*dx(i)
            dx(i + 1_ip) = da*dx(i + 1_ip)
            dx(i + 2_ip) = da*dx(i + 2_ip)
            dx(i + 3_ip) = da*dx(i + 3_ip)
            dx(i + 4_ip) = da*dx(i + 4_ip)
         end do
      else
!
!        code for increment not equal to 1
!
         nincx = n*incx
         do i = 1_ip, nincx, incx
            dx(i) = da*dx(i)
         end do
      end if

   end subroutine dscal

   subroutine dswap(n, dx, incx, dy, incy)
         !! DSWAP interchanges two vectors.
         !! uses unrolled loops for increments equal to 1.

      integer(ip) :: incx, incy, n
      real(wp) :: dx(*), dy(*)

      real(wp) :: dtemp
      integer(ip) :: i, ix, iy, m, mp1

      if (n <= 0_ip) return
      if (incx == 1_ip .and. incy == 1_ip) then
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
         m = mod(n, 3_ip)
         if (m /= 0_ip) then
            do i = 1_ip, m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            end do
            if (n < 3_ip) return
         end if
         mp1 = m + 1_ip
         do i = mp1, n, 3_ip
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i + 1_ip)
            dx(i + 1_ip) = dy(i + 1_ip)
            dy(i + 1_ip) = dtemp
            dtemp = dx(i + 2_ip)
            dx(i + 2_ip) = dy(i + 2_ip)
            dy(i + 2_ip) = dtemp
         end do
      else
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         if (incx < 0_ip) ix = (-n + 1_ip)*incx + 1_ip
         if (incy < 0_ip) iy = (-n + 1_ip)*incy + 1_ip
         do i = 1_ip, n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         end do
      end if

   end subroutine dswap

   real(wp) function ddot(n, dx, incx, dy, incy)
         !! ddot forms the dot product of two vectors.
         !! uses unrolled loops for increments equal to one.

      integer(ip) :: incx, incy, n
      real(wp) :: dx(*), dy(*)

      real(wp) :: dtemp
      integer(ip) :: i, ix, iy, m, mp1

      ddot = 0.0_wp
      dtemp = 0.0_wp
      if (n <= 0_ip) return
      if (incx == 1_ip .and. incy == 1_ip) then
         ! code for both increments equal to 1
         ! clean-up loop
         m = mod(n, 5_ip)
         if (m /= 0_ip) then
            do i = 1_ip, m
               dtemp = dtemp + dx(i)*dy(i)
            end do
            if (n < 5_ip) then
               ddot = dtemp
               return
            end if
         end if
         mp1 = m + 1_ip
         do i = mp1, n, 5_ip
            dtemp = dtemp + dx(i)*dy(i) + &
                    dx(i + 1_ip)*dy(i + 1_ip) + dx(i + 2_ip)*dy(i + 2_ip) + &
                    dx(i + 3_ip)*dy(i + 3_ip) + dx(i + 4_ip)*dy(i + 4_ip)
         end do
      else
         ! code for unequal increments or equal increments
         ! not equal to 1
         ix = 1_ip
         iy = 1_ip
         if (incx < 0_ip) ix = (-n + 1_ip)*incx + 1_ip
         if (incy < 0_ip) iy = (-n + 1_ip)*incy + 1_ip
         do i = 1_ip, n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         end do
      end if
      ddot = dtemp

      end function ddot

#endif

   end module bspline_blas_module
