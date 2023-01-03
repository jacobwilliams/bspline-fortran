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
!  univ. of california berkeley, univ. of colorado denver and nag ltd.
!```
!
!### See also
!  * [BLAS Sourcecode](https://github.com/Reference-LAPACK/lapack/tree/master/BLAS/SRC)

module bspline_blas_module

#ifndef HAS_BLAS

   use bspline_kinds_module, only: wp, ip

   implicit none
   private

   public :: daxpy,dcopy,dscal,dswap,ddot,dnrm2,dasum,idamax,drotm,drotmg

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
         ! code for both increments equal to 1
         ! clean-up loop
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
         ! code for unequal increments or equal increments
         ! not equal to 1
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
         ! code for both increments equal to 1
         ! clean-up loop
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
         ! code for unequal increments or equal increments
         ! not equal to 1
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
         ! code for increment equal to 1
         ! clean-up loop
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
         ! code for increment not equal to 1
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
         ! code for both increments equal to 1
         ! clean-up loop
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
         ! code for unequal increments or equal increments not equal
         ! to 1
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

   function dnrm2( n, x, incx )
      !! returns the euclidean norm of a vector
         real(wp) :: dnrm2

         real(wp), parameter :: zero = 0.0_wp
         real(wp), parameter :: one  = 1.0_wp
         real(wp), parameter :: maxN = huge(0.0_wp)

         real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
             (minexponent(0._wp) - 1) * 0.5_wp)
         real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
             (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
         real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
             (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
         real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
             (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))

         integer(ip) :: incx, n
         real(wp) :: x(*)

         integer(ip) :: i, ix
         logical :: notbig
         real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
      !
      !  Quick return if possible
      !
         DNRM2 = zero
         if( n <= 0 ) return
      !
         scl = one
         sumsq = zero
      !
      !  Compute the sum of squares in 3 accumulators:
      !     abig -- sums of squares scaled down to avoid overflow
      !     asml -- sums of squares scaled up to avoid underflow
      !     amed -- sums of squares that do not require scaling
      !  The thresholds and multipliers are
      !     tbig -- values bigger than this are scaled down by sbig
      !     tsml -- values smaller than this are scaled up by ssml
      !
         notbig = .true.
         asml = zero
         amed = zero
         abig = zero
         ix = 1
         if( incx < 0 ) ix = 1 - (n-1)*incx
         do i = 1, n
            ax = abs(x(ix))
            if (ax > tbig) then
               abig = abig + (ax*sbig)**2
               notbig = .false.
            else if (ax < tsml) then
               if (notbig) asml = asml + (ax*ssml)**2
            else
               amed = amed + ax**2
            end if
            ix = ix + incx
         end do
      !
      !  Combine abig and amed or amed and asml if more than one
      !  accumulator was used.
      !
         if (abig > zero) then
      !
      !     Combine abig and amed if abig > 0.
      !
            if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
               abig = abig + (amed*sbig)*sbig
            end if
            scl = one / sbig
            sumsq = abig
         else if (asml > zero) then
      !
      !     Combine amed and asml if asml > 0.
      !
            if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
               amed = sqrt(amed)
               asml = sqrt(asml) / ssml
               if (asml > amed) then
                  ymin = amed
                  ymax = asml
               else
                  ymin = asml
                  ymax = amed
               end if
               scl = one
               sumsq = ymax**2*( one + (ymin/ymax)**2 )
            else
               scl = one / ssml
               sumsq = asml
            end if
         else
      !
      !     Otherwise all values are mid-range
      !
            scl = one
            sumsq = amed
         end if
         DNRM2 = scl*sqrt( sumsq )
         return
   end function

   real(wp) function dasum(n,dx,incx)
         !! dasum takes the sum of the absolute values.

         integer(ip) :: incx,n
         real(wp) :: dx(*)

         real(wp) dtemp
         integer(ip) i,m,mp1,nincx

         dasum = 0.0_wp
         dtemp = 0.0_wp
         if (n<=0 .or. incx<=0) return
         if (incx==1) then
            ! code for increment equal to 1
            ! clean-up loop
            m = mod(n,6)
            if (m/=0) then
               do i = 1,m
                  dtemp = dtemp + abs(dx(i))
               end do
               if (n<6) then
                  dasum = dtemp
                  return
               end if
            end if
            mp1 = m + 1
            do i = mp1,n,6
               dtemp = dtemp + abs(dx(i)) + abs(dx(i+1)) + &
                       abs(dx(i+2)) + abs(dx(i+3)) + &
                       abs(dx(i+4)) + abs(dx(i+5))
            end do
         else
            ! code for increment not equal to 1
            nincx = n*incx
            do i = 1,nincx,incx
               dtemp = dtemp + abs(dx(i))
            end do
         end if
         dasum = dtemp

   end function dasum

   integer function idamax(n,dx,incx)
         !! idamax finds the index of the first element having maximum absolute value.

         integer(ip) :: incx,n
         real(wp) :: dx(*)

         real(wp) :: dmax
         integer(ip) :: i,ix

         idamax = 0
         if (n<1 .or. incx<=0) return
         idamax = 1
         if (n==1) return
         if (incx==1) then
            ! code for increment equal to 1
            dmax = abs(dx(1))
            do i = 2,n
               if (abs(dx(i))>dmax) then
                  idamax = i
                  dmax = abs(dx(i))
               end if
            end do
         else
            ! code for increment not equal to 1
            ix = 1
            dmax = abs(dx(1))
            ix = ix + incx
            do i = 2,n
               if (abs(dx(ix))>dmax) then
                  idamax = i
                  dmax = abs(dx(ix))
               end if
               ix = ix + incx
            end do
         end if

   end function idamax

   subroutine drotm(n,dx,incx,dy,incy,dparam)
      !! apply the modified givens transformation, H, to the 2 by n matrix

      integer(ip) :: incx,incy,n
      real(wp) :: dparam(5),dx(*),dy(*)

      real(wp) :: dflag,dh11,dh12,dh21,dh22,w,z
      integer(ip) :: i,kx,ky,nsteps

      real(wp),parameter :: zero = 0.0_wp
      real(wp),parameter :: two = 2.0_wp

      dflag = dparam(1)
      if (n<=0 .or. (dflag+two==zero)) return
      if (incx==incy.and.incx>0) then

         nsteps = n*incx
         if (dflag<zero) then
            dh11 = dparam(2)
            dh12 = dparam(4)
            dh21 = dparam(3)
            dh22 = dparam(5)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w*dh11 + z*dh12
               dy(i) = w*dh21 + z*dh22
            end do
         else if (dflag==zero) then
            dh12 = dparam(4)
            dh21 = dparam(3)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w + z*dh12
               dy(i) = w*dh21 + z
            end do
         else
            dh11 = dparam(2)
            dh22 = dparam(5)
            do i = 1,nsteps,incx
               w = dx(i)
               z = dy(i)
               dx(i) = w*dh11 + z
               dy(i) = -w + dh22*z
            end do
         end if
      else
         kx = 1
         ky = 1
         if (incx<0) kx = 1 + (1-n)*incx
         if (incy<0) ky = 1 + (1-n)*incy
         if (dflag<zero) then
            dh11 = dparam(2)
            dh12 = dparam(4)
            dh21 = dparam(3)
            dh22 = dparam(5)
            do i = 1,n
               w = dx(kx)
               z = dy(ky)
               dx(kx) = w*dh11 + z*dh12
               dy(ky) = w*dh21 + z*dh22
               kx = kx + incx
               ky = ky + incy
            end do
         else if (dflag==zero) then
            dh12 = dparam(4)
            dh21 = dparam(3)
            do i = 1,n
               w = dx(kx)
               z = dy(ky)
               dx(kx) = w + z*dh12
               dy(ky) = w*dh21 + z
               kx = kx + incx
               ky = ky + incy
            end do
         else
               dh11 = dparam(2)
               dh22 = dparam(5)
               do i = 1,n
                  w = dx(kx)
                  z = dy(ky)
                  dx(kx) = w*dh11 + z
                  dy(ky) = -w + dh22*z
                  kx = kx + incx
                  ky = ky + incy
            end do
         end if
      end if

   end subroutine drotm

   subroutine drotmg(dd1,dd2,dx1,dy1,dparam)
      !! construct the modified givens transformation matrix H

      real(wp) :: dd1,dd2,dx1,dy1
      real(wp) :: dparam(5)

      real(wp) :: dflag,dh11,dh12,dh21,dh22,dp1,dp2,dq1,dq2,dtemp,&
                  du

      real(wp),parameter :: zero = 0.0_wp
      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: two = 2.0_wp
      real(wp),parameter :: gam    = 4096.0_wp
      real(wp),parameter :: gamsq  = gam*gam    !! 16777216.0_wp
      real(wp),parameter :: rgamsq = one/gamsq  !! 5.9604645e-8_wp

         if (dd1<zero) then
            ! go zero-h-d-and-dx1..
            dflag = -one
            dh11 = zero
            dh12 = zero
            dh21 = zero
            dh22 = zero
            dd1 = zero
            dd2 = zero
            dx1 = zero
         else
            ! case-dd1-nonnegative
            dp2 = dd2*dy1
            if (dp2==zero) then
               dflag = -two
               dparam(1) = dflag
               return
            end if
            ! regular-case..
            dp1 = dd1*dx1
            dq2 = dp2*dy1
            dq1 = dp1*dx1
            if (abs(dq1)>abs(dq2)) then
               dh21 = -dy1/dx1
               dh12 = dp2/dp1
               du = one - dh12*dh21
               if (du>zero) then
                  dflag = zero
                  dd1 = dd1/du
                  dd2 = dd2/du
                  dx1 = dx1*du
               else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  dflag = -one
                  dh11 = zero
                  dh12 = zero
                  dh21 = zero
                  dh22 = zero
                  dd1 = zero
                  dd2 = zero
                  dx1 = zero
               end if
            else

               if (dq2<zero) then
                  ! go zero-h-d-and-dx1..
                  dflag = -one
                  dh11 = zero
                  dh12 = zero
                  dh21 = zero
                  dh22 = zero
                  dd1 = zero
                  dd2 = zero
                  dx1 = zero
               else
                  dflag = one
                  dh11 = dp1/dp2
                  dh22 = dx1/dy1
                  du = one + dh11*dh22
                  dtemp = dd2/du
                  dd2 = dd1/du
                  dd1 = dtemp
                  dx1 = dy1*du
               end if
            end if

            ! procedure..scale-check
            if (dd1/=zero) then
               do while ((dd1<=rgamsq) .or. (dd1>=gamsq))
                  if (dflag==zero) then
                     dh11 = one
                     dh22 = one
                     dflag = -one
                  else
                     dh21 = -one
                     dh12 = one
                     dflag = -one
                  end if
                  if (dd1<=rgamsq) then
                     dd1 = dd1*gam**2
                     dx1 = dx1/gam
                     dh11 = dh11/gam
                     dh12 = dh12/gam
                  else
                     dd1 = dd1/gam**2
                     dx1 = dx1*gam
                     dh11 = dh11*gam
                     dh12 = dh12*gam
                  end if
               enddo
            end if

            if (dd2/=zero) then
               do while ( (abs(dd2)<=rgamsq) .or. (abs(dd2)>=gamsq) )
                  if (dflag==zero) then
                     dh11 = one
                     dh22 = one
                     dflag = -one
                  else
                     dh21 = -one
                     dh12 = one
                     dflag = -one
                  end if
                  if (abs(dd2)<=rgamsq) then
                     dd2 = dd2*gam**2
                     dh21 = dh21/gam
                     dh22 = dh22/gam
                  else
                     dd2 = dd2/gam**2
                     dh21 = dh21*gam
                     dh22 = dh22*gam
                  end if
               end do
            end if

         end if

         if (dflag<zero) then
            dparam(2) = dh11
            dparam(3) = dh21
            dparam(4) = dh12
            dparam(5) = dh22
         else if (dflag==zero) then
            dparam(3) = dh21
            dparam(4) = dh12
         else
            dparam(2) = dh11
            dparam(5) = dh22
         end if

         dparam(1) = dflag

      end subroutine drotmg

#endif

   end module bspline_blas_module
