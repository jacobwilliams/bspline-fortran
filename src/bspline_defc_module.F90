!*****************************************************************************************
!>
!  [[defc]] and [[dfc]] procedures and support routines from [SLATEC](https://netlib.org/slatec/src/).
!  For fitting B-splines polynomials to discrete 1D data.
!
!### References
!
!  For a description of the B-splines and usage instructions to
!  evaluate them, see:
!
!   * C. W. de Boor, Package for Calculating with B-Splines.
!     SIAM J. Numer. Anal., p. 441, (June, 1977).
!
!  For further discussion of (constrained) curve fitting using
!    B-splines, see reference 2.
!
!   * R. J. Hanson, Constrained least squares curve fitting
!     to discrete data using B-splines, a users guide,
!     Report SAND78-1291, Sandia Laboratories, December
!     1978.
!
!### History
!  * Dec 2022 (Jacob Williams) : Cleanup and modernization of the SLATEC routines.
!
!@note This module does not support the user-defined `ip` integer kind.
!      It only uses the default integer kind.
!
!@todo add `iflag` outputs to be consistent with the rest of the library.

module bspline_defc_module

   use bspline_kinds_module, only: wp !, ip
#ifndef HAS_BLAS
   use bspline_blas_module
#endif

   implicit none

   private

   real(wp),parameter :: drelpr = epsilon(1.0_wp) !! machine precision (`d1mach(4)`)

#ifdef HAS_BLAS
   ! user is linking against an external BLAS library
   double precision,external :: ddot, dasum
   integer,external :: idamax
   external :: daxpy,dcopy,dscal,dswap,dnrm2,drotm,drotmg
#endif

   public :: defc, dfc, dcv

   contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subprogram fits a piecewise polynomial curve
!  to discrete data.  The piecewise polynomials are
!  represented as B-splines.
!  The fitting is done in a weighted least squares sense.
!
!  The data can be processed in groups of modest size.
!  The size of the group is chosen by the user.  This feature
!  may be necessary for purposes of using constrained curve fitting
!  with subprogram [[DFC]] on a very large data set.
!
!### Evaluating the Fitted Curve
!
!  To evaluate derivative number `IDER` at `XVAL`,
!  use the function subprogram [[DBVALU]].
!
!```fortran
! f = dbvalu(bkpt,coeff,nbkpt-nord,nord,ider,xval,inbv,workb)
!```
!
!  The output of this subprogram will not be
!  defined unless an output value of `MDEOUT=1`
!  was obtained from [[DEFC]], `XVAL` is in the data
!  interval, and `IDER` is nonnegative and `< NORD`.
!
!  The first time [[DBVALU]] is called, `INBV=1`
!  must be specified.  This value of `INBV` is the
!  overwritten by [[DBVALU]].  The array `WORKB(*)`
!  must be of length at least `3*NORD`, and must
!  not be the same as the `W(*)` array used in the
!  call to [[DEFC]].
!
!  [[DBVALU]] expects the breakpoint array `BKPT(*)`
!  to be sorted.
!
!### Revision history
!
!   * 800801  DATE WRITTEN.
!     WRITTEN BY R. HANSON, SANDIA NATL. LABS.,
!     ALB., N. M., AUGUST-SEPTEMBER, 1980.
!   * 890531  Changed all specific intrinsics to generic.  (WRB)
!   * 890531  REVISION DATE from Version 3.2
!   * 891214  Prologue converted to Version 4.0 format.  (BAB)
!   * 900510  Change Prologue comments to refer to XERMSG.  (RWC)
!   * 900607  Editorial changes to Prologue to make Prologues for EFC,
!     DEFC, FC, and DFC look as much the same as possible.  (RWC)
!   * 920501  Reformatted the REFERENCES section.  (WRB)
!   * Jacob Williams, 2022 : modernized

   subroutine defc(Ndata, Xdata, Ydata, Sddata, Nord, Nbkpt, Bkpt, Mdein, &
                   Mdeout, Coeff, Lw, w)

      integer,intent(in) :: Ndata !! number of points (size of `xdata` and `ydata`).
                                  !! Any non-negative value of `NDATA` is allowed.
                                  !! A negative value of `NDATA` is an error.
      real(wp),dimension(ndata), intent(in) :: Xdata !! X data array. No sorting of `XDATA(*)` is required.
      real(wp),dimension(ndata), intent(in) :: Ydata !! Y data array.
      real(wp),dimension(ndata), intent(in) :: Sddata !! Y value standard deviation or uncertainty.
                                                      !! A zero value for any entry of
                                                      !! `SDDATA(*)` will weight that data point as 1.
                                                      !! Otherwise the weight of that data point is
                                                      !! the reciprocal of this entry.
      integer,intent(in) :: Nord !! B-spline order.
                                 !! (The order of the spline is one more than the
                                 !! degree of the piecewise polynomial defined on
                                 !! each interval.  This is consistent with the
                                 !! B-spline package convention.  For example,
                                 !! `NORD=4` when we are using piecewise cubics.)
                                 !! `NORD` must be in the range `1 <= NORD <= 20`.
      integer,intent(in) :: Nbkpt !! The value of `NBKPT` must satisfy the condition `NBKPT >= 2*NORD`.
      real(wp),dimension(:),intent(in) :: Bkpt !! `NBKPT` knots of the B-spline.
                                     !! Normally the
                                     !! problem data interval will be included between
                                     !! the limits `BKPT(NORD)` and `BKPT(NBKPT-NORD+1)`.
                                     !! The additional end knots `BKPT(I),I=1,...,NORD-1`
                                     !! and `I=NBKPT-NORD+2,...,NBKPT`, are
                                     !! required to compute the functions used to fit
                                     !! the data.  No sorting of `BKPT(*)` is required.
                                     !! Internal to [[DEFC]] the extreme end knots may
                                     !! be reduced and increased respectively to
                                     !! accommodate any data values that are exterior
                                     !! to the given knot values.  The contents of
                                     !! `BKPT(*)` is not changed.
      integer,intent(in) :: Mdein !! An integer flag, with one of two possible
                                  !! values (1 or 2), that directs the subprogram
                                  !! action with regard to new data points provided
                                  !! by the user:
                                  !!
                                  !! * `= 1`  The first time that [[DEFC]] has been
                                  !!   entered.  There are NDATA points to process.
                                  !! * `= 2`  This is another entry to DEFC().  The
                                  !!   subprogram [[DEFC]] has been entered with MDEIN=1
                                  !!   exactly once before for this problem.  There
                                  !!   are NDATA new additional points to merge and
                                  !!   process with any previous points.
                                  !!   (When using [[DEFC]] with MDEIN=2 it is
                                  !!   important that the set of knots remain fixed at the
                                  !!   same values for all entries to [[DEFC]].)
      integer,intent(out) :: Mdeout !! An output flag that indicates the status
                                    !! of the curve fit:
                                    !!
                                    !!  * `=-1`  A usage error of [[DEFC]] occurred.  The
                                    !!    offending condition is noted with the SLATEC
                                    !!    library error processor, `XERMSG( )`.  In case
                                    !!    the working array `W(*)` is not long enough, the
                                    !!    minimal acceptable length is printed.
                                    !!
                                    !!  * `=1`  The B-spline coefficients for the fitted
                                    !!    curve have been returned in array `COEFF(*)`.
                                    !!
                                    !!  * `=2`  Not enough data has been processed to
                                    !!    determine the B-spline coefficients.
                                    !!    The user has one of two options.  Continue
                                    !!    to process more data until a unique set
                                    !!    of coefficients is obtained, or use the
                                    !!    subprogram [[DFC]] to obtain a specific
                                    !!    set of coefficients.  The user should read
                                    !!    the usage instructions for [[DFC]] for further
                                    !!    details if this second option is chosen.
      real(wp),intent(out) :: Coeff(*) !! If the output value of `MDEOUT=1`, this array
                                       !! contains the unknowns obtained from the least
                                       !! squares fitting process.  These `N=NBKPT-NORD`
                                       !! parameters are the B-spline coefficients.
                                       !! For `MDEOUT=2`, not enough data was processed to
                                       !! uniquely determine the B-spline coefficients.
                                       !! In this case, and also when `MDEOUT=-1`, all
                                       !! values of `COEFF(*)` are set to zero.
                                       !!
                                       !! If the user is not satisfied with the fitted
                                       !! curve returned by [[DEFC]], the constrained
                                       !! least squares curve fitting subprogram [[DFC]]
                                       !! may be required.  The work done within [[DEFC]]
                                       !! to accumulate the data can be utilized by
                                       !! the user, if so desired.  This involves
                                       !! saving the first `(NBKPT-NORD+3)*(NORD+1)`
                                       !! entries of `W(*)` and providing this data
                                       !! to [[DFC]] with the "old problem" designation.
                                       !! The user should read the usage instructions
                                       !! for subprogram [[DFC]] for further details.
      integer,intent(in) :: Lw !! The amount of working storage actually
                               !! allocated for the working array `W(*)`.
                               !! This quantity is compared with the
                               !! actual amount of storage needed in [[DEFC]].
                               !! Insufficient storage allocated for `W(*)` is
                               !! an error.  This feature was included in [[DEFC]]
                               !! because misreading the storage formula
                               !! for `W(*)` might very well lead to subtle
                               !! and hard-to-find programming bugs.
                               !!
                               !! The length of the array `W(*)` must satisfy
                               !!```
                               !! LW >= (NBKPT-NORD+3)*(NORD+1)+
                               !!         (NBKPT+1)*(NORD+1)+
                               !!       2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
                               !!```
      real(wp) :: w(*) !! Working Array.
                       !! Its length is specified as an input parameter
                       !! in `LW` as noted above. The contents of `W(*)`
                       !! must not be modified by the user between calls
                       !! to [[DEFC]] with values of `MDEIN=1,2,2,...` .
                       !! The first `(NBKPT-NORD+3)*(NORD+1)` entries of
                       !! `W(*)` are acceptable as direct input to [[DFC]]
                       !! for an "old problem" only when `MDEOUT=1` or `2`.

      integer :: lbf, lbkpt, lg, lptemp, lww, lxtemp, mdg, mdw

      ! LWW=1               USAGE IN DEFCMN( ) OF W(*)..
      ! LWW,...,LG-1        W(*,*)
      ! LG,...,LXTEMP-1     G(*,*)
      ! LXTEMP,...,LPTEMP-1 XTEMP(*)
      ! LPTEMP,...,LBKPT-1  PTEMP(*)
      ! LBKPT,...,LBF       BKPT(*) (LOCAL TO DEFCMN( ))
      ! LBF,...,LBF+NORD**2 BF(*,*)

      mdg = Nbkpt + 1
      mdw = Nbkpt - Nord + 3
      lww = 1
      lg = lww + mdw*(Nord + 1)
      lxtemp = lg + mdg*(Nord + 1)
      lptemp = lxtemp + max(Ndata, Nbkpt)
      lbkpt = lptemp + max(Ndata, Nbkpt)
      lbf = lbkpt + Nbkpt

      call defcmn(Ndata, Xdata, Ydata, Sddata, Nord, Nbkpt, Bkpt, Mdein, Mdeout, &
                  Coeff, w(lbf), w(lxtemp), w(lptemp), w(lbkpt), w(lg), mdg, &
                  w(lww), mdw, Lw)

   end subroutine defc
!*****************************************************************************************

!*****************************************************************************************
!>
!  This is a companion subprogram to [[DEFC]].
!  This subprogram does weighted least squares fitting of data by
!  B-spline curves.
!  The documentation for [[DEFC]] has complete usage instructions.
!
!### Revision history
!  * 800801  DATE WRITTEN. Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and extensively revised (WRB & RWC)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900328  Added TYPE section.  (WRB)
!  * 900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!  * 900604  DP version created from SP version.  (RWC)

   subroutine defcmn(Ndata, Xdata, Ydata, Sddata, Nord, Nbkpt, Bkptin, &
                     Mdein, Mdeout, Coeff, Bf, Xtemp, Ptemp, Bkpt, g, Mdg, w, &
                     Mdw, Lw)

      integer :: Lw, Mdein, Mdeout, Mdg, Mdw, Nbkpt, Ndata, Nord
      real(wp) :: Bf(Nord, *), Bkpt(*), Bkptin(*), Coeff(*), &
                  g(Mdg, *), Ptemp(*), Sddata(*), w(Mdw, *), &
                  Xdata(*), Xtemp(*), Ydata(*)

      real(wp) :: rnorm, xmax, xmin, xval
      integer :: i, idata, ileft, intseq, ip, ir, irow, l, mt, n, &
                 nb, nordm1, nordp1, np1
      character(len=8) :: xern1, xern2
      integer :: dfspvn_j
      real(wp), dimension(20) :: dfspvn_deltam, dfspvn_deltap

      ! Initialize variables and analyze input.

      n = Nbkpt - Nord
      np1 = n + 1

      ! Initially set all output coefficients to zero.

      call dcopy(n, [0.0_wp], 0, Coeff, 1)

      Mdeout = -1
      if (Nord < 1 .or. Nord > 20) then
         write (*, *) 'IN DEFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.'
         return
      end if

      if (Nbkpt < 2*Nord) then
         write (*, *) 'IN DEFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.'
         return
      end if

      if (Ndata < 0) then
         write (*, *) 'IN DEFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.'
         return
      end if

      nb = (Nbkpt - Nord + 3)*(Nord + 1) + (Nbkpt + 1)*(Nord + 1) &
           + 2*max(Nbkpt, Ndata) + Nbkpt + Nord**2
      if (Lw < nb) then
         write (xern1, '(I8)') nb
         write (xern2, '(I8)') Lw
         write (*, *) 'IN DEFC, INSUFFICIENT STORAGE FOR W(*).  CHECK FORMULA '// &
            'THAT READS LW>= ... .  NEED = '//xern1// &
            ' GIVEN = '//xern2
         Mdeout = -1
         return
      end if

      if (Mdein /= 1 .and. Mdein /= 2) then
         write (*, *) 'IN DEFC, INPUT VALUE OF MDEIN MUST BE 1-2.'
         return
      end if

      ! Sort the breakpoints.

      call dcopy(Nbkpt, Bkptin, 1, Bkpt, 1)
      call dsort(Nbkpt, 1, Bkpt)

      ! Save interval containing knots.

      xmin = Bkpt(Nord)
      xmax = Bkpt(np1)
      nordm1 = Nord - 1
      nordp1 = Nord + 1

      ! Process least squares equations.

      ! Sort data and an array of pointers.

      call dcopy(Ndata, Xdata, 1, Xtemp, 1)
      do i = 1, Ndata
         Ptemp(i) = i
      end do
      ! JW : really Ptemp should be an integer array.
      !      it is real because they are stuffing it in
      !      a real work array and also using dsort on it.

      if (Ndata > 0) then
         call dsort(Ndata, 2, Xtemp, Ptemp)
         xmin = min(xmin, Xtemp(1))
         xmax = max(xmax, Xtemp(Ndata))
      end if

      ! Fix breakpoint array if needed. This should only involve very
      ! minor differences with the input array of breakpoints.

      do i = 1, Nord
         Bkpt(i) = min(Bkpt(i), xmin)
      end do

      do i = np1, Nbkpt
         Bkpt(i) = max(Bkpt(i), xmax)
      end do

      ! set up variables for dfspvn
      dfspvn_j = 1
      dfspvn_deltam = 0.0_wp
      dfspvn_deltap = 0.0_wp

      ! Initialize parameters of banded matrix processor, DBNDAC( ).

      mt = 0
      ip = 1
      ir = 1
      ileft = Nord
      intseq = 1
      do idata = 1, Ndata

         ! Sorted indices are in PTEMP(*).

         l = int(Ptemp(idata))
         xval = Xdata(l)

         ! When interval changes, process equations in the last block.

         if (xval >= Bkpt(ileft + 1)) then
            call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)
            mt = 0

            ! Move pointer up to have BKPT(ILEFT)<=XVAL, ILEFT<=N.

            do ileft = ileft, n
               if (xval < Bkpt(ileft + 1)) exit
               if (Mdein == 2) then
                  !  Data is being sequentially accumulated.
                  !  Transfer previously accumulated rows from W(*,*) to
                  !  G(*,*) and process them.
                  call dcopy(nordp1, w(intseq, 1), Mdw, g(ir, 1), Mdg)
                  call dbndac(g, Mdg, Nord, ip, ir, 1, intseq)
                  intseq = intseq + 1
               end if
            end do
         end if

         ! Obtain B-spline function value.

         call dfspvn(Bkpt, Nord, 1, xval, ileft, Bf, &
                     dfspvn_j, dfspvn_deltam, dfspvn_deltap)

         ! Move row into place.

         irow = ir + mt
         mt = mt + 1
         call dcopy(Nord, Bf, 1, g(irow, 1), Mdg)
         g(irow, nordp1) = Ydata(l)

         ! Scale data if uncertainty is nonzero.

         if (Sddata(l) /= 0.0_wp) call dscal(nordp1, 1.0_wp/Sddata(l), g(irow, 1), Mdg)

         ! When staging work area is exhausted, process rows.

         if (irow == Mdg - 1) then
            call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)
            mt = 0
         end if
      end do

      ! Process last block of equations.

      call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)

      ! Finish processing any previously accumulated rows from W(*,*)
      ! to G(*,*).

      if (Mdein == 2) then
         do i = intseq, np1
            call dcopy(nordp1, w(i, 1), Mdw, g(ir, 1), Mdg)
            call dbndac(g, Mdg, Nord, ip, ir, 1, min(n, i))
         end do
      end if

      ! Last call to adjust block positioning.

      call dcopy(nordp1, [0.0_wp], 0, g(ir, 1), Mdg)
      call dbndac(g, Mdg, Nord, ip, ir, 1, np1)

      ! Transfer accumulated rows from G(*,*) to W(*,*) for
      ! possible later sequential accumulation.

      do i = 1, np1
         call dcopy(nordp1, g(i, 1), Mdg, w(i, 1), Mdw)
      end do

     ! Solve for coefficients when possible.

      do i = 1, n
         if (g(i, 1) == 0.0_wp) then
            Mdeout = 2
            return
         end if
      end do

      ! All the diagonal terms in the accumulated triangular
      ! matrix are nonzero.  The solution can be computed but
      ! it may be unsuitable for further use due to poor
      ! conditioning or the lack of constraints.  No checking
      ! for either of these is done here.

      call dbndsl(1, g, Mdg, Nord, ip, ir, Coeff, n, rnorm)
      Mdeout = 1

   end subroutine defcmn
!*****************************************************************************************

!*****************************************************************************************
!>
!  These subroutines solve the least squares problem `Ax = b` for
!  banded matrices A using sequential accumulation of rows of the
!  data matrix.  Exactly one right-hand side vector is permitted.
!
!  These subroutines are intended for the type of least squares
!  systems that arise in applications such as curve or surface
!  fitting of data.  The least squares equations are accumulated and
!  processed using only part of the data.  This requires a certain
!  user interaction during the solution of Ax = b.
!
!  Specifically, suppose the data matrix (A B) is row partitioned
!  into Q submatrices.  Let (E F) be the T-th one of these
!  submatrices where E = (0 C 0).  Here the dimension of E is MT by N
!  and the dimension of C is MT by NB.  The value of NB is the
!  bandwidth of A.  The dimensions of the leading block of zeros in E
!  are MT by JT-1.
!
!  The user of the subroutine DBNDAC provides MT,JT,C and F for
!  T=1,...,Q.  Not all of this data must be supplied at once.
!
!  Following the processing of the various blocks (E F), the matrix
!  (A B) has been transformed to the form (R D) where R is upper
!  triangular and banded with bandwidth NB.  The least squares
!  system Rx = d is then easily solved using back substitution by
!  executing the statement CALL DBNDSL(1,...). The sequence of
!  values for JT must be nondecreasing.  This may require some
!  preliminary interchanges of rows and columns of the matrix A.
!
!  The primary reason for these subroutines is that the total
!  processing can take place in a working array of dimension MU by
!  NB+1.  An acceptable value for MU is
!
!                    MU = MAX(MT + N + 1),
!
!  where N is the number of unknowns.
!
!  Here the maximum is taken over all values of MT for T=1,...,Q.
!  Notice that MT can be taken to be a small as one, showing that
!  MU can be as small as N+2.  The subprogram DBNDAC processes the
!  rows more efficiently if MU is large enough so that each new
!  block (C F) has a distinct value of JT.
!
!  The four principle parts of these algorithms are obtained by the
!  following call statements:
!
!   * `CALL [[DBNDAC]](...)`  Introduce new blocks of data.
!   * `CALL [[DBNDSL]](1,...)` Compute solution vector and length of
!     residual vector.
!   * `CALL [[DBNDSL]](2,...)` Given any row vector H solve YR = H for the
!     row vector Y.
!   * `CALL [[DBNDSL]](3,...)` Given any column vector W solve RZ = W for
!     the column vector Z.
!
!### Remarks
!
!  To obtain the upper triangular matrix and transformed right-hand
!  side vector D so that the super diagonals of R form the columns
!  of G(*,*), execute the following Fortran statements.
!
!```fortran
!     nbp1=nb+1
!     do j=1, nbp1
!       g(ir,j) = 0.0
!     end do
!     mt=1
!     jt=n+1
!     call dbndac(g,mdg,nb,ip,ir,mt,jt)
!```
!
!### References
!
!  * C. L. Lawson and R. J. Hanson, Solving Least Squares
!    Problems, Prentice-Hall, Inc., 1974, Chapter 27.
!
!### Revision history
!  * 790101  DATE WRITTEN. Lawson, C. L., (JPL), Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dbndac(g, Mdg, Nb, Ip, Ir, Mt, Jt)

      implicit none

      integer,intent(in) :: Mdg !! The number of rows in the working array
                                !! `G(*,*)`.  The value of MDG should be `>= MU`.
                                !! The value of `MU` is defined in the abstract
                                !! of these subprograms.
      real(wp),intent(inout) :: g(Mdg, *) !! `G(MDG,NB+1)`
                                          !!
                                          !! *Input*
                                          !! The working array into which the user will
                                          !! place the `MT` by `NB+1` block `(C F)` in rows `IR`
                                          !! through `IR+MT-1`, columns 1 through `NB+1`.
                                          !! See descriptions of `IR` and `MT` below.
                                          !!
                                          !! *Output*
                                          !! The working array which will contain the
                                          !! processed rows of that part of the data
                                          !! matrix which has been passed to [[DBNDAC]].
      integer,intent(in) :: Nb !! The bandwidth of the data matrix `A`.
      integer,intent(inout) :: Ip !! *Input*
                                  !! Set by the user to the value 1 before the
                                  !! first call to [[DBNDAC]].  Its subsequent value
                                  !! is controlled by [[DBNDAC]] to set up for the
                                  !! next call to [[DBNDAC]].
                                  !!
                                  !! *Output*
                                  !! The value of this argument is advanced by
                                  !! [[DBNDAC]] to be ready for storing and processing
                                  !! a new block of data in `G(*,*)`.
      integer,intent(inout) :: Ir !! *Input*
                                  !! Index of the row of `G(*,*)` where the user is
                                  !! to place the new block of data `(C F)`.  Set by
                                  !! the user to the value 1 before the first call
                                  !! to [[DBNDAC]].  Its subsequent value is controlled
                                  !! by [[DBNDAC]]. A value of `IR > MDG` is considered
                                  !! an error.
                                  !!
                                  !! *Output*
                                  !! The value of this argument is advanced by
                                  !! [[DBNDAC]] to be ready for storing and processing
                                  !! a new block of data in `G(*,*)`.
      integer,intent(in) :: Mt !! Set by the user to indicate the
                               !! number of new rows of data in the block
      integer,intent(in) :: Jt !! Set by the user to indicate
                               !! the index of the first nonzero column in that
                               !! set of rows `(E F) = (0 C 0 F)` being processed.

      real(wp) :: rho
      integer :: i, ie, ig, ig1, ig2, iopt, j, jg, &
                 k, kh, l, lp1, mh, mu, nbp1, nerr

      real(wp), parameter :: zero = 0.0_wp

      ! ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.

      nbp1 = Nb + 1
      if (Mt <= 0 .or. Nb <= 0) return
      if (.not. Mdg < Ir) then
         if (Jt /= Ip) then
            if (Jt > Ir) then
               do i = 1, Mt
                  ig1 = Jt + Mt - i
                  ig2 = Ir + Mt - i
                  do j = 1, nbp1
                     g(ig1, j) = g(ig2, j)
                  end do
               end do
               ie = Jt - Ir
               do i = 1, ie
                  ig = Ir + i - 1
                  do j = 1, nbp1
                     g(ig, j) = zero
                  end do
               end do
               Ir = Jt
            end if
            mu = min(Nb - 1, Ir - Ip - 1)
            if (mu /= 0) then
               do l = 1, mu
                  k = min(l, Jt - Ip)
                  lp1 = l + 1
                  ig = Ip + l
                  do i = lp1, Nb
                     jg = i - k
                     g(ig, jg) = g(ig, i)
                  end do
                  do i = 1, k
                     jg = nbp1 - i
                     g(ig, jg) = zero
                  end do
               end do
            end if
            Ip = Jt
         end if
         mh = Ir + Mt - Ip
         kh = min(nbp1, mh)
         do i = 1, kh
            call dh12(1, i, max(i + 1, Ir - Ip + 1), mh, g(Ip, i), 1, &
                      rho, g(Ip, i + 1), 1, Mdg, nbp1 - i)
         end do
         Ir = Ip + kh
         if (kh >= nbp1) then
            do i = 1, Nb
               g(Ir - 1, i) = zero
            end do
         end if
      else
         nerr = 1
         iopt = 2
         write (*, *) 'MDG<IR, Probable error.'
      end if

   end subroutine dbndac
!*****************************************************************************************

!*****************************************************************************************
!>
!  These subroutines solve the least squares problem `Ax = b` for
!  banded matrices A using sequential accumulation of rows of the
!  data matrix.  Exactly one right-hand side vector is permitted.
!
!  See [[dbndac]] for a full description of how to use them.
!
!### References
!
!  * C. L. Lawson and R. J. Hanson, Solving Least Squares
!    Problems, Prentice-Hall, Inc., 1974, Chapter 27.
!
!### Revision history
!  * 790101  DATE WRITTEN. Lawson, C. L., (JPL), Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

    subroutine dbndsl(Mode, g, Mdg, Nb, Ip, Ir, x, n, Rnorm)

      integer,intent(in) :: Mode !! Set by the user to one of the values 1, 2, or
                                 !! 3. These values respectively indicate that
                                 !! the solution of `AX = B`, `YR = H` or `RZ = W` is
                                 !! required.
      integer,intent(in) :: Mdg !! The number of rows in the working array
                                !! `G(*,*)`.  The value of `MDG` should be `>= MU`.
                                !! The value of `MU` is defined in the abstract
                                !! of these subprograms.
                                !!
                                !! This argument has the same meaning and
                                !! contents as following the last call to [[DBNDAC]].
      real(wp),intent(in) :: g(Mdg, *) !! `G(MDG,NB+1)`
                                       !!
                                       !! This argument has the same meaning and
                                       !! contents as following the last call to [[DBNDAC]].
      integer,intent(in) :: Nb !! This argument has the same meaning and
                               !! contents as following the last call to [[DBNDAC]].
      integer,intent(in) :: Ip !! This argument has the same meaning and
                               !! contents as following the last call to [[DBNDAC]].
      integer,intent(in) :: Ir !! This argument has the same meaning and
                               !! contents as following the last call to [[DBNDAC]].
      real(wp),intent(inout) :: x(*) !! `X(N)`
                                     !!
                                     !! *Input* With mode=2 or 3 this array contains,
                                     !! respectively, the right-side vectors H or W of
                                     !! the systems YR = H or RZ = W.
                                     !!
                                     !! *Output* This array contains the solution vectors `X`,
                                     !! `Y` or `Z` of the systems `AX = B`, `YR = H` or
                                     !! `RZ = W` depending on the value of `MODE`=1,
                                     !! 2 or 3.
      integer,intent(in) :: n !! The number of variables in the solution
                              !! vector.  If any of the `N` diagonal terms are
                              !! zero the subroutine [[DBNDSL]] prints an
                              !! appropriate message.  This condition is
                              !! considered an error.
      real(wp),intent(out) :: Rnorm !! If `MODE=1`, `RNORM` is the Euclidean length of the
                                    !! residual vector `AX-B`.  When `MODE=2` or `3` RNORM`
                                    !! is set to zero.

      real(wp) :: rsq, s
      integer :: i, i1, i2, ie, ii, iopt, irm1, ix, j, &
                 jg, l, nerr, np1

      real(wp), parameter :: zero = 0.0_wp

      main: block

         Rnorm = zero
         select case (Mode)
         case (1)
            ! ALG. STEP 26
            do j = 1, n
               x(j) = g(j, Nb + 1)
            end do
            rsq = zero
            np1 = n + 1
            irm1 = Ir - 1
            if (np1 <= irm1) then
               do j = np1, irm1
                  rsq = rsq + g(j, Nb + 1)**2
               end do
               Rnorm = sqrt(rsq)
            end if
         case (2)
            do j = 1, n
               s = zero
               if (j /= 1) then
                  i1 = max(1, j - Nb + 1)
                  i2 = j - 1
                  do i = i1, i2
                     l = j - i + 1 + max(0, i - Ip)
                     s = s + x(i)*g(i, l)
                  end do
               end if
               l = max(0, j - Ip)
               if (g(j, l + 1) == 0) exit main
               x(j) = (x(j) - s)/g(j, l + 1)
            end do
            return
         end select

         ! MODE = 3
         do ii = 1, n
            i = n + 1 - ii
            s = zero
            l = max(0, i - Ip)
            if (i /= n) then
               ie = min(n + 1 - i, Nb)
               do j = 2, ie
                  jg = j + l
                  ix = i - 1 + j
                  s = s + g(i, jg)*x(ix)
               end do
            end if
            if (g(i, l + 1) == 0) exit main
            x(i) = (x(i) - s)/g(i, l + 1)
         end do

         return

      end block main

      ! error handling
      nerr = 1
      iopt = 2
      write (*, *) 'A zero diagonal term is in the n by n upper triangular matrix.'

   end subroutine dbndsl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculates the value of all possibly nonzero B-splines at `X` of
!  order `MAX(JHIGH,(J+1)(INDEX-1))` on `T`.
!
!### Revision history
!  * 780801  DATE WRITTEN
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)
!  * JW : made threadsafe. See also [[dbspvn]]

   subroutine dfspvn(t, Jhigh, Index, x, Ileft, Vnikx, j, deltam, deltap)

      real(wp),intent(in) :: t(*)
      integer,intent(in) :: Jhigh
      integer,intent(in) :: Index
      real(wp),intent(in) :: x
      integer,intent(in) :: Ileft
      real(wp) :: Vnikx(*)
      integer, intent(inout) :: j !! JW : added
      real(wp), dimension(20), intent(inout) :: deltam, deltap !! JW : added

      real(wp) :: vm, vmprev
      integer :: imjp1, ipj, jp1, jp1ml, l

      if (Index /= 2) then
         j = 1
         Vnikx(1) = 1.0_wp
         if (j >= Jhigh) return
      end if

      do
         ipj = Ileft + j
         deltap(j) = t(ipj) - x
         imjp1 = Ileft - j + 1
         deltam(j) = x - t(imjp1)
         vmprev = 0.0_wp
         jp1 = j + 1
         do l = 1, j
            jp1ml = jp1 - l
            vm = Vnikx(l)/(deltap(l) + deltam(jp1ml))
            Vnikx(l) = vm*deltap(l) + vmprev
            vmprev = vm*deltam(jp1ml)
         end do
         Vnikx(jp1) = vmprev
         j = jp1
         if (j >= Jhigh) exit
      end do

   end subroutine dfspvn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Construction and/or application of a single
!  Householder transformation. `Q = I + U*(U**T)/B`
!
!### Reference
!
!  * C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
!    to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
!
!### Revision history
!  * 790101  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)
!  * 900911  Added DDOT to real(wp) statement.  (WRB)

    subroutine dh12(Mode, Lpivot, l1, m, u, Iue, Up, c, Ice, Icv, Ncv)

    integer,intent(in) :: Mode !! 1 or 2   to select algorithm  H1  or  H2 .
    integer,intent(in) :: Lpivot !! the index of the pivot element.
    integer,intent(in) :: l1 !! If `L1 <= M` the transformation will be constructed to
                             !! zero elements indexed from `L1` through `M`. If `L1 > M`
                             !! the subroutine does an identity transformation.
    integer,intent(in) :: m !! see `l1`
    integer,intent(in) :: Iue !! the storage increment between elements of `U`.
    real(wp),intent(inout) :: u(Iue, *) !! On entry to H1 `U()` contains the pivot vector.
                                        !! On exit from H1 `U()` and `UP`
                                        !! contain quantities defining the vector `U` of the
                                        !! Householder transformation.   On entry to H2 `U()`
                                        !! and `UP` should contain quantities previously computed
                                        !! by H1.  These will not be modified by H2.
    real(wp),intent(inout)  :: Up !! see `u`
    real(wp),intent(inout) :: c(*) !! On entry to H1 or H2 `C()` contains a matrix which will be
                                   !! regarded as a set of vectors to which the Householder
                                   !! transformation is to be applied.  On exit `C()` contains the
                                   !! set of transformed vectors.
    integer,intent(in) :: Ice !! Storage increment between elements of vectors in `C()`.
    integer,intent(in) :: Icv !! Storage increment between vectors in `C()`.
    integer,intent(in) :: Ncv !! Number of vectors in `C()` to be transformed. If `NCV <= 0`
                              !! no operations will be done on `C()`.

      integer :: i, i2, i3, i4, incr, j, kl1, &
                 kl2, klp, l1m1, mml1p2
      real(wp) :: b, cl, clinv, ul1m1, sm

      real(wp), parameter :: one = 1.0_wp

      if (0 < Lpivot .and. Lpivot < l1 .and. l1 <= m) then
         cl = abs(u(1, Lpivot))
         if (Mode /= 2) then
            ! ****** CONSTRUCT THE TRANSFORMATION. ******
            do j = l1, m
               cl = max(abs(u(1, j)), cl)
            end do
            if (cl <= 0.0_wp) return
            clinv = one/cl
            sm = (u(1, Lpivot)*clinv)**2
            do j = l1, m
               sm = sm + (u(1, j)*clinv)**2
            end do
            cl = cl*sqrt(sm)
            if (u(1, Lpivot) > 0.0_wp) cl = -cl
            Up = u(1, Lpivot) - cl
            u(1, Lpivot) = cl
            ! ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
         elseif (cl <= 0.0_wp) then
            return
         end if

         if (Ncv > 0) then
            b = Up*u(1, Lpivot)
            ! B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
            if (b < 0.0_wp) then
               b = one/b
               mml1p2 = m - l1 + 2
               if (mml1p2 <= 20) then
                  i2 = 1 - Icv + Ice*(Lpivot - 1)
                  incr = Ice*(l1 - Lpivot)
                  do j = 1, Ncv
                     i2 = i2 + Icv
                     i3 = i2 + incr
                     i4 = i3
                     sm = c(i2)*Up
                     do i = l1, m
                        sm = sm + c(i3)*u(1, i)
                        i3 = i3 + Ice
                     end do
                     if (sm /= 0.0_wp) then
                        sm = sm*b
                        c(i2) = c(i2) + sm*Up
                        do i = l1, m
                           c(i4) = c(i4) + sm*u(1, i)
                           i4 = i4 + Ice
                        end do
                     end if
                  end do
               else
                  l1m1 = l1 - 1
                  kl1 = 1 + (l1m1 - 1)*Ice
                  kl2 = kl1
                  klp = 1 + (Lpivot - 1)*Ice
                  ul1m1 = u(1, l1m1)
                  u(1, l1m1) = Up
                  if (Lpivot /= l1m1) call dswap(Ncv, c(kl1), Icv, c(klp), Icv)
                  do j = 1, Ncv
                     sm = ddot(mml1p2, u(1, l1m1), Iue, c(kl1), Ice)
                     sm = sm*b
                     call daxpy(mml1p2, sm, u(1, l1m1), Iue, c(kl1), Ice)
                     kl1 = kl1 + Icv
                  end do
                  u(1, l1m1) = ul1m1
                  if (Lpivot /= l1m1) then
                     kl1 = kl2
                     call dswap(Ncv, c(kl1), Icv, c(klp), Icv)
                  end if
               end if
            end if
         end if
      end if

   end subroutine dh12
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sort an array and optionally make the same interchanges in
!  an auxiliary array.  The array may be sorted in increasing
!  or decreasing order.
!
!### History
!  * 29-dec-2022 : Replaced original routines.
!    Now just a wraper for [[sort_ascending]] recursive quicksort (JW)

   subroutine dsort(n, Kflag, Dx, Dy)
      implicit none

      integer, intent(in) :: n !! number of values in array DX to be sorted
      integer, intent(in) :: Kflag !! control parameter:
                                   !!  * Kflag < 0 : sort DX in decreasing order and optionally carry DY along.
                                   !!  * Kflag > 0 : sort DX in increasing order and optionally carry DY along.
      real(wp), dimension(*), intent(inout) :: Dx !! array of values to be sorted   (usually abscissas)
      real(wp), dimension(*), intent(inout), optional :: Dy !! array to be (optionally) carried along

      if (n < 1) then
         write (*, *) 'The number of values to be sorted is not positive.'
         return
      end if

      if (abs(Kflag) == 0) then
         write (*, *) 'The sort control parameter, K, cannot be 0.'
         return
      end if

      ! Alter array DX to get decreasing order if needed
      if (Kflag < 0) Dx(1:n) = -Dx(1:n)
      call sort_ascending(n, Dx, Dy)
      if (Kflag < 0) Dx(1:n) = -Dx(1:n)

   end subroutine dsort
!*****************************************************************************************

!*****************************************************************************************
!>
!  Recursive quicksoft.
!  Modified to also carry along a second array.
!
!### Author
!  * Jacob Williams

   subroutine sort_ascending(n, dx, dy)

      integer, intent(in) :: n
      real(wp), dimension(*), intent(inout) :: dx !! array of values to be sorted
      real(wp), dimension(*), intent(inout), optional :: dy !! array to be (optionally) carried along

      logical :: carry_dy !! if `dy` is to be also sorted

      integer, parameter :: max_size_for_insertion_sort = 20 !! max size for using insertion sort.
                                                             !! (otherwise, use quicksort)


      carry_dy = present(dy)
      call quicksort(1, n)

   contains

      recursive subroutine quicksort(ilow, ihigh)

         !! Sort the array (ascending order).

         integer, intent(in) :: ilow
         integer, intent(in) :: ihigh

         integer :: ipivot !! pivot element
         integer :: i      !! counter
         integer :: j      !! counter

         if (ihigh - ilow <= max_size_for_insertion_sort .and. ihigh > ilow) then

            ! do insertion sort:
            do i = ilow + 1, ihigh
               do j = i, ilow + 1, -1
                  if (dx(j) < dx(j - 1)) then
                     call swap(dx(j), dx(j - 1))
                     if (carry_dy) call swap(dy(j), dy(j - 1))
                  else
                     exit
                  end if
               end do
            end do

         else if (ihigh - ilow > max_size_for_insertion_sort) then

            ! do the normal quicksort:
            call partition(ilow, ihigh, ipivot)
            call quicksort(ilow, ipivot - 1)
            call quicksort(ipivot + 1, ihigh)

         end if

      end subroutine quicksort

      subroutine partition(ilow, ihigh, ipivot)

         !! Partition the array

         integer, intent(in)  :: ilow
         integer, intent(in)  :: ihigh
         integer, intent(out) :: ipivot

         integer :: i, ip, im

         im = (ilow + ihigh)/2
         call swap(dx(ilow), dx(im))
         if (carry_dy) call swap(dy(ilow), dy(im))
         ip = ilow
         do i = ilow + 1, ihigh
            if (dx(i) < dx(ilow)) then
               ip = ip + 1
               call swap(dx(ip), dx(i))
               if (carry_dy) call swap(dy(ip), dy(i))
            end if
         end do
         call swap(dx(ilow), dx(ip))
         if (carry_dy) call swap(dy(ilow), dy(ip))
         ipivot = ip

      end subroutine partition

      subroutine swap(v1, v2)
         !! swap two real values
         real(wp), intent(inout) :: v1
         real(wp), intent(inout) :: v2
         real(wp) :: tmp
         tmp = v1
         v1 = v2
         v2 = tmp
      end subroutine swap

   end subroutine sort_ascending
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subprogram fits a piecewise polynomial curve
!  to discrete data.  The piecewise polynomials are
!  represented as B-splines.
!  The fitting is done in a weighted least squares sense.
!  Equality and inequality constraints can be imposed on the
!  fitted curve.
!
!### Evaluating the Variance Function
!
!  To evaluate the variance function (assuming
!  that the uncertainties of the Y values were
!  provided to [[DFC]] and an input value of
!  MODE=2 or 4 was used), use the function
!  subprogram [[DCV]]
!
!```fortran
!    var = dcv(xval,ndata,nconst,nord,nbkpt, bkpt,w)
!```
!
!  Here XVAL is the point where the variance is
!  desired.  The other arguments have the same
!  meaning as in the usage of [[DFC]].
!
!  For those users employing the old problem
!  designation, let MDATA be the number of data
!  points in the problem.  (This may be different
!  from NDATA if the old problem designation
!  feature was used.)  The value, VAR, should be
!  multiplied by the quantity
!
!  `DBLE(MAX(NDATA-N,1))/DBLE(MAX(MDATA-N,1))`
!
!  The output of this subprogram is not defined
!  if an input value of MODE=1 or 3 was used in
!  FC( ) or if an output value of MODE=-1, 2, or
!  3 was obtained.  The variance function, except
!  for the scaling factor noted above, is given
!  by
!
!  `VAR=(transpose of B(XVAL))*C*B(XVAL)`
!
!  The vector B(XVAL) is the B-spline basis
!  function values at X=XVAL.
!  The covariance matrix, C, of the solution
!  coefficients accounts only for the least
!  squares equations and the explicitly stated
!  equality constraints.  This fact must be
!  considered when interpreting the variance
!  function from a data fitting problem that has
!  inequality constraints on the fitted curve.
!
!### Evaluating the Fitted Curve
!
!  * Refer to the [[defc]] header
!
!### Revision history
!  * 780801  DATE WRITTEN. Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900510  Convert references to XERRWV to references to XERMSG.  (RWC)
!  * 900607  Editorial changes to Prologue to make Prologues for EFC,
!    DEFC, FC, and DFC look as much the same as possible.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dfc (ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, &
                   nconst, xconst, yconst, nderiv, mode, coeff, w, iw)

   integer, intent(in) :: ndata !! number of points (size of `xdata` and `ydata`).
                                !! Any non-negative value of `NDATA` is allowed.
                                !! A negative value of `NDATA` is an error.
   real(wp), intent(in) :: xdata(*) !! X data array. No sorting of `XDATA(*)` is required.
   real(wp), intent(in) :: ydata(*) !! Y data array.
   real(wp), intent(in) :: sddata(*) !! Y value standard deviation or uncertainty.
                                     !! A zero value for any entry of
                                     !! `SDDATA(*)` will weight that data point as 1.
                                     !! Otherwise the weight of that data point is
                                     !! the reciprocal of this entry.
   integer, intent(in) :: nord !! B-spline order.
                               !! (The order of the spline is one more than the
                               !! degree of the piecewise polynomial defined on
                               !! each interval.  This is consistent with the
                               !! B-spline package convention.  For example,
                               !! `NORD=4` when we are using piecewise cubics.)
                               !! `NORD` must be in the range `1 <= NORD <= 20`.
   integer,intent(in) :: Nbkpt !! The value of `NBKPT` must satisfy the condition `NBKPT >= 2*NORD`.
   real(wp),dimension(*),intent(in) :: Bkpt !! `NBKPT` knots of the B-spline.
                                            !! Normally the
                                            !! problem data interval will be included between
                                            !! the limits `BKPT(NORD)` and `BKPT(NBKPT-NORD+1)`.
                                            !! The additional end knots `BKPT(I),I=1,...,NORD-1`
                                            !! and `I=NBKPT-NORD+2,...,NBKPT`, are
                                            !! required to compute the functions used to fit
                                            !! the data.  No sorting of `BKPT(*)` is required.
                                            !! Internal to [[DEFC]] the extreme end knots may
                                            !! be reduced and increased respectively to
                                            !! accommodate any data values that are exterior
                                            !! to the given knot values.  The contents of
                                            !! `BKPT(*)` is not changed.
   integer, intent(in) :: nconst !! The number of conditions that constrain the
                                 !! B-spline is NCONST.  A constraint is specified
                                 !! by an (X,Y) pair in the arrays XCONST(*) and
                                 !! YCONST(*), and by the type of constraint and
                                 !! derivative value encoded in the array
                                 !! NDERIV(*).
   real(wp), intent(in) :: xconst(*) !! X value of constraint.
                                     !! No sorting of XCONST(*) is required.
   real(wp), intent(in) :: yconst(*) !! Y value of constraint
   integer, intent(in) :: nderiv(*) !! The value of NDERIV(*) is
                                    !! determined as follows.  Suppose the I-th
                                    !! constraint applies to the J-th derivative
                                    !! of the B-spline.  (Any non-negative value of
                                    !! J < NORD is permitted.  In particular the
                                    !! value J=0 refers to the B-spline itself.)
                                    !! For this I-th constraint, set
                                    !!```
                                    !!  XCONST(I)=X,
                                    !!  YCONST(I)=Y, and
                                    !!  NDERIV(I)=ITYPE+4*J, where
                                    !!
                                    !!  ITYPE = 0,      if (J-th deriv. at X) <= Y.
                                    !!        = 1,      if (J-th deriv. at X) >= Y.
                                    !!        = 2,      if (J-th deriv. at X) == Y.
                                    !!        = 3,      if (J-th deriv. at X) ==
                                    !!                     (J-th deriv. at Y).
                                    !!```
                                    !! (A value of NDERIV(I)=-1 will cause this
                                    !! constraint to be ignored.  This subprogram
                                    !! feature is often useful when temporarily
                                    !! suppressing a constraint while still
                                    !! retaining the source code of the calling
                                    !! program.)
   integer, intent(inout) :: mode !! *Input*
                                  !!
                                  !! An input flag that directs the least squares
                                  !! solution method used by [[DFC]].
                                  !!
                                  !! The variance function, referred to below,
                                  !! defines the square of the probable error of
                                  !! the fitted curve at any point, XVAL.
                                  !! This feature of [[DFC]] allows one to use the
                                  !! square root of this variance function to
                                  !! determine a probable error band around the
                                  !! fitted curve.
                                  !!
                                  !!  * `=1`  a new problem.  No variance function.
                                  !!  * `=2`  a new problem.  Want variance function.
                                  !!  * `=3`  an old problem.  No variance function.
                                  !!  * `=4`  an old problem.  Want variance function.
                                  !!
                                  !! Any value of MODE other than 1-4 is an error.
                                  !!
                                  !! The user with a new problem can skip directly
                                  !! to the description of the input parameters
                                  !! IW(1), IW(2).
                                  !!
                                  !! If the user correctly specifies the new or old
                                  !! problem status, the subprogram [[DFC]] will
                                  !! perform more efficiently.
                                  !! By an old problem it is meant that subprogram
                                  !! [[DFC]] was last called with this same set of
                                  !! knots, data points and weights.
                                  !!
                                  !! Another often useful deployment of this old
                                  !! problem designation can occur when one has
                                  !! previously obtained a Q-R orthogonal
                                  !! decomposition of the matrix resulting from
                                  !! B-spline fitting of data (without constraints)
                                  !! at the breakpoints BKPT(I), I=1,...,NBKPT.
                                  !! For example, this matrix could be the result
                                  !! of sequential accumulation of the least
                                  !! squares equations for a very large data set.
                                  !! The user writes this code in a manner
                                  !! convenient for the application.  For the
                                  !! discussion here let
                                  !!
                                  !! `N=NBKPT-NORD, and K=N+3`
                                  !!
                                  !! Let us assume that an equivalent least squares
                                  !! system
                                  !!
                                  !! `RC=D`
                                  !!
                                  !! has been obtained.  Here R is an N+1 by N
                                  !! matrix and D is a vector with N+1 components.
                                  !! The last row of R is zero.  The matrix R is
                                  !! upper triangular and banded.  At most NORD of
                                  !! the diagonals are nonzero.
                                  !! The contents of R and D can be copied to the
                                  !! working array W(*) as follows.
                                  !!
                                  !! The I-th diagonal of R, which has N-I+1
                                  !! elements, is copied to W(*) starting at
                                  !!
                                  !! `W((I-1)*K+1),`
                                  !!
                                  !! for I=1,...,NORD.
                                  !! The vector D is copied to W(*) starting at
                                  !!
                                  !! `W(NORD*K+1)`
                                  !!
                                  !! The input value used for NDATA is arbitrary
                                  !! when an old problem is designated.  Because
                                  !! of the feature of [[DFC]] that checks the
                                  !! working storage array lengths, a value not
                                  !! exceeding NBKPT should be used.  For example,
                                  !! use NDATA=0.
                                  !!
                                  !! (The constraints or variance function request
                                  !! can change in each call to [[DFC]].)  A new
                                  !! problem is anything other than an old problem.
                                  !!
                                  !! *Output*
                                  !!
                                  !! An output flag that indicates the status
                                  !! of the constrained curve fit.
                                  !!
                                  !!  * `=-1`  a usage error of [[DFC]] occurred.  The
                                  !!    offending condition is noted with the
                                  !!    SLATEC library error processor, XERMSG.
                                  !!    In case the working arrays W(*) or IW(*)
                                  !!    are not long enough, the minimal
                                  !!    acceptable length is printed.
                                  !! * `= 0`  successful constrained curve fit.
                                  !! * `= 1`  the requested equality constraints
                                  !!   are contradictory.
                                  !! * `= 2`  the requested inequality constraints
                                  !!    are contradictory.
                                  !! * `= 3`  both equality and inequality constraints
                                  !!   are contradictory.
   real(wp), intent(out) :: coeff(*) !! If the output value of MODE=0 or 1, this array
                                     !! contains the unknowns obtained from the least
                                     !! squares fitting process.  These N=NBKPT-NORD
                                     !! parameters are the B-spline coefficients.
                                     !! For MODE=1, the equality constraints are
                                     !! contradictory.  To make the fitting process
                                     !! more robust, the equality constraints are
                                     !! satisfied in a least squares sense.  In this
                                     !! case the array COEFF(*) contains B-spline
                                     !! coefficients for this extended concept of a
                                     !! solution.  If MODE=-1,2 or 3 on output, the
                                     !! array COEFF(*) is undefined.
   real(wp) :: w(*) !! real work array of length `IW(1)`. The
                    !! contents of `W(*)` must not be modified by the
                    !! user if the variance function is desired.
                    !!
                    !! The length of W(*) must be at least
                    !!```
                    !!   NB=(NBKPT-NORD+3)*(NORD+1)+
                    !!       2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
                    !!```
                    !! Whenever possible the code uses banded matrix
                    !! processors DBNDAC( ) and DBNDSL( ).  These
                    !! are utilized if there are no constraints,
                    !! no variance function is required, and there
                    !! is sufficient data to uniquely determine the
                    !! B-spline coefficients.  If the band processors
                    !! cannot be used to determine the solution,
                    !! then the constrained least squares code DLSEI
                    !! is used.  In this case the subprogram requires
                    !! an additional block of storage in W(*).  For
                    !! the discussion here define the integers NEQCON
                    !! and NINCON respectively as the number of
                    !! equality (ITYPE=2,3) and inequality
                    !! (ITYPE=0,1) constraints imposed on the fitted
                    !! curve.  Define
                    !!
                    !! `L = NBKPT-NORD+1`
                    !!
                    !! and note that
                    !!
                    !! `NCONST = NEQCON+NINCON`
                    !!
                    !! When the subprogram [[DFC]] uses [[DLSEI]] the
                    !! length of the working array W(*) must be at
                    !! least
                    !!
                    !! `LW = NB+(L+NCONST)*L+2*(NEQCON+L)+(NINCON+L)+(NINCON+2)*(L+6)`
   integer  :: iw(*) !! integer work array of length `IW(2)`
                     !!
                     !! `IW(1),IW(2)` are the amounts of working storage actually
                     !! allocated for the working arrays W(*) and
                     !! IW(*).  These quantities are compared with the
                     !! actual amounts of storage needed in [[DFC]].
                     !! Insufficient storage allocated for either
                     !! W(*) or IW(*) is an error.  This feature was
                     !! included in [[DFC]] because misreading the
                     !! storage formulas for W(*) and IW(*) might very
                     !! well lead to subtle and hard-to-find
                     !! programming bugs.
                     !!
                     !! The length of the array IW(*) must be at least
                     !!
                     !! `IW1 = NINCON+2*L`
                     !!
                     !! in any case.

   integer :: i1, i2, i3, i4, i5, i6, i7, mdg, mdw

   mdg = nbkpt - nord + 3
   mdw = nbkpt - nord + 1 + nconst

   ! USAGE IN DFCMN( ) OF W(*)..
   !     I1,...,I2-1      G(*,*)
   !     I2,...,I3-1      XTEMP(*)
   !     I3,...,I4-1      PTEMP(*)
   !     I4,...,I5-1      BKPT(*) (LOCAL TO [[DFCMN]])
   !     I5,...,I6-1      BF(*,*)
   !     I6,...,I7-1      W(*,*)
   !     I7,...           WORK(*) FOR [[DLSEI]]

   i1 = 1
   i2 = i1 + mdg*(nord+1)
   i3 = i2 + max(ndata,nbkpt)
   i4 = i3 + max(ndata,nbkpt)
   i5 = i4 + nbkpt
   i6 = i5 + nord*nord
   i7 = i6 + mdw*(nbkpt-nord+1)
   call dfcmn(ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, nconst, &
              xconst, yconst, nderiv, mode, coeff, w(i5), w(i2), w(i3), &
              w(i4), w(i1), mdg, w(i6), mdw, w(i7), iw)

   end subroutine dfc
!*****************************************************************************************

!*****************************************************************************************
!>
!  This is a companion subprogram to [[DFC]].
!  The documentation for [[DFC]] has complete usage instructions.
!
!### Revision history
!  * 780801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and extensively revised (WRB & RWC)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900328  Added TYPE section.  (WRB)
!  * 900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!  * 900604  DP version created from SP version.  (RWC)

   subroutine dfcmn (ndata, xdata, ydata, sddata, nord, nbkpt, &
                     bkptin, nconst, xconst, yconst, nderiv, mode, coeff, bf, xtemp, &
                     ptemp, bkpt, g, mdg, w, mdw, work, iwork)

   integer :: iwork(*), mdg, mdw, mode, nbkpt, nconst, ndata, nderiv(*), &
              nord
   real(wp) :: bf(nord,*), bkpt(*), bkptin(*), coeff(*), &
               g(mdg,*), ptemp(*), sddata(*), w(mdw,*), work(*), &
               xconst(*), xdata(*), xtemp(*), yconst(*), ydata(*)

   real(wp) :: prgopt(10), rnorm, rnorme, rnorml, xmax, &
               xmin, xval, yval
   integer :: i, idata, ideriv, ileft, intrvl, intw1, ip, ir, irow, &
              itype, iw1, iw2, l, lw, mt, n, nb, neqcon, nincon, nordm1, &
              nordp1, np1
   logical :: band, new, var
   character(len=8) :: xern1
   integer :: dfspvn_j
   real(wp), dimension(20) :: dfspvn_deltam, dfspvn_deltap

   ! Analyze input.

   if (nord<1 .or. nord>20) then
      write(*,*) 'IN DFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.'
      mode = -1
      return
   elseif (nbkpt<2*nord) then
      write(*,*) 'IN DFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' // &
                 'THE B-SPLINE ORDER.'
      mode = -1
      return
   endif
   if (ndata<0) then
      write(*,*) 'IN DFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.'
      mode = -1
      return
   endif

   ! Amount of storage allocated for W(*), IW(*).

   iw1 = iwork(1)
   iw2 = iwork(2)
   nb = (nbkpt-nord+3)*(nord+1) + 2*max(ndata,nbkpt) + nbkpt + &
        nord**2

   ! See if sufficient storage has been allocated.

   if (iw1<nb) then
      write (xern1, '(I8)') nb
      write(*,*) 'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = ' // xern1
      mode = -1
      return
   endif

   select case (mode)
   case (1)
      band = .true.
      var = .false.
      new = .true.
   case (2)
      band = .false.
      var = .true.
      new = .true.
   case (3)
      band = .true.
      var = .false.
      new = .false.
   case (4)
      band = .false.
      var = .true.
      new = .false.
   case default
      write(*,*) 'IN DFC, INPUT VALUE OF MODE MUST BE 1-4.'
      mode = -1
      return
   end select
   mode = 0

   ! Sort the breakpoints.

   call dcopy (nbkpt, bkptin, 1, bkpt, 1)
   call dsort (nbkpt, 1, bkpt)

   ! Initialize variables.

   neqcon = 0
   nincon = 0
   do i = 1,nconst
      l = nderiv(i)
      itype = mod(l,4)
      if (itype<2) then
         nincon = nincon + 1
      else
         neqcon = neqcon + 1
      endif
   end do
   ! set up variables for dfspvn
   dfspvn_j = 1
   dfspvn_deltam = 0.0_wp
   dfspvn_deltap = 0.0_wp

   ! Compute the number of variables.

   n = nbkpt - nord
   np1 = n + 1
   lw = nb + (np1+nconst)*np1 + 2*(neqcon+np1) + (nincon+np1) + &
        (nincon+2)*(np1+6)
   intw1 = nincon + 2*np1

   ! Save interval containing knots.

   xmin = bkpt(nord)
   xmax = bkpt(np1)

   ! Find the smallest referenced independent variable value in any
   ! constraint.

   do i = 1,nconst
      xmin = min(xmin,xconst(i))
      xmax = max(xmax,xconst(i))
   end do
   nordm1 = nord - 1
   nordp1 = nord + 1

   ! Define the option vector PRGOPT(1-10) for use in [[DLSEI]].

   prgopt(1) = 4

   ! Set the covariance matrix computation flag.

   prgopt(2) = 1
   if (var) then
      prgopt(3) = 1
   else
      prgopt(3) = 0
   endif

   ! Increase the rank determination tolerances for both equality
   ! constraint equations and least squares equations.

   prgopt(4) = 7
   prgopt(5) = 4
   prgopt(6) = 1.0e-4_wp
   prgopt(7) = 10
   prgopt(8) = 5
   prgopt(9) = 1.0e-4_wp
   prgopt(10) = 1

   ! Turn off work array length checking in [[DLSEI]].

   iwork(1) = 0
   iwork(2) = 0

   ! Initialize variables and analyze input.

   if (new) then

      ! To process least squares equations sort data and an array of
      ! pointers.

      call dcopy (ndata, xdata, 1, xtemp, 1)
      do i = 1,ndata
         ptemp(i) = i
      end do

      if (ndata>0) then
         call dsort (ndata, 2, xtemp, ptemp)
         xmin = min(xmin,xtemp(1))
         xmax = max(xmax,xtemp(ndata))
      endif

      ! Fix breakpoint array if needed.

      do i = 1,nord
         bkpt(i) = min(bkpt(i),xmin)
      end do

      do i = np1,nbkpt
         bkpt(i) = max(bkpt(i),xmax)
      end do

      ! Initialize parameters of banded matrix processor, DBNDAC( ).

      mt = 0
      ip = 1
      ir = 1
      ileft = nord
      do idata = 1,ndata

         ! Sorted indices are in PTEMP(*).

         l = ptemp(idata)
         xval = xdata(l)

         ! When interval changes, process equations in the last block.

         if (xval>=bkpt(ileft+1)) then
            call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)
            mt = 0

            ! Move pointer up to have BKPT(ILEFT)<=XVAL,
            ! ILEFT<NP1.

            do
               if (xval>=bkpt(ileft+1) .and. ileft<n) then
                  ileft = ileft + 1
               else
                  exit
               endif
            end do

         endif

         ! Obtain B-spline function value.

         call dfspvn (bkpt, nord, 1, xval, ileft, bf, &
                      dfspvn_j, dfspvn_deltam, dfspvn_deltap)

         ! Move row into place.

         irow = ir + mt
         mt = mt + 1
         call dcopy (nord, bf, 1, g(irow,1), mdg)
         g(irow,nordp1) = ydata(l)

         ! Scale data if uncertainty is nonzero.

         if (sddata(l)/=0.0_wp) call dscal (nordp1, 1.0_wp/sddata(l), &
                                            g(irow,1), mdg)

         ! When staging work area is exhausted, process rows.

         if (irow==mdg-1) then
            call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)
            mt = 0
         endif
      end do

      ! Process last block of equations.

      call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)

      ! Last call to adjust block positioning.

      call dcopy (nordp1, [0.0_wp], 0, g(ir,1), mdg)
      call dbndac (g, mdg, nord, ip, ir, 1, np1)
   endif

   band = band .and. nconst==0
   do i = 1,n
      band = band .and. g(i,1)/=0.0_wp
   end do

   ! Process banded least squares equations.

   if (band) then
      call dbndsl (1, g, mdg, nord, ip, ir, coeff, n, rnorm)
      return
   endif

   ! Check further for sufficient storage in working arrays.

   if (iw1<lw) then
      write (xern1, '(I8)') lw
      write(*,*) 'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = ' // xern1
      mode = -1
      return
   endif

   if (iw2<intw1) then
      write (xern1, '(I8)') intw1
      write(*,*) 'IN DFC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = ' // xern1
      mode = -1
      return
   endif

   ! Write equality constraints.
   ! Analyze constraint indicators for an equality constraint.

   neqcon = 0
   do idata = 1,nconst
      l = nderiv(idata)
      itype = mod(l,4)
      if (itype>1) then
         ideriv = l/4
         neqcon = neqcon + 1
         ileft = nord
         xval = xconst(idata)

         do
            if (xval<bkpt(ileft+1) .or. ileft>=n) exit
            ileft = ileft + 1
         end do

         call dfspvd (bkpt, nord, xval, ileft, bf, ideriv+1)
         call dcopy (np1, [0.0_wp], 0, w(neqcon,1), mdw)
         call dcopy (nord, bf(1,ideriv+1), 1, w(neqcon,ileft-nordm1), &
                     mdw)

         if (itype==2) then
            w(neqcon,np1) = yconst(idata)
         else
            ileft = nord
            yval = yconst(idata)

            do
               if (yval<bkpt(ileft+1) .or. ileft>=n) exit
               ileft = ileft + 1
            end do

            call dfspvd (bkpt, nord, yval, ileft, bf, ideriv+1)
            call daxpy (nord, -1.0_wp, bf(1, ideriv+1), 1, &
                        w(neqcon, ileft-nordm1), mdw)
         endif
      endif
   end do

   ! Transfer least squares data.

   do i = 1,np1
      irow = i + neqcon
      call dcopy (n, [0.0_wp], 0, w(irow,1), mdw)
      call dcopy (min(np1-i, nord), g(i,1), mdg, w(irow,i), mdw)
      w(irow,np1) = g(i,nordp1)
   end do

   ! Write inequality constraints.
   ! Analyze constraint indicators for inequality constraints.

   nincon = 0
   do idata = 1,nconst
      l = nderiv(idata)
      itype = mod(l,4)
      if (itype<2) then
         ideriv = l/4
         nincon = nincon + 1
         ileft = nord
         xval = xconst(idata)

         do
            if (xval<bkpt(ileft+1) .or. ileft>=n) exit
            ileft = ileft + 1
         end do

         call dfspvd (bkpt, nord, xval, ileft, bf, ideriv+1)
         irow = neqcon + np1 + nincon
         call dcopy (n, [0.0_wp], 0, w(irow,1), mdw)
         intrvl = ileft - nordm1
         call dcopy (nord, bf(1, ideriv+1), 1, w(irow, intrvl), mdw)

         if (itype==1) then
            w(irow,np1) = yconst(idata)
         else
            w(irow,np1) = -yconst(idata)
            call dscal (nord, -1.0_wp, w(irow, intrvl), mdw)
         endif
      endif
   end do

   ! Solve constrained least squares equations.

   call dlsei(w, mdw, neqcon, np1, nincon, n, prgopt, coeff, rnorme, &
              rnorml, mode, work, iwork)

   end subroutine dfcmn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculates value and derivs of all B-splines which do not vanish at `X`
!
!  Fill `VNIKX(J,IDERIV), J=IDERIV, ... ,K`  with nonzero values of
!  B-splines of order `K+1-IDERIV , IDERIV=NDERIV, ... ,1`, by repeated
!  calls to [[DFSPVN]]
!
!### Revision history
!  * 780801  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 890911  Removed unnecessary intrinsics.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)

   subroutine dfspvd (t, k, x, ileft, vnikx, nderiv)

   real(wp) :: t(*)
   integer :: k
   real(wp) :: x
   integer :: ileft
   real(wp) :: vnikx(k,*)
   integer :: nderiv

   real(wp) :: a(20,20)
   integer :: ideriv,idervm,i,j,kmd,m,jm1,ipkmd,l,jlow
   real(wp) :: fkmd,diff,v
   integer :: dfspvn_j
   real(wp), dimension(20) :: dfspvn_deltam, dfspvn_deltap

   ! set up variables for dfspvn
   dfspvn_j = 1
   dfspvn_deltam = 0.0_wp
   dfspvn_deltap = 0.0_wp

   call dfspvn(t,k+1-nderiv,1,x,ileft,vnikx(nderiv,nderiv),&
               dfspvn_j,dfspvn_deltam,dfspvn_deltap)
   if (nderiv <= 1) return

   ideriv = nderiv
   do i=2,nderiv
      idervm = ideriv-1
      do j=ideriv,k
          vnikx(j-1,idervm) = vnikx(j,ideriv)
      end do
      ideriv = idervm
      call dfspvn(t,0,2,x,ileft,vnikx(ideriv,ideriv),&
                  dfspvn_j,dfspvn_deltam,dfspvn_deltap)
   end do

   do i=1,k
      do j=1,k
         a(i,j) = 0.0_wp
      end do
      a(i,i) = 1.0_wp
   end do
   kmd = k
   do m=2,nderiv
      kmd = kmd-1
      fkmd = kmd
      i = ileft
      j = k
      do
         jm1 = j-1
         ipkmd = i + kmd
         diff = t(ipkmd) - t(i)
         if (jm1 == 0) exit
         if (diff /= 0.0_wp) then
            do l=1,j
               a(l,j) = (a(l,j) - a(l,j-1))/diff*fkmd
            end do
         end if
         j = jm1
         i = i - 1
      end do
      if (diff /= 0.0_wp) then
         a(1,1) = a(1,1)/diff*fkmd
      end if
      do i=1,k
         v = 0.0_wp
         jlow = max(i,m)
         do j=jlow,k
            v = a(i,j)*vnikx(j,m) + v
         end do
         vnikx(i,m) = v
      end do
   end do

   end subroutine dfspvd
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve a least squares problem for banded matrices using
!  sequential accumulation of rows of the data matrix.
!  Exactly one right-hand side vector is permitted.
!
!  This subroutine solves a linear least squares problem or a set of
!  linear least squares problems having the same matrix but different
!  right-side vectors.  The problem data consists of an M by N matrix
!  A, an M by NB matrix B, and an absolute tolerance parameter TAU
!  whose usage is described below.  The NB column vectors of B
!  represent right-side vectors for NB distinct linear least squares
!  problems.
!
!  This set of problems can also be written as the matrix least
!  squares problem
!
!  `A = B`,
!
!  where X is the N by NB solution matrix.
!
!  Note that if B is the M by M identity matrix, then X will be the
!  pseudo-inverse of A.
!
!  This subroutine first transforms the augmented matrix (A B) to a
!  matrix (R C) using premultiplying Householder transformations with
!  column interchanges.  All subdiagonal elements in the matrix R are
!  zero and its diagonal elements satisfy
!
!```
!  abs(r(i,i))>=abs(r(i+1,i+1)),
!  i = 1,...,l-1, where
!  l = min(m,n).
!```
!
!  The subroutine will compute an integer, KRANK, equal to the number
!  of diagonal terms of R that exceed TAU in magnitude. Then a
!  solution of minimum Euclidean length is computed using the first
!  KRANK rows of (R C).
!
!  To be specific we suggest that the user consider an easily
!  computable matrix norm, such as, the maximum of all column sums of
!  magnitudes.
!
!  Now if the relative uncertainty of B is EPS, (norm of uncertainty/
!  norm of B), it is suggested that TAU be set approximately equal to
!  EPS*(norm of A).
!
!### References
!  * C. L. Lawson and R. J. Hanson, Solving Least Squares
!    Problems, Prentice-Hall, Inc., 1974, Chapter 14.
!
!### Revision history
!  * 790101  DATE WRITTEN. Lawson, C. L., (JPL), Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dhfti (a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h, g, ip)

   integer,intent(in) :: mda !! actual leading dimension of `a`
   integer,intent(in) :: mdb !! actual leading dimension of `b`
   real(wp),intent(inout) :: a(mda,*) !! `A(MDA,N)`.
                                      !! The array A(*,*) initially contains the M by N
                                      !! matrix A of the least squares problem AX = B.
                                      !! The first dimensioning parameter of the array
                                      !! A(*,*) is MDA, which must satisfy MDA>=M
                                      !! Either M>=N or M<N is permitted.  There
                                      !! is no restriction on the rank of A.  The
                                      !! condition MDA<M is considered an error.
                                      !!
                                      !! The contents of the array A(*,*) will be
                                      !! modified by the subroutine. These contents
                                      !! are not generally required by the user.
   integer,intent(in) :: m
   integer,intent(in) :: n
   real(wp),intent(inout) :: b(mdb,*) !! `(B(MDB,NB) or B(M))`.
                                      !! If NB = 0 the subroutine will perform the
                                      !! orthogonal decomposition but will make no
                                      !! references to the array B(*).  If NB>0
                                      !! the array B(*) must initially contain the M by
                                      !! NB matrix B of the least squares problem AX =
                                      !! B.  If NB>=2 the array B(*) must be doubly
                                      !! subscripted with first dimensioning parameter
                                      !! MDB>=MAX(M,N).  If NB = 1 the array B(*) may
                                      !! be either doubly or singly subscripted.  In
                                      !! the latter case the value of MDB is arbitrary
                                      !! but it should be set to some valid integer
                                      !! value such as MDB = M.
                                      !!
                                      !! The condition of NB>1.AND.MDB< MAX(M,N)
                                      !! is considered an error.
                                      !!
                                      !! On return the array B(*) will contain the N by
                                      !! NB solution matrix X.
   integer,intent(in) :: nb
   real(wp),intent(in) :: tau !! Absolute tolerance parameter provided by user
                              !! for pseudorank determination.
   integer,intent(out) :: krank !! Set by the subroutine to indicate the
                                !! pseudorank of A.
   real(wp),intent(out)  :: rnorm(*) !! `RNORM(NB)`.
                                     !! On return, RNORM(J) will contain the Euclidean
                                     !! norm of the residual vector for the problem
                                     !! defined by the J-th column vector of the array
                                     !! B(*,*) for J = 1,...,NB.
   real(wp) :: h(*) !! `H(N)`. Array of working space used by DHFTI.
                    !! On return, contains
                    !! elements of the pre-multiplying
                    !! Householder transformations used to compute
                    !! the minimum Euclidean length solution.
                    !! not generally required by the user.
   real(wp) :: g(*) !! `G(N)`. Array of working space used by DHFTI.
                    !! On return, contain
                    !! elements of the post-multiplying
                    !! Householder transformations used to compute
                    !! the minimum Euclidean length solution.
                    !! not generally required by the user.
   integer :: ip(*) !! `IP(N)`. Array of working space used by DHFTI.
                    !! Array in which the subroutine records indices
                    !! describing the permutation of column vectors.
                    !! not generally required by the user.

   integer :: i, ii, iopt, ip1, j, jb, jj, k, kp1, l, ldiag, lmax, nerr
   real(wp) :: dzero, factor, hmax, sm, sm1, szero, tmp
   logical :: lmax_found

   szero = 0.0_wp
   dzero = 0.0_wp
   factor = 0.001_wp

   k = 0
   ldiag = min(m,n)
   if (ldiag > 0) then

      if (mda < m) then
         nerr = 1
         iopt = 2
         write(*,*) 'MDA<M, PROBABLE ERROR.'
         return
      end if
      if (nb > 1 .and. max(m,n) > mdb) then
         nerr = 2
         iopt = 2
         write(*,*) 'MDB<MAX(M,N).AND.NB>1. PROBABLE ERROR.'
         return
      end if

      do j = 1, ldiag
         lmax_found = .false.
         if (j /= 1) then
            ! UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
            lmax = j
            do l = j, n
               h(l) = h(l) - a(j-1,l)**2
               if (h(l) > h(lmax)) lmax = l
            end do
            lmax_found = (factor*h(lmax) > hmax*drelpr)
         end if
         if (.not. lmax_found) then
            ! COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
            lmax = j
            do l = j, n
               h(l) = 0.0_wp
               do i = j, m
                  h(l) = h(l) + a(i,l)**2
               end do
               if (h(l) > h(lmax)) lmax = l
            end do
            hmax = h(lmax)
         end if
         ! LMAX HAS BEEN DETERMINED
         ! DO COLUMN INTERCHANGES IF NEEDED.
         ip(j) = lmax
         if (ip(j) /= j) then
            do i = 1, m
               tmp = a(i,j)
               a(i,j) = a(i,lmax)
               a(i,lmax) = tmp
            end do
            h(lmax) = h(j)
         end if
         ! COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
         ! AND B.
         call dh12(1,j,j+1,m,a(1,j),1,h(j),a(1,j+1),1,mda,n-j)
         call dh12(2,j,j+1,m,a(1,j),1,h(j),b,1,mdb,nb)
      end do

      ! DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
      do j = 1, ldiag
         if (abs(a(j,j)) <= tau) then
            k = j - 1
            exit
         else
            if (j==ldiag) k = ldiag
         end if
      end do

         kp1 = k + 1
         ! COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
         if (nb >= 1) then
            do jb = 1, nb
               tmp = szero
               if (m >= kp1) then
                  do i = kp1, m
                     tmp = tmp + b(i,jb)**2
                  end do
               end if
               rnorm(jb) = sqrt(tmp)
            end do
         end if
         ! SPECIAL FOR PSEUDORANK = 0
         if (k > 0) then
            ! IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
            ! DECOMPOSITION OF FIRST K ROWS.

            if (k /= n) then
               do ii = 1, k
                  i = kp1 - ii
                  call dh12(1,i,kp1,n,a(i,1),mda,g(i),a,mda,1,i-1)
               end do
            end if

            if (nb >= 1) then
               do jb = 1, nb
                  ! SOLVE THE K BY K TRIANGULAR SYSTEM.
                  do l = 1, k
                     sm = dzero
                     i = kp1 - l
                     ip1 = i + 1
                     if (k >= ip1) then
                        do j = ip1, k
                           sm = sm + a(i,j)*b(j,jb)
                        end do
                     end if
                     sm1 = sm
                     b(i,jb) = (b(i,jb) - sm1)/a(i,i)
                  end do

                  ! COMPLETE COMPUTATION OF SOLUTION VECTOR.
                  if (k /= n) then
                     do j = kp1, n
                        b(j,jb) = szero
                     end do
                     do i = 1, k
                        call dh12(2,i,kp1,n,a(i,1),mda,g(i),b(1,jb),1,mdb,1)
                     end do
                  end if

                  ! RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
                  ! COLUMN INTERCHANGES.

                  do jj = 1, ldiag
                     j = ldiag + 1 - jj
                     if (ip(j) /= j) then
                        l = ip(j)
                        tmp = b(l,jb)
                        b(l,jb) = b(j,jb)
                        b(j,jb) = tmp
                     end if
                  end do

               end do

            end if

         elseif (nb >= 1) then
           do jb = 1, nb
              do i = 1, n
                 b(i,jb) = szero
             end do
           end do
         end if

   end if

   ! THE SOLUTION VECTORS, X, ARE NOW
   ! IN THE FIRST N ROWS OF THE ARRAY B(,).

   krank = k

   end subroutine dhfti
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determine an N1-vector W, and
!            an N2-vector Z
!  which minimizes the Euclidean length of W
!  subject to G*W+H*Z >= Y.
!  This is the least projected distance problem, LPDP.
!  The matrices G and H are of respective
!  dimensions M by N1 and M by N2.
!
!  Called by subprogram [[DLSI]].
!
!```
!  The matrix
!             (G H Y)
!
!  occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
!
!  The solution (W) is returned in X(*).
!               (Z)
!```
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)
!  * 910408  Updated the AUTHOR section.  (WRB)

   subroutine dlpdp (a, mda, m, n1, n2, prgopt, x, wnorm, mode, ws, is)

   integer,intent(in) :: mda
   integer :: m
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   real(wp) :: a(mda,*) !! `A(MDA,N+1)`, where `N=N1+N2`.
   real(wp) :: prgopt(*)
   real(wp) :: x(*) !! `X(N)`, where `N=N1+N2`.
   real(wp) :: wnorm
   integer,intent(out) :: mode !! The value of MODE indicates the status of
                               !! the computation after returning to the user.
                               !!
                               !!  * `MODE=1` The solution was successfully obtained.
                               !!  * `MODE=2` The inequalities are inconsistent.
   real(wp) :: ws(*) !! `WS((M+2)*(N+7))`, where `N=N1+N2`. This is a slight overestimate for WS(*).
   integer :: is(*) !! `IS(M+N+1)`, where `N=N1+N2`.

   integer :: i, iw, ix, j, l, modew, n, np1
   real(wp) :: rnorm, sc, ynorm

   real(wp),parameter :: zero = 0.0_wp
   real(wp),parameter :: one = 1.0_wp
   real(wp),parameter :: fac = 0.1_wp

   n = n1 + n2
   mode = 1
   if (m <= 0) then
      if (n > 0) then
         x(1) = zero
         call dcopy(n,x,0,x,1)
      end if
      wnorm = zero
      return
   end if

      np1 = n + 1

      ! SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
      do i = 1, m
         sc = dnrm2(n,a(i,1),mda)
         if (sc /= zero) then
            sc = one/sc
            call dscal(np1,sc,a(i,1),mda)
         end if
      end do

      ! SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
      ynorm = dnrm2(m,a(1,np1),1)
      if (ynorm /= zero) then
         sc = one/ynorm
         call dscal(m,sc,a(1,np1),1)
      end if

      ! SCALE COLS OF MATRIX H.
      j = n1 + 1
      do
         if (j > n) exit
         sc = dnrm2(m,a(1,j),1)
         if (sc /= zero) sc = one/sc
         call dscal(m,sc,a(1,j),1)
         x(j) = sc
         j = j + 1
      end do

      if (n1 > 0) then

         ! COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
         iw = 0
         do i = 1, m

            ! MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
            call dcopy(n2,a(i,n1+1),mda,ws(iw+1),1)
            iw = iw + n2

            ! MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
            call dcopy(n1,a(i,1),mda,ws(iw+1),1)
            iw = iw + n1

            ! MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
            ws(iw+1) = a(i,np1)
            iw = iw + 1
         end do
         ws(iw+1) = zero
         call dcopy(n,ws(iw+1),0,ws(iw+1),1)
         iw = iw + n
         ws(iw+1) = one
         iw = iw + 1

         ! SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U>=0.  THE
         ! MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
         ! F = TRANSPOSE OF (0,...,0,1).
         ix = iw + 1
         iw = iw + m

         ! DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
         ! DWNNLS( ).
         is(1) = 0
         is(2) = 0
         call dwnnls(ws,np1,n2,np1-n2,m,0,prgopt,ws(ix),rnorm, &
                     modew,is,ws(iw+1))

         ! COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
         sc = one - ddot(m,a(1,np1),1,ws(ix),1)
         if (one + fac*abs(sc) == one .or. rnorm <= zero) then
            mode = 2
            return
         end if
         sc = one/sc
         do j = 1, n1
            x(j) = sc*ddot(m,a(1,j),1,ws(ix),1)
         end do
         ! COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
         ! VECTOR.
         do i = 1, m
            a(i,np1) = a(i,np1) - ddot(n1,a(i,1),mda,x,1)
         end do

      end if

      if (n2 > 0) then

         ! COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
         iw = 0
         do i = 1, m
            call dcopy(n2,a(i,n1+1),mda,ws(iw+1),1)
            iw = iw + n2
            ws(iw+1) = a(i,np1)
            iw = iw + 1
         end do
         ws(iw+1) = zero
         call dcopy(n2,ws(iw+1),0,ws(iw+1),1)
         iw = iw + n2
         ws(iw+1) = one
         iw = iw + 1
         ix = iw + 1
         iw = iw + m

         ! SOLVE RV=S SUBJECT TO V>=0.  THE MATRIX R =(TRANSPOSE
         ! OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
         ! OF (0,...,0,1)).
         !
         ! DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
         ! DWNNLS( ).
         is(1) = 0
         is(2) = 0
         call dwnnls(ws,n2+1,0,n2+1,m,0,prgopt,ws(ix),rnorm,modew, &
                     is,ws(iw+1))

         ! COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
         sc = one - ddot(m,a(1,np1),1,ws(ix),1)
         if (one + fac*abs(sc) == one .or. rnorm <= zero) then
            mode = 2
            return
         end if
         sc = one/sc
         do j = 1, n2
            l = n1 + j
            x(l) = sc*ddot(m,a(1,l),1,ws(ix),1)*x(l)
         end do
      end if

      ! ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
      call dscal(n,ynorm,x,1)
      wnorm = dnrm2(n1,x,1)

   end subroutine dlpdp
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subprogram solves a linearly constrained least squares
!  problem with both equality and inequality constraints, and, if the
!  user requests, obtains a covariance matrix of the solution
!  parameters.
!
!  Suppose there are given matrices E, A and G of respective
!  dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
!  respective lengths ME, MA and MG.  This subroutine solves the
!  linearly constrained least squares problem
!
!  * `EX = F, (E ME by N)` (equations to be exactly satisfied)
!  * `AX = B, (A MA by N)` (equations to be approximately satisfied, least squares sense)
!  * `GX >= H,(G MG by N)` (inequality constraints)
!
!  The inequalities GX >= H mean that every component of the
!  product GX must be >= the corresponding component of H.
!
!  In case the equality constraints cannot be satisfied, a
!  generalized inverse solution residual vector length is obtained
!  for F-EX.  This is the minimal length possible for F-EX.
!
!  Any values ME >= 0, MA >= 0, or MG >= 0 are permitted.  The
!  rank of the matrix E is estimated during the computation.  We call
!  this value KRANKE.  It is an output parameter in IP(1) defined
!  below.  Using a generalized inverse solution of EX=F, a reduced
!  least squares problem with inequality constraints is obtained.
!  The tolerances used in these tests for determining the rank
!  of E and the rank of the reduced least squares problem are
!  given in Sandia Tech. Rept. SAND-78-1290.  They can be
!  modified by the user if new values are provided in
!  the option list of the array PRGOPT(*).
!
!  The user must dimension all arrays appearing in the call list..
!  W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
!  where K=MAX(MA+MG,N).  This allows for a solution of a range of
!  problems in the given working space.  The dimension of WS(*)
!  given is a necessary overestimate.  Once a particular problem
!  has been run, the output parameter IP(3) gives the actual
!  dimension required for that problem.
!
!  The parameters for [[DLSEI]] are
!
!```
!  Input.. All TYPE REAL variables are DOUBLE PRECISION
!
!  W(*,*),MDW,   The array W(*,*) is doubly subscripted with
!  ME,MA,MG,N    first dimensioning parameter equal to MDW.
!                For this discussion let us call M = ME+MA+MG.  Then
!                MDW must satisfy MDW >= M.  The condition
!                MDW < M is an error.
!
!                The array W(*,*) contains the matrices and vectors
!
!                               (E  F)
!                               (A  B)
!                               (G  H)
!
!                in rows and columns 1,...,M and 1,...,N+1
!                respectively.
!
!                The integers ME, MA, and MG are the
!                respective matrix row dimensions
!                of E, A and G.  Each matrix has N columns.
!
!  PRGOPT(*)    This real-valued array is the option vector.
!               If the user is satisfied with the nominal
!               subprogram features set
!
!               PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!               Otherwise PRGOPT(*) is a linked list consisting of
!               groups of data of the following form
!
!               LINK
!               KEY
!               DATA SET
!
!               The parameters LINK and KEY are each one word.
!               The DATA SET can be comprised of several words.
!               The number of items depends on the value of KEY.
!               The value of LINK points to the first
!               entry of the next group of data within
!               PRGOPT(*).  The exception is when there are
!               no more options to change.  In that
!               case, LINK=1 and the values KEY and DATA SET
!               are not referenced.  The general layout of
!               PRGOPT(*) is as follows.
!
!            ...PRGOPT(1) = LINK1 (link to first entry of next group)
!            .  PRGOPT(2) = KEY1 (key to the option change)
!            .  PRGOPT(3) = data value (data value for this change)
!            .       .
!            .       .
!            .       .
!            ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
!            .                       next group)
!            .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
!            .  PRGOPT(LINK1+2) = data value
!            ...     .
!            .       .
!            .       .
!            ...PRGOPT(LINK) = 1 (no more options to change)
!
!               Values of LINK that are nonpositive are errors.
!               A value of LINK > NLINK=100000 is also an error.
!               This helps prevent using invalid but positive
!               values of LINK that will probably extend
!               beyond the program limits of PRGOPT(*).
!               Unrecognized values of KEY are ignored.  The
!               order of the options is arbitrary and any number
!               of options can be changed with the following
!               restriction.  To prevent cycling in the
!               processing of the option array, a count of the
!               number of options changed is maintained.
!               Whenever this count exceeds NOPT=1000, an error
!               message is printed and the subprogram returns.
!
!               Options..
!
!               KEY=1
!                      Compute in W(*,*) the N by N
!               covariance matrix of the solution variables
!               as an output parameter.  Nominally the
!               covariance matrix will not be computed.
!               (This requires no user input.)
!               The data set for this option is a single value.
!               It must be nonzero when the covariance matrix
!               is desired.  If it is zero, the covariance
!               matrix is not computed.  When the covariance matrix
!               is computed, the first dimensioning parameter
!               of the array W(*,*) must satisfy MDW >= MAX(M,N).
!
!               KEY=10
!                      Suppress scaling of the inverse of the
!               normal matrix by the scale factor RNORM**2/
!               MAX(1, no. of degrees of freedom).  This option
!               only applies when the option for computing the
!               covariance matrix (KEY=1) is used.  With KEY=1 and
!               KEY=10 used as options the unscaled inverse of the
!               normal matrix is returned in W(*,*).
!               The data set for this option is a single value.
!               When it is nonzero no scaling is done.  When it is
!               zero scaling is done.  The nominal case is to do
!               scaling so if option (KEY=1) is used alone, the
!               matrix will be scaled on output.
!
!               KEY=2
!                      Scale the nonzero columns of the
!                      entire data matrix.
!               (E)
!               (A)
!               (G)
!
!               to have length one.  The data set for this
!               option is a single value.  It must be
!               nonzero if unit length column scaling
!               is desired.
!
!               KEY=3
!                      Scale columns of the entire data matrix
!               (E)
!               (A)
!               (G)
!
!               with a user-provided diagonal matrix.
!               The data set for this option consists
!               of the N diagonal scaling factors, one for
!               each matrix column.
!
!               KEY=4
!                      Change the rank determination tolerance for
!               the equality constraint equations from
!               the nominal value of SQRT(DRELPR).  This quantity can
!               be no smaller than DRELPR, the arithmetic-
!               storage precision.  The quantity DRELPR is the
!               largest positive number such that T=1.+DRELPR
!               satisfies T == 1.  The quantity used
!               here is internally restricted to be at
!               least DRELPR.  The data set for this option
!               is the new tolerance.
!
!               KEY=5
!                      Change the rank determination tolerance for
!               the reduced least squares equations from
!               the nominal value of SQRT(DRELPR).  This quantity can
!               be no smaller than DRELPR, the arithmetic-
!               storage precision.  The quantity used
!               here is internally restricted to be at
!               least DRELPR.  The data set for this option
!               is the new tolerance.
!
!               For example, suppose we want to change
!               the tolerance for the reduced least squares
!               problem, compute the covariance matrix of
!               the solution parameters, and provide
!               column scaling for the data matrix.  For
!               these options the dimension of PRGOPT(*)
!               must be at least N+9.  The Fortran statements
!               defining these options would be as follows:
!
!               PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
!               PRGOPT(2)=1 (covariance matrix key)
!               PRGOPT(3)=1 (covariance matrix wanted)
!
!               PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
!               PRGOPT(5)=5 (least squares equas.  tolerance key)
!               PRGOPT(6)=... (new value of the tolerance)
!
!               PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
!               PRGOPT(8)=3 (user-provided column scaling key)
!
!               CALL DCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
!                 scaling factors from the user array D(*)
!                 to PRGOPT(9)-PRGOPT(N+8))
!
!               PRGOPT(N+9)=1 (no more options to change)
!
!               The contents of PRGOPT(*) are not modified
!               by the subprogram.
!               The options for WNNLS( ) can also be included
!               in this array.  The values of KEY recognized
!               by WNNLS( ) are 6, 7 and 8.  Their functions
!               are documented in the usage instructions for
!               subroutine WNNLS( ).  Normally these options
!               do not need to be modified when using [[DLSEI]].
!
!  IP(1),       The amounts of working storage actually
!  IP(2)        allocated for the working arrays WS(*) and
!               IP(*), respectively.  These quantities are
!               compared with the actual amounts of storage
!               needed by [[DLSEI]].  Insufficient storage
!               allocated for either WS(*) or IP(*) is an
!               error.  This feature was included in [[DLSEI]]
!               because miscalculating the storage formulas
!               for WS(*) and IP(*) might very well lead to
!               subtle and hard-to-find execution errors.
!
!               The length of WS(*) must be at least
!
!               LW = 2*(ME+N)+K+(MG+2)*(N+7)
!
!               where K = max(MA+MG,N)
!               This test will not be made if IP(1)<=0.
!
!               The length of IP(*) must be at least
!
!               LIP = MG+2*N+2
!               This test will not be made if IP(2)<=0.
!
!  Output.. All TYPE REAL variables are DOUBLE PRECISION
!
!  X(*),RNORME,  The array X(*) contains the solution parameters
!  RNORML        if the integer output flag MODE = 0 or 1.
!                The definition of MODE is given directly below.
!                When MODE = 0 or 1, RNORME and RNORML
!                respectively contain the residual vector
!                Euclidean lengths of F - EX and B - AX.  When
!                MODE=1 the equality constraint equations EX=F
!                are contradictory, so RNORME /= 0.  The residual
!                vector F-EX has minimal Euclidean length.  For
!                MODE >= 2, none of these parameters is defined.
!
!  MODE          Integer flag that indicates the subprogram
!                status after completion.  If MODE >= 2, no
!                solution has been computed.
!
!                MODE =
!
!                0  Both equality and inequality constraints
!                   are compatible and have been satisfied.
!
!                1  Equality constraints are contradictory.
!                   A generalized inverse solution of EX=F was used
!                   to minimize the residual vector length F-EX.
!                   In this sense, the solution is still meaningful.
!
!                2  Inequality constraints are contradictory.
!
!                3  Both equality and inequality constraints
!                   are contradictory.
!
!                The following interpretation of
!                MODE=1,2 or 3 must be made.  The
!                sets consisting of all solutions
!                of the equality constraints EX=F
!                and all vectors satisfying GX >= H
!                have no points in common.  (In
!                particular this does not say that
!                each individual set has no points
!                at all, although this could be the
!                case.)
!
!                4  Usage error occurred.  The value
!                   of MDW is < ME+MA+MG, MDW is
!                   < N and a covariance matrix is
!                   requested, or the option vector
!                   PRGOPT(*) is not properly defined,
!                   or the lengths of the working arrays
!                   WS(*) and IP(*), when specified in
!                   IP(1) and IP(2) respectively, are not
!                   long enough.
!
!  W(*,*)        The array W(*,*) contains the N by N symmetric
!                covariance matrix of the solution parameters,
!                provided this was requested on input with
!                the option vector PRGOPT(*) and the output
!                flag is returned with MODE = 0 or 1.
!
!  IP(*)         The integer working array has three entries
!                that provide rank and working array length
!                information after completion.
!
!                   IP(1) = rank of equality constraint
!                           matrix.  Define this quantity
!                           as KRANKE.
!
!                   IP(2) = rank of reduced least squares
!                           problem.
!
!                   IP(3) = the amount of storage in the
!                           working array WS(*) that was
!                           actually used by the subprogram.
!                           The formula given above for the length
!                           of WS(*) is a necessary overestimate.
!                           If exactly the same problem matrices
!                           are used in subsequent executions,
!                           the declared dimension of WS(*) can
!                           be reduced to this output value.
!  User Designated
!  Working Arrays..
!
!  WS(*),IP(*)              These are respectively type real
!                           and type integer working arrays.
!                           Their required minimal lengths are
!                           given above.
!```
!
!### References
!  * K. H. Haskell and R. J. Hanson, An algorithm for
!    linear least squares problems with equality and
!    nonnegativity constraints, Report SAND77-0552, Sandia
!    Laboratories, June 1978.
!  * K. H. Haskell and R. J. Hanson, Selected algorithms for
!    the linearly constrained least squares problem - a
!    users guide, Report SAND78-1290, Sandia Laboratories,
!    August 1979.
!  * K. H. Haskell and R. J. Hanson, An algorithm for
!    linear least squares problems with equality and
!    nonnegativity constraints, Mathematical Programming
!    21 (1981), pp. 98-118.
!  * R. J. Hanson and K. H. Haskell, Two algorithms for the
!    linearly constrained least squares problem, ACM
!    Transactions on Mathematical Software, September 1982.
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and extensively revised (WRB & RWC)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!  * 900604  DP version created from SP version.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dlsei (w, mdw, me, ma, mg, n, prgopt, x, rnorme, &
                     rnorml, mode, ws, ip)

   integer,intent(in) :: mdw
   real(wp) :: w(mdw,*)
   integer :: me
   integer :: ma
   integer :: mg
   integer :: n
   real(wp) :: prgopt(*)
   real(wp) :: x(*)
   real(wp) :: rnorme
   real(wp) :: rnorml
   integer :: mode
   real(wp) :: ws(*)
   integer :: ip(3)

   real(wp) :: enorm, fnorm, gam, rb, rn, rnmax, size, &
               sn, snmax, t, tau, uj, up, vj, xnorm, xnrme
   integer :: i, imax, j, jp1, k, key, kranke, last, lchk, link, m, &
              mapke1, mdeqc, mend, mep1, n1, n2, next, nlink, nopt, np1, &
              ntimes
   logical :: cov, done
   character(len=8) :: xern1, xern2, xern3, xern4

   ! Set the nominal tolerance used in the code for the equality
   ! constraint equations.
   tau = sqrt(drelpr)

   ! Check that enough storage was allocated in WS(*) and IP(*).
   mode = 4
   if (min(n,me,ma,mg) < 0) then
      write (xern1, '(I8)') n
      write (xern2, '(I8)') me
      write (xern3, '(I8)') ma
      write (xern4, '(I8)') mg
      write(*,*) 'ALL OF THE VARIABLES N, ME,' // &
                 ' MA, MG MUST BE >= 0. ENTERED ROUTINE WITH: ' // &
                    'N = ' // trim(adjustl(xern1)) // &
                 ', ME = ' // trim(adjustl(xern2)) // &
                 ', MA = ' // trim(adjustl(xern3)) // &
                 ', MG = ' // trim(adjustl(xern4))
      return
   endif

   if (ip(1)>0) then
      lchk = 2*(me+n) + max(ma+mg,n) + (mg+2)*(n+7)
      if (ip(1)<lchk) then
         write (xern1, '(I8)') lchk
         write(*,*) 'INSUFFICIENT STORAGE ' // &
                    'ALLOCATED FOR WS(*), NEED LW = ' // xern1
         return
      endif
   endif

   if (ip(2)>0) then
      lchk = mg + 2*n + 2
      if (ip(2)<lchk) then
         write (xern1, '(I8)') lchk
         write(*,*) 'INSUFFICIENT STORAGE ' // &
                    'ALLOCATED FOR IP(*), NEED LIP = ' // xern1
         return
      endif
   endif

   ! Compute number of possible right multiplying Householder
   ! transformations.

   m = me + ma + mg
   if (n<=0 .or. m<=0) then
      mode = 0
      rnorme = 0
      rnorml = 0
      return
   endif

   if (mdw<m) then
      write(*,*) 'MDW < ME+MA+MG IS AN ERROR'
      return
   endif

   np1 = n + 1
   kranke = min(me,n)
   n1 = 2*kranke + 1
   n2 = n1 + n

   ! Set nominal values.
   !
   ! The nominal column scaling used in the code is
   ! the identity scaling.

   call dcopy (n, [1.0_wp], 0, ws(n1), 1)

   ! No covariance matrix is nominally computed.

   cov = .false.

   ! Process option vector.
   ! Define bound for number of options to change.

   nopt = 1000
   ntimes = 0

   ! Define bound for positive values of LINK.

   nlink = 100000
   last = 1
   link = prgopt(1)
   if (link==0 .or. link>nlink) then
      write(*,*) 'THE OPTION VECTOR IS UNDEFINED'
      return
   endif

   do
      if (link<=1) exit
      ntimes = ntimes + 1
      if (ntimes>nopt) then
         write(*,*) 'THE LINKS IN THE OPTION VECTOR ARE CYCLING.'
         return
      endif

      key = prgopt(last+1)
      if (key==1) then
         cov = prgopt(last+2) /= 0.0_wp
      elseif (key==2 .and. prgopt(last+2)/=0.0_wp) then
         do j = 1,n
            t = dnrm2(m,w(1,j),1)
            if (t/=0.0_wp) t = 1.0_wp/t
            ws(j+n1-1) = t
         end do
      elseif (key==3) then
         call dcopy (n, prgopt(last+2), 1, ws(n1), 1)
      elseif (key==4) then
         tau = max(drelpr,prgopt(last+2))
      endif

      next = prgopt(link)
      if (next<=0 .or. next>nlink) then
         write(*,*) 'THE OPTION VECTOR IS UNDEFINED'
         return
      endif

      last = link
      link = next
   end do

   do j = 1,n
      call dscal (m, ws(n1+j-1), w(1,j), 1)
   end do

   if (cov .and. mdw<n) then
      write(*,*) 'MDW < N WHEN COV MATRIX NEEDED, IS AN ERROR'
      return
   endif

   ! Problem definition and option vector OK.

   mode = 0

   ! Compute norm of equality constraint matrix and right side.

   enorm = 0.0_wp
   do j = 1,n
      enorm = max(enorm,dasum(me,w(1,j),1))
   end do

   fnorm = dasum(me,w(1,np1),1)
   snmax = 0.0_wp
   rnmax = 0.0_wp
   do i = 1,kranke

      ! Compute maximum ratio of vector lengths. Partition is at
      ! column I.

      do k = i,me
         sn = ddot(n-i+1,w(k,i),mdw,w(k,i),mdw)
         rn = ddot(i-1,w(k,1),mdw,w(k,1),mdw)
         if (rn==0.0_wp .and. sn>snmax) then
            snmax = sn
            imax = k
         elseif (k==i .or. sn*rnmax>rn*snmax) then
            snmax = sn
            rnmax = rn
            imax = k
         endif
      end do

      ! Interchange rows if necessary.

      if (i/=imax) call dswap (np1, w(i,1), mdw, w(imax,1), mdw)
      if (snmax>rnmax*tau**2) then
         ! Eliminate elements I+1,...,N in row I.
         call dh12 (1, i, i+1, n, w(i,1), mdw, ws(i), w(i+1,1), mdw, 1, m-i)
      else
         kranke = i - 1
         exit
      endif
   end do

   ! Save diagonal terms of lower trapezoidal matrix.

   call dcopy (kranke, w, mdw+1, ws(kranke+1), 1)

   ! Use Householder transformation from left to achieve
   ! KRANKE by KRANKE upper triangular form.

   if (kranke<me) then
      do k = kranke,1,-1
         ! Apply transformation to matrix cols. 1,...,K-1.
         call dh12 (1, k, kranke+1, me, w(1,k), 1, up, w, 1, mdw, k-1)
         ! Apply to rt side vector.
         call dh12 (2, k, kranke+1, me, w(1,k), 1, up, w(1,np1), 1, 1, 1)
      end do
   endif

   ! Solve for variables 1,...,KRANKE in new coordinates.

   call dcopy (kranke, w(1, np1), 1, x, 1)
   do i = 1,kranke
      x(i) = (x(i)-ddot(i-1,w(i,1),mdw,x,1))/w(i,i)
   end do

   ! Compute residuals for reduced problem.

   mep1 = me + 1
   rnorml = 0.0_wp
   do i = mep1,m
      w(i,np1) = w(i,np1) - ddot(kranke,w(i,1),mdw,x,1)
      sn = ddot(kranke,w(i,1),mdw,w(i,1),mdw)
      rn = ddot(n-kranke,w(i,kranke+1),mdw,w(i,kranke+1),mdw)
      if (rn<=sn*tau**2 .and. kranke<n) &
         call dcopy (n-kranke, [0.0_wp], 0, w(i,kranke+1), mdw)
   end do

   ! Compute equality constraint equations residual length.

   rnorme = dnrm2(me-kranke,w(kranke+1,np1),1)

   ! Move reduced problem data upward if KRANKE<ME.

   if (kranke<me) then
      do j = 1,np1
         call dcopy (m-me, w(me+1,j), 1, w(kranke+1,j), 1)
      end do
   endif

   ! Compute solution of reduced problem.

   call dlsi(w(kranke+1, kranke+1), mdw, ma, mg, n-kranke, prgopt, &
             x(kranke+1), rnorml, mode, ws(n2), ip(2))

   ! Test for consistency of equality constraints.
   done = .false.

   if (me>0) then
      mdeqc = 0
      xnrme = dasum(kranke,w(1,np1),1)
      if (rnorme>tau*(enorm*xnrme+fnorm)) mdeqc = 1
      mode = mode + mdeqc

      ! Check if solution to equality constraints satisfies inequality
      ! constraints when there are no degrees of freedom left.

      if (kranke==n .and. mg>0) then
         xnorm = dasum(n,x,1)
         mapke1 = ma + kranke + 1
         mend = ma + kranke + mg
         do i = mapke1,mend
            size = dasum(n,w(i,1),mdw)*xnorm + abs(w(i,np1))
            if (w(i,np1)>tau*size) then
               mode = mode + 2
               done = .true.
               exit
            endif
         end do
      endif
   endif

   if (.not. done) then
      ! Replace diagonal terms of lower trapezoidal matrix.

      if (kranke>0) then
         call dcopy (kranke, ws(kranke+1), 1, w, mdw+1)
         ! Reapply transformation to put solution in original coordinates.
         do i = kranke,1,-1
            call dh12 (2, i, i+1, n, w(i,1), mdw, ws(i), x, 1, 1, 1)
         end do

         ! Compute covariance matrix of equality constrained problem.

         if (cov) then
            do j = min(kranke,n-1),1,-1
               rb = ws(j)*w(j,j)
               if (rb/=0.0_wp) rb = 1.0_wp/rb
               jp1 = j + 1
               do i = jp1,n
                  w(i,j) = rb*ddot(n-j,w(i,jp1),mdw,w(j,jp1),mdw)
               end do
               gam = 0.5_wp*rb*ddot(n-j,w(jp1,j),1,w(j,jp1),mdw)
               call daxpy (n-j, gam, w(j,jp1), mdw, w(jp1,j), 1)
               do i = jp1,n
                  do k = i,n
                     w(i,k) = w(i,k) + w(j,i)*w(k,j) + w(i,j)*w(j,k)
                     w(k,i) = w(i,k)
                  end do
               end do
               uj = ws(j)
               vj = gam*uj
               w(j,j) = uj*vj + uj*vj
               do i = jp1,n
                  w(j,i) = uj*w(i,j) + vj*w(j,i)
               end do
               call dcopy (n-j, w(j, jp1), mdw, w(jp1,j), 1)
            end do
         endif
      endif

      ! Apply the scaling to the covariance matrix.

      if (cov) then
         do i = 1,n
            call dscal (n, ws(i+n1-1), w(i,1), mdw)
            call dscal (n, ws(i+n1-1), w(1,i), 1)
         end do
      endif

   end if

   ! Rescale solution vector.

   if (mode<=1) then
      do j = 1,n
         x(j) = x(j)*ws(n1+j-1)
      end do
   endif

   ip(1) = kranke
   ip(3) = ip(3) + 2*kranke + n

   end subroutine dlsei
!*****************************************************************************************

!*****************************************************************************************
!>
!  This is a companion subprogram to [[DLSEI]].  The documentation for
!  [[DLSEI]] has complete usage instructions.
!
!  Solve:
!```
!    AX = B,  A  MA by N  (least squares equations)
!```
!
!  subject to:
!```
!    GX>=H, G  MG by N  (inequality constraints)
!```
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and extensively revised (WRB & RWC)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)
!  * 900604  DP version created from SP version.  (RWC)
!  * 920422  Changed CALL to DHFTI to include variable MA.  (WRB)

   subroutine dlsi (w, mdw, ma, mg, n, prgopt, x, rnorm, mode, ws, ip)

   integer,intent(in) :: mdw !! contain (resp) var. dimension of `W(*,*)`, and matrix dimensions.
   integer,intent(in) :: ma  !! contain (resp) var. dimension of `W(*,*)`, and matrix dimensions.
   integer,intent(in) :: mg  !! contain (resp) var. dimension of `W(*,*)`, and matrix dimensions.
   integer,intent(in) :: n   !! contain (resp) var. dimension of `W(*,*)`, and matrix dimensions.
   real(wp) :: w(mdw,*) !! `W(*,*)` contains:
                        !!
                        !!```
                        !! (A B)
                        !! (G H)
                        !!```
                        !!
                        !! in rows `1,...,MA+MG`,
                        !! cols `1,...,N+1`.
   real(wp),intent(in) :: prgopt(*) !! Program option vector.
   real(wp),intent(out) :: x(*) !! Solution vector(unless MODE=2)
   real(wp),intent(out) :: rnorm !! length of AX-B.
   integer,intent(out) :: mode !! * `=0`   Inequality constraints are compatible.
                               !! * `=2`   Inequality constraints contradictory.
   real(wp) :: ws(*) !! Working storage of dimension `K+N+(MG+2)*(N+7)`,
                     !! where `K=MAX(MA+MG,N)`.
   integer :: ip(*) !! `IP(MG+2*N+1)` Integer working storage

   real(wp) :: anorm, fac, gam, rb, tau, tol, xnorm
   integer :: i, j, k, key, krank, krm1, krp1, l, last, link, m, map1, &
              mdlpdp, minman, n1, n2, n3, next, np1
   logical :: cov, sclcov
   real(wp) :: rnorm_(1) !! JW added for call to [[dhfti]]

   ! Set the nominal tolerance used in the code.
   tol = sqrt(drelpr)

   mode = 0
   rnorm = 0.0_wp
   m = ma + mg
   np1 = n + 1
   krank = 0

   main : block

      if (n<=0 .or. m<=0) exit main

      ! To process option vector.

      cov = .false.
      sclcov = .true.
      last = 1
      link = prgopt(1)

      do
         if (link<=1) exit
         key = prgopt(last+1)
         if (key==1) cov = prgopt(last+2) /= 0.0_wp
         if (key==10) sclcov = prgopt(last+2) == 0.0_wp
         if (key==5) tol = max(drelpr,prgopt(last+2))
         next = prgopt(link)
         last = link
         link = next
      end do

      ! Compute matrix norm of least squares equations.

      anorm = 0.0_wp
      do j = 1,n
         anorm = max(anorm,dasum(ma,w(1,j),1))
      end do

      ! Set tolerance for DHFTI( ) rank test.

      tau = tol*anorm

      ! Compute Householder orthogonal decomposition of matrix.

      call dcopy (n, [0.0_wp], 0, ws, 1)
      call dcopy (ma, w(1, np1), 1, ws, 1)
      k = max(m,n)
      minman = min(ma,n)
      n1 = k + 1
      n2 = n1 + n
      rnorm_(1) = rnorm ! JW
      call dhfti (w, mdw, ma, n, ws, ma, 1, tau, krank, rnorm_, ws(n2), &
                  ws(n1), ip)
      rnorm = rnorm_(1) ! JW
      fac = 1.0_wp
      gam = ma - krank
      if (krank<ma .and. sclcov) fac = rnorm**2/gam

      ! Reduce to DLPDP and solve.

      map1 = ma + 1

      ! Compute inequality rt-hand side for DLPDP.

      if (ma<m) then
         if (minman>0) then
            do i = map1,m
               w(i,np1) = w(i,np1) - ddot(n,w(i,1),mdw,ws,1)
            end do

            ! Apply permutations to col. of inequality constraint matrix.

            do i = 1,minman
               call dswap (mg, w(map1,i), 1, w(map1,ip(i)), 1)
            end do

            ! Apply Householder transformations to constraint matrix.

            if (krank>0 .and. krank<n) then
               do i = krank,1,-1
                  call dh12 (2, i, krank+1, n, w(i,1), mdw, ws(n1+i-1), &
                           w(map1,1), mdw, 1, mg)
               end do
            endif

            ! Compute permuted inequality constraint matrix times r-inv.

            do i = map1,m
               do j = 1,krank
                  w(i,j) = (w(i,j)-ddot(j-1,w(1,j),1,w(i,1),mdw))/w(j,j)
               end do
         end do
         endif

         ! Solve the reduced problem with DLPDP algorithm,
         ! the least projected distance problem.

         call dlpdp(w(map1,1), mdw, mg, krank, n-krank, prgopt, x, &
                  xnorm, mdlpdp, ws(n2), ip(n+1))

         ! Compute solution in original coordinates.

         if (mdlpdp==1) then
            do i = krank,1,-1
               x(i) = (x(i)-ddot(krank-i,w(i,i+1),mdw,x(i+1),1))/w(i,i)
            end do

            ! Apply Householder transformation to solution vector.

            if (krank<n) then
               do i = 1,krank
                  call dh12 (2, i, krank+1, n, w(i,1), mdw, ws(n1+i-1), &
                           x, 1, 1, 1)
               end do
            endif

            ! Repermute variables to their input order.

            if (minman>0) then
               do i = minman,1,-1
                  call dswap (1, x(i), 1, x(ip(i)), 1)
               end do

               ! Variables are now in original coordinates.
               ! Add solution of unconstrained problem.

               do i = 1,n
                  x(i) = x(i) + ws(i)
               end do

               ! Compute the residual vector norm.

               rnorm = sqrt(rnorm**2+xnorm**2)
            endif
         else
            mode = 2
         endif
      else
         call dcopy (n, ws, 1, x, 1)
      endif

      ! Compute covariance matrix based on the orthogonal decomposition
      ! from DHFTI( ).

      if (.not.cov .or. krank<=0) exit main
      krm1 = krank - 1
      krp1 = krank + 1

      ! Copy diagonal terms to working array.

      call dcopy (krank, w, mdw+1, ws(n2), 1)

      ! Reciprocate diagonal terms.

      do j = 1,krank
         w(j,j) = 1.0_wp/w(j,j)
      end do

      ! Invert the upper triangular QR factor on itself.

      if (krank>1) then
         do i = 1,krm1
            do j = i+1,krank
               w(i,j) = -ddot(j-i,w(i,i),mdw,w(i,j),1)*w(j,j)
            end do
         end do
      endif

      ! Compute the inverted factor times its transpose.

      do i = 1,krank
         do j = i,krank
            w(i,j) = ddot(krank+1-j,w(i,j),mdw,w(j,j),mdw)
         end do
      end do

      ! Zero out lower trapezoidal part.
      ! Copy upper triangular to lower triangular part.

      if (krank<n) then
         do j = 1,krank
            call dcopy (j, w(1,j), 1, w(j,1), mdw)
         end do

         do i = krp1,n
            call dcopy (i, [0.0_wp], 0, w(i,1), mdw)
         end do

         ! Apply right side transformations to lower triangle.

         n3 = n2 + krp1
         do i = 1,krank
            l = n1 + i
            k = n2 + i
            rb = ws(l-1)*ws(k-1)

            ! If RB>=0.0_wp, transformation can be regarded as zero.

            if (rb<0.0_wp) then
               rb = 1.0_wp/rb

               ! Store unscaled rank one Householder update in work array.

               call dcopy (n, [0.0_wp], 0, ws(n3), 1)
               l = n1 + i
               k = n3 + i
               ws(k-1) = ws(l-1)

               do j = krp1,n
                  ws(n3+j-1) = w(i,j)
               end do

               do j = 1,n
                  ws(j) = rb*(ddot(j-i,w(j,i),mdw,ws(n3+i-1),1)+ &
                        ddot(n-j+1,w(j,j),1,ws(n3+j-1),1))
               end do

               l = n3 + i
               gam = 0.5_wp*rb*ddot(n-i+1,ws(l-1),1,ws(i),1)
               call daxpy (n-i+1, gam, ws(l-1), 1, ws(i), 1)
               do j = i,n
                  do l = 1,i-1
                     w(j,l) = w(j,l) + ws(n3+j-1)*ws(l)
                  end do
                  do l = i,j
                     w(j,l) = w(j,l) + ws(j)*ws(n3+l-1)+ws(l)*ws(n3+j-1)
                  end do
               end do
            endif
         end do

         ! Copy lower triangle to upper triangle to symmetrize the
         ! covariance matrix.

         do i = 1,n
            call dcopy (i, w(i,1), mdw, w(1,i), 1)
         end do
      endif

      ! Repermute rows and columns.

      do i = minman,1,-1
         k = ip(i)
         if (i/=k) then
            call dswap (1, w(i,i), 1, w(k,k), 1)
            call dswap (i-1, w(1,i), 1, w(1,k), 1)
            call dswap (k-i-1, w(i,i+1), mdw, w(i+1,k), 1)
            call dswap (n-k, w(i, k+1), mdw, w(k, k+1), mdw)
         endif
      end do

      ! Put in normalized residual sum of squares scale factor
      ! and symmetrize the resulting covariance matrix.

      do j = 1,n
         call dscal (j, fac, w(1,j), 1)
         call dcopy (j, w(1,j), 1, w(j,1), mdw)
      end do

   end block main

   ip(1) = krank
   ip(2) = n + max(m,n) + (mg+2)*(n+7)

   end subroutine dlsi
!*****************************************************************************************

!*****************************************************************************************
!>
!  This is a companion subprogram to [[DWNNLS]].
!  The documentation for [[DWNNLS]] has complete usage instructions.
!
!  Note: The `M` by `(N+1)` matrix `W( , )` contains the rt. hand side
!        `B` as the `(N+1)`st col.
!
!  Triangularize `L1` by `L1` subsystem, where `L1=MIN(M,L)`, with
!  col interchanges.
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and revised.  (WRB & RWC)
!  * 890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900328  Added TYPE section.  (WRB)
!  * 900604  DP version created from SP version. .  (RWC)

   subroutine dwnlit (w, mdw, m, n, l, ipivot, itype, h, scale, &
                      rnorm, idope, dope, done)

   integer :: idope(*), ipivot(*), itype(*), l, m, mdw, n
   real(wp) :: dope(*), h(*), rnorm, scale(*), w(mdw,*)
   logical :: done

   real(wp) :: alsq, amax, eanorm, factor, hbar, rn, sparam(5), &
               t, tau
   integer :: i, i1, imax, ir, j, j1, jj, jp, krank, l1, lb, lend, me, &
              mend, niv, nsoln
   logical :: indep, recalc

   me     = idope(1)
   nsoln  = idope(2)
   l1     = idope(3)
   alsq   = dope(1)
   eanorm = dope(2)
   tau    = dope(3)
   lb     = min(m-1,l)
   recalc = .true.
   rnorm  = 0.0_wp
   krank  = 0

   ! We set FACTOR=1.0 so that the heavy weight ALAMDA will be
   ! included in the test for column independence.
   factor = 1.0_wp
   lend = l

   main : block

      do i=1,lb

         ! Set IR to point to the I-th row.
         ir = i
         mend = m
         call dwnlt1 (i, lend, m, ir, mdw, recalc, imax, hbar, h, scale, w)

         ! Update column SS and find pivot column.
         call dwnlt3 (i, imax, m, mdw, ipivot, h, w)

         do
            ! Perform column interchange.
            ! Test independence of incoming column.
            if (dwnlt2(me, mend, ir, factor, tau, scale, w(1,i))) then

               ! Eliminate I-th column below diagonal using modified Givens
               ! transformations applied to (A B).
               !
               ! When operating near the ME line, use the largest element
               ! above it as the pivot.
               do j=m,i+1,-1
                  jp = j-1
                  if (j==me+1) then
                     imax = me
                     amax = scale(me)*w(me,i)**2
                     do jp=j-1,i,-1
                        t = scale(jp)*w(jp,i)**2
                        if (t>amax) then
                           imax = jp
                           amax = t
                        endif
                     end do
                     jp = imax
                  endif
                  if (w(j,i)/=0.0_wp) then
                     call drotmg (scale(jp), scale(j), w(jp,i), w(j,i), &
                                 sparam)
                     w(j,i) = 0.0_wp
                     call drotm (n+1-i, w(jp,i+1), mdw, w(j,i+1), mdw, &
                                 sparam)
                  endif
               end do
               exit
            else if (lend>i) then
               ! Column I is dependent.  Swap with column LEND.
               ! Perform column interchange,
               ! and find column in remaining set with largest SS.
               call dwnlt3 (i, lend, m, mdw, ipivot, h, w)
               lend = lend - 1
               imax = idamax(lend-i+1, h(i), 1) + i - 1
               hbar = h(imax)
            else
               krank = i - 1
               exit main
            endif
         end do

      end do
      krank = l1

   end block main

   if (krank<me) then
      factor = alsq
      do i=krank+1,me
         call dcopy (l, [0.0_wp], 0, w(i,1), mdw)
      end do

      ! Determine the rank of the remaining equality constraint
      ! equations by eliminating within the block of constrained
      ! variables.  Remove any redundant constraints.

      recalc = .true.
      lb = min(l+me-krank, n)
      do i=l+1,lb
         ir = krank + i - l
         lend = n
         mend = me
         call dwnlt1 (i, lend, me, ir, mdw, recalc, imax, hbar, h, &
                      scale, w)

         ! Update col ss and find pivot col
         call dwnlt3 (i, imax, m, mdw, ipivot, h, w)

         ! Perform column interchange
         ! Eliminate elements in the I-th col.
         do j=me,ir+1,-1
            if (w(j,i)/=0.0_wp) then
              call drotmg (scale(j-1), scale(j), w(j-1,i), w(j,i), &
                           sparam)
               w(j,i) = 0.0_wp
               call drotm (n+1-i, w(j-1,i+1), mdw,w(j,i+1), mdw, &
                           sparam)
            endif
         end do

         ! I=column being eliminated.
         ! Test independence of incoming column.
         ! Remove any redundant or dependent equality constraints.
         if (.not.dwnlt2(me, mend, ir, factor,tau,scale,w(1,i))) then
            jj = ir
            do ir=jj,me
               call dcopy (n, [0.0_wp], 0, w(ir,1), mdw)
               rnorm = rnorm + (scale(ir)*w(ir,n+1)/alsq)*w(ir,n+1)
               w(ir,n+1) = 0.0_wp
               scale(ir) = 1.0_wp
               ! Reclassify the zeroed row as a least squares equation.
               itype(ir) = 1
            end do
            ! Reduce ME to reflect any discovered dependent equality
            ! constraints.
            me = jj - 1
            exit
         endif
      end do
   endif

   ! Try to determine the variables KRANK+1 through L1 from the
   ! least squares equations.  Continue the triangularization with
   ! pivot element W(ME+1,I).
   if (krank<l1) then
      recalc = .true.

      ! Set FACTOR=ALSQ to remove effect of heavy weight from
      ! test for column independence.
      factor = alsq
      do i=krank+1,l1

         ! Set IR to point to the ME+1-st row.
         ir = me+1
         lend = l
         mend = m
         call dwnlt1 (i, l, m, ir, mdw, recalc, imax, hbar, h, scale, w)

         ! Update column SS and find pivot column.
         call dwnlt3 (i, imax, m, mdw, ipivot, h, w)

         ! Perform column interchange.
         ! Eliminate I-th column below the IR-th element.
         do j=m,ir+1,-1
            if (w(j,i)/=0.0_wp) then
               call drotmg (scale(j-1), scale(j), w(j-1,i), w(j,i), sparam)
               w(j,i) = 0.0_wp
               call drotm (n+1-i, w(j-1,i+1),  mdw, w(j,i+1), mdw, sparam)
            endif
         end do

         ! Test if new pivot element is near zero.
         ! If so, the column is dependent.
         ! Then check row norm test to be classified as independent.
         t = scale(ir)*w(ir,i)**2
         indep = t > (tau*eanorm)**2
         if (indep) then
            rn = 0.0_wp
            do i1=ir,m
               do j1=i+1,n
                  rn = max(rn, scale(i1)*w(i1,j1)**2)
               end do
            end do
            indep = t > rn*tau**2
         endif

         ! If independent, swap the IR-th and KRANK+1-th rows to
         ! maintain the triangular form.  Update the rank indicator
         ! KRANK and the equality constraint pointer ME.
         if (.not.indep) exit
         call dswap(n+1, w(krank+1,1), mdw, w(ir,1), mdw)
         call dswap(1, scale(krank+1), 1, scale(ir), 1)

         ! Reclassify the least square equation as an equality
         ! constraint and rescale it.
         itype(ir) = 0
         t = sqrt(scale(krank+1))
         call dscal(n+1, t, w(krank+1,1), mdw)
         scale(krank+1) = alsq
         me = me+1
         krank = krank+1
      end do
   endif

   ! If pseudorank is less than L, apply Householder transformation.
   ! from right.
   if (krank<l) then
      do j=krank,1,-1
         call dh12 (1, j, krank+1, l, w(j,1), mdw, h(j), w, mdw, 1, &
                   j-1)
      end do
   endif

   niv = krank + nsoln - l
   if (l==n) done = .true.

   ! End of initial triangularization.
   idope(1) = me
   idope(2) = krank
   idope(3) = niv

   end subroutine dwnlit
!*****************************************************************************************

!*****************************************************************************************
!>
!  This is a companion subprogram to [[DWNNLS]].
!  The documentation for [[DWNNLS]] has complete usage instructions.
!
!  In addition to the parameters discussed in the prologue to
!  subroutine [[DWNNLS]], the following work arrays are used in
!  subroutine [[DWNLSM]]  (they are passed through the calling
!  sequence from [[DWNNLS]] for purposes of variable dimensioning).
!  Their contents will in general be of no interest to the user.
!
!      IPIVOT(*)
!         An array of length N.  Upon completion it contains the
!      pivoting information for the cols of W(*,*).
!
!      ITYPE(*)
!         An array of length M which is used to keep track
!      of the classification of the equations.  ITYPE(I)=0
!      denotes equation I as an equality constraint.
!      ITYPE(I)=1 denotes equation I as a least squares
!      equation.
!
!      WD(*)
!         An array of length N.  Upon completion it contains the
!      dual solution vector.
!
!      H(*)
!         An array of length N.  Upon completion it contains the
!      pivot scalars of the Householder transformations performed
!      in the case KRANK<L.
!
!      SCALE(*)
!         An array of length M which is used by the subroutine
!      to store the diagonal matrix of weights.
!      These are used to apply the modified Givens
!      transformations.
!
!      Z(*),TEMP(*)
!         Working arrays of length N.
!
!      D(*)
!         An array of length N that contains the
!      column scaling for the matrix (E).
!                                    (A)
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and revised.  (WRB & RWC)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900328  Added TYPE section.  (WRB)
!  * 900510  Fixed an error message.  (RWC)
!  * 900604  DP version created from SP version.  (RWC)
!  * 900911  Restriction on value of ALAMDA included.  (WRB)

   subroutine dwnlsm (w, mdw, mme, ma, n, l, prgopt, x, rnorm, mode, &
                      ipivot, itype, wd, h, scale, z, temp, d)

   integer :: ipivot(*), itype(*), l, ma, mdw, mme, mode, n
   real(wp) :: d(*), h(*), prgopt(*), rnorm, scale(*), temp(*), &
               w(mdw,*), wd(*), x(*), z(*)

   real(wp) :: alamda, alpha, alsq, amax, blowup, bnorm, &
               dope(3), eanorm, fac, sm, sparam(5), t, tau, wmax, z2, &
               zz
   integer :: i, idope(3), imax, isol, itemp, iter, itmax, iwmax, j, &
              jcon, jp, key, krank, l1, last, link, m, me, next, niv, nlink, &
              nopt, nsoln, ntimes
   logical :: done, feasbl, hitcon, pos

   ! Set the nominal tolerance used in the code.
   tau = sqrt(drelpr)

   m = ma + mme
   me = mme
   mode = 2

   ! To process option vector
   fac = 1.0e-4_wp

   ! Set the nominal blow up factor used in the code.
   blowup = tau

   ! The nominal column scaling used in the code is
   ! the identity scaling.
   call dcopy (n, [1.0_wp], 0, d, 1)

   ! Define bound for number of options to change.
   nopt = 1000

   ! Define bound for positive value of LINK.
   nlink = 100000
   ntimes = 0
   last = 1
   link = prgopt(1)
   if (link<=0 .or. link>nlink) then
      write(*,*) 'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED'
      return
   endif

   do
      if (link<=1) exit
      ntimes = ntimes + 1
      if (ntimes>nopt) then
         write(*,*) 'IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.'
         return
      endif

      key = prgopt(last+1)
      if (key==6 .and. prgopt(last+2)/=0.0_wp) then
         do j = 1,n
            t = dnrm2(m,w(1,j),1)
            if (t/=0.0_wp) t = 1.0_wp/t
            d(j) = t
         end do
      endif

      if (key==7) call dcopy (n, prgopt(last+2), 1, d, 1)
      if (key==8) tau = max(drelpr,prgopt(last+2))
      if (key==9) blowup = max(drelpr,prgopt(last+2))

      next = prgopt(link)
      if (next<=0 .or. next>nlink) then
         write(*,*) 'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED'
         return
      endif

      last = link
      link = next
   end do

   do j = 1,n
      call dscal (m, d(j), w(1,j), 1)
   end do

   ! Process option vector
   done = .false.
   iter = 0
   itmax = 3*(n-l)
   mode = 0
   nsoln = l
   l1 = min(m,l)

   ! Compute scale factor to apply to equality constraint equations.
   do j = 1,n
      wd(j) = dasum(m,w(1,j),1)
   end do

   imax = idamax(n,wd,1)
   eanorm = wd(imax)
   bnorm = dasum(m,w(1,n+1),1)
   alamda = eanorm/(drelpr*fac)

   ! On machines, such as the VAXes using D floating, with a very
   ! limited exponent range for double precision values, the previously
   ! computed value of ALAMDA may cause an overflow condition.
   ! Therefore, this code further limits the value of ALAMDA.
   alamda = min(alamda,sqrt(huge(1.0_wp)))

   ! Define scaling diagonal matrix for modified Givens usage and
   ! classify equation types.
   alsq = alamda**2
   do i = 1,m
      ! When equation I is heavily weighted ITYPE(I)=0,
      ! else ITYPE(I)=1.
      if (i<=me) then
         t = alsq
         itemp = 0
      else
         t = 1.0_wp
         itemp = 1
      endif
      scale(i) = t
      itype(i) = itemp
   end do

   ! Set the solution vector X(*) to zero and the column interchange
   ! matrix to the identity.
   call dcopy (n, [0.0_wp], 0, x, 1)
   do i = 1,n
      ipivot(i) = i
   end do

   ! Perform initial triangularization in the submatrix
   ! corresponding to the unconstrained variables.
   ! Set first L components of dual vector to zero because
   ! these correspond to the unconstrained variables.
   call dcopy (l, [0.0_wp], 0, wd, 1)

   ! The arrays IDOPE(*) and DOPE(*) are used to pass
   ! information to DWNLIT().  This was done to avoid
   ! a long calling sequence or the use of COMMON.
   idope(1) = me
   idope(2) = nsoln
   idope(3) = l1
   dope(1) = alsq
   dope(2) = eanorm
   dope(3) = tau
   call dwnlit (w, mdw, m, n, l, ipivot, itype, h, scale, rnorm, &
                idope, dope, done)
   me    = idope(1)
   krank = idope(2)
   niv   = idope(3)

   main : do

      ! Perform WNNLS algorithm using the following steps.
      !
      ! Until(DONE)
      !    compute search direction and feasible point
      !    when (HITCON) add constraints
      !    else perform multiplier test and drop a constraint
      !    fin
      ! Compute-Final-Solution
      !
      ! To compute search direction and feasible point,
      ! solve the triangular system of currently non-active
      ! variables and store the solution in Z(*).
      !
      ! To solve system
      ! Copy right hand side into TEMP vector to use overwriting method.

      if (done) exit main
      isol = l + 1
      if (nsoln>=isol) then
         call dcopy (niv, w(1,n+1), 1, temp, 1)
         do j = nsoln,isol,-1
            if (j>krank) then
               i = niv - nsoln + j
            else
               i = j
            endif
            if (j>krank .and. j<=l) then
               z(j) = 0.0_wp
            else
               z(j) = temp(i)/w(i,j)
               call daxpy (i-1, -z(j), w(1,j), 1, temp, 1)
            endif
         end do
      endif

      ! Increment iteration counter and check against maximum number
      ! of iterations.
      iter = iter + 1
      if (iter>itmax) then
         mode = 1
         done = .true.
      endif

      ! Check to see if any constraints have become active.
      ! If so, calculate an interpolation factor so that all
      ! active constraints are removed from the basis.
      alpha = 2.0_wp
      hitcon = .false.
      do j = l+1,nsoln
         zz = z(j)
         if (zz<=0.0_wp) then
            t = x(j)/(x(j)-zz)
            if (t<alpha) then
               alpha = t
               jcon = j
            endif
            hitcon = .true.
         endif
      end do

      ! Compute search direction and feasible point
      if (hitcon) then

         ! To add constraints, use computed ALPHA to interpolate between
         ! last feasible solution X(*) and current unconstrained (and
         ! infeasible) solution Z(*).
         do j = l+1,nsoln
            x(j) = x(j) + alpha*(z(j)-x(j))
         end do
         feasbl = .false.

         do
            ! Remove column JCON and shift columns JCON+1 through N to the
            ! left.  Swap column JCON into the N th position.  This achieves
            ! upper Hessenberg form for the nonactive constraints and
            ! leaves an upper Hessenberg matrix to retriangularize.
            do i = 1,m
               t = w(i,jcon)
               call dcopy (n-jcon, w(i, jcon+1), mdw, w(i, jcon), mdw)
               w(i,n) = t
            end do

            ! Update permuted index vector to reflect this shift and swap.
            itemp = ipivot(jcon)
            do i = jcon,n - 1
               ipivot(i) = ipivot(i+1)
            end do
            ipivot(n) = itemp

            ! Similarly permute X(*) vector.
            call dcopy (n-jcon, x(jcon+1), 1, x(jcon), 1)
            x(n) = 0.0_wp
            nsoln = nsoln - 1
            niv = niv - 1

            ! Retriangularize upper Hessenberg matrix after adding
            ! constraints.
            i = krank + jcon - l
            do j = jcon,nsoln
               if (itype(i)==0 .and. itype(i+1)==0) then
                  ! Zero IP1 to I in column J
                  if (w(i+1,j)/=0.0_wp) then
                     call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j), &
                                 sparam)
                     w(i+1,j) = 0.0_wp
                     call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw, &
                                 sparam)
                  endif
               elseif (itype(i)==1 .and. itype(i+1)==1) then
                  ! Zero IP1 to I in column J
                  if (w(i+1,j)/=0.0_wp) then
                     call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j), &
                                 sparam)
                     w(i+1,j) = 0.0_wp
                     call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw, &
                                 sparam)
                  endif
               elseif (itype(i)==1 .and. itype(i+1)==0) then
                  call dswap (n+1, w(i,1), mdw, w(i+1,1), mdw)
                  call dswap (1, scale(i), 1, scale(i+1), 1)
                  itemp = itype(i+1)
                  itype(i+1) = itype(i)
                  itype(i) = itemp
                  ! Swapped row was formerly a pivot element, so it will
                  ! be large enough to perform elimination.
                  ! Zero IP1 to I in column J.
                  if (w(i+1,j)/=0.0_wp) then
                     call drotmg (scale(i), scale(i+1), w(i,j), w(i+1,j), &
                                 sparam)
                     w(i+1,j) = 0.0_wp
                     call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw, &
                                 sparam)
                  endif
               elseif (itype(i)==0 .and. itype(i+1)==1) then
                  if (scale(i)*w(i,j)**2/alsq>(tau*eanorm)**2) then
                     ! Zero IP1 to I in column J
                     if (w(i+1,j)/=0.0_wp) then
                        call drotmg (scale(i), scale(i+1), w(i,j), &
                                    w(i+1,j), sparam)
                        w(i+1,j) = 0.0_wp
                        call drotm (n+1-j, w(i,j+1), mdw, w(i+1,j+1), mdw, &
                                    sparam)
                     endif
                  else
                     call dswap (n+1, w(i,1), mdw, w(i+1,1), mdw)
                     call dswap (1, scale(i), 1, scale(i+1), 1)
                     itemp = itype(i+1)
                     itype(i+1) = itype(i)
                     itype(i) = itemp
                     w(i+1,j) = 0.0_wp
                  endif
               endif
               i = i + 1
            end do

            ! See if the remaining coefficients in the solution set are
            ! feasible.  They should be because of the way ALPHA was
            ! determined.  If any are infeasible, it is due to roundoff
            ! error.  Any that are non-positive will be set to zero and
            ! removed from the solution set.
            do jcon = l+1,nsoln
               if (x(jcon)<=0.0_wp) then
                  exit
               else
                  if (jcon==nsoln) feasbl = .true.
               end if
            end do
            if (feasbl) exit
         end do

      else

         ! To perform multiplier test and drop a constraint.
         call dcopy (nsoln, z, 1, x, 1)
         if (nsoln<n) call dcopy (n-nsoln, [0.0_wp], 0, x(nsoln+1), 1)

         ! Reclassify least squares equations as equalities as necessary.
         i = niv + 1
         do
            if (i>me) exit
            if (itype(i)==0) then
               i = i + 1
            else
               call dswap (n+1, w(i,1), mdw, w(me,1), mdw)
               call dswap (1, scale(i), 1, scale(me), 1)
               itemp = itype(i)
               itype(i) = itype(me)
               itype(me) = itemp
               me = me - 1
            endif
         end do

         ! Form inner product vector WD(*) of dual coefficients.
         do j = nsoln+1,n
            sm = 0.0_wp
            do i = nsoln+1,m
               sm = sm + scale(i)*w(i,j)*w(i,n+1)
            end do
            wd(j) = sm
         end do

         do
            ! Find J such that WD(J)=WMAX is maximum.  This determines
            ! that the incoming column J will reduce the residual vector
            ! and be positive.
            wmax = 0.0_wp
            iwmax = nsoln + 1
            do j = nsoln+1,n
               if (wd(j)>wmax) then
                  wmax = wd(j)
                  iwmax = j
               endif
            end do
            if (wmax<=0.0_wp) exit main

            ! Set dual coefficients to zero for incoming column.
            wd(iwmax) = 0.0_wp

            ! WMAX > 0.0_wp, so okay to move column IWMAX to solution set.
            ! Perform transformation to retriangularize, and test for near
            ! linear dependence.
            !
            ! Swap column IWMAX into NSOLN-th position to maintain upper
            ! Hessenberg form of adjacent columns, and add new column to
            ! triangular decomposition.
            nsoln = nsoln + 1
            niv = niv + 1
            if (nsoln/=iwmax) then
               call dswap (m, w(1,nsoln), 1, w(1,iwmax), 1)
               wd(iwmax) = wd(nsoln)
               wd(nsoln) = 0.0_wp
               itemp = ipivot(nsoln)
               ipivot(nsoln) = ipivot(iwmax)
               ipivot(iwmax) = itemp
            endif

            ! Reduce column NSOLN so that the matrix of nonactive constraints
            ! variables is triangular.
            do j = m,niv+1,-1
               jp = j - 1

               ! When operating near the ME line, test to see if the pivot
               ! element is near zero.  If so, use the largest element above
               ! it as the pivot.  This is to maintain the sharp interface
               ! between weighted and non-weighted rows in all cases.
               if (j==me+1) then
                  imax = me
                  amax = scale(me)*w(me,nsoln)**2
                  do jp = j - 1,niv,-1
                     t = scale(jp)*w(jp,nsoln)**2
                     if (t>amax) then
                        imax = jp
                        amax = t
                     endif
                  end do
                  jp = imax
               endif

               if (w(j,nsoln)/=0.0_wp) then
                  call drotmg (scale(jp), scale(j), w(jp,nsoln), w(j,nsoln), sparam)
                  w(j,nsoln) = 0.0_wp
                  call drotm (n+1-nsoln, w(jp,nsoln+1), mdw, w(j,nsoln+1), mdw, sparam)
               endif
            end do

            ! Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
            ! this is nonpositive or too large.  If this was true or if the
            ! pivot term was zero, reject the column as dependent.
            if (w(niv,nsoln)/=0.0_wp) then
               isol = niv
               z2 = w(isol,n+1)/w(isol,nsoln)
               z(nsoln) = z2
               pos = z2 > 0.0_wp
               if (z2*eanorm>=bnorm .and. pos) then
                  pos = .not. (blowup*z2*eanorm>=bnorm)
               endif

            elseif (niv<=me .and. w(me+1,nsoln)/=0.0_wp) then
               ! Try to add row ME+1 as an additional equality constraint.
               ! Check size of proposed new solution component.
               ! Reject it if it is too large.
               isol = me + 1
               if (pos) then
                  ! Swap rows ME+1 and NIV, and scale factors for these rows.
                  call dswap (n+1, w(me+1,1), mdw, w(niv,1), mdw)
                  call dswap (1, scale(me+1), 1, scale(niv), 1)
                  itemp = itype(me+1)
                  itype(me+1) = itype(niv)
                  itype(niv) = itemp
                  me = me + 1
               endif
            else
               pos = .false.
            endif

            if (.not.pos) then
               nsoln = nsoln - 1
               niv = niv - 1
            endif
            if (pos.or.done) exit
         end do

      endif

   end do main

   ! Else perform multiplier test and drop a constraint.  To compute
   ! final solution.  Solve system, store results in X(*).
   !
   ! Copy right hand side into TEMP vector to use overwriting method.
   isol = 1
   if (nsoln>=isol) then
      call dcopy (niv, w(1,n+1), 1, temp, 1)
      do j = nsoln,isol,-1
         if (j>krank) then
            i = niv - nsoln + j
         else
            i = j
         endif
         if (j>krank .and. j<=l) then
            z(j) = 0.0_wp
         else
            z(j) = temp(i)/w(i,j)
            call daxpy (i-1, -z(j), w(1,j), 1, temp, 1)
         endif
      end do
   endif

   ! Solve system.
   call dcopy (nsoln, z, 1, x, 1)

   ! Apply Householder transformations to X(*) if KRANK<L
   if (krank<l) then
      do i = 1,krank
         call dh12 (2, i, krank+1, l, w(i,1), mdw, h(i), x, 1, 1, 1)
      end do
   endif

   ! Fill in trailing zeroes for constrained variables not in solution.
   if (nsoln<n) call dcopy (n-nsoln, [0.0_wp], 0, x(nsoln+1), 1)

   ! Permute solution vector to natural order.
   do i = 1,n
      j = i
      do
         if (ipivot(j)==i) exit
         j = j + 1
      end do
      ipivot(j) = ipivot(i)
      ipivot(i) = j
      call dswap (1, x(j), 1, x(i), 1)
   end do

   ! Rescale the solution using the column scaling.
   do j = 1,n
      x(j) = x(j)*d(j)
   end do

   do i = nsoln+1,m
      t = w(i,n+1)
      if (i<=me) t = t/alamda
      t = (scale(i)*t)*t
      rnorm = rnorm + t
   end do

   rnorm = sqrt(rnorm)

   end subroutine dwnlsm
!*****************************************************************************************

!*****************************************************************************************
!>
!  To update the column Sum Of Squares and find the pivot column.
!  The column Sum of Squares Vector will be updated at each step.
!  When numerically necessary, these values will be recomputed.
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!  * 900604  DP version created from SP version.  (RWC)

   subroutine dwnlt1 (i, lend, mend, ir, mdw, recalc, imax, hbar, h, &
                      scale, w)

   integer :: i, imax, ir, lend, mdw, mend
   real(wp) :: h(*), hbar, scale(*), w(mdw,*)
   logical :: recalc

   integer :: j, k

   if (ir/=1 .and. (.not.recalc)) then
      ! Update column SS=sum of squares.
      do j=i,lend
         h(j) = h(j) - scale(ir-1)*w(ir-1,j)**2
      end do
      ! Test for numerical accuracy.
      imax = idamax(lend-i+1, h(i), 1) + i - 1
      recalc = (hbar+1.e-3*h(imax)) == hbar
   endif

   ! If required, recalculate column SS, using rows IR through MEND.
   if (recalc) then
      do j=i,lend
         h(j) = 0.0_wp
         do k=ir,mend
            h(j) = h(j) + scale(k)*w(k,j)**2
         end do
      end do
      ! Find column with largest SS.
      imax = idamax(lend-i+1, h(i), 1) + i - 1
      hbar = h(imax)
   endif

   end subroutine dwnlt1
!*****************************************************************************************

!*****************************************************************************************
!>
!  To test independence of incoming column.
!
!  Test the column IC to determine if it is linearly independent
!  of the columns already in the basis.  In the initial tri. step,
!  we usually want the heavy weight ALAMDA to be included in the
!  test for independence.  In this case, the value of `FACTOR` will
!  have been set to 1.0 before this procedure is invoked.
!  In the potentially rank deficient problem, the value of FACTOR
!  will have been set to `ALSQ=ALAMDA**2` to remove the effect of the
!  heavy weight from the test for independence.
!
!  Write new column as partitioned vector
!
!   * `(A1)`  number of components in solution so far `= NIV`
!   * `(A2)`  `M-NIV` components
!
!  And compute
!
!   * `SN` = inverse weighted length of `A1`
!   * `RN` = inverse weighted length of `A2`
!
!  Call the column independent when `RN > TAU*SN`
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!  * 900604  DP version created from SP version.  (RWC)

   logical function dwnlt2 (me, mend, ir, factor, tau, scale, wic)

   real(wp) :: factor, scale(*), tau, wic(*)
   integer :: ir, me, mend

   real(wp) :: rn, sn, t
   integer :: j

   sn = 0.0_wp
   rn = 0.0_wp
   do j=1,mend
      t = scale(j)
      if (j<=me) t = t/factor
      t = t*wic(j)**2
      if (j<ir) then
         sn = sn + t
      else
         rn = rn + t
      endif
   end do
   dwnlt2 = rn > sn*tau**2

   end function dwnlt2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Perform column interchange.
!  Exchange elements of permuted index vector and perform column
!  interchanges.
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890620  Code extracted from WNLIT and made a subroutine.  (RWC))
!  * 900604  DP version created from SP version.  (RWC)

   subroutine dwnlt3 (i, imax, m, mdw, ipivot, h, w)

   integer,intent(in) :: i
   integer,intent(in) :: imax
   integer,intent(inout) :: ipivot(*)
   integer ,intent(in):: m
   integer,intent(in) :: mdw
   real(wp),intent(inout) :: h(*)
   real(wp),intent(inout) :: w(mdw,*)

   real(wp) :: t
   integer :: itemp

   if (imax/=i) then
      itemp        = ipivot(i)
      ipivot(i)    = ipivot(imax)
      ipivot(imax) = itemp
      call dswap(m, w(1,imax), 1, w(1,i), 1)
      t       = h(imax)
      h(imax) = h(i)
      h(i)    = t
   endif

   end subroutine dwnlt3
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subprogram solves a linearly constrained least squares
!  problem.  Suppose there are given matrices `E` and `A` of
!  respective dimensions `ME` by `N` and `MA` by `N`, and vectors `F`
!  and `B` of respective lengths `ME` and `MA`.  This subroutine
!  solves the problem
!
!   * `EX = F`, (equations to be exactly satisfied)
!   * `AX = B`, (equations to be approximately satisfied, in the least squares sense)
!
!  subject to components `L+1,...,N` nonnegative
!
!  Any values `ME>=0`, `MA>=0` and `0<= L <=N` are permitted.
!
!  The problem is reposed as problem [[DWNNLS]]
!
!```
!            (WT*E)X = (WT*F)
!            (   A)    (   B), (least squares)
!            subject to components L+1,...,N nonnegative.
!```
!
!  The subprogram chooses the heavy weight (or penalty parameter) `WT`.
!
!  The parameters for [[DWNNLS]] are
!
!```
!  INPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!  W(*,*),MDW,  The array W(*,*) is double subscripted with first
!  ME,MA,N,L    dimensioning parameter equal to MDW.  For this
!               discussion let us call M = ME + MA.  Then MDW
!               must satisfy MDW>=M.  The condition MDW<M
!               is an error.
!
!               The array W(*,*) contains the matrices and vectors
!
!                    (E  F)
!                    (A  B)
!
!               in rows and columns 1,...,M and 1,...,N+1
!               respectively.  Columns 1,...,L correspond to
!               unconstrained variables X(1),...,X(L).  The
!               remaining variables are constrained to be
!               nonnegative. The condition L<0 or L>N is
!               an error.
!
!  PRGOPT(*)    This double precision array is the option vector.
!               If the user is satisfied with the nominal
!               subprogram features set
!
!               PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!               Otherwise PRGOPT(*) is a linked list consisting of
!               groups of data of the following form
!
!               LINK
!               KEY
!               DATA SET
!
!               The parameters LINK and KEY are each one word.
!               The DATA SET can be comprised of several words.
!               The number of items depends on the value of KEY.
!               The value of LINK points to the first
!               entry of the next group of data within
!               PRGOPT(*).  The exception is when there are
!               no more options to change.  In that
!               case LINK=1 and the values KEY and DATA SET
!               are not referenced. The general layout of
!               PRGOPT(*) is as follows.
!
!            ...PRGOPT(1)=LINK1 (link to first entry of next group)
!            .  PRGOPT(2)=KEY1 (key to the option change)
!            .  PRGOPT(3)=DATA VALUE (data value for this change)
!            .       .
!            .       .
!            .       .
!            ...PRGOPT(LINK1)=LINK2 (link to the first entry of
!            .                       next group)
!            .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
!            .  PRGOPT(LINK1+2)=DATA VALUE
!            ...     .
!            .       .
!            .       .
!            ...PRGOPT(LINK)=1 (no more options to change)
!
!               Values of LINK that are nonpositive are errors.
!               A value of LINK>NLINK=100000 is also an error.
!               This helps prevent using invalid but positive
!               values of LINK that will probably extend
!               beyond the program limits of PRGOPT(*).
!               Unrecognized values of KEY are ignored.  The
!               order of the options is arbitrary and any number
!               of options can be changed with the following
!               restriction.  To prevent cycling in the
!               processing of the option array a count of the
!               number of options changed is maintained.
!               Whenever this count exceeds NOPT=1000 an error
!               message is printed and the subprogram returns.
!
!               OPTIONS..
!
!               KEY=6
!                      Scale the nonzero columns of the
!               entire data matrix
!               (E)
!               (A)
!               to have length one. The DATA SET for
!               this option is a single value.  It must
!               be nonzero if unit length column scaling is
!               desired.
!
!               KEY=7
!                      Scale columns of the entire data matrix
!               (E)
!               (A)
!               with a user-provided diagonal matrix.
!               The DATA SET for this option consists
!               of the N diagonal scaling factors, one for
!               each matrix column.
!
!               KEY=8
!                      Change the rank determination tolerance from
!               the nominal value of SQRT(SRELPR).  This quantity
!               can be no smaller than SRELPR, The arithmetic-
!               storage precision.  The quantity used
!               here is internally restricted to be at
!               least SRELPR.  The DATA SET for this option
!               is the new tolerance.
!
!               KEY=9
!                      Change the blow-up parameter from the
!               nominal value of SQRT(SRELPR).  The reciprocal of
!               this parameter is used in rejecting solution
!               components as too large when a variable is
!               first brought into the active set.  Too large
!               means that the proposed component times the
!               reciprocal of the parameter is not less than
!               the ratio of the norms of the right-side
!               vector and the data matrix.
!               This parameter can be no smaller than SRELPR,
!               the arithmetic-storage precision.
!
!               For example, suppose we want to provide
!               a diagonal matrix to scale the problem
!               matrix and change the tolerance used for
!               determining linear dependence of dropped col
!               vectors.  For these options the dimensions of
!               PRGOPT(*) must be at least N+6.  The FORTRAN
!               statements defining these options would
!               be as follows.
!
!               PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
!               PRGOPT(2)=7 (user-provided scaling key)
!
!               CALL DCOPY(N,D,1,PRGOPT(3),1) (copy the N
!               scaling factors from a user array called D(*)
!               into PRGOPT(3)-PRGOPT(N+2))
!
!               PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
!               PRGOPT(N+4)=8 (linear dependence tolerance key)
!               PRGOPT(N+5)=... (new value of the tolerance)
!
!               PRGOPT(N+6)=1 (no more options to change)
!
!
!  IWORK(1),    The amounts of working storage actually allocated
!  IWORK(2)     for the working arrays WORK(*) and IWORK(*),
!               respectively.  These quantities are compared with
!               the actual amounts of storage needed for DWNNLS( ).
!               Insufficient storage allocated for either WORK(*)
!               or IWORK(*) is considered an error.  This feature
!               was included in DWNNLS( ) because miscalculating
!               the storage formulas for WORK(*) and IWORK(*)
!               might very well lead to subtle and hard-to-find
!               execution errors.
!
!               The length of WORK(*) must be at least
!
!               LW = ME+MA+5*N
!               This test will not be made if IWORK(1)<=0.
!
!               The length of IWORK(*) must be at least
!
!               LIW = ME+MA+N
!               This test will not be made if IWORK(2)<=0.
!
!  OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!  X(*)         An array dimensioned at least N, which will
!               contain the N components of the solution vector
!               on output.
!
!  RNORM        The residual norm of the solution.  The value of
!               RNORM contains the residual vector length of the
!               equality constraints and least squares equations.
!
!  MODE         The value of MODE indicates the success or failure
!               of the subprogram.
!
!               MODE = 0  Subprogram completed successfully.
!
!                    = 1  Max. number of iterations (equal to
!                         3*(N-L)) exceeded. Nearly all problems
!                         should complete in fewer than this
!                         number of iterations. An approximate
!                         solution and its corresponding residual
!                         vector length are in X(*) and RNORM.
!
!                    = 2  Usage error occurred.  The offending
!                         condition is noted with the error
!                         processing subprogram, XERMSG( ).
!
!  User-designated
!  Working arrays..
!
!  WORK(*)      A double precision working array of length at least
!               M + 5*N.
!
!  IWORK(*)     An integer-valued working array of length at least
!               M+N.
!```
!
!### References
!  * K. H. Haskell and R. J. Hanson, An algorithm for
!    linear least squares problems with equality and
!    nonnegativity constraints, Report SAND77-0552, Sandia
!    Laboratories, June 1978.
!  * K. H. Haskell and R. J. Hanson, Selected algorithms for
!    the linearly constrained least squares problem - a
!    users guide, Report SAND78-1290, Sandia Laboratories,
!    August 1979.
!  * K. H. Haskell and R. J. Hanson, An algorithm for
!    linear least squares problems with equality and
!    nonnegativity constraints, Mathematical Programming
!    21 (1981), pp. 98-118.
!  * R. J. Hanson and K. H. Haskell, Two algorithms for the
!    linearly constrained least squares problem, ACM
!    Transactions on Mathematical Software, September 1982.
!  * C. L. Lawson and R. J. Hanson, Solving Least Squares
!    Problems, Prentice-Hall, Inc., 1974.
!
!### Revision history
!  * 790701  DATE WRITTEN. Hanson, R. J., (SNLA), Haskell, K. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890618  Completely restructured and revised.  (WRB & RWC)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900510  Convert XERRWV calls to XERMSG calls, change Prologue
!           comments to agree with WNNLS.  (RWC)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dwnnls (w, mdw, me, ma, n, l, prgopt, x, rnorm, mode, &
                      iwork, work)

   integer :: iwork(*), l, l1, l2, l3, l4, l5, liw, lw, ma, mdw, me, &
              mode, n
   real(wp) :: prgopt(*), rnorm, w(mdw,*), work(*), x(*)
   character(len=8) :: xern1

   mode = 0
   if (ma+me<=0 .or. n<=0) return

   if (iwork(1)>0) then
      lw = me + ma + 5*n
      if (iwork(1)<lw) then
         write (xern1, '(I8)') lw
         write(*,*) 'INSUFFICIENT STORAGE ' // &
                    'ALLOCATED FOR WORK(*), NEED LW = ' // xern1
         mode = 2
         return
      endif
   endif

   if (iwork(2)>0) then
      liw = me + ma + n
      if (iwork(2)<liw) then
         write (xern1, '(I8)') liw
         write(*,*) 'INSUFFICIENT STORAGE ' // &
                    'ALLOCATED FOR IWORK(*), NEED LIW = ' // xern1
         mode = 2
         return
      endif
   endif

   if (mdw<me+ma) then
      write(*,*) 'THE VALUE MDW<ME+MA IS AN ERROR'
      mode = 2
      return
   endif

   if (l<0 .or. l>n) then
      write(*,*) 'L>=0 .AND. L<=N IS REQUIRED'
      mode = 2
      return
   endif

   ! THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
   ! WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
   ! REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).

   l1 = n + 1
   l2 = l1 + n
   l3 = l2 + me + ma
   l4 = l3 + n
   l5 = l4 + n

   call dwnlsm(w, mdw, me, ma, n, l, prgopt, x, rnorm, mode, iwork, &
               iwork(l1), work(1), work(l1), work(l2), work(l3), &
               work(l4), work(l5))

   end subroutine dwnnls
!*****************************************************************************************

!*****************************************************************************************
!>
!  [[dcv]] is a companion function subprogram for [[dfc]].  The
!  documentation for [[dfc]] has complete usage instructions.
!
!  [[dcv]] is used to evaluate the variance function of the curve
!  obtained by the constrained B-spline fitting subprogram, [[dfc]].
!  The variance function defines the square of the probable error
!  of the fitted curve at any point, XVAL.  One can use the square
!  root of this variance function to determine a probable error band
!  around the fitted curve.
!
!  [[dcv]] is used after a call to [[dfc]].  MODE, an input variable to
!  [[dfc]], is used to indicate if the variance function is desired.
!  In order to use [[dcv]], MODE must equal 2 or 4 on input to [[dfc]].
!  MODE is also used as an output flag from [[dfc]].  Check to make
!  sure that MODE = 0 after calling [[dfc]], indicating a successful
!  constrained curve fit.  The array SDDATA, as input to [[dfc]], must
!  also be defined with the standard deviation or uncertainty of the
!  Y values to use [[dcv]].
!
!  To evaluate the variance function after calling [[dfc]] as stated
!  above, use [[dcv]] as shown here
!
!  `VAR = DCV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W)`
!
!  The variance function is given by
!
!  `VAR = (transpose of B(XVAL))*C*B(XVAL)/DBLE(MAX(NDATA-N,1))`
!
!  where `N = NBKPT - NORD`.
!
!  The vector B(XVAL) is the B-spline basis function values at
!  X=XVAL.  The covariance matrix, C, of the solution coefficients
!  accounts only for the least squares equations and the explicitly
!  stated equality constraints.  This fact must be considered when
!  interpreting the variance function from a data fitting problem
!  that has inequality constraints on the fitted curve.
!
!  All the variables in the calling sequence for [[dcv]] are used in
!  [[dfc]] except the variable XVAL.  Do not change the values of
!  these variables between the call to [[dfc]] and the use of [[dcv]].
!
!### Reference
!  * R. J. Hanson, Constrained least squares curve fitting
!    to discrete data using B-splines, a users guide,
!    Report SAND78-1291, Sandia Laboratories, December
!    1978.
!
!### Revision history
!  * 780801  DATE WRITTEN. Hanson, R. J., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 890911  Removed unnecessary intrinsics.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   real(wp) function dcv (xval, ndata, nconst, nord, nbkpt, bkpt, w)

      real(wp),intent(in) :: xval !! The point where the variance is desired
      integer,intent(in) :: nbkpt !! The number of knots in the array BKPT(*).
                                  !! The value of NBKPT must satisfy NBKPT .GE. 2*NORD.
      integer,intent(in) :: nconst !! The number of conditions that constrained the B-spline in
                                   !! [[dfc]].
      integer,intent(in) :: ndata !! The number of discrete (X,Y) pairs for which [[dfc]]
                                  !! calculated a piece-wise polynomial curve.
      integer,intent(in) :: nord !! The order of the B-spline used in [[dfc]].
                                 !! The value of NORD must satisfy 1 < NORD < 20 .
                                 !!
                                 !! (The order of the spline is one more than the degree of
                                 !! the piece-wise polynomial defined on each interval.  This
                                 !! is consistent with the B-spline package convention.  For
                                 !! example, NORD=4 when we are using piece-wise cubics.)
      real(wp),intent(in) :: bkpt(*) !! The array of knots.  Normally the problem
                                     !! data interval will be included between the limits
                                     !! BKPT(NORD) and BKPT(NBKPT-NORD+1).  The additional end
                                     !! knots BKPT(I),I=1,...,NORD-1 and I=NBKPT-NORD+2,...,NBKPT,
                                     !! are required by [[dfc]] to compute the functions used to
                                     !! fit the data.
      real(wp) :: w(*) !! Real work array as used in [[dfc]].  See [[dfc]]
                       !! for the required length of W(*).  The contents of W(*)
                       !! must not be modified by the user if the variance function
                       !! is desired.

      real(wp) :: v(40)
      integer :: i, ileft, ip, is, last, mdg, mdw, n
      integer :: dfspvn_j
      real(wp), dimension(20) :: dfspvn_deltam, dfspvn_deltap

      real(wp),parameter :: zero = 0.0_wp

      ! set up variables for dfspvn
      dfspvn_j = 1
      dfspvn_deltam = zero
      dfspvn_deltap = zero

      mdg = nbkpt - nord + 3
      mdw = nbkpt - nord + 1 + nconst
      is = mdg*(nord + 1) + 2*max(ndata,nbkpt) + nbkpt + nord**2
      last = nbkpt - nord + 1
      ileft = nord
      do
         if (xval < bkpt(ileft+1) .or. ileft >= last - 1) exit
         ileft = ileft + 1
      end do
      call dfspvn(bkpt,nord,1,xval,ileft,v(nord+1), &
                  dfspvn_j, dfspvn_deltam, dfspvn_deltap)
      ileft = ileft - nord + 1
      ip = mdw*(ileft - 1) + ileft + is
      n = nbkpt - nord
      do i = 1, nord
         v(i) = ddot(nord,w(ip),1,v(nord+1),1)
         ip = ip + mdw
      end do
      dcv = max(ddot(nord,v,1,v(nord+1),1),zero)

      ! scale the variance so it is an unbiased estimate.
      dcv = dcv/max(ndata-n,1)

   end function dcv
!*****************************************************************************************

!*****************************************************************************************
   end module bspline_defc_module
!*****************************************************************************************