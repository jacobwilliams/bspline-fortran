!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!
!  [[DEFC]] procedure and support routines from SLATEC.

module bspline_defc_module

   use bspline_kinds_module, only: wp !, ip
   use bspline_blas_module

   implicit none

   private

   public :: defc

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
!  with subprogram DFC( ) on a very large data set.
!
!### References
!
!  1. C. W. de Boor, Package for Calculating with B-Splines.
!     SIAM J. Numer. Anal., p. 441, (June, 1977).
!  2. R. J. Hanson, Constrained least squares curve fitting
!     to discrete data using B-splines, a users guide,
!     Report SAND78-1291, Sandia Laboratories, December
!     1978.
!
!### Notes
!
!  * For a description of the B-splines and usage instructions to
!    evaluate them, see reference 1.
!  * For further discussion of (constrained) curve fitting using
!    B-splines, see reference 2.
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
                                     !! Internal to `DEFC( )` the extreme end knots may
                                     !! be reduced and increased respectively to
                                     !! accommodate any data values that are exterior
                                     !! to the given knot values.  The contents of
                                     !! `BKPT(*)` is not changed.
      integer,intent(in) :: Mdein !! An integer flag, with one of two possible
                                  !! values (1 or 2), that directs the subprogram
                                  !! action with regard to new data points provided
                                  !! by the user:
                                  !!
                                  !! * `= 1`  The first time that DEFC( ) has been
                                  !!   entered.  There are NDATA points to process.
                                  !! * `= 2`  This is another entry to DEFC().  The
                                  !!   subprogram DEFC( ) has been entered with MDEIN=1
                                  !!   exactly once before for this problem.  There
                                  !!   are NDATA new additional points to merge and
                                  !!   process with any previous points.
                                  !!   (When using DEFC( ) with MDEIN=2 it is
                                  !!   important that the set of knots remain fixed at the
                                  !!   same values for all entries to DEFC( ).)
      integer,intent(out) :: Mdeout !! An output flag that indicates the status
                                    !! of the curve fit:
                                    !!
                                    !!  * `=-1`  A usage error of `DEFC( )` occurred.  The
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
                                    !!    subprogram `DFC( )` to obtain a specific
                                    !!    set of coefficients.  The user should read
                                    !!    the usage instructions for `DFC( )` for further
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
                                       !! curve returned by `DEFC( )`, the constrained
                                       !! least squares curve fitting subprogram `DFC( )`
                                       !! may be required.  The work done within `DEFC( )`
                                       !! to accumulate the data can be utilized by
                                       !! the user, if so desired.  This involves
                                       !! saving the first `(NBKPT-NORD+3)*(NORD+1)`
                                       !! entries of `W(*)` and providing this data
                                       !! to `DFC( )` with the "old problem" designation.
                                       !! The user should read the usage instructions
                                       !! for subprogram `DFC( )` for further details.
      integer,intent(in) :: Lw !! The amount of working storage actually
                               !! allocated for the working array `W(*)`.
                               !! This quantity is compared with the
                               !! actual amount of storage needed in `DEFC( )`.
                               !! Insufficient storage allocated for `W(*)` is
                               !! an error.  This feature was included in `DEFC`
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
                       !! to `DEFC( )` with values of `MDEIN=1,2,2,...` .
                       !! The first `(NBKPT-NORD+3)*(NORD+1)` entries of
                       !! `W(*)` are acceptable as direct input to `DFC( )`
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

!
!     Initialize variables and analyze input.
!
      n = Nbkpt - Nord
      np1 = n + 1
!
!     Initially set all output coefficients to zero.
!
      call dcopy(n, [0.0_wp], 0, Coeff, 1)

      !.... JW : add iflag outputs to be consistent with the rest of the library... start with 4000
      Mdeout = -1
      if (Nord < 1 .or. Nord > 20) then
         write (*, *) 'IN DEFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.'
         return
      end if
!
      if (Nbkpt < 2*Nord) then
         write (*, *) 'IN DEFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.'
         return
      end if
!
      if (Ndata < 0) then
         write (*, *) 'IN DEFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.'
         return
      end if
!
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
!
      if (Mdein /= 1 .and. Mdein /= 2) then
         write (*, *) 'IN DEFC, INPUT VALUE OF MDEIN MUST BE 1-2.'
         return
      end if
!
!     Sort the breakpoints.
!
      call dcopy(Nbkpt, Bkptin, 1, Bkpt, 1)
      call dsort(Nbkpt, 1, Bkpt)
!
!     Save interval containing knots.
!
      xmin = Bkpt(Nord)
      xmax = Bkpt(np1)
      nordm1 = Nord - 1
      nordp1 = Nord + 1
!
!     Process least squares equations.
!
!     Sort data and an array of pointers.
!
      call dcopy(Ndata, Xdata, 1, Xtemp, 1)
      do i = 1, Ndata
         Ptemp(i) = i
      end do
      ! JW : really Ptemp should be an integer array.
      !      it's real because they are stuffing it in
      !      a real work array and also using dsort on it.
!
      if (Ndata > 0) then
         call dsort(Ndata, 2, Xtemp, Ptemp)
         xmin = min(xmin, Xtemp(1))
         xmax = max(xmax, Xtemp(Ndata))
      end if
!
!     Fix breakpoint array if needed. This should only involve very
!     minor differences with the input array of breakpoints.
!
      do i = 1, Nord
         Bkpt(i) = min(Bkpt(i), xmin)
      end do
!
      do i = np1, Nbkpt
         Bkpt(i) = max(Bkpt(i), xmax)
      end do

      ! set up variables for dfspvn
      dfspvn_j = 1
      dfspvn_deltam = 0.0_wp
      dfspvn_deltap = 0.0_wp

!
!     Initialize parameters of banded matrix processor, DBNDAC( ).
!
      mt = 0
      ip = 1
      ir = 1
      ileft = Nord
      intseq = 1
      do idata = 1, Ndata
!
!        Sorted indices are in PTEMP(*).
!
         l = int(Ptemp(idata))
         xval = Xdata(l)
!
!        When interval changes, process equations in the last block.
!
         if (xval >= Bkpt(ileft + 1)) then
            call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)
            mt = 0
!
!           Move pointer up to have BKPT(ILEFT)<=XVAL, ILEFT<=N.
!
            do ileft = ileft, n
               if (xval < Bkpt(ileft + 1)) exit
               if (Mdein == 2) then
!
!                 Data is being sequentially accumulated.
!                 Transfer previously accumulated rows from W(*,*) to
!                 G(*,*) and process them.
!
                  call dcopy(nordp1, w(intseq, 1), Mdw, g(ir, 1), Mdg)
                  call dbndac(g, Mdg, Nord, ip, ir, 1, intseq)
                  intseq = intseq + 1
               end if
            end do
         end if
!
!        Obtain B-spline function value.
!
         call dfspvn(Bkpt, Nord, 1, xval, ileft, Bf, &
                     dfspvn_j, dfspvn_deltam, dfspvn_deltap)
!
!        Move row into place.
!
         irow = ir + mt
         mt = mt + 1
         call dcopy(Nord, Bf, 1, g(irow, 1), Mdg)
         g(irow, nordp1) = Ydata(l)
!
!        Scale data if uncertainty is nonzero.
!
         if (Sddata(l) /= 0.0_wp) call dscal(nordp1, 1.0_wp/Sddata(l), g(irow, 1), Mdg)
!
!        When staging work area is exhausted, process rows.
!
         if (irow == Mdg - 1) then
            call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)
            mt = 0
         end if
      end do
!
!     Process last block of equations.
!
      call dbndac(g, Mdg, Nord, ip, ir, mt, ileft - nordm1)
!
!     Finish processing any previously accumulated rows from W(*,*)
!     to G(*,*).
!
      if (Mdein == 2) then
         do i = intseq, np1
            call dcopy(nordp1, w(i, 1), Mdw, g(ir, 1), Mdg)
            call dbndac(g, Mdg, Nord, ip, ir, 1, min(n, i))
         end do
      end if
!
!     Last call to adjust block positioning.
!
      call dcopy(nordp1, [0.0_wp], 0, g(ir, 1), Mdg)
      call dbndac(g, Mdg, Nord, ip, ir, 1, np1)
!
!     Transfer accumulated rows from G(*,*) to W(*,*) for
!     possible later sequential accumulation.
!
      do i = 1, np1
         call dcopy(nordp1, g(i, 1), Mdg, w(i, 1), Mdw)
      end do
!
!     Solve for coefficients when possible.
!
      do i = 1, n
         if (g(i, 1) == 0.0_wp) then
            Mdeout = 2
            return
         end if
      end do
!
!     All the diagonal terms in the accumulated triangular
!     matrix are nonzero.  The solution can be computed but
!     it may be unsuitable for further use due to poor
!     conditioning or the lack of constraints.  No checking
!     for either of these is done here.
!
      call dbndsl(1, g, Mdg, Nord, ip, ir, Coeff, n, rnorm)
      Mdeout = 1

   end subroutine defcmn
!*****************************************************************************************

!*****************************************************************************************
!>
!  These subroutines solve the least squares problem Ax = b for
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
!  following call statements
!
!  CALL DBNDAC(...)  Introduce new blocks of data.
!
!  CALL DBNDSL(1,...)Compute solution vector and length of
!                    residual vector.
!
!  CALL DBNDSL(2,...)Given any row vector H solve YR = H for the
!                    row vector Y.
!
!  CALL DBNDSL(3,...)Given any column vector W solve RZ = W for
!                    the column vector Z.
!
!  The dots in the above call statements indicate additional
!  arguments that will be specified in the following paragraphs.
!
!  The user must dimension the array appearing in the call list..
!  G(MDG,NB+1)
!
!  Description of calling sequence for DBNDAC..
!
!  The entire set of parameters for DBNDAC are
!
!  Input.. All Type REAL variables are real(wp)
!
!  G(*,*)            The working array into which the user will
!                    place the MT by NB+1 block (C F) in rows IR
!                    through IR+MT-1, columns 1 through NB+1.
!                    See descriptions of IR and MT below.
!
!  MDG               The number of rows in the working array
!                    G(*,*).  The value of MDG should be >= MU.
!                    The value of MU is defined in the abstract
!                    of these subprograms.
!
!  NB                The bandwidth of the data matrix A.
!
!  IP                Set by the user to the value 1 before the
!                    first call to DBNDAC.  Its subsequent value
!                    is controlled by DBNDAC to set up for the
!                    next call to DBNDAC.
!
!  IR                Index of the row of G(*,*) where the user is
!                    to place the new block of data (C F).  Set by
!                    the user to the value 1 before the first call
!                    to DBNDAC.  Its subsequent value is controlled
!                    by DBNDAC. A value of IR .GT. MDG is considered
!                    an error.
!
!  MT,JT             Set by the user to indicate respectively the
!                    number of new rows of data in the block and
!                    the index of the first nonzero column in that
!                    set of rows (E F) = (0 C 0 F) being processed.
!
!  Output.. All Type REAL variables are real(wp)
!
!  G(*,*)            The working array which will contain the
!                    processed rows of that part of the data
!                    matrix which has been passed to DBNDAC.
!
!  IP,IR             The values of these arguments are advanced by
!                    DBNDAC to be ready for storing and processing
!                    a new block of data in G(*,*).
!
!  Description of calling sequence for DBNDSL..
!
!  The user must dimension the arrays appearing in the call list..
!
!  G(MDG,NB+1), X(N)
!
!  The entire set of parameters for DBNDSL are
!
!  Input.. All Type REAL variables are real(wp)
!
!  MODE              Set by the user to one of the values 1, 2, or
!                    3.  These values respectively indicate that
!                    the solution of AX = B, YR = H or RZ = W is
!                    required.
!
!  G(*,*),MDG,       These arguments all have the same meaning and
!   NB,IP,IR         contents as following the last call to DBNDAC.
!
!  X(*)              With mode=2 or 3 this array contains,
!                    respectively, the right-side vectors H or W of
!                    the systems YR = H or RZ = W.
!
!  N                 The number of variables in the solution
!                    vector.  If any of the N diagonal terms are
!                    zero the subroutine DBNDSL prints an
!                    appropriate message.  This condition is
!                    considered an error.
!
!  Output.. All Type REAL variables are real(wp)
!
!  X(*)              This array contains the solution vectors X,
!                    Y or Z of the systems AX = B, YR = H or
!                    RZ = W depending on the value of MODE=1,
!                    2 or 3.
!
!  RNORM             If MODE=1 RNORM is the Euclidean length of the
!                    residual vector AX-B.  When MODE=2 or 3 RNORM
!                    is set to zero.
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

      real(wp) :: g(Mdg, *), rho
      integer :: i, ie, ig, ig1, ig2, iopt, Ip, Ir, j, jg, Jt, &
                 k, kh, l, lp1, Mdg, mh, Mt, mu, Nb
      integer :: nbp1, nerr

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
         write (*, *) 'MDG<IR, PROBABLE ERROR.'
      end if

   end subroutine dbndac
!*****************************************************************************************

!*****************************************************************************************
!>
!  These subroutines solve the least squares problem Ax = b for
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
!  following call statements
!
!  CALL DBNDAC(...)  Introduce new blocks of data.
!
!  CALL DBNDSL(1,...)Compute solution vector and length of
!                    residual vector.
!
!  CALL DBNDSL(2,...)Given any row vector H solve YR = H for the
!                    row vector Y.
!
!  CALL DBNDSL(3,...)Given any column vector W solve RZ = W for
!                    the column vector Z.
!
!  The dots in the above call statements indicate additional
!  arguments that will be specified in the following paragraphs.
!
!  The user must dimension the array appearing in the call list..
!  G(MDG,NB+1)
!
!  Description of calling sequence for DBNDAC..
!
!  The entire set of parameters for DBNDAC are
!
!  Input.. All Type REAL variables are real(wp)
!
!  G(*,*)            The working array into which the user will
!                    place the MT by NB+1 block (C F) in rows IR
!                    through IR+MT-1, columns 1 through NB+1.
!                    See descriptions of IR and MT below.
!
!  MDG               The number of rows in the working array
!                    G(*,*).  The value of MDG should be >= MU.
!                    The value of MU is defined in the abstract
!                    of these subprograms.
!
!  NB                The bandwidth of the data matrix A.
!
!  IP                Set by the user to the value 1 before the
!                    first call to DBNDAC.  Its subsequent value
!                    is controlled by DBNDAC to set up for the
!                    next call to DBNDAC.
!
!  IR                Index of the row of G(*,*) where the user is
!                    the user to the value 1 before the first call
!                    to DBNDAC.  Its subsequent value is controlled
!                    by DBNDAC. A value of IR .GT. MDG is considered
!                    an error.
!
!  MT,JT             Set by the user to indicate respectively the
!                    number of new rows of data in the block and
!                    the index of the first nonzero column in that
!                    set of rows (E F) = (0 C 0 F) being processed.
!  Output.. All Type REAL variables are real(wp)
!
!  G(*,*)            The working array which will contain the
!                    processed rows of that part of the data
!                    matrix which has been passed to DBNDAC.
!
!  IP,IR             The values of these arguments are advanced by
!                    DBNDAC to be ready for storing and processing
!                    a new block of data in G(*,*).
!
!  Description of calling sequence for DBNDSL..
!
!  The user must dimension the arrays appearing in the call list..
!
!  G(MDG,NB+1), X(N)
!
!  The entire set of parameters for DBNDSL are
!
!  Input..
!
!  MODE              Set by the user to one of the values 1, 2, or
!                    3.  These values respectively indicate that
!                    the solution of AX = B, YR = H or RZ = W is
!                    required.
!
!  G(*,*),MDG,       These arguments all have the same meaning and
!   NB,IP,IR         contents as following the last call to DBNDAC.
!
!  X(*)              With mode=2 or 3 this array contains,
!                    respectively, the right-side vectors H or W of
!                    the systems YR = H or RZ = W.
!
!  N                 The number of variables in the solution
!                    vector.  If any of the N diagonal terms are
!                    zero the subroutine DBNDSL prints an
!                    appropriate message.  This condition is
!                    considered an error.
!
!  Output..
!
!  X(*)              This array contains the solution vectors X,
!                    Y or Z of the systems AX = B, YR = H or
!                    RZ = W depending on the value of MODE=1,
!                    2 or 3.
!
!  RNORM             If MODE=1 RNORM is the Euclidean length of the
!                    residual vector AX-B.  When MODE=2 or 3 RNORM
!                    is set to zero.
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
!  * 890831  Modified array declarations.  (WRB)
!  * 891006  Cosmetic changes to prologue.  (WRB)
!  * 891006  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dbndsl(Mode, g, Mdg, Nb, Ip, Ir, x, n, Rnorm)

      real(wp) :: g(Mdg, *), Rnorm, rsq, s, x(*)
      integer :: i, i1, i2, ie, ii, iopt, Ip, Ir, irm1, ix, j, &
                 jg, l, Mdg, Mode, n, Nb, nerr, np1

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
      write (*, *) 'A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR MATRIX.'

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
!  Householder transformation. Q = I + U*(U**T)/B
!
!     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
!     LPIVOT is the index of the pivot element.
!     L1,M   If L1 <= M   the transformation will be constructed to
!            zero elements indexed from L1 through M.   If L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
!                   IUE is the storage increment between elements.
!                                       On exit from H1 U() and UP
!                   contain quantities defining the vector U of the
!                   Householder transformation.   On entry to H2 U()
!                   and UP should contain quantities previously computed
!                   by H1.  These will not be modified by H2.
!     C()    On entry to H1 or H2 C() contains a matrix which will be
!            regarded as a set of vectors to which the Householder
!            transformation is to be applied.  On exit C() contains the
!            set of transformed vectors.
!     ICE    Storage increment between elements of vectors in C().
!     ICV    Storage increment between vectors in C().
!     NCV    Number of vectors in C() to be transformed. If NCV <= 0
!            no operations will be done on C().
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

      integer :: i, i2, i3, i4, Ice, Icv, incr, Iue, j, kl1, &
                 kl2, klp, l1, l1m1, Lpivot, m, mml1p2, Mode, Ncv
      real(wp) :: b, c(*), cl, clinv, ul1m1, sm, u(Iue, *), Up

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
   end module bspline_defc_module
!*****************************************************************************************