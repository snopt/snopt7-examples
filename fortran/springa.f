!     ------------------------------------------------------------------
!     File springA.f for snOptA
!
!     This is a main program to generate an optimal control problem
!     of arbitrary size and solve it by calling snOptA.
!
!     The problem size depends on a parameter T.  There are
!     2T constraints and 3T + 2 variables, as well as bounds
!     on the variables.  The first T constraints are quadratic in
!     T + 1 variables, and the objective function is quadratic in
!     T + 1 other variables.
!
!     The control problem models a spring, mass and damper system.
!     It is of the form
!
!   --------------------------------------------------------------------
!   | minimize    1/2 sum x(t)**2   (t = 0 to T)                       |
!   |                                                                  |
!   | subject to                                                       |
!   |     y(t+1)  =  y(t)  -  0.01 y(t)**2  -  0.004 x(t)  +  0.2 u(t) |
!   |                                                                  |
!   |     x(t+1)  =  x(t)  +  0.2  y(t),                               |
!   |                                                                  |
!   |     y(t)   >=  -1,     -0.2  <=  u(t)  <=  0.2,                  |
!   |                                                                  |
!   |                (all for t = 0 to T-1)                            |
!   | and                                                              |
!   |     y(0)    =   0,      y(T)  =  0,       x(0) = 10.             |
!   --------------------------------------------------------------------
!
!     For large enough T (e.g. T >= 90), the optimal objective value
!     is approximately 1186.382.
!
!     This model with T = 100 was used as test problem 5.11 in
!     B. A. Murtagh and M. A. Saunders (1982), A projected Lagrangian
!     algorithm and its implementation for sparse nonlinear constraints,
!     Mathematical Programming Study 16, 84--117.
!
!     14 Nov 1994: First version of spring.f, derived from manne.f.
!     25 Apr 2001: Updated for SNOPT 6.1.
!     ------------------------------------------------------------------
      program
     &     springa
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG, lencw, leniw, lenrw
      parameter
     &   ( maxF   = 1000,  maxn   = 1500,
     &     lenA   = 10000, lenG   = 5000,
     &     nxname = 1,     nFname =    1 )
      integer
     &     iAfun(lenA), jAvar(lenA),
     &     iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character*8
     &     Prob, xnames(nxname), Fnames(nFname)
      double precision
     &     objAdd, sInf,
     &     A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)

!     SNOPT workspace

      parameter
     &     (lenrw = 50000,
     &      leniw = 50000,
     &      lencw =   500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character*8
     &     cw(lencw)
      integer
     &     DerOpt, Errors, neA, neG, ObjRow, INFO, iPrint, iPrt, iSpecs,
     &     iSumm, iSum, itnlim, mincw, miniw, minrw,
     &     nF, n, nInf, nOut, nS, T
      logical
     &     byname
      character*20
     &     lfile
      external
     &     userf, userfg
!     ------------------------------------------------------------------
      integer            Cold
      parameter         (Cold   = 0)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by sntest.

      iSpecs =  4  ! equivalenced to springa.spc
      iPrint =  9  ! equivalenced to springa.out
      iSumm  =  6  ! summary file goes to standard output...
      nOut   =  6  ! ... as do messages from this program.

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'springa.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'springa.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     1. Solve springa with only F specified.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '1. springa with  only F specified.'
      write(nOut, *) '-------------------------------------------'

!     Generate a  T-period problem.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.

      Errors = 0

      T      = 100
      iPrt   =   0
      iSum   =   0
      call snSeti
     &   ( 'Problem number', T, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call springData0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) then
         go to 910
      end if

!     snOptA will compute the Jacobian by finite-differences.
!     The user has the option of calling  snJac  to define the
!     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).

      call snJac
     &   ( INFO, nF, n, userf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) go to 900

!     ------------------------------------------------------------------
!     Warn snOptA that userf does not set the derivatives.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, userf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'springa (f) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .gt. 30) go to 900

!     ------------------------------------------------------------------
!     2. Solve springa with F and G specified.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '2. with F and G specified.'
      write(nOut, *) '-------------------------------------------'

!     Read a Specs file.  This must include "Problem number  T"
!     for some integer T.

      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

!     Generate an T-period problem.

      call springData1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 910

!     Specify options that were not set in the Specs file.

      DerOpt = 1
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      itnlim = 20000
      call snSeti
     &   ( 'Iterations       ', itnlim, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, userfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'springa (fg) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900

!     ------------------------------------------------------------------
!     3. Solve springa with linear constraint scaling.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '3. springa with linear constraint scaling.'
      write(nOut, *) '-------------------------------------------'

      call springData1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 910

      call snSeti
     &   ( 'Scale option     ', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, userfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'springa (LC scaling) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900

      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
  910 stop

 4000 format(/  a, 2x, a  )

      end ! program springa

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine springData0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, nF, n, ObjRow, lencw, leniw, lenrw,
     &     xstate(maxn), iw(leniw)
      double precision
     &     objAdd,
     &     xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     springData0  generates data for the test problem springa
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     On exit,
!     Errors    is 0 if there is enough storage, 1 otherwise.
!     nF        is the number of problem functions
!               (objective and constraints, linear and nonlinear).
!     n         is the number of variables.
!     xlow      holds the lower bounds on x.
!     xupp      holds the upper bounds on x.
!     Flow      holds the lower bounds on F.
!     Fupp      holds the upper bounds on F.

!     xstate(1:n)  is a set of initial states for each x  (0,1,2,3,4,5).
!     Fstate(1:nF) is a set of initial states for each F  (0,1,2,3,4,5).
!     x (1:n)      is a set of initial values for x.
!     Fmul(1:nF)   is a set of initial values for the dual variables.
!
!     19 Jul 2000: First version of springA based on SNOPT 5.3 spring.
!     03 Jun 2001: Current version.
!     ==================================================================
      integer
     &     i, nOut, jt, ju, ju0, jx0, jx, jy0, jy, nCon, Obj, T
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
!     ------------------------------------------------------------------
      nOut = 6
!     ------------------------------------------------------------------
!     The following call fetches T, the number of nonlinear constraints.
!     ------------------------------------------------------------------
      Errors = 0
      call snGeti
     &   ( 'Problem number', T,   Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     Check if there is enough storage.

      if ( T    .le. 1         .or.
     &     maxF .lt. 2*T       .or.
     &     maxn .lt. 3*T + 2       ) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  T =', T
         Errors = 1
      end if

      if (Errors .ge. 1) go to 910

!     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  100) then
         prob(1:6) = 'Spring'
      else if (T .lt. 1000) then
         prob(1:5) = 'Sprin'
      else
         prob(1:3) = 'Spr'
      end if

      write(nOut, *) 'Problem springa.    T =', T

      n      = 3*T  + 2
      nCon   = 2*T
      nF     = nCon + 1

      Obj    = nF
      ObjRow = Obj

!     objAdd = 0.0 means there is no constant to be added to the
!            objective.

      objAdd = zero

!     The variables are ordered as follows:
!     state   variables y:     1: T+1
!     state   variables x:   T+2:2T+2
!     control variables u:  2T+3:3T+2
!     jy, jx, ju are the base indices for the corresponding variables

      jy0  = 1                  ! points to state   y(0)  in x.
      jx0  = jy0 + T + 1        ! points to state   x(0)  in x.
      ju0  = jx0 + T + 1        ! points to control u(0)  in x

!     ------------------------------------------------------------------
!     Initialize the bounds
!     ------------------------------------------------------------------
      do jt = 0, T
         jy = jy0 + jt
         jx = jx0 + jt

!        Initialize the bounds and values of the states x.

         xlow(jx)   =  bminus
         xupp(jx)   =  bplus
         x(jx)      =  zero
         xstate(jx) =  3

!        Initialize the bounds and values of the states y.

         xlow(jy)   = -one
         xupp(jy)   =  bplus
         x(jy)      = -one
         xstate(jy) =  0
      end do

      do jt = 0, T-1
         ju = ju0 + jt

!        Initialize the bounds and values of the constrols u.

         xlow(ju)   = -0.2d0
         xupp(ju)   =  0.2d0
         x(ju)      =  zero
         xstate(ju) =  3
      end do

!     Fix the boundary conditions.

      xlow(jy0)   = zero
      xupp(jy0)   = zero

      xlow(jy0+T) = zero
      xupp(jy0+T) = zero
      x(jy0+T)    = zero

      xlow(jx0)   = 10.0d0
      xupp(jx0)   = 10.0d0
      x(jx0)      = 10.0d0

!     Bounds on F

      do i = 1, nCon
         Flow(i) = zero
         Fupp(i) = zero
      end do

      do i = 1, nCon
         Fmul(i) = zero
      end do

!     Set the objective and its bounds.

      Fmul(ObjRow) =  zero
      Flow(ObjRow) =  bminus
      Fupp(ObjRow) =  bplus

  910 return

      end ! subroutine springData0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userf
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, n, nF, lenG,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     This is userf for problem springa.
!
!     ==================================================================
      double precision
     &     FObj, ut, xt, xtp1, yt, ytp1
      integer
     &     jt, ju, ju0, jx0, jx, jy0, jy, lin, nln, ObjRow, T
!     ------------------------------------------------------------------
      T      = (n - 2)/3
      ObjRow = 2*T + 1

      lin    = 1                ! points to linear    constraints in F.
      nln    = lin + T          ! points to nonlinear constraints in F.

      jy0    = 1                ! points to state   y(0)  in x.
      jx0    = jy0 + T + 1      ! points to state   x(0)  in x.
      ju0    = jx0 + T + 1      ! points to control u(0)  in x

      FObj = 0.0d0

      do jt = 0, T-1
         jy = jy0 + jt
         jx = jx0 + jt
         ju = ju0 + jt

         yt   = x(jy)
         ytp1 = x(jy+1)

         xt   = x(jx)
         xtp1 = x(jx+1)

         ut   = x(ju)

         F(nln) = 1.0d-2*yt**2 - yt + ytp1 + 4.0d-3*xt       - 0.2d0*ut
         F(lin) = -0.2d0*yt                -        xt + xtp1
         FObj   = FObj +  xt**2

         nln    = nln + 1
         lin    = lin + 1
      end do

!     Set the objective row.

      F(ObjRow) = (FObj + xtp1**2)/2.0d+0

      end ! subroutine userf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine springData1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     objAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     springData1  generates data for the test problem springa
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     The  nF by n  Jacobian is written as the sum of
!     two  nF by n  sparse matrices G and A, i.e.,  J = A + G,  where
!     A  and  G  contain contain the constant and nonlinear
!     elements of J respectively.
!
!     The nonzero pattern of G and A is specified by listing the
!     coordinates (i.e., row and column indices) of each nonzero.
!     Note that the coordinates specify the overall STRUCTURE of the
!     sparsity, not the Jacobian entries that happen to be zero at
!     the initial point.
!
!     The coordinates of the kth nonzero of  G  are defined
!     by iGfun(k) and jGvar(k)  (i.e., if i=iGfun(k) and j=jGvar(k),
!     then G(k) is the ijth element of G.)  Any known values of G(k)
!     must be assigned by the user in the routine userfg.
!
!     The coordinates of the kth nonzero of  A  are defined by
!     iAfun(k) and jAvar(k)  (i.e., if i=iAfun(k) and j=jAvar(k),
!     then A(k) is the ijth element of A.)  All values of A must be
!     assigned before the call to SNOPT.
!
!     The elements of A and G can be stored in any order, (e.g., by rows
!     or columns or mixed).
!
!     RESTRICTIONS:
!     1.  A nonzero entry of J must be specified as either an element
!         of A, or an element of G, but NOT BOTH (i.e.,  coordinates of
!         A  and  G  must not overlap.  Elements that are a sum of a
!         constant and varying part must be included in G and loaded
!         by userfg.
!
!     2.  If the computed value of an element of G happens to be zero
!         at a given point, it must still be loaded in userfg. (The
!         order of the coordinates is meaningful in SNOPT.)
!
!     On exit,
!     Errors    is 0 if there is enough storage, 1 otherwise.
!     nF        is the number of problem functions
!               (objective and constraints, linear and nonlinear).
!     n         is the number of variables.
!     neG       is the number of nonzeros in Jn.
!     neA       is the number of nonzeros in Jc.
!     xlow      holds the lower bounds on x.
!     xupp      holds the upper bounds on x.
!     Flow      holds the lower bounds on F.
!     Fupp      holds the upper bounds on F.

!     xstate(1:n)  are the initial states for each x  (0,1,2,3,4,5).
!     Fstate(1:nF) are the initial states for each F  (0,1,2,3,4,5).
!     x (1:n)      are the initial values for x.
!     Fmul(1:nF)   are the initial values for the dual variables.
!
!     19 Jul 2000: First version of springA based on SNOPT 5.3 spring.
!     03 Jun 2001: Current version.
!     ==================================================================
      integer
     &     i, nOut, jt, ju, ju0, jx0, jx, jy0, jy, lin, nCon,
     &     nln, Obj, T
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
!     ------------------------------------------------------------------
      nOut = 6
!     ------------------------------------------------------------------
!     The following call fetches T, the number of nonlinear constraints.
!     It is specified at runtime in the SPECS file.
!     ------------------------------------------------------------------
      Errors = 0
      call snGeti
     &   ( 'Problem number', T,
     &     Errors, cw, lencw, iw, leniw, rw, lenrw )

!     Check if there is enough storage.

      if ( T    .le. 1         .or.
     &     maxF .lt. 2*T       .or.
     &     maxn .lt. 3*T + 2   .or.
     &     lenA .lt. 6*T       .or.
     &     lenG .lt. 2*T + 1       ) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  T =', T
         Errors = 1
      end if

      if (Errors .ge. 1) go to 910

!     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  100) then
         prob(1:6) = 'Spring'
      else if (T .lt. 1000) then
         prob(1:5) = 'Sprin'
      else
         prob(1:3) = 'Spr'
      end if

      write(nOut, *) 'Problem springA.    T =', T

      n      = 3*T  + 2
      nCon   = 2*T
      nF     = nCon + 1

      Obj    = nF
      ObjRow = Obj

!     objAdd = 0.0 means there is no constant to be added to the
!            objective.

      objAdd = zero

!     The variables are ordered as follows:
!     state   variables y:     1: T+1
!     state   variables x:   T+2:2T+2
!     control variables u:  2T+3:3T+2
!     jx, jy, ju are the base indices for the corresponding variables

      neA  = 0
      neG  = 0

      lin  = 1                  ! points to linear    constraints in F.
      nln  = lin + T            ! points to nonlinear constraints in F.

      jy0  = 1                  ! points to state   y(0)  in x.
      jx0  = jy0 + T + 1        ! points to state   x(0)  in x.
      ju0  = jx0 + T + 1        ! points to control u(0)  in x

      do jt = 0, T-1
         jx = jx0 + jt
         jy = jy0 + jt
         ju = ju0 + jt

!        x states.

         neA        =  neA + 1
         iAfun(neA) =  nln
         jAvar(neA) =  jx
         A(neA)     =  0.004D0

         neA        =  neA + 1
         iAfun(neA) =  lin
         jAvar(neA) =  jx
         A(neA)     = -one

         neA        =  neA + 1
         iAfun(neA) =  lin
         jAvar(neA) =  jx  + 1
         A(neA)     =  one

!        y states

         neG        =  neG + 1
         iGfun(neG) =  nln
         jGvar(neG) =  jy
!        G(neG)     =  0.02d0*yt - 1

         neA        =  neA + 1
         iAfun(neA) =  nln
         jAvar(neA) =  jy  + 1
         A(neA)     =  one

         neA        =  neA + 1
         iAfun(neA) =  lin
         jAvar(neA) =  jy
         A(neA)     = -0.2d0

!        u controls.

         neA        = neA + 1
         iAfun(neA) = nln
         jAvar(neA) = ju
         A(neA)     = -0.2d0

!        Objective

         neG        =  neG + 1
         iGfun(neG) =  ObjRow
         jGvar(neG) =  jx
!        G(neG)     =  xt

         nln        = nln + 1
         lin        = lin + 1
      end do

!     One last gradient element.

      neG        =  neG + 1
      iGfun(neG) =  ObjRow
      jGvar(neG) =  jx0 + T
!     G(neG)     =  xt

!     ------------------------------------------------------------------
!     Initialize the bounds
!     ------------------------------------------------------------------
      do jt = 0, T
         jy = jy0 + jt
         jx = jx0 + jt

!        Initialize the bounds and values of the states y.

         xlow(jy)   = -one
         xupp(jy)   =  bplus
         x(jy)      = -one
         xstate(jy) =  0

!        Initialize the bounds and values of the states x.

         xlow(jx)   =  bminus
         xupp(jx)   =  bplus
         x(jx)      =  zero
         xstate(jx) =  3
      end do

      do jt = 0, T-1
         ju = ju0 + jt

!        Initialize the bounds and values of the constrols u.

         xlow(ju)   = -0.2d0
         xupp(ju)   =  0.2d0
         x(ju)      =  zero
         xstate(ju) =  3
      end do

!     Fix the boundary conditions.

      xlow(jy0)   = zero
      xupp(jy0)   = zero

      xlow(jy0+T) = zero
      xupp(jy0+T) = zero
      x(jy0+T)    = zero

      xlow(jx0)   = 10.0d0
      xupp(jx0)   = 10.0d0
      x(jx0)      = 10.0d0

!     Bounds on F

      do i = 1, nCon
         Flow(i) = zero
         Fupp(i) = zero
      end do

      do i = 1, nCon
         Fmul(i) = zero
      end do

!     Set the objective and its bounds.

      Fmul(ObjRow) =  zero
      Flow(ObjRow) =  bminus
      Fupp(ObjRow) =  bplus

  910 return

      end ! subroutine springData1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, n, nF, lenG,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     This is userfg for problem springa.
!
!     ==================================================================
      double precision
     &     FObj, xt, yt
      integer
     &     jt, jx0, jx, jy0, jy, lin, neG, nln, ObjRow, T
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero  = 0.0d+0,   one    = 1.0d+0)
!     ------------------------------------------------------------------

      T      = (n - 2)/3
      ObjRow = 2*T + 1

!     The variables are ordered as follows:
!     state   variables y:     1: T+1
!     state   variables x:   T+2:2T+2
!     control variables u:  2T+3:3T+2
!     jy, jx, ju are the base indices for the corresponding variables

      lin  = 1                  ! points to linear    constraints in F.
      nln  = lin + T            ! points to nonlinear constraints in F.
      neG  = 0

      jy0  = 1                  ! points to state   y(0)  in x.
      jx0  = jy0 + T + 1        ! points to state   x(0)  in x.

      Fobj = zero

      do jt = 0, T-1
         jy = jy0 + jt
         jx = jx0 + jt

!        y states

         yt   = x(jy)

!        x states.

         xt   = x(jx)

         if (needF .gt. 0) then
            F(nln) =  0.01d0*yt**2 - yt
            FObj   = FObj +  xt**2
         end if

         if (needG .gt. 0) then
            neG        =  neG + 1
!           iGfun(neG) =  nln
!           jGvar(neG) =  jy
            G(neG)     =  0.02d0*yt - one

!           Objective

            neG        =  neG + 1
!           iGfun(neG) =  ObjRow
!           jGvar(neG) =  jx
            G(neG)     =  xt
         end if

         nln        = nln + 1
      end do

!     One last gradient element.

      if (needF .gt. 0) then
          F(ObjRow)    = (FObj +  x(jx0+T)**2)/2.0d+0
      end if

      if (needG .gt. 0) then
         neG        =  neG + 1
!        iGfun(neG) =  ObjRow
!        jGvar(neG) =  jx0 + T
         G(neG)     =  x(jx0+T)
      end if

      end ! of userfg
