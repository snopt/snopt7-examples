!     ------------------------------------------------------------------
!     File t4manneA.f
!     This is a main program to test subroutine snOptA, which is
!     part of the SNOPT 6 package.
!     It generates the problem called MANNE and asks snOptA to solve it.
!
!     19 Jul 2000: First version of t4manneA.
!     27 Oct 2002: Current version.
!     ------------------------------------------------------------------
      program
     &     t4main

      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG,
     &     lencw, leniw, lenrw
      parameter
     &   ( maxF   = 1000,
     &     maxn   = 1500,
     &     lenA   = 5000, lenG   = 5000,
     &     nxname =    1, nFname =    1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character
     &     Prob*8, xnames(nxname)*8, Fnames(nFname)*8
      double precision
     &     ObjAdd, sInf, A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)

!     SNOPT workspace

      parameter
     &     (lenrw = 100000,
     &      leniw =  50000,
     &      lencw =    500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character
     &     cw(lencw)*8
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
      integer            Cold,       Basis,      Warm
      parameter         (Cold   = 0, Basis  = 1, Warm  = 2)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by sntest.

      iSpecs =  4  ! equivalenced to manne.spc
      iPrint =  9  ! equivalenced to manne.out
      iSumm  =  6  ! summary file goes to standard output...
      nOut   =  6  ! ... as do messages from this program.

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't4mannea.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't4mannea.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Generate a  T-period problem.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors = 0

      T      = 30
      iPrt   =  0
      iSum   =  0
      call snSeti( 'Problem number', T, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call t4dat0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ge. 1) go to 910

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
     &     ObjAdd, ObjRow, Prob, userf,
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
      write(nOut, *) 't4mannea (FD) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900

!     ------------------------------------------------------------------
!     Solve the problem again, using a Cold start (Start = 0).
!
!     Read a Specs file (optional).
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

!     ------------------------------------------------------------------
!     Generate an T-period problem.
!     ------------------------------------------------------------------
      call t4dat1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ge. 1) then
         go to 910
      end if

!     ------------------------------------------------------------------
!     Specify options that were not set in the Specs file.
!
!     As the problem has not changed,  we set the call status to
!     indicate that there is nothing special about the first call to
!     the user-defined functions.
!     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      itnlim = 20000
      call snSeti( 'Iterations',        itnlim, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Call status       ',             0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Solve the problem again.
!     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, userfg,
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
      write(nOut, *) 't4mannea (g, J) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 900

!     ------------------------------------------------------------------
!     Solve the problem again, using a Warm start at the solution.
!     ------------------------------------------------------------------
      call snSeti
     &   ( 'Call status       ',             0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Warm, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, userfg,
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
      write(nOut, *) 't4mannea (warm start) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
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

      end ! program t4mainA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t4dat0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, nF, n, ObjRow, lencw, leniw, lenrw,
     &     xstate(maxn), iw(leniw)
      double precision
     &     ObjAdd,
     &     xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     t4dat0  generates data for the test problem t4manne
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
!     19 Jul 2000: First version of t4manne based on SNOPT 5.3 t4manne.
!     27 Oct 2002: Current version.
!     ==================================================================
      integer
     &     i, Out, iPrt, iSum, jC, jI, jM, jY, k, Obj, T
      double precision
     &     scale
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
!     ------------------------------------------------------------------
      Out = 6
!     ------------------------------------------------------------------
!     The following call fetches T, the number of nonlinear constraints.
!     It is specified at runtime in the SPECS file.
!     ------------------------------------------------------------------
      call snGeti
     &   ( 'Problem number', T, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (T .le. 1  .or.  T .gt. maxF/2) then
         write(Out, *) 'Invalid  T  specified:', T
         go to 910
      end if

!     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  1000) then
         prob(1:5) = 'Manne'
      else if (T .lt. 10000) then
         prob(1:4) = 'Mann'
      else
         prob(1:3) = 'Man'
      end if

      write(Out, *) 'Problem MANNE.    T =', T

      nF     = T*2 + 1
      n      = T*3

!     Check if there is enough storage.

      Errors = 0
      if (nF    .gt. maxF) Errors = 1
      if (n     .gt. maxn) Errors = 1
      if (Errors  .ge.    1) then
         write(Out, *) 'Not enough storage to generate a problem ',
     &                  'with  T =', T
         go to 910
      end if

      iPrt   = 0
      iSum   = 0
      call snSet
     &   ( 'Maximize', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      ObjAdd = zero             ! no additive obj. constant
      Obj    = nF
      ObjRow = Obj

!     The variables are ordered as follows:
!     Kt Variables    1: T    represent Kapital  (nonlinear)
!     Ct            T+1:2T    represent Consumption
!     It           2T+1:3T    represent Investment
!
!     The columns for Kapital are nonlinear in the first T rows.
!     The columns for Consumption form -I in the first T rows.
!     The columns for Investment  form -I in the first T rows and -I
!        in the last T rows.
!
!     The first T constraints (money)    are nonlinear, and
!     the next  T constraints (capacity) are linear.
!     The nonlinear constraint part of the Jacobian is a T x T diagonal.
!     We generate the sparsity pattern here.
!     Real values for the derivatives are computed by userfg.

!     jC and jI are base indices for the Ct and It variables.

      jC    =   T
      jI    = 2*T

!     The first T components are nonlinear, and the next T are linear.

!     Set lower and upper bounds for Kt, Ct, It.
!     Also initial values and initial states for all variables.
!     The nonlinear variables are the most important.
!     We make them all superbasic.
!     The rest are ok nonbasic.
!     For test purposes, we want the initial x to be infeasible
!     with respect to the linear constraints.
!     Try setting the last Kapital too high.

      do k = 1, T
         xlow(   k) = 3.05d+0
         xupp(   k) = bplus
         xlow(jC+k) = 0.95d+0
         xupp(jC+k) = bplus
         xlow(jI+k) = 0.05d+0
         xupp(jI+k) = bplus

         x(   k)  = 3.0d+0 + (k - 1)/10.0d+0
         x(jC+k)  = xlow(jC+k)
         x(jI+k)  = xlow(jI+k)

!-->     xstate(   k) = 2
         xstate(   k) = 0
         xstate(jC+k) = 0
         xstate(jI+k) = 0

         if (k .eq. T) then
            x(k)      = 1.0d+3
            xstate(k) = 2
         end if
      end do

!     The first Capital is fixed.
!     The last three Investments are bounded.
!     Fudge them to be the normal ones for T = 10.

      scale        = T / 10.0d+0
      xupp(1)      = xlow(1)
      x(1)         = xlow(1)
      xstate(1)    = 0
      xupp(jI+T-2) = 0.112d+0 * scale
      xupp(jI+T-1) = 0.114d+0 * scale
      xupp(jI+T  ) = 0.116d+0 * scale

!     Set bounds on F.
!     The T nonlinear (Money)    components are >=.
!     The T    linear (Capacity) components are <=.
!     We no longer need to set initial values and states for slacks
!     (assuming SNOPT does a cold start).

      jM     = 0
      jY     = T

      do k = 1, T
         Flow(jM+k) = zero
         Fupp(jM+k) = bplus
         Flow(jY+k) = bminus
         Fupp(jY+k) = zero

!-       x (jM+k) = zero
!-       x (jY+k) = zero
!-       Fstate(jM+k) = 0
!-       Fstate(jY+k) = 0
      end do

!     The last Money and Capacity components have a Range.

      Fupp(jM+T) =   10.0d+0
      Flow(jY+T) = - 20.0d+0

!     Set the objective and its bounds.

      Fmul(Obj) = zero
      Flow(Obj) = bminus
      Fupp(Obj) = bplus

!     Initialize pi.

      do i = 1, T
         Fmul(i)   = - one
         Fmul(T+i) = + one
      end do

  910 return

      end ! subroutine t4dat0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userf
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     This is userf for problem t4manneA.
!
!     The data bt(*) is computed by userf  on first entry.
!
!     ==================================================================
      integer
     &     i, Out, j, jC, jI, k, Obj, T
      double precision
     &     a, beta, Fi, FObj, gfac, grow, xc0, xCon, xi0, xk0
!     ------------------------------------------------------------------
      double precision   growth
      parameter         (growth = .03d+0)
      double precision   zero,              one
      parameter         (zero   = 0.0d+0,   one    = 1.0d+0)
      double precision   b, at(365), bt(365)
      common    /manne / b, at     , bt
!     ------------------------------------------------------------------
      T      = n/3
      Obj    = 2*T + 1
      Out    = 6

!     ---------------------------------------
!     First entry.  Define b, at(*) and bt(*)
!     for this and all subsequent entries.
!     ---------------------------------------
      if (Status .eq. 1) then
         grow   = 0.03d+0
         beta   = 0.95d+0
         xk0    = 3.0d+0
         xc0    = 0.95d+0
         xi0    = 0.05d+0
         b      = 0.25d+0
         if (Out .gt. 0) write(Out, 1000) T, b

         a      = (xc0 + xi0) / xk0**b
         gfac   = (one + grow)**(one - b)
         at(1)  = a*gfac
         bt(1)  = beta

         do j  = 2, T
            at(j) = at(j-1)*gfac
            bt(j) = bt(j-1)*beta
         end do

         bt(T) = bt(T) / (one - beta)
      end if

!     -------------
!     Normal entry.
!     -------------
!     Constraints first.
!     Only the first T components of F are nonlinear.

!     jC and jI are base indices for the Ct and It variables.

      jC    =   T
      jI    = 2*T

      do    i = 1, T
         j    = i
         F(i) = at(i)*x(j)**b - x(jC+j) - x(jI+j)
      end do

      do k = 1, T
         i = T + k

         if (k .lt. T) then
            F(i) =      - x(k) + x(k+1) - x(jI+k)
         else
            F(i) = growth*x(k)          - x(jI+k)
         end if
      end do

!     Set the objective element.

      FObj = zero
      do j = 1, T
         xCon = x(T+j)
         FObj = FObj  +  bt(j) * log(xCon)
!Min     FObj = FObj  -  bt(j) * log(xCon)
      end do

      F(Obj) = FObj

!     ------------
!     Final entry.
!     ------------
      if (Status .ge. 2) then
         if (Out .gt. 0) write(Out, 2000) (F(j), j = 1, T)
      end if
      return

 1000 format(// ' This is problem  Manne.   T =', i4, '   b =', f8.3/)
 2000 format(// ' Final nonlinear function values' / (5f12.5))

      end ! subroutine userf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t4dat1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     t4dat1  generates data for the test problem t4manne
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
!     neH       is the number of nonzeros in H(neH).
!     xlow      holds the lower bounds on x.
!     xupp      holds the upper bounds on x.
!     Flow      holds the lower bounds on F.
!     Fupp      holds the upper bounds on F.

!     xstate(1:n)  is a set of initial states for each x  (0,1,2,3,4,5).
!     Fstate(1:nF) is a set of initial states for each F  (0,1,2,3,4,5).
!     x (1:n)      is a set of initial values for x.
!     Fmul(1:nF)   is a set of initial values for the dual variables.
!
!     19 Jul 2000: First version of t4manne based on SNOPT 5.3 t4manne.
!     03 Dec 2000: Current version.
!     ==================================================================
      integer
     &     i, Out, iPrt, iSum, j, jC, jI, jM, jY, k, Obj, T
      double precision
     &     scale
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
      double precision   growth
      parameter         (growth = .03d+0)
!     ------------------------------------------------------------------
      Out = 6
!     ------------------------------------------------------------------
!     The following call fetches T, the number of nonlinear constraints.
!     It is specified at runtime in the SPECS file.
!     ------------------------------------------------------------------
      call snGeti
     &   ( 'Problem number', T, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (T .le. 1  .or.  T .gt. maxF/2) then
         write(Out, *) 'Invalid  T  specified:', T
         go to 910
      end if

!     Write T into the problem name.

      write(prob, '(i8)') T
      if      (T .lt.  1000) then
         prob(1:5) = 'Manne'
      else if (T .lt. 10000) then
         prob(1:4) = 'Mann'
      else
         prob(1:3) = 'Man'
      end if

      write(Out, *) 'Problem MANNE.    T =', T

      nF     = T*2 + 1
      n      = T*3

!     Check if there is enough storage.

      Errors = 0
      if (nF    .gt. maxF ) Errors = 1
      if (n     .gt. maxn ) Errors = 1
      if (Errors  .ge. 1) then
         write(Out, *) 'Not enough storage to generate a problem ',
     &                 'with  T =', T
         go to 910
      end if

      iPrt   = 0
      iSum   = 0
      call snSet
     &   ( 'Maximize', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      ObjAdd = zero             ! no additive obj. constant
      Obj    = nF
      ObjRow = Obj

!     The variables are ordered as follows:
!     Kt Variables    1: T    represent Kapital  (nonlinear)
!     Ct            T+1:2T    represent Consumption
!     It           2T+1:3T    represent Investment
!
!     The columns for Kapital are nonlinear in the first T rows.
!     The columns for Consumption form -I in the first T rows.
!     The columns for Investment  form -I in the first T rows and -I
!        in the last T rows.
!
!     The first T constraints (money)    are nonlinear, and
!     the next  T constraints (capacity) are linear.
!     The nonlinear constraint part of the Jacobian is a T x T diagonal.
!     We generate the sparsity pattern here.
!     Real values for the derivatives are computed by userfg.

!     jC and jI are base indices for the Ct and It variables.

      jC    =   T
      jI    = 2*T

!     The first T components are nonlinear, and the next T are linear.

      neG    = 0
      neA    = 0

      do i = 1, T
         j          = i

!        F(i)       = at(i)*x(j)**b - x(jC+j) - x(jI+j)

         neG        = neG + 1
         iGfun(neG) = i
         jGvar(neG) = j
!        G(neG)     = b*F(i)/xKap

         neA        = neA + 1
         iAfun(neA) = i
         jAvar(neA) = jC + j
         A(neA)     = -one

         neA        = neA + 1
         iAfun(neA) = i
         jAvar(neA) = jI + j
         A(neA)     = -one
      end do

!     Linear components of F.

      do k = 1, T
         i = T + k

         if (k .lt. T) then

!           F(i) = - x(k) + x(k+1) - x(jI+k)

            neA        = neA + 1
            iAfun(neA) = i
            jAvar(neA) = k
            A(neA)     = -one

            neA        = neA + 1
            iAfun(neA) = i
            jAvar(neA) = k + 1
            A(neA)     =  one
         else

!           F(T+T) = growth*x(T)   - x(jI+T)

            neA        = neA + 1
            iAfun(neA) = i
            jAvar(neA) = k
            A(neA)     = growth
         end if

         neA        = neA + 1
         iAfun(neA) = i
         jAvar(neA) = jI + k
         A(neA)     = -one
      end do

!     Set the objective element.

      do j = 1, T

!Max     FObj = FObj  +  bt(j) * log(x(T+j))
!Min     FObj = FObj  -  bt(j) * log(x(T+j))

         neG        = neG + 1
         iGfun(neG) = Obj
         jGvar(neG) = T + j
!Max     G(neG)     = + bt(j) / xCon
!Min     G(neG)     = - bt(j) / xCon

      end do

!     Set lower and upper bounds for Kt, Ct, It.
!     Also initial values and initial states for all variables.
!     The nonlinear variables are the most important.
!     We make them all superbasic.
!     The rest are ok nonbasic.
!     For test purposes, we want the initial x to be infeasible
!     with respect to the linear constraints.
!     Try setting the last Kapital too high.

      do k = 1, T
         xlow(   k) = 3.05d+0
         xupp(   k) = bplus
         xlow(jC+k) = 0.95d+0
         xupp(jC+k) = bplus
         xlow(jI+k) = 0.05d+0
         xupp(jI+k) = bplus

         x(   k)  = 3.0d+0 + (k - 1)/10.0d+0
         x(jC+k)  = xlow(jC+k)
         x(jI+k)  = xlow(jI+k)

!-->     xstate(   k) = 2
         xstate(   k) = 0
         xstate(jC+k) = 0
         xstate(jI+k) = 0

         if (k .eq. T) then
            x(k)      = 1.0d+3
            xstate(k) = 2
         end if
      end do

!     The first Capital is fixed.
!     The last three Investments are bounded.
!     Fudge them to be the normal ones for T = 10.

      scale        = T / 10.0d+0
      xupp(1)      = xlow(1)
      x(1)         = xlow(1)
      xstate(1)    = 0
      xupp(jI+T-2) = 0.112d+0 * scale
      xupp(jI+T-1) = 0.114d+0 * scale
      xupp(jI+T  ) = 0.116d+0 * scale

!     Set bounds on F.
!     The T nonlinear (Money)    components are >=.
!     The T    linear (Capacity) components are <=.
!     We no longer need to set initial values and states for slacks
!     (assuming SNOPT does a cold start).

      jM     = 0
      jY     = T

      do k = 1, T
         Flow(jM+k) = zero
         Fupp(jM+k) = bplus
         Flow(jY+k) = bminus
         Fupp(jY+k) = zero

!-       x (jM+k) = zero
!-       x (jY+k) = zero
!-       Fstate(jM+k) = 0
!-       Fstate(jY+k) = 0
      end do

!     The last Money and Capacity components have a Range.

      Fupp(jM+T) =   10.0d+0
      Flow(jY+T) = - 20.0d+0

!     Set the objective and its bounds.

      Fmul(Obj) = zero
      Flow(Obj) = bminus
      Fupp(Obj) = bplus

!     Initialize pi.

      do i = 1, T
         Fmul(i)   = - one
         Fmul(T+i) = + one
      end do

  910 return

      end ! of t4dat1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     This is userfg for problem t4manne.
!
!     The data bt(*) is computed by userfg  on first entry.
!
!     ==================================================================
      integer
     &     i, Out, j, neG, Obj, T
      double precision
     &     a, beta, Fi, FObj, gfac, Gi, grow,
     &     xc0, xCon, xKap, xi0, xk0
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero  = 0.0d+0,   one    = 1.0d+0)
      double precision   b, at(365), bt(365)
      common    /manne / b, at     , bt
!     ------------------------------------------------------------------
      T      = n/3
      Obj    = 2*T + 1

      Out    = 6

!     ---------------------------------------
!     First entry.  Define b, at(*) and bt(*)
!     for this and all subsequent entries.
!     ---------------------------------------
      if (Status .eq. 1) then
         grow   = 0.03d+0
         beta   = 0.95d+0
         xk0    = 3.0d+0
         xc0    = 0.95d+0
         xi0    = 0.05d+0
         b      = 0.25d+0
         if (Out .gt. 0) write(Out, 1000) T, b

         a      = (xc0 + xi0) / xk0**b
         gfac   = (one + grow)**(one - b)
         at(1)  = a*gfac
         bt(1)  = beta

         do j  = 2, T
            at(j) = at(j-1)*gfac
            bt(j) = bt(j-1)*beta
         end do

         bt(T) = bt(T) / (one - beta)
      end if

!     -------------
!     Normal entry.
!     -------------
!     Constraints first.
!     Only the first T components of F are nonlinear.

      neG = 0
      do i = 1, T
         j    = i
         xKap = x(j)
         Fi   = at(i)*xKap**b
         Gi   =  b*Fi/xKap

         if (needF .gt. 0) F(i) = Fi
         if (needG .gt. 0) then
            neG        = neG + 1
!           iGfun(neG) = i
!           jGvar(neG) = j
            G(neG)     = Gi
         end if
      end do

!     Set the objective element.

      FObj = zero
      do j = 1, T
         xCon = x(T+j)
         if (needF .gt. 0) then
            FObj = FObj  +  bt(j) * log(xCon)
!Min        FObj = FObj  -  bt(j) * log(xCon)
         end if

         if (needG .gt. 0) then
            neG    = neG + 1

!           iGfun(neG) = Obj
!           jGvar(neG) = T + j
            G(neG) = + bt(j) / xCon
!Min        G(neG) = - bt(j) / xCon
         end if
      end do

      F(Obj) = FObj

!     ------------
!     Final entry.
!     ------------
      if (Status .ge. 2) then
         if (Out .gt. 0) write(Out, 2000) (F(j), j = 1, T)
      end if
      return

 1000 format(// ' Starting  problem  Manne.   T =', i4, '   b =', f8.3/)
 2000 format(// ' Finishing problem  Manne.'/
     &          ' Final nonlinear function values' / (5f12.5))

      end ! of userfg
