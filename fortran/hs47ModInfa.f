!     ------------------------------------------------------------------
!     File hs47ModInfa.f : solves a modified version of hs47.
!     This is problem is hs47Modc with the bounds changed to give
!     an infeasible problem:
!
!     (1) the lower bounds on the two general linear constraints are
!         switched.
!     (2) lower bounds of 1.0 are imposed on x(1:n).
!
!     12 October 2014: First  version.
!     ------------------------------------------------------------------
      program
     &     hs47ModInf
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 30, maxn   =  20,
     &     lenA   = 50, lenG   = 100,
     &     nxname =  1, nFname =   1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character*8
     &     Prob, xnames(nxname), Fnames(nFname)
      double precision
     &     objAdd, sInf,
     &     A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)
      integer
     &     DerOpt, Errors, neA, neG, ObjRow, INFO, iPrt, iPrint, iSpecs,
     &     iSum, iSumm, Mjrlim, mincw, miniw, minrw, nF, n, nInf,
     &     nOut, nS
      logical
     &     byname
      integer
     &     lunit
      character*20
     &     lfile
      external
     &     hs47Modf, hs47Modfg, hs47Modfg2, hs47ModExf, hs47ModExfg
!     ------------------------------------------------------------------
!     SNOPT workspace

      integer               lenrw
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      integer               leniw
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      integer               lencw
      parameter          (  lencw =   600)
      character*8        cw(lencw)
!     ------------------------------------------------------------------
      integer             Cold,       Basis,      Warm
      parameter          (Cold   = 0, Basis  = 1, Warm  = 2)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!
!     nOut    is an output file used here by hs47Mod

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      nOut   =  6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lunit = iSpecs
         lfile = 'hs47ModInfa.spc'
         open( lunit, file=lfile, status='OLD',     err=800 )
         lunit = iPrint
         lfile = 'hs47ModInfa.out'
         open( lunit, file=lfile, status='UNKNOWN', err=800 )
      else

!        VMS  systems.  Define units for the Specs and print files.

         lunit = iSpecs
         open( lunit, status='OLD',     err=900 )
         lunit = iPrint
         open( lunit, status='UNKNOWN', err=900 )
      end if

!     ------------------------------------------------------------------
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      Errors = 0
      iPrt   = 0
      iSum   = 0

!     ------------------------------------------------------------------
!     1. Solve the modified hs47 with an infeasible constraint.
!        This is a test of implicit elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------'
      write(nOut, *) '1. Infeasible Modified hs47   '
      write(nOut, *) '   using implicit elastic mode'
      write(nOut, *) '   No derivatives supplied    '
      write(nOut, *) '------------------------------'

!     Set up the problem to be solved.
!     No derivatives are set in this case.

      call hs47fData
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     snOptA will compute the Jacobian by finite-differences.
!     The user has the option of calling  snJac  to define the
!     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).

      call snJac
     &   ( INFO, nF, n, hs47Modf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) then
         go to 920
      end if

!     ------------------------------------------------------------------
!     Warn snOptA that userf does not set the derivatives.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      DerOpt = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, hs47Modf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfa (f only) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .gt. 30) go to 920

!     ------------------------------------------------------------------
!     2. Solve the modified hs47 with an infeasible constraint.
!        This is a test of implicit elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------'
      write(nOut, *) '2. Infeasible Modified hs47   '
      write(nOut, *) '   using implicit elastic mode'
      write(nOut, *) '   Derivatives are supplied   '
      write(nOut, *) '------------------------------'

!     ------------------------------------------------------------------
!     Read a Specs file (Optional).
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

!     Set up the data structure for the sparse Jacobian.

      call hs47fgData
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      MjrLim = 250
      call snSeti
     &   ( 'Major Iterations', MjrLim, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Solve the problem again, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, hs47Modfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfa (f and g) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     3. Solve the same problem with the variables and  constraints in
!        nonlinear/linear order.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '--------------------------------------'
      write(nOut, *) '3. Infeasible Modified hs47           '
      write(nOut, *) '   Everything in nonlnear/linear order'
      write(nOut, *) '   Implicit elastic mode              '
      write(nOut, *) '   Derivatives  supplied              '
      write(nOut, *) '--------------------------------------'

      call hs47fgData2
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

      DerOpt = 1
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      MjrLim = 250
      call snSeti
     &   ( 'Major Iterations', MjrLim, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, hs47Modfg2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfa (f and g) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)

      if (INFO .ge. 30) go to 920

!     ------------------------------------------------------------------
!     4. Solve the modified hs47 with an infeasible constraint.
!        This is a test of implicit elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------'
      write(nOut, *) '4. Infeasible Modified hs47   '
      write(nOut, *) '   using explicit elastic mode'
      write(nOut, *) '   No derivatives supplied    '
      write(nOut, *) '------------------------------'

!     Set up the problem to be solved.
!     No derivatives are set in this case.

      call hs47ExfData
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     snOptA will compute the Jacobian by finite-differences.
!     The user has the option of calling  snJac  to define the
!     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).

      call snJac
     &   ( INFO, nF, n, hs47ModExf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) then
         go to 920
      end if

      DerOpt = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, hs47ModExf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfa Explicit (f only) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .gt. 30) go to 920

!     ------------------------------------------------------------------
!     5. Solve the modified hs47 with an infeasible constraint.
!        This is a test of implicit elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------'
      write(nOut, *) '5. Infeasible Modified hs47   '
      write(nOut, *) '   Explicit elastic mode      '
      write(nOut, *) '   Derivatives supplied       '
      write(nOut, *) '------------------------------'

!     Set up the data structure for the sparse Jacobian.

      call hs47ExfgData
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snSet
     &   ( 'Verify level 3', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Solution   Yes', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Solve the problem again, using a Cold start (Start = 0).
!     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     objAdd, ObjRow, Prob, hs47ModExfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfa Explicit (f and g) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, 4010) 'Error while opening unit', lunit
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'Insufficient space to hold the problem'
      stop

  920 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! program hsmain

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47fData
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
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     hs47fData defines a modified form of the problem HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 1
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     On exit:
!     Errors      is 0 if there is enough storage, 1 otherwise.
!     nF          is the number of problem functions
!                 (objective and constraints, linear and nonlinear).
!     n           is the number of variables.
!     xlow        holds the lower bounds on x.
!     xupp        holds the upper bounds on x.
!     Flow        holds the lower bounds on F.
!     Fupp        holds the upper bounds on F.

!     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n)     is a set of initial values for x.
!     Fmul(1:nF)  is a set of initial values for the dual variables.
!
!     24 Dec 1997: First version of hs47fData.
!     30 Oct 2002: Current version.
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'hs47f   '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6
      Obj    = 6                ! HS47 objective row
      ObjRow = 6                ! Can be 0
      n      = 5

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do i = 1, n
         xLow(i)   =  1.0d+0       ! was  - plInfy
         xupp(i)   =  plInfy
      end do

!     Nonlinear constraints first.
!     Followed by the linear constraints.
!     The Objective row is free.

      Flow(1)   =  1.0d+0       ! Linear constraint
      Fupp(1)   =  plInfy

      Flow(2)   =  3.0d+0
      Fupp(2)   =  3.0d+0

      Flow(3)   =  3.0d+0       ! Linear constraint
      Fupp(3)   =  plInfy

      Flow(4)   =  1.0d+0
      Fupp(4)   =  1.0d+0

      Flow(5)   =  1.0d+0
      Fupp(5)   =  1.0d+0

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      objAdd    = 0.0d+0

!     ------------------------------------------------------------------
!     Initialize x, xstate and Fmul.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  2.0d+0
      x(2)   =  sqrt(2.0d+0)
      x(3)   = -1.0d+0
      x(4)   =  2.0d+0 - sqrt(2.0d+0)
      x(5)   =  0.5d+0

      do i = 1, n
         xstate(i) =  0
      end do

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      end ! subroutine hs47fData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47Modf
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG,
     &     lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n),
     &     ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for a
!     modified version of HS47.
!     nF = 6, n = 5.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     ==================================================================
      integer
     &     Obj
!     ------------------------------------------------------------------
      Obj = 6

!     Constraints.

      F(1) =                           x(3)     + x(4) + x(5)
      F(2) =     x(1)      + x(2)**2 + x(3)**3
      F(3) =     x(1)      + x(2)
      F(4) =               + x(2)    - x(3)**2  + x(4)
      F(5) =     x(1)*x(5)

!     Objective

      F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3 +
     &          (x(3)-x(4))**4 + (x(4)-x(5))**4

      end ! subroutine hs47Modf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47fgData
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
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     hs47fgData defines a modified form of the problem HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     On exit:
!        nF   is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'HS47 (1)'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6
      Obj    = 6                ! HS47 objective row
      ObjRow = 6                ! Can be 0
      n      = 5

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF) Errors = 1
      if (n      .gt. maxn) Errors = 1
      if (Errors .gt.    0) return

!     ------------------------------------------------------------------
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The pattern of the Jacobian is as follows, where
!     L = constant element, N = nonlinear element.
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |          L    L    L
!         2  | L   N    N
!         3  | L   L
!         4  |     L    N    L
!         5  | N                  N
!     row 6  | N   N    N    N    N    Objective row
!
!
!     First we assign the list of varying derivative entries.
!     The pattern of nonlinear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |
!         2  |     6    7
!         3  |
!         4  |          8
!         5  | 9                 10
!     row 6  | 1   2    3    4    5    Objective row
!
!
!     ------------------------------------------------------------------
!     Nonlinear Objective derivatives

      neG        =  1
!     G(neG)     =  2.0d+0   * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * (x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(3)-x(4))**3 - 3.0d+0 * (x(2)-x(3))**2
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(4)-x(5))**3 -  4.0d+0 * (x(3)-x(4))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     = -4.0d+0  * (x(4)-x(5))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  5

!     Nonlinear constraints.

      neG        =  neG + 1
!     G(neG)     =  2.0d+0   * x(2)
      iGfun(neG) =  2
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * x(3)**2
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     = -2.0d+0   * x(3)
      iGfun(neG) =  4
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  x(5)
      iGfun(neG) =  5
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  x(1)
      iGfun(neG) =  5
      jGvar(neG) =  5

!     neG        = 10

!     ------------------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     The pattern of linear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |          4    5    6
!         2  | 1
!         3  | 7   8
!         4  |     2         3
!         5  |
!     row 6  |                         Objective row
!
!     ------------------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  4
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  3
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  4
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  5
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

!     neA        =  8

!     ----------------
!     Initial x.
!     ----------------
      objAdd = 0.0d+0

      x(1)   =  2.0d+0
      x(2)   =  sqrt(2.0d+0)
      x(3)   = -1.0d+0
      x(4)   =  2.0d+0 - sqrt(2.0d+0)
      x(5)   =  0.5d+0

      do i = 1, n
         xLow(i)   =  1.0d+0       ! was  -infBnd
         xUpp(i)   =  plInfy
      end do

      Flow(1)   =  1.0d+0       ! Linear constraint
      Fupp(1)   =  plInfy

      Flow(2)   =  3.0d+0
      Fupp(2)   =  3.0d+0

      Flow(3)   =  3.0d+0       ! Linear constraint
      Fupp(3)   =  plInfy

      Flow(4)   =  1.0d+0
      Fupp(4)   =  1.0d+0

      Flow(5)   =  1.0d+0
      Fupp(5)   =  1.0d+0

      Flow(Obj) = -plInfy       ! Objective row
      Fupp(Obj) =  plInfy

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      do i = 1, n
         xstate(i) =  0
      end do

      end ! subroutine hs47fgData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47Modfg
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
!     Computes the nonlinear objective and constraint terms for a
!     modified (and infeasible version of HS47.
!     nF = 6, n = 5.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define
!     the sparsity pattern and values of the nonlinear elements
!     of the Jacobian.
!     ==================================================================
      integer
     &     neG, Obj
!     ------------------------------------------------------------------
      Obj = 6

      if (needF .gt. 0) then

!        Constraints.  Only the nonlinear terms are required.
!        The linear terms are defined elsewhere  via the triples
!          (A(k),iAfun(k),jAvar(k)), k = 1,2,...,neA,

         F(2) =     x(2)*x(2) + x(3)*x(3)*x(3)
         F(4) = -   x(3)*x(3)
         F(5) =     x(1)*x(5)

!        Objective.  Only the nonlinear terms are required.

         F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3
     &          +  (x(3)-x(4))**4 + (x(4)-x(5))**4
      end if

      neG =  0

      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  2.0d+0  *(x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  3.0d+0*(x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(3)-x(4))**3 - 3.0d+0 * (x(2)-x(3))**2
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(4)-x(5))**3 -  4.0d+0 * (x(3)-x(4))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     = -4.0d+0 *(x(4)-x(5))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  5

!        Nonlinear constraints.

         neG        =  neG + 1
         G(neG)     =  2.0d+0  *x(2)
!        iGfun(neG) =  2
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  3.0d+0 * x(3)**2
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     = -2.0d+0   * x(3)
!        iGfun(neG) =  4
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  x(5)
!        iGfun(neG) =  5
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  x(1)
!        iGfun(neG) =  5
!        jGvar(neG) =  5

!        neG        = 10
      end if

      end ! subroutine hs47Modfg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47fgData2
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
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     hs47fgData2 defines a modified form of the problem HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1)      + x(2)                            >= 3
!                                       x(3)    + x(4) + x(5) >= 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     On exit:
!        nF   is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'HS47 (2)'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6
      Obj    = 6                ! HS47 objective row
      ObjRow = 6                ! Can be 0
      n      = 5

      objAdd = 0.0d+0

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF) Errors = 1
      if (n      .gt. maxn) Errors = 1
      if (Errors .gt.    0) return

!     ------------------------------------------------------------------
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(5))**4
!                                                 + (x(5)-x(4))**4
!
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(5)         = 1
!                 x(1)*x(4)                                    = 1
!                 x(1)      + x(2)                            >= 1
!                                       x(3)    + x(4) + x(5) >= 3
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The pattern of the Jacobian is as follows, where
!     L = constant element, N = nonlinear element.
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  | L   N    N
!         2  |     L    N         L
!         3  | N             N
!         4  | L   L
!         5  |          L    L    L
!     row 6  | N   N    N    N    N    Objective row
!
!
!     First we assign the list of varying derivative entries.
!     The pattern of nonlinear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |
!         2  |     6    7
!         3  |
!         4  |          8
!         5  | 9                 10
!     row 6  | 1   2    3    4    5    Objective row
!
!
!     ------------------------------------------------------------------
!     Nonlinear Objective derivatives

      neG        =  1
!     G(neG)     =  2.0d+0 * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * (x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(3)-x(5))**3 - 3.0d+0 * (x(2)-x(3))**2
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     = -4.0d+0  * (x(5)-x(4))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(5)-x(4))**3 -  4.0d+0 * (x(3)-x(5))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  5

!     Nonlinear constraints.

      neG        =  neG + 1
!     G(neG)     =  2.0d+0   * x(2)
      iGfun(neG) =  1
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * x(3)**2
      iGfun(neG) =  1
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     = -2.0d+0   * x(3)
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  x(4)
      iGfun(neG) =  3
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  x(1)
      iGfun(neG) =  3
      jGvar(neG) =  4

!     neG        = 10

!     ------------------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     The pattern of linear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |          4    5    6
!         2  | 1
!         3  | 7   8
!         4  |     2         3
!         5  |
!     row 6  |                         Objective row
!
!     ------------------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  5
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  5
      jAvar(neA) =  3
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  5
      jAvar(neA) =  4
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  5
      jAvar(neA) =  5
      A(neA)     =  1.0d+0

!     neA        =  8

!     ----------------
!     Initial x.
!     ----------------
      x(1)   =  2.0d+0
      x(2)   =  sqrt(2.0d+0)
      x(3)   = -1.0d+0
      x(4)   =  2.0d+0 - sqrt(2.0d+0)
      x(5)   =  0.5d+0

      do i = 1, n
         xLow(i)   =  1.0d+0       ! was  -infBnd
         xUpp(i)   =  plInfy
      end do

      Flow(1)   =  3.0d+0
      Fupp(1)   =  3.0d+0

      Flow(2)   =  1.0d+0
      Fupp(2)   =  1.0d+0

      Flow(3)   =  1.0d+0
      Fupp(3)   =  1.0d+0

      Flow(4)   =  3.0d+0       ! Linear constraint
      Fupp(4)   =  plInfy

      Flow(5)   =  1.0d+0       ! Linear constraint
      Fupp(5)   =  plInfy

      Flow(Obj) = -plInfy       ! Objective row
      Fupp(Obj) =  plInfy

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      do i = 1, n
         xstate(i) =  0
      end do

      end ! subroutine hs47fgData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47Modfg2
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
!     Computes the nonlinear objective and constraint terms for a
!     modified (and infeasible) version of HS47.
!     nF = 6, n = 5.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(5))**4
!                                                 + (x(5)-x(4))**4
!
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(5)         = 1
!                 x(1)*x(4)                                    = 1
!                 x(1)      + x(2)                            >= 3
!                                       x(3)    + x(4) + x(5) >= 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define
!     the sparsity pattern and values of the nonlinear elements
!     of the Jacobian.
!     ==================================================================
      integer
     &     neG, Obj
!     ------------------------------------------------------------------
      Obj = 6

      if (needF .gt. 0) then

!        Constraints.  Only the nonlinear terms are required.
!        The linear terms are defined elsewhere  via the triples
!          (A(k),iAfun(k),jAvar(k)), k = 1,2,...,neA,

         F(1)   =  x(2)*x(2) + x(3)*x(3)*x(3)
         F(2)   = -x(3)*x(3)
         F(3)   =  x(1)*x(4)

!        Objective.  Only the nonlinear terms are required.

         F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3
     &          +  (x(3)-x(5))**4 + (x(5)-x(4))**4
      end if

      neG =  0

      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  2.0d+0  *(x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  3.0d+0*(x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(3)-x(5))**3 - 3.0d+0 * (x(2)-x(3))**2
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     = -4.0d+0 *(x(5)-x(4))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(5)-x(4))**3 -  4.0d+0 * (x(3)-x(5))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  5

!        Nonlinear constraints.

         neG        =  neG + 1
         G(neG)     =  2.0d+0  *x(2)
!        iGfun(neG) =  1
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  3.0d+0 * x(3)**2
!        iGfun(neG) =  1
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     = -2.0d+0   * x(3)
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  x(4)
!        iGfun(neG) =  3
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  x(1)
!        iGfun(neG) =  3
!        jGvar(neG) =  4

!        neG        = 10
      end if

      end ! subroutine hs47Modfg2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ExfData
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
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     hs47fData defines a modified form of the problem HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 1
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     Given a penalty parameter penParm, explicit elastic variables
!     u(1:3), v(1:3) are added to give the elastic problem:
!
!     Minimize    f(x) + penParm*(u + v)
!
!     subject to
!                   (             x         )
!                   (             u         )
!            bl le  (             v         ) le bu
!                   ( fCon(x) + J x - u + v )
!                   (            Ax         )
!
!     where the Jacobian for fCon(x) + J x - u + v is stored in
!     JCon(ldJ,*), with with dimensions  mNCon by n.  The elements of
!     the Jacobian of fCon comprise the first nnJac columns JCon.
!
!     The elastic variables u and v are held in x(6:11).
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     On exit:
!     Errors      is 0 if there is enough storage, 1 otherwise.
!     nF          is the number of problem functions
!                 (objective and constraints, linear and nonlinear).
!     n           is the number of variables.
!     xlow        holds the lower bounds on x.
!     xupp        holds the upper bounds on x.
!     Flow        holds the lower bounds on F.
!     Fupp        holds the upper bounds on F.

!     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n)     is a set of initial values for x.
!     Fmul(1:nF)  is a set of initial values for the dual variables.
!
!     24 Dec 1997: First version of hs47fData.
!     30 Oct 2002: Current version.
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'hs47f Ex'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6
      Obj    = 6                ! objective row
      ObjRow = 6                ! Can be 0
      n      = 11

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do i = 1, 5
         xLow(i)   =  1.0d+0       ! was  - plInfy
         xupp(i)   =  plInfy
      end do

      do i = 6, n
         xLow(i) =  0.0d+0
         xUpp(i) =  plInfy
      end do

!     Nonlinear constraints first.
!     Followed by the linear constraints.
!     The Objective row is free.

      Flow(1)   =  1.0d+0       ! Linear constraint
      Fupp(1)   =  plInfy

      Flow(2)   =  3.0d+0
      Fupp(2)   =  3.0d+0

      Flow(3)   =  3.0d+0       ! Linear constraint
      Fupp(3)   =  plInfy

      Flow(4)   =  1.0d+0
      Fupp(4)   =  1.0d+0

      Flow(5)   =  1.0d+0
      Fupp(5)   =  1.0d+0

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      objAdd    = 0.0d+0

!     ------------------------------------------------------------------
!     Initialize x, xstate and Fmul.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x( 1)   =  2.0d+0
      x( 2)   =  sqrt(2.0d+0)
      x( 3)   = -1.0d+0
      x( 4)   =  2.0d+0 - sqrt(2.0d+0)
      x( 5)   =  0.5d+0
      x( 6)   =  0.0d+0
      x( 7)   =  0.0d+0
      x( 8)   =  0.0d+0
      x( 9)   =  0.0d+0
      x(10)   =  0.0d+0
      x(11)   =  0.0d+0

      do i = 1, n
         xstate(i) =  0
      end do

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      end ! subroutine hs47ExfData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ModExf
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG,
     &     lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n),
     &     ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for a
!     modified version of HS47.
!     nF = 6, n = 11.
!
!     Minimize    f(x) + penParm*(u + v)
!
!     subject to
!                   (             x         )
!                   (             u         )
!            bl le  (             v         ) le bu
!                   ( fCon(x) + J x - u + v )
!                   (            Ax         )
!
!     where the Jacobian for fCon(x) + J x - u + v is stored in
!     JCon(ldJ,*), with with dimensions  mNCon by n.  The elements of
!     the Jacobian of fCon comprise the first nnJac columns JCon.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     ==================================================================
      integer
     &     Obj
!     ------------------------------------------------------------------
      double precision   wtInf
      parameter         (wtInf  = 1.0d+4)
!     ------------------------------------------------------------------
      Obj = 6

!     Constraints.

      F(1) =                           x(3)     + x(4) + x(5)
      F(2) =     x(1)      + x(2)**2 + x(3)**3
     &         - x(6)                + x(9)
      F(3) =     x(1)      + x(2)
      F(4) =               + x(2)    - x(3)**2  + x(4)
     &                - x(7)                + x(10)
      F(5) =     x(1)*x(5)
     &                     - x(8)              + x(11)

!     Objective

      F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3 +
     &          (x(3)-x(4))**4 + (x(4)-x(5))**4
     &         + wtInf*(x(6) + x(7) + x(8) + x(9) + x(10) + x(11))

      end ! subroutine hs47ModExf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ExfgData
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
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     hs47ExfgData defines a modified form of the problem HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     On exit:
!        nF   is the number of objective and constraint functions
!               (including linear and nonlinear)
!        n    is the number of variables.
!
!        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
!             of the nonzero problem derivatives.
!             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
!             of the problem vector F(i), i = 0,1,2,...,nF,  with
!             objective function in position 0 and constraint functions
!             in positions  1  through  m.
!
!        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
!             of the nonzero constant problem derivatives.
!
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision   wtInf
      parameter         (wtInf  = 1.0d+4)
      double precision  plInfy
      parameter         (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'hs47Exfg'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 6
      Obj    = 6                ! objective row
      ObjRow = 6                ! Can be 0
      n      = 11

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

!     ------------------------------------------------------------------
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The pattern of the Jacobian is as follows, where
!     L = constant element, N = nonlinear element.
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |          L    L    L
!         2  | L   N    N
!         3  | L   L
!         4  |     L    N    L
!         5  | N                  N
!     row 6  | N   N    N    N    N    Objective row
!
!
!     First we assign the list of varying derivative entries.
!     The pattern of nonlinear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |
!         2  |     6    7
!         3  |
!         4  |          8
!         5  | 9                 10
!     row 6  | 1   2    3    4    5    Objective row
!
!
!     ------------------------------------------------------------------
!     Nonlinear Objective derivatives

      neG        =  1
!     G(neG)     =  2.0d+0   * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * (x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
      iGfun(neG) =  Obj
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(3)-x(4))**3 - 3.0d+0 * (x(2)-x(3))**2
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  4.0d+0  * (x(4)-x(5))**3 -  4.0d+0 * (x(3)-x(4))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  4

      neG        =  neG + 1
!     G(neG)     = -4.0d+0  * (x(4)-x(5))**3
      iGfun(neG) =  Obj
      jGvar(neG) =  5

!     Nonlinear constraints.

      neG        =  neG + 1
!     G(neG)     =  2.0d+0   * x(2)
      iGfun(neG) =  2
      jGvar(neG) =  2

      neG        =  neG + 1
!     G(neG)     =  3.0d+0 * x(3)**2
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     = -2.0d+0   * x(3)
      iGfun(neG) =  4
      jGvar(neG) =  3

      neG        =  neG + 1
!     G(neG)     =  x(5)
      iGfun(neG) =  5
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  x(1)
      iGfun(neG) =  5
      jGvar(neG) =  5

!     neG        = 10

!     ------------------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     The pattern of linear elements is as follows:
!
!              Column
!            | 1   2    3    4    5   6   7   8   9   10   11
!            +----------------------
!         1  |          4    5    6
!         2  | 1                      9          10
!         3  | 7   8
!         4  |     2         3           11           12
!         5  |                               13            14
!     row 6  |                       15  16  17  18   19   20
!
!     ------------------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  4
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  3
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  4
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  1
      jAvar(neA) =  5
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  1
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  2
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  6
      A(neA)     = -1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  9
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  7
      A(neA)     = -1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) = 10
      A(neA)     =  1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  5
      jAvar(neA) =  8
      A(neA)     = -1.0d+0

      neA        =  neA + 1
      iAfun(neA) =  5
      jAvar(neA) = 11
      A(neA)     =  1.0d+0

!!
      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =   6
      A(neA)     =  wtInf

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =   7
      A(neA)     =  wtInf

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =   8
      A(neA)     =  wtInf

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =   9
      A(neA)     =  wtInf

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =  10
      A(neA)     =  wtInf

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =  11
      A(neA)     =  wtInf

!     neA        = 20

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do i = 1, 5
         xLow(i)   =  1.0d+0       ! was  - plInfy
         xupp(i)   =  plInfy
      end do

      do i = 6, n
         xLow(i) =  0.0d+0
         xUpp(i) =  plInfy
      end do

!     Nonlinear constraints first.
!     Followed by the linear constraints.
!     The Objective row is free.

      Flow(1)   =  1.0d+0       ! Linear constraint
      Fupp(1)   =  plInfy

      Flow(2)   =  3.0d+0
      Fupp(2)   =  3.0d+0

      Flow(3)   =  3.0d+0       ! Linear constraint
      Fupp(3)   =  plInfy

      Flow(4)   =  1.0d+0
      Fupp(4)   =  1.0d+0

      Flow(5)   =  1.0d+0
      Fupp(5)   =  1.0d+0

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      objAdd    = 0.0d+0

!     ------------------------------------------------------------------
!     Initialize x, xstate and Fmul.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x( 1)   =  2.0d+0
      x( 2)   =  sqrt(2.0d+0)
      x( 3)   = -1.0d+0
      x( 4)   =  2.0d+0 - sqrt(2.0d+0)
      x( 5)   =  0.5d+0
      x( 6)   =  0.0d+0
      x( 7)   =  0.0d+0
      x( 8)   =  0.0d+0
      x( 9)   =  0.0d+0
      x(10)   =  0.0d+0
      x(11)   =  0.0d+0

      do i = 1, n
         xstate(i) =  0
      end do

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      end ! subroutine hs47ExfData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ModExfg
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG,
     &     lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n),
     &     ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for a
!     modified version of HS47.
!     nF = 6, n = 11.
!
!     Minimize    f(x) + penParm*(u + v)
!
!     subject to
!                   (             x         )
!                   (             u         )
!            bl le  (             v         ) le bu
!                   ( fCon(x) + J x - u + v )
!                   (            Ax         )
!
!     where the Jacobian for fCon(x) + J x - u + v is stored in
!     JCon(ldJ,*), with with dimensions  mNCon by n.  The elements of
!     the Jacobian of fCon comprise the first nnJac columns JCon.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to                        x(3)    + x(4) + x(5) >= 3
!                 x(1)      + x(2)**2 + x(3)**3                = 3
!                 x(1)      + x(2)                            >= 1
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!                 x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     ==================================================================
      integer
     &     neG, Obj
!     ------------------------------------------------------------------
      double precision   wtInf
      parameter         (wtInf  = 1.0d+4)
!     ------------------------------------------------------------------
      Obj = 6

!     Constraints.

      if (needF .gt. 0) then

!        Constraints.  Only the nonlinear terms are required.
!        The linear terms are defined elsewhere  via the triples
!          (A(k),iAfun(k),jAvar(k)), k = 1,2,...,neA,

         F(2) =     x(2)*x(2) + x(3)*x(3)*x(3)
         F(4) = -   x(3)*x(3)
         F(5) =     x(1)*x(5)

!        Objective.  Only the nonlinear terms are required.

         F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3
     &          +  (x(3)-x(4))**4 + (x(4)-x(5))**4
      end if

      neG =  0

      if (needG .gt. 0) then
         neG        =  neG + 1
         G(neG)     =  2.0d+0  *(x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  3.0d+0*(x(2)-x(3))**2 - 2.0d+0   * (x(1)-x(2))
!        iGfun(neG) =  Obj
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(3)-x(4))**3 - 3.0d+0 * (x(2)-x(3))**2
!        iGfun(neG) =  Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  4.0d+0 *(x(4)-x(5))**3 -  4.0d+0 * (x(3)-x(4))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     = -4.0d+0 *(x(4)-x(5))**3
!        iGfun(neG) =  Obj
!        jGvar(neG) =  5

!        Nonlinear constraints.

         neG        =  neG + 1
         G(neG)     =  2.0d+0  *x(2)
!        iGfun(neG) =  2
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  3.0d+0 * x(3)**2
!        iGfun(neG) =  2
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     = -2.0d+0   * x(3)
!        iGfun(neG) =  4
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  x(5)
!        iGfun(neG) =  5
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  x(1)
!        iGfun(neG) =  5
!        jGvar(neG) =  5

!        neG        = 10
      end if

      end ! subroutine hs47ModExfg

