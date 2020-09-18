!     ------------------------------------------------------------------
!     File hs47a.f
!     This is a main program to illustrate the snOptA,  which is part
!     of the SNOPT 6 package.
!
!     Includes hs47.
!
!     31 Jul 1996: First   version.
!     25 Sep 1999: Updated for SNOPT 6
!     30 Oct 2002: Current version.
!     ------------------------------------------------------------------
      program
     &     hsmain
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 30, maxn   =  10,
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
     &     hs47f, hs47fg
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
!     nOut    is an output file used here by hs47

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      nOut   =  6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lunit = iSpecs
         lfile = 'hs47a.spc'
         open( lunit, file=lfile, status='OLD',     err=800 )
         lunit = iPrint
         lfile = 'hs47a.out'
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

!     Set up the problem to be solved.
!     No derivatives are set in this case.

      Errors = 0
      call hs47fData
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     objAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     snOptA will compute the Jacobian by finite-differences.
!     The user has the option of calling  snJac  to define the
!     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).

      call snJac
     &   ( INFO, nF, n, hs47f,
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
     &     objAdd, ObjRow, Prob, hs47f,
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
      write(nOut, *) 'hs47a (f only) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', objAdd + F(ObjRow)
      if (INFO .gt. 30) go to 920

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
     &     objAdd, ObjRow, Prob, hs47fg,
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
      write(nOut, *) 'hs47a (f and g) finished.'
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
!     hs47fData defines the problem hs47a.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
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
!     24 Dec 1997: First version of hs47d0.
!     30 Oct 2002: Current version.
!     ==================================================================
      integer
     &     i, Obj
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'hs47a   '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      Obj    = 4                ! HS47 objective row
      ObjRow = 4                ! Can be 0
      n      = 5

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Nonlinear constraints first.
!     Followed by the linear constraints.
!     The Objective row is free.

      Flow(1)   =  3.0d+0
      Fupp(1)   =  3.0d+0
      Flow(2)   =  1.0d+0
      Fupp(2)   =  1.0d+0
      Flow(3)   =  1.0d+0
      Fupp(3)   =  1.0d+0
      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      objAdd = 0.0d+0

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do i = 1, n
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
      end do

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

      subroutine hs47f
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
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!
!     ==================================================================
      integer
     &     Obj
!     ------------------------------------------------------------------
      Obj = 4

!     Constraints.

      F(1) =     x(1)      + x(2)**2 + x(3)**3
      F(2) =               + x(2)    - x(3)**2  + x(4)
      F(3) =     x(1)*x(5)

!     Objective

      F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3
     &       +  (x(3)-x(4))**4 + (x(4)-x(5))**4

      end ! subroutine hs47f

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
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
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

      Prob   = 'hs47a   '

!     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      Obj    = 4                ! HS47 objective row
      ObjRow = 4                ! Can be 0
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
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!
!     The pattern of the Jacobian is as follows, where
!     L = constant element, N = nonlinear element.
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  | L   N    N
!         2  |     L    N    L
!         3  | N                  N
!     row 4  | N   N    N    N    N    Objective row
!
!
!     First we assign the list of varying derivative entries.
!     The pattern of nonlinear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  |     6    7
!         2  |          8
!         3  | 9                 10
!     row 4  | 1   2    3    4    5    Objective row
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

!     ----------------------
!     Nonlinear constraints.
!     ----------------------
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
!     G(neG)     =  x(5)
      iGfun(neG) =  3
      jGvar(neG) =  1

      neG        =  neG + 1
!     G(neG)     =  x(1)
      iGfun(neG) =  3
      jGvar(neG) =  5

!     neG        = 10

!     ------------------------------------------------------------------
!     Next we assign the list of constant derivative entries.
!     The pattern of linear elements is as follows:
!
!              Column
!            | 1   2    3    4    5
!            +----------------------
!         1  | 1
!         2  |     2         3
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
      jAvar(neA) =  4
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
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      Flow(1)   =  3.0d+0
      Fupp(1)   =  3.0d+0
      Flow(2)   =  1.0d+0
      Fupp(2)   =  1.0d+0
      Flow(3)   =  1.0d+0
      Fupp(3)   =  1.0d+0
      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      do i = 1, nF
         Fmul(i) = 0.0d+0
      end do

      end ! subroutine hs47fgData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47fg
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
!     modified version of HS47.
!     nF = 6, n = 5.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to  x(1)      + x(2)**2 + x(3)**3                = 3
!                             x(2)    - x(3)**2 + x(4)         = 1
!                 x(1)*x(5)                                    = 1
!
!     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define
!     the sparsity pattern and values of the nonlinear elements
!     of the Jacobian.
!     ==================================================================
      integer
     &     neG, Obj
!     ------------------------------------------------------------------
      Obj = 4

      if (needF .gt. 0) then

!        Constraints.  Only the nonlinear terms are required.
!        The linear terms are defined elsewhere  via the triples
!          (A(k),iAfun(k),jAvar(k)), k = 1,2,...,neA,

         F(1) =     x(2)*x(2) + x(3)*x(3)*x(3)
         F(2) = -   x(3)*x(3)
         F(3) =     x(1)*x(5)

!        Objective.  Only the nonlinear terms are required.

         F(Obj) =  (x(1)-x(2))**2 + (x(2)-x(3))**3 +
     &             (x(3)-x(4))**4 + (x(4)-x(5))**4
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

!        ----------------------
!        Nonlinear constraints.
!        ----------------------
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
         G(neG)     =  x(5)
!        iGfun(neG) =  3
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  x(1)
!        iGfun(neG) =  3
!        jGvar(neG) =  5

!        neG        = 10
      end if

      end ! subroutine hs47fg

