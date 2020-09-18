!     ------------------------------------------------------------------
!     File sntoyA.f
!     This is a main program to illustrate the interface snOptA,
!     part of the SNOPT 7 package.
!
!     19 Oct 2002: First version for SNOPT 6
!     19 Oct 2003: Updated for SNOPT 7
!     19 Oct 2003: Current version.
!     ------------------------------------------------------------------
      program
     &     sntoyA

      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 10,
     &     maxn   = 10,
     &     lenA   = 10, lenG   = 10,
     &     nxname =  1, nFname =  1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character*8
     &     Prob, xnames(nxname), Fnames(nFname)
      double precision
     &     ObjAdd, sInf, A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)
      integer
     &     DerOpt, Errors, neA, neG, ObjRow, INFO, iPrt, iPrint, iSpecs,
     &     iSum, iSumm, Major, mincw, miniw, minrw, nF, n, nInf, nOut,
     &     nS
      logical
     &     byname
      integer
     &     lunit
      character*20
     &     lfile
      external
     &     userf, userfg
!     ------------------------------------------------------------------
!     SNOPT workspace

      integer               lenrw
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)

      integer               leniw
      parameter          (  leniw = 10000)
      integer            iw(leniw)

      integer               lencw
      parameter          (  lencw =   500)
      character*8        cw(lencw)

      integer             Cold,       Basis,      Warm
      parameter          (Cold   = 0, Basis  = 1, Warm  = 2)
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!
!     nOut    is an output file used here by hs106.

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      nOut   =  6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lunit = iSpecs
         lfile = 'sntoya.spc'
         open( lunit, file=lfile, status='OLD',     err=800 )
         lunit = iPrint
         lfile = 'sntoya.out'
         open( lunit, file=lfile, status='UNKNOWN', err=800 )
      else

!        VMS  systems.  Define units for the Specs and print files.

         lunit = iSpecs
         open( lunit, status='OLD',     err=900 )
         lunit = iPrint
         open( lunit, status='UNKNOWN', err=900 )
      end if

!     ==================================================================
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!     ==================================================================
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     1. Solve sntoya.
!        No derivatives are supplied.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------------'
      write(nOut, *) ' 1. Solve sntoya (f). No derivatives supplied'
      write(nOut, *) '---------------------------------------------'

      call Toy0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     ObjAdd, ObjRow, x, xstate, xlow, xupp,
     &     F, Flow, Fupp, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     snOptA will compute the Jacobian by finite-differences.
!     The user has the option of calling  snJac  to define the
!     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
!
      call snJac
     &   ( INFO, nF, n, userf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) go to 920

!     ------------------------------------------------------------------
!     Warn snOptA that userf does not set the derivatives.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors = 0

      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
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
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'sntoya (f) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd +  F(objRow)
      if (INFO .gt. 30) go to 920

!     ------------------------------------------------------------------
!     2. Solve sntoya with derivatives specified.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '--------------------------------------------'
      write(nOut, *) ' 2. Solve sntoya (fg) derivatives specified.'
      write(nOut, *) '--------------------------------------------'

!     ------------------------------------------------------------------
!     Read a Specs file (Optional).
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      call Toy1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, x, xstate, xlow, xupp,
     &     F, Flow, Fupp, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      Major  = 250
      call snSeti
     &   ( 'Major Iteration limit', Major, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Solve the problem again, this time with derivatives specified.
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
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'sntoya (fg) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd +  F(objRow)
      if (INFO .gt. 30) go to 920
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

      end ! main program

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Toy0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     ObjAdd, ObjRow, x, xstate, xlow, xupp,
     &     F, Flow, Fupp, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, nF, n, ObjRow, lencw, leniw, lenrw,
     &     Fstate(maxF), xstate(maxn), iw(leniw)
      double precision
     &     ObjAdd,
     &     xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     Toy0   defines input data for the toy problem discussed in the
!     snOptA Users Guide.
!
!        Minimize                      x(2)
!
!        subject to   x(1)**2      + 4 x(2)**2  <= 4,
!                    (x(1) - 2)**2 +   x(2)**2  <= 5,
!                     x(1) >= 0.
!
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
!     ==================================================================
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'sntoya f'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 3
      ObjRow = 1                !  = 0 for a feasible point.
      n      = 2

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.    0 ) return

      ObjAdd = 0.0d+0

!     -------------------------------
!     Set the upper and lower bounds.
!     -------------------------------
      do j = 1, n
         xlow(j)   = -plInfy
         xupp(j)   =  plInfy
         xstate(j) =  0
      end do

      xlow(1) =  0.0d+0

      do i = 1, nF
         Fstate(i) = 0
         Flow(i)   = -plInfy
         Fupp(i)   =  plInfy
         Fmul(i)   =  0.0d+0
         F(i)      =  0.0d+0
      end do

      Fupp(2)   =  4.0d+0
      Fupp(3)   =  5.0d+0

!     ----------------
!     Initialize  x.
!     ----------------
      x(1)   = 1.0d+0
      x(2)   = 1.0d+0

      end ! subroutine Toy0

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
!     Computes the nonlinear objective and constraint terms for the toy
!     problem featured in the snOptA users guide.
!     nF = 3, n = 2.
!
!        Minimize                      x(2)
!
!        subject to   x(1)**2      + 4 x(2)**2  <= 4,
!                    (x(1) - 2)**2 +   x(2)**2  <= 5,
!                     x(1) >= 0.
!
!     ==================================================================
      integer
     &     Out
!     ------------------------------------------------------------------
      Out = 6
!     --------------------------------------------
!     Print something on the first and last entry.
!     -------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a/)') ' Starting  sntoya (f)'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a/)') ' Finishing sntoya (f)'
         return
      end if

      F(1) =                           x(2)     ! The objective row
      F(2) =   x(1)**2        + 4.0d+0*x(2)**2
      F(3) =  (x(1) - 2.0d+0)**2 +     x(2)**2

      end ! subroutine userf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Toy1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, x, xstate, xlow, xupp,
     &     F, Flow, Fupp, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, Fstate(maxF), xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

!     ==================================================================
!     Toy1   defines input data for the toy problem discussed in the
!     snOptA Users Guide.
!
!        Minimize                      x(2)
!
!        subject to   x(1)**2      + 4 x(2)**2  <= 4,
!                    (x(1) - 2)**2 +   x(2)**2  <= 5,
!                     x(1) >= 0.
!
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
!             To keep things simple, no constant elements are set here.
!
!     ==================================================================
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'sntoya g'

!     Assign the dimensions of the constraint Jacobian.

      nF     = 3
      ObjRow = 1                !  = 0 for a feasible point.
      n      = 2

!     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF) Errors = 1
      if (n      .gt. maxn) Errors = 1
      if (Errors .gt.    0) return

!     ------------------------------------------------------------------
      neG        = 0

      neG        = neG + 1
      iGfun(neG) = 1
      jGvar(neG) = 1
!     g(neG)     = 0.0d+0

      neG        = neG + 1
      iGfun(neG) = 1
      jGvar(neG) = 2
!     g(neG)     = 1.0d+0

      neG        = neG + 1
      iGfun(neG) = 2
      jGvar(neG) = 1
!     g(neG)     = 2.0d+0* x(1)

      neG        = neG + 1
      iGfun(neG) = 2
      jGvar(neG) = 2
!     g(neG)     = 8.0d+0*x(2)

      neG        = neG + 1
      iGfun(neG) = 3
      jGvar(neG) = 1
!     g(neG)     = 2.0d+0*(x(1) - 2.0d+0)

      neG        = neG + 1
      iGfun(neG) = 3
      jGvar(neG) = 2
!     g(neG)     = 2.0d+0*x(2)

!     neG        = 6 derivatives in all

      neA        =  0

      ObjAdd = 0.0d+0

!     -------------------------------
!     Set the upper and lower bounds.
!     -------------------------------
      do j = 1, n
         xlow(j)   = -plInfy
         xupp(j)   =  plInfy
         xstate(j) =  0
      end do

      xlow(1) =  0.0d+0

      do i = 1, nF
         Fstate(i) = 0
         Flow(i)   = -plInfy
         Fupp(i)   =  plInfy
         Fmul(i)   =  0.0d+0
         F(i)      =  0.0d+0
      end do

      Fupp(2)   =  4.0d+0
      Fupp(3)   =  5.0d+0

!     ----------------
!     Initialize  x.
!     ----------------
      x(1)   = 1.0d+0
      x(2)   = 1.0d+0

      end ! subroutine Toy1

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
!     Computes the nonlinear objective and constraint terms for the toy
!     problem featured in the snOptA users guide.
!     nF = 3, n = 2.
!
!        Minimize                      x(2)
!
!        subject to   x(1)**2      + 4 x(2)**2  <= 4,
!                    (x(1) - 2)**2 +   x(2)**2  <= 5,
!                     x(1) >= 0.
!
!     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define
!     the sparsity pattern and values of the nonlinear elements
!     of the Jacobian.
!     ==================================================================
      integer
     &     neG, nOut
!     ------------------------------------------------------------------
      nOut = 6

!     --------------------------------------------
!     Print something on the first and last entry.
!     -------------------------------------------
      if (Status .eq. 1) then       ! First
         if (nOut .gt. 0) write(nOut, '(/a)') ' Starting  sntoya (fg)'
      else  if (Status .ge. 2) then ! Last
         if (nOut .gt. 0) write(nOut, '(/a)') ' Finishing sntoya (fg)'
         return
      end if

      if (needF .gt. 0) then
         F(1) =                           x(2) ! The objective row
         F(2) =   x(1)**2        + 4.0d+0*x(2)**2
         F(3) =  (x(1) - 2.0d+0)**2 +     x(2)**2
      end if

      if (needG .gt. 0) then
         neG        =  0

         neG        = neG + 1
!        iGfun(neG) = 1
!        jGvar(neG) = 1
         g(neG)     = 0.0d+0

         neG        = neG + 1
!        iGfun(neG) = 1
!        jGvar(neG) = 2
         g(neG)     =  1.0d+0

         neG        = neG + 1
!        iGfun(neG) = 2
!        jGvar(neG) = 1
         g(neG)     =  2.0d+0* x(1)

         neG        = neG + 1
!        iGfun(neG) = 2
!        jGvar(neG) = 2
         g(neG)     = 8.0d+0*x(2)

         neG        = neG + 1
!        iGfun(neG) = 3
!        jGvar(neG) = 1
         g(neG)     = 2.0d+0*(x(1) - 2.0d+0)

         neG        = neG + 1
!        iGfun(neG) = 3
!        jGvar(neG) = 2
         g(neG)     = 2.0d+0*x(2)
      end if

      end ! subroutine userfg
