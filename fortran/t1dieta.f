*     ------------------------------------------------------------------
*     File t1dieta.f
*     This is program illustrates the use of the interface snOptA,
*     which is part of the SNOPT 7 package.
*
*     04 Dec 2002: First version for SNOPT 6
*     04 Dec 2002: Current version.
*     ------------------------------------------------------------------
      program
     &     t1diet
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 10,
     &     maxn   = 10,
     &     lenA   = 30, lenG   = 1,
     &     nxname =  1, nFname = 1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character
     &     Prob*8, xnames(nxname)*8, Fnames(nFname)*8
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
*     ------------------------------------------------------------------
*     SNOPT workspace

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
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by hs106.

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lunit = iSpecs
         lfile = 't1dieta.spc'
         open( lunit, file=lfile, status='OLD',     err=800 )
         lunit = iPrint
         lfile = 't1dieta.out'
         open( lunit, file=lfile, status='UNKNOWN', err=800 )
      else

*        VMS  systems.  Define units for the Specs and print files.

         lunit = iSpecs
         open( lunit, status='OLD',     err=900 )
         lunit = iPrint
         open( lunit, status='UNKNOWN', err=900 )
      end if

*     ==================================================================
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ==================================================================
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     Set up the problem to be solved.
*     No derivatives are set in this case.

      Errors = 0

      call Diet0
     &   ( Errors, maxF, maxn, Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

*     snOptA will compute the Jacobian by finite-differences.
*     The user has the option of calling  snJac  to define the
*     coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
*
      call snJac
     &   ( INFO, nF, n, userf,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 102) go to 920

*     ------------------------------------------------------------------
*     Warn snOptA that userf does not set the derivatives.
*     The parameters iPrt and iSum may refer to the Print and Summary
*     file respectively.  Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option',    DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( ' Print   frequency 1', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( ' Summary frequency 1', DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
*     ------------------------------------------------------------------
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
      write(nOut, *) 't1dieta (1) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

*     Set up the data structure for the sparse linear constraints.

      call Diet1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

*     ------------------------------------------------------------------
*     Specify any options not set in the Specs file.
*     ------------------------------------------------------------------
      DerOpt = 1
      call snSeti
     &   ( 'Derivative option',    DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      Major  = 250
      call snSeti
     &   ( 'Major Iteration limit', Major, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Solve the problem again, this time with derivatives specified.
*     ------------------------------------------------------------------
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
      write(nOut, *) 't1dieta (2) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920
      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
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

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Diet0
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

*     ==================================================================
*     Diet0   defines input data for the Diet problem of Chvatal, 1983.
*
*
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 ) ( = objective row c')

*     Errors      is 0 if there is enough storage, 1 otherwise.
*     nF          is the number of problem functions
*                 (objective and constraints).
*     n           is the number of variables.
*     xlow        holds the lower bounds on x.
*     xupp        holds the upper bounds on x.
*     Flow        holds the lower bounds on F = Ax.
*     Fupp        holds the upper bounds on F = Ax.

*     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
*     x (1:n)     is a set of initial values for x.
*     Fmul(1:nF)  is a set of initial values for the dual variables.
*
*     ==================================================================
      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob = 'Diet LP.'

*     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      ObjRow = 4
      n      = 6

*     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      ObjAdd = 0.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on x.
*     ------------------------------------------------------------------
      do j = 1, n
         xlow(j) = 0.0d+0
      end do

      xupp(1) = 4.0d+0
      xupp(2) = 3.0d+0
      xupp(3) = 2.0d+0
      xupp(4) = 8.0d+0
      xupp(5) = 2.0d+0
      xupp(6) = 2.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on  Ax.
*     The objective row is free (i.e., infinite upper and lower bounds).
*     ------------------------------------------------------------------
      Flow( 1) =  2000.0d+0
      Flow( 2) =    55.0d+0
      Flow( 3) =   800.0d+0
      Flow( 4) = - plInfy

      do i = 1, nF
         Fupp(i) =  plInfy
         Fmul(i) =  0.0d+0
      end do

*     ----------------
*     Initialize  x.
*     ----------------
      do j = 1, n
         x(j) = 1.0d+0
      end do

      do j = 1, n
         xstate(j) = 0
      end do

      end ! subroutine Diet0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

*     ==================================================================
*     Computes the objective and constraint terms for the Diet problem
*     of Chvatal, 1983.
*
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 ) ( = objective row c')
*
*     ==================================================================
      integer
     &     Out
*     ------------------------------------------------------------------
      Out = 0
*     --------------------------------------------
*     Print something on the first and last entry.
*     -------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is problem  Diet0'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished problem Diet0'
         return
      end if

      F(1) = 110.0d+0*x(1) + 205.0d+0*x(2) + 160.0d+0*x(3)
     &     + 160.0d+0*x(4) + 420.0d+0*x(5) + 260.0d+0*x(6)
      F(2) =   4.0d+0*x(1) +  32.0d+0*x(2) +  13.0d+0*x(3)
     &     +   8.0d+0*x(4) +   4.0d+0*x(5) +  14.0d+0*x(6)
      F(3) =   2.0d+0*x(1) +  12.0d+0*x(2) +  54.0d+0*x(3)
     &     + 285.0d+0*x(4) +  22.0d+0*x(5) +  80.0d+0*x(6)
      F(4) =   3.0d+0*x(1) +  24.0d+0*x(2) +  13.0d+0*x(3)
     &     +   9.0d+0*x(4) +  20.0d+0*x(5) +  19.0d+0*x(6) ! Objective

      end ! subroutine userf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Diet1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character*8
     &     Prob, cw(lencw)

*     ==================================================================
*     Diet0   defines input data for the Diet problem of Chvatal, 1983.
*
*
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 ) ( = objective row c')

*     Errors      is 0 if there is enough storage, 1 otherwise.
*     nF          is the number of problem functions
*                 (objective and constraints).
*     n           is the number of variables.
*     xlow        holds the lower bounds on x.
*     xupp        holds the upper bounds on x.
*     Flow        holds the lower bounds on F = Ax.
*     Fupp        holds the upper bounds on F = Ax.

*     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
*     x (1:n)     is a set of initial values for x.
*     Fmul(1:nF)  is a set of initial values for the dual variables.
*
*     ==================================================================
      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision     plInfy
      parameter           (plInfy = 1.0d+20)
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob = 'Diet LP.'

*     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      ObjRow = 4
      n      = 6

*     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      ObjAdd = 0.0d+0

      neA    = 0

*     Column 1

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 1
      A(neA)     = 110.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 1
      A(neA)     = 4.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 1
      A(neA)     = 2.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 1
      A(neA)     = 3.0d+0

*     Column 2.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 2
      A(neA)     = 205.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 2
      A(neA)     = 32.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 2
      A(neA)     = 12.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 2
      A(neA)     = 24.0d+0

*     Column 3.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 3
      A(neA)     = 160.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 3
      A(neA)     = 13.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 3
      A(neA)     = 54.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 3
      A(neA)     = 13.0d+0

*     Column 4.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 4
      A(neA)     = 160.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 4
      A(neA)     = 8.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 4
      A(neA)     = 285.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 4
      A(neA)     = 9.0d+0

*     Column 5.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 5
      A(neA)     = 420.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 5
      A(neA)     = 4.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 5
      A(neA)     = 22.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 5
      A(neA)     = 20.0d+0

*     Column 6.

      neA        = neA + 1
      iAfun(neA) = 1
      jAvar(neA) = 6
      A(neA)     = 260.0d+0

      neA        = neA + 1
      iAfun(neA) = 2
      jAvar(neA) = 6
      A(neA)     = 14.0d+0

      neA        = neA + 1
      iAfun(neA) = 3
      jAvar(neA) = 6
      A(neA)     = 80.0d+0

      neA        = neA + 1
      iAfun(neA) = 4
      jAvar(neA) = 6
      A(neA)     = 19.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         xlow(j) = 0.0d+0
      end do

      xupp(1) = 4.0d+0
      xupp(2) = 3.0d+0
      xupp(3) = 2.0d+0
      xupp(4) = 8.0d+0
      xupp(5) = 2.0d+0
      xupp(6) = 2.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on  Ax.
*     The objective row is free (i.e., infinite upper and lower bounds).
*     ------------------------------------------------------------------
      Flow( 1) =  2000.0d+0
      Flow( 2) =    55.0d+0
      Flow( 3) =   800.0d+0
      Flow( 4) = - plInfy

      do i = 1, nF
         Fupp(i) =  plInfy
         Fmul(i) =  0.0d+0
      end do

*     ----------------
*     Initialize  x.
*     ----------------
      do j = 1, n
         x(j) = 1.0d+0
      end do

      do j = 1, n
         xstate(j) = 0
      end do

      end ! subroutine Diet1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

*     ==================================================================
*     Dummy  objective and constraint function for the Diet problem
*     of Chvatal, 1983.
*
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 ) ( = objective row c')
*
*     ==================================================================
*     Relax, A*x  is computed by snOptA from A.

      end ! subroutine userfg
