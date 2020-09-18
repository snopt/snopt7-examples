!     ------------------------------------------------------------------
!     File snmemtesta.f
!     This is a main program to illustrate how the routine snMemA
!     can be used to estimate the storage needed for snOptA.
!
!     27 Dec 2002: First   version.
!     11 Mar 2008: Current version.
!     ------------------------------------------------------------------
      program
     &     sntest

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
      character
     &     Prob*8, xnames(nxname)*8, Fnames(nFname)*8
      double precision
     &     ObjAdd, sInf, A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)
      integer
     &     Errors, neA, neG, ObjRow, INFO, i, iPrt, iPrint,
     &     iSum, iSumm, lencw, leniw, lenrw,
     &     Highcw, Highiw, Highrw,
     &     mincw, miniw, minrw, nF, n, nInf,
     &     nS, lenct, lenit, lenrt
      logical
     &     byname
      integer
     &     lunit
      character
     &     lfile*20
      external
     &     userfg

!     SNOPT workspace
      parameter          (  Highrw = 50000)
      double precision   rw(Highrw)
      parameter          (  Highiw = 50000)
      integer            iw(Highiw)
      parameter          (  Highcw =   500)
      character          cw(Highcw)*8

      parameter          (lenrt = 500)
      double precision   rt(lenrt)
      parameter          (lenit = 500)
      integer            it(lenit)
      parameter          (lenct = 500)
      character          ct(lenct)*8

      integer            Cold
      parameter         (Cold   = 0)

!     Specify some SNOPT files.

      iSumm  =  6               ! Summary file
      iPrint =  9               ! Print   file
      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the print file.

         lunit = iPrint
         lfile = 'snmemtesta.out'
         open( lunit, file=lfile, status='UNKNOWN', err=800 )
      else

!        VMS  systems.  Define a unit for the print file.

         lunit = iPrint
         open( lunit, status='UNKNOWN', err=900 )
      end if

!     ==================================================================
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!
!     Note that  ltmpcw, ltmpiw and ltmprw are installed as the default
!     upper limits on the SNOPT workspace.  These defaults will need
!     to be increased when better values have been estimated.
!     ==================================================================
      call snInit
     &   ( iPrint, iSumm, ct, lenct, it, lenit, rt, lenrt )

!     Set up the problem to be solved.

      Errors = 0

      call Hex1
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     ct, lenct, it, lenit, rt, lenrt )

!     ------------------------------------------------------------------
!     Compute an estimate of the storage needed by snOptA.
!     The storage estimation procedure is not strictly required in this
!     program because an overestimate of the required memory has been
!     defined above, as is required for f77.  Storage allocation is
!     included here as a guide for users wishing to allocate storage
!     using f90 or C.
!
!     The required values are mincw, miniw and minrw.
!     The default upper limits on the SNOPT workspace must be updated
!     with these values.
!     ------------------------------------------------------------------
      call snMemA
     &   ( INFO, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     ct, lenct, it, lenit, rt, lenrt )

      lencw = mincw
      leniw = miniw
      lenrw = minrw

!     allocate( cw(lencw), iw(leniw), rw(lenrw) )
!     Copy the first 500 elements of ct, it and rt into cw, iw and rw.

      do i = 1, 500
         cw(i) = ct(i)
         iw(i) = it(i)
         rw(i) = rt(i)
      end do

      iPrt   = 0                ! Suppress print   output
      iSum   = 0                ! Suppress summary output
      call snSeti
     &   ( 'Total character workspace', lencw, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Total integer   workspace', leniw, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Total real      workspace', lenrw, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start (Start = 0).
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

      write(iSumm, *) ' '
      write(iSumm, *) 'snmemtesta finished.'
      write(iSumm, *) 'snOptA INFO =', INFO
      write(iSumm, *) 'nInf        =', nInf
      write(iSumm, *) 'sInf        =', sInf
      if (ObjRow .gt. 0)
     &write(iSumm, *) 'Obj         =', F(ObjRow)
      if (INFO .ge. 30) go to 910
      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(iSumm, 4000) 'Error while opening file', lfile
      stop

  900 write(iSumm, 4010) 'Error while opening unit', lunit
      stop

  910 write(iSumm, *) ' '
      write(iSumm, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! program sntest

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Hex1
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
      character
     &     Prob*8, cw(lencw)*8

!     ==================================================================
!     Hex1 generates data for the Hexagon problem, first derivatives
!     only.
!
!     On entry,
!     maxF, maxn are upper limits on nF and n.
!
!     The  nF by n  Jacobian is written as the sum of
!     two  nF by n  sparse matrices G and A, i.e.,  J = A + G,  where
!     A  and  G  contain contain the constant and nonlinear
!     elements of J respectively.
!
!     The (fixed) coordinates of the kth nonzero of  G  are defined
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
!         order of the coordinates is meaningful in snOptA.)
!
!     On exit,
!     Errors      is 0 if there is enough storage, 1 otherwise.
!     nF          is the number of problem functions
!                 (objective and constraints, linear and nonlinear).
!     n           is the number of variables.
!     neG         is the number of nonzeros in Jn.
!     neA         is the number of nonzeros in Jc.
!     xlow        holds the lower bounds on x.
!     xupp        holds the upper bounds on x.
!     Flow        holds the lower bounds on F.
!     Fupp        holds the upper bounds on F.

!     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n)     is a set of initial values for x.
!     Fmul(1:nF)  is a set of initial values for the dual variables.
!
!     24 Dec 1997: First version of Hex1.
!     27 Oct 2002: Current version.
!     ==================================================================
      integer
     &     i, iPrt, iSum, j, Obj
!     ------------------------------------------------------------------
      double precision   InfBnd
      parameter         (InfBnd = 1.0d+20)
      double precision   zero,             one
      parameter         (zero  = 0.0d+0,   one    = 1.0d+0)
!     ------------------------------------------------------------------
!     Give the problem a name.

      Prob   = 'Hex1    '

      nF     = 19
      n      =  9
      Obj    = 19               ! Hexagon objective row
      ObjRow = 19               ! =0 if just a feasible point is needed

!     Check if there is enough storage.

      Errors = 0
      if (nF      .gt. maxF ) Errors = 1
      if (n       .gt. maxn ) Errors = 1
      if (Errors  .gt.   0  ) return

      iPrt   = 0
      iSum   = 0
      call snSet ( 'Maximize   ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     F( 1)  =    x_1^2 + x_6^2
!     F( 2)  =   (x_2   - x_1)^2  +  (x_7 - x_6)^2
!     F( 3)  =   (x_3   - x_1)^2  +   x_6^2
!     F( 4)  =   (x_1   - x_4)^2  +  (x_6 - x_8)^2
!     F( 5)  =   (x_1   - x_5)^2  +  (x_6 - x_9)^2
!     F( 6)  =    x_2^2 + x_7^2
!     F( 7)  =   (x_3   - x_2)^2  +   x_7^2
!     F( 8)  =   (x_4   - x_2)^2  +  (x_8 - x_7)^2
!     F( 9)  =   (x_2   - x_5)^2  +  (x_7 - x_9)^2
!     F(10)  =   (x_4   - x_3)^2  +   x_8^2
!     F(11)  =   (x_5   - x_3)^2  +   x_9^2
!     F(12)  =    x_4^2 +  x_8^2
!     F(13)  =   (x_4   - x_5)^2 + (x_9 - x_8)^2
!     F(14)  =    x_5^2 + x_9^2
!     F(15)  =  -x_1 + x_2
!     F(16)  =        -x_2 + x_3
!     F(17)  =               x_3 - x_4
!     F(18)  =                     x_4 - x_5
!     F(Obj) =  x_2 x_6 - x_1 x_7 + x_3 x_7 + x_5 x_8 - x_4 x_9 - x_3 x_8
!
!     The pattern of the Jacobian is as follows:
!
!               Column
!             | 1   2    3    4    5   6   7   8   9
!             +--------------------------------------
!     Row  1  | N                      N
!          2  | N   N                  N   N
!          3  | N        N             N
!          4  | N             N        N       N
!          5  | N                  N   N           N
!          6  |     N                      N
!          7  |     N    N                 N
!          8  |     N         N            N   N
!          9  |     N              N       N       N
!         10  |          N    N                N
!         11  |          N         N               N
!         12  |               N                N
!         13  |               N    N           N   N
!         14  |                    N               N
!         15  | L   L
!         16  |     L    L
!         17  |          L    L
!         18  |               L    L
!         19  | N   N    N    N    N   N   N   N   N  Objective row
!

!     --------------------------
!     Linear Jacobian  elements
!     --------------------------
      neA        =  0

      neA        = neA + 1
      iAfun(neA) = 15
      jAvar(neA) =  1
      A(neA)     = -one

      neA        = neA + 1
      iAfun(neA) = 15
      jAvar(neA) =  2
      A(neA)     =  one

      neA        = neA + 1
      iAfun(neA) = 16
      jAvar(neA) =  2
      A(neA)     = -one

      neA        = neA + 1
      iAfun(neA) = 16
      jAvar(neA) =  3
      A(neA)     =  one

      neA        = neA + 1
      iAfun(neA) = 17
      jAvar(neA) =  3
      A(neA)     =  one

      neA        = neA + 1
      iAfun(neA) = 17
      jAvar(neA) =  4
      A(neA)     = -one

      neA        = neA + 1
      iAfun(neA) = 18
      jAvar(neA) =  4
      A(neA)     =  one

      neA        = neA + 1
      iAfun(neA) = 18
      jAvar(neA) =  5
      A(neA)     = -one

!     neA        =  8

!     ------------------------------
!     Nonlinear Jacobian elements
!     ------------------------------
      neG        =  0

      neG        =  neG + 1
      iGfun(neG) =  1
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) =  1
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) =  2
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) =  2
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) =  2
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) =  2
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) =  3
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) =  3
      jGvar(neG) =  3

      neG        =  neG + 1
      iGfun(neG) =  3
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) =  4
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) =  4
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) =  4
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) =  4
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) =  5
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) =  5
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) =  5
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) =  5
      jGvar(neG) =  9

      neG        =  neG + 1
      iGfun(neG) =  6
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) =  6
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) =  7
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) =  7
      jGvar(neG) =  3

      neG        =  neG + 1
      iGfun(neG) =  7
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) =  8
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) =  8
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) =  8
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) =  8
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) =  9
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) =  9
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) =  9
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) =  9
      jGvar(neG) =  9

      neG        =  neG + 1
      iGfun(neG) = 10
      jGvar(neG) =  3

      neG        =  neG + 1
      iGfun(neG) = 10
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) = 10
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) = 11
      jGvar(neG) =  3

      neG        =  neG + 1
      iGfun(neG) = 11
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) = 11
      jGvar(neG) =  9

      neG        =  neG + 1
      iGfun(neG) = 12
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) = 12
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) = 13
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) = 13
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) = 13
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) = 13
      jGvar(neG) =  9

      neG        =  neG + 1
      iGfun(neG) = 14
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) = 14
      jGvar(neG) =  9

!     --------------------------
!     Objective gradient
!     --------------------------
      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  1

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  2

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  3

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  4

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  5

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  6

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  7

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  8

      neG        =  neG + 1
      iGfun(neG) = Obj
      jGvar(neG) =  9

!     Final neG is 53

!     ------------------------------------------------------------------
!     Ranges for the problem functions.
!     ------------------------------------------------------------------
!     Nonlinear constraints first.
!     Followed by the linear constraints.
!     The Objective row is free.

      do i = 1, 14
         Flow(i) = -InfBnd
         Fupp(i) =  one
      end do


      do i = 15, 18
         Flow(i) =  zero
         Fupp(i) =  InfBnd
      end do

      Flow(Obj) = -InfBnd
      Fupp(Obj) =  InfBnd

      ObjAdd  = zero

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, n
         xlow(j) = -InfBnd
         xupp(j) =  InfBnd
      end do

      xlow(1) =  zero
      xlow(3) = -one
      xlow(5) =  zero
      xlow(6) =  zero
      xlow(7) =  zero

      xupp(3) =  one
      xupp(8) =  zero
      xupp(9) =  zero

!     ------------------------------------------------------------------
!     Initialize x, xstate and Fmul.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  .1d+0
      x(2)   =  .125d+0
      x(3)   =  .666666d+0
      x(4)   =  .142857d+0
      x(5)   =  .111111d+0
      x(6)   =  .2d+0
      x(7)   =  .25d+0
      x(8)   = -.2d+0
      x(9)   = -.25d+0

      do j = 1, n
         xstate(j) = 0
      end do

      do i = 1, nF
         Fmul(i) = zero
!        Fmul(i) = one
      end do

      end ! subroutine fgdata

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg
     &   ( Status,    n, x,
     &      needF,   nF, F,
     &      needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Computes the nonlinear objective and constraint terms for the
!     problem Hexagon.
!
!     No user-defined storage is needed.
!     ==================================================================
      integer
     &     neG, Obj
!     ------------------------------------------------------------------
      double precision   two
      parameter         (two = 2.0d+0)
      integer            Out
      parameter         (Out  = 6)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is  Hex1 (Maximized)'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished Hex1'
         return
      end if

      Obj = 19

      if (needF .gt. 0) then
         F( 1) =    x(1)**2          +  x(6)**2
         F( 2) =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         F( 3) =   (x(3) - x(1))**2  +  x(6)**2
         F( 4) =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         F( 5) =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         F( 6) =    x(2)**2          +  x(7)**2
         F( 7) =   (x(3) - x(2))**2  +  x(7)**2
         F( 8) =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         F( 9) =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         F(10) =   (x(4) - x(3))**2  +  x(8)**2
         F(11) =   (x(5) - x(3))**2  +  x(9)**2
         F(12) =    x(4)**2          +  x(8)**2
         F(13) =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         F(14) =    x(5)**2          +  x(9)**2

!        Objective (maximized).

         F(Obj) =   x(2)*x(6)        - x(1)*x(7) + x(3)*x(7) + x(5)*x(8)
     &                               - x(4)*x(9) - x(3)*x(8)
      end if

      if (needG .gt. 0) then

         neG        = 0

!        Constraint gradients

         neG        =  neG + 1
         G(neG)     =  two*x(1)
!        iGfun(neG) =  1
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  two*x(6)
!        iGfun(neG) =  1
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     = -two*(x(2) - x(1))
!        iGfun(neG) =  2
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =   two*(x(2) - x(1))
!        iGfun(neG) =  2
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  -two*(x(7) - x(6))
!        iGfun(neG) =  2
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     =   two*(x(7) - x(6))
!        iGfun(neG) =  2
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =  -two*(x(3) - x(1))
!        iGfun(neG) =  3
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =   two*(x(3) - x(1))
!        iGfun(neG) =  3
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =   two*x(6)
!        iGfun(neG) =  3
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     =   two*(x(1) - x(4))
!        iGfun(neG) =  4
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  -two*(x(1) - x(4))
!        iGfun(neG) =  4
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =   two*(x(6) - x(8))
!        iGfun(neG) =  4
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     =  -two*(x(6) - x(8))
!        iGfun(neG) =  4
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     =   two*(x(1) - x(5))
!        iGfun(neG) =  5
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  -two*(x(1) - x(5))
!        iGfun(neG) =  5
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =   two*(x(6) - x(9))
!        iGfun(neG) =  5
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     =  -two*(x(6) - x(9))
!        iGfun(neG) =  5
!        jGvar(neG) =  9

         neG        =  neG + 1
         G(neG)     =   two*x(2)
!        iGfun(neG) =  6
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =   two*x(7)
!        iGfun(neG) =  6
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =  -two*(x(3) - x(2))
!        iGfun(neG) =  7
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =   two*(x(3) - x(2))
!        iGfun(neG) =  7
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =   two*x(7)
!        iGfun(neG) =  7
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =  -two*(x(4) - x(2))
!        iGfun(neG) =  8
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =   two*(x(4) - x(2))
!        iGfun(neG) =  8
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  -two*(x(8) - x(7))
!        iGfun(neG) =  8
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =   two*(x(8) - x(7))
!        iGfun(neG) =  8
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     =   two*(x(2) - x(5))
!        iGfun(neG) =  9
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =  -two*(x(2) - x(5))
!        iGfun(neG) =  9
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =   two*(x(7) - x(9))
!        iGfun(neG) =  9
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =  -two*(x(7) - x(9))
!        iGfun(neG) =  9
!        jGvar(neG) =  9

         neG        =  neG + 1
         G(neG)     =  -two*(x(4) - x(3))
!        iGfun(neG) = 10
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =   two*(x(4) - x(3))
!        iGfun(neG) = 10
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =   two*x(8)
!        iGfun(neG) = 10
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     =  -two*(x(5) - x(3))
!        iGfun(neG) = 11
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =   two*(x(5) - x(3))
!        iGfun(neG) = 11
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =   two*x(9)
!        iGfun(neG) = 11
!        jGvar(neG) =  9

         neG        =  neG + 1
         G(neG)     =   two*x(4)
!        iGfun(neG) = 12
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =   two*x(8)
!        iGfun(neG) = 12
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     =   two*(x(4) - x(5))
!        iGfun(neG) = 13
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  -two*(x(4) - x(5))
!        iGfun(neG) = 13
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =  -two*(x(9) - x(8))
!        iGfun(neG) = 13
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     =   two*(x(9) - x(8))
!        iGfun(neG) = 13
!        jGvar(neG) =  9

         neG        =  neG + 1
         G(neG)     =   two*x(5)
!        iGfun(neG) = 14
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =   two*x(9)
!        iGfun(neG) = 14
!        jGvar(neG) =  9

!        --------------------------
!        Objective gradient
!        --------------------------
         neG        =  neG + 1
         G(neG)     = - x(7)
!        iGfun(neG) = Obj
!        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =   x(6)
!        iGfun(neG) = Obj
!        jGvar(neG) =  2

         neG        =  neG + 1
         G(neG)     =   x(7) - x(8)
!        iGfun(neG) = Obj
!        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     = - x(9)
!        iGfun(neG) = Obj
!        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =   x(8)
!        iGfun(neG) = Obj
!        jGvar(neG) =  5

         neG        =  neG + 1
         G(neG)     =   x(2)
!        iGfun(neG) = Obj
!        jGvar(neG) =  6

         neG        =  neG + 1
         G(neG)     =   x(3) - x(1)
!        iGfun(neG) = Obj
!        jGvar(neG) =  7

         neG        =  neG + 1
         G(neG)     =   x(5) - x(3)
!        iGfun(neG) = Obj
!        jGvar(neG) =  8

         neG        =  neG + 1
         G(neG)     = - x(4)
!        iGfun(neG) = Obj
!        jGvar(neG) =  9
      end if

      end ! subroutine userfg
