!     ------------------------------------------------------------------
!     File snmemtestb.f
!     This is a main program to illustrate how the routine snMemB
!     can be used to estimate the storage needed for snOptB.
!
!     27 Dec 2002: First   version.
!     11 Mar 2008: Current version.
!     ------------------------------------------------------------------
      program
     &     snmemtestb
      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxne  = 3000,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)
      integer
     &     lenrw, leniw, lencw,
     &     Highcw, Highiw, Highrw,
     &     lenct, lenit, lenrt
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  Highrw = 50000)
      double precision   rw(Highrw)
      parameter          (  Highiw = 50000)
      integer            iw(Highiw)
      parameter          (  Highcw =  1000)
      character          cw(Highcw)*8

      parameter          (lenrt = 500)
      double precision   rt(lenrt)
      parameter          (lenit = 500)
      integer            it(lenit)
      parameter          (lenct = 500)
      character          ct(lenct)*8
!     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      external
     &     HexCon, HexObj
      integer
     &     Errors, INFO, i, iObj, iPrt, iPrint, iSum, iSumm, lunit, m,
     &     mincw, miniw, minrw, n, ne, neGcon, nInf, nnCon, nnJac,
     &     nnObj, nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
!     ------------------------------------------------------------------
!     Specify some snOptB files.
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).

      iSumm  =  6               ! Summary file
      iPrint =  9               ! Print   file
      nOut   =  6
      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the print file.

         lunit = iPrint
         lfile = 'snmemtestb.out'
         open( lunit, file=lfile, status='UNKNOWN', err=800 )
      else

!        VMS  systems.  Define a unit for the print file.

         lunit = iPrint
         open( lunit, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     The storage estimation procedure is not strictly required in this
!     program because an overestimate of the required memory has been
!     defined above, as is required for f77.  Storage allocation is
!     included here as a guide for users wishing to allocate storage
!     using f90 or C.
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, ct, lenct, it, lenit, rt, lenrt )

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      Errors = 0

      call HexDat
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi,
     &     ct, lenct, it, lenit, rt, lenrt )

!     ------------------------------------------------------------------
!     Compute an estimate of the storage needed by snOptB.
!     Copy the first 500 elements of ct, it and rt into cw, iw and rw.
!     The required values are mincw, miniw and minrw.
!     The default upper limits on the SNOPT workspace must be updated
!     with these values.
!     ------------------------------------------------------------------
!     Everything is known except neGcon, the number of nonlinear
!     Jacobian elements.  We give an upper bound  nnCon*nnJac  based
!     on the assumption that the nonlinear Jacobian is dense.

      neGcon = nnCon*nnJac

      call snMemB
     &   ( INFO, m, n, ne, neGcon,
     &     nnCon, nnJac, nnObj,
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
!     Go for it, using a Cold start.
!     ------------------------------------------------------------------
      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     HexCon, HexObj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'snmemtestb finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910

      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program snmemtestb

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObj
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem Hexagon.
!     No user-defined storage is used.
!     ==================================================================

      fObj    =   x(2)*x(6) - x(1)*x(7) + x(3)*x(7) + x(5)*x(8)
     &          - x(4)*x(9) - x(3)*x(8)

      gObj(1) = - x(7)
      gObj(2) =   x(6)
      gObj(3) =   x(7) - x(8)
      gObj(4) = - x(9)
      gObj(5) =   x(8)
      gObj(6) =   x(2)
      gObj(7) =   x(3) - x(1)
      gObj(8) =   x(5) - x(3)
      gObj(9) = - x(4)

      end ! subroutine HexObj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexCon
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem Hexagon.
!
!     No user-defined storage is used.
!     ==================================================================
      double precision   two
      parameter         (two   = 2.0d+0)
!     ------------------------------------------------------------------
      fCon( 1) =    x(1)**2          +  x(6)**2
      fCon( 2) =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
      fCon( 3) =   (x(3) - x(1))**2  +  x(6)**2
      fCon( 4) =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
      fCon( 5) =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
      fCon( 6) =    x(2)**2          +  x(7)**2
      fCon( 7) =   (x(3) - x(2))**2  +  x(7)**2
      fCon( 8) =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
      fCon( 9) =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
      fCon(10) =   (x(4) - x(3))**2  +  x(8)**2
      fCon(11) =   (x(5) - x(3))**2  +  x(9)**2
      fCon(12) =    x(4)**2          +  x(8)**2
      fCon(13) =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
      fCon(14) =    x(5)**2          +  x(9)**2

!     Nonlinear Jacobian elements for column 1.
!     rows = (1,2,3,4,5).

      gCon( 1) =   two*x(1)
      gCon( 2) = - two*(x(2) - x(1))
      gCon( 3) = - two*(x(3) - x(1))
      gCon( 4) =   two*(x(1) - x(4))
      gCon( 5) =   two*(x(1) - x(5))

!     Nonlinear Jacobian elements for column 2.
!     Rows = (2,6,7,8,9).

      gCon( 6) =   two*(x(2) - x(1))
      gCon( 7) =   two*x(2)
      gCon( 8) = - two*(x(3) - x(2))
      gCon( 9) = - two*(x(4) - x(2))
      gCon(10) =   two*(x(2) - x(5))

!     Nonlinear Jacobian elements for column 3.
!     Rows = (3,7,10,11).

      gCon(11) =   two*(x(3) - x(1))
      gCon(12) =   two*(x(3) - x(2))
      gCon(13) = - two*(x(4) - x(3))
      gCon(14) = - two*(x(5) - x(3))

!     Nonlinear Jacobian elements for column 4.
!     Rows = (4,8,10,12,13).

      gCon(15) = - two*(x(1) - x(4))
      gCon(16) =   two*(x(4) - x(2))
      gCon(17) =   two*(x(4) - x(3))
      gCon(18) =   two*x(4)
      gCon(19) =   two*(x(4) - x(5))

!     Nonlinear Jacobian elements for column 5.
!     Rows = (5,9,11,13,14).

      gCon(20) = - two*(x(1) - x(5))
      gCon(21) = - two*(x(2) - x(5))
      gCon(22) =   two*(x(5) - x(3))
      gCon(23) = - two*(x(4) - x(5))
      gCon(24) =   two*x(5)

!     Nonlinear Jacobian elements for column 6.
!     Rows = (1,2,3,4,5).

      gCon(25) =   two*x(6)
      gCon(26) = - two*(x(7) - x(6))
      gCon(27) =   two*x(6)
      gCon(28) =   two*(x(6) - x(8))
      gCon(29) =   two*(x(6) - x(9))

!     Nonlinear Jacobian elements for column 7.
!     Rows = (2,6,7,8,9).

      gCon(30) =   two*(x(7) - x(6))
      gCon(31) =   two*x(7)
      gCon(32) =   two*x(7)
      gCon(33) = - two*(x(8) - x(7))
      gCon(34) =   two*(x(7) - x(9))

!     Nonlinear Jacobian elements for column 8.
!     Rows = (4,8,10,12,13).

      gCon(35) = - two*(x(6) - x(8))
      gCon(36) =   two*(x(8) - x(7))
      gCon(37) =   two*x(8)
      gCon(38) =   two*x(8)
      gCon(39) = - two*(x(9) - x(8))

!     Nonlinear Jacobian elements for column 9.
!     Rows = (5,9,11,13,14).

      gCon(40) = - two*(x(6) - x(9))
      gCon(41) = - two*(x(7) - x(9))
      gCon(42) =   two*x(9)
      gCon(43) =   two*(x(9) - x(8))
      gCon(44) =   two*x(9)

      end ! subroutine HexCon

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexDat
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencw, leniw, lenrw, maxm, maxn, maxne, Errors, m, n,
     &     ne, nnCon, nnObj, nnJac, iObj, indA(maxne), hs(maxn+maxm),
     &     locA(maxn+1), iw(leniw)
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

!     ------------------------------------------------------------------
!     HexDat generates data for the Hexagon problem.
!     The constraints take the form
!              c(x) + A*x - s = 0,
!     where the Jacobian for c(x) + Ax is stored in Acol(*), and any
!     terms coming from c(x) are in the TOP LEFT-HAND CORNER of Acol(*),
!     with dimensions  nnCon x nnJac.
!     Note that the right-hand side is zero.
!     s is a set of slack variables whose bounds contain any constants
!     that might have formed a right-hand side.
!
!     The objective function is
!             f(x) + d'x
!     where d would be row iobj of A (but there is no such row in
!     this example).  f(x) involves only the FIRST nnObj variables.
!
!     On entry,
!     maxm, maxn, maxne are upper limits on m, n, ne.
!
!     On exit,
!     Errors  is 0 if there is enough storage, 1 otherwise.
!     m       is the number of nonlinear and linear constraints.
!     n       is the number of variables.
!     ne      is the number of nonzeros in Acol(*).
!     nnCon   is the number of nonlinear constraints (they come first).
!     nnObj   is the number of nonlinear objective variables.
!     nnJac   is the number of nonlinear Jacobian variables.
!     Acol    is the constraint matrix (Jacobian), stored column-wise.
!     indA    is the list of row indices for each nonzero in Acol(*).
!     locA    is a set of pointers to the start of each column of Acol.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n) is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!
!     24 Dec 1997: First version of HexDat.
!     ------------------------------------------------------------------
      integer
     &     i, iPrt, iSum, j
!     ------------------------------------------------------------------
      double precision   bplus
      parameter         (bplus   = 1.0d+20)
      double precision   zero,               one
      parameter         (zero    = 0.0d+0,   one    = 1.0d+0)
!     ------------------------------------------------------------------

!     Give a name to the Problem.

      Prob   = 'Hexagon '

      ne     = 52
      n      =  9
      m      = 18

      nnCon  = 14
      nnJac  =  n
      nnObj  =  n

!     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

      iPrt   = 0
      iSum   = 0
      call snSet
     &   ( 'Maximize   ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     --------------------------------------
!     Set up the list of row indices in indA.
!     --------------------------------------
!     Column  1
!     Nonlinear elements in rows (1, 2, 3, 4, 5)  first.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) =  2
      indA( 3) =  3
      indA( 4) =  4
      indA( 5) =  5

      Acol( 1) =  zero
      Acol( 2) =  zero
      Acol( 3) =  zero
      Acol( 4) =  zero
      Acol( 5) =  zero

!     Column 1.
!     Linear element in row 6 next.

      indA( 6) = 15

      Acol( 6)  = -one

!     Column 2.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).

      locA( 2) =  7

      indA( 7) =  2
      indA( 8) =  6
      indA( 9) =  7
      indA(10) =  8
      indA(11) =  9

      Acol( 7) =  zero
      Acol( 8) =  zero
      Acol( 9) =  zero
      Acol(10) =  zero
      Acol(11) =  zero

!     Column 2.
!     Linear elements in rows (15,16).

      indA(12) = 15
      indA(13) = 16

      Acol(12) =  one
      Acol(13) = -one

!     Column 3.
!     Nonlinear elements in rows (3, 7, 10, 11).

      locA( 3) =  14

      indA(14) =  3
      indA(15) =  7
      indA(16) = 10
      indA(17) = 11

      Acol(14) =  zero
      Acol(15) =  zero
      Acol(16) =  zero
      Acol(17) =  zero

!     Column 3.
!     Linear elements in rows (16, 17).

      indA(18) = 16
      indA(19) = 17

      Acol(18) =  one
      Acol(19) =  one

!     Column 4.
!     Nonlinear elements in rows (20, 21, 22, 23, 24).

      locA( 4) = 20

      indA(20) =  4
      indA(21) =  8
      indA(22) = 10
      indA(23) = 12
      indA(24) = 13

      Acol(20) =  zero
      Acol(21) =  zero
      Acol(22) =  zero
      Acol(23) =  zero
      Acol(24) =  zero

!     Column 4.
!     Linear elements in rows (17, 18).

      indA(25) = 17
      indA(26) = 18

      Acol(25) = -one
      Acol(26) =  one

!     Column 5.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).

      locA( 5) = 27

      indA(27) =  5
      indA(28) =  9
      indA(29) = 11
      indA(30) = 13
      indA(31) = 14

      Acol(27) =  zero
      Acol(28) =  zero
      Acol(29) =  zero
      Acol(30) =  zero
      Acol(31) =  zero

!     Column 5.
!     Linear element in row 18.

      indA(32) = 18

      Acol(32) = -one

!     Column 6.
!     Nonlinear elements in rows (1, 2, 3, 4, 5, 6).

      locA(6)  = 33

      indA(33) =  1
      indA(34) =  2
      indA(35) =  3
      indA(36) =  4
      indA(37) =  5

      Acol(33) =  zero
      Acol(34) =  zero
      Acol(35) =  zero
      Acol(36) =  zero
      Acol(37) =  zero

!     Column 7.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).

      locA(7)  =  38

      indA(38) =  2
      indA(39) =  6
      indA(40) =  7
      indA(41) =  8
      indA(42) =  9

      Acol(38) =  zero
      Acol(39) =  zero
      Acol(40) =  zero
      Acol(41) =  zero
      Acol(42) =  zero

!     Column 8.
!     Nonlinear elements in rows (4, 8, 10, 12, 13).

      locA(8)  =  43

      indA(43) =  4
      indA(44) =  8
      indA(45) = 10
      indA(46) = 12
      indA(47) = 13

      Acol(43) =  zero
      Acol(44) =  zero
      Acol(45) =  zero
      Acol(46) =  zero
      Acol(47) =  zero

!     Column 9.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).

      locA(9)  =  48

      indA(48) =  5
      indA(49) =  9
      indA(50) = 11
      indA(51) = 13
      indA(52) = 14

      Acol(48) =  zero
      Acol(49) =  zero
      Acol(50) =  zero
      Acol(51) =  zero
      Acol(52) =  zero

!     Don't forget to finish off  locA.
!     This is crucial.

      locA(10) =  ne + 1

!     ------------------------------------------------------------------
!     Constraint ranges
!     ------------------------------------------------------------------
!     Nonlinear constraints first.

      do i = 1, nnCon
         bl(i+n) = -bplus
         bu(i+n) =  one
      end do

!     Followed by the linear constraints.

      do i = nnCon+1, m
         bl(i+n) =  zero
         bu(i+n) =  bplus
      end do

!     No linear objective term for this problem.

      iObj    = 0
      ObjAdd  = zero

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -bplus
         bu(j) =  bplus
      end do

      bl(1) =  zero
      bl(3) = -one
      bl(5) =  zero
      bl(6) =  zero
      bl(7) =  zero

      bu(3) =  one
      bu(8) =  zero
      bu(9) =  zero

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
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
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = zero
      end do

      end ! subroutine HexDat
