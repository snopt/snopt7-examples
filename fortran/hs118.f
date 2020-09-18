!     ------------------------------------------------------------------
!     File hs118.f
!     This is an example of solving a QP using a call to SQOPT, which is
!     part of the SNOPT package.
!
!     04 Oct 1994: First   version.
!     09 Jul 2005: Current version.
!     ------------------------------------------------------------------
      program
     &     hsmain
      implicit
     &     none
      integer
     &     maxm, maxn, maxne
      parameter
     &   ( maxm   = 10000,
     &     maxn   = 15000,
     &     maxne  = 30000 )

      character
     &     Prob*8, Names(maxn+maxm)*8
      integer
     &     indA(maxne) , eType(maxn+maxm), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm), pi(maxm), rc(maxn+maxm)

!     SQOPT workspace---------------------------------------------------
      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lencw = 500, leniw = 5000, lenrw = 10000 )
      character
     &     cw(lencw)*8
      integer
     &     iw(leniw)
      double precision
     &     rw(lenrw)
!     ------------------------------------------------------------------
      logical
     &     byname
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSumm, iSpecs, j,
     &     lencObj, m, maxS, mincw, miniw, minrw, n, ncolH, ne, nInf,
     &     nName, nOut, nS
      double precision
     &     infBnd, Obj, ObjAdd, sInf
      character
     &     lfile*20
      external
     &     myHx

!     ------------------------------------------------------------------
!     Specify some of the SQOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!
!     nOut    is an output file used here by main.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'hs118.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'hs118.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     First,  sqInit MUST be called to initialize optional parameters
!     to their default values.
!     ------------------------------------------------------------------
      call sqInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file (Optional).
!     ------------------------------------------------------------------
      call sqSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

!     ------------------------------------------------------------------
!     Define the problem.
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants ObjAdd and iObj so that the objective to
!         be minimized is     ObjAdd + x(n+iObj) + half x'*H*x
!     ------------------------------------------------------------------
      call hs118
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, x )

!     ------------------------------------------------------------------
!     Fix the column variables to be non-elastic and the row  variables
!     to be elastic.
!     ------------------------------------------------------------------
      do j = 1, n
         eType(j) = 0
      end do

      do j = n+1, n+m
         eType(j) = 3
      end do

!     ------------------------------------------------------------------
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      do j = 1, n
         hs(j) = 0
      end do

!     ------------------------------------------------------------------
!     Specify options that were not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors = 0

      maxS   = 10
      i1     =  0
      i2     =  0
      call sqseti
     &   ( 'Superbasics Limit ', maxS  , i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset   ! Required for following Warm start call
     &   ( 'Sticky parameters yes', i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
!     hs     need not be set if a basis file is to be input.
!            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
!            The values are used by the Crash procedure s2crsh
!            to choose an initial basis B.
!            If hs(j) = 0 or 1, column j is eligible for B.
!            If hs(j) = 2, column j is initially superbasic (not in B).
!            If hs(j) = 3, column j is eligible for B and is given
!                          preference over columns with hs(j) = 0 or 1.
!            If hs(j) = 4 or 5, column j is initially nonbasic.
!     ------------------------------------------------------------------
      call sqopt
     &   ( 'Cold', myHx, m,
     &     n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 990
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs118 (Cold start) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     Alter some options and call sqopt again, testing the Warm start.
!     The following illustrates the use of sqset, sqseti and sqsetr
!     to set specific options.  We can ensure that all unspecified
!     options take default values by first calling
!     sqset ( 'Defaults', ... ).
!     Beware that certain parameters would then need to be redefined.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) 'Alter options and test Warm start:'

      Errors  =  0
      infBnd  =  1.0d+21
      call sqset
     &   ( '                  ',         iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Defaults          ',         iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Scale option      ',      0, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print frequency   ',      1, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Summary frequency ',      1, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqsetr
     &   ( 'Infinite Bound    ', infBnd, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (Errors .gt. 0) then
         write(nOut, *) 'NOTE: Some of the options were not recognized'
      end if

!     Test the Warm start.
!     hs(*) specifies a complete basis from the previous call.
!     A Warm start uses hs(*) directly, without calling Crash.
!
!     Warm  starts are normally used after sqopt has solved a
!     problem with the SAME DIMENSIONS but perhaps altered data.
!     Here we have not altered the data, so very few iterations
!     should be required.

      call sqopt
     &   ( 'Warm', myHx, m,
     &     n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '
      write(nOut, *) 'hs118 (warm start) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

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

      end ! program hsmain

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine myHx
     &   ( ncolH, x, Hx, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     ncolH, nState, lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(ncolH), Hx(ncolH), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     This is myHx   for problem Hock-Schittkowski 118.
!     ==================================================================
      integer
     &     i, iPrint
!     ------------------------------------------------------------------
      iPrint = 9

!     ---------------------------------------
!     First entry.  Print something.
!     ---------------------------------------
      if (nState .eq. 1) then
         if (iPrint .gt. 0) write(iPrint, 1000) ncolH
      end if

!     -------------
!     Normal entry.
!     -------------
      do i = 1, 5
         Hx(3*i-2) =  2.0d-4 * x(3*i-2)
         Hx(3*i-1) =  2.0d-4 * x(3*i-1)
         Hx(3*i)   =  3.0d-4 * x(3*i)
      end do

!     ------------
!     Last entry.
!     ------------
      if (nState .ge. 2) then
         if (iPrint .gt. 0) write(iPrint, 2000)
      end if

      return

 1000 format(// ' This is problem  hs118.   ncolH =', i4)
 2000 format(// ' Finished         hs118.')

      end ! subroutine myHx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs118
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nName, lencObj, ncolH, iObj,
     &     indA(maxne), locA(maxn+1)
      character
     &     Prob*8, Names(maxn+maxm)*8
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

!     ==================================================================
!     hs118   sets the constraints and bounds for the
!     quadratic program (Problem 118 of Hock and Schittkowski).
!     Note that the linear objective term is the last row of A.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     infBnd
!     ------------------------------------------------------------------
!     Give a name to the problem.

      Prob   = 'HS 118..'

      n      = 15
      m      = 18      ! Includes the objective row
      ne     = 54

      lencObj = 0
      ncolH   = 15
      nName   = 1

      infBnd =  1.0d+20

!     ------------------------------------------------------------------
!     Assign the constraint nonzeros to Acol, column by column.
!     indA(i) gives the row index of element Acol(i).
!     locA(j) gives the index in a of the start of column j.
!     ------------------------------------------------------------------
!     Column 1.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) = 13
      indA( 3) = 18

      Acol( 1) = -1.0d+0
      Acol( 2) =  1.0d+0
      Acol( 3) =  2.3d+0

!     Column 2.

      locA( 2) =  4

      indA( 4) =  5
      indA( 5) = 13
      indA( 6) = 18

      Acol( 4) = -1.0d+0
      Acol( 5) =  1.0d+0
      Acol( 6) =  1.7d+0

!     Column 3.

      locA( 3) =  7

      indA( 7) =  9
      indA( 8) = 13
      indA( 9) = 18

      Acol( 7) = -1.0d+0
      Acol( 8) =  1.0d+0
      Acol( 9) =  2.2D+0

!     Column 4.

      locA( 4) = 10

      indA(10) =  1
      indA(11) =  2
      indA(12) = 14
      indA(13) = 18

      Acol(10) =  1.0d+0
      Acol(11) = -1.0d+0
      Acol(12) =  1.0d+0
      Acol(13) =  2.3d+0

!     Column 5.

      locA( 5) = 14

      indA(14) =  5
      indA(15) =  6
      indA(16) = 14
      indA(17) = 18

      Acol(14)  =  1.0d+0
      Acol(15)  = -1.0d+0
      Acol(16)  =  1.0d+0
      Acol(17)  =  1.7d+0

!     Column 6.

      locA(6)  = 18

      indA(18) =  9
      indA(19) = 10
      indA(20) = 14
      indA(21) = 18

      Acol(18) =  1.0d+0
      Acol(19) = -1.0d+0
      Acol(20) =  1.0d+0
      Acol(21) =  2.2D+0

!     Column 7.

      locA(7)  = 22

      indA(22) =  2
      indA(23) =  3
      indA(24) = 15
      indA(25) = 18

      Acol(22) =  1.0d+0
      Acol(23) = -1.0d+0
      Acol(24) =  1.0d+0
      Acol(25) =  2.3d+0

!     Column 8.

      locA(8)  = 26

      indA(26) =  6
      indA(27) =  7
      indA(28) = 15
      indA(29) = 18

      Acol(26) =  1.0d+0
      Acol(27) = -1.0d+0
      Acol(28) =  1.0d+0
      Acol(29) =  1.7d+0

!     Column 9.

      locA(9)  = 30

      indA(30) = 10
      indA(31) = 11
      indA(32) = 15
      indA(33) = 18

      Acol(30) =  1.0d+0
      Acol(31) = -1.0d+0
      Acol(32) =  1.0d+0
      Acol(33) =  2.2D+0

!     Column 10.

      locA(10) = 34

      indA(34) =  3
      indA(35) =  4
      indA(36) = 16
      indA(37) = 18

      Acol(34) =  1.0d+0
      Acol(35) = -1.0d+0
      Acol(36) =  1.0d+0
      Acol(37) =  2.3d+0

!     Column 11.

      locA(11) = 38

      indA(38) =  7
      indA(39) =  8
      indA(40) = 16
      indA(41) = 18

      Acol(38) =  1.0d+0
      Acol(39) = -1.0d+0
      Acol(40) =  1.0d+0
      Acol(41) =  1.7d+0

!     Column 12.

      locA(12) = 42

      indA(42) = 11
      indA(43) = 12
      indA(44) = 16
      indA(45) = 18

      Acol(42) =  1.0d+0
      Acol(43) = -1.0d+0
      Acol(44) =  1.0d+0
      Acol(45) =  2.2D+0

!     Column 13.

      locA(13) = 46

      indA(46) =  4
      indA(47) = 17
      indA(48) = 18

      Acol(46) =  1.0d+0
      Acol(47) =  1.0d+0
      Acol(48) =  2.3d+0

!     Column 14.

      locA(14) = 49

      indA(49) =  8
      indA(50) = 17
      indA(51) = 18

      Acol(49) =  1.0d+0
      Acol(50) =  1.0d+0
      Acol(51) =  1.7d+0

!     Column 15.

      locA(15) = 52

      indA(52) = 12
      indA(53) = 17
      indA(54) = 18

      Acol(52) =  1.0d+0
      Acol(53) =  1.0d+0
      Acol(54) =  2.2D+0

!     locA(n+1)-1 points to the last nonzero of the nth column.

      locA(16) = 55

!     ------------------------------------------------------------------
!     Define l and u such that l <=  Ax  <= u.
!     Temporarily store  l  and  u  in  bl(n+1:n+m)  and  bu(n+1:n+m).
!     Set the default l and u.
!     ------------------------------------------------------------------
      do i = 1, m
         j     = n + i
         bl(j) = 0.0d+0
         bu(j) = infBnd
      end do

!     iObj   = 18  means the linear objective is row 18 in Acol(*).
!     The objective row is free.

      iObj   = 18
      bl(n+iObj) = -infBnd

      bl(n+1)  =  -7.0d+0
      bu(n+1)  =   6.0d+0

      bl(n+2)  =  -7.0d+0
      bu(n+2)  =   6.0d+0

      bl(n+3)  =  -7.0d+0
      bu(n+3)  =   6.0d+0

      bl(n+4)  =  -7.0d+0
      bu(n+4)  =   6.0d+0

      bl(n+5)  =  -7.0d+0
      bu(n+5)  =   7.0d+0

      bl(n+6)  =  -7.0d+0
      bu(n+6)  =   7.0d+0

      bl(n+7)  =  -7.0d+0
      bu(n+7)  =   7.0d+0

      bl(n+8)  =  -7.0d+0
      bu(n+8)  =   7.0d+0

      bl(n+9)  =  -7.0d+0
      bu(n+9)  =   6.0d+0

      bl(n+10) =  -7.0d+0
      bu(n+10) =   6.0d+0

      bl(n+11) =  -7.0d+0
      bu(n+11) =   6.0d+0

      bl(n+12) =  -7.0d+0
      bu(n+12) =   6.0d+0

      bl(n+13) =  60.0d+0
      bl(n+14) =  50.0d+0
      bl(n+15) =  70.0d+0
      bl(n+16) =  85.0d+0
      bl(n+17) = 100.0d+0

!     ----------------------------------------------------------------
!     Set the default upper and lower bounds on the variables x(1:n).
!     ----------------------------------------------------------------
      do j = 1, n
         bl(j) = 0.0d+0
         bu(j) = infBnd
      end do

      bl( 1) =   8.0d+0
      bu( 1) =  21.0d+0

      bl( 2) =  43.0d+0
      bu( 2) =  57.0d+0

      bl( 3) =   3.0d+0
      bu( 3) =  16.0d+0

      bu( 4) =  90.0d+0
      bu( 5) = 120.0d+0
      bu( 6) =  60.0d+0
      bu( 7) =  90.0d+0
      bu( 8) = 120.0d+0
      bu( 9) =  60.0d+0
      bu(10) =  90.0d+0
      bu(11) = 120.0d+0
      bu(12) =  60.0d+0
      bu(13) =  90.0d+0
      bu(14) = 120.0d+0
      bu(15) =  60.0d+0

!     ObjAdd = 0.0 means there is no constant to be added to the QP
!                  objective.

      ObjAdd =   0.0d+0

!     ------------------------------------------------------------------
!     Set the initial estimate of the solution.
!     ------------------------------------------------------------------
      x( 1)  =  20.0d+0
      x( 2)  =  55.0d+0
      x( 3)  =  15.0d+0
      x( 4)  =  20.0d+0
      x( 5)  =  60.0d+0
      x( 6)  =  20.0d+0
      x( 7)  =  20.0d+0
      x( 8)  =  60.0d+0
      x( 9)  =  20.0d+0
      x(10)  =  20.0d+0
      x(11)  =  60.0d+0
      x(12)  =  20.0d+0
      x(13)  =  20.0d+0
      x(14)  =  60.0d+0
      x(15)  =  20.0d+0

      end ! subroutine hs118
