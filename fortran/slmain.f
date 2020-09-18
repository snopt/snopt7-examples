*     ------------------------------------------------------------------
*     File slmain.f
*     This is an example of solving an LP using a call to SQOPT, which is
*     part of the SNOPT package.
*
*     04 Oct 1994: First   version.
*     09 Jul 2005: Current version.
*     ------------------------------------------------------------------
      program
     &     slmain
      implicit
     &     none
      integer
     &     maxm, maxn, maxne
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1500,
     &     maxne  = 3000 )

      character
     &     Prob*8, Names(maxm+maxn)*8
      integer
     &     indA(maxne) , eType(maxn+maxm), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm), cObj(maxn),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)

*     SQOPT workspace---------------------------------------------------
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
*     ------------------------------------------------------------------
      logical
     &     byname
      integer
     &     Errors, INFO, iObj, iPrint, iSumm, iSpecs, itnlim, i1, i2, j,
     &     lencObj, m, mincw, miniw, minrw, n, ncolH, ne, nInf, nName,
     &     nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
      character
     &     lfile*20
      external
     &     nullHx
*     ------------------------------------------------------------------
*     Specify some of the SQOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by main.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'slmain.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'slmain.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call sqInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     1. Solve an LP problem with an explicit linear objective.
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and cObj so that the explicit
*         objective is
*             ObjAdd + cObj'*x
*     ------------------------------------------------------------------
      call sldat1
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, x )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call sqSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

*     ------------------------------------------------------------------
*     Fix the column variables to be non-elastic and the row  variables
*     to be elastic.
*     ------------------------------------------------------------------
      do j = 1, n
         eType(j) = 0
      end do

      do j = n+1, n+m
         eType(j) = 3
      end do

*     ------------------------------------------------------------------
*     Set the initial value and status of each variable.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         hs(j) = 0
      end do

*     ------------------------------------------------------------------
*     Specify options that were not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      Errors = 0

      itnlim = 40
      i1     =  0
      i2     =  0
      call sqseti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
*     hs     need not be set if a basis file is to be input.
*            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
*            The values are used by the Crash procedure s2crsh
*            to choose an initial basis B.
*            If hs(j) = 0 or 1, column j is eligible for B.
*            If hs(j) = 2, column j is initially superbasic (not in B).
*            If hs(j) = 3, column j is eligible for B and is given
*                          preference over columns with hs(j) = 0 or 1.
*            If hs(j) = 4 or 5, column j is initially nonbasic.
*
*     No user-workspace, so we use cw, iw, rw.
*     ------------------------------------------------------------------
      call sqopt
     &   ( 'Cold', nullHx, m,
     &     n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'diet problem (explicit obj) finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     2. Solve the same problem but with the objective row as part of
*        the constraint matrix A.
*        Use a Cold Start because the dimensions are different.
*     ------------------------------------------------------------------
      call sldat2
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, x )

*     ------------------------------------------------------------------
*     Set the initial status of each variable.
*     For fun, we use the hs values from the previous run.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      call sqset      ! Required for following Warm start call
     &   ( 'Sticky parameters yes', iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqopt
     &   ( 'Cold', nullHx, m,
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
      write(nOut, *) 'diet problem (implicit obj) finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj        =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj        =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     Alter some options and call sqopt again, testing the Warm start.
*     The following illustrates the use of snset, snseti and snsetr
*     to set specific options.  If necessary, we could ensure that
*     all unspecified options take default values
*     by first calling snset ( 'Defaults', ... ).
*     Beware that certain parameters would then need to be redefined.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) 'Alter options and test Warm start:'

      Errors = 0
      itnlim = 500
      call sqset
     &   ( ' ',                          iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Print  level     0',         iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Scale option     0',         iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Iterations',         itnlim, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (Errors .gt. 0) then
         write(nOut, *) 'NOTE: Some of the options were not recognized'
      end if

*     Test the Warm start.
*     hs(*) specifies a complete basis from the previous call.
*     A Warm start uses hs(*) directly, without calling Crash.
*
*     Warm starts are normally used after sqopt has solved a
*     problem with the SAME DIMENSIONS but perhaps altered data.
*     Here we have not altered the data, so very few iterations
*     should be required.

      call sqopt
     &   ( 'Warm', nullHx, m,
     &     n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '
      write(nOut, *) 'diet problem (warm start) finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'sqOpt  INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *)
     &               'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *)
     &               'Obj           =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program slmain

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sldat1
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
     &     x(maxn+maxm), cObj(maxn)

*     ------------------------------------------------------------------
*     Diet problem with explicit linear objective.
*
*     Define the problem.
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and cObj so that the objective
*         to be minimized is   ObjAdd + cObj'*x
*     ------------------------------------------------------------------
      integer
     &     i, j
      double precision
     &     infBnd
*     ------------------------------------------------------------------
      double precision    zero
      parameter         ( zero   =0.0d+0)
*     ------------------------------------------------------------------

*     Give the problem a name.

      Prob   = 'Diet 1..'
      infBnd =  1.0d+20

*     ------------------------------------------------------------------
*     This is the Diet problem of Chvatal, 1983.
*     Assign the constraint nonzeros to a, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 )
*
*     c' =  (   3   24   13    9   20   19 )
*

      n       =  6
      m       =  4     ! Includes the objective row
      ne      = 24     ! n*m for a dense A.

      nName   =  1

      lencObj =  6
      ncolH   =  0

*     Column 1.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) =  2
      indA( 3) =  3
      indA( 4) =  4

      Acol( 1) =   110.0d+0
      Acol( 2) =     4.0d+0
      Acol( 3) =     2.0d+0
      Acol( 4) =     3.0d+0

*     Column 2.

      locA( 2) =  5

      indA( 5) =  1
      indA( 6) =  2
      indA( 7) =  3
      indA( 8) =  4

      Acol( 5) =   205.0d+0
      Acol( 6) =    32.0d+0
      Acol( 7) =    12.0d+0
      Acol( 8) =    24.0d+0

*     Column 3.

      locA( 3) =  9

      indA( 9) =  1
      indA(10) =  2
      indA(11) =  3
      indA(12) =  4

      Acol( 9)  =   160.0d+0
      Acol(10)  =    13.0d+0
      Acol(11)  =    54.0d+0
      Acol(12)  =    13.0d+0

*     Column 4.

      locA( 4) = 13

      indA(13) =  1
      indA(14) =  2
      indA(15) =  3
      indA(16) =  4

      Acol(13) =   160.0d+0
      Acol(14) =     8.0d+0
      Acol(15) =   285.0d+0
      Acol(16) =     9.0d+0

*     Column 5.

      locA( 5) = 17

      indA(17) =  1
      indA(18) =  2
      indA(19) =  3
      indA(20) =  4

      Acol(17) =   420.0d+0
      Acol(18) =     4.0d+0
      Acol(19) =    22.0d+0
      Acol(20) =    20.0d+0

*     Column 6.

      locA(6)  = 21

      indA(21) =  1
      indA(22) =  2
      indA(23) =  3
      indA(24) =  4

      Acol(21) =   260.0d+0
      Acol(22) =    14.0d+0
      Acol(23) =    80.0d+0
      Acol(24) =    19.0d+0

*     Don't forget to finish off  locA(n+1) = ne+1
*     This is crucial.

      locA( 7) = 25

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = zero
      end do

      bu(1)  =     4.0d+0
      bu(2)  =     3.0d+0
      bu(3)  =     2.0d+0
      bu(4)  =     8.0d+0
      bu(5)  =     2.0d+0
      bu(6)  =     2.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on  Ax.
*     The last row is free (i.e., infinite upper and lower bounds).
*     ------------------------------------------------------------------
      bl( 7) =  2000.0d+0
      bl( 8) =    55.0d+0
      bl( 9) =   800.0d+0
      bl(10) = - infBnd

      do i = 1, m
         j     = n + i
         bu(j) = infBnd
      end do

      cObj(1) =   3.0d+0
      cObj(2) =  24.0d+0
      cObj(3) =  13.0d+0
      cObj(4) =   9.0d+0
      cObj(5) =  20.0d+0
      cObj(6) =  19.0d+0
      iObj    =   0
      ObjAdd  = zero

*     ------------------------------------------------------------------
*     Set the initial estimate of the solution.
*     ------------------------------------------------------------------
      x(1) = 1.0d+0
      x(2) = 1.0d+0
      x(3) = 1.0d+0
      x(4) = 1.0d+0
      x(5) = 1.0d+0
      x(6) = 1.0d+0

      end ! subroutine sldat1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sldat2
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
     &     x(maxn+maxm), cObj(maxn)

*     ------------------------------------------------------------------
*     Define an LP problem with objective row as row iObj of Acol.
*
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and iObj so that the objective
*         to be minimized is   ObjAdd + x(n+iObj)
*     ------------------------------------------------------------------
      integer
     &     i, j
      double precision
     &     infBnd
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   =0.0d+0)
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'Diet 2..'

      infBnd =  1.0d+20

*     ------------------------------------------------------------------
*     This is the Diet problem of Chvatal, 1983.
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 )
*

      n       =  6
      m       =  4     ! Includes the objective row
      ne      = 24     ! n*m for a dense A.

      nName   =  1

      lencObj =  0
      ncolH   =  0

*     Column 1.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) =  2
      indA( 3) =  3
      indA( 4) =  4

      Acol( 1)  =   110.0d+0
      Acol( 2)  =     4.0d+0
      Acol( 3)  =     2.0d+0
      Acol( 4)  =     3.0d+0

*     Column 2.

      locA( 2) =  5

      indA( 5) =  1
      indA( 6) =  2
      indA( 7) =  3
      indA( 8) =  4

      Acol( 5)  =   205.0d+0
      Acol( 6)  =    32.0d+0
      Acol( 7)  =    12.0d+0
      Acol( 8)  =    24.0d+0

*     Column 3.

      locA( 3) =  9

      indA( 9) =  1
      indA(10) =  2
      indA(11) =  3
      indA(12) =  4

      Acol( 9)  =   160.0d+0
      Acol(10)  =    13.0d+0
      Acol(11)  =    54.0d+0
      Acol(12)  =    13.0d+0

*     Column 4.

      locA( 4) = 13

      indA(13) =  1
      indA(14) =  2
      indA(15) =  3
      indA(16) =  4

      Acol(13)  =   160.0d+0
      Acol(14)  =     8.0d+0
      Acol(15)  =   285.0d+0
      Acol(16)  =     9.0d+0

*     Column 5.

      locA( 5) = 17

      indA(17) =  1
      indA(18) =  2
      indA(19) =  3
      indA(20) =  4

      Acol(17)  =   420.0d+0
      Acol(18)  =     4.0d+0
      Acol(19)  =    22.0d+0
      Acol(20)  =    20.0d+0

*     Column 6.

      locA(6)  = 21

      indA(21) =  1
      indA(22) =  2
      indA(23) =  3
      indA(24) =  4

      Acol(21)  =   260.0d+0
      Acol(22)  =    14.0d+0
      Acol(23)  =    80.0d+0
      Acol(24)  =    19.0d+0

*     Don't forget to finish off  locA(n+1) = ne+1
*     This is crucial.

      locA( 7) = 25

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j)  = zero
      end do

      bu(1)  =     4.0d+0
      bu(2)  =     3.0d+0
      bu(3)  =     2.0d+0
      bu(4)  =     8.0d+0
      bu(5)  =     2.0d+0
      bu(6)  =     2.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on  Ax.
*     The last row is free (i.e., infinite upper and lower bounds).
*     ------------------------------------------------------------------
      bl( 7) =  2000.0d+0
      bl( 8) =    55.0d+0
      bl( 9) =   800.0d+0

      do i = 1, m
         j     = n + i
         bu(j) = infBnd
      end do

      bl(10) = - infBnd
      bu(10) =   infBnd

      iObj   = 4
      ObjAdd = zero

*     ------------------------------------------------------------------
*     Set the initial estimate of the solution.
*     ------------------------------------------------------------------
      x(1) = 0.0d+0
      x(2) = 0.0d+0
      x(3) = 0.0d+0
      x(4) = 0.0d+0
      x(5) = 0.0d+0
      x(6) = 0.0d+0

      end ! subroutine sldat2

