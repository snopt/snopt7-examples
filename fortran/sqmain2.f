*     ------------------------------------------------------------------
*     File sqmain2.f
*     More examples of calls to subroutine SQOPT.
*
*     04 Oct 1994: First   version.
*     27 Oct 2003: Updated for SNOPT 7.
*     23 Apr 2011: Indefinite QP example added.
*     ------------------------------------------------------------------
      program
     &     sqmain2

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
     &     indA(maxne), eType(maxm+maxn), hs(maxm+maxn), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxm+maxn), bu(maxm+maxn), cObj(maxn),
     &     x(maxm+maxn), pi(maxm), rc(maxm+maxn)

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
     &     Errors, INFO, iObj, iPrint, iSumm, iSpecs, itnlim,
     &     iPrt, iSum, j, lencObj, m, maxS, mincw, miniw, minrw,
     &     n, ncolH, ne, nInf, nName, nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
      character
     &     lfile*20
      external
     &     userHx, userHx3

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

         lfile = 'sqmain2.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'sqmain2.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

      Errors = 0
      iPrt   = 0
      iSum   = 0

*     ------------------------------------------------------------------
*     First,  sqInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call sqInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call sqSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

*     ------------------------------------------------------------------
*     1. Solve a QP problem with an explicit linear objective.
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and c so that the explicit
*         objective is
*             ObjAdd + c'*x + half*x'*H*x
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) '1. QP with an explicit linear objective.'
      write(nOut, *) '----------------------------------------'

      call sqdata1
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

*     ------------------------------------------------------------------
*     Specify options not set in the Specs file.
*     iPrt and iSum refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      maxS   = ncolH + 1
      call sqseti
     &   ( 'Superbasics Limit', maxS  , iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      itnlim = 200
      call sqseti
     &   ( 'Iterations',        itnlim, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Solve the QP using a Cold start.
*     hs     need not be set if a basis file is to be input.
*            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
*            The values are used by the Crash procedure s2crsh
*            to choose an initial basis B.
*            If hs(j) = 0 or 1, column j is eligible for B.
*            If hs(j) = 2, column j is initially superbasic (not in B).
*            If hs(j) = 3, column j is eligible for B and is given
*                          preference over columns with hs(j) = 0 or 1.
*            If hs(j) = 4 or 5, column j is initially nonbasic.
*     ------------------------------------------------------------------
      call sqopt
     &   ( 'Cold', userHx, m,
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
      write(nOut, *) 'sqmain2 (1) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     2. Solve the same problem but with the objective row as part of
*     the constraint matrix A.
*     Use a Cold Start because the dimensions are different.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) '2. Same QP with objective as row of A.  '
      write(nOut, *) '----------------------------------------'

      call sqset
     &   ( 'Solution no    ',        iPrt,    iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset    ! Required for following Warm start
     &   ( 'Sticky parameters yes',  iPrt,    iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqdata2
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      call sqopt
     &   ( 'Cold', userHx, m,
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
      write(nOut, *) ' '
      write(nOut, *) 'sqmain2 (2) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     3. Alter some options and call sqopt again, with a Warm start.
*     The following illustrates the use of sqset, sqseti and sqsetr
*     to set specific options.  We could ensure that all unspecified
*     options take default values by first calling
*     sqset ( 'Defaults', ... ).
*     Beware that certain parameters would then need to be redefined.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) '3. Same QP with Warm start.             '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) ' '

      Errors = 0
      itnlim = 500
      call sqset
     &   ( 'Defaults',             iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print   frequency', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Summary frequency', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( ' ',                    iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Scale option 0',       iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Iterations', itnlim,   iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Solution no    ',      iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset    ! Required for following Hot start
     &   ( 'Sticky parameters yes',iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (Errors .gt. 0) then
         write(nOut, *) 'NOTE: Some of the options were not recognized'
      end if

*     ------------------------------------------------------------------
*     Test the Warm start.
*     hs(*) specifies a complete basis from the previous call.
*     A Warm start uses hs(*) directly, without calling Crash.
*
*     Warm starts are normally used after sqopt has solved a
*     problem with the SAME DIMENSIONS but perhaps altered data.
*     Here we have not altered the data, so very few iterations
*     should be required.
*     ------------------------------------------------------------------
      call sqopt
     &   ( 'Warm', userHx, m,
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
      write(nOut, *) 'sqmain2 (3) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     4. Call sqopt with the quasi-Newton solver and a hot start.
*
*     The following illustrates the use of sqset, sqseti and sqsetr
*     to set specific options.  We could ensure that all unspecified
*     options take default values by first calling
*     sqset ( 'Defaults', ... ).
*     Beware that certain parameters would then need to be redefined.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) '4. QP with QN solver and hot start.     '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) ' '

      INFO = 0
      call sqset
     &   ( 'QPsolver QN',     iPrt, iSum, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqopt
     &   ( 'Hot H', userHx, m,
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
      write(nOut, *) 'sqmain2 (4) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)

*     ------------------------------------------------------------------
*     5. Attempt to solve a nonconvex QP.
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and c so that the explicit
*         objective is
*             ObjAdd + c'*x + half*x'*H*x
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) '5. Nonconvex QP (to test error exit)    '
      write(nOut, *) '----------------------------------------'
      write(nOut, *) ' '

      call sqdata3
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

*     ------------------------------------------------------------------
*     Specify options not set in the Specs file.
*     iPrt and iSum refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      call sqset
     &   ( 'Defaults',           iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print   level',     1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print   frequency', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Summary frequency', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Solve the QP using a Cold start.
*     hs     need not be set if a basis file is to be input.
*            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
*            The values are used by the Crash procedure s2crsh
*            to choose an initial basis B.
*            If hs(j) = 0 or 1, column j is eligible for B.
*            If hs(j) = 2, column j is initially superbasic (not in B).
*            If hs(j) = 3, column j is eligible for B and is given
*                          preference over columns with hs(j) = 0 or 1.
*            If hs(j) = 4 or 5, column j is initially nonbasic.
*     ------------------------------------------------------------------
      call sqopt
     &   ( 'Cold', userHx3, m,
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
      write(nOut, *) 'sqmain2 (5) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj
      if (INFO .gt. 30) go to 910

      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

      stop

 4000 format(/  a, 2x, a  )

      end ! program sqmain2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqdata1
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nName, lencObj, ncolH, iObj,
     &      eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     Prob*8, Names(maxn+maxm)*8
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

*     ==================================================================
*     sqdata1   defines the problem discusses in the SQOPT Users Guide.
*
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and c so that the explicit
*         objective is
*             ObjAdd + cObj'*x + half*x'*H*x
*         If necessary, include an additional linear objective terms
*         as row iObj of A.
*     ==================================================================
      integer
     &     i, j, neA
      double precision
     &     infBnd
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
*     Give the problem a name

      Prob   = 'sqProb  '

*     ------------------------------------------------------------------
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*
*     -inf     ( 1  -1                           )    0
*     -inf     (     1  -1                       )    0
*     -inf     (         1  -1                   )    0
*     -inf le  (             1  -1               ) le 0
*     -inf     (                 1  -1           )    0
*     -inf     (                     1  -1       )    0
*     -inf     (                         1  -1   )    0
*        1 =   ( 1   1   1   1   1   1   1   1   ) =  1
*     ------------------------------------------------------------------

      n       = 30
      m       = n               ! Does not include an objective row
      ne      = n + 2*(n-1)

      ObjAdd  = one
      lencObj = n
      iObj    = 0

      ncolH   = n
      nName   = 1

      infBnd  =  1.0d+20

      neA = 0

      do j = 1, n

         ! Set the elements of column j

         locA( j) =  neA + 1

         if (j .gt. 1) then
            neA       =  neA + 1
            indA(neA) =  j   - 1
            Acol(neA) = -one
         endif

         if (j .lt. n) then
            neA       =  neA + 1
            indA(neA) =  j
            Acol(neA) =  one
         end if

         neA       =  neA + 1
         indA(neA) =  m
         Acol(neA) =  one

      end do

      ! Don't forget to finish off  locA.
      ! This is crucial.

      locA(n+1) =  neA + 1

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do i = 1, n-1

         ! Bounds on  x

         bl(i) = zero
         bu(i) = infBnd

         ! Bounds on Ax

         j     = n + i
         bl(j) = -infBnd
         bu(j) =  zero
      end do

      bl(n) = zero
      bu(n) = infBnd

      bl(n+m) =  one
      bu(n+m) =  one

*     ------------------------------------------------------------------
*     Set the objective terms.
*     The objective linear term is explicit.
*     ------------------------------------------------------------------
      cObj(1) = one
      do j = 2, n
         cObj(j) = - cObj(j-1)
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
*     Set the initial estimate of the solution.
*     ------------------------------------------------------------------
      x(1) = -one
      do j = 2, n
         x(j) = -x(j-1)
      end do

      end ! subroutine sqdata1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqdata2
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nName, lencObj, ncolH, iObj,
     &      eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     Prob*8, Names(maxn+maxm)*8
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

*     ==================================================================
*     sqdata2  defines the problem discusses in the SQOPT Users Guide.
*
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and cObj so that the explicit
*         objective is
*             ObjAdd + (row iObj)*x + half*x'*H*x
*     ==================================================================
      integer
     &     i, j, neA
      double precision
     &     infBnd, cj
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
*     Give the problem a name

      Prob   = 'sqProb 2'

*     ------------------------------------------------------------------
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*
*     -inf     ( 1  -1                           )    0
*     -inf     (     1  -1                       )    0
*     -inf     (         1  -1                   )    0
*     -inf le  (             1  -1               ) le 0
*     -inf     (                 1  -1           )    0
*     -inf     (                     1  -1       )    0
*     -inf     (                         1  -1   )    0
*     -inf le  (-1   1  -1   1  -1   1  -1   1   ) le inf
*        1 =   ( 1   1   1   1   1   1   1   1   ) =  1
*     ------------------------------------------------------------------

      n       = 30
      m       = n + 1           ! Includes the objective row
      ne      = 2*n + 2*(n-1)

      ObjAdd  = one
      iObj    = n
      lencObj = 0

      ncolH   = n
      nName   = 1

      infBnd  =  1.0d+20

      neA     = 0

      cj      = one

      do j = 1, n

         ! Set the elements of column j

         locA( j) =  neA + 1

         if (j .gt. 1) then
            neA       =  neA + 1
            indA(neA) =  j   - 1
            Acol(neA) = -one
         endif

         if (j .lt. n) then
            neA       =  neA + 1
            indA(neA) =  j
            Acol(neA) =  one
         end if

         neA       =  neA + 1
         indA(neA) =  iObj
         Acol(neA) =   cj
         cj        =  -cj

         neA       =  neA + 1
         indA(neA) =  m
         Acol(neA) =  one

      end do

*     Don't forget to finish off  locA.
*     This is crucial.

      locA(n+1) =  neA + 1

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do i = 1, n-1

         ! Bounds on  x

         bl(i) = zero
         bu(i) = infBnd

         ! Bounds on Ax

         j     = n + i
         bl(j) = -infBnd
         bu(j) =  zero
      end do

      bl(n) = zero
      bu(n) = infBnd

      bl(n+m) =  one
      bu(n+m) =  one

      bl(n+iObj) = -infBnd
      bu(n+iObj) =  infBnd


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
*     Set the initial estimate of the solution.
*     ------------------------------------------------------------------
      x(1) = -one
      do j = 2, n
         x(j) = -x(j-1)
      end do

      end ! subroutine sqdata2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqdata3
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nName, lencObj, ncolH, iObj,
     &      eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     Prob*8, Names(maxn+maxm)*8
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

*     ==================================================================
*     sqdata3  defines the problem discusses in the SQOPT Users Guide.
*
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and cObj so that the explicit
*         objective is
*             ObjAdd + (row iObj)*x + half*x'*H*x
*     ==================================================================
      integer
     &     i, j, neA
      double precision
     &     infBnd, cj
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
*     Give the problem a name

      Prob   = 'sqProb 3'

*     ------------------------------------------------------------------
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*
*     -inf     ( 1  -1                           )    0
*     -inf     (     1  -1                       )    0
*     -inf     (         1  -1                   )    0
*     -inf le  (             1  -1               ) le 0
*     -inf     (                 1  -1           )    0
*     -inf     (                     1  -1       )    0
*     -inf     (                         1  -1   )    0
*     ------------------------------------------------------------------

      n       = 8
      m       = 7        ! Acol does not include the objective row
      ne      = 14

      ObjAdd  = one
      iObj    = n
      lencObj = 0

      ncolH   = n
      nName   = 1

      infBnd  =  1.0d+20

      neA     = 0

      do j = 1, n

         ! Set the elements of column j

         locA( j) =  neA + 1

         if (j .gt. 1) then
            neA       =  neA + 1
            indA(neA) =  j   - 1
            Acol(neA) = -one
         endif

         if (j .lt. n) then
            neA       =  neA + 1
            indA(neA) =  j
            Acol(neA) =  one
         end if

      end do

      ! Don't forget to finish off  locA.
      ! This is crucial.

      locA(n+1) =  neA + 1

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------

      do  j = 1, n
       ! bl(m+j) = - j - 0.1_rp*(j - 1)
       ! bu(m+j) =   j
         bl(m+j) = - infBnd
         bu(m+j) =   infBnd
         cObj(j) =   8 - j
         x(j)    = - j
      end do

      do i = 1, m
         bl(i)  = - one - 0.05d+1*(i - 1)
         bu(i)  =   infBnd
      end do

      ! iObj   = 0   => exoplicit  linear objective.
      ! ObjAdd = 0.0 => no constant added to the QP objective.

      iObj   = 0
      ObjAdd = zero

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
*     Fix the column variables to be non-elastic and the row  variables
*     to be elastic.
*     ------------------------------------------------------------------
      do j = 1, n
         eType(j) = 0
      end do

      do j = n+1, n+m
         eType(j) = 3
      end do

      end ! subroutine sqdata3

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userHx
     &   ( ncolH, x, Hx, State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     ncolH, State, lencu, leniu, lenru, iu(leniu)
      character
     &     cu(lencu)*8
      double precision
     &     x(ncolH), Hx(ncolH), ru(lenru)

*     ==================================================================
*     This is the user-defined Hessian-vector product for the
*     sqopt example program.
*
*          (  2  -1   0   0   0  0  0 . 0  0 )
*          ( -1   2  -1   0   0  0  0 . 0  0 )
*          (  0  -1   2  -1   0  0  0 . 0  0 )
*          (  0   0   .   .   .  0  0 . 0  0 )
*          (  0   0   0   .   2 -1  0 . 0  0 )
*     H =  (  0   0   0   0  -1  2  0 . 0  0 )   rank(H) = ncolH
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  .   .   .   .   .  .  . . .  . )
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  0   0   0   0   0  0  0 . 0  0 )

*     ==================================================================
      integer
     &     j
      double precision
     &     rho
*     ------------------------------------------------------------------
      double precision   two,            five
      parameter         (two   = 2.0d+0, five  = 5.0d+0)
      integer            nOut
      parameter         (nOut  = 6)
*     ------------------------------------------------------------------

*     First entry.  Print something on the standard output.

      if (State .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000) ncolH
      end if

      rho = five                ! Smoothing parameter

*     -------------
*     Normal entry.
*     -------------
      Hx(1) = - x(2) + two*x(1)
      do j = 2, ncolH-1
         Hx(j) = - x(j+1) + two*x(j) - x(j-1)
      end do
      Hx(ncolH) =  two*x(ncolH) - x(ncolH-1)

      do j = 1, ncolH
         Hx(j) = rho*Hx(j)
      end do

*     ------------
*     Final entry.
*     ------------
      if (State .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is problem  sqmain.   ncolH =', i4)
 2000 format(/ ' Finished         sqmain.')

      end ! subroutine userHx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userHx3
     &   ( ncolH, x, Hx, State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     ncolH, State, lencu, leniu, lenru, iu(leniu)
      character
     &     cu(lencu)*8
      double precision
     &     x(ncolH), Hx(ncolH), ru(lenru)

*     ==================================================================
*     This is the user-defined Hessian-vector product for the
*     sqopt example program.
*
*          (  2  -1   0   0   0  0  0 . 0  0 )
*          ( -1   2  -1   0   0  0  0 . 0  0 )
*          (  0  -1   2  -1   0  0  0 . 0  0 )
*          (  0   0   .   .   .  0  0 . 0  0 )
*          (  0   0   0   .   2 -1  0 . 0  0 )
*     H =  (  0   0   0   0  -1  2  0 . 0  0 )   rank(H) = ncolH
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  .   .   .   .   .  .  . . .  . )
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  0   0   0   0   0  0  0 . 0  0 )
*          (  0   0   0   0   0  0  0 . 0  0 )

*     ==================================================================
      integer
     &     i, j
      double precision
     &     sum
*     ------------------------------------------------------------------
      double precision   two,            five
      parameter         (two   = 2.0d+0, five  = 5.0d+0)
      integer            nOut
      parameter         (nOut  = 6)
*     ------------------------------------------------------------------

*     First entry.  Print something on the standard output.

      if (State .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000) ncolH
      end if

*     -------------
*     Normal entry.
*     -------------
      do i = 1, ncolH
         sum    = 1.69*x(i)
         do  j = 1, ncolH
            sum = sum + abs(i-j)*x(j)
         end do
         Hx(i)  = sum
      end do

*     ------------
*     Final entry.
*     ------------
      if (State .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is problem  sqmain.   ncolH =', i4)
 2000 format(/ ' Finished         sqmain.')

      end ! subroutine userHx
