*     ------------------------------------------------------------------
*     File sqmain.f
*     Example of a call to subroutine SQOPT, which is
*     part of the SNOPT package.
*
*     04 Oct 1994: First   version.
*     27 Oct 2003: Updated for SNOPT 7.
*     03 Jul 2010: Indefinite case added.
*     ------------------------------------------------------------------
      program
     &     sqmain

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
     &     Errors, INFO, iObj, iPrint, iSumm, iSpecs, itnlim, iP, iS,
     &     lencObj, m, maxS, mincw, miniw, minrw, n, ncolH, ne, nInf,
     &     nName, nOut, nS
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

         lfile = 'sqmain.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'sqmain.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

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
      write(nOut, 1000)

      call sqdata
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

*     ------------------------------------------------------------------
*     Specify options not set in the Specs file.
*     iP and iS refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      Errors = 0

      maxS   = ncolH + 1
      itnlim = 200
      iP     =  0
      iS     =  0
      call sqseti
     &   ( 'Superbasics Limit', maxS  , iP, iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Iterations',        itnlim, iP, iS, Errors,
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
      write(nOut, *) 'sqmain (explicit obj) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     2. Solve the same problem but with the objective row as part of
*     the constraint matrix A.
*     Use a Cold Start because the dimensions are different.
*     ------------------------------------------------------------------
      write(nOut, 2000)

      Errors = 0                ! Reset count of input errors.

      call sqset
     &   ( 'Solution no    ',        iP,    iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset                ! Required for the following Warm start call
     &   ( 'Sticky parameters yes',  iP,    iS, Errors,
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
      write(nOut, *) 'sqmain (implicit obj) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
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

*     ------------------------------------------------------------------
*     3. Alter some options and call sqopt again, with a Warm start.
*     The following illustrates the use of sqset, sqseti and sqsetr
*     to set specific options.  We could ensure that all unspecified
*     options take default values by first calling
*     sqset ( 'Defaults', ... ).
*     Beware that certain parameters would then need to be redefined.
*     ------------------------------------------------------------------
      write(nOut, 3000)

      Errors = 0

      itnlim = 500
      call sqset
     &   ( 'Defaults',             iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print   frequency', 1, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Summary frequency', 1, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( ' ',                    iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Scale option 0',       iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset
     &   ( 'Solution no    ',      iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Iterations', itnlim,   iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqset                ! Required for the following HOT start call
     &   ( 'Sticky parameters yes',iPrint, iSumm, Errors,
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
      write(nOut, *) 'sqmain (warm start) finished.'
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
      write(nOut, 4000)

      Errors = 0

      call sqset
     &   ( 'QPsolver QN',     iPrint, iSumm, Errors,
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
      write(nOut, *) 'sqmain (QN) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj + x(n+iObj)

*     ------------------------------------------------------------------
*     5. Attempt to solve an indefinite  QP problem.
*        As sqpopt is intended for convex QPs it may or may not find a
*        solution, depending on the starting point.
*     ------------------------------------------------------------------
      call sqdata3
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nName, lencObj, ncolH,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

*     ------------------------------------------------------------------
*     Specify options not set in the Specs file.
*     iP and iS refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      write(nOut, 5000)

      Errors = 0

      maxS   = ncolH + 1
      itnlim = 200
      iP     =  0
      iS     =  0
      call sqseti
     &   ( 'Superbasics Limit', maxS  , iP, iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Iterations',        itnlim, iP, iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print level',            1, iP, iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Print   frequency',      1, iP, iS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Summary frequency',      1, iP, iS, Errors,
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
      write(nOut, *) 'sqmain (indefinite) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'sqOpt INFO    =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =',  ObjAdd + Obj + x(n+iObj)
      if (INFO .gt. 30) go to 910

      stop

 1000 format(/ ' ---------------------------------------------'
     &       //' sqOpt with an explicit linear objective term')
 2000 format(/ ' ---------------------------------------------'
     &       //' sqOpt with an implicit linear objective term')
 3000 format(/ ' ---------------------------------------------'
     &       //' sqOpt with a warm start and different options')
 4000 format(/ ' ---------------------------------------------'
     &       //' sqOpt with a hot start and QN solver')
 5000 format(/ ' ---------------------------------------------'
     &       //' sqOpt applied to a nonconvex QP')

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 9000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

      stop

 9000 format(/  a, 2x, a  )

      end ! program sqmain

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqdata
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
*     sqdat0   defines the problem discusses in the SQOPT Users Guide.
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
     &     infBnd, x0j
*     ------------------------------------------------------------------
      double precision   zero,          one,          half
      parameter         (zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0)
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

      ncolH   = n

      ObjAdd  = zero
      do j = 1, ncolH
         x0j     = half
         ObjAdd  = ObjAdd + x0j*x0j
         cObj(j) = - x0j
      end do

      ObjAdd  = half*ObjAdd
      lencObj = ncolH
      iObj    = 0

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

      end ! subroutine sqdata

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
*     The matrix H is the just the identity.
*
*     ==================================================================
      integer
     &     j
*     ------------------------------------------------------------------
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
      do j = 1, ncolH
         Hx(j) = x(j)
      end do

*     ------------
*     Final entry.
*     ------------
      if (State .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is  sqmain problem 1.   ncolH =', i4)
 2000 format(/ ' Finished sqmain problem 1.')

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
*     sqopt  indefinite example program.
*
*         (1.69  1    2    3    4    5    6    7   )
*         (1     1.69 1    2    3    4    5    6   )
*         (2     1    1.69 1    2    3    4    5   )
*         (3     2    1    1.69 1    2    3    4   )
*     H = (4     3    2    1    1.69 1    2    3   )
*         (5     4    3    2    1    1.69 1    2   )
*         (6     5    4    3    2    1    1.69 1   )
*         (7     6    5    4    3    2    1    1.69)
*
*     ==================================================================
      integer
     &     i, j
      double precision
     &     sum, rij
*     ------------------------------------------------------------------
      integer            nOut
      parameter         (nOut  = 6)
*     ------------------------------------------------------------------

*     First entry.  Print something on the standard output.

      if (State .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000) ncolH
      end if

*     -------------------
*     Normal entry. n = 8
*     -------------------
      do i = 1, ncolH
         sum    = 1.69d+0*x(i)
         do j = 1, ncolH
            rij = abs(i-j)
            sum = sum + rij*x(j)
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

 1000 format(/ ' This is  sqmain (indefinite) problem 2.   ncolH =', i4)
 2000 format(/ ' Finished sqmain (indefinite) problem 2.')

      end ! subroutine userHx3

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
*     sqdat0   defines the problem discussed in the SQOPT Users Guide.
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
     &     infBnd, x0j
*     ------------------------------------------------------------------
      double precision   zero,          one,          half
      parameter         (zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0)
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
*     -inf le  (-h  -h  -h  -h  -h  -h  -h  -h   1   ) le inf
*        1 =   ( 1   1   1   1   1   1   1   1   ) =  1
*     with h = 1/2.
*     ------------------------------------------------------------------
      n       = 30
      m       = n + 1           ! Includes the objective row
      ne      = 2*n + 2*(n-1)

      ncolH   = n

      ObjAdd  = zero
      do j = 1, ncolH
         x0j     = half
         ObjAdd  = ObjAdd + x0j*x0j
      end do
      ObjAdd  = half*ObjAdd

      iObj    = n
      lencObj = 0

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

         neA       =  neA + 1
         indA(neA) =  iObj
         Acol(neA) = -half

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
*     sqdata3   defines an indeinite QP.
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
     &     infBnd
*     ------------------------------------------------------------------
      double precision   one
      parameter         (one = 1.0d+0)
*     ------------------------------------------------------------------
*     Give the problem a name

      Prob   = 'sqProb 3'

*     ------------------------------------------------------------------
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*
*     -1         ( -1   1                           )     inf
*     -1.05      (     -1   1                       )     inf
*     -1.1       (         -1   1                   )     inf
*     -1.15  le  (             -1   1               ) le  inf
*     -1.2       (                 -1   1           )     inf
*     -1.25      (                     -1   1       )     inf
*     -1.3       (                         -1   1   )     inf
*     ------------------------------------------------------------------
      n       = 8
      m       = n               ! Includes the objective row
      ne      = 2*(n-1) + n

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
            Acol(neA) =  one
         endif

         if (j .lt. n) then
            neA       =  neA + 1
            indA(neA) =  j
            Acol(neA) = -one
         end if

         neA       =  neA + 1
         indA(neA) =  iObj
         Acol(neA) =  8 - j

      end do

*     Don't forget to finish off  locA.
*     This is crucial.

      locA(n+1) =  neA + 1

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do i = 1, n

         ! Upper bounds on  x

         bu(i) = i

         ! Upper bounds on Ax

         j     = n + i
         bu(j) = infBnd
      end do

      bl(1)   = -1.0d+0
      bl(2)   = -2.1d+0
      bl(3)   = -3.2d+0
      bl(4)   = -4.3d+0
      bl(5)   = -5.4d+0
      bl(6)   = -6.5d+0
      bl(7)   = -7.6d+0
      bl(8)   = -8.7d+0
      bl(n+1) = -1.0d+0
      bl(n+2) = -1.05d+0
      bl(n+3) = -1.1d+0
      bl(n+4) = -1.15d+0
      bl(n+5) = -1.2d+0
      bl(n+6) = -1.25d+0
      bl(n+7) = -1.3d+0

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
*     x(j) = 0.0d+0 gives "nonconvex QP" error
      do j = 1, n
         x(j) = -j
      end do

      end ! subroutine sqdata3

