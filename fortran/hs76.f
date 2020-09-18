!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File hs76.f    Example program for SQOPT
!
!     Various examples related to problem  Hock-Schittkowski 76.
!
!     Sample program for SQOPT Version 2.1-0   Sep 2014.
!
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      program hs76

      implicit
     &     none
      integer
     &     maxm, maxn, maxne
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1500,
     &     maxne  = 3000 )

      character
     &     ProblemName*8, Names(maxm+maxn)*8
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
     &     Errors, INFO, iObjA, iPrint, iSumm, iSpecs, itnlim,
     &     iPrt, iSum, lencObj, m, maxS, mincw, miniw, minrw, n, nb,
     &     ncolH, ne, nInf, nNames, nOut, nS
      double precision
     &     infBnd, objQP, objAdd, sInf
      character
     &     lfile*20
      external
     &     userHx

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

         lfile = 'hs76.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'hs76.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     First,  sqInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call sqInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     1. Solve hs76, a QP with an explicit linear objective.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '---------------------'
      write(nOut, *) ' 1. Solve the QP hs76'
      write(nOut, *) '---------------------'

!     ------------------------------------------------------------------
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants objAdd and c so that the explicit
!         objective is
!             objAdd + c'*x + half*x'*H*x
!     ------------------------------------------------------------------
      call hs76Data
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

      Errors = 0
      iPrt   = 0               ! Suppress printing while reading options
      iSum   = 0

      infBnd = 1.0d+20
      call sqsetr
     &   ( 'Infinite Bound size =', infBnd, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqSet
     &   ( 'Sticky parameters =  yes',      iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     Read a specs file.

      call sqSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

      call sqopt
     &   ( 'Cold', userHx, m,
     &     n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, objQP,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs76 finished.  '
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'hs76   INFO    =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', objQP
      write(nOut, *) ' '
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     2. Solve hs76 in elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------'
      write(nOut, *) ' 2. hs76 starting in elastic mode'
      write(nOut, *) '---------------------------------'

      call hs76EData
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

      call sqseti
     &   ( 'Elastic mode =',      2, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Elastic objective =', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqopt
     &   ( 'Cold', userHx, m,
     &     n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, objQP,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs76 (in elatic mode) finished. '
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'hs76   INFO    =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', objQP
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     3. Solve hs76 with the quasi-Newton solver.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) ' 3. Solve hs76 with the quasi-Newton solver'
      write(nOut, *) '-------------------------------------------'

      call hs76Data
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

      call sqseti
     &   ( 'Elastic mode =',      1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Elastic objective =', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqset
     &   ( 'QPsolver     QN',        iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqopt
     &   ( 'Cold', userHx, m,
     &     n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, objQP,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs76 (QN option) finished.  '
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'hs76   INFO    =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', objQP
      write(nOut, *) ' '
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     4. Solve hs76 with the quasi-Newton solver.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) ' 4. Solve hs76 with the quasi-Newton solver'
      write(nOut, *) '    Starting in elastic mode               '
      write(nOut, *) '-------------------------------------------'

      call hs76Data
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

      call sqseti
     &   ( 'Elastic mode =',      2, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Elastic objective =', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqopt
     &   ( 'Cold', userHx, m,
     &     n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj, Names,
     &     eType, hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, objQP,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs76 finished.  '
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'hs76   INFO    =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', objQP
      write(nOut, *) ' '
      if (INFO .gt. 30) go to 910


      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 9000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

      stop

 9000 format(/  a, 2x, a  )

      end ! program hs76ModMain

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs76Data
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nNames, lencObj, ncolH, iObjA,
     &      eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     ProblemName*8, Names(maxn+maxm)*8
      double precision
     &     objAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

!     ==================================================================
!     hs76Data   defines the problem HS76.
!
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants objAdd and c so that the explicit
!         objective is
!             objAdd + cObj'*x + half*x'*H*x
!     ==================================================================
      integer
     &     i, j, nb
      double precision
     &     infBnd
!     ------------------------------------------------------------------
      double precision   zero,          one,          half
      parameter         (zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0)
!     ------------------------------------------------------------------
!     Give the problem a name

      ProblemName = 'HS76    '

      n       = 4
      m       = 3
      ncolH   = n
      nb      = n + m

      iObjA   = 0
      objAdd  = 0.0d+0
      lencObj = ncolH

      nNames  = 1

      infBnd  =  1.0d+20

      cObj(1) = -1.0d+0
      cObj(2) = -3.0d+0
      cObj(3) =  1.0d+0
      cObj(4) = -1.0d+0

!     -------------------------------------
!     General linear constraints
!
!     Set up the list of row indices in indA.
!     -------------------------------------
!     Column  1.

      ne = 0
      locA( 1) =  1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  3.0d+0

!     -------------------------------------------
!     Column 2.
!     -------------------------------------------
      locA( 2) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) =  1.0d+0

!     ------------------------------------------
!     Column 3.
!     ------------------------------------------
      locA( 3) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) =  4.0d+0

!     -------------------------------------------
!     Column 4.
!     ------------------------------------------
      locA( 4) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) = -1.0d+0

!     -----------------------------------------
!     Don't forget to finish off locA.
!     This is crucial.
!     -----------------------------------------
      locA(n+1)= ne + 1

!     ------------------------------------------------------------------
!     Set the upper and lower bounds on the variables
!     ------------------------------------------------------------------
      do i = 1, n
         bl(i) = zero
         bu(i) = infBnd
      end do

      bl(n+1) = -infBnd
      bl(n+2) = -infBnd
      bl(n+3) =  1.5d+0

      bu(n+1) =  5.0d+0
      bu(n+2) =  4.0d+0
      bu(n+3) =  infBnd

!     ------------------------------------------------------------------
!     Set the initial value and state of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash will set the rest.
!     ------------------------------------------------------------------
      do j = 1, n
         hs(j) = 0
      end do

!     ------------------------------------------------------------------
!     Fix all constraints to be nonelastic.
!     ------------------------------------------------------------------
      do j = 1, nb
         eType(j) = 0
      end do

!     ------------------------------------------------------------------
!     Set the initial estimate of the solution.
!     ------------------------------------------------------------------
      do j = 1, n
         x(j) = 5.0d-1
      end do

      end ! subroutine hs76Data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs76EData
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nNames, lencObj, ncolH, iObjA,
     &      eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     ProblemName*8, Names(maxn+maxm)*8
      double precision
     &     objAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

!     ==================================================================
!     hs76EData   defines the problem HS76.
!
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants objAdd and c so that the explicit
!         objective is
!             objAdd + cObj'*x + half*x'*H*x
!     ==================================================================
      integer
     &     i, j, nb
      double precision
     &     infBnd
!     ------------------------------------------------------------------
      double precision   zero,          one,          half
      parameter         (zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0)
!     ------------------------------------------------------------------
!     Give the problem a name

      ProblemName = 'HS76 E  '

      n       = 4
      m       = 3
      ncolH   = n
      nb      = n + m

      iObjA   = 0
      objAdd  = 0.0d+0
      lencObj = ncolH

      nNames  = 1

      infBnd  =  1.0d+20

      cObj(1) = -1.0d+0
      cObj(2) = -3.0d+0
      cObj(3) =  1.0d+0
      cObj(4) = -1.0d+0

!     ---------------------------
!     General linear constraints.
!
!     Set up the list of row indices in indA.
!     -------------------------------------
!     Column  1.

      ne = 0
      locA( 1) =  1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  3.0d+0

!     -------------------------------------------
!     Column 2.
!     -------------------------------------------
      locA( 2) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) =  1.0d+0

!     ------------------------------------------
!     Column 3.
!     ------------------------------------------
      locA( 3) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  2.0d+0

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) =  4.0d+0

!     -------------------------------------------
!     Column 4.
!     ------------------------------------------
      locA( 4) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  1.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) = -1.0d+0

!     -----------------------------------------
!     Don't forget to finish off locA.
!     This is crucial.
!     -----------------------------------------
      locA(n+1)= ne + 1

!     ------------------------------------------------------------------
!     Set the upper and lower bounds on the variables
!     ------------------------------------------------------------------
      do i = 1, n
         bl(i) = zero
         bu(i) = infBnd
      end do

      bl(n+1) = -infBnd
      bl(n+2) = -infBnd
      bl(n+3) =  1.5d+0

      bu(n+1) =  5.0d+0
      bu(n+2) =  4.0d+0
      bu(n+3) =  infBnd

!     ------------------------------------------------------------------
!     Set the initial value and state of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash will set the rest.
!     ------------------------------------------------------------------
      do j = 1, n
         hs(j) = 0
      end do

!     ------------------------------------------------------------------
!     Fix the variable bounds to be non-elastic and the constraints
!     to be elastic.
!     ------------------------------------------------------------------
      do j = 1, n
         eType(j) = 0
      end do

      do j = n+1, nb
         eType(j) = 3
      end do

!     ------------------------------------------------------------------
!     Set the initial estimate of the solution.
!     ------------------------------------------------------------------
      do j = 1, n
         x(j) = 5.0d-1
      end do

      end ! subroutine hs76EData

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userHx
     &   ( ncolH, x, Hx, qpState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     ncolH, qpState, lencu, leniu, lenru, iu(leniu)
      character
     &     cu(lencu)*8
      double precision
     &     x(ncolH), Hx(ncolH), ru(lenru)

!     ==================================================================
!     This is the user-defined Hessian-vector product for hs76.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     H(4,4)
!     ------------------------------------------------------------------
      integer            nOut
      parameter         (nOut  = 6)
!     ------------------------------------------------------------------

!     First entry.  Print something on the standard output.

      if (qpState .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000) ncolH
      end if

      H(1,1)   =   2.0d+0
      H(2,1)   =   0.0d+0
      H(3,1)   =  -1.0d+0
      H(4,1)   =   0.0d+0

      H(1,2)   =   0.0d+0
      H(2,2)   =   1.0d+0
      H(3,2)   =   0.0d+0
      H(4,2)   =   0.0d+0

      H(1,3)   =  -1.0d+0
      H(2,3)   =   0.0d+0
      H(3,3)   =   2.0d+0
      H(4,3)   =   1.0d+0

      H(1,4)   =   0.0d+0
      H(2,4)   =   0.0d+0
      H(3,4)   =   1.0d+0
      H(4,4)   =   1.0d+0

!     -------------
!     Normal entry.
!     -------------
      do i = 1, ncolH
         hx(i) = 0.0d+0
         do j = 1, ncolH
            Hx(i) = Hx(i) + H(i,j)*x(j)
         end do
      end do

!     ------------
!     Final entry.
!     ------------
      if (qpState .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is  hs76.   ncolH =', i4)
 2000 format(/ ' Finished hs76.')

      end ! subroutine userHx

