!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File hs76ModInf.f    Example program for SQOPT
!
!     Various examples related to problem  Hock-Schittkowski 76.
!
!     Sample program for SQOPT Version 2.1-0   Sep 2014.
!
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      program hs76ModInf

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

         lfile = 'hs76ModInf.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'hs76ModInf.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     First,  sqInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call sqInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      Errors = 0
      iPrt   = 0               ! Suppress printing while reading options
      iSum   = 0

!     ------------------------------------------------------------------
!     1. Solve hs76 modified to have infeasible constraints.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------------'
      write(nOut, *) ' 1. hs76ModInf: hs76 with infeasible constraints '
      write(nOut, *) '    Solved in implicit elastic mode              '
      write(nOut, *) '-------------------------------------------------'

!     Read a specs file.

      call sqSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )
      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

      call hs76ModData
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

      infBnd = 1.0d+20
      call sqseti
     &   ( 'Elastic mode =',      2,      iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqseti
     &   ( 'Elastic objective =', 1,      iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqsetr
     &   ( 'Elastic weight =   ', 1.0d+4, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call sqsetr
     &   ( 'Infinite Bound size =', infBnd, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call sqSet
     &   ( 'Sticky parameters =  yes',      iPrt, iSum, Errors,
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
      write(nOut, *) 'hs76ModInf finished. '
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'sqopt INFO     =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      if (nInf .eq. 0)
     &write(nOut, *) 'Obj            =', objQP
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     2. Solve hs76 modified to have explicit elastic variables.
!        The explicit objective elastic term is 1.0d+5. This value
!        almost minimizes the one-norm of the constraint violations.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------------'
      write(nOut, *) ' 2. hs76ModInf: hs76 with infeasible constraints.'
      write(nOut, *) '    Explicit elastic mode with cObj obj. term    '
      write(nOut, *) '-------------------------------------------------'
      write(nOut, *) ' '

      call hs76ModExData
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

!     No need for elastic mode. The problem is always feasible.
!     Elastic mode      = 1 => enter elastic mode only if infeasible.

      call sqseti
     &   ( 'Elastic mode =',      1, iPrt, iSum, Errors,
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
      write(nOut, *) 'hs76ModInf (Ex) finished.'
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'sqopt INFO     =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', ObjQP
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     3. Solve the modified hs76 with explicit elastic variables.
!        Include the elastic objective elements as part of A.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------------'
      write(nOut, *) ' 3. hs76ModInf: hs76 with infeasible constraints.'
      write(nOut, *) '    Explicit elastic mode with impl. obj. term   '
      write(nOut, *) '-------------------------------------------------'

      call hs76ModEx2Data
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      nb     = n + m

!     No need for elastic mode. The problem is always feasible.

      call sqseti
     &   ( 'Elastic mode =',      1, iPrt, iSum, Errors,
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
      write(nOut, *) 'hs76ModInf (Ex2) finished.'
      write(nOut, *) 'Input errors   =', Errors
      write(nOut, *) 'sqopt INFO     =', INFO
      write(nOut, *) 'nInf           =', nInf
      write(nOut, *) 'sInf           =', sInf
      write(nOut, *) 'Obj            =', ObjQP
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

      end ! program hs76ModInf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs76ModData
     &   ( maxm, maxn, maxne,
     &     m, n, ne, nNames, lencObj, ncolH,
     &     iObjA, objAdd, ProblemName,
     &     Acol, indA, locA, bl, bu, cObj,
     &     Names, eType, hs, x )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, m, n, ne, nNames, lencObj, ncolH, iObjA,
     &     eType(maxn+maxm), hs(maxn+maxm), indA(maxne), locA(maxn+1)
      character
     &     ProblemName*8, Names(maxn+maxm)*8
      double precision
     &     objAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     cObj(maxn), x(maxn+maxm)

!     ==================================================================
!     hs76ModData   defines the problem HS76 modified to have infeasible
!     constraints.
!
!     The 3rd and 4th constraints are inconsistent.
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

      ProblemName = 'HS76Mod '

      n       = 4
      m       = 4               ! Includes 1 infeasible constraint
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

!     A(1,1)  =  1.0d+0
!     A(2,1)  =  3.0d+0
!     A(3,1)  =  0.0d+0
!     A(4,1)  =  0.0d+0
!
!     A(1,2)  =  2.0d+0
!     A(2,2)  =  1.0d+0
!     A(3,2)  =  1.0d+0
!     A(4,2)  = -1.0d+0
!
!     A(1,3)  =  1.0d+0
!     A(2,3)  =  2.0d+0
!     A(3,3)  =  4.0d+0
!     A(4,3)  = -4.0d+0
!
!     A(1,4)  =  1.0d+0
!     A(2,4)  = -1.0d+0
!     A(3,4)  =  0.0d+0
!     A(4,4)  =  0.0d+0

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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -1.0d+0

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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -4.0d+0

!
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

!     bl(n+1) = -infBnd         ! hs76 bounds
!     bu(n+1) =  5.0d+0

!     bl(n+2) = -infBnd
!     bu(n+2) =  4.0d+0

!     bl(n+3) =  1.5d+0
!     bu(n+3) =  infBnd

      bl(n+1) =   8.0+0         ! Modified bounds
      bu(n+1) =  10.0d+0

      bl(n+2) =   4.0+0
      bu(n+2) =   8.0d+0

      bl(n+3) =   2.0d+0
      bu(n+3) =   3.0d+0

      bl(n+4) =   2.0d+0
      bu(n+4) =   3.0d+0

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
      x(1) = 1.5d+0             !Solution
      x(2) = 2.0d+0
      x(3) = 0.0d+0
      x(4) = 2.5d+0

      do j = 1, n
         x(j) = 5.0d-1
      end do

      end ! subroutine hs76ModData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs76ModExData
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
!     hs76ExData   defines the problem HS76 with elastic variables.
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

      double precision   wtInf
      parameter         (wtInf = 1.0d+4)
!     ------------------------------------------------------------------
!     Give the problem a name

      ProblemName = 'HS76Ex '

      n       = 12              ! x(1:4), u(1:4) v(1:4)
      m       =  4              ! Includes 1 infeasible constraint
      ncolH   =  4
      nb      =  n + m

      iObjA   = 0
      objAdd  = 0.0d+0
      lencObj = n

      nNames  = 1

      infBnd  =  1.0d+20

      cObj( 1) =  -1.0d+0
      cObj( 2) =  -3.0d+0
      cObj( 3) =   1.0d+0
      cObj( 4) =  -1.0d+0

      cObj( 5) =   wtInf
      cObj( 6) =   wtInf
      cObj( 7) =   wtInf
      cObj( 8) =   wtInf
      cObj( 9) =   wtInf
      cObj(10) =   wtInf
      cObj(11) =   wtInf
      cObj(12) =   wtInf

!      A(1, 1)   =  1.0d+0
!      A(2, 1)   =  3.0d+0
!
!      A(1, 2)   =  2.0d+0
!      A(2, 2)   =  1.0d+0
!      A(3, 2)   =  1.0d+0
!      A(4, 2)   = -1.0d+0
!
!      A(1, 3)   =  1.0d+0
!      A(2, 3)   =  2.0d+0
!      A(3, 3)   =  4.0d+0
!      A(4, 3)   = -4.0d+0
!
!      A(1, 4)   =  1.0d+0
!      A(2, 4)   = -1.0d+0
!
!      A(1, 5)   = -1.0d+0       ! w(1),  violation of bu(1)
!
!      A(2, 6)   = -1.0d+0       ! w(2),  violation of bu(2)
!
!      A(3, 7)   = -1.0d+0       ! w(3),  violation of bu(3)
!
!      A(4, 8)   = -1.0d+0       ! w(4),  violation of bu(4)
!
!      A(1, 9)   = +1.0d+0       ! v(1),  violation of bl(1)
!
!      A(2,10)   = +1.0d+0       ! v(2),  violation of bl(2)
!
!      A(3,11)   = +1.0d+0       ! v(3),  violation of bl(3)
!
!      A(4,12)   = +1.0d+0       ! v(4),  violation of bl(4)

!     ------------------------------------------
!     General linear constraints  Ax - w + v = 0
!
!     Set up the list of row indices in indA.
!     ------------------------------------------
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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -1.0d+0

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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -4.0d+0

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

!     -------------------------------------------
!     Column 5.
!     ------------------------------------------
      locA( 5) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) = -1.0d+0        ! w(1),  violation of bu(1)

!     -------------------------------------------
!     Column 6.
!     ------------------------------------------
      locA( 6) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) = -1.0d+0        ! w(2),  violation of bu(2)

!     -------------------------------------------
!     Column 7.
!     ------------------------------------------
      locA( 7) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) = -1.0d+0        ! w(3),  violation of bu(3)

!     -------------------------------------------
!     Column 8.
!     ------------------------------------------
      locA( 8) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -1.0d+0        ! w(4),  violation of bu(4)

!     -------------------------------------------
!     Column 9.
!     ------------------------------------------
      locA( 9) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) = +1.0d+0        ! v(1),  violation of bl(1)

!     -------------------------------------------
!     Column 10.
!     ------------------------------------------
      locA(10) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) = +1.0d+0        ! v(2),  violation of bl(2)

!     -------------------------------------------
!     Column 11.
!     ------------------------------------------
      locA(11) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) = +1.0d+0        ! v(3),  violation of bl(3)

!     -------------------------------------------
!     Column 12.
!     ------------------------------------------
      locA(12) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = +1.0d+0        ! v(4),  violation of bl(4)

!!     -----------------------------------------
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

      bl(n+1) =   8.0+0
      bu(n+1) =  10.0d+0
      bl(n+2) =   4.0+0
      bu(n+2) =   8.0d+0
      bl(n+3) =   2.0d+0
      bu(n+3) =   3.0d+0
      bl(n+4) =   2.0d+0
      bu(n+4) =   3.0d+0

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

      end ! subroutine hs76ExData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs76ModEx2Data
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
!     hs76Ex2Data  defines problem HS76 with explicit elastic variables.
!
!     Constant linear objective is stored as row iObjA of A.
!
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants objAdd and c so that the explicit
!         objective is
!             objAdd + A(iObjA,:)'*x + cObj'*x + half*x'*H*x
!     ==================================================================
      integer
     &     i, j, nb
      double precision
     &     infBnd
!     ------------------------------------------------------------------
      double precision   zero,          one,          half
      parameter         (zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0)

      double precision   wtInf
      parameter         (wtInf = 1.0d+5)
!     ------------------------------------------------------------------
!     Give the problem a name

      ProblemName = 'HS76Ex2'

      n       = 12              ! x(1:4), u(1:4) v(1:4)
      m       =  5              ! Includes the linear objective.

      ncolH   =  4
      nb      =  n + m

      iObjA   = 5
      objAdd  = 0.0d+0
      lencObj = 4

      nNames  = 1

      infBnd   =   1.0d+20

      cObj( 1) =  -1.0d+0
      cObj( 2) =  -3.0d+0
      cObj( 3) =   1.0d+0
      cObj( 4) =  -1.0d+0

!      A(1, 1)   =  1.0d+0
!      A(2, 1)   =  3.0d+0
!
!
!      A(1, 2)   =  2.0d+0
!      A(2, 2)   =  1.0d+0
!      A(3, 2)   =  1.0d+0
!      A(4, 2)   = -1.0d+0
!
!      A(1, 3)   =  1.0d+0
!      A(2, 3)   =  2.0d+0
!      A(3, 3)   =  4.0d+0
!      A(4, 3)   = -4.0d+0
!
!      A(1, 4)   =  1.0d+0
!      A(2, 4)   = -1.0d+0
!
!      A(1, 5)   = -1.0d+0
!      A(5, 5)   =  wtInf
!
!      A(2, 6)   = -1.0d+0
!      A(5, 6)   =  wtInf
!
!      A(3, 7)   = -1.0d+0
!      A(5, 7)   =  wtInf
!
!      A(4, 8)   = -1.0d+0
!      A(5, 8)   =  wtInf
!
!      A(1, 9)   = +1.0d+0
!      A(5, 9)   =  wtInf
!
!      A(2,10)   = +1.0d+0
!      A(5,10)   =  wtInf
!
!      A(3,11)   = +1.0d+0
!      A(5,11)   =  wtInf
!
!      A(4,12)   = +1.0d+0
!      A(5,12)   =  wtInf


!     --------------------------------------
!     General linear constraints  Ax - u + v
!
!     Set up the list of row indices in indA.
!     ---------------------------------------
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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -1.0d+0

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

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -4.0d+0

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

!     -------------------------------------------
!     Column 5.
!     ------------------------------------------
      locA( 5) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) = -1.0d+0

      ne       =  ne + 1
      indA(ne) =  5             ! Ojective row
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 6.
!     ------------------------------------------
      locA( 6) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) = -1.0d+0

      ne       =  ne + 1
      indA(ne) =  5             ! Ojective row
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 7.
!     ------------------------------------------
      locA( 7) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) = -1.0d+0

      ne       =  ne + 1
      indA(ne) =  5             ! Ojective row
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 8.
!     ------------------------------------------
      locA( 8) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = -1.0d+0

      ne       =  ne + 1
      indA(ne) =  5             ! Ojective row
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 9.
!     ------------------------------------------
      locA( 9) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) = +1.0d+0

      ne       =  ne + 1
      indA(ne) =  5             ! Ojective row
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 10.
!     ------------------------------------------
      locA(10) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) = +1.0d+0

      ne       =  ne + 1
      indA(ne) =  5
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 11.
!     ------------------------------------------
      locA(11) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  3
      Acol(ne) = +1.0d+0

      ne       =  ne + 1
      indA(ne) =  5
      Acol(ne) =  wtInf

!     -------------------------------------------
!     Column 12.
!     ------------------------------------------
      locA(12) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  4
      Acol(ne) = +1.0d+0

      ne       =  ne + 1
      indA(ne) =  5
      Acol(ne) =  wtInf

!!     -----------------------------------------
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

      bl(n+1) =   8.0+0
      bu(n+1) =  10.0d+0
      bl(n+2) =   4.0+0
      bu(n+2) =   8.0d+0
      bl(n+3) =   2.0d+0
      bu(n+3) =   3.0d+0
      bl(n+4) =   2.0d+0
      bu(n+4) =   3.0d+0
      bl(n+5) =  -infBnd        ! Objective row
      bu(n+5) =   infBnd        ! Objective row

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

      end ! subroutine hs76Ex2Data

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

 1000 format(/ ' This is  hs76ModInf.   ncolH =', i4)
 2000 format(/ ' Finished hs76ModInf.')

      end ! subroutine userHx

