!***********************************************************************
!                                                                      *
!     File  lpsolver.f                                                 *
!                                                                      *
!     Generic main program for stand-alone SQOPT.                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!                               S Q O P T                              *
!                                                                      *
!    Sparse Quadratic Optimization                                     *
!                                                                      *
!                      Version 7.2                     Jul 1, 2005     *
!                                                                      *
!    Philip E. Gill    Walter  Murray          Michael A. Saunders     *
!    UC San Diego      Stanford University     Stanford University     *
!                                                                      *
!                                                                      *
!    (C) 1992--2005  Regents of the University of California           *
!                    and the Trustees of Stanford University           *
!                                                                      *
!     This software is NOT in the public domain. Its use is governed   *
!     by a license agreement with either the University of California  *
!     or Stanford University.  It is a breach of copyright to make     *
!     copies except as authorized by the license agreement.            *
!                                                                      *
!     This material is based upon work partially supported by the      *
!     National Science Foundation under Grants DMI-9204208 and         *
!     DMI-9204547; and the Office of Naval Research Grant              *
!     N00014-90-J-1242.                                                *
!***********************************************************************
!
!  SQOPT Fortran source files:
!
!  1. lpsolver   Main program
!  2. sq02lib    SQOPT routines and auxiliaries
!  3. sn02lib    SNOPT routines and auxiliaries
!  4. sn03prnt   Print routine
!  5. sn10mach   Machine-dependent routines
!  6. sn15blas   Level-1 Basic Linear Algebra Subprograms (a subset)
!  7. sn17util   linear algebra subprograms
!  8. sn20amat   Core allocation and manipulation of the ( A -I )
!  9. sn25bfac   Basis factorization routines
! 10. sn27LU     LU factorization routines
! 11. sn30spec   SPECS file routines
! 12. sn35mps    MPS file routines
! 13. sn37wrap   Interface and argument checking routines
! 14. sn40bfil   Basis file and solution output routines
! 15. sn50lp     Routines for the primal simplex method
! 16. sn55qp     Routines for quadratic programming
! 17. sn56qncg   QN and CG routines
! 18. sn57qopt   QP and Memory allocation routines called by SQOPT
! 19. sn65rmod   For maintaining R, the approximate reduced Hessian
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      program            sqmps

!     ------------------------------------------------------------------
!     This is the default main program for SQOPT
!     It provides all of the necessary workspace.
!     ------------------------------------------------------------------
      integer            lencw,          leniw,          lenrw
      parameter         (lencw = 150000, leniw = 400000, lenrw = 600000)
      character          cw(lencw)*8
      integer            iw(leniw)
      double precision   rw(lenrw)

      call sqmps1
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      end ! program sqmps

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqmps1
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqmps1 is used for the stand-alone version of the optimizer.
!     It is called by the main program (or equivalent driver).
!     It repeatedly looks for a new problem in the SPECS file
!     and asks for it to be solved, until s3file returns inform gt 1,
!     which means an ENDRUN card was found in the SPECS file,
!     or end-of-file was encountered.
!
!     15 Nov 1991: First version based on Minos 5.4 routine minos1.
!     11 Nov 2000: Updated for SNOPT 6.1.
!     26 Oct 2003: Updated for SNOPT 7.1.
!     05 Oct 2005: lvlTim set correctly.
!     17 Jun 2007: s3unsetOpt used to reset the optional parameters.
!     ==================================================================
      character
     &     cStart*8, Solver*6, str*80, str2*80, title*30
      logical
     &     GotR, PrintMem
      integer
     &     calls, Errors, flMax, InfBnd, INFO, inform, iObj,
     &     iPrint, iSpecs, iSumm, lAcol, lbl, lbu, lenR, lenrhs, lenx0,
     &     lkx, lgObj, leType, lhs, liwEst, lrwEst, lindA, llocA,
     &     lNames, loop, lpi, lrc,
     &     lvlStart, lvlTim, lx, maxn, mincw, miniw, minrw, n, neA,
     &     nextcw, nextiw, nextrw, nnCon, nnJac, nnH, nnL, nnObj, nb,
     &     ngObj, ngObj0, ngQP, nkx, nInf, nInfE, nlocA, nName,
     &     nnH0, nrhs, nS, nx0, m, maxm, maxR, maxS, maxneA, maxru,
     &     maxiu, maxcu, maxrw, maxiw, maxcw, startType
      integer
     &     neH, nlocH, indH(1), locH(1)
      double precision
     &     Hcol(1)
      double precision
     &     objAdd, objQP, objTrue, rhs(1), sInf, sInfE, x0(1),
     &     s1flmx
      external
     &     s1flmx, s3opt, qpHx, sqHx, sqlog
!     ------------------------------------------------------------------
      integer            DefltF,     OpenF
      parameter         (DefltF = 0, OpenF  = 1)

      parameter         (flMax     =   8) ! est. of the largest pos. real
      parameter         (InfBnd    =  70) ! definition of an infinite bound

      parameter         (maxru     =   2) ! maxru+1  starts SNOPT  rw
      parameter         (maxiu     =   4) ! maxiu+1  starts SNOPT  iw
      parameter         (maxcu     =   6) ! maxcu+1  starts SNOPT  cw

      parameter         (maxrw     =   3) ! end of SNOPT part of rw
      parameter         (maxiw     =   5) ! end of SNOPT part of iw
      parameter         (maxcw     =   7) ! end of SNOPT part of cw
      parameter         (iSpecs    =  11) ! Specs (options) file
      parameter         (iPrint    =  12) ! Print file
      parameter         (iSumm     =  13) ! Summary file
      parameter         (nnJac     =  21) ! # nonlinear Jacobian vars
      parameter         (nnObj     =  22) ! # variables in gObj
      parameter         (nnCon     =  23) ! # of nonlinear constraints
      parameter         (nnL       =  24) !   max( nnObj, nnJac )
      parameter         (ngObj     =  26) ! length of QP constant vector
      parameter         (nnH       =  27) ! # QP Hessian columns
      parameter         (lvlStart  =  69) ! = 0(1) => cold(warm) start
      parameter         (lvlTim    = 182) ! Timing level

      parameter         (maxm      = 133) ! Est. number of rows
      parameter         (maxn      = 134) ! Est. number of columns
      parameter         (maxneA    = 135) ! Est. number of elements
!     -------------------------------------------------------------------
      Solver = 'SQOPT '
      INFO   = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         inform = 81       ! Work arrays must have at least 500 elements
         call snWRAP( inform, Solver, str, str2, iw, leniw )
         go to 999
      end if

      call s3unsetAll
     &   ( cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Define global files (reader, printer, etc.)
!     ------------------------------------------------------------------
      iw(iSpecs) = 4
      iw(iPrint) = 9
      iw(iSumm ) = 6
      call s1file( DefltF, iw, leniw )

!     ==================================================================
!     Loop through each problem in the SPECS file.
!     ==================================================================
      do loop  = 1, 100000
         calls = loop

         call s3unsetPrm
     &      ( cw, lencw, iw, leniw, rw, lenrw )

!        ---------------------------------------------------------------
!        Initialize some global values.
!        ---------------------------------------------------------------
         rw(flmax)  = s1flmx( )
         rw(InfBnd) = 1.0d+20
         iw(lvlTim) = 3


         iw(maxru) = 500        ! rw(1:500) contains sqopt variables
         iw(maxiu) = 500        ! iw(1:500) contains sqopt variables
         iw(maxcu) = 500        ! cw(1:500) contains sqopt variables
         iw(maxrw) = lenrw
         iw(maxiw) = leniw
         iw(maxcw) = lencw

!        Initialize timers.

         iw(lvlTim) = 1
         call s1time( 0, 0, iw, leniw, rw, lenrw )

!        Initialize dimensions that can be set in the specs file.

         iw(nnCon )  =   0
         iw(nnJac )  =   0
         iw(nnObj )  =   0
         iw(ngObj )  =   0
         iw(nnL   )  =   0
         iw(nnH   )  =   0

         iw(maxm  )  =   0
         iw(maxn  )  =   0
         iw(maxneA)  =   0

!        ---------------------------------------------------------------
!        Define the SQOPT title and read the Specs file.
!        ---------------------------------------------------------------
         call sqTitle( title )
         call s1init ( title, iw, leniw, rw, lenrw )
         call s3file
     &      ( INFO,
     &        calls, iw(iSpecs), s3opt,
     &        title, iw(iPrint), iw(iSumm), Errors,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (INFO .ge. 2) then
            INFO = 100 + INFO
            go to 999
         end if

         call s1file
     &      ( OpenF, iw, leniw  )

!        Set undefined MPS options to their default values.

         call s3dflt
     &      ( cw, lencw, iw, leniw, rw, lenrw )

!        ---------------------------------------------------------------
!        Check memory limits and fetch the workspace starting positions.
!        ---------------------------------------------------------------
         call s2Mem0
     &      ( INFO,
     &        Solver, lencw, leniw, lenrw, iw,
     &        mincw, miniw, minrw, iw(maxcw), iw(maxiw), iw(maxrw),
     &        nextcw, nextiw, nextrw )
         if (INFO .gt. 0) go to 999 ! Exit without printing

!        ---------------------------------------------------------------
!        Input data in MPS format.
!        Get values of the problem dimensions:
!           maxm , maxn , maxneA
!           lenR , maxS
!           ngObj
!        Initialize  locA, indA, Acol,
!                    bl, bu, iObj, and  objAdd.
!        Compute the array pointers accordingly.
!        ---------------------------------------------------------------
         call s1time
     &      ( 1, 0, iw, leniw, rw, lenrw )
         call s3inpt
     &      ( INFO,
     &        iw(maxm), iw(maxn), iw(maxneA),
     &        iw(nnCon), iw(nnJac), iw(nnObj),
     &        m, n, neA, iObj, objAdd,
     &        nextcw, nextiw, nextrw,
     &        mincw, miniw, minrw,
     &        cw, lencw, iw, leniw, rw, lenrw )
         call s1time
     &      (-1, 0, iw, leniw, rw, lenrw )
         if (INFO .ne. 0) go to 800

!        Fetch the addresses of the problem arrays (set in s3inpt).

         lAcol   = iw(256) ! Acol(neA)   = Constraints by columns
         llocA   = iw(257) ! locA(n+1)   = column pointers for indA
         lindA   = iw(258) ! indA(neA) holds the row indices for Aij
         lbl     = iw(271) ! bl(nb)      = lower bounds
         lbu     = iw(272) ! bu(nb)      = upper bounds
         lx      = iw(299) ! x(nb)       = the solution (x,s)
         lpi     = iw(279) ! pi(m)       = the pi-vector
         lhs     = iw(282) ! the column state vector
         leType  = iw(283) ! eType(nb) definition of elastic vars
         lNames  = iw(359) ! Names(nName)

         call iload ( n, 0, iw(leType)  , 1 )
         call iload ( m, 3, iw(leType+n), 1 )

         nb        = n   + m
         nlocA     = n   + 1
         nName     = nb
         lrc       = lpi + m
         lgObj     = lbl

         nnH0      = max( iw(nnH)  , 1 )
         ngQP      = max( iw(nnH)  , iw(ngObj) )
         ngObj0    = max( iw(ngObj), 1   )

         nrhs      = 0
         lenrhs    = max( nrhs , 1   )
         nx0       = 0
         lenx0     = max( nx0  , 1   )

         cStart    = 'Cold'       ! Preempted by lvlStart

         startType = iw(lvlStart)
         call s3chkArgsQ
     &      ( INFO,
     &        cStart, m, n, neA, nName, nS,
     &        iw(ngObj), iObj, iw(nnH),
     &        iw(lindA), iw(llocA), rw(lbl), rw(lbu), cw(lNames),
     &        iw(lhs), rw(lpi), startType, Errors,
     &        iw, leniw, rw, lenrw )

!        Record n, m, neA and iObj for s5Defaults.

         iw( 15)    = n    ! copy of the number of columns
         iw( 16)    = m    ! copy of the number of rows
         iw( 17)    = neA  ! copy of the number of nonzeros in Acol
         iw(204)    = iObj ! position of the objective row in A

!        ---------------------------------------------------------------
!        Check options.
!        Open any files needed for this problem.
!        ---------------------------------------------------------------
         call s5Defaults
     &      ( m, n, iw(ngObj), iw(nnH),
     &        cw, lencw, iw, leniw, rw, lenrw )
         call s3printQ
     &      ( m, n, iw(ngObj), iw(nnH), startType, iw, leniw, rw,lenrw )
         call s1file
     &      ( OpenF , iw, leniw  )

!        ---------------------------------------------------------------
!        Compute the storage requirements for SQOPT  from the following
!        variables:
!           m    ,  n   , neA
!           lenR , maxS , maxR
!           ngObj, nnH
!        All are now known.
!        ---------------------------------------------------------------
         maxR    = iw( 52) ! max columns of R.
         maxS    = iw( 53) ! max # of superbasics

         nkx     = nb

         lenR    = maxR*(maxR + 1)/2 +  (maxS - maxR)
         iw( 28) = lenR         ! R(lenR) is the reduced Hessian factor

!        ---------------------------------------------------------------
!        Allocate the local arrays for snOptQ.
!        s5Map  maps snOptQ integer and double arrays.
!        s2BMap maps the arrays for the LU routines.
!        s2Mem  checks what space is available and prints any messages.
!        ---------------------------------------------------------------
         call s5Map
     &      ( m, n, nkx, iw(ngObj), iw(nnH),
     &        lenR, maxR, maxS,
     &        nextcw, nextiw, nextrw, iw, leniw )
         call s2Bmap
     &      ( m, n, neA, maxS,
     &        nextiw, nextrw, iw(maxiw), iw(maxrw), liwEst, lrwEst,
     &        iw, leniw )
         PrintMem = .true.        ! Print all messages in s2Mem
         call s2Mem
     &      ( inform, PrintMem, liwEst, lrwEst,
     &        nextcw, nextiw, nextrw,
     &        iw(maxcw), iw(maxiw), iw(maxrw), lencw, leniw, lenrw,
     &        mincw, miniw, minrw, iw )
         if (inform .ne. 0) then
            INFO = inform
            go to 800
         end if

!        Define the row and column ordering for J.
!        SQOPT  uses natural order throughout, so kx = kxN.

         iw(247) = nkx     ! dimension of kx and its inverse, kxN
         lkx     = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
         iw(252) = lkx     ! jN = kxN(j ) => col j of Jcol is variable jN

         call s1perm( n, iw(lkx) )
         call s1perm( m, iw(lkx+n) )

!        ------------------------------------------------------------------
!        Solve the problem.
!        ------------------------------------------------------------------
         call s5solve
     &      ( INFO,
     &        Solver, startType,
     &        sqlog, sqHx, qpHx, GotR,
     &        m, n, nb, nnH0, iw(nnH), nName, ngQP, ngObj0, iw(ngObj),
     &        iObj, objAdd, objQP, objTrue,
     &        nInf, sInf, nInfE, sInfE,
     &        neA, nlocA, iw(llocA), iw(lindA), rw(lAcol),
     &        neH, nlocH,     locH,     indH,       Hcol,
     &        rw(lbl), rw(lbu), rw(lgObj), cw(lNames),
     &        lenrhs, nrhs, rhs, lenx0, nx0, x0,
     &        iw(leType), iw(lhs), rw(lx), rw(lpi), rw(lrc), nS,
     &        cw, lencw, iw, leniw, rw, lenrw,
     &        cw, lencw, iw, leniw, rw, lenrw )

!        Print times for all clocks (if lvlTim > 0).

         call s1time
     &      ( 0, 2, iw, leniw, rw, lenrw )

      end do
!     ==================================================================
!     End of loop through SPECS file.
!     ==================================================================

  800 call snWRAP( INFO, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine sqmps1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine qpHx
     &   ( nnH, x, Hx, Status, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     lencu, leniu, lenru, nnH, Status, iu(leniu)
      double precision
     &     Hx(nnH), x(nnH), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     This version of qpHx is a dummy routine used for solving
!     LP's with the stand-alone version of sqopt.
!     It should never be called by SQOPT.
!
!     Warn the user (via the standard output) that it has been called.
!     ==================================================================
      external
     &     s1outpt
      integer
     &     nOut, s1outpt
!     ------------------------------------------------------------------
      nOut = s1outpt( )
      if (nOut .gt. 0) write(nOut, 9000)
      return

 9000 format(/ ' XXX dummy qpHx has been called in error.')

      end ! subroutine qpHx
