!     ------------------------------------------------------------------
!     File meyer3b.f
!     Illustrates using SNOPT on the unconstrained problem meyer3.
!
!     15 May 1998: First   version.
!     04 Sep 2006: Current version.
!     ------------------------------------------------------------------
      program
     &     meyer3main

      implicit
     &     none

      integer
     &     maxm, maxn, maxne, nName
      parameter
     &   ( maxm   =    1,
     &     maxn   = 1000,
     &     maxne  = 3000,
     &     nName  = 1 )

      integer
     &     indA(maxne), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)
      character
     &     Prob*8, Names(nName)*8

*     ------------------------------------------------------------------
*     SNOPT workspace

      integer
     &     lenrw, leniw, lencw
      parameter
     &     (lenrw = 10000,
     &      leniw =  5000,
     &      lencw =   500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character
     &     cw(lencw)*8
*     ------------------------------------------------------------------
*     USER workspace (used by calvar routines)

      integer
     &     lenru, leniu, lencu
      parameter
     &     (lenru = 16,
     &      leniu =  1,
     &      lencu =  1)
      integer
     &     iu(leniu)
      double precision
     &     ru(lenru)
      character
     &     cu(lencu)*8
!     ------------------------------------------------------------------
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim,
     &     m, mincw, miniw, minrw, n, ne, nInf, nnCon,
     &     nnJac, nnObj, nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
      logical
     &     byname
      character
     &     lfile*20
      external
     &     dummy, meyer3obj
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!
!     nOut    is an output file used here by meyer3.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'meyer3b.spc'
         open( iSpecs, file=lfile, status='OLD',     err=900 )

         lfile = 'meyer3b.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=900 )
      end if

!     ------------------------------------------------------------------
!     First,  snInit MUST be called to initialize optional parameters
!     to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file (Optional).
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 800
      end if

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      Errors = 0

      call meyer3Data
     &   ( maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      itnlim = 2500
      i1     =    0
      i2     =    0
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Hessian full      ',         i1, i2, Errors,
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
      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     dummy, meyer3obj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '
      write(nOut, *) 'Meyer3 finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      write(nOut, *) 'Obj           =', ObjAdd + Obj

  800 if ( byname ) then
         close(iSpecs)
         close(iPrint)
      end if

      stop

!     ------------------------------------------------------------------
!     File Error.
!     ------------------------------------------------------------------
  900 write(nOut, 4000) 'Error while opening file', lfile
      stop

 4000 format(/  a, 2x, a  )

      end ! program meyer3main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine meyer3obj
     &   ( mode, n, x, fObj, gObj, State,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, n, State, lencu, leniu, lenru, iu(leniu)
      double precision
     &     fObj, x(n), gObj(n)
      character
     &     cu(lencu)*8
      double precision
     &     ru(lenru)

!     ------------------------------------------------------------------
!     This is funobj for problem meyer33
!     ------------------------------------------------------------------
      integer
     &     i, mData
      double precision
     &     e1, g1, g2, g3, sum, s1, s2, ti, x1, x2, x3
*     ------------------------------------------------------------------
      double precision   zero,          two
      parameter         (zero = 0.0d+0, two = 2.0d+0)
!     ------------------------------------------------------------------

*     ---------------------------------------
*     First entry.  Define b, at(*) and bt(*)
*     for this and all subsequent entries.
*     ---------------------------------------
      mData = 16

      if (State .eq. 1) then
         ru( 1) = 34780.0d+0
         ru( 2) = 28610.0d+0
         ru( 3) = 23650.0d+0
         ru( 4) = 19630.0d+0
         ru( 5) = 16370.0d+0
         ru( 6) = 13720.0d+0
         ru( 7) = 11540.0d+0
         ru( 8) =  9744.0d+0
         ru( 9) =  8261.0d+0
         ru(10) =  7030.0d+0
         ru(11) =  6005.0d+0
         ru(12) =  5147.0d+0
         ru(13) =  4427.0d+0
         ru(14) =  3820.0d+0
         ru(15) =  3307.0d+0
         ru(16) =  2872.0d+0
      end if

      x1  = x(1)
      x2  = x(2)
      x3  = x(3)

      sum = zero
      g1  = zero
      g2  = zero
      g3  = zero

      ti  = 45
      do i = 1, mData
         ti  = ti   + 5
         e1  = exp(x2/(ti + x3))
         s1  = x1*e1 - ru(i)
         s2  = ti + x3
         sum = sum +  s1*s1

         if (mode .gt. 0) then
            g1 = g1 + two*s1*e1
            g2 = g2 + two*s1*e1*x1/s2
            g3 = g3 - two*s1*e1*x1*x2/(s2*s2)
         end if
      end do

      fObj    = sum
      gObj(1) = g1
      gObj(2) = g2
      gObj(3) = g3


      end ! subroutine meyer3obj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dummy
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnCon, nnJac, neJac, State, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem meyer33.
!     No nonlinear constraints.
!     ==================================================================

!     Relax

      end ! subroutine dummy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine meyer3Data
     &   ( maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none

      integer
     &     Errors, maxm, maxn, maxne, m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character
     &     Prob*8

!     ------------------------------------------------------------------
!     Define the problem.
!     (1) Compute l, u, and A so that the constraints are ranges of the
!         form  l <= Ax <= u.
!         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
!
!     (2) Set up the constants ObjAdd and iObj so that the explicit
!         objective is
!             ObjAdd + (row iObj of A)'*x + f(x)
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
!     a       is the constraint matrix (Jacobian), stored column-wise.
!     indA      is the list of row indices for each nonzero in Acol(*).
!     locA      is a set of pointers to the beginning of each column of a.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n) is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!     ------------------------------------------------------------------
      integer
     &     i, j
      double precision
     &     bplus, zero, ten
      parameter         (bplus   =  1.0d+21)
      parameter         (zero    =  0.0d+0,  ten  = 10.0d+0)
!     ------------------------------------------------------------------
!     Name the Problem.

      Prob = 'meyer3'

      n     =  3
      m     =  1
      ne    =  1                ! Dummy row

!     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

      nnCon   = 0
      nnJac   = 0
      nnObj   = n

!     ==================================================================
!     SNOPT requires at least one row in a, so we have to
!     provide a dummy linear constraint.
!
!     Define a (sparse) free row that has just one element.
!     Set this element to be zero.
!     Note that locA(j) = locA(j+1) for all columns j that have no entries.
!     ------------------------------------------------------------------
      locA(1) =  1
      indA(1) =  1
      Acol(1) =  zero

      do j = 2, n
         locA(j) = 2
      end do
      locA(n+1) =  ne + 1

!     Make the linear row a free row.

      bl(n+1) =  -bplus
      bu(n+1) =   bplus

!     ------------------------------------------------------------------
!     Set the upper and lower bounds on the variables
!     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -bplus
         bu(j) =  bplus
      end do

!     No linear objective term for this problem.

      iObj    = 0
      ObjAdd  = zero

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1) =     0.0056096365d+0 ! Optimal solution
      x(2) =  6181.3463462519d+0
      x(3) =   345.2236346230d+0

      x(1) =    2.0d-2 ! Standard start.
      x(2) = 4000.0d+0
      x(3) =  250.0d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = zero
      end do

      end ! subroutine meyer3Data
