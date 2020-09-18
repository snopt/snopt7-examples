*     ------------------------------------------------------------------
*     File t4mannec.f for snOptC
*     Illustrates using snOptC on a constrained problem.
*     It generates the problem called MANNE on Pages 98-108 of the
*     MINOS 5.1 User's Guide, then asks snOptC to solve it.
*
*     20 Jun 2005: First version of t4mannec.f.
*     22 Jun 2005: Current version.
*     ------------------------------------------------------------------
      program
     &     t4mannec

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1500,
     &       maxne  = 3000,
     &       nName  = 1 )

      character
     &     Names(nName)*8
      integer
     &     indA(maxne) , hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenru, leniu, lencu, lenrw, leniw, lencw
*     ------------------------------------------------------------------
*     SNOPT workspace

      parameter          (  lenrw = 10000)
      double precision   rw(lenrw)
      parameter          (  leniw = 5000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character*8        cw(lencw)
*     ------------------------------------------------------------------
*     User workspace

      parameter          (  lenru = 1)
      double precision   ru(lenru)
      parameter          (  leniu = 1)
      integer            iu(leniu)
      parameter          (  lencu = 1)
      character          cu(lencu)*8
*     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20, ProbNm*8
      external
     &     usrfun
      integer
     &     i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim, m, maxS,
     &     mincw, miniw, minrw, n, ne, nInf, nnCon, nnJac, nnObj, nOut,
     &     nS, nT
      double precision
     &     Obj, ObjAdd, sInf
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*     nOut    is an output file used here by snmain.

      iSpecs = 4
      iPrint = 15
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't4mannec.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't4mannec.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file.  This must include "Nonlinear constraints T"
*     for some integer T.
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101) then
         go to 990
      end if

!     ------------------------------------------------------------------
!     1. Solve manne with default options.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '1. Manne with default options'
      write(nOut, *) '-------------------------------------------'

*     The following call fetches nnCon, which defines T, the number of
*     nonlinear constraints.  It is specified at runtime in the
*     SPECS file

      call snGeti
     &   ( 'Nonlinear constraints', nnCon,
     &     INFO, cw, lencw, iw, leniw, rw, lenrw)

      nT = nnCon

      if (nT .le. 1 .or. nT .gt. maxm/2) then
         write(nout,*) 'Invalid number of Nonlinear constraints:', nT
         stop
      end if

*     Write nT into the problem name.

      write(ProbNm, '(i8)') nT
      if (nT .lt. 1000) then
         ProbNm(1:5) = 'Manne'
      else if (nT .lt. 10000) then
         ProbNm(1:4) = 'Mann'
      else
         ProbNm(1:3) = 'Man'
      end if
      write(nout,*) 'Problem MANNE. T = ', nT

*     Set up the data structure for the sparse Jacobian.
*     Assign dummy values for the nonlinear elements.

      call t4ManneDat
     &   ( nT, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      if (INFO .ge. 1) then
	write(nout, *) 'Not enough storage to generate a problem ',
     &		      'with nT =', nT
	go to 990
      end if

*     ------------------------------------------------------------------
*     Specify any options not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      maxS   =  2*nT
      itnlim = 40*nT
      i1     =  0
      i2     =  0

      call snSeti
     &	( ' Superbasics Limit ', maxS  , i1, i2, INFO,
     &	  cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
*     hs     need not be set if a basis file is to be input.
*            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
*            The values are used by the Crash procedure m2crsh
*            to choose an initial basis B.
*            If hs(j) = 0 or 1, column j is eligible for B.
*            If hs(j) = 2, column j is initially superbasic (not in B).
*            If hs(j) = 3, column j is eligible for B and is given
*                          preference over columns with hs(j) = 0 or 1.
*            If hs(j) = 4 or 5, column j is initially nonbasic.
*     ------------------------------------------------------------------
      call snOptC
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, ProbNm,
     &     usrfun,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 't4mannec finished.'
      write(nOut, *) 'snOptC INFO  =', INFO
      write(nOut, *) 'nInf  =', nInf
      write(nOut, *) 'sInf  =', sInf
      write(nOut, *) 'Obj   =', Obj
      if (INFO .ge. 10) go to 910

!     ------------------------------------------------------------------
!     2. Solve manne with some scaled linear constraints and variables.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '2. Manne with scaling'
      write(nOut, *) '-------------------------------------------'

*     Set up the data structure for the sparse Jacobian again.

      call t4ManneDat
     &   ( nT, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      if (INFO .ge. 1) then
	write(nout, *) 'Not enough storage to generate a problem ',
     &		      'with nT =', nT
	go to 990
      end if

      call snSeti
     &   ( 'Scale option', 1, i1, i2, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptC
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, ProbNm,
     &     usrfun,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 't4mannec (scaled) finished.'
      write(nOut, *) 'snOptC INFO  =', INFO
      write(nOut, *) 'nInf  =', nInf
      write(nOut, *) 'sInf  =', sInf
      write(nOut, *) 'Obj   =', Obj
      if (INFO .ge. 10) go to 910

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

      end ! program snoptc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun
     &   ( mode, nnObj, nnCon, nnJac, nnL, neJac,
     &     x, fObj, gObj, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnObj, nnCon, nnJac, nnL, neJac, State,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnL)
      double precision
     &     fObj, gObj(nnObj)
      double precision
     &     fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     This is usrfun for problem t4mannec.
*
*     The data bt(*) is computed by funcon on its first entry.
*
*     For test purposes, we look at    Derivative level
*     and sometimes pretend that we don't know the first
*     three elements of the gradient.
*     ==================================================================
      logical
     &     needf, needg
      integer
     &     nT, j, jKap, jCon, nOut
      double precision
     &     a, axKap, beta, gfac, grow, xCon, xKap, xKap0, xCon0,
     &     xInt0
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      double precision   b, at(365), bt(365)
      common    /manne / b, at,      bt
*     ------------------------------------------------------------------
      nT   = nnCon              ! Size of the problem.
      nOut = 15

*     ---------------------------------------
*     First entry.  Define b, at(*) and bt(*)
*     for this and all subsequent entries.
*     ---------------------------------------
      if (State .eq. 1) then
         grow  = 0.03d+0
         beta  = 0.95d+0
         xKap0 = 3.0d+0
         xCon0 = 0.95d+0
         xInt0 = 0.05d+0
         b     = 0.25d+0
         if (nOut .gt. 0) write(nOut, 1000) nt, b

         a     = (xCon0 + xInt0) / xKap0**b
         gfac  = (one + grow)**(one - b)
         at(1) = a*gfac
         bt(1) = beta

         do j  = 2, nT
            at(j) = at(j-1)*gfac
            bt(j) = bt(j-1)*beta
         end do

         bt(nT) = bt(nT) / (one - beta)
      end if

      needf = mode .eq. 0  .or.  mode .eq. 2
      needg = mode .eq. 1  .or.  mode .eq. 2

      if (needf) fObj = zero

      jKap   = 0                  ! counts the K(t) variables.
      jCon   = nT                 ! counts the C(t) variables.

      do j = 1, nT
         jKap  = jKap + 1
         jCon  = jCon + 1


         xCon  = x(jCon)
         xKap  = x(jKap)
         axKap = at(jKap)*xKap**b

         if (needf) then
            fCon(jKap) = axKap
            fObj       = fObj + bt(jKap)*log(xCon)
         end if

         if (needg) then
            gCon(jKap) = b*axKap/xKap

            gObj(jKap) = zero
            gObj(jCon) = bt(jKap)/xCon
         end if
      end do

*     ------------
*     Final entry.
*     ------------
      if (State .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000) (fCon(j), j = 1, nnCon)
      end if

      return

 1000 format(// ' This is problem  t4manne.   nT =', i4, '   b =', f8.3)
 2000 format(// ' Final nonlinear function values' / (5f12.5))

      end ! subroutine usrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t4ManneDat
     &   ( nT, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, INFO, m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne) , hs(maxn+maxm), locA(maxn+1), nT
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

*     ------------------------------------------------------------------
*     t4data  generates data for the test problem t4manne
*     (called problem MANNE in the SNOPT 4.3 User's Guide).
*     The constraints take the form
*              f(x) + A*x - s = 0,
*     where the Jacobian for f(x) + Ax is stored in Acol(*), and any
*     terms coming from f(x) are in the TOP LEFT-HAND CORNER of Acol(*),
*     with dimensions  nnCon x nnJac.
*     Note that the right-hand side is zero.
*     s is a set of slack variables whose bounds contain any constants
*     that might have formed a right-hand side.
*
*     The objective function is
*             f(x) + c'x
*     where c would be row iObj of A (but there is no such row in
*     this example).  f(x) involves only the FIRST nnObj variables.
*
*     On entry,
*     nT      is T, the number of time periods.
*     maxm, maxn, maxne are upper limits on m, n, ne.
*
*     On exit,
*     INFO    is 0 if there is enough storage, 1 otherwise.
*     m       is the number of nonlinear and linear constraints.
*     n       is the number of variables.
*     ne      is the number of nonzeros in Acol(*).
*     nnCon   is the number of nonlinear constraints (they come first).
*     nnObj   is the number of nonlinear objective variables.
*     nnJac   is the number of nonlinear Jacobian variables.
*     a       is the constraint matrix (Jacobian), stored column-wise.
*     indA      is the list of row indices for each nonzero in Acol(*).
*     locA      is a set of pointers to the beginning of each column of a.
*     bl      is the lower bounds on x and s.
*     bu      is the upper bounds on x and s.
*     hs(1:n) is a set of initial states for each x (0,1,2,3,4,5).
*     x(1:n)  is a set of initial values for x.
*     pi(1:m) is a set of initial values for the dual variables pi.
*
*     09 Jul 1992: No need to initialize x and hs for the slacks.
*     15 Oct 1993: pi is now an output parameter.  (Should have been
*                  all along.)
*     ------------------------------------------------------------------
      integer
     &     i, jI, jC, jM, jY, k, T
      double precision
     &     dummy, growth, bplus, bminus, scale
*     ------------------------------------------------------------------
      double precision   zero,              one
      parameter        ( zero   = 0.0d+0,   one    = 1.0d+0  )
      parameter        ( dummy  = 0.1d+0,   growth = .03d+0,
     &                   bplus  = 1.0d+20,  bminus = - bplus )
*     ------------------------------------------------------------------
*     nT defines the dimension of the problem.

      m      = nT*2
      n      = nT*3
      nnCon  = nT
      nnObj  = nT*2
      nnJac  = nT
      ne     = nT*6 - 1
      T      = nT

*     No linear component of the objective function.
      iObj = 0
      ObjAdd = zero

*     Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm ) INFO = 1
      if (n     .gt. maxn ) INFO = 1
      if (ne    .gt. maxne) INFO = 1
      if (INFO  .gt.   0  ) return

*     Generate columns for Capital (Kt, t = 1 to nT).
*     The first nT rows are nonlinear, and the next nT are linear.
*     The Jacobian is an nT x nT diagonal.
*     We generate the sparsity pattern here.
*     We put in dummy numerical values of 0.1 for the gradients.
*     Real values for the gradients are computed by funcon.

      ne     = 0

      do k = 1, nT

*        There is one Jacobian nonzero per column.

         ne       = ne + 1
         locA(k)  = ne
         indA(ne) = k
         Acol(ne) = dummy

*        The linear constraints form an upper bidiagonal pattern.

         if (k .gt. 1) then
            ne       = ne + 1
            indA(ne) = nT + k - 1
            Acol(ne) = one
         end if

         ne       = ne + 1
         indA(ne) = nT + k
         Acol(ne) = - one
      end do

      Acol(ne)  = growth   ! The last nonzero is special.

*     Generate columns for Consumption (Ct for t = 1 to nT).
*     They form -I in the first nT rows.
*     jC and jI are base indices for the Ct and It variables.

      jC    = nT
      jI    = nT*2

      do k = 1, nT
         ne         = ne + 1
         locA(jC+k) = ne
         indA(ne)   = k
         Acol(ne)   = - one
      end do

*     Generate columns for Investment (It for t = 1 to nT).
*     They form -I in the first nT rows and -I in the last nT rows.

      do k = 1, nT
         ne         = ne + 1
         locA(jI+k) = ne
         indA(ne)   = k
         Acol(ne)   = - one
         ne         = ne + 1
         Acol(ne)   = - one
         indA(ne)   = nT + k
      end do

*     locA(*) has one extra element.

      locA(n+1) = ne + 1

*     Set lower and upper bounds for Kt, Ct, It.
*     Also initial values and initial states for all variables.
*     The Jacobian variables are the most important.
*     We make them all superbasic.
*     The others are ok nonbasic.
*     For test purposes, we want the initial x to be infeasible
*     with respect to the linear constraints.
*     Try setting the last Kapital too high.


      do k = 1, nT
         bl(   k) = 3.05d+0
         bu(   k) = bplus
         bl(jC+k) = 0.95d+0
         bu(jC+k) = bplus
         bl(jI+k) = 0.05d+0
         bu(jI+k) = bplus

         x(   k)  = 3.0d+0 + (k - 1)/10.0d+0
         x(jC+k)  = bl(jC+k)
         x(jI+k)  = bl(jI+k)

*-->     hs(   k) = 2
         hs(   k) = 0

         hs(jC+k) = 0
         hs(jI+k) = 0

         if (k .eq. nT) then
            x(k)  = 1.0d+3
            hs(k) = 2
         end if
      end do

*     The first Capital is fixed.
*     The last three Investments are bounded.
*     Fudge them to be the normal ones for T = 10.

      scale       = T / 10.0d+0
      bu(1)       = bl(1)
      x(1)        = bl(1)
      hs(1)       = 0
      bu(jI+nT-2) = 0.112d+0 * scale
      bu(jI+nT-1) = 0.114d+0 * scale
      bu(jI+nT  ) = 0.116d+0 * scale

*     Set bounds on the slacks.
*     The nT nonlinear (Money)    rows are >=.
*     The nT    linear (Capacity) rows are <=.
*     We no longer need to set initial values and states for slacks
*     (assuming SNOPT does a cold start).

      jM     = n
      jY     = n + nT

      do k = 1, nT
         bl(jM+k) = zero
         bu(jM+k) = bplus
         bl(jY+k) = bminus
         bu(jY+k) = zero

*-       x (jM+k) = zero
*-       x (jY+k) = zero
*-       hs(jM+k) = 0
*-       hs(jY+k) = 0
      end do

*     The last Money and Capacity rows have a Range.

      bu(jM+nT) =   10.0d+0
      bl(jY+nT) = - 20.0d+0

*     Initialize pi.
*     SNOPT requires only pi(1:nnCon) to be initialized.
*     We initialize all of pi just in case.

      do i = 1, nT
         pi(i)    = - one
         pi(nT+i) = + one
      end do

      end ! subroutine t4ManneDat
