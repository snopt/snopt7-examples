!     ------------------------------------------------------------------
!     File t4manneb.f
!     Illustrates using snOptB on a constrained problem.
!     It generates the problem called MANNE on Pages 98-108 of the
!     MINOS 5.1 User's Guide, then asks snOptB to solve it.
!
!     16 May 1998: First   version.
!     11 Apr 2005: Current version.
!     ------------------------------------------------------------------
      program
     &     t4main

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1500,
     &     maxne  = 3000,
     &     nName  = 1 )

      character
     &     ProbNm*8, Names(nName)*8
      integer
     &     indA(maxne), hs(maxm+maxn), locA(maxn+1)
      double precision
     &     Acol(maxne) , bl(maxm+maxn), bu(maxm+maxn),
     &     x(maxm+maxn), pi(maxm)     , rc(maxm+maxn)
      integer
     &     lenrw, leniw, lencw
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  lenrw = 10000)
      double precision   rw(lenrw)
      parameter          (  leniw =  5000)
      integer            iw(leniw)
      parameter          (  lencw =    500)
      character          cw(lencw)*8
!     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      external
     &     ManCon, ManObj
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim,
     &     m, maxS, mincw, miniw, minrw, n, ne, nInf, nnCon, nnJac,
     &     nnObj, nOut, nS, nT
      double precision
     &     Obj, ObjAdd, sInf
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by sntest.

      iSpecs = 4   ! equivalenced to t4manne.spc
      iPrint = 9   ! equivalenced to t4manne.out
      iSumm  = 6   ! summary file goes to standard output...
      nOut   = 6   ! ... as do messages from this program.

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't4manneb.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't4manneb.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file.  This must include "Nonlinear constraints  T"
!     for some integer T.
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

!     ------------------------------------------------------------------
!     The following assignment allows access to nnCon,
!     which defines T, the number of nonlinear constraints.
!     It is specified at runtime in the SPECS file.
!     ------------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Nonlinear constraints', nnCon, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      nT = nnCon
      if (nT .le. 1  .or.  nT .gt. maxm/2) then
         write(nOut, *) 'Invalid  nT  specified:', nT
         go to 990
      end if

!     Write nT into the problem name.

      write(probnm, '(i8)') nT
      if      (nT .lt.  1000) then
         probnm(1:5) = 'Manne'
      else if (nT .lt. 10000) then
         probnm(1:4) = 'Mann'
      else
         probnm(1:3) = 'Man'
      end if

      write(nOut, *) 'Problem MANNE.    T =', nT

!     ------------------------------------------------------------------
!     Generate an nT-period problem.
!     ------------------------------------------------------------------
      call t4data
     &   ( nT, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      if (Errors .ge. 1) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  nT =', nT
         go to 990
      end if

!     ------------------------------------------------------------------
!     Specify options that were not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      maxS   =  2*nT
      itnlim = 40*nT
      i1     = 0
      i2     = 0
      call snSeti
     &   ( 'Superbasics Limit ', maxS  , i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
!     iObj   = 0 means there is no linear objective row in Acol(*).
!     ObjAdd = 0.0 means there is no constant to be added to the
!            objective.
!     hs     need not be set if a basis file is to be input.
!            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
!            The values are used by the Crash procedure s2crsh
!            to choose an initial basis B.
!            If hs(j) = 0 or 1, column j is eligible for B.
!            If hs(j) = 2, column j is initially superbasic (not in B).
!            If hs(j) = 3, column j is eligible for B and is given
!                          preference over columns with hs(j) = 0 or 1.
!            If hs(j) = 4 or 5, column j is initially nonbasic.
!
!     SNOPT is called with iw and rw used for USER workspace.
!     This allows access to SNOPT variables in  ManCon and ManObj.
!     ------------------------------------------------------------------
      iObj   = 0
      ObjAdd = 0.0d+0

      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, ProbNm,
     &     ManCon, ManObj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '
      write(nOut, *) 't4manneb finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *)
     &               'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *)
     &               'Obj           =', ObjAdd + Obj
      end if

      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  990 stop

 4000 format(/  a, 2x, a  )

!     end of main program to test subroutine snOpt
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t4data
     &   ( nT, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     nT, maxm, maxn, maxne, Errors, m, n, ne, nnCon, nnObj, nnJac,
     &     indA(maxne), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

!     ------------------------------------------------------------------
!     t4data  generates data for the test problem t4manne
!     (called problem MANNE in the SNOPT 4.3 User's Guide).
!     The constraints take the form
!              f(x) + A*x - s = 0,
!     where the Jacobian for f(x) + Ax is stored in Acol(*), and any
!     terms coming from f(x) are in the TOP LEFT-HAND CORNER of Acol(*),
!     with dimensions  nnCon x nnJac.
!     Note that the right-hand side is zero.
!     s is a set of slack variables whose bounds contain any constants
!     that might have formed a right-hand side.
!
!     The objective function is
!             f(x) + c'x
!     where c would be row iObj of A (but there is no such row in
!     this example).  f(x) involves only the FIRST nnObj variables.
!
!     On entry,
!     nT      is T, the number of time periods.
!     maxm, maxn, maxne are upper limits on m, n, ne.
!
!     On exit,
!     Errors    is 0 if there is enough storage, 1 otherwise.
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
!     hs(1:n) is a set of initial states for each x (0,1,2,3,4,5).
!     x(1:n)  is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!
!     09 Jul 1992: No need to initialize x and hs for the slacks.
!     15 Oct 1993: pi is now an output parameter.  (Should have been
!                  all along.)
!     ------------------------------------------------------------------
      integer
     &     i, jI, jC, jM, jY, k, T
      double precision
     &     dummy, growth, bplus, bminus, scale
!     ------------------------------------------------------------------
      double precision   zero,              one
      parameter        ( zero   = 0.0d+0,   one    = 1.0d+0  )
      parameter        ( dummy  = 0.1d+0,   growth = .03d+0,
     &                   bplus  = 1.0d+20,  bminus = - bplus )
!     ------------------------------------------------------------------
!     nT defines the dimension of the problem.

      m      = nT*2
      n      = nT*3
      nnCon  = nT
      nnObj  = nT*2
      nnJac  = nT
      ne     = nT*6 - 1
      T      = nT

!     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

!     Generate columns for Capital (Kt, t = 1 to nT).
!     The first nT rows are nonlinear, and the next nT are linear.
!     The Jacobian is an nT x nT diagonal.
!     We generate the sparsity pattern here.
!     We put in dummy numerical values of 0.1 for the gradients.
!     Real values for the gradients are computed by funcon.

      ne     = 0
      do k = 1, nT

!        There is one Jacobian nonzero per column.

         ne       = ne + 1
         locA(k)  = ne
         indA(ne) = k
         Acol(ne) = dummy

!        The linear constraints form an upper bidiagonal pattern.

         if (k .gt. 1) then
            ne       = ne + 1
            indA(ne) = nT + k - 1
            Acol(ne) = one
         end if

         ne       = ne + 1
         indA(ne) = nT + k
         Acol(ne) = - one
      end do

!     The last nonzero is special.

      Acol(ne)  = growth

!     Generate columns for Consumption (Ct for t = 1 to nT).
!     They form -I in the first nT rows.
!     jC and jI are base indices for the Ct and It variables.

      jC    = nT
      jI    = nT*2

      do k = 1, nT
         ne         = ne + 1
         locA(jC+k) = ne
         indA(ne)   = k
         Acol(ne)   = - one
      end do

!     Generate columns for Investment (It for t = 1 to nT).
!     They form -I in the first nT rows and -I in the last nT rows.

      do k = 1, nT
         ne         = ne + 1
         locA(jI+k) = ne
         indA(ne)   = k
         Acol(ne)   = - one
         ne         = ne + 1
         Acol(ne)   = - one
         indA(ne)   = nT + k
      end do

!     locA(*) has one extra element.

      locA(n+1) = ne + 1

!     Set lower and upper bounds for Kt, Ct, It.
!     Also initial values and initial states for all variables.
!     The Jacobian variables are the most important.
!     We make them all superbasic.
!     The others are ok nonbasic.
!     For test purposes, we want the initial x to be infeasible
!     with respect to the linear constraints.
!     Try setting the last Kapital too high.


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

!-->     hs(   k) = 2
         hs(   k) = 0

         hs(jC+k) = 0
         hs(jI+k) = 0

         if (k .eq. nT) then
            x(k)  = 1.0d+3
            hs(k) = 2
         end if
      end do

!     The first Capital is fixed.
!     The last three Investments are bounded.
!     Fudge them to be the normal ones for T = 10.

      scale       = T / 10.0d+0
      bu(1)       = bl(1)
      x(1)        = bl(1)
      hs(1)       = 0
      bu(jI+nT-2) = 0.112d+0 * scale
      bu(jI+nT-1) = 0.114d+0 * scale
      bu(jI+nT  ) = 0.116d+0 * scale

!     Set bounds on the slacks.
!     The nT nonlinear (Money)    rows are >=.
!     The nT    linear (Capacity) rows are <=.
!     We no longer need to set initial values and states for slacks
!     (assuming SNOPT does a cold start).

      jM     = n
      jY     = n + nT

      do k = 1, nT
         bl(jM+k) = zero
         bu(jM+k) = bplus
         bl(jY+k) = bminus
         bu(jY+k) = zero

!-       x (jM+k) = zero
!-       x (jY+k) = zero
!-       hs(jM+k) = 0
!-       hs(jY+k) = 0
      end do

!     The last Money and Capacity rows have a Range.

      bu(jM+nT) =   10.0d+0
      bl(jY+nT) = - 20.0d+0

!     Initialize pi.
!     SNOPT requires only pi(1:nnCon) to be initialized.
!     We initialize all of pi just in case.

      do i = 1, nT
         pi(i)    = - one
         pi(nT+i) = + one
      end do

      end ! subroutine t4data

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ManObj
     &   ( mode, nnObj, x, fObj, gObj, State,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     mode, nnObj, State, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), rw(lenrw)
      character
     &     cw(lencw)*8

!     ------------------------------------------------------------------
!     This is funobj for problem t4manne.
!
!     The data bt(*) is computed by funcon  on its first entry.
!
!     For test purposes, we look at    Derivative level
!     and sometimes pretend that we don't know the first
!     three elements of the gradient.
!
!     12-Nov-93: Changed from Maximize to Minimize to steer around bugs.
!     13-Jan-95: Reverted to Maximize.
!     ==================================================================
      intrinsic
     &     log
      integer
     &     j, lvlDer, nT
      logical
     &     gknown
      double precision
     &     xcon
!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   b, at(365), bt(365)
      common    /manne / b, at,      bt
!     ------------------------------------------------------------------
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level

      gknown = lvlDer .eq. 1  .or.  lvlDer .eq. 3
      nT     = nnObj/2
      fObj   = zero

      do j = 1, nT
         xcon = x(nT+j)
         fObj    = fObj  +  bt(j) * log(xcon)
!Min     fObj    = fObj  -  bt(j) * log(xcon)
         if (mode .gt. 0) then
            gObj(j) = zero
            if (gknown                ) gObj(nT+j) = + bt(j) / xcon
!           if (gknown  .or.  j .gt. 3) gObj(nT+j) = + bt(j) / xcon
!Min        if (gknown  .or.  j .gt. 3) gObj(nT+j) = - bt(j) / xcon
         end if
      end do

      end ! subroutine ManObj (funobj for t4manne)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ManCon
     &   ( mode, nnCon, nnJac, njac, x, fCon, gCon, State,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, njac, State,
     &     lencw, leniw, lenrw, iw(leniw)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(njac), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     This is funcon for problem t4manne.
!
!     For test purposes, we look at    Derivative level
!     and sometimes pretend that we don't know the first
!     three elements of the gradient.
!     ==================================================================
      logical
     &     gknown
      integer
     &     iPrint, j, lvlDer, nT
      double precision
     &     a, beta, gfac, grow, xk0, xKap, xc0, xi0
!     ------------------------------------------------------------------
      double precision    one
      parameter         ( one = 1.0d+0 )
      double precision    b, at(365), bt(365)
      common     /manne / b, at,      bt
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level

      gknown = lvlDer .ge. 2
      nT     = nnJac

!     ---------------------------------------
!     First entry.  Define b, at(*) and bt(*)
!     for this and all subsequent entries.
!     ---------------------------------------
      if (State .eq. 1) then
         grow  = 0.03d+0
         beta  = 0.95d+0
         xk0   = 3.0d+0
         xc0   = 0.95d+0
         xi0   = 0.05d+0
         b     = 0.25d+0
         if (iPrint .gt. 0) write(iPrint, 1000) nt, b

         a     = (xc0 + xi0) / xk0**b
         gfac  = (one + grow)**(one - b)
         at(1) = a*gfac
         bt(1) = beta

         do j  = 2, nT
            at(j) = at(j-1)*gfac
            bt(j) = bt(j-1)*beta
         end do

         bt(nT) = bt(nT) / (one - beta)
      end if

!     -------------
!     Normal entry.
!     -------------
      do j = 1, nT
         xkap    = x(j)
         fCon(j) = at(j) * xkap**b
         if (mode .gt. 0) then
            if (gknown                ) gCon(j) = b*fCon(j) / xkap
!           if (gknown  .or.  j .gt. 3) gCon(j) = b*fCon(j) / xkap
         end if
      end do

!     ------------
!     Final entry.
!     ------------
      if (State .ge. 2) then
         if (iPrint .gt. 0) write(iPrint, 2000) (fCon(j), j = 1, nT)
      end if

      return

 1000 format(// ' This is problem  t4manne.   nT =', i4, '   b =', f8.3)
 2000 format(// ' Final nonlinear function values' / (5f12.5))

      end ! subroutine ManCon (funcon for t4manne)
