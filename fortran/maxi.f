!     ------------------------------------------------------------------
!     File maxi.f
!
!     maxi     mxdata   funcon   funmax
!
!     This is a program to generate and solve optimization problems
!     of the form
!
!         Given P, maximize d, the minimum diameter of P nonoverlapping
!         circles that will fit inside a unit square.
!     i.e.
!         maximize   d2
!         such that  (x(p) - x(q))**2  +  (y(p) - y(q))**2  >=  d2,
!                                           p = 1..P, q = 1..p-1,
!                    0 <= x(p), y(p) <= 1,  p = 1..P.
!         where d2 = d**2.
!
!     maxi    is the main program.
!     mxdata  is the problem generator.
!     funcon  defines the constraints and Jacobian.
!     funmax  does the work for funcon.
!     snOptB   gets called to solve the problem.
!
!
!     28 Aug 1992: First version, based on AMPL formulation sent to
!                  MAS by David Gay.  MINOS 5.4 (Jul 1992) had trouble.
!     30 May 1994: Modified to call SNOPT by PEG.
!     27 Nov 1996: Modified for SNOPT versions 4.3 and higher.
!     26 Sep 1997: Modified for SNOPT versions 5.3 and higher.
!     14 Mar 2004: Modified for SNOPT versions 7.1 and higher.
!     ------------------------------------------------------------------
      program
     &     maxi
      implicit
     &     none
      integer
     &     maxP, maxm, maxn, maxne, nname
      parameter
     &   ( maxP   = 50,
     &     maxm   = (maxP - 1)*maxP / 2 + 1,
     &     maxn   = 2*maxP + 1,
     &     maxne  = 2*(maxP - 1)*maxP + maxm,
     &     nname  = 1 )
      character
     &     Prob*8, Names(nname)*8
      integer
     &     indA(maxne), hs(maxm+maxn), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxm+maxn), bu(maxm+maxn),
     &     x(maxm+maxn), pi(maxm), rc(maxm+maxn)

      integer
     &     leniw, lenrw, lencw
      parameter
     &     (  leniw = 100000 +  50*maxP)
      parameter
     &     (  lenrw = 100000 + 100*maxP)
      parameter
     &     (  lencw = 500              )
      character
     &     cw(lencw)*8
      integer
     &     iw(leniw)
      double precision
     &     rw(lenrw)
!     ------------------------------------------------------------------
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm,
     &     ltime, m, mincw, miniw, minrw, n, ne, nInf, nnCon,
     &     nnJac, nnObj, nOut, nProb, nS, P
      logical
     &     byname
      double precision
     &     Obj, ObjAdd, sInf
      character
     &     charP*4, lfile*20
      external
     &     funcon, funobj
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'maxi.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'maxi.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file.
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         stop
      end if

!     ------------------------------------------------------------------
!     Set P to be the number of circles.  We assume the Specs file
!     does this via the line
!        Problem number    P
!     The following call fetches nProb, which defines P.
!     ------------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Problem number', nProb, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      P     = nProb
      write(charP, '(i4)') P

!     Give names to the Problem, Objective, Rhs, Ranges and Bounds.

      Prob  = 'maxi' // charP

      call mxdata
     &   ( P, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, x, pi, hs )

      if (INFO .ge. 1) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  P =', P
         stop
      end if

!     ------------------------------------------------------------------
!     Specify options that were not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      i1     = 0
      i2     = 0
      ltime  = 2
      call snSeti
     &   ( 'Timing level      ', ltime, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
!     iObj   = m   specifies the linear objective row in Acol(*).
!     Objadd = 0.0 means there is no constant to be added to the
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
!     ------------------------------------------------------------------
      iObj   = m
      ObjAdd = 0.0d0

      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     funcon, funobj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 900
      end if

      write(nOut, *) ' '
      write(nOut, *) 'maxi finished.'
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
      if (INFO .gt. 30) go to 900
      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
      stop

 4000 format(/  a, 2x, a  )

      end ! program maxi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mxdata
     &   ( P, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, x, pi, hs )

      implicit
     &     none
      integer
     &     P, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     indA(maxne), hs(maxm+maxn), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxm+maxn), bu(maxm+maxn), x(maxm+maxn),
     &     pi(maxm)

!     ------------------------------------------------------------------
!     mxdata  generates data for the test problems maxi.
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
!     where c is row iobj of A (but there is no f(x) in this example).
!
!     On entry,
!     P       is the number of time circles.
!     maxm, maxn, maxne are upper limits on m, n, ne.
!
!     On exit,
!     INFO    is 0 if there is enough storage, 1 otherwise.
!     m       is the number of nonlinear and linear constraints.
!     n       is the number of variables.
!     ne      is the number of nonzeros in Acol(*).
!     nnCon   is the number of nonlinear constraints (they come first).
!     nnObj   is the number of nonlinear objective variables.
!     nnJac   is the number of nonlinear Jacobian variables.
!     a       is the constraint matrix (Jacobian), stored column-wise.
!     indA    is the list of row indices for each nonzero in Acol(*).
!     locA    is a set of pointers to the beginning of each column of a.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     x (1:n) is a set of initial values for x.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!
!     28 Aug 1992: First version, based on t4data.
!     ------------------------------------------------------------------
      integer
     &     i, id, iz, j, k, l, nb
      double precision
     &     d2min, realP
      double precision
     &     zero, one, dummy, bplus, bminus
      parameter
     &   ( zero   = 0.0d+0,   one    = 1.0d+0,
     &     dummy  = 0.1d+0,
     &     bplus  = 1.0d+20,  bminus = - bplus )
!     ------------------------------------------------------------------
!     P defines the dimension of the problem.

      m      = (P - 1)*P / 2 + 1
      n      = 2*P + 1
      nb     = n + m
      nnCon  = m - 1
      nnObj  = 0
      nnJac  = n - 1
      ne     = 2*(P - 1)*P + m

!     Check if there is enough storage.

      INFO   = 0
      if (m      .gt. maxm ) INFO = 1
      if (n      .gt. maxn ) INFO = 1
      if (ne     .gt. maxne) INFO = 1
      if (INFO   .gt.   0  ) return

!     ------------------------------------------------------------------
!     The variables are P pairs (x(p),y(p)) and d2.
!     The x's come first, then the y's.
!     For P = 5, the sparsity pattern of the Jacobian is
!        a z       a z       -1
!        b   z     b   z     -1
!          b z       b z     -1
!        c     z   c     z   -1
!          c   z     c   z   -1
!            c z       c z   -1
!        d       z d       z -1
!          d     z   d     z -1
!            d   z     d   z -1
!              d z       d z -1
!                             1
!     where in general, a, b, c, d,... indicate P-1 diagonal matrices,
!     and there is a set of consecutive z's in each column.
!     The first P columns have the same structure as the next P.
!     The last column corresponds to the linear variable d2,
!     and the last row corresponds to the objective "maximize d2".
!     We put in dummy numerical values of 0.1 for a, b, c, d and z.
!     Real values for the gradients are computed by funcon.
!     ------------------------------------------------------------------

!     Generate the structure of the first P columns of the Jacobian.

      ne     = 0   ! total number of Jacobian elements
      iz     = 0   ! row index of each z element in current column
      id     = 0   ! row index of first diagonal element a or b or c...
      realP  = P

      do k = 1, P
         locA(k) = ne + 1

!        Column k has k-1 "z" entries.

         do  l = 1, k-1
            ne = ne + 1
            iz = iz + 1
            indA(ne) = iz
            Acol(ne) = dummy
         end do

!        Column k has P-k entries from the diagonal matrices.

         id    = iz + k
         do  l = k, P-1
            ne = ne + 1
            indA(ne) = id
            Acol(ne) = dummy
            id       = id + l
         end do

!        Set bounds and initial value and state.

         bl(k)  = zero
         bu(k)  = one
         x (k)  = k / (realP + one)
         hs(k)  = 0
      end do

!     The next P columns look exactly the same.

      call dcopy ( ne, Acol, 1, Acol(ne+1), 1 )
      call icopy ( ne, indA, 1, indA(ne+1), 1 )
      do k = 1, P
         locA(P+k) = locA(k) + ne
      end do
      call dcopy (  P, bl, 1, bl( P+1), 1 )
      call dcopy (  P, bu, 1, bu( P+1), 1 )
      call dcopy (  P, x , 1,  x( P+1), 1 )
      call icopy (  P, hs, 1, hs( P+1), 1 )
      ne      = 2*ne
      locA(n) = ne + 1

!     The last column is the linear variable d2.
!     It is -1 in all rows except the objective (the last row).
!     The last nonzero is 1 for the linear objective.

      do i  = 1, m
         ne =   ne + 1
         indA(ne) =   i
         Acol(ne) = - one
      end do

      Acol(ne)  = one
      locA(n+1) = ne + 1

!     ------------------------------------------------------------------
!     Set bounds for d2 (variable n) and for each slack.
!     The last (objective) slack is a free variable.
!
!     d2min is an obvious lower bound on d2, corresponding to
!     all circles centered along one edge of the square.
!     We initialize d2 = d2min
!     and set the lower bound to be half this
!     just so d2 is away from the lower bound.
!     ------------------------------------------------------------------
      d2min   = (2.0d+0 / (realP - one))**2
      x(n)    = d2min
      bl(n)   = d2min * 0.5d+0
      bu(n)   = bplus
      hs(n)   = 0

      do j = n+1, nb
         bl(j) = zero
         bu(j) = bplus
      end do

      bl(nb)   = bminus

      call dload ( nnCon, zero, pi, 1 )

      end ! subroutine mxdata

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj
     &   ( mode, n, x, f, g, State,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none

      integer
     &     mode, n, State, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     f, x(n), g(n), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     Nonlinear objective routines for maxi.
!     maxi does not have a nonlinear objective.
!     ==================================================================
      integer
     &     iPrint, iSumm
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      if (iPrint .gt. 0) write(iPrint, 9000)
      if (iSumm  .gt. 0) write(iSumm , 9000)
      mode   = -1
      return

 9000 format(/ ' XXX maxi should not call funobj.')

      end ! subroutine funobj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funcon
     &   ( mode, m, n, neJac, x, f, g, State,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none

      integer
     &     mode, m, n, neJac, State, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     x(n), f(m), g(neJac), rw(lenrw)
      character
     &     cw(lencw)*8

!     ------------------------------------------------------------------
!     This funcon defines nonlinear constraints and Jacobian entries
!     for problem maxi.  See main program "maxi" and generator "mxdata".
!
!     There are  P = n/2        circles
!           and  m = (P-1)*P/2  constraints of the form
!         (x(p) - x(q))**2  +  (y(p) - y(q))**2  >=  diameter.
!     For P = 4, the constraint functions are
!
!         (x2 - x1)**2   +   (y2 - y1)**2
!         (x3 - x1)**2   +   (y3 - y1)**2
!         (x3 - x2)**2   +   (y3 - y2)**2
!         (x4 - x1)**2   +   (y4 - y1)**2
!         (x4 - x2)**2   +   (y4 - y2)**2
!         (x4 - x3)**2   +   (y4 - y3)**2
!
!     and the first P columns of the Jacobian look like this:
!
!        -2(x2 - x1)   2(x2 - x1)
!        -2(x3 - x1)                2(x3 - x1)
!                     -2(x3 - x2)   2(x3 - x2)
!        -2(x4 - x1)                             2(x4 - x1)
!                     -2(x4 - x2)                2(x4 - x2)
!                                  -2(x4 - x3)   2(x4 - x3)
!
!     The next P columns look the same in the y variables,
!     so we make use of the auxiliary routine funmax.
!
!     28 Aug 1992: First version.  Boy, this is a lot tougher than AMPL.
!     ------------------------------------------------------------------
      integer
     &     neJacP, P
      double precision
     &     zero
      parameter
     &   ( zero   = 0.0d+0 )
!     ------------------------------------------------------------------
      P      = n/2
      neJacP = neJac/2
      call dload
     &   ( m, zero, f, 1 )
      call funmax
     &   ( m, P, neJacP, x     , f, g           )
      call funmax
     &   ( m, P, neJacP, x(P+1), f, g(neJacP+1) )

      end ! subroutine funcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funmax
     &   ( m, P, neJacP, x, f, g )

      implicit
     &     none
      integer
     &     m, P, neJacP
      double precision
     &     x(P), f(m), g(neJacP)

!     ------------------------------------------------------------------
!     funmax is an auxiliary routine for funcon for problem maxi.
!     Its structure parallels that of mxdata.
!     ------------------------------------------------------------------
      integer
     &     iz, k, l, ne
      double precision
     &     t, xk

      double precision
     &     two
      parameter
     &   ( two    = 2.0d+0 )
!     ------------------------------------------------------------------
      ne    = 0
      iz    = 0

      do  k = 1, P
         xk = x(k)

!        Column k has k-1 "z" entries.

         do l = 1, k-1
            ne    = ne + 1
            iz    = iz + 1
            t     = xk - x(l)
            f(iz) = f(iz) + t**2
            g(ne) = two * t
         end do

!        Column k has P-k entries from the diagonal matrices.

         do l = k, P-1
            ne    =   ne + 1
            t     =   x(l+1) - xk
            g(ne) = - two * t
         end do
      end do

      end ! subroutine funmax
