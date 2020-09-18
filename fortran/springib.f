!     ------------------------------------------------------------------
!     File springib.f for snOptB
!
!     This is a main program to generate an optimal control problem
!     of arbitrary size and solve it by calling SNOPT as a subroutine.
!
!     The problem is identical to springb.f except that the constraints
!     are infeasible.
!
!     The problem size depends on a parameter T.  There are
!     2T constraints and 3T + 2 variables, as well as bounds
!     on the variables.  The first T constraints are quadratic in
!     T + 1 variables, and the objective function is quadratic in
!     T + 1 other variables.
!
!     The control problem models a spring, mass and damper system.
!     It is of the form
!
!   --------------------------------------------------------------------
!   | minimize    1/2 sum x(t)**2   (t = 0 to T)                       |
!   |                                                                  |
!   | subject to                                                       |
!   |     y(t+1)  =  y(t)  -  0.01 y(t)**2  -  0.004 x(t)  +  0.2 u(t) |
!   |                                                                  |
!   |     x(t+1)  =  x(t)  +  0.2  y(t),                               |
!   |                                                                  |
!   |     y(t)   >=  -1,     -0.2  <=  u(t)  <=  0.2,                  |
!   |                                                                  |
!   |                (all for t = 0 to T-1)                            |
!   | and                                                              |
!   |     y(0)    =   0,      y(T)  =  0,       x(0) = 10.             |
!   --------------------------------------------------------------------
!
!     For large enough T (e.g. T >= 90), the optimal objective value
!     is about 1186.382.
!
!     This model with T = 100 was used as test problem 5.11 in
!     B. A. Murtagh and M. A. Saunders (1982), A projected Lagrangian
!     algorithm and its implementation for sparse nonlinear constraints,
!     Mathematical Programming Study 16, 84--117.
!
!     14 Nov 1994: First version of spring.f, derived from manne.f.
!     26 Sep 1997: Updated for SNOPT 5.3.
!     03 Sep 2015: Fixed some typos in the comments.
!     ------------------------------------------------------------------
      program            springib

      implicit
     &     none
      integer
     &     maxT, maxm, maxn, maxne, nName
      parameter
     &   ( maxT   = 2000,
     &     maxm   = 2*maxT,
     &     maxn   = 3*maxT + 2,
     &     maxne  = 7*maxT,
     &     nName  = 1 )

      character
     &     ProbNm*8, Names(nName)*8
      integer
     &     indA(maxne) , hs(maxm+maxn), locA(maxn+1)
      double precision
     &     Acol(maxne) , bl(maxm+maxn), bu(maxm+maxn),
     &     x(maxm+maxn), pi(maxm)     , rc(maxm+maxn)

!     ------------------------------------------------------------------
!     USER workspace (none required)

      integer
     &     lenru, leniu, lencu
      parameter
     &   (  lenru = 1,
     &      leniu = 1,
     &      lencu = 1)
      integer
     &     iu(leniu)
      double precision
     &     ru(lenru)
      character
     &     cu(lencu)*8
!     ------------------------------------------------------------------
!     SNOPT workspace

      integer
     &     lenrw, leniw, lencw
      parameter
     &     (lenrw = 50000,
     &      leniw = 50000,
     &      lencw =   500)
      integer
     &     iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8
!     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      integer
     &     Errors, INFO, iObj, iSpecs, iPrint, iSumm, m, mincw, miniw,
     &     minrw, n, ne, nnCon, nnObj, nnJac, nOut, nS, ninf, T
      double precision
     &     ObjAdd, sinf, Obj
      external
     &     SprCon, SprObj
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by the main program.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'springib.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'springib.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

!     ------------------------------------------------------------------
!     Set options to their default values.
!     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Read a Specs file.  This must include "Problem number  T"
!     for some integer T.
!     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         stop
      end if

!     ------------------------------------------------------------------
!     The following call fetches T, which defines the number of
!     nonlinear constraints.
!     It is specified at runtime in the SPECS file.
!     ------------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Problem number', T, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (T .le. 1  .or.  T .gt. maxm/2) then
         write(nOut, *) 'Invalid no. of Nonlinear constraints:', T
         stop
      end if

!     Write T into the problem name.

      write(ProbNm, '(i8)') T
      if      (T .lt.   100) then
         ProbNm(1:6) = 'Spring'
      else if (T .lt.  1000) then
         ProbNm(1:5) = 'Sprng'
      else if (T .lt. 10000) then
         ProbNm(1:4) = 'Spri'
      else
         ProbNm(1:3) = 'Spr'
      end if

      write(nOut, *) 'Problem springib.   T =', T

!     ------------------------------------------------------------------
!     Generate an T-period problem.
!     ------------------------------------------------------------------
      call spdata
     &   ( T, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      if (Errors .gt. 0) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  Nonlinear constraints =', T
         stop
      end if

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
!     iObj   = 0 means there is no linear objective row in Acol(*).
!     Objadd = 0.0 means there is no constant to be added to the
!            objective.
!     hs     need not be set if a basis file is to be input.
!            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
!            The values are used by the Crash procedure
!            to choose an initial basis B.
!            If hs(j) = 0 or 1, column j is eligible for B.
!            If hs(j) = 2, column j is initially superbasic (not in B).
!            If hs(j) = 3, column j is eligible for B and is given
!                          preference over columns with hs(j) = 0 or 1.
!            If hs(j) = 4 or 5, column j is initially nonbasic.
!     ------------------------------------------------------------------
      iObj   = 0
      ObjAdd = 0.0d+0

      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, ProbNm,
     &     SprCon, SprObj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, ninf, sinf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 900
      end if

      write(nOut, *) ' '
      write(nOut, *) 'springib (fg) finished.'
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
      if (INFO .ge. 30) go to 900
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

      end ! program springib

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine spdata
     &   ( T, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     Errors, m, maxm, maxn, maxne, n, ne, nnCon, nnObj, nnJac, T,
     &     indA(maxne), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne)    , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

!     ------------------------------------------------------------------
!     spdata generates data for the "Spring" optimal control problem.
!     The constraints take the form
!              c(x) + A*x - s = 0,
!     where the Jacobian for c(x) + Ax is stored in Acol(*), and any
!     terms coming from c(x) are in the TOP LEFT-HAND CORNER of Acol(*),
!     with dimensions  nnCon x nnJac.
!     Note that the right-hand side is zero.
!     s is a set of slack variables whose bounds contain any constants
!     that might have formed a right-hand side.
!
!     The objective function is
!             f(x) + d'x
!     where d would be row iobj of A (but there is no such row in
!     this example).  f(x) involves only the FIRST nnObj variables.
!
!     On entry,
!     T       is the number of time periods.
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
!     Acol    is the constraint matrix (Jacobian), stored column-wise.
!     indA    is the list of row indices for each nonzero in Acol(*).
!     locA    is a set of pointers to the beginning of each column of a.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n) is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!
!     14 Nov 1994: First version of spdata.
!     ------------------------------------------------------------------
      integer
     &     i, j, k, nb
!     ------------------------------------------------------------------
      double precision zero,              one
      parameter      ( zero   = 0.0d+0,   one    = 1.0d+0 )
      double precision bplus,             dummy
      parameter      ( bplus  = 1.0d+20,  dummy  = 0.111111d+0 )
!     ------------------------------------------------------------------
!     T defines the dimension of the problem.

      m      = T*2
      n      = T*3 + 2  ! y(0:T), x(0:T) and u(0:T-1)
      nb     = n   + m
      nnCon  = T
      nnObj  = T*2 + 2  ! y(0:T) and x(0:T)
      nnJac  = T   + 1  ! y(0:T)
      ne     = T*7

!     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

!     ------------------------------------------------------------------
!     Generate columns for y(t), t = 0 to T.
!     The first T rows are nonlinear, and the next T are linear.
!     The Jacobian is T x (T+1) upper bidiagonal.
!     We generate the sparsity pattern here.
!     We put in dummy numerical values for the nonlinear gradients.
!     The true non-constant values are computed by funcon.
!     ------------------------------------------------------------------
      j      = 0   ! counts the variables
      ne     = 0   ! counts the Jacobian and linear constraint entries

      do   k = 0, T
         j       =   j  + 1
         locA(j) =   ne + 1
         bl(j)   = - one
         bu(j)   =   bplus
         x (j)   = - one
         hs(j)   =   0      ! Make the y(t) nonbasic.

!        There are two Jacobian nonzeros per column,
!        except in the first and last column.

         if (k .gt. 0) then    !  Aij = 1
            ne       =   ne + 1
            indA(ne) =   k
            Acol(ne) =   one
         end if

         if (k .lt. T) then    !  Aij = .02y - 1  (nonlinear)
            ne       =   ne + 1
            indA(ne) =   k  + 1
            Acol(ne) =   dummy
         end if

!        Below the Jacobian the linear constraints are diagonal.

         if (k .lt. T) then
            ne       =   ne + 1
            indA(ne) =   T  + k + 1
            Acol(ne) = - 0.2d+0
         end if
      end do

!     ------------------------------------------------------------------
!     Generate columns for x(t), t = 0 to T.
!     They form 0.004*I in the first T rows,
!     and an upper-bidiagonal in the last T rows.
!     ------------------------------------------------------------------
      do       k = 0, T
         j       =   j  + 1
         locA(j) =   ne + 1
         bl(j)   = - bplus
         bu(j)   =   bplus
         x (j)   =   zero
         hs(j)   =   3     ! Make the x(t) basic.

!        Part of 0.004*I.

         if (k .lt. T) then
            ne       =   ne + 1
            indA(ne) =   k  + 1
            Acol(ne) =   0.004d+0
         end if

!        The bidiagonal parts have two entries
!        except in the first and last columns.

         if (k .gt. 0) then    !  Aij = 1
            ne       =   ne + 1
            indA(ne) =   T  + k
            Acol(ne) =   one
         end if

         if (k .lt. T) then    !  Aij = - 1
            ne       =   ne + 1
            indA(ne) =   T  + k + 1
            Acol(ne) = - one
         end if
      end do

!     ------------------------------------------------------------------
!     Generate columns for u(t), t = 0 to T-1.
!     They form -0.2I in the first T rows.
!     ------------------------------------------------------------------
      do k = 0, T - 1
         j       =   j  + 1
         locA(j) =   ne + 1
         bl(j)   = - 0.2d+0
         bu(j)   =   0.2d+0
         x (j)   =   zero
         hs(j)   =   3     ! Make the u(t) basic.

         ne       =   ne + 1
         indA(ne) =   k  + 1
         Acol(ne) = - 0.2d+0
      end do

!     locA(*) has one extra element.
!     Some of the variables are fixed.

      locA(n+1) = ne + 1

      bl(1)   = zero      ! y(0) = 0
      bu(1)   = zero

      bl(T+1) = zero      ! y(T) = 0
      bu(T+1) = zero

      bl(T+2) = 10.0d+0   ! x(0) = 10
      bu(T+2) = 10.0d+0

!     Added 11/9/96 to make springib infeasible

      bl(2*T+3) = -10.0d+0   ! u(0) = -10
      bu(2*T+3) = -10.0d+0

!     ------------------------------------------------------------------
!     Set bounds on the slacks.
!     We don't need to set initial values and states for slacks
!     (assuming SNOPT does a cold start).
!     ------------------------------------------------------------------
      do      j = n + 1, nb
         bl(j)  = zero
         bu(j)  = zero
      end do

!     Initialize pi.
!     SNOPT requires only pi(1:nnCon) to be initialized.
!     We initialize all of pi just in case.

      do       i = 1, T
         pi(i)   =   zero
         pi(T+i) =   zero
      end do

      end ! subroutine spdata

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SprObj
     &   ( mode, n, x, f, g, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     lencu, leniu, lenru, mode, n, nState, iu(leniu)
      double precision
     &     f, x(n), g(n), ru(lenru)
      character
     &     cu(lencu)*8
!     ==================================================================
!     This is funobj for problem Springib  (an optimal control problem).
!     ==================================================================
      integer
     &     jy, jx, k, T
      double precision
     &     u
!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------
      T      = (n - 2)/2
      f      = zero
      jy     = 0
      jx     = T + 1

      do k = 0, T
         jy    = jy + 1
         jx    = jx + 1
         u     = x(jx)
         f     = f  +  u**2
         g(jy) = zero
         g(jx) = u
      end do

      f = f / 2.0d+0

      end ! subroutine SprObj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SprCon
     &   ( mode, m, n, njac, x, f, g, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     lencu, leniu, lenru, mode, m, n, njac, nState, iu(leniu)
      double precision
     &     x(n), f(m), g(njac), ru(lenru)
      character
     &     cu(lencu)*8

!     ------------------------------------------------------------------
!     This is funcon for problem Springib  (Optimal Control Problem).
!     The Jacobian is upper bidiagonal,
!     and only the diagonal terms are nonlinear.
!     The constant 1's in the Jacobian are not regenerated here.
!     ------------------------------------------------------------------
      integer
     &     i, jg, jy, T
      double precision
     &     yt, ytp1
!     ------------------------------------------------------------------
      double precision   one
      parameter        ( one = 1.0d+0 )
!     ------------------------------------------------------------------
      T     = n - 1
      jy    =     0    ! Counts y(t) variables
      jg    =   - 1    ! Counts nonlinear Jacobian elements

      do i = 1, T
         jy     = jy + 1
         jg     = jg + 2
         yt     = x(jy)
         ytp1   = x(jy + 1)
         f(i)   = 0.01d+0 * yt**2  +  (ytp1  -  yt)
         g(jg)  = 0.02d+0 * yt               -  one
!--      g(jg+1)= one      ! Constant element set by spdata.
      end do

      end ! subroutine SprCon

