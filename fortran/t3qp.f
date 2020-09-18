*     ------------------------------------------------------------------
*     File t3qp.f  (Unix version)
*     Illustrates using SNOPT on a quadratic program.
*
*     15 May 1998: First   version.
*     17 Jul 2005: Current version.
*     ------------------------------------------------------------------
      program
     &     t3main

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1000,
     &     maxne  = 3000,
     &     nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne)  , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)   , rc(maxn+maxm)

*     SNOPT workspace---------------------------------------------------
      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lencw = 500, leniw = 10000, lenrw = 20000 )
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
     &     Errors, iSpecs, iPrint, iSumm, i1, i2, INFO,
     &     iObj, itnlim, m, mincw, miniw, minrw, n, ne, nInf,
     &     nnCon, nnJac, nOut, nnObj, nS
      double precision
     &     Obj, ObjAdd, sInf
      character
     &     lfile*20
      external
     &     t3obj, dumcon

*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*     nOut    is an output file used here by t3qp.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't3qp.spc'
         open( iSpecs, file=lfile, status='OLD',     err=900 )

         lfile = 't3qp.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=900 )
      end if

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 800
      end if

*     Set up the data structure for the sparse Jacobian.
*     Assign dummy values for the nonlinear elements.

      Errors = 0

      call t3data
     &   ( maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

*     ------------------------------------------------------------------
*     Specify any options not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      itnlim = 250
      i1     =   0
      i2     =   0
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
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
      call snOpt
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     dumcon, t3obj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '

      if (INFO .gt. 30) then
         write(nOut, *) ' '
         write(nOut, *) 'STOPPING because of error condition'
      else
         write(nOut, *) ' '
         write(nOut, *) 't3qp finished.'
         write(nOut, *) 'Input  errors =', Errors
         write(nOut, *) 'snOptB INFO   =', INFO
         write(nOut, *) 'nInf          =', nInf
         write(nOut, *) 'sInf          =', sInf
         if (iObj .gt. 0) then
            write(nOut, *)
     &                  'Obj           =', ObjAdd + x(n+iObj) + Obj
         else
            write(nOut, *)
     &                  'Obj           =', ObjAdd + Obj
         end if
      end if

  800 if ( byname ) then
         close(iSpecs)
         close(iPrint)
      end if

      stop

*     ------------------------------------------------------------------
*     File Error.
*     ------------------------------------------------------------------
  900 write(nOut, 4000) 'Error while opening file', lfile
      stop

 4000 format(/  a, 2x, a  )

      end ! program t3main

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t3obj
     &   ( mode, n, x, f, g, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, n, nState, lencu, leniu, lenru, iu(leniu)
      character
     &     cu(lencu)*8
      double precision
     &     f, x(n), g(n), ru(lenru)
*     ------------------------------------------------------------------
*     This is funobj for problem t3qp.
*     f(x) = 1/2 x'Qx,   g(x) = Qx.
*     ------------------------------------------------------------------
      integer
     &     i,j
      double precision
     &     ddot, xj
      external
     &     ddot
*     ------------------------------------------------------------------
      double precision   zero,           half
      parameter        ( zero = 0.0d+0,  half = 0.5d+0 )
      double precision   two,            four
      parameter        ( two  = 2.0d+0,  four = 4.0d+0 )
      integer            maxn
      parameter        ( maxn = 10 )
*     ------------------------------------------------------------------
      double precision   Q(maxn,maxn)
      save               Q
*     ------------------------------------------------------------------
      if (nState .eq. 1) then

*        Define Q on the first entry.
*        Here we assume n = 3.

         Q(1,1) = four
         Q(1,2) = two
         Q(1,3) = two
         Q(2,2) = four
         Q(2,3) = zero
         Q(3,3) = two

*        Make Q symmetric.

         do j = 1, n-1
            do i = 2, n
               Q(i,j) = Q(j,i)
            end do
         end do
      end if

*     Compute f and g on all entries.
*     We first compute g = Qx, then f = 1/2 x'g.
*     In Fortran it is best to run down the columns of Q.

      do i = 1, n
         g(i)   = zero
      end do

      do j = 1, n
         xj     = x(j)
         do i = 1, n
            g(i)  = g(i)  +  Q(i,j) * xj
         end do
      end do

      f   = half  *  ddot  ( n, x, 1, g, 1 )

      end ! subroutine t3obj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dumcon
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     Problem t3qp.
*     No nonlinear constraints.
*     ==================================================================

*     Relax

      end ! subroutine dummy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t3data
     &   ( maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      character
     &     Prob*8
      integer
     &     Errors, iObj, maxm, maxn, maxne, m, n, ne,
     &     nnCon, nnObj, nnJac, indA(maxne), hs(maxn+maxm),
     &     locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

*     ------------------------------------------------------------------
*     Define the problem.
*     (1) Compute l, u, and A so that the constraints are ranges of the
*         form  l <= Ax <= u.
*         Store l and u in bl(n+1:n+m) and bu(n+1:n+m).
*
*     (2) Set up the constants ObjAdd and iObj so that the explicit
*         objective is
*             ObjAdd + (row iObj of A)'*x + f(x)
*
*     On entry,
*     maxm, maxn, maxne are upper limits on m, n, ne.
*
*     On exit,
*     Errors  is 0 if there is enough storage, 1 otherwise.
*     m       is the number of nonlinear and linear constraints.
*     n       is the number of variables.
*     ne      is the number of nonzeros in Acol(*).
*     nnCon   is the number of nonlinear constraints (they come first).
*     nnObj   is the number of nonlinear objective variables.
*     nnJac   is the number of nonlinear Jacobian variables.
*     Acol    is the constraint matrix (Jacobian), stored column-wise.
*     indA    is the list of row indices for each nonzero in Acol(*).
*     locA    is a set of pointers to the beginning of each column of a.
*     bl      is the lower bounds on x and s.
*     bu      is the upper bounds on x and s.
*     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
*     x (1:n) is a set of initial values for x.
*     pi(1:m) is a set of initial values for the dual variables pi.
*     ------------------------------------------------------------------
      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision   infBnd
      parameter         (infBnd   =  1.0d+21)
*     ------------------------------------------------------------------

*     Name the Problem.

      Prob = 't3qp....'

      n    =  3
      m    =  2
      ne   =  6               ! n*m for a dense A.

*     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

      nnCon   = 0
      nnJac   = 0
      nnObj   = n

*     The linear objective is row 2 of a.

      iObj    = 2
      ObjAdd  = 0.0d+0

*     Column 1.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) =  2

      Acol( 1) =  1.0d+0
      Acol( 2) = -8.0d+0

*     Column 2.

      locA( 2) =  3

      indA( 3) =  1
      indA( 4) =  2

      Acol( 3) =  1.0d+0
      Acol( 4) = -6.0d+0

*     Column 3.

      locA( 3) =  5

      indA( 5) =  1
      indA( 6) =  2

      Acol( 5) =  2.0d+0
      Acol( 6) = -4.0d+0

*     Don't forget to finish off  locA(n+1) = ne+1
*     This is crucial.

      locA(n+1) =  ne + 1

*     ------------------------------------------------------------------
*     Set the upper and lower bound on the row.
*     Make the linear row a free row.
*     ------------------------------------------------------------------
      bl(n+1) = -infBnd
      bu(n+1) =  3.0d+0

      bl(n+2) =  -infBnd
      bu(n+2) =   infBnd

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j) =  0.0d+0
         bu(j) =  infBnd
      end do

*     ------------------------------------------------------------------
*     Initialize x, hs and pi.
*     Set the initial value and status of each variable.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         x (j) = 0.0d+0
         hs(j) = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine t3data

