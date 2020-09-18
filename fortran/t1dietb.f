*     ------------------------------------------------------------------
*     File t1diet.f  (Unix version)
*     Illustrates the use of subroutine SNOPT on a linear program,
*
*     15 May 1998: First   version.
*     27 Oct 2003: Current version.
*     ------------------------------------------------------------------
      program            t1main

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1000,
     &     maxne  = 3000,
     &     nName  =    1 )
      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)

*     SNOPT workspace

      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lenrw = 20000)
      parameter
     &   ( leniw = 10000)
      parameter
     &   ( lencw =   500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character
     &     cw(lencw)*8
*     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      integer
     &     Errors, iSpecs, iPrint, iSumm, i1, i2, INFO, iObj,
     &     itnlim, m, mincw, miniw, minrw, n, ne, nInf, nnCon,
     &     nnJac, nOut, nnObj, nS
      double precision
     &     obj, ObjAdd, sInf
      external
     &     funcon, funobj       !dummy subroutines
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by t1diet.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't1dietb.spc'
         open( iSpecs, file=lfile, status='OLD',     err=900 )

         lfile = 't1dietb.out'
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

      call t1data
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

      write(nOut, *) ' '
      write(nOut, *) 't1dietb finished.'
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

      end ! subroutine t1diet

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnObj, nState, lencu, leniu, lenru, iu(leniu)
      character
     &     cu(lencu)*8
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)

*     ==================================================================
*     Dummy funobj for t1diet.
*     No nonlinear objective
*     ==================================================================

*     Relax

      end ! subroutine funobj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funcon
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
*     Dummy funcon for t1diet.
*     No nonlinear constraints.
*     ==================================================================

*     Relax

      end ! subroutine funcon

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t1data
     &   ( maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      character
     &     Prob*8
      integer
     &     maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne), hs(maxn+maxm), locA(maxn+1)
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
*             ObjAdd + (row iObj of A)'*x
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
*     Acol    is the constraint matrix (Jacobian) A, stored column-wise.
*     indA    is the list of row indices for each nonzero in Acol(*).
*     locA    is a set of pointers to the beginning of each column of A.
*     bl      is the lower bounds on x and s.
*     bu      is the upper bounds on x and s.
*     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
*     x (1:n) is a set of initial values for x.
*     pi(1:m) is a set of initial values for the dual variables pi.
*     ------------------------------------------------------------------
      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision   bplus
      parameter         (bplus   = 1.0d+21)
      double precision   zero,               one
      parameter         (zero    = 0.0d+0,   one    = 1.0d+0)
*     ------------------------------------------------------------------
*     Name the Problem.

      Prob = 'Diet LP.'

*     ------------------------------------------------------------------
*     This is the Diet problem of Chvatal, 1983.
*     Assign the constraint nonzeros to Acol, column by column.
*     indA(i) gives the row index of element Acol(i).
*     locA(j) gives the index in a of the start of column j.
*     ------------------------------------------------------------------
*           ( 110  205  160  160  420  260 )
*     A  =  (   4   32   13    8    4   14 )
*           (   2   12   54  285   22   80 )
*           (   3   24   13    9   20   19 ) ( = objective row c')

      n      =  6
      m      =  4               ! Includes the objective row
      ne     = 24               ! n*m for a dense A.

*     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.   0  ) return

      nnCon  =  0
      nnObj  =  0
      nnJac  =  0

*     Column 1.

      locA( 1) =  1

      indA( 1) =  1
      indA( 2) =  2
      indA( 3) =  3
      indA( 4) =  4

      Acol( 1) =   110.0d+0
      Acol( 2) =     4.0d+0
      Acol( 3) =     2.0d+0
      Acol( 4) =     3.0d+0

*     Column 2.

      locA( 2) =  5

      indA( 5) =  1
      indA( 6) =  2
      indA( 7) =  3
      indA( 8) =  4

      Acol( 5) =   205.0d+0
      Acol( 6) =    32.0d+0
      Acol( 7) =    12.0d+0
      Acol( 8) =    24.0d+0

*     Column 3.

      locA( 3) =  9

      indA( 9) =  1
      indA(10) =  2
      indA(11) =  3
      indA(12) =  4

      Acol( 9) =   160.0d+0
      Acol(10) =    13.0d+0
      Acol(11) =    54.0d+0
      Acol(12) =    13.0d+0

*     Column 4.

      locA( 4) = 13

      indA(13) =  1
      indA(14) =  2
      indA(15) =  3
      indA(16) =  4

      Acol(13) =   160.0d+0
      Acol(14) =     8.0d+0
      Acol(15) =   285.0d+0
      Acol(16) =     9.0d+0

*     Column 5.

      locA( 5) = 17

      indA(17) =  1
      indA(18) =  2
      indA(19) =  3
      indA(20) =  4

      Acol(17) =   420.0d+0
      Acol(18) =     4.0d+0
      Acol(19) =    22.0d+0
      Acol(20) =    20.0d+0

*     Column 6.

      locA(6)  = 21

      indA(21) =  1
      indA(22) =  2
      indA(23) =  3
      indA(24) =  4

      Acol(21) =   260.0d+0
      Acol(22) =    14.0d+0
      Acol(23) =    80.0d+0
      Acol(24) =    19.0d+0

*     Don't forget to finish off  locA(n+1) = ne+1
*     This is crucial.

      locA( 7) = 25

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on the variables
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j)  = zero
      end do

      bu(1)  =     4.0d+0
      bu(2)  =     3.0d+0
      bu(3)  =     2.0d+0
      bu(4)  =     8.0d+0
      bu(5)  =     2.0d+0
      bu(6)  =     2.0d+0

*     ------------------------------------------------------------------
*     Set the upper and lower bounds on  Ax.
*     The last row is free (i.e., infinite upper and lower bounds).
*     ------------------------------------------------------------------
      bl( 7) =  2000.0d+0
      bl( 8) =    55.0d+0
      bl( 9) =   800.0d+0
      bl(10) =  -bplus

      do i = 1, m
         j     = n + i
         bu(j) = bplus
      end do

      iObj   = 4
      ObjAdd = zero

*     ------------------------------------------------------------------
*     Initialize x, hs and pi.
*     Set the initial value and status of each variable.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         x(j)   = one
      end do

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = zero
      end do

      end ! subroutine t1data

