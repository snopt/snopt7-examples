*     ------------------------------------------------------------------
*     sntoyb.f implements the sample nonlinear problem defined in the
*     snOptB User's Guide.
*
*     snOptB is an interface from the SNOPT package.
*
*     31 Jul 1996: First   version.
*     19 Oct 2003: Updated for SNOPT 7
*     23 Apr 2007: Case with sparse Jacobian added.
*     ------------------------------------------------------------------
      program
     &     sntoyb

      implicit
     &     none

      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm   = 100,
     &       maxn   = 100,
     &       maxne  = 300,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim,
     &     lencw, leniw, lenrw, m, mincw, miniw, minrw, n, ne, nInf,
     &     nnCon, nnJac, nnObj, nOut, nS
      double precision
     &     infBnd, Obj, ObjAdd, sInf
*     ------------------------------------------------------------------
*     SNOPT workspace
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character          cw(lencw)*8
*     ------------------------------------------------------------------
      external
     &     funcon, funcon2, funobj
      logical
     &     byname
      character
     &     lfile*20
*     ------------------------------------------------------------------
*     Give a name to the Problem.

      Prob = 'Toy NLP '

*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by sntoy.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'sntoyb.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'sntoyb.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (optional).
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

*     ------------------------------------------------------------------
*     Define what we mean by an infinite bound.
*     ------------------------------------------------------------------
      Errors = 0
      infBnd = 1.0d+20
      call snSetr
     &   ( 'Infinite Bound', infBnd, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Set up the problem constraints and bounds.
*     Assign dummy values for the nonlinear Jacobian elements.
*     ------------------------------------------------------------------
      call toyDat
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     ObjAdd, infBnd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

*     ------------------------------------------------------------------
*     Specify options that were not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      Errors = 0

      itnlim = 2500
      i1     =    0
      i2     =    0
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

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'sntoyb (1) finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

*     ------------------------------------------------------------------
*     Solve the problem again. This time the Jacobian is stored as a
*     sparse matrix.  The objective function is the same.
*     ------------------------------------------------------------------
      call toyDat2
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     ObjAdd, infBnd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     funcon2, funobj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'sntoyb (2) finished again.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .gt. 30) go to 910

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

      end ! program snoptb

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funcon
     &   ( mode, nnCon, nnJac, neJac,
     &     x, fCon, gCon, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(nnCon,nnJac), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     Toy NLP problem from the SNOPT User's Guide.
*     ==================================================================
      integer
     &     nOut
*     ------------------------------------------------------------------
      nOut   = 9

      if (nState .eq. 1) then    ! First entry
         if (nOut .gt. 0) write(nOut, '(/a)') ' This is problem  Toy'
      end if

      if (mode .eq. 0  .or.  mode .eq. 2) then
         fCon( 1)  = x(1)**2 +  x(2)**2
         fCon( 2)  =            x(2)**4
      end if

      if (mode .ge. 1) then

*        Jacobian elements for column 1.

         gCon(1,1) = 2.0d+0*x(1) ! Jacobian elements for column 1
         gCon(2,1) = 0.0d+0      ! Can't be omitted

*        Jacobian elements for column 2.

         gCon(1,2) = 2.0d+0*x(2) ! Jacobian elements for column 2
         gCon(2,2) = 4.0d+0*x(2)**3
      end if

      if (nState .ge. 2) then    ! Last entry
         if (nOut .gt. 0) write(nOut, '(/a)') ' Finished problem  Toy'
      end if

      end ! subroutine funcon

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj
     &   ( mode, nnObj,
     &     x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu, lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     Toy NLP problem from the SNOPT User's Guide.
*     ==================================================================
      double precision
     &     sum
*     ------------------------------------------------------------------
      sum    = x(1) + x(2) + x(3)

      if (mode .eq. 0  .or.  mode .eq. 2) then
         fObj    = sum*sum
      end if

      if (mode .eq. 1  .or.  mode .eq. 2) then
         sum     = 2.0d+0*sum
         gObj(1) = sum
         gObj(2) = sum
         gObj(3) = sum
      end if

      end ! subroutine funobj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine toyDat
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     ObjAdd, infBnd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      character
     &     Prob*8
      integer
     &     maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, infBnd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

*     ------------------------------------------------------------------
*     toyDat generates data for the snOptB example problem described in
*     the SNOPT User's Guide.
*
*      minimize   (x(1)    +  x(2)    + x(3))**2  + 3*x(3)  +  5*x(4)
*      subject to  x(1)**2 +  x(2)**2 + x(3)                =   2
*                             x(2)**4             + x(4)    =   4
*                2*x(1)   + 4*x(2)                         >=   0.
*
*     The constraints must be supplied in the form
*        low <=  c(x) + Ax <= upp
*     where  low  and  upp  are m-vectors of constant lower and upper
*     bounds on the constraints,
*
*     The lower bounds low are stored in bl(n+1), bl(n+2), ..., bl(n+m).
*     The upper bounds upp are stored in bu(n+1), bu(n+2), ..., bu(n+m).
*
*     One of the rows of A may define the elements of a constant term
*     in the objective (see the definition of the objective function
*     below).  An objective row of A must be defined as being "free",
*     i.e., it must have lower and upper bounds of -Infty and +Infty,
*     where Infty is a large positive number (see the definitional of
*     the optional parameter "Infinite bound" in the User's Guide.
*
*     The Jacobian for c(x) + Ax is stored in Acol(*), and any
*     terms coming from c(x) are in the TOP LEFT-HAND CORNER of Acol(*),
*     with dimensions  nnCon by nnJac.
*
*     In toyDat, the Jacobian of c(x) is held as a dense matrix.
*
*     snOptB converts the constraints into the form
*              c(x) + A*x - s = 0,
*     where  s is a set of slack variables whose lower and upper bounds
*     are defined from the elements of bl and bu.  Note that the
*     right-hand side of this internal form of the constraints is zero.
*
*     The objective function is
*             f(x) + d'x
*     where d is row  iObj  of Acol. (iObj = 4 in this example.)
*     The function f(x) involves  only the FIRST nnObj variables.
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
*
*     24 Dec 1997: First version of toyDat.
*     24 Apr 2007: Comments updated.
*     ------------------------------------------------------------------
      integer
     &     i, j
*     ------------------------------------------------------------------
*     Give a name to the Problem.

      Prob   = 'Toy NLP '

      ne     = 10
      n      =  4
      m      =  4

      nnCon  =  2
      nnJac  =  2
      nnObj  =  3

*     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.     0) return

*     -------------------------------------
*     Set up the list of row indices in indA.
*     -------------------------------------
*     Column  1
*     Nonlinear elements in rows (1, 2)  first.

      locA( 1) =  1

      indA( 1) =  1
      Acol( 1) =  0.0d+0

      indA( 2) =  2
      Acol( 2) =  0.0d+0

*     Linear element in row 3 next.

      indA( 3) =  3
      Acol( 3) =  2.0d+0

*     Column 2.
*     Nonlinear elements in rows (1, 2).

      locA( 2) =  4

      indA( 4) =  1
      Acol( 4) =  0.0d+0

      indA( 5) =  2
      Acol( 5) =  0.0d+0

*     Linear element in row 3.

      indA( 6) =  3
      Acol( 6) =  4.0d+0

*     Column 3.
*     Linear element in row 1.

      locA( 3) =  7

      indA( 7) =  1
      Acol( 7) =  1.0d+0

*     Objective row element in row 4.

      indA( 8) =  4
      Acol( 8) =  3.0d+0

*     Column 4.
*     Linear element in row 2

      locA( 4) =  9

      indA( 9) =  2
      Acol( 9) =  1.0d+0

*     Objective row element in row 4

      indA(10) =  4
      Acol(10) =  5.0d+0

*     Don't forget to finish off  locA.
*     This is crucial.

      locA(5) =  ne + 1

*     ------------------------------------------------------------------
*     Constraint ranges
*     ------------------------------------------------------------------
*     Nonlinear constraints first.

      bl(n+1) =  2.0d+0
      bu(n+1) =  2.0d+0

      bl(n+2) =  4.0d+0
      bu(n+2) =  4.0d+0

*     Followed by the linear constraints.

      bl(n+3) =  0.0d+0
      bu(n+3) =  infBnd

*     The linear objective term is row 4.

      iObj    = 4

*     The objective row is a free row.

      bl(n+4) = -infBnd
      bu(n+4) =  infBnd

      ObjAdd  = 0.0d+0

*     ------------------------------------------------------------------
*     Variable ranges
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -infBnd
         bu(j) =  infBnd
      end do

      bl(3) =  0.0d+0
      bl(4) =  0.0d+0

*     ------------------------------------------------------------------
*     Initialize x, hs and pi.
*     Set the initial value and status of each variable.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         x(j)  = 0.0d+0
      end do

      do j = 1, n
         hs(j) = 0
      end do

      do i = 1, m
         pi(i) = 0.0d+0
      end do

      end ! subroutine toyDat

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine toyDat2
     &   ( Prob, maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     ObjAdd, infBnd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      character
     &     Prob*8
      integer
     &     maxm, maxn, maxne, Errors,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj,
     &     indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, infBnd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)

*     ------------------------------------------------------------------
*     toyDat2 generates data for the snOptB example problem described in
*     the SNOPT User's Guide.
*
*     In toyDat2, the Jacobian of c(x) is held as a sparse matrix.
*
*     24 Apr 2007: First version of toyDat2.
*     ------------------------------------------------------------------
      integer
     &     i, j
*     ------------------------------------------------------------------
*     Give a name to the Problem.

      Prob   = 'Toy NLP2'

      ne     =  9
      n      =  4
      m      =  4

      nnCon  =  2
      nnJac  =  2
      nnObj  =  3

*     Check if there is enough storage.

      Errors = 0
      if (m      .gt. maxm ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (ne     .gt. maxne) Errors = 1
      if (Errors .gt.     0) return

*     -------------------------------------
*     Set up the list of row indices in indA.
*     -------------------------------------
*     Column  1
*     Nonlinear elements in rows (1, 2)  first.

      locA(1) =  1

      indA(1) =  1
      Acol(1) =  0.0d+0

*     Linear element in row 3 next.

      indA(2) =  3
      Acol(2) =  2.0d+0

*     Column 2.
*     Nonlinear elements in rows (1, 2).

      locA(2) =  3

      indA(3) =  1
      Acol(3) =  0.0d+0

      indA(4) =  2
      Acol(4) =  0.0d+0

*     Linear element in row 3.

      indA(5) =  3
      Acol(5) =  4.0d+0

*     Column 3.
*     Linear element in row 1.

      locA(3) =  6

      indA(6) =  1
      Acol(6) =  1.0d+0

*     Objective row element in row 4.

      indA(7) =  4
      Acol(7) =  3.0d+0

*     Column 4.
*     Linear element in row 2

      locA(4) =  8

      indA(8) =  2
      Acol(8) =  1.0d+0

*     Objective row element in row 4

      indA(9) =  4
      Acol(9) =  5.0d+0

*     Don't forget to finish off  locA.
*     This is crucial.

      locA(5) =  ne + 1

*     ------------------------------------------------------------------
*     Constraint ranges
*     ------------------------------------------------------------------
*     Nonlinear constraints first.

      bl(n+1) =  2.0d+0
      bu(n+1) =  2.0d+0

      bl(n+2) =  4.0d+0
      bu(n+2) =  4.0d+0

*     Followed by the linear constraints.

      bl(n+3) =  0.0d+0
      bu(n+3) =  infBnd

*     The linear objective term is row 4.

      iObj    = 4

*     The objective row is a free row.

      bl(n+4) = -infBnd
      bu(n+4) =  infBnd

      ObjAdd  = 0.0d+0

*     ------------------------------------------------------------------
*     Variable ranges
*     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -infBnd
         bu(j) =  infBnd
      end do

      bl(3) =  0.0d+0
      bl(4) =  0.0d+0

*     ------------------------------------------------------------------
*     Initialize x, hs and pi.
*     Set the initial value and status of each variable.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         x(j)  = 0.0d+0
      end do

      do j = 1, n
         hs(j) = 0
      end do

      do i = 1, m
         pi(i) = 0.0d+0
      end do

      end ! subroutine toyDat2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funcon2
     &   ( mode, nnCon, nnJac, neJac,
     &     x, fCon, gCon, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

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
*     Toy NLP problem from the SNOPT User's Guide.
*     ==================================================================
      integer
     &     nOut
*     ------------------------------------------------------------------
      nOut   = 6

      if (nState .eq. 1) then    ! First entry
         if (nOut .gt. 0) write(nOut, '(/a)') ' This is problem  Toy2'
      end if

      if (mode .eq. 0  .or.  mode .eq. 2) then
         fCon( 1)  = x(1)**2 +  x(2)**2
         fCon( 2)  =            x(2)**4
      end if

      if (mode .ge. 1) then
         gCon(1) = 2.0d+0*x(1)  ! Jacobian elements for column 1

         gCon(2) = 2.0d+0*x(2)  ! Jacobian elements for column 2
         gCon(3) = 4.0d+0*x(2)**3
      end if

      if (nState .ge. 2) then    ! Last entry
         if (nOut .gt. 0) write(nOut, '(/a)') ' Finished problem Toy2'
      end if

      end ! subroutine funcon2

