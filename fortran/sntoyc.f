!     ------------------------------------------------------------------
!     File sntoyc.f
!     This is a main program to illustrate the use of subroutine snOptC,
!     which is part of the SNOPT package.
!
!     08 Jun 2005: First version.
!     14 Nov 2010: Fixed bug in definition of indA.
!     ------------------------------------------------------------------
      program
     &     sntoyc

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxne  = 3000,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne) , hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenrw, leniw, lencw
      !-----------------------------------------------------------------
      ! SNOPT workspace
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character*8        cw(lencw)
      !-----------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      external
     &     usrfun
      integer
     &     Errors, i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim, m,
     &     mincw, miniw, minrw, n, ne, nInf, nnCon, nnJac, nnObj, nOut,
     &     nS
      double precision
     &     Obj, ObjAdd, sInf
      !-----------------------------------------------------------------
      ! Specify some of the SNOPT files.
      ! iSpecs  is the Specs file   (0 if none).
      ! iPrint  is the Print file   (0 if none).
      ! iSumm   is the Summary file (0 if none).
      ! nOut    is an output file used here by snmain.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

      byname = .true.

      if ( byname ) then

         ! Unix and DOS systems.  Open the Specs and print files.

         lfile = 'sntoyc.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'sntoyc.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

      !-----------------------------------------------------------------
      ! First,  snInit MUST be called to initialize optional parameters
      ! to their default values.
      !-----------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Read a Specs file (Optional).
      !-----------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

      ! Set up the data structure for the sparse Jacobian.
      ! Assign dummy values for the nonlinear elements.

      call ToyDat
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      !-----------------------------------------------------------------
      ! Specify any options not set in the Specs file.
      ! i1 and i2 may refer to the Print and Summary file respectively.
      ! Setting them to 0 suppresses printing.
      !-----------------------------------------------------------------
      Errors = 0

      itnlim = 250
      i1     =   0
      i2     =   0
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Go for it, using a Cold start.
      ! hs   need not be set if a basis file is to be input.
      !      Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
      !      The values are used by the Crash procedure s2crsh
      !      to choose an initial basis B.
      !      If hs(j) = 0 or 1, column j is eligible for B.
      !      If hs(j) = 2, column j is initially superbasic (not in B).
      !      If hs(j) = 3, column j is eligible for B and is given
      !                    preference over columns with hs(j) = 0 or 1.
      !      If hs(j) = 4 or 5, column j is initially nonbasic.
      !-----------------------------------------------------------------
      call snOptC
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     usrfun,
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
      write(nOut, *) 'sntoyc finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910
      stop

      !-----------------------------------------------------------------
      ! Error exit.
      !-----------------------------------------------------------------
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
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nnCon, nnJac, nnL, neJac, nState,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnL), sum
      double precision
     &     fObj, gObj(nnObj)
      double precision
     &     fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

      !=================================================================
      ! Problem Toy.
      !    minimize (x1 + x2 +x3)^2 + 3*x3 + 5*x4
      !
      !     subject to    x1^2 +   x2^2 + x3       = 2
      !                            x2^4      + x4  = 4
      !                 2*x1   + 4*x2             >= 0.
      !
      !=================================================================
      !-----------------------------------------------------------------
      ! Objective.
      !-----------------------------------------------------------------
      sum = x(1) + x(2) + x(3)

      fObj = sum * sum

      gObj(1) = 2.0d+0 * sum
      gObj(2) = 2.0d+0 * sum
      gObj(3) = 2.0d+0 * sum

      !-----------------------------------------------------------------
      ! Constraints.
      !-----------------------------------------------------------------
      fCon(1) = x(1)**2 + x(2)**2
      fCon(2) = x(2)**4

      ! Nonlinear Jacobian elements for column 1.
      ! row  (1).

      gCon(1) = 2.0d+0*x(1)

      ! Nonlinear Jacobian elements for column 2.
      ! rows (1, 2).

      gCon(2) = 2.0d+0*x(2)
      gCon(3) = 4.0d+0*x(2)**3

      end ! subroutine usrfun

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ToyDat
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, ne, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, INFO, m, n, ne, nnCon, nnObj, nnJac,
     &     iObj, indA(maxne) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character*8
     &     Prob
      !-----------------------------------------------------------------
      ! ToyDat generates data for the Toy problem from the users manual.
      ! The constraints take the form
      !          c(x) + A*x - s = 0,
      ! where the Jacobian matrix of derivatives for c(x) + Ax is stored
      ! in Acol(*), and any terms coming from the derivatives of c(x)
      ! are stored in the TOP LEFT-HAND CORNER of Acol(*), with
      ! dimensions  nnCon x nnJac.

      ! Note that the right-hand side is zero.  s denotes the vector of
      ! slack variables whose bounds contain any constants that might
      ! have formed a right-hand side.
      !
      ! Problem Toy:
      !    minimize (x1 + x2 +x3)^2 + 3*x3 + 5*x4
      !
      !     subject to    x1^2 +   x2^2 + x3       = 2
      !                            x2^4      + x4  = 4
      !                 2*x1   + 4*x2             >= 0.
      !
      ! Implicitly, the objective function is written as
      !        f(x) + d'x
      ! where d is row iobj of A. In this case,  f(x) involves only
      ! the first nnObj variables.
      !
      !
      ! On entry,
      !  maxm, maxn, maxne are upper limits on m, n, ne.
      !
      ! On exit,
      !  INFO    is 0 if there is enough storage, 1 otherwise.
      !  m       is the number of nonlinear and linear constraints.
      !  n       is the number of variables.
      !  ne      is the number of nonzeros in Acol(*).
      !  nnCon   is the number of nonlinear constraints (listed first).
      !  nnObj   is the number of nonlinear objective variables.
      !  nnJac   is the number of nonlinear Jacobian variables.
      !  Acol    is the constraint Jacobian, stored column-wise.
      !  indA    is the list of row indices for each nonzero in Acol(*).
      !  locA    consists of the pointers to the start of each column
      !          of Acol.
      !  bl      is the lower bounds on x and s.
      !  bu      is the upper bounds on x and s.
      !  hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
      !  x (1:n) is a set of initial values for x.
      !  pi(1:m) is a set of initial values for the dual variables pi.
      !
      ! 08 Jun 2005: First version of ToyDat.
      ! 14 Nov 2010: Spiffed up for snopt8.
      !-----------------------------------------------------------------
      integer
     &     i, j
      !-----------------------------------------------------------------
      double precision   bplus
      parameter         (bplus   = 1.0d+20)
      !-----------------------------------------------------------------
      ! Define a name for the Problem.

      Prob   = 'Toy 1     '

      ne     = 9
      n      = 4
      m      = 4

      nnCon  = 2
      nnJac  = 2
      nnObj  = 3

      ! Linear objective term is in row 4.

      iObj    = 4
      ObjAdd  = 0.0d+0

      ! Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm ) INFO = 1
      if (n     .gt. maxn ) INFO = 1
      if (ne    .gt. maxne) INFO = 1
      if (INFO  .gt.   0  ) return

      !--------------------------------------------
      ! Set up the data structure for the Jacobian.
      ! The nonzeros are stored columnwise.
      ! Nonlinear elements must be set (0 is okay).
      !--------------------------------------------

      !===========================================
      ! Column  1
      !===========================================
      locA(1) = 1

      ! Nonlinear element in row 1.

      indA(1) = 1
      Acol(1) = 0.0d+0

      ! Linear element in row 3.

      indA(2) = 3
      Acol(2) = 2.0d+0

      !===========================================
      ! Column  2
      !===========================================
      locA(2) = 3

      ! Nonlinear elements in rows (1, 2).

      indA(3) = 1
      Acol(3) = 0.0d+0

      indA(4) = 2
      Acol(4) = 0.0d+0

      ! Linear element in row 3.

      indA(5) = 3
      Acol(5) = 4.0d+0

      !===========================================
      ! Column  3
      !===========================================
      locA(3) =  6

      ! No Nonlinear elements (as nnJac = 2)
      ! Linear elements in rows (1, 4).

      indA(6) = 1
      Acol(6) = 1.0d+0

      indA(7) = iObj            ! Objective row
      Acol(7) = 3.0d+0

      !===========================================
      ! Column  4
      !===========================================
      locA(4) = 8

      ! No Nonlinear elements.
      ! Linear elements in rows (2,4).

      indA(8) = 2
      Acol(8) = 1.0d+0

      indA(9) = iObj            ! Objective row
      Acol(9) = 5.0d+0

      !===========================================
      ! Don't forget to finish off  locA.
      ! This is crucial.
      !===========================================
      locA(5) =  ne + 1 ! = 10

      !-----------------------------------------------------------------
      ! Constraint ranges
      !-----------------------------------------------------------------
      ! Nonlinear constraints first.

      bl(n+1) = 2.0d+0
      bu(n+1) = 2.0d+0
      bl(n+2) = 4.0d+0
      bu(n+2) = 4.0d+0

      ! Followed by the linear constraints.

      bl(n+3) =  0.0d+0
      bu(n+3) =  bplus
      bl(n+4) = -bplus
      bu(n+4) =  bplus

      !-----------------------------------------------------------------
      ! Variable ranges
      !-----------------------------------------------------------------
      do j = 1, n
         bl(j) = -bplus
         bu(j) =  bplus
      end do

      bl(3) = 0.0d+0
      bl(4) = 0.0d+0

      !-----------------------------------------------------------------
      ! Initialize x, hs and pi.
      ! Set the initial value and status of each variable.
      ! For want of something better to do, make the variables x(1:n)
      ! temporarily fixed at their current values.
      ! The crash can set the rest.
      !-----------------------------------------------------------------
      x(1)   =  .1d+0
      x(2)   =  .125d+0
      x(3)   =  .666666d+0
      x(4)   =  .142857d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine toyDat
