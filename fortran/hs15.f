!     ------------------------------------------------------------------
!     File hs15.f
!
!     This is a main program to illustrate the use of snOptB,
!     which is part of the SNOPT 7 package.
!
!     22 Apr 2020: First version.
!     ------------------------------------------------------------------
      program
     &     hs15

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, nName
      parameter
     &     ( maxm  = 2,
     &       maxn  = 2,
     &       maxne = 4,
     &       nName = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxne), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenru, leniu, lencu, lenrw, leniw, lencw
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  lenrw = 2000)
      double precision   rw(lenrw)
      parameter          (  leniw = 1000)
      integer            iw(leniw)
      parameter          (  lencw =  500)
      character          cw(lencw)*8
!     ------------------------------------------------------------------
!     No user workspace (cw, iw, and rw could be used instead).

      parameter          (  lenru = 1)
      double precision   ru(lenru)
      parameter          (  leniu = 1)
      integer            iu(leniu)
      parameter          (  lencu = 1)
      character          cu(lencu)*8
!     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      integer
     &     i1, i2, INFO, iObj, iPrint, iSpecs, iSumm, itnlim, lunit,
     &     m, MjrLim, maxiu, maxru, mincw, miniw, minrw,
     &     n, ne, nInf, nlocA, nnCon, nnJac, nnObj,
     &     nOut, nS
      double precision
     &     Obj, objAdd, sInf
      external
     &     hs15Con, hs15Obj
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).
!     nOut    is an output file used here by snmain.

      iSpecs = 4
      iPrint = 9
      iSumm  = 6
      nOut   = 6

!     Open the Specs and print files.
      lfile = 'hs15.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'hs15.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

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

      if (INFO .ne. 101) then
         go to 990
      end if

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      call hs15Data
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     ne, Acol, nlocA, indA, locA,
     &     bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     i1 and i2 may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      i1     = 0
      i2     = 0

      MjrLim = 5
      call snSeti
     &   ( 'Major Iterations', MjrLim, i1, i2, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     Sticky parameters yes is required for HOT start
      call snSet
     &   ( 'Sticky parameters yes', i1, i2, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Go for it, using a Cold start.
!     hs     need not be set if a basis file is to be input.
!            Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
!            The values are used by the Crash procedure m2crsh
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
     &     iObj, objAdd, Prob,
     &     hs15Con, hs15Obj,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Warm', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     hs15Con, hs15Obj,
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
      write(nOut, *) 'hs15  finished.'
      write(nOut, *) 'INFO  =', INFO
      write(nOut, *) 'nInf  =', nInf
      write(nOut, *) 'sInf  =', sInf
      write(nOut, *) 'Obj   =', Obj
      if (INFO .ge. 10) go to 910
      stop

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, 4010) 'Error while opening unit', lunit
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! program hs15

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs15Data
     &   ( Prob, maxm, maxn, maxne, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     ne, Acol, nlocA, indA, locA,
     &     bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxne, INFO, m, n, ne, nlocA,
     &     nnCon, nnObj, nnJac, iObj,
     &     indA(maxne), hs(maxn+maxm), locA(maxn+1)
      double precision
     &     objAdd, Acol(maxne) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character
     &     Prob*8

!     ------------------------------------------------------------------
!     hs15Data generates data for problem hsmain.
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
!     where d would be row iobj of A.
!     f(x) involves only the FIRST nnObj variables.
!
!     On entry,
!     maxm, maxn, maxne are upper limits on m, n, ne.
!
!     On exit,
!     INFO    is 0 if there is enough storage, 1 otherwise.
!     m       is the number of nonlinear and linear constraints.
!     n       is the number of variables.
!     ne     is the number of nonzeros in Acol(*).
!     nnCon   is the number of nonlinear constraints (they come first).
!     nnObj   is the number of nonlinear objective variables.
!     nnJac   is the number of nonlinear Jacobian variables.
!     Acol    is the constraint matrix (Jacobian), stored column-wise.
!     indA    is the list of row indices for each nonzero in Acol(*).
!     locA    is a set of pointers to the start of each column of Acol.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n) is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!
!     12 Aug 2005: First version of hs15Data.
!     ------------------------------------------------------------------
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision   bplus
      parameter         (bplus   = 1.0d+20)
!     ------------------------------------------------------------------
!     Give a name to the Problem.

      Prob  = 'hs15    '

      ne    =  4
      n     =  2
      m     =  2

      nnCon =  2
      nnJac =  2
      nnObj =  2

      nlocA = n   + 1

!     Check if there is enough storage.

      INFO  = 0
      if (m    .gt. maxm ) INFO = 1
      if (n    .gt. maxn ) INFO = 1
      if (ne   .gt. maxne) INFO = 1
      if (INFO .gt.   0  ) return

!     -------------------------------------
!     Set up the list of row indices in indA.
!     -------------------------------------
!     Column  1.

      ne = 0
      locA( 1) =  1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  0.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  1.0d+0

!     -------------------------------------------
!     Column 2.
!     Nonlinear elements in rows (1, 2).
!     -------------------------------------------
      locA( 2) =  ne + 1

      ne       =  ne + 1
      indA(ne) =  1
      Acol(ne) =  0.0d+0

      ne       =  ne + 1
      indA(ne) =  2
      Acol(ne) =  0.0d+0

!     -----------------------------------------
!     Don't forget to finish off locA.
!     This is crucial.
!     -----------------------------------------
      locA(n+1)= ne + 1

!     ------------------------------------------------------------------
!     Constraint ranges
!     ------------------------------------------------------------------
!     Nonlinear constraints first.

      bl(n+1) = 1.0d+0
      bu(n+1) = bplus

      bl(n+2) = 0.0d+0
      bu(n+2) = bplus

!     No linear objective term for this problem.

      iObj    = 0
      objAdd  = 0.0d+0

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -bplus
         bu(j) =  bplus
      end do

      bu(1) = 0.5d+0

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  -2.0d+0
      x(2)   =   1.0d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine hs15Dat

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs15Obj
     &   ( mode, nnObj, x, fObj, gObj, State,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, State, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     This is the objective for hs15.  The variables and constraints
!     have been reordered so the the nonlinear variables and constraints
!     appear first.
!     ==================================================================

      fObj    = 100.0d0 * (x(2) - x(1)**2)**2 + (1-x(1))**2

      gObj(1) = -400.0d0*(x(2) - x(1)**2) * x(1) - 2.0d0*(1.0d0 - x(1))
      gObj(2) = 200.0d0*(x(2) - x(1)**2)

      end ! subroutine hs15Obj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs15Con
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, State, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character*8
     &     cu(lencu)

!     ------------------------------------------------------------------
!     This is userfn for interface snoptb.
!
!     No user-defined storage is used.
!     ------------------------------------------------------------------
      integer
     &     neG
!     ------------------------------------------------------------------
      integer            Out
      parameter         (Out = 6)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (      State .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  hs15'
      else  if (State .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing hs15'
      end if

!     ------------------------------------------------------------------
!     Constraints and Objective.
!     ------------------------------------------------------------------
      fCon( 1) = x(1)*x(2)
      fCon( 2) = x(1) + x(2)**2

      neG = 0

!     -------------------------------------------
!     Nonlinear elements for column 1 (g=df/dx1).
!     g rows = (1).
!     -------------------------------------------
      neG = neG + 1
      gCon(neG) =  x(2)

      neG = neG + 1
      gCon(neG) =  1.0d0

!     -------------------------------------------
!     Nonlinear elements for column 2 (g=df/dx2).
!     g Rows = (1,2).
!     -------------------------------------------
      neG = neG + 1
      gCon(neG) =  x(1)

      neG = neG + 1
      gCon(neG) = 2.0d+0*x(2)

      end ! subroutine hs15Con
