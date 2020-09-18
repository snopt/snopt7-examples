!     ------------------------------------------------------------------
!     File snmainc.f
!     This is a main program to illustrate the use of snOptC,
!     which is part of the SNOPT 7 package.
!
!     31 Jul 1996: First   version.
!     19 Oct 2003: Updated for SNOPT 7
!     ------------------------------------------------------------------
      program
     &     snmain

      implicit
     &     none
      integer
     &     maxm, maxn, maxneA, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxneA = 3000,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indA(maxneA), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxneA), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenru, leniu, lencu, lenrw, leniw, lencw
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character          cw(lencw)*8
!     ------------------------------------------------------------------
!     User workspace (not used)

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
     &     Errors, iPrt, iSum, INFO, iObj, iPrint, iSpecs, iSumm,
     &     itnlim, j, m, mincw, miniw, minrw, n, neA, nInf, nnCon,
     &     nnJac, nnObj, nOut, nS
      double precision
     &     Obj, ObjAdd, sInf
      external
     &     hexfun, hexfun1, hexfun2

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

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'snmainc.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'snmainc.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
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
         go to 990
      end if

!     ------------------------------------------------------------------
!     1. Solve hexagon as a minimization problem with a cold start.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-----------------------------------------------'
      write(nOut, *) '1. Hexagon (minimize) problem with a cold start'
      write(nOut, *) '-----------------------------------------------'

      call hexData
     &   ( Prob, maxm, maxn, maxneA, INFO,
     &     m, n, neA, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     iPrt and iSum may refer to the Print and Summary files.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors = 0

      itnlim = 2500
      iPrt   =    0
      iSum   =    0
      call snSeti
     &   ( 'Iterations        ', itnlim, iPrt, iSum, Errors,
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
      call snOptC
     &   ( 'Cold', m, n, neA, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     hexfun,
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
      write(nOut, *) 'snmainc (1) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910

!     ------------------------------------------------------------------
!     2. Solve hexagon with a war start.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------------'
      write(nOut, *) '2. Hexagon problem with a warm start'
      write(nOut, *) '------------------------------------'

!     Set the option for a warm start.
!     Setting the warm start as an optional parameter takes precedence
!     over the cold/warm start argument in the call to snOptC.

      call snSet ( 'Warm start ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snSeti
     &   ( 'Minor print level', 0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Summary frequency', 0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSeti
     &   ( 'Print   frequency', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptC
     &   ( 'Cold', m, n, neA, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     hexfun,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *) ' '
      write(nOut, *) 'snmainc (2) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910

!     ------------------------------------------------------------------
!     3. Solve hexagon with some missing derivatives.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------'
      write(nOut, *) '3. Hexagon with three missing obj gradients'
      write(nOut, *) '-------------------------------------------'

      call snSet
     &   ( 'Defaults', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call hexData
     &   ( Prob, maxm, maxn, maxneA, INFO,
     &     m, n, neA, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      call snSeti
     &   ( 'Derivative level', 2, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptC
     &   ( 'Cold', m, n, neA, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     hexfun2,
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
      write(nOut, *) 'snmainc (3) (missing obj gradients) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910

!     ------------------------------------------------------------------
!     4. Solve hexagon with some missing derivatives.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------'
      write(nOut, *) '4. Hexagon with no constraint gradients'
      write(nOut, *) '---------------------------------------'

      call hexData
     &   ( Prob, maxm, maxn, maxneA, INFO,
     &     m, n, neA, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      call snSeti
     &   ( 'Derivative level', 1, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptC
     &   ( 'Cold', m, n, neA, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     hexfun1,
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
      write(nOut, *) 'snmainc (4) (no constraint gradients) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptC INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *) 'Obj           =', ObjAdd + Obj
      end if
      if (INFO .ge. 30) go to 910

!     ------------------------------------------------------------------
!     5. Solve hexagon with some missing derivatives.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-------------------------------------------------'
      write(nOut, *) '5. Find a feasible point for the Hexagon problem.'
      write(nOut, *) '   No constraint derivatives specified'
      write(nOut, *) '------------------------------------------------'

      call hexData
     &   ( Prob, maxm, maxn, maxneA, INFO,
     &     m, n, neA, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      do j = 1, n
         x(j) = 10.0d+0
      end do

      call snSeti
     &   ( 'Derivative level', 0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Feasible point', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptC
     &   ( 'Cold', m, n, neA, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, Prob,
     &     hexfun1,
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
      write(nOut, *) 'snmainc (5) (no constraint gradients) finished.'
      write(nOut, *) 'Input errors  =', Errors
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

!     ------------------------------------------------------------------
!     Error exit.
!     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program snoptc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hexData
     &   ( Prob, maxm, maxn, maxneA, INFO,
     &     m, n, neA, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Acol, indA, locA, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxneA, INFO, m, n, neA, nnCon, nnObj, nnJac,
     &     iObj, indA(maxneA) , hs(maxn+maxm), locA(maxn+1)
      double precision
     &     ObjAdd, Acol(maxneA) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character
     &     Prob*8

!     ------------------------------------------------------------------
!     hexData generates data for the Hexagon problem.
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
!     maxm, maxn, maxneA are upper limits on m, n, neA.
!
!     On exit,
!     INFO    is 0 if there is enough storage, 1 otherwise.
!     m       is the number of nonlinear and linear constraints.
!     n       is the number of variables.
!     neA     is the number of nonzeros in Acol(*).
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
!     24 Dec 1997: First version of hexData.
!     ------------------------------------------------------------------
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision   bplus
      parameter         (bplus   = 1.0d+20)
      double precision   zero,               one
      parameter         (zero    = 0.0d+0,   one    = 1.0d+0)
!     ------------------------------------------------------------------
!     Give a name to the Problem.

      Prob   = 'Hexagon '

      neA    = 52
      n      =  9
      m      = 18

      nnCon  = 14
      nnJac  =  n
      nnObj  =  n

!     Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm  ) INFO = 1
      if (n     .gt. maxn  ) INFO = 1
      if (neA   .gt. maxneA) INFO = 1
      if (INFO  .gt. 0     ) return

!     ---------------------------------------
!     Set up the list of row indices in indA.
!     ---------------------------------------
      neA       =  0

!     ---------------------------------------------------
!     Column  1
!     First,  nonlinear elements in rows (1, 2, 3, 4, 5).
!     ---------------------------------------------------
      locA( 1)  =  1

      neA       =  neA + 1
      indA(neA) =  1
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  2
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  3
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  4
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  5
      Acol(neA) =  zero

!     Column 1.
!     Next,  linear element in row 15.

      neA       =  neA + 1
      indA(neA) =  15
      Acol(neA) = -one

!     -------------------------------------------
!     Column 2.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).
!     -------------------------------------------
      locA( 2)  =  neA + 1

      neA       =  neA + 1
      indA(neA) =  2
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  6
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  7
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  8
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  9
      Acol(neA) =  zero

!     Column 2.
!     Linear elements in rows (15, 16).

      neA       =  neA + 1
      indA(neA) =  15
      Acol(neA) =  one

      neA       =  neA + 1
      indA(neA) =  16
      Acol(neA) = -one

!     ------------------------------------------
!     Column 3.
!     Nonlinear elements in rows (3, 7, 10, 11).
!     ------------------------------------------
      locA( 3)  =  neA + 1

      neA       =  neA + 1
      indA(neA) =  3
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  7
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  10
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  11
      Acol(neA) =  zero

!     Column 3.
!     Linear elements in rows (16, 17).

      neA       =  neA + 1
      indA(neA) =  16
      Acol(neA) =  one

      neA       =  neA + 1
      indA(neA) =  17
      Acol(neA) =  one

!     -------------------------------------------
!     Column 4.
!     Nonlinear elements in rows (20, 21, 22, 23, 24).
!     -------------------------------------------
      locA( 4)  =  neA + 1

      neA       =  neA + 1
      indA(neA) =  4
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  8
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  10
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  12
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  13
      Acol(neA) =  zero

!     Column 4.
!     Linear elements in rows (17, 18).

      neA       =  neA + 1
      indA(neA) =  17
      Acol(neA) = -one

      neA       =  neA + 1
      indA(neA) =  18
      Acol(neA) =  one

!     -------------------------------------------
!     Column 5.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).
!     -------------------------------------------
      locA( 5)  =  neA + 1

      neA       =  neA + 1
      indA(neA) =  5
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  9
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  11
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  13
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  14
      Acol(neA) =  zero

!     Column 5.
!     Linear element in row 18.

      neA       =  neA + 1
      indA(neA) =  18
      Acol(neA) = -one

!     ------------------------------------------
!     Column 6.
!     Nonlinear elements in rows (1, 2, 3, 4, 5, 6).
!     ------------------------------------------
      locA(6)   =  neA + 1

      neA       =  neA + 1
      indA(neA) =  1
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  2
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  3
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  4
      Acol(neA) =  zero

      neA       = neA + 1
      indA(neA) =  5
      Acol(neA) =  zero

!     -------------------------------------------
!     Column 7.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).
!     -------------------------------------------
      locA(7)   =  neA + 1

      neA       =  neA + 1
      indA(neA) =  2
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  6
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  7
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  8
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  9
      Acol(neA) =  zero

!     -------------------------------------------
!     Column 8.
!     Nonlinear elements in rows (4, 8, 10, 12, 13).
!     -------------------------------------------
      locA(8)   =   neA + 1

      neA       =  neA + 1
      indA(neA) =  4
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  8
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  10
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  12
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  13
      Acol(neA) =  zero

!     -------------------------------------------
!     Column 9.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).
!     -------------------------------------------
      locA(9)   =  neA + 1

      neA       =  neA + 1
      indA(neA) =  5
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  9
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  11
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  13
      Acol(neA) =  zero

      neA       =  neA + 1
      indA(neA) =  14
      Acol(neA) =  zero

!     -----------------------------------------
!     Don't forget to finish off locA and locH.
!     This is crucial.
!     -----------------------------------------
      locA(10)  =  neA + 1

!     ------------------------------------------------------------------
!     Constraint ranges
!     ------------------------------------------------------------------
!     Nonlinear constraints first.

      do i = 1, nnCon
         bl(i+n) = -bplus
         bu(i+n) =  one
      end do

!     Followed by the linear constraints.

      do i = nnCon+1, m
         bl(i+n) =  zero
         bu(i+n) =  bplus
      end do

!     No linear objective term for this problem.

      iObj    = 0
      ObjAdd  = zero

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, n
         bl(j) = -bplus
         bu(j) =  bplus
      end do

      bl(1) =  zero
      bl(3) = -one
      bl(5) =  zero
      bl(6) =  zero
      bl(7) =  zero

      bu(3) =  one
      bu(8) =  zero
      bu(9) =  zero

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  .1d+0
      x(2)   =  .125d+0
      x(3)   =  .666666d+0
      x(4)   =  .142857d+0
      x(5)   =  .111111d+0
      x(6)   =  .2d+0
      x(7)   =  .25d+0
      x(8)   = -.2d+0
      x(9)   = -.25d+0

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = zero
      end do

      end ! subroutine hexData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hexfun
     &   ( mode, nnObj, nnCon, nnJac, nnL, negCon,
     &     x, fObj, gObj, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnObj, nnCon, nnJac, nnL, negCon, State,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     fObj, fCon(nnCon), gCon(negCon), gObj(nnObj),
     &     x(nnL), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem Hexagon.
!
!     mode
!     ----
!      0  fObj, fCon
!      1  gObj, gCon
!      2  fObj, fCon, gObj, gCon
!
!     No user-defined storage is used.
!     ==================================================================
      logical
     &     needf, needg
      integer
     &     neG
!     ------------------------------------------------------------------
      integer            Out
      parameter         (Out = 6)
      double precision   two
      parameter         (two = 2.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (      State .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  Hex (Minimize)'
      else  if (State .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing Hex'
      end if

!     ------------------------------------------------------------------
!     Constraints and Objective   f = [ fCon fObj ].
!     ------------------------------------------------------------------
      needf = mode .eq. 0  .or.  mode .eq. 2
      needg = mode .eq. 1  .or.  mode .eq. 2

      if (needf) then
         fCon( 1) =    x(1)**2          +   x(6)**2
         fCon( 2) =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         fCon( 3) =   (x(3) - x(1))**2  +   x(6)**2
         fCon( 4) =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         fCon( 5) =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         fCon( 6) =    x(2)**2          +   x(7)**2
         fCon( 7) =   (x(3) - x(2))**2  +   x(7)**2
         fCon( 8) =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         fCon( 9) =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         fCon(10) =   (x(4) - x(3))**2  +   x(8)**2
         fCon(11) =   (x(5) - x(3))**2  +   x(9)**2
         fCon(12) =    x(4)**2          +   x(8)**2
         fCon(13) =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         fCon(14) =    x(5)**2          +   x(9)**2
         fObj     =  - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &                                       + x(4)*x(9) + x(3)*x(8)
      end if

      neG = 0

      if (needg) then
!        -------------------------------------------
!        Nonlinear elements for column 1 (g=df/dx1).
!        g rows = (1,2,3,4,5,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two* x(1)         ! row  1
         neG       =  neG + 1
         gCon(neG) = -two*(x(2) - x(1)) ! row  2
         neG       =  neG + 1
         gCon(neG) = -two*(x(3) - x(1)) ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(1) - x(4)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(1) - x(5)) ! row  5

         gObj(  1) =       x(7)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 2 (g=df/dx2).
!        g Rows = (2,6,7,8,9,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(2) - x(1)) ! row 2
         neG       =  neG + 1
         gCon(neG) =  two* x(2)         ! row 6
         neG       =  neG + 1
         gCon(neG) = -two*(x(3) - x(2)) ! row 7
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(2)) ! row 8
         neG       =  neG + 1
         gCon(neG) =  two*(x(2) - x(5)) ! row 9

         gObj( 2)  =      -x(6)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 3 (g=df/dx3).
!        g Rows = (3,7,10,11,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(3) - x(1)) ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(3) - x(2)) ! row  7
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(3)) ! row 10
         neG       =  neG + 1
         gCon(neG) = -two*(x(5) - x(3)) ! row 11

         gObj(  3) =     - x(7) + x(8)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 4 (g=df/dx4).
!        g Rows = (4,8,10,12,13,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(1) - x(4)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(2)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(3)) ! row 10
         neG       =  neG + 1
         gCon(neG) =  two* x(4)         ! row 12
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(5)) ! row 13

         gObj(  4) =       x(9)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 5 (g=df/dx5).
!        g Rows = (5,9,11,13,14,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(1) - x(5)) ! row  5
         neG       =  neG + 1
         gCon(neG) = -two*(x(2) - x(5)) ! row  9
         neG       =  neG + 1
         gCon(neG) =  two*(x(5) - x(3)) ! row 11
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(5)) ! row 13
         neG       =  neG + 1
         gCon(neG) =  two* x(5)         ! row 14

         gObj( 5)  = -      x(8)        ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 6 (g=df/dx6).
!        g Rows = (1,2,3,4,5,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two* x(6)         ! row  1
         neG       =  neG + 1
         gCon(neG) = -two*(x(7) - x(6)) ! row  2
         neG       =  neG + 1
         gCon(neG) =  two* x(6)         ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(6) - x(8)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(6) - x(9)) ! row  5

         gObj(  6) =      -x(2)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 7 (g=df/dx7).
!        g rows = (2,6,7,8,9,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(7) - x(6)) ! row  2
         neG       =  neG + 1
         gCon(neG) =  two* x(7)         ! row  6
         neG       =  neG + 1
         gCon(neG) =  two* x(7)         ! row  7
         neG       =  neG + 1
         gCon(neG) = -two*(x(8) - x(7)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*(x(7) - x(9)) ! row  9

         gObj(  7) =     - x(3) + x(1)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 8 (g=df/dx8).
!        g Rows = (4,8,10,12,13,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(6) - x(8)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(8) - x(7)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*x(8)          ! row 10
         neG       =  neG + 1
         gCon(neG) =  two*x(8)          ! row 12
         neG       =  neG + 1
         gCon(neG) = -two*(x(9) - x(8)) ! row 13

         gObj(  8) =     - x(5) + x(3)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 9 (g=df/dx9).
!        g Rows = (5,9,11,13,14,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(6) - x(9)) ! row  5
         neG       =  neG + 1
         gCon(neG) = -two*(x(7) - x(9)) ! row  9
         neG       =  neG + 1
         gCon(neG) =  two* x(9)         ! row 11
         neG       =  neG + 1
         gCon(neG) =  two*(x(9) - x(8)) ! row 13
         neG       =  neG + 1
         gCon(neG) =  two* x(9)         ! row 14
         gObj(  9) =       x(4)         ! Obj 15

      end if

      end ! subroutine hexfun

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hexfun1
     &   ( mode, nnObj, nnCon, nnJac, nnL, negCon,
     &     x, fObj, gObj, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnObj, nnCon, nnJac, nnL, negCon, State,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     fObj, fCon(nnCon), gCon(negCon), gObj(nnObj),
     &     x(nnL), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem Hexagon1
!
!     mode
!     ----
!      0  fObj, fCon
!      1  gObj, gCon
!      2  fObj, fCon, gObj, gCon
!
!     No constraint gradients.
!     ==================================================================
      logical
     &     needf, needg
!     ------------------------------------------------------------------
      integer            Out
      parameter         (Out = 6)
      double precision   two
      parameter         (two = 2.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (      State .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  Hex (Minimize)'
      else  if (State .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing Hex'
      end if

!     ------------------------------------------------------------------
!     Constraints and Objective   f = [ fCon fObj ].
!     ------------------------------------------------------------------
      needf = mode .eq. 0  .or.  mode .eq. 2
      needg = mode .eq. 1  .or.  mode .eq. 2

      if (needf) then
         fCon( 1) =    x(1)**2          +   x(6)**2
         fCon( 2) =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         fCon( 3) =   (x(3) - x(1))**2  +   x(6)**2
         fCon( 4) =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         fCon( 5) =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         fCon( 6) =    x(2)**2          +   x(7)**2
         fCon( 7) =   (x(3) - x(2))**2  +   x(7)**2
         fCon( 8) =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         fCon( 9) =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         fCon(10) =   (x(4) - x(3))**2  +   x(8)**2
         fCon(11) =   (x(5) - x(3))**2  +   x(9)**2
         fCon(12) =    x(4)**2          +   x(8)**2
         fCon(13) =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         fCon(14) =    x(5)**2          +   x(9)**2
         fObj     =  - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &                                       + x(4)*x(9) + x(3)*x(8)
      end if

      if (needg) then
         gObj( 1) =       x(7)
         gObj( 2) =      -x(6)
         gObj( 3) =     - x(7) + x(8)
         gObj( 4) =       x(9)
         gObj( 5) = -     x(8)
         gObj( 6) =      -x(2)
         gObj( 7) =     - x(3) + x(1)
         gObj( 8) =     - x(5) + x(3)
         gObj( 9) =       x(4)
      end if

      end ! subroutine hexfun1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hexfun2
     &   ( mode, nnObj, nnCon, nnJac, nnL, negCon,
     &     x, fObj, gObj, fCon, gCon,
     &     State, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnObj, nnCon, nnJac, nnL, negCon, State,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     fObj, fCon(nnCon), gCon(negCon), gObj(nnObj),
     &     x(nnL), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Problem Hexagon2
!
!     mode
!     ----
!      0  fObj, fCon
!      1  gObj, gCon
!      2  fObj, fCon, gObj, gCon
!
!     Three objective gradients are missing.
!     ==================================================================
      logical
     &     needf, needg
      integer
     &     neG
!     ------------------------------------------------------------------
      integer            Out
      parameter         (Out = 6)
      double precision   two
      parameter         (two = 2.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (      State .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  Hex (Minimize)'
      else  if (State .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing Hex'
      end if

!     ------------------------------------------------------------------
!     Constraints and Objective   f = [ fCon fObj ].
!     ------------------------------------------------------------------
      needf = mode .eq. 0  .or.  mode .eq. 2
      needg = mode .eq. 1  .or.  mode .eq. 2

      if (needf) then
         fCon( 1) =    x(1)**2          +   x(6)**2
         fCon( 2) =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         fCon( 3) =   (x(3) - x(1))**2  +   x(6)**2
         fCon( 4) =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         fCon( 5) =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         fCon( 6) =    x(2)**2          +   x(7)**2
         fCon( 7) =   (x(3) - x(2))**2  +   x(7)**2
         fCon( 8) =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         fCon( 9) =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         fCon(10) =   (x(4) - x(3))**2  +   x(8)**2
         fCon(11) =   (x(5) - x(3))**2  +   x(9)**2
         fCon(12) =    x(4)**2          +   x(8)**2
         fCon(13) =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         fCon(14) =    x(5)**2          +   x(9)**2
         fObj     =  - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &                                       + x(4)*x(9) + x(3)*x(8)
      end if

      neG = 0

      if (needg) then
!        -------------------------------------------
!        Nonlinear elements for column 1 (g=df/dx1).
!        g rows = (1,2,3,4,5,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two* x(1)         ! row  1
         neG       =  neG + 1
         gCon(neG) = -two*(x(2) - x(1)) ! row  2
         neG       =  neG + 1
         gCon(neG) = -two*(x(3) - x(1)) ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(1) - x(4)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(1) - x(5)) ! row  5

!        gObj(  1) =       x(7)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 2 (g=df/dx2).
!        g Rows = (2,6,7,8,9,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(2) - x(1)) ! row 2
         neG       =  neG + 1
         gCon(neG) =  two* x(2)         ! row 6
         neG       =  neG + 1
         gCon(neG) = -two*(x(3) - x(2)) ! row 7
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(2)) ! row 8
         neG       =  neG + 1
         gCon(neG) =  two*(x(2) - x(5)) ! row 9

!        gObj( 2)  =      -x(6)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 3 (g=df/dx3).
!        g Rows = (3,7,10,11,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(3) - x(1)) ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(3) - x(2)) ! row  7
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(3)) ! row 10
         neG       =  neG + 1
         gCon(neG) = -two*(x(5) - x(3)) ! row 11

!        gObj(  3) =     - x(7) + x(8)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 4 (g=df/dx4).
!        g Rows = (4,8,10,12,13,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(1) - x(4)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(2)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(3)) ! row 10
         neG       =  neG + 1
         gCon(neG) =  two* x(4)         ! row 12
         neG       =  neG + 1
         gCon(neG) =  two*(x(4) - x(5)) ! row 13

         gObj(  4) =       x(9)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 5 (g=df/dx5).
!        g Rows = (5,9,11,13,14,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(1) - x(5)) ! row  5
         neG       =  neG + 1
         gCon(neG) = -two*(x(2) - x(5)) ! row  9
         neG       =  neG + 1
         gCon(neG) =  two*(x(5) - x(3)) ! row 11
         neG       =  neG + 1
         gCon(neG) = -two*(x(4) - x(5)) ! row 13
         neG       =  neG + 1
         gCon(neG) =  two* x(5)         ! row 14

         gObj( 5)  = -      x(8)        ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 6 (g=df/dx6).
!        g Rows = (1,2,3,4,5,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two* x(6)         ! row  1
         neG       =  neG + 1
         gCon(neG) = -two*(x(7) - x(6)) ! row  2
         neG       =  neG + 1
         gCon(neG) =  two* x(6)         ! row  3
         neG       =  neG + 1
         gCon(neG) =  two*(x(6) - x(8)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(6) - x(9)) ! row  5

         gObj(  6) =      -x(2)         ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 7 (g=df/dx7).
!        g rows = (2,6,7,8,9,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) =  two*(x(7) - x(6)) ! row  2
         neG       =  neG + 1
         gCon(neG) =  two* x(7)         ! row  6
         neG       =  neG + 1
         gCon(neG) =  two* x(7)         ! row  7
         neG       =  neG + 1
         gCon(neG) = -two*(x(8) - x(7)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*(x(7) - x(9)) ! row  9

         gObj(  7) =     - x(3) + x(1)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 8 (g=df/dx8).
!        g Rows = (4,8,10,12,13,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(6) - x(8)) ! row  4
         neG       =  neG + 1
         gCon(neG) =  two*(x(8) - x(7)) ! row  8
         neG       =  neG + 1
         gCon(neG) =  two*x(8)          ! row 10
         neG       =  neG + 1
         gCon(neG) =  two*x(8)          ! row 12
         neG       =  neG + 1
         gCon(neG) = -two*(x(9) - x(8)) ! row 13

         gObj(  8) =     - x(5) + x(3)  ! Obj 15

!        -------------------------------------------
!        Nonlinear elements for column 9 (g=df/dx9).
!        g Rows = (5,9,11,13,14,15).
!        -------------------------------------------
         neG       =  neG + 1
         gCon(neG) = -two*(x(6) - x(9)) ! row  5
         neG       =  neG + 1
         gCon(neG) = -two*(x(7) - x(9)) ! row  9
         neG       =  neG + 1
         gCon(neG) =  two* x(9)         ! row 11
         neG       =  neG + 1
         gCon(neG) =  two*(x(9) - x(8)) ! row 13
         neG       =  neG + 1
         gCon(neG) =  two* x(9)         ! row 14
         gObj(  9) =       x(4)         ! Obj 15

      end if

      end ! subroutine hexfun2
