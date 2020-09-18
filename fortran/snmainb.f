!     ------------------------------------------------------------------
!     File snmainb.f  (Unix version)
!     This is a main program to illustrate the use of subroutine snOptB,
!     which is part of the SNOPT 7 package.
!
!     31 Jul 1996: First   version.
!     19 Oct 2003: Updated for SNOPT 7
!     19 Oct 2003: Current version.
!     ------------------------------------------------------------------
      program
     &     snmainb

      implicit
     &     none
      integer
     &     maxm, maxn, maxneJ, nName
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxneJ = 3000,
     &       nName  = 1 )

      character
     &     Prob*8, Names(nName)*8
      integer
     &     indJ(maxneJ) , hs(maxn+maxm)
      integer
     &     locJ(maxn+1)
      double precision
     &     Jcol(maxneJ), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenrw, leniw, lencw
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   500)
      character*8        cw(lencw)
!     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      external
     &     HexCon,   HexObjMin,   HexObjMax, HexObjMax2, HexCon1,
     &     HexConFD, HexObjMinFD, HexObjMaxFD
      integer
     &     Errors, iPrt, iSum, INFO, iObj, iPrint, iSpecs, iSumm,
     &     itns, j, lvlDer, m, mincw, miniw, minrw,
     &     n, neJ, nInf, nnCon, nnJac, nnObj, nOut, nS
      double precision
     &     obj, objAdd, sInf
!     ------------------------------------------------------------------
!     Specify some of the SNOPT files.
!     iSpecs  is the Specs file   (0 if none).
!     iPrint  is the Print file   (0 if none).
!     iSumm   is the Summary file (0 if none).

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

!        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'snmainb.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'snmainb.out'
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

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     iPrt and iSum may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors =   0
      iPrt   =   0
      iSum   =   0

      call snSet
     &   ( 'Sticky parameters =  yes  ', iPrt, iSum, Errors,
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
      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon, HexObjMin,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (1) (minimize) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     2. Solve hexagon as a minimization problem. No derivatives.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------------'
      write(nOut, *) '2. Hexagon (minimize) problem. No derivatives'
      write(nOut, *) '---------------------------------------------'

!     Set up the data structure for the sparse Jacobian.
!     Assign dummy values for the nonlinear elements.

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     iPrt and iSum may refer to the Print and Summary file respectively.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      Errors =   0
      iPrt   =   0
      iSum   =   0

      lvlDer = 0
      call snSeti
     &   ( 'Derivative level          ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexConFD, HexObjMinFD,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (2) (minimize) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     3. Solve hexagon as a maximization problem.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-----------------------------------------------'
      write(nOut, *) '3. Hexagon (maximize) problem with a cold start'
      write(nOut, *) '-----------------------------------------------'

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      call snSet ( 'Maximize   ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      lvlDer = 3
      call snSeti
     &   ( 'Derivative level          ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon, HexObjMax,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (3) (maximize) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     4. Solve hexagon as a maximization problem. No derivatives.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------------'
      write(nOut, *) '4. Hexagon (maximize) problem. No derivatives'
      write(nOut, *) '---------------------------------------------'

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      lvlDer = 0
      call snSeti
     &   ( 'Derivative level          ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexConFD, HexObjMaxFD,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (4) (maximize, fd) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     5. Solve hexagon with some constraints gradients missing.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------------'
      write(nOut, *) '5. Hexagon with missing constraint derivatives'
      write(nOut, *) '---------------------------------------------'

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      call snSet
     &   ( 'Defaults', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snSet ( 'Maximize   ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      lvlDer = 1
      call snSeti
     &   ( 'Derivative level   ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon1, HexObjMax,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (5) (missing constr. grads) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     6. Solve hexagon with some obj gradients missing.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------------'
      write(nOut, *) '6. Hexagon with missing objective derivatives'
      write(nOut, *) '---------------------------------------------'

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      lvlDer = 2
      call snSeti
     &   ( 'Derivative level   ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon, HexObjMax2,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (6) (missing obj. grads) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     7. Find a feasible point for Hexagon. (maximize)
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-----------------------------------------------'
      write(nOut, *) '7. Hexagon (feasible point, maximized objective'
      write(nOut, *) '-----------------------------------------------'

      call HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, objAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      do j = 1, n
         x(j) = 10.0d+0
      end do

      lvlDer = 3
      call snSeti
     &   ( 'Derivative level   ', lvlDer, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Maximize           ',         iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Feasible point     ',         iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon, HexObjMax,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      call snGeti
     &   ( 'iw 421        ', itns, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (7) (feasible point) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

!     ------------------------------------------------------------------
!     8. Maximize Hexagon starting at another feasible point.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-----------------------------------------------'
      write(nOut, *) '8. Hexagon (feasible point, maximized objective'
      write(nOut, *) '-----------------------------------------------'

      call snSet
     &   ( 'Maximize           ',         iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     HexCon, HexObjMax,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *)    ' '
      write(nOut, *)    'snmainb (8) (maximize) finished.'
      write(nOut, *)    'Input  errors =', Errors
      write(nOut, *)    'snOptB INFO   =', INFO
      write(nOut, *)    'Iterations    =', itns
      write(nOut, *)    'nInf          =', nInf
      write(nOut, *)    'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *) 'Obj           =', objAdd + x(n+iObj) + obj
      else
         write(nOut, *) 'Obj           =', objAdd + obj
      end if
      if (INFO .gt. 30) go to 910

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

      end ! program snmain

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObjMin
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     No user-defined storage is used.
!     ==================================================================

      fObj    = - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &          + x(4)*x(9) + x(3)*x(8)

      gObj(1) =   x(7)
      gObj(2) = - x(6)
      gObj(3) = - x(7) + x(8)
      gObj(4) =   x(9)
      gObj(5) = - x(8)
      gObj(6) = - x(2)
      gObj(7) = - x(3) + x(1)
      gObj(8) = - x(5) + x(3)
      gObj(9) =   x(4)

      end ! subroutine HexObjMin

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObjMinFD
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     No derivatives.
!     ==================================================================
      fObj    = - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &          + x(4)*x(9) + x(3)*x(8)

      end ! subroutine HexObjMinFD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObjMax
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     No user-defined storage is used.
!     ==================================================================

      fObj    =   x(2)*x(6) - x(1)*x(7) + x(3)*x(7) + x(5)*x(8)
     &          - x(4)*x(9) - x(3)*x(8)

      gObj(1) = - x(7)
      gObj(2) =   x(6)
      gObj(3) =   x(7) - x(8)
      gObj(4) = - x(9)
      gObj(5) =   x(8)
      gObj(6) =   x(2)
      gObj(7) =   x(3) - x(1)
      gObj(8) =   x(5) - x(3)
      gObj(9) = - x(4)

      end ! subroutine HexObjMax

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObjMaxFD
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     No derivatives supplied.
!     ==================================================================

      fObj    =   x(2)*x(6) - x(1)*x(7) + x(3)*x(7) + x(5)*x(8)
     &          - x(4)*x(9) - x(3)*x(8)

      end ! subroutine HexObjMaxFD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexObjMax2
     &   ( mode, nnObj, x, fObj, gObj, nState,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnObj, nState, lencu, leniu,  lenru, iu(leniu)
      double precision
     &     fObj, x(nnObj), gObj(nnObj), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     THree missing derivatives.
!     ==================================================================

      fObj    =   x(2)*x(6) - x(1)*x(7) + x(3)*x(7) + x(5)*x(8)
     &          - x(4)*x(9) - x(3)*x(8)

      gObj(1) = - x(7)
      gObj(2) =   x(6)
      gObj(3) =   x(7) - x(8)
      gObj(4) = - x(9)
      gObj(5) =   x(8)
      gObj(6) =   x(2)
!     gObj(7) =   x(3) - x(1)
!     gObj(8) =   x(5) - x(3)
!     gObj(9) = - x(4)

      end ! subroutine HexObjMax2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexConFD
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!     No derivatives supplied.
!     ==================================================================
      integer            Out
      parameter         (Out  = 6)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (nState .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  snmainb'
      else  if (nState .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing snmainb'
         return
      end if

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

      end ! subroutine HexConFD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexCon
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon.
!
!     No user-defined storage is used.
!     ==================================================================
      integer            Out
      parameter         (Out  = 6)
      double precision   two
      parameter         (two  = 2.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (nState .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  snmainb'
      else  if (nState .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing snmainb'
         return
      end if

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

!     Nonlinear Jacobian elements for column 1.
!     rows = (1,2,3,4,5).

      gCon( 1) =   two*x(1)
      gCon( 2) = - two*(x(2) - x(1))
      gCon( 3) = - two*(x(3) - x(1))
      gCon( 4) =   two*(x(1) - x(4))
      gCon( 5) =   two*(x(1) - x(5))

!     Nonlinear Jacobian elements for column 2.
!     Rows = (2,6,7,8,9).

      gCon( 6) =   two*(x(2) - x(1))
      gCon( 7) =   two*x(2)
      gCon( 8) = - two*(x(3) - x(2))
      gCon( 9) = - two*(x(4) - x(2))
      gCon(10) =   two*(x(2) - x(5))

!     Nonlinear Jacobian elements for column 3.
!     Rows = (3,7,10,11).

      gCon(11) =   two*(x(3) - x(1))
      gCon(12) =   two*(x(3) - x(2))
      gCon(13) = - two*(x(4) - x(3))
      gCon(14) = - two*(x(5) - x(3))

!     Nonlinear Jacobian elements for column 4.
!     Rows = (4,8,10,12,13).

      gCon(15) = - two*(x(1) - x(4))
      gCon(16) =   two*(x(4) - x(2))
      gCon(17) =   two*(x(4) - x(3))
      gCon(18) =   two*x(4)
      gCon(19) =   two*(x(4) - x(5))

!     Nonlinear Jacobian elements for column 5.
!     Rows = (5,9,11,13,14).

      gCon(20) = - two*(x(1) - x(5))
      gCon(21) = - two*(x(2) - x(5))
      gCon(22) =   two*(x(5) - x(3))
      gCon(23) = - two*(x(4) - x(5))
      gCon(24) =   two*x(5)

!     Nonlinear Jacobian elements for column 6.
!     Rows = (1,2,3,4,5).

      gCon(25) =   two*x(6)
      gCon(26) = - two*(x(7) - x(6))
      gCon(27) =   two*x(6)
      gCon(28) =   two*(x(6) - x(8))
      gCon(29) =   two*(x(6) - x(9))

!     Nonlinear Jacobian elements for column 7.
!     Rows = (2,6,7,8,9).

      gCon(30) =   two*(x(7) - x(6))
      gCon(31) =   two*x(7)
      gCon(32) =   two*x(7)
      gCon(33) = - two*(x(8) - x(7))
      gCon(34) =   two*(x(7) - x(9))

!     Nonlinear Jacobian elements for column 8.
!     Rows = (4,8,10,12,13).

      gCon(35) = - two*(x(6) - x(8))
      gCon(36) =   two*(x(8) - x(7))
      gCon(37) =   two*x(8)
      gCon(38) =   two*x(8)
      gCon(39) = - two*(x(9) - x(8))

!     Nonlinear Jacobian elements for column 9.
!     Rows = (5,9,11,13,14).

      gCon(40) = - two*(x(6) - x(9))
      gCon(41) = - two*(x(7) - x(9))
      gCon(42) =   two*x(9)
      gCon(43) =   two*(x(9) - x(8))
      gCon(44) =   two*x(9)

      end ! subroutine HexCon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexCon1
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character*8
     &     cu(lencu)

!     ==================================================================
!     Problem Hexagon1.
!
!     Five missing constraint gradients.
!     ==================================================================
      integer            Out
      parameter         (Out  = 6)
      double precision   two
      parameter         (two  = 2.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (nState .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  snmainb'
      else  if (nState .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing snmainb'
         return
      end if

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

!     Nonlinear Jacobian elements for column 1.
!     rows = (1,2,3,4,5).

      gCon( 1) =   two*x(1)
      gCon( 2) = - two*(x(2) - x(1))
      gCon( 3) = - two*(x(3) - x(1))
      gCon( 4) =   two*(x(1) - x(4))
      gCon( 5) =   two*(x(1) - x(5))

!     Nonlinear Jacobian elements for column 2.
!     Rows = (2,6,7,8,9).

      gCon( 6) =   two*(x(2) - x(1))
      gCon( 7) =   two*x(2)
      gCon( 8) = - two*(x(3) - x(2))
      gCon( 9) = - two*(x(4) - x(2))
      gCon(10) =   two*(x(2) - x(5))

!     Nonlinear Jacobian elements for column 3.
!     Rows = (3,7,10,11).

!     gCon(11) =   two*(x(3) - x(1))
      gCon(12) =   two*(x(3) - x(2))
      gCon(13) = - two*(x(4) - x(3))
      gCon(14) = - two*(x(5) - x(3))

!     Nonlinear Jacobian elements for column 4.
!     Rows = (4,8,10,12,13).

      gCon(15) = - two*(x(1) - x(4))
      gCon(16) =   two*(x(4) - x(2))
      gCon(17) =   two*(x(4) - x(3))
      gCon(18) =   two*x(4)
      gCon(19) =   two*(x(4) - x(5))

!     Nonlinear Jacobian elements for column 5.
!     Rows = (5,9,11,13,14).

!     gCon(20) = - two*(x(1) - x(5))
      gCon(21) = - two*(x(2) - x(5))
      gCon(22) =   two*(x(5) - x(3))
      gCon(23) = - two*(x(4) - x(5))
      gCon(24) =   two*x(5)

!     Nonlinear Jacobian elements for column 6.
!     Rows = (1,2,3,4,5).

!     gCon(25) =   two*x(6)
!     gCon(26) = - two*(x(7) - x(6))
!     gCon(27) =   two*x(6)
      gCon(28) =   two*(x(6) - x(8))
      gCon(29) =   two*(x(6) - x(9))

!     Nonlinear Jacobian elements for column 7.
!     Rows = (2,6,7,8,9).

      gCon(30) =   two*(x(7) - x(6))
      gCon(31) =   two*x(7)
      gCon(32) =   two*x(7)
      gCon(33) = - two*(x(8) - x(7))
      gCon(34) =   two*(x(7) - x(9))

!     Nonlinear Jacobian elements for column 8.
!     Rows = (4,8,10,12,13).

      gCon(35) = - two*(x(6) - x(8))
      gCon(36) =   two*(x(8) - x(7))
      gCon(37) =   two*x(8)
      gCon(38) =   two*x(8)
      gCon(39) = - two*(x(9) - x(8))

!     Nonlinear Jacobian elements for column 9.
!     Rows = (5,9,11,13,14).

      gCon(40) = - two*(x(6) - x(9))
      gCon(41) = - two*(x(7) - x(9))
      gCon(42) =   two*x(9)
      gCon(43) =   two*(x(9) - x(8))
      gCon(44) =   two*x(9)

      end ! subroutine HexCon1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine HexDat
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, neJ, nnCon, nnObj, nnJac, iObj, ObjAdd,
     &     Jcol, indJ, locJ, bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxneJ, INFO, m, n, neJ, nnCon, nnObj, nnJac,
     &     iObj, indJ(maxneJ) , hs(maxn+maxm), locJ(maxn+1)
      double precision
     &     ObjAdd, Jcol(maxneJ) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character*8
     &     Prob

!     ------------------------------------------------------------------
!     HexDat generates data for the Hexagon problem.
!     The constraints take the form
!              c(x) + A*x - s = 0,
!     where the Jacobian for c(x) + Ax is stored in Jcol(*), and any
!     terms coming from c(x) are in the TOP LEFT-HAND CORNER of Jcol(*),
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
!     maxm, maxn, maxneJ are upper limits on m, n, ne.
!
!     On exit,
!     INFO    is 0 if there is enough storage, 1 otherwise.
!     m       is the number of nonlinear and linear constraints.
!     n       is the number of variables.
!     neJ     is the number of nonzeros in Jcol(*).
!     nnCon   is the number of nonlinear constraints (they come first).
!     nnObj   is the number of nonlinear objective variables.
!     nnJac   is the number of nonlinear Jacobian variables.
!     Jcol    is the constraint matrix (Jacobian), stored column-wise.
!     indJ    is the list of row indices for each nonzero in Jcol(*).
!     locJ    is a set of pointers to the start of each column of Jcol.
!     bl      is the lower bounds on x and s.
!     bu      is the upper bounds on x and s.
!     hs(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
!     x (1:n) is a set of initial values for x.
!     pi(1:m) is a set of initial values for the dual variables pi.
!
!     24 Dec 1997: First version of HexDat.
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

      Prob   = 'snmainb '

      neJ    = 52
      n      =  9
      m      = 18

      nnCon  = 14
      nnJac  =  n
      nnObj  =  n

!     Check if there is enough storage.

      INFO = 0
      if (m     .gt.   maxm) INFO = 1
      if (n     .gt.   maxn) INFO = 1
      if (neJ   .gt. maxneJ) INFO = 1
      if (INFO  .gt.      0) return

!     --------------------------------------
!     Set up the list of row indices in indJ.
!     --------------------------------------
!     Column  1
!     Nonlinear elements in rows (1, 2, 3, 4, 5)  first.

      locJ( 1) =  1

      indJ( 1) =  1
      indJ( 2) =  2
      indJ( 3) =  3
      indJ( 4) =  4
      indJ( 5) =  5

      Jcol( 1) =  zero
      Jcol( 2) =  zero
      Jcol( 3) =  zero
      Jcol( 4) =  zero
      Jcol( 5) =  zero

!     Column 1.
!     Linear element in row 6 next.

      indJ( 6) = 15

      Jcol( 6)  = -one

!     Column 2.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).

      locJ( 2) =  7

      indJ( 7) =  2
      indJ( 8) =  6
      indJ( 9) =  7
      indJ(10) =  8
      indJ(11) =  9

      Jcol( 7) =  zero
      Jcol( 8) =  zero
      Jcol( 9) =  zero
      Jcol(10) =  zero
      Jcol(11) =  zero

!     Column 2.
!     Linear elements in rows (15,16).

      indJ(12) = 15
      indJ(13) = 16

      Jcol(12) =  one
      Jcol(13) = -one

!     Column 3.
!     Nonlinear elements in rows (3, 7, 10, 11).

      locJ( 3) =  14

      indJ(14) =  3
      indJ(15) =  7
      indJ(16) = 10
      indJ(17) = 11

      Jcol(14) =  zero
      Jcol(15) =  zero
      Jcol(16) =  zero
      Jcol(17) =  zero

!     Column 3.
!     Linear elements in rows (16, 17).

      indJ(18) = 16
      indJ(19) = 17

      Jcol(18) =  one
      Jcol(19) =  one

!     Column 4.
!     Nonlinear elements in rows (20, 21, 22, 23, 24).

      locJ( 4) = 20

      indJ(20) =  4
      indJ(21) =  8
      indJ(22) = 10
      indJ(23) = 12
      indJ(24) = 13

      Jcol(20) =  zero
      Jcol(21) =  zero
      Jcol(22) =  zero
      Jcol(23) =  zero
      Jcol(24) =  zero

!     Column 4.
!     Linear elements in rows (17, 18).

      indJ(25) = 17
      indJ(26) = 18

      Jcol(25) = -one
      Jcol(26) =  one

!     Column 5.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).

      locJ( 5) = 27

      indJ(27) =  5
      indJ(28) =  9
      indJ(29) = 11
      indJ(30) = 13
      indJ(31) = 14

      Jcol(27) =  zero
      Jcol(28) =  zero
      Jcol(29) =  zero
      Jcol(30) =  zero
      Jcol(31) =  zero

!     Column 5.
!     Linear element in row 18.

      indJ(32) = 18

      Jcol(32) = -one

!     Column 6.
!     Nonlinear elements in rows (1, 2, 3, 4, 5, 6).

      locJ(6)  = 33

      indJ(33) =  1
      indJ(34) =  2
      indJ(35) =  3
      indJ(36) =  4
      indJ(37) =  5

      Jcol(33) =  zero
      Jcol(34) =  zero
      Jcol(35) =  zero
      Jcol(36) =  zero
      Jcol(37) =  zero

!     Column 7.
!     Nonlinear elements in rows (2, 6, 7, 8, 9).

      locJ(7)  =  38

      indJ(38) =  2
      indJ(39) =  6
      indJ(40) =  7
      indJ(41) =  8
      indJ(42) =  9

      Jcol(38) =  zero
      Jcol(39) =  zero
      Jcol(40) =  zero
      Jcol(41) =  zero
      Jcol(42) =  zero

!     Column 8.
!     Nonlinear elements in rows (4, 8, 10, 12, 13).

      locJ(8)  =  43

      indJ(43) =  4
      indJ(44) =  8
      indJ(45) = 10
      indJ(46) = 12
      indJ(47) = 13

      Jcol(43) =  zero
      Jcol(44) =  zero
      Jcol(45) =  zero
      Jcol(46) =  zero
      Jcol(47) =  zero

!     Column 9.
!     Nonlinear elements in rows (5, 9, 11, 13, 14).

      locJ(9)  =  48

      indJ(48) =  5
      indJ(49) =  9
      indJ(50) = 11
      indJ(51) = 13
      indJ(52) = 14

      Jcol(48) =  zero
      Jcol(49) =  zero
      Jcol(50) =  zero
      Jcol(51) =  zero
      Jcol(52) =  zero

!     Don't forget to finish off  locJ.
!     This is crucial.

      locJ(10) =  neJ + 1

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

      end ! subroutine HexDat
