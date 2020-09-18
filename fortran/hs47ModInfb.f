!     ------------------------------------------------------------------
!     File hs47ModInfb.f : solves a modified version of hs47.
!     This is problem is hs47Modc with the bounds changed to give
!     an infeasible problem:
!
!     (1) the lower bounds on the two general linear constraints are
!         switched.
!     (2) lower bounds of 1.0 are imposed on x(1:n).
!
!     12 October 2014: First  version.
!     ------------------------------------------------------------------
      program
     &     hs47ModInf

      implicit
     &     none
      integer
     &     maxm, maxn, maxneJ, nNames
      parameter
     &     ( maxm   = 1000,
     &       maxn   = 1000,
     &       maxneJ = 3000,
     &       nNames = 1 )

      character
     &     Prob*8, Names(nNames)*8
      integer
     &     indJ(maxneJ), hs(maxn+maxm)
      integer
     &     locJ(maxn+1)
      double precision
     &     Jcol(maxneJ), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)     , rc(maxn+maxm)
      integer
     &     lenru, leniu, lencu, lenrw, leniw, lencw
!     ------------------------------------------------------------------
!     SNOPT workspace

      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      parameter          (  lencw =   600)
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
     &     errors, iPrt, iSum, INFO, iObj, iPrint, iSpecs, iSumm,
     &     itnlim, lunit, m, MjrLim, maxiu, maxru, mincw, miniw, minrw,
     &     n, neJ, nInf, nlocJ, nnCon, nnJac, nnObj,
     &     nOut, nS
      double precision
     &     Obj, objAdd, sInf
      external
     &     hs47Con, hs47Obj
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

         lfile = 'hs47ModInfb.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'hs47ModInfb.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      else

!        VMS  systems.  Define units for the Specs and print files.

         lunit = iSpecs
         open( lunit, status='OLD',     err=900 )
         lunit = iPrint
         open( lunit, status='UNKNOWN', err=900 )
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

      if (INFO .ne. 101) then
         go to 990
      end if

!     ------------------------------------------------------------------
!     1. Solve the modified hs47 with an infeasible constraint.
!        This is a test of implicit elastic mode.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------'
      write(nOut, *) '1. Infeasible Modified hs47'
      write(nOut, *) '   using implicit elastic mode'
      write(nOut, *) '------------------------------'

      call hs47ModInfData
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     neJ, Jcol, nlocJ, indJ, locJ,
     &     bl, bu, hs, x, pi )

!     ------------------------------------------------------------------
!     Specify any options not set in the Specs file.
!     iPrt and iSum may refer to the print and summary files.
!     Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      errors = 0
      iPrt   = 0
      iSum   = 0

      itnlim = 2500
      call snSeti
     &   ( 'Iterations        ', itnlim, iPrt, iSum, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      MjrLim = 250
      call snSeti
     &   ( 'Major Iterations', MjrLim, iPrt, iSum, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &   ( 'Cold', m, n, neJ, nNames,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     hs47Con, hs47Obj,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfb finished.'
      write(nOut, *) 'INFO  =', INFO
      write(nOut, *) 'nInf  =', nInf
      write(nOut, *) 'sInf  =', sInf
      write(nOut, *) 'Obj   =', Obj
      if (INFO .ge. 30) go to 910

!     ------------------------------------------------------------------
!     2. Solve modified hs47 with explicit elastic variables.
!        The elastic  objective term is defined by a free row of Jcon.
!        Number of nonlinear objective variables = 5.
!     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '------------------------------------------------'
      write(nOut, *) '2. Infeas Modified hs47 with explicit elastics'
      write(nOut, *) '   Number of nonlinear objective variables = 5  '
      write(nOut, *) '------------------------------------------------'

      call hs47ModInfExData
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     neJ, Jcol, nlocJ, indJ, locJ,
     &     bl, bu, hs, x, pi )

      call snseti
     &   ( 'Verify level = ',      3, iPrt, iSum, errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptB
     &  ( 'Cold', m, n, neJ, nNames,
     &     nnCon, nnObj, nnJac,
     &     iObj, objAdd, Prob,
     &     hs47Con, hs47Obj,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'hs47ModInfc (Ex) finished.'
      write(nOut, *) 'INFO  =', INFO
      write(nOut, *) 'nInf  =', nInf
      write(nOut, *) 'sInf  =', sInf
      write(nOut, *) 'Obj   =', Obj
      if (INFO .ge. 30) go to 910

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

      end ! program hs47ModInf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ModInfData
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     neJ, Jcol, nlocJ, indJ, locJ,
     &     bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxneJ, INFO, m, n, neJ, nlocJ,
     &     nnCon, nnObj, nnJac, iObj,
     &     indJ(maxneJ), hs(maxn+maxm), locJ(maxn+1)
      double precision
     &     objAdd, Jcol(maxneJ) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character
     &     Prob*8

!     ------------------------------------------------------------------
!     hs47ModInfData generates data for a modified version of HS47.
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to
!                  x(1)      + x(2)**2 + x(3)**3                = 3
!                              x(2)    - x(3)**2        + x(5)  = 1
!                  x(1)*x(4)                                    = 1
!                  x(1)      + x(2)                            >= 3
!                                        x(3)    + x(4) + x(5) >= 1
!                  x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     The problem is input in generic form:
!
!     Minimize               f(x)
!
!     subject to
!                       (      x     )
!                bl le  ( fCon(x)+Jx ) le bu
!                       (     Ax     )
!
!     where the Jacobian for fCon(x) + J x is stored in JCon(ldJ,*),
!     with with dimensions  mNCon by n.  The elements of
!     the Jacobian of fCon comprise the first nnJac columns of JCon.
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
!     19 Sep 2014: First version of hs47ModInfData.
!     ------------------------------------------------------------------
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision   infBnd
      parameter         (infBnd   = 1.0d+20)
!     ------------------------------------------------------------------
!     Give a name to the Problem.

      Prob   = 'HS47MOD '

      neJ    =  13
      n      =   5
      m      =   5

      nnCon  =   3
      nnJac  =   4
      nnObj  =   5

      iObj   = 0               ! No linear objective term.
      objAdd = 0.0d+0

      nlocJ  = n   + 1

!     Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm  ) INFO = 1
      if (n     .gt. maxn  ) INFO = 1
      if (neJ   .gt. maxneJ) INFO = 1
      if (INFO  .gt.   0  ) return

!     -------------------------------------
!     Set up the list of row indices in indJ.
!     -------------------------------------
!     Column  1.
!     Elements in nonlinear rows 1,3 first.

      neJ = 0
      locJ( 1) =  1

      neJ       = neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  1.0d+0          ! Constant

      neJ       = neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) =  0.0d+0

!     Column 1.
!     Linear element in row 4 next.

      neJ       = neJ + 1
      indJ(neJ) = 4
      Jcol(neJ) = 1.0d+0

!     -------------------------------------------
!     Column 2.
!     Elements in nonlinear rows  (1, 2).
!     -------------------------------------------
      locJ( 2)  =  neJ + 1

      neJ       = neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  0.0d+0

      neJ       = neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  1.0d+0          ! Constant

!     Column 2.
!     Linear element in row (4).

      neJ       = neJ + 1
      indJ(neJ) = 4
      Jcol(neJ) = 1.0d+0

!     ------------------------------------------
!     Column 3.
!     Elements in nonlinear rows  (1, 2).
!     ------------------------------------------
      locJ( 3) = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  0.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  0.0d+0

!     Column3.
!     Linear element in rows (5).

      neJ       =  neJ + 1
      indJ(neJ) =  5
      Jcol(neJ) =  1.0d+0

!     -------------------------------------------
!     Column 4.
!     Nonlinear element in row 3.

      locJ( 4) = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) =  0.0d+0

!     Column 4.
!     Linear element in row (5).

      neJ       = neJ + 1
      indJ(neJ) = 5
      Jcol(neJ) = 1.0d+0

!     -------------------------------------------
!     Column 5.
!     Linear elements in rows (2,5).

      locJ(5) = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  5
      Jcol(neJ) =  1.0d+0

!     -----------------------------------------
!     Don't forget to finish off locJ.
!     This is crucial.
!     -----------------------------------------
      locJ(6) = neJ + 1

!     ------------------------------------------------------------------
!     Constraint ranges
!     ------------------------------------------------------------------
!     Nonlinear constraints first.

      bl(n+1) = 3.0d+0
      bu(n+1) = 3.0d+0

      bl(n+2) = 1.0d+0
      bu(n+2) = 1.0d+0

      bl(n+3) = 1.0d+0
      bu(n+3) = 1.0d+0

!     Followed by the linear constraints.

      bl(n+4) = 3.0d+0
      bu(n+4) = infBnd

      bl(n+5) = 1.0d+0
      bu(n+5) = infBnd

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, n
         bl(j) =  1.0d+0        ! was  -infBnd
         bu(j) =  infBnd
      end do

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  2.0d+0
      x(2)   =  sqrt(2.0d+0)
      x(3)   = -1.0d+0
      x(4)   =  0.5d+0
      x(5)   =  2.0d+0 - sqrt(2.0d+0)

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine hs47ModInfData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47ModInfExData
     &   ( Prob, maxm, maxn, maxneJ, INFO,
     &     m, n, nnCon, nnObj, nnJac, iObj, objAdd,
     &     neJ, Jcol, nlocJ, indJ, locJ,
     &     bl, bu, hs, x, pi )

      implicit
     &     none
      integer
     &     maxm, maxn, maxneJ, INFO, m, n, neJ, nlocJ,
     &     nnCon, nnObj, nnJac, iObj,
     &     indJ(maxneJ), hs(maxn+maxm), locJ(maxn+1)
      double precision
     &     objAdd, Jcol(maxneJ) , bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm)
      character
     &     Prob*8

!     ------------------------------------------------------------------
!     hs47ModInfExData generates the data for the modified HS47 with
!     explicit elastic variables. The constant elastic graient term is
!     held as a free row of A.
!
!     The modified version of hs47 is:
!
!     Minimize    (x(1)-x(2))**2 + (x(2)-x(3))**3 + (x(3)-x(4))**4
!                                                 + (x(4)-x(5))**4
!
!     subject to
!                  x(1)      + x(2)**2 + x(3)**3                = 3
!                              x(2)    - x(3)**2        + x(5)  = 1
!                  x(1)*x(4)                                    = 1
!                  x(1)      + x(2)                            >= 3
!                                        x(3)    + x(4) + x(5) >= 1
!                  x(1),       x(2),     x(3),    x(4),   x(5) >= 1
!
!     Given a paenalty parameter penParm, explicit elastic variables
!     u(1:3), v(1:3) are added to give the elastic problem:
!
!     Minimize    f(x) + penParm*(u + v)
!
!     subject to
!                   (             x         )
!                   (             u         )
!            bl le  (             v         ) le bu
!                   ( fCon(x) + J x - u + v )
!                   (            Ax         )
!
!     where the Jacobian for fCon(x) + J x - u + v is stored in
!     JCon(ldJ,*), with with dimensions  mNCon by n.  The elements of
!     the Jacobian of fCon comprise the first nnJac columns JCon.
!
!     The elastic variables u and v are held in x(6:11).
!
!
!     On entry,
!     maxm, maxn, maxne are upper limits on m, n, ne.
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
!     12 Aug 2005: First version of hs47ModExDat.
!     ------------------------------------------------------------------
      integer
     &     i, j
!     ------------------------------------------------------------------
      double precision   infBnd
      parameter         (infBnd = 1.0d+20)

      double precision   wtInf
      parameter         (wtInf  = 1.0d+4)
!     ------------------------------------------------------------------
!     Give a name to the Problem.

      Prob   = 'hs47In'

      neJ    = 25               ! 13 + 12 elastic grads
      n      = 11               ! x(1:5), u(1:3) v(1:3)
      m      =  6               ! m = 5 + free row

      nnCon  =  3
      nnJac  =  4
      nnObj  =  5
      iObj   =  6

      nlocJ  = n   + 1

      objAdd = 0.0d+0           ! No constant objective term.

!     Check if there is enough storage.

      INFO   = 0
      if (m     .gt. maxm  ) INFO = 1
      if (n     .gt. maxn  ) INFO = 1
      if (neJ   .gt. maxneJ) INFO = 1
      if (INFO  .gt.   0   ) return

!     -------------------------------------
!     Set up the list of row indices in indJ.
!     -------------------------------------
!     Column  1.
!     Elements in nonlinear rows 1,3 first.

      neJ = 0
      locJ( 1)  =  1

      neJ       = neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  1.0d+0          ! Constant

      neJ       = neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) =  0.0d+0

!     Linear element in row 4 next.

      neJ       = neJ + 1
      indJ(neJ) =  4
      Jcol(neJ) =  1.0d+0

!     -------------------------------------------
!     Column 2.
!     Elements in nonlinear rows  (1, 2).
!     -------------------------------------------
      locJ( 2)  =  neJ + 1

      neJ       = neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  0.0d+0

      neJ       = neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  1.0d+0          ! Constant

!     Linear element in row (4).

      neJ       = neJ + 1
      indJ(neJ) = 4
      Jcol(neJ) = 1.0d+0

!     ------------------------------------------
!     Column 3.
!     Elements in nonlinear rows  (1, 2).
!     ------------------------------------------
      locJ( 3)  = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) =  0.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  0.0d+0

!     Linear element in rows (5).

      neJ       =  neJ + 1
      indJ(neJ) =  5
      Jcol(neJ) =  1.0d+0

!     -------------------------------------------
!     Column 4.
!     Nonlinear element in row 3.

      locJ( 4) = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) =  0.0d+0

!     Column 4.

      neJ       = neJ + 1
      indJ(neJ) = 5
      Jcol(neJ) = 1.0d+0

!     -------------------------------------------
!     Column 5.
!     Linear elements in rows (2,5).

      locJ(5) = neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) =  1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  5
      Jcol(neJ) =  1.0d+0

!     -------------------------------------------
!     Column 6.

      locJ(6)   =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) = -1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -------------------------------------------
!     Column 7.

      locJ(7)   =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) = -1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -------------------------------------------
!     Column 8.

      locJ(8)   =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) = -1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -------------------------------------------
!     Column 9.

      locJ(9)   =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  1
      Jcol(neJ) = +1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -------------------------------------------
!     Column 10.

      locJ(10)  =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  2
      Jcol(neJ) = +1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -------------------------------------------
!     Column 11.

      locJ(11)  =  neJ + 1

      neJ       =  neJ + 1
      indJ(neJ) =  3
      Jcol(neJ) = +1.0d+0

      neJ       =  neJ + 1
      indJ(neJ) =  iObj
      Jcol(neJ) =  wtInf

!     -----------------------------------------
!     Don't forget to finish off locJ.
!     This is crucial.
!     -----------------------------------------
      locJ(n+1) = neJ + 1

!     ------------------------------------------------------------------
!     Constraint ranges
!     ------------------------------------------------------------------
!     Nonlinear constraints first.

      bl(n+1) =  3.0d+0
      bu(n+1) =  3.0d+0

      bl(n+2) =  1.0d+0
      bu(n+2) =  1.0d+0

      bl(n+3) =  1.0d+0
      bu(n+3) =  1.0d+0

!     Linear constraints next.

      bl(n+4) =  3.0d+0
      bu(n+4) =  infBnd

      bl(n+5) =  1.0d+0
      bu(n+5) = +infBnd

!     Objective free row

      bl(n+6) = -infBnd
      bu(n+6) = +infBnd

!     ------------------------------------------------------------------
!     Variable ranges
!     ------------------------------------------------------------------
      do j = 1, nnObj
         bl(j) =  1.0d+0        ! was  -infBnd
         bu(j) =  infBnd
      end do

      do j = nnObj+1, n
         bl(j) =  0.0d+0
         bu(j) =  infBnd
      end do

!     ------------------------------------------------------------------
!     Initialize x, hs and pi.
!     Set the initial value and status of each variable.
!     For want of something better to do, make the variables x(1:n)
!     temporarily fixed at their current values.
!     The crash can set the rest.
!     ------------------------------------------------------------------
      x(1)   =  2.0d+0
      x(2)   =  sqrt(2.0d+0)
      x(3)   = -1.0d+0
      x(4)   =  0.5d+0
      x(5)   =  2.0d+0 - sqrt(2.0d+0)

      do j = nnObj+1, n
         x(j) = 0.0d+0
      end do

      do j = 1, n
         hs(j)  = 0
      end do

      do i = 1, m
         pi(i)  = 0.0d+0
      end do

      end ! subroutine hs47ModInfExData

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47Obj
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
!     Problem Hexagon.
!     No user-defined storage is used.
!     ==================================================================
      double precision   two
      parameter         (two    = 2.0d+0)
      double precision   three
      parameter         (three  = 3.0d+0)
      double precision   four
      parameter         (four   = 4.0d+0)
!     ------------------------------------------------------------------
      fObj    = (x(1)-x(2))**2 + (x(2)-x(3))**3
     &            + (x(3)-x(5))**4 + (x(5)-x(4))**4

      gObj(1) =  two*( x(1)-x(2) )
      gObj(2) = -two*(x(1)-x(2)) + three*(x(2)-x(3))**2
      gObj(3) =  four*(x(3)-x(5))**3 - three*(x(2)-x(3))**2
      gObj(4) = -four*(x(5)-x(4))**3
      gObj(5) =  four*( (x(5)-x(4))**3 - (x(3)-x(5))**3 )

      end ! subroutine hs47Obj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hs47Con
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
      double precision   two
      parameter         (two    = 2.0d+0)
      double precision   one
      parameter         (one    = 1.0d+0)
      double precision   three
      parameter         (three  = 3.0d+0)
      double precision   four
      parameter         (four   = 4.0d+0)
!     ------------------------------------------------------------------
!     Print something on the first and last entry.

      if (      State .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  hs47ModInfb'
      else  if (State .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing hs47ModInfb'
      end if

!     ------------------------------------------------------------------
!     Constraints and Objective.
!     ------------------------------------------------------------------
      fCon( 1) = x(1) + x(2)**2 + x(3)**3
      fCon( 2) = x(2) - x(3)**2
      fCon( 3) = x(1)*x(4)

      neG = 0

!     -------------------------------------------
!     Nonlinear elements for column 1 (g=df/dx1).
!     g rows = (1).
!     -------------------------------------------
      neG = neG + 1
      gCon(neG) =  one          ! row 1
      neG = neG + 1
      gCon(neG) =  x(4)         ! row 3

!     -------------------------------------------
!     Nonlinear elements for column 2 (g=df/dx2).
!     g Rows = (1,2).
!     -------------------------------------------
      neG = neG + 1
      gCon(neG) =  two*x(2)     ! row 1
      neG = neG + 1
      gCon(neG) =  one          ! row 2

!     -------------------------------------------
!     Nonlinear elements for column 3 (g=df/dx3).
!     g Rows = (3,7,10,11,15).
!     -----------------------------------------
      neG = neG + 1
      gCon(neG) =  three*x(3)**2 ! row 1
      neG = neG + 1
      gCon(neG) =  -two*x(3)     ! row 2

!     -------------------------------------------
!     Nonlinear elements for column 4 (g=df/dx4).
!     g Rows = (3,7,10,11,15).
!     -----------------------------------------
      neG = neG + 1
      gCon(neG) =  x(1)         ! row 3


      end ! subroutine hs47Con

