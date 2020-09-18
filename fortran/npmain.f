*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File npmain.f
*
*     Sample program for NPOPT Version 7.1  December 2003.
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      program
     &     npmain
      implicit
     &     none

*     ==================================================================
*     Set the declared array dimensions.
*     ldA    = the declared row dimension of  a.
*     ldg    = the declared row dimension of  gCon.
*     ldH    = the declared row dimension of  Hess.
*     maxn   = maximum no. of variables allowed for.
*     maxbnd = maximum no. of variables + linear & nonlinear constrnts.
*     liwork = the length of the integer work array.
*     lwork  = the length of the double precision work array.
*     ==================================================================
      integer
     &     ldA, ldg, ldH, liwork, lwork, maxbnd, maxn
      parameter
     &    (ldA    =  5, ldg    =   20, ldH   =   10,
     &     maxn   =  9, liwork = 2500, lwork = 2500,
     &     maxbnd =  maxn + ldA + ldg)

      integer
     &     Errors, i, INFO, iOptns, iP, iPrint, iS, iSumm, iter, j,
     &     lunit, n, nbnd, nclin, ncnln, nOut,
     &     istate(maxbnd), iwork(liwork)
      double precision
     &     fObj, InfBnd, A(ldA,maxn), bl(maxbnd), bu(maxbnd),
     &     fCon(ldg), gCon(ldg,maxn), cMul(maxbnd),
     &     gObj(maxn), Hess(ldH,maxn), x(maxn), work(lwork)
      logical
     &     byname, byunit
      character
     &     lFile*20
      external
     &     fnobj1, fnobj2, fncon1, fncon2

*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
*     Assign file numbers and open files by various means.
*     (Some systems don't need explicit open statements.)
*     iOptns = unit number for the Options file.
*     iPrint = unit number for the Print file.
*     iSumm  = unit number for the Summary file.
*     ------------------------------------------------------------------
      iOptns =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6
      byname = .true.
      byunit = .false.

      if ( byname ) then
         lFile = 'npmain.spc'
         open( iOptns, file=lFile, status='OLD',     err=800 )

         lFile = 'npmain.out'
         open( iPrint, file=lFile, status='UNKNOWN', err=800 )

      else if ( byunit ) then
         lUnit = iOptns
         open( lUnit, status='OLD',     err=900 )

         lUnit = iPrint
         open( lUnit, status='UNKNOWN', err=900 )
      end if

*     =============================================================
*     Set the actual problem dimensions.
*     n      = the number of variables.
*     nclin  = the number of general linear constraints (may be 0).
*     ncnln  = the number of nonlinear constraints (may be 0).
*     =============================================================
      n      = 9
      nclin  = 4
      ncnln  = 14
      nbnd   = n + nclin + ncnln

*     ------------------------------------------------------------------
*     Assign the data arrays.
*     iPrint = the unit number for printing.
*     iOptns = the unit number for reading the options file.
*     bounds  .ge.    infBnd  will be treated as plus  infinity.
*     bounds  .le.  - infBnd  will be treated as minus infinity.
*     A      = the linear constraint matrix.
*     bl     = the lower bounds on  x,  A'x  and  fCon(x).
*     bu     = the upper bounds on  x,  A'x  and  fCon(x).
*     x      = the initial estimate of the solution.
*     ------------------------------------------------------------------

*     Set the matrix  A.

      do j = 1, n
         do  i = 1, nclin
            A(i,j) = zero
         end do
      end do
      A(1,1) = -one
      A(1,2) =  one
      A(2,2) = -one
      A(2,3) =  one
      A(3,3) =  one
      A(3,4) = -one
      A(4,4) =  one
      A(4,5) = -one

*     Set the bounds.

      infBnd =  1.0d+21

      do j =  1, nbnd
         bl(j) = -infBnd
         bu(j) =  infBnd
      end do
      bl(1)  =  zero
      bl(3)  = -one
      bl(5)  =  zero
      bl(6)  =  zero
      bl(7)  =  zero

      bu(3)  =  one
      bu(8)  =  zero
      bu(9)  =  zero

*     Set lower bounds of zero for all four linear constraints.

      do j =  n+1, n+nclin
         bl(j) =  zero
      end do

*     Set upper bounds of one for all 14 nonlinear constraints.

      do j =  n+nclin+1, nbnd
         bu(j) =  one
      end do

*     Set the initial estimate of  x.

      x(1)   =  .1d+0
      x(2)   =  .125d+0
      x(3)   =  .666666d+0
      x(4)   =  .142857d+0
      x(5)   =  .111111d+0
      x(6)   =  .2d+0
      x(7)   =  .25d+0
      x(8)   = -.2d+0
      x(9)   = -.25d+0

*     ------------------------------------------------------------------
*     First,  npInit MUST be called to initialize optional parameters
*     to their default values.
*     The Print file   will be on unit iPrint.
*     The Summary file will be on unit iSumm (typically the screen).
*     ------------------------------------------------------------------
      call npInit
     &   ( iPrint, iSumm, iwork, liwork, work, lwork )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      call npSpec
     &   ( iOptns, INFO, iwork, liwork, work, lwork )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         stop
      end if

*     ------------------------------------------------------------------
*     Set a few options in-line.
*     iP and iS may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      Errors = 0

      iP     =  0
      iS     =  0
      call npSetr
     &   ( 'Infinite Bound =', infBnd, iP, iS, Errors,
     &     iwork, liwork, work, lwork )

*     ------------------------------------------------------------------
*     Solve the problem.
*     ------------------------------------------------------------------
      call npOpt
     &   ( n, nclin, ncnln, ldA, ldg, ldH,
     &     A, bl, bu,
     &     fncon1, fnobj1,
     &     INFO, iter, istate,
     &     fCon, gCon, cMul, fObj, gObj, Hess, x,
     &     iwork, liwork, work, lwork )

      if (INFO .gt. 30) go to 999

*     ------------------------------------------------------------------
*     The following is for illustrative purposes only.
*     A second run solves the same problem,  but defines the objective
*     and constraints via the subroutines fnobj2 and fncon2.  Some
*     objective derivatives and the constant Jacobian elements are not
*     supplied.
*     We do a warm start using
*              istate    (the working set)
*              cMul    (the Lagrange multipliers)
*              Hess      (the Hessian approximation)
*     from the previous run, but with a slightly perturbed starting
*     point.
*     ------------------------------------------------------------------
      do j = 1, n
         x(j) = x(j) + 0.01d+0
      end do

      call npSet
     &   ( ' ',                    iP, iS, Errors,
     &     iwork, liwork, work, lwork )
      call npSet
     &   ( 'Defaults',             iP, iS, Errors,
     &     iwork, liwork, work, lwork )
      call npSeti
     &   ( 'Derivative Level',  0, iP, iS, Errors,
     &     iwork, liwork, work, lwork )
      call npSet
     &   ( 'Verify no',            iP, iS, Errors,
     &     iwork, liwork, work, lwork )
      call npSet
     &   ( 'Warm Start',           iP, iS, Errors,
     &     iwork, liwork, work, lwork )
      call npSeti
     &   ( 'Major Iterations', 20, iP, iS, Errors,
     &     iwork, liwork, work, lwork )

*     ------------------------------------------------------------------
*     Go for it...
*     ------------------------------------------------------------------
      call npOpt
     &   ( n, nclin, ncnln, ldA, ldg, ldH,
     &     A, bl, bu,
     &     fncon2, fnobj2,
     &     INFO, iter, istate,
     &     fCon, gCon, cMul, fObj, gObj, Hess, x,
     &     iwork, liwork, work, lwork )

      write(nOut, *) ' '
      write(nOut, *) 'npmain finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'npOpt INFO    =', INFO
      write(nOut, *) 'Obj           =', fObj
      if (INFO .gt. 30) go to 999

      stop

*     Error conditions.

  800 write(nOut , 4000) 'Error while opening file', lfile
      stop

  900 write(nOut , 4010) 'Error while opening unit', lunit
      stop

  999 write(iPrint, 3010)
      stop

 3010 format(/ ' npOpt  terminated with error condition')
 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! program npmain

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine fnobj1
     &     ( mode, n, x, fObj, gObj, nState )

      implicit
     &     none
      integer
     &     mode, n, nState
      double precision
     &     fObj, x(n), gObj(n)

*     ==================================================================
*     fnobj1  computes the value and first derivatives of the nonlinear
*     objective function.
*     ==================================================================
      fObj   = - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &         + x(4)*x(9) + x(3)*x(8)

      gObj(1) =   x(7)
      gObj(2) = - x(6)
      gObj(3) = - x(7) + x(8)
      gObj(4) =   x(9)
      gObj(5) = - x(8)
      gObj(6) = - x(2)
      gObj(7) = - x(3) + x(1)
      gObj(8) = - x(5) + x(3)
      gObj(9) =   x(4)

      end ! subroutine fnobj1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine fncon1
     &   ( mode, ncnln, n, ldg, needc, x, fCon, gCon, nState )

      implicit
     &     none
      integer
     &     mode, ncnln, n, ldg, nState, needc(*)
      double precision
     &     x(n), fCon(*), gCon(ldg,*)

*     ==================================================================
*     fncon1  computes the values and first derivatives of the nonlinear
*     constraints.
*
*     The zero elements of Jacobian matrix are set only once.  This
*     occurs during the first call to fncon1  (nState = 1).
*     ==================================================================
      integer            Out
      parameter         (Out  = 6)

      integer
     &     i, j
*     ------------------------------------------------------------------
      double precision   zero,          two
      parameter         (zero = 0.0d+0, two = 2.0d+0)
*     ------------------------------------------------------------------
      if (nState .eq. 1) then

*        First call to fncon1.  set all jacobian elements to zero.
*        Note: this will only work with `derivative level = 3'.

         do j = 1, n
            do i = 1, ncnln
               gCon(i,j) = zero
            end do
         end do

         if (Out .gt. 0) write(Out, '(/a)') ' Starting  Hex1 (Maximize)'

      else  if (nState .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing Hex1'
         return

      end if

      if (needc(1) .gt. 0) then
         fCon(1)    =   x(1)**2  +  x(6)**2
         gCon(1,1)  =   two*x(1)
         gCon(1,6)  =   two*x(6)
      end if

      if (needc(2) .gt. 0) then
         fCon(2)    =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         gCon(2,1)  = - two*(x(2) - x(1))
         gCon(2,2)  =   two*(x(2) - x(1))
         gCon(2,6)  = - two*(x(7) - x(6))
         gCon(2,7)  =   two*(x(7) - x(6))
      end if

      if (needc(3) .gt. 0) then
         fCon(3)    =   (x(3) - x(1))**2  +  x(6)**2
         gCon(3,1)  = - two*(x(3) - x(1))
         gCon(3,3)  =   two*(x(3) - x(1))
         gCon(3,6)  =   two*x(6)
      end if

      if (needc(4) .gt. 0) then
         fCon(4)    =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         gCon(4,1)  =   two*(x(1) - x(4))
         gCon(4,4)  = - two*(x(1) - x(4))
         gCon(4,6)  =   two*(x(6) - x(8))
         gCon(4,8)  = - two*(x(6) - x(8))
      end if

      if (needc(5) .gt. 0) then
         fCon(5)    =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         gCon(5,1)  =   two*(x(1) - x(5))
         gCon(5,5)  = - two*(x(1) - x(5))
         gCon(5,6)  =   two*(x(6) - x(9))
         gCon(5,9)  = - two*(x(6) - x(9))
      end if

      if (needc(6) .gt. 0) then
         fCon(6)    =   x(2)**2  +  x(7)**2
         gCon(6,2)  =   two*x(2)
         gCon(6,7)  =   two*x(7)
      end if

      if (needc(7) .gt. 0) then
         fCon(7)    =   (x(3) - x(2))**2  +  x(7)**2
         gCon(7,2)  = - two*(x(3) - x(2))
         gCon(7,3)  =   two*(x(3) - x(2))
         gCon(7,7)  =   two*x(7)
      end if

      if (needc(8) .gt. 0) then
         fCon(8)    =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         gCon(8,2)  = - two*(x(4) - x(2))
         gCon(8,4)  =   two*(x(4) - x(2))
         gCon(8,7)  = - two*(x(8) - x(7))
         gCon(8,8)  =   two*(x(8) - x(7))
      end if

      if (needc(9) .gt. 0) then
         fCon(9)    =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         gCon(9,2)  =   two*(x(2) - x(5))
         gCon(9,5)  = - two*(x(2) - x(5))
         gCon(9,7)  =   two*(x(7) - x(9))
         gCon(9,9)  = - two*(x(7) - x(9))
      end if

      if (needc(10) .gt. 0) then
         fCon(10)   =   (x(4) - x(3))**2  +  x(8)**2
         gCon(10,3) = - two*(x(4) - x(3))
         gCon(10,4) =   two*(x(4) - x(3))
         gCon(10,8) =   two*x(8)
      end if

      if (needc(11) .gt. 0) then
         fCon(11)   =   (x(5) - x(3))**2  +  x(9)**2
         gCon(11,3) = - two*(x(5) - x(3))
         gCon(11,5) =   two*(x(5) - x(3))
         gCon(11,9) =   two*x(9)
      end if

      if (needc(12) .gt. 0) then
         fCon(12)   =   x(4)**2  +  x(8)**2
         gCon(12,4) =   two*x(4)
         gCon(12,8) =   two*x(8)
      end if

      if (needc(13) .gt. 0) then
         fCon(13)   =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         gCon(13,4) =   two*(x(4) - x(5))
         gCon(13,5) = - two*(x(4) - x(5))
         gCon(13,8) = - two*(x(9) - x(8))
         gCon(13,9) =   two*(x(9) - x(8))
      end if

      if (needc(14) .gt. 0) then
         fCon(14)   =   x(5)**2  +  x(9)**2
         gCon(14,5) =   two*x(5)
         gCon(14,9) =   two*x(9)
      end if

      end ! subroutine fncon1

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine fnobj2
     &   ( mode, n, x, fObj, gObj, nState )

      implicit
     &     none
      integer
     &     mode, n, nstate
      double precision
     &     fObj, x(n), gObj(n)
*     ==================================================================
*     fnobj2  computes the value and some first derivatives of the
*     nonlinear objective function.
*     ==================================================================
      fObj    = - x(2)*x(6) + x(1)*x(7) - x(3)*x(7) - x(5)*x(8)
     &          + x(4)*x(9) + x(3)*x(8)

      gObj(3) = - x(7) + x(8)
      gObj(7) = - x(3) + x(1)
      gObj(8) = - x(5) + x(3)

      end ! subroutine fnobj2

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine fncon2
     &   ( mode, ncnln, n, ldg, needc, x, fCon, gCon, nState )

      implicit
     &     none
      integer
     &     mode, ncnln, n, ldg, nstate, needc(*)
      double precision
     &     x(n), fCon(*), gCon(ldg,*)
*     ==================================================================
*     fncon2  computes the values and the non-constant derivatives of
*     the nonlinear constraints.
*     ==================================================================
      integer            Out
      parameter         (Out  = 6)

      double precision   two
      parameter         (two = 2.0d+0)
*     ------------------------------------------------------------------
*     Print something on the first and last entry.

      if (      nState .eq. 1) then ! First
         if (Out .gt. 0) write(Out, '(/a)') ' Starting  Hex2 (Maximize)'
      else  if (nState .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finishing Hex2'
         return
      end if

      if (needc(1) .gt. 0) then
         fCon(1)    =   x(1)**2  +  x(6)**2
         gCon(1,1)  =   two*x(1)
         gCon(1,6)  =   two*x(6)
      end if

      if (needc(2) .gt. 0) then
         fCon(2)    =   (x(2) - x(1))**2  +  (x(7) - x(6))**2
         gCon(2,1)  = - two*(x(2) - x(1))
         gCon(2,2)  =   two*(x(2) - x(1))
         gCon(2,6)  = - two*(x(7) - x(6))
         gCon(2,7)  =   two*(x(7) - x(6))
      end if

      if (needc(3) .gt. 0) then
         fCon(3)    =   (x(3) - x(1))**2  +  x(6)**2
         gCon(3,1)  = - two*(x(3) - x(1))
         gCon(3,3)  =   two*(x(3) - x(1))
         gCon(3,6)  =   two*x(6)
      end if

      if (needc(4) .gt. 0) then
         fCon(4)    =   (x(1) - x(4))**2  +  (x(6) - x(8))**2
         gCon(4,1)  =   two*(x(1) - x(4))
         gCon(4,4)  = - two*(x(1) - x(4))
         gCon(4,6)  =   two*(x(6) - x(8))
         gCon(4,8)  = - two*(x(6) - x(8))
      end if

      if (needc(5) .gt. 0) then
         fCon(5)    =   (x(1) - x(5))**2  +  (x(6) - x(9))**2
         gCon(5,1)  =   two*(x(1) - x(5))
         gCon(5,5)  = - two*(x(1) - x(5))
         gCon(5,6)  =   two*(x(6) - x(9))
         gCon(5,9)  = - two*(x(6) - x(9))
      end if

      if (needc(6) .gt. 0) then
         fCon(6)    =   x(2)**2  +  x(7)**2
         gCon(6,2)  =   two*x(2)
         gCon(6,7)  =   two*x(7)
      end if

      if (needc(7) .gt. 0) then
         fCon(7)    =   (x(3) - x(2))**2  +  x(7)**2
         gCon(7,2)  = - two*(x(3) - x(2))
         gCon(7,3)  =   two*(x(3) - x(2))
         gCon(7,7)  =   two*x(7)
      end if

      if (needc(8) .gt. 0) then
         fCon(8)    =   (x(4) - x(2))**2  +  (x(8) - x(7))**2
         gCon(8,2)  = - two*(x(4) - x(2))
         gCon(8,4)  =   two*(x(4) - x(2))
         gCon(8,7)  = - two*(x(8) - x(7))
         gCon(8,8)  =   two*(x(8) - x(7))
      end if

      if (needc(9) .gt. 0) then
         fCon(9)    =   (x(2) - x(5))**2  +  (x(7) - x(9))**2
         gCon(9,2)  =   two*(x(2) - x(5))
         gCon(9,5)  = - two*(x(2) - x(5))
         gCon(9,7)  =   two*(x(7) - x(9))
         gCon(9,9)  = - two*(x(7) - x(9))
      end if

      if (needc(10) .gt. 0) then
         fCon(10)   =   (x(4) - x(3))**2  +  x(8)**2
         gCon(10,3) = - two*(x(4) - x(3))
         gCon(10,4) =   two*(x(4) - x(3))
         gCon(10,8) =   two*x(8)
      end if

      if (needc(11) .gt. 0) then
         fCon(11)   =   (x(5) - x(3))**2  +  x(9)**2
         gCon(11,3) = - two*(x(5) - x(3))
         gCon(11,5) =   two*(x(5) - x(3))
         gCon(11,9) =   two*x(9)
      end if

      if (needc(12) .gt. 0) then
         fCon(12)   =   x(4)**2  +  x(8)**2
         gCon(12,4) =   two*x(4)
         gCon(12,8) =   two*x(8)
      end if

      if (needc(13) .gt. 0) then
         fCon(13)   =   (x(4) - x(5))**2  +  (x(9) - x(8))**2
         gCon(13,4) =   two*(x(4) - x(5))
         gCon(13,5) = - two*(x(4) - x(5))
         gCon(13,8) = - two*(x(9) - x(8))
         gCon(13,9) =   two*(x(9) - x(8))
      end if

      if (needc(14) .gt. 0) then
         fCon(14)   =   x(5)**2  +  x(9)**2
         gCon(14,5) =   two*x(5)
         gCon(14,9) =   two*x(9)
      end if

      end ! subroutine fncon2
