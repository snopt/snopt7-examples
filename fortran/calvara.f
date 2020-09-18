*     ------------------------------------------------------------------
*     File calvarA.f
*     This is a main program to test subroutine snOptA, which is part
*     of the SNOPT package.
*     It generates two calvar problems and calls snOptA to solve them.
*
*     11 Jun 2006: First version of calvara.
*     ------------------------------------------------------------------
      program
     &     calvarmain
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG,
     &     lencw, leniw, lenrw
      parameter
     &   ( maxF   = 1,    maxn  = 5000,
     &     lenA   = 1,    lenG  = 5000,
     &     nxname = 1, nFname = 1 )
      integer
     &     iAfun(lenA), jAvar(lenA),
     &     iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character*8
     &     Prob, xnames(nxname), Fnames(nFname)
      double precision
     &     ObjAdd, sInf,
     &     A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)

*     ------------------------------------------------------------------
*     SNOPT workspace

      parameter
     &     (lenrw = 20000,
     &      leniw =  5000,
     &      lencw =   500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
*     ------------------------------------------------------------------
*     USER workspace (used by calvar routines)

      integer
     &     lenru, leniu, lencu
      parameter
     &     (lenru = 5000,
     &      leniu = 1,
     &      lencu = 1)
      integer
     &     iu(leniu)
      double precision
     &     ru(lenru)
      character
     &     cu(lencu)*8
*     ------------------------------------------------------------------
      integer
     &     DerOpt, Errors, neA, neG, ObjRow, INFO, iPrint, iPrt, iSpecs,
     &     iSumm, iSum, j, mincw, miniw, minrw, nF, n, nInf,
     &     nOut, nS
      logical
     &     byname
      character
     &     lfile*20, cw(lencw)*8
      external
     &     calvar1, calvar2, calvar3
*     ------------------------------------------------------------------
      integer            Cold,       Basis,      Warm
      parameter         (Cold   = 0, Basis  = 1, Warm  = 2)
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*     nOut    is an output file used here by sntest.

      iSpecs =  4  ! equivalenced to banana.spc
      iPrint =  9  ! equivalenced to banana.out
      iSumm  =  6  ! summary file goes to standard output...
      nOut   =  6  ! ... as do messages from this program.

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 'calvara.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 'calvara.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

*     ------------------------------------------------------------------
*     Set options to default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file.  This must include "Problem number nn"
*     for some integer nn.  This defines the number of variables.
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 910
      end if

*     ------------------------------------------------------------------
*     Set up the problem to be solved.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '----------------------------------------------'
      write(nOut, *) '1. First calvar problem with default QP solver'
      write(nOut, *) '----------------------------------------------'

      Errors = 0
      call caldat
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) then
         go to 910
      end if

*     ------------------------------------------------------------------
*     Specify options that were not set in the Specs file.
*     The parameters iPrt and iSum may refer to the Print and Summary
*     file respectively.  Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      DerOpt = 1
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option',  DerOpt, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Solve problem calvar1, using a Cold start (Start = 0).
*     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, calvar1,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'calvar1 (QP) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 40) go to 900

*     ------------------------------------------------------------------
*     Solve calvar1 again.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '--------------------------------------'
      write(nOut, *) '2. First calvar problem with QN option'
      write(nOut, *) '--------------------------------------'

      call snSeti
     &   ( 'Verify level',  0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'QPSolver QN', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      do j = 1, n
         xstate(j) = 0
*        x(j)      = zero
         x(j)      = float(j)/float(n)
      end do

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, calvar1,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'calvar1 (QN) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 40) go to 900

*     ------------------------------------------------------------------
*     Solve problem calvar2, using a Cold start (Start = 0).
*     ------------------------------------------------------------------
*     Set the initial value and status of each variable.
*     Set the range for the objective.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '-----------------------------------------------'
      write(nOut, *) '3. Second calvar problem with default QP solver'
      write(nOut, *) '-----------------------------------------------'

      call snSeti
     &   ( 'Verify level',  3, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'QPSolver Cholesky', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      do j = 1, n
         xstate(j) = 0
*        x(j)      = zero
         x(j)      = float(j)/float(n)
      end do

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, calvar2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'calvar2 (QP) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 40) go to 900

*     ------------------------------------------------------------------
*     Solve problem calvar3, using a Cold start (Start = 0).
*     ------------------------------------------------------------------
      write(nOut, *) ' '
      write(nOut, *) '---------------------------------------'
      write(nOut, *) '4. Second calvar problem with QN option'
      write(nOut, *) '---------------------------------------'

      call snSeti
     &   ( 'Verify level',  0, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'QPSolver QN', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      do j = 1, n
         xstate(j) = 0
*        x(j)      = zero
         x(j)      = float(j)/float(n)
      end do

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, calvar2,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 910
      end if

      write(nOut, *) ' '
      write(nOut, *) 'calvar2 (QN) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'snOptA INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 40) go to 900
      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
  910 stop

 4000 format(/  a, 2x, a  )

      end ! program calvarmain

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine caldat
     &   ( Errors, Prob, maxF, maxn, nF, n,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp, x, xstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

*     ==================================================================
*     caldat  generates data for the test problem calvar1
*
*     On entry,
*     maxF, maxn are upper limits on nF and n.
*
*     The  nF by n  Jacobian is written as the sum of
*     two  nF by n  sparse matrices G and A, i.e.,  J = A + G,  where
*     A  and  G  contain contain the constant and nonlinear
*     elements of J respectively.
*
*     The nonzero pattern of G and A is specified by listing the
*     coordinates (i.e., row and column indices) of each nonzero.
*     Note that the coordinates specify the overall STRUCTURE of the
*     sparsity, not the Jacobian entries that happen to be zero at
*     the initial point.
*
*     The coordinates of the kth nonzero of  G  are defined
*     by iGfun(k) and jGvar(k)  (i.e., if i=iGfun(k) and j=jGvar(k),
*     then G(k) is the ijth element of G.)  Any known values of G(k)
*     must be assigned by the user in the routine userfg.
*
*     The coordinates of the kth nonzero of  A  are defined by
*     iAfun(k) and jAvar(k)  (i.e., if i=iAfun(k) and j=jAvar(k),
*     then A(k) is the ijth element of A.)  All values of A must be
*     assigned before the call to SNOPT.
*
*     The elements of A and G can be stored in any order, (e.g., by rows
*     or columns or mixed).
*
*     RESTRICTIONS:
*     1.  A nonzero entry of J must be specified as either an element
*         of A, or an element of G, but NOT BOTH (i.e.,  coordinates of
*         A  and  G  must not overlap.  Elements that are a sum of a
*         constant and varying part must be included in G and loaded
*         by userfg.
*
*     2.  If the computed value of an element of G happens to be zero
*         at a given point, it must still be loaded in userfg. (The
*         order of the coordinates is meaningful in SNOPT.)
*
*     On exit,
*     Errors    is 0 if there is enough storage, 1 otherwise.
*     nF        is the number of problem functions
*               (objective and constraints, linear and nonlinear).
*     n         is the number of variables.
*     neG       is the number of nonzeros in Jn.
*     neA       is the number of nonzeros in Jc.
*     xlow      holds the lower bounds on x.
*     xupp      holds the upper bounds on x.
*     Flow      holds the lower bounds on F.
*     Fupp      holds the upper bounds on F.

*     xstate(1:n) is a set of initial states for each x  (0,1,2,3,4,5).
*     Fstate(1:m) is a set of initial states for each F  (0,1,2,3,4,5).
*     x (1:n)     is a set of initial values for x.
*     Fmul(1:m)   is a set of initial values for the dual variables.
*
*     11 Jun 2006: First version of calvara based on SNOPT 7.2 t2banana.
*     ==================================================================
      integer
     &     j, nOut, Obj
*     ------------------------------------------------------------------
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
*     Name the Problem.

      Prob = 'calvar'

      nOut = 6
*     ------------------------------------------------------------------
*     The following call fetches n, which defines the number of
*     variables.  It is specified at runtime in the SPECS file.
*     ------------------------------------------------------------------
      Errors = 0
      call sngeti
     &   ( 'Problem number', n, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     Check if there is enough storage.

      if ( n    .le. 1 .or.
     &     maxn .lt. n .or.
     &     lenG .lt. n     ) then
         write(nOut, *) 'Not enough storage to generate calvar',
     &                  ' problem with  n =', n
         Errors = 1
      end if

      if (Errors .ne. 0) go to 910

      write(nOut, *) 'Problem calvar.   n =', n

      nF     = 1
      Obj    = nF
      ObjRow = Obj

      ObjAdd = zero             ! no additive obj. constant

*     ------------------------------------------------------------------
*     Set the initial value and status of each variable.
*     Set the range for the objective.
*     For want of something better to do, make the variables x(1:n)
*     temporarily fixed at their current values.
*     The crash can set the rest.
*     ------------------------------------------------------------------
      do j = 1, n
         iGfun(j)  = Obj
         jGvar(j)  = j
         xlow(j)   = bminus
         xupp(j)   = bplus
         xstate(j) = 0
*        x(j)      = zero
         x(j)      = float(j)/float(n)
      end do

      neG    = n
      neA    = 0

      Flow(1) = bminus
      Fupp(1) = bplus

  910 return

      end ! subroutine caldat

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine calvar1
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     Calculus of variations function 1
!     n must be an even integer.

!     int_0^1  [y^2 + ydot atan(ydot) - log(sqrt(1+ydot^2))] dt
!
!     If n = 50,  fmin = 2.1387374355879d0
!     If n =100,  fmin = 2.138675071500110d0
!
!     Starting point:  x(i) = 0.0
!
!     ==================================================================
      intrinsic
     &     atan, log
      integer
     &     halfn, i, j, k, nn, nOut, ObjRow
      double precision
     &     absc(4), atan, yL, ypL, yR, ypR, fy, fyp,
     &     gaL, gbL, gaR, gbR, gaLR,
     &     gbLR, h, log, rn, aL, apL, bL, bpL, aR, apR,
     &     bR, bpR, Obj, qsum, sum, t, dL, dR, u, v, vL, vR,
     &     weight(4), wi, tL, tR, y, yp
!     ------------------------------------------------------------------
!     parameter        ( nmax = 200 )
!     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   two
      parameter         (two  = 2.0d+0)
!     ------------------------------------------------------------------
      nOut   = 0

      halfn  = n/2

!     First entry.  Print something on the standard output.
!     Set up the panels.

      if (Status .eq. 1) then
         if (nOut .gt. 0) then
            write(nOut, 1000)
         end if

         if (lenru .lt.  2*n+2) then
            write(6, 9999)
            Status = - 2 ! Terminates the run
            return
         end if

         rn    = halfn
         rn    = one/rn
         sum   = zero
         nn    = halfn + 1
         do i = 1, nn
            ru(i) = sum
            sum   = sum + rn
         end do
      end if

*     -------------
*     Normal entry.
*     -------------
*     Set the objective element.

      ObjRow = 1

!     Set the abscissas and weights for the quadrature over the panel.

      absc(1) = - 0.861136311594053d+0
      absc(2) = - 0.339981043584856d+0
      absc(3) = - absc(2)
      absc(4) = - absc(1)
      weight(1) = 0.347854845137454d+0
      weight(2) = 0.652145154862546d+0
      weight(3) = weight(2)
      weight(4) = weight(1)

      gaLR = zero
      gbLR = zero
      Obj  = zero

!     Integrate over the j-th panel [tL,tR].

      do j = 1, halfn
         tL  = ru(j)
         tR  = ru(j+1)
         h   = tR - tL
         u   = half*(tR + tL)
         v   = half*h

         if (j .eq. 1) then
            yL = one
         else
            yL = x(2*j-2)
         end if

         ypL = x(2*j-1)

         if (j .eq. halfn) then
            yR  = two
            ypR = x(n)
         else
            yR  = x(2*j)
            ypR = x(2*j+1)
         end if

         qsum = zero
         gaL  = zero
         gbL  = zero
         gaR  = zero
         gbR  = zero

         do i = 1, 4
            wi   = weight(i)
            t    = u + v*absc(i)
            dL   = t - tL
            vL   = dL/h
            dR   = t - tR
            vR   = dR/h
            aL   = vR*vR*(one + vL + vL)
            apL  = two*vR*(one + vL + vL + vR)/h
            bL   = vR*vR*dL
            bpL  = vR*(vL + vL + vR)
            aR   = vL*vL*(one - vR - vR)
            apR  = two*vL*(one - vR - vR - vL)/h
            bR   = vL*vL*dR
            bpR  = vL*(vR + vR + vL)
            y    = yL*aL  + ypL*bL  + yR*aR  + ypR*bR
            yp   = yL*apL + ypL*bpL + yR*apR + ypR*bpR

            fyp  = atan(yp)
            fy   = two*y
            qsum = qsum + (y*y + yp*fyp - half*log(one + yp*yp))*wi

            gaL  = gaL + (fy*aL + fyp*apL)*wi
            gbL  = gbL + (fy*bL + fyp*bpL)*wi
            gaR  = gaR + (fy*aR + fyp*apR)*wi
            gbR  = gbR + (fy*bR + fyp*bpR)*wi
         end do
         gaL = v*gaL
         gbL = v*gbL
         k   = 2*j - 2

         if (needG .gt. 0) then
            if (j .eq. 1) then
               G(1)   = gbL
            else
               G(k)   = gaL + gaLR
               G(k+1) = gbL + gbLR
            end if
         end if

         gaLR = v*gaR
         gbLR = v*gbR
         Obj  = Obj + v*qsum
      end do

      if (needG .gt. 0) then
         G(n) = gbLR
      end if

      if (needF .gt. 0) F(ObjRow) = Obj

*     ------------
*     Final entry.
*     ------------
      if (Status .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is problem  calvar1.'/)
 2000 format(/ ' Finished         calvar1.'/)
 9999 format(/ ' Not enough workspace to define calvar1.'/)

      end ! subroutine calvar1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine calvar2
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     calvar2 is the second calculus of variation problem
!
!     int_0^1  [100(y - ydot^2)^2 + (1 - ydot)^2] dt
!
!     if n = 50,  fmin = 5.2179888453973d1
!     if n =100,  fmin = 5.217988845250500d1
!
!     Starting point:  x(i) = zero
!
!     ==================================================================
      integer
     &     halfn, i, j, k, nn, nOut, ObjRow
      double precision
     &     absc(4), yL, ypL, yR, ypR, fy, fyp,
     &     gaL, gbL, gaR, gbR, gaLR, gbLR,
     &     h, rn, aL, apL, bL, bpL, aR, apR, bR,
     &     bpR, Obj, qsum, sum, t, p1, p2, dL, dR, u, v, vL, vR,
     &     weight(4), wi, tL, tR, y, yp
!     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one    = 1.0d+0)
      double precision   two,           ten,           twenty
      parameter         (two  = 2.0d+0, ten  =10.0d+0, twenty =20.0d+0)
!     ------------------------------------------------------------------
      nOut   = 0

      halfn  = n/2

*     First entry.  Print something on the standard output.

      if (Status .eq. 1) then
         if (nOut .gt. 0) then
            write(nOut, 1000)
         end if

         if (lenru .lt.  2*n+2) then
            write(6, 9999)
            Status = - 2 ! Terminates the run
            return
         end if
         rn    = halfn
         rn    = one/rn
         sum   = zero
         nn    = halfn + 1
         do i = 1, nn
            ru(i) = sum
            sum   = sum + rn
         end do
      end if

*     -------------
*     Normal entry.
*     -------------
      ObjRow = 1

!     Set the abscissas and weights for the quadrature over the panel

      absc(1) = - 0.861136311594053d+0
      absc(2) = - 0.339981043584856d+0
      absc(3) = - absc(2)
      absc(4) = - absc(1)
      weight(1) = 0.347854845137454d+0
      weight(2) = 0.652145154862546d+0
      weight(3) = weight(2)
      weight(4) = weight(1)

      gaLR = zero
      gbLR = zero
      Obj  = zero

!     Integrate over the j-th panel

      do j = 1, halfn
         tL   = ru(j)
         tR   = ru(j + 1)
         h    = tR - tL
         u    = half*(tR + tL)
         v    = half* h
         if (j .eq. 1) then
            yL = zero
         else
            yL = x(2*j-2)
         end if

         ypL = x(2*j - 1)

         if (j .eq. halfn) then
            yR  = one
            ypR = x(n)
         else
            yR  = x(2*j)
            ypR = x(2*j+1)
         end if

         qsum = zero
         gaL  = zero
         gbL  = zero
         gaR  = zero
         gbR  = zero

         do i = 1, 4
            wi   = weight(i)
            t    = u + v*absc(i)
            dL   = t - tL
            vL   = dL/h
            dR   = t - tR
            vR   = dR/h
            aL   = vR*vR*(one + vL + vL)
            apL  = two*vR*(one + vL + vL + vR)/h
            bL   = vR*vR*dL
            bpL  = vR*(vL + vL + vR)
            aR   = vL*vL*(one - vR - vR)
            apR  = two*vL*(one - vR - vR - vL)/h
            bR   = vL*vL*dR
            bpR  = vL*(vR + vR + vL)
            y    = yL*aL  + ypL*bL  + yR*aR  + ypR*bR
            yp   = yL*apL + ypL*bpL + yR*apR + ypR*bpR

            p1   =  ten*(yp - y*y)
            p2   =  one - y
            fyp  =  twenty*p1
            fy   = -two*(y*fyp + p2)
            qsum = qsum + (p1*p1 + p2*p2)*wi

            gaL  = gaL  + (fy*aL + fyp*apL)*wi
            gbL  = gbL  + (fy*bL + fyp*bpL)*wi
            gaR  = gaR  + (fy*aR + fyp*apR)*wi
            gbR  = gbR  + (fy*bR + fyp*bpR)*wi
         end do
         gaL = v*gaL
         gbL = v*gbL
         k   = 2*j - 2
         if (needG .gt. 0) then
            if (j .eq. 1) then
               G(1) = gbL
            else
               G(k)   = gaL + gaLR
               G(k+1) = gbL + gbLR
            end if
         end if

         gaLR = v*gaR
         gbLR = v*gbR
         Obj  = Obj + v*qsum
      end do

      if (needG .gt. 0) then
         G(n) = gbLR
      end if

      if (needF .gt. 0) F(ObjRow) = Obj

*     ------------
*     Final entry.
*     ------------
      if (Status .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is problem  calvar2.'/)
 2000 format(/ ' Finished         calvar2.'/)
 9999 format(/ ' Not enough workspace to define calvar2.'/)

      end ! subroutine calvar2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine calvar3
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, nF, n, lenG, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     calvar3 is the third calculus of variation problem
!
!     int_0^1  [exp(-2y^2)*(ydot^2 - 1)] dt
!
!     If n = 50,  fmin = -1.471978186169d-1
!     If n =100,  fmin = -1.471979018668d-1
!
!     starting point:  x(i) = zero
!
!     ==================================================================
      intrinsic
     &     exp
      integer
     &     halfn, i, j, k, nn, nOut, ObjRow
      double precision
     &     absc(4), yL, ypL, yR, ypR, ex, eprd, exp, fx, fy,
     &     gaL, gbL, gaR, gbR, gaLR, gbLR,
     &     h, rn, aL, apL, bL, bpL, aR, apR, bR,
     &     bpR, Obj, qsum, sum, t, dL, dR, u, v, vL, vR, weight(4), wi,
     &     tL, tR, y, yp
!     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   two,           four
      parameter         (two  = 2.0d+0, four = 4.0d+0)
!     ------------------------------------------------------------------
      nOut   = 0

      halfn  = n/2

*     First entry.  Print something on the standard output.

      if (Status .eq. 1) then
         if (nOut .gt. 0) then
            write(nOut, 1000)
         end if

         if (lenru .lt.  2*n+2) then
            write(6, 9999)
            Status = - 2 ! Terminates the run
            return
         end if

         rn    = halfn
         rn    = one/rn
         sum   = zero
         nn    = halfn + 1
         do i = 1, nn
            ru(i) = sum
            sum   = sum + rn
         end do
      end if

*     -------------
*     Normal entry.
*     -------------
*     Set the objective element.

      ObjRow = 1

!     Set the abscissas and weights for the quadrature over the panel.

      absc(1) = - 0.861136311594053d+0
      absc(2) = - 0.339981043584856d+0
      absc(3) = - absc(2)
      absc(4) = - absc(1)
      weight(1) = 0.347854845137454d+0
      weight(2) = 0.652145154862546d+0
      weight(3) = weight(2)
      weight(4) = weight(1)

      gaLR = zero
      gbLR = zero
      Obj  = zero

!     Integrate over the j-th panel

      do j = 1, halfn
         tL   = ru(j)
         tR   = ru(j + 1)
         h    = tR - tL
         u    = half*(tR + tL)
         v    = half*h

         if (j .eq. 1) then
            yL  = zero
         else
            yL  = x(2*j-2)
         end if

         ypL = x(2*j-1)
         if (j .eq. halfn) then
            yR  = one
            ypR = x(n)
         else
            yR  = x(2*j)
            ypR = x(2*j+1)
         end if

         qsum = zero
         gaL   = zero
         gbL   = zero
         gaR   = zero
         gbR   = zero

         do i = 1, 4
            wi   = weight(i)
            t    = u + v*absc(i)
            dL   = t - tL
            vL   = dL/h
            dR   = t - tR
            vR   = dR/h

            aL   = vR*vR*(one + vL + vL)
            apL  = two*vR*(one + vL + vL + vR)/h
            bL   = vR*vR*dL
            bpL  = vR*(vL + vL + vR)
            aR   = vL*vL*(one - vR - vR)
            apR  = two*vL*(one - vR - vR - vL)/h
            bR   = vL*vL*dR
            bpR  = vL*(vR + vR + vL)

            y    = yL*aL  + ypL*bL  + yR*aR  + ypR*bR
            yp   = yL*apL + ypL*bpL + yR*apR + ypR*bpR
            ex   = exp(-two*y*y)
            eprd = ex*(yp*yp - one)
            fy   = two*yp*ex
            fx   = -four*y*eprd
            qsum = qsum + eprd*wi

            gaL = gaL + (fx*aL + fy*apL)*wi
            gbL = gbL + (fx*bL + fy*bpL)*wi
            gaR = gaR + (fx*aR + fy*apR)*wi
            gbR = gbR + (fx*bR + fy*bpR)*wi
         end do
         gaL = v*gaL
         gbL = v*gbL

         if (needG .gt. 0) then
            k  = 2*j - 2
            if (j .eq. 1) then
               G(1) = gbL
            else
               G(k)   = gaL + gaLR
               G(k+1) = gbL + gbLR
            end if
         end if

         gaLR = v*gaR
         gbLR = v*gbR
         Obj  = Obj + v*qsum
      end do

      if (needG .gt. 0) then
         G(n) = gbLR
      end if

      if (needF .gt. 0) F(ObjRow) = Obj

*     ------------
*     Final entry.
*     ------------
      if (Status .ge. 2) then
         if (nOut .gt. 0) write(nOut, 2000)
      end if
      return

 1000 format(/ ' This is problem  calvar3.'/)
 2000 format(/ ' Finished         calvar3.'/)
 9999 format(/ ' Not enough workspace to define calvar3.'/)

      end ! subroutine calvar3

