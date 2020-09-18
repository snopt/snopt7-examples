*     ------------------------------------------------------------------
*     File catmixa.f
*
*     The AMPL model is:
*
*     param tf := 1;      # Final time
*     param nh;           # Number of subintervals
*     param x1_0;         # Initial condition for x1
*     param x2_0;         # Initial condition for x2
*
*     param alpha;        # smoothing parameter;
*
*     param h := tf/nh;
*
*     var  u {0..nh} <= 1, >= 0;
*     var x1 {0..nh};
*     var x2 {0..nh};
*
*     minimize objective:  -1 + x1[nh] + x2[nh]
*                             + alpha*h*sum{i in 0..nh-1} (u[i+1] - u[i])^2 ;
*
*     subject to ode1 {i in 0..(nh-1)}:
*     x1[i+1] = x1[i] + (h/2)*(u[i]*(10*x2[i]-x1[i]) + u[i+1]*(10*x2[i+1]-x1[i+1]));
*
*     subject to ode2 {i in 0..(nh-1)}:
*     x2[i+1] = x2[i] + (h/2)*(u[i]*(x1[i]-10*x2[i]) - (1-u[i])*x2[i] +
*                              u[i+1]*(x1[i+1]-10*x2[i+1]) - (1-u[i+1])*x2[i+1]);
*
*     subject to ic1:
*     x1[0] = x1_0;
*
*     subject to ic2:
*     x2[0] = x2_0;
*
*     Data:
*     param nh := 800;
*     param x1_0 := 1;
*     param x2_0:= 0;
*     param alpha := 0.0;;
*
*     let {i in 0..nh}  u[i] := 0;
*     let {i in 0..nh} x1[i] := 1;
*     let {i in 0..nh} x2[i] := 0;
*
*     29 Jul 2001: First version of catmixa.f, derived from spring.f.
*     16 Mar 2004: Current version.
*     ------------------------------------------------------------------
      program
     &     catmix

      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG, lencw, leniw, lenrw
      parameter
     &   ( maxF   = 2010,
     &     maxn   = 3010,
     &     lenA   = 10000, lenG  = 30000,
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

*     SNOPT workspace

      parameter
     &     (lenrw = 600000,
     &      leniw = 350000,
     &      lencw =     501)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character*8
     &     cw(lencw)
      integer
     &     Errors, neA, neG, ObjRow, INFO, iPrint, iPrt, iSpecs, iSumm,
     &     iSum, mincw, miniw, minrw, nF, n, nInf, nOut, nS
      character*20
     &     lfile
      external
     &     usrfun
*     ------------------------------------------------------------------
      integer            Cold
      parameter         (Cold   = 0)
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*     nOut    is an output file used here by sntest.

      iSpecs =  4  ! equivalenced to catmixa.spc
      iPrint =  9  ! equivalenced to catmixa.out
      iSumm  =  6  ! summary file goes to standard output...
      nOut   =  6  ! ... as do messages from this program.

*     Open the Specs and print files.

      lfile = 'catmixa.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )

      lfile = 'catmixa.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

*     ------------------------------------------------------------------
*     Set options to default values.
*     ------------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      Errors = 0

*     ------------------------------------------------------------------
*     Read a Specs file.  This must include "Problem number nh"
*     for some integer nh.  This defines 2*nh nonlinear constraints.
*     ------------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 900
      end if

*     ------------------------------------------------------------------
*     Generate a nh-period problem.
*     ------------------------------------------------------------------
      call catdat
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 900

*     ------------------------------------------------------------------
*     Go for it, using a Cold start.
*     ------------------------------------------------------------------
      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 900
      end if

      write(nOut, *) ' '
      write(nOut, *) 'catmixa finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .gt. 30) go to 900

*     ------------------------------------------------------------------
*     Repeat with the SNOPT CG option
*     ------------------------------------------------------------------
      call catdat
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 900

      iPrt   = 0
      iSum   = 0

      call snSet
     &   ( 'QPsolver CG         ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'CG preconditioning 1', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 900
      end if

      write(nOut, *) ' '
      write(nOut, *) 'catmixa (CG) finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .gt. 30) go to 900

*     ------------------------------------------------------------------
*     Repeat with the SNOPT quasi-Newton option
*     ------------------------------------------------------------------
      call catdat
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .ne. 0) go to 900

      iPrt   = 0
      iSum   = 0

      call snSet
     &   ( 'QPsolver QN         ', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snSet
     &   ( 'Solution         Yes', iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snOptA
     &   ( Cold, nF, n, nxname, nFname,
     &     ObjAdd, ObjRow, Prob, usrfun,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .eq. 82 .or. INFO .eq. 83 .or. INFO .eq. 84) then
         go to 900
      end if

      write(nOut, *) ' '
      write(nOut, *) 'catmixa finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .gt. 30) go to 900

      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  900 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'
      stop

 4000 format(/  a, 2x, a  )

      end ! program catmixa

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine catdat
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, n, nF,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn), Fstate(maxF),
     &     iAfun(lenA), jAvar(lenA),
     &     iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), xmul(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

*     ==================================================================
*     catdat  generates data for the test problem catmix.
*     The AMPL model is:
*
*     param tf := 1;      # Final time
*     param nh;           # Number of subintervals
*     param x1_0;         # Initial condition for x1
*     param x2_0;         # Initial condition for x2
*
*     param alpha;        # smoothing parameter;
*
*     param h := tf/nh;
*
*     var u {0..nh} <= 1, >= 0;
*     var x1 {0..nh};
*     var x2 {0..nh};
*
*     minimize objective:  -1 + x1[nh] + x2[nh]
*                             + alpha*h*sum{i in 0..nh-1} (u[i+1] - u[i])^2 ;
*
*     subject to ode1 {i in 0..(nh-1)}:
*     x1[i+1] = x1[i] + (h/2)*(u[i]*(10*x2[i]-x1[i]) + u[i+1]*(10*x2[i+1]-x1[i+1]));
*
*     subject to ode2 {i in 0..(nh-1)}:
*     x2[i+1] = x2[i] + (h/2)*(u[i]*(x1[i]-10*x2[i]) - (1-u[i])*x2[i] +
*                              u[i+1]*(x1[i+1]-10*x2[i+1]) - (1-u[i+1])*x2[i+1]);
*
*     subject to ic1:
*     x1[0] = x1_0;
*
*     subject to ic2:
*     x2[0] = x2_0;
*
*     Data:
*     param nh := 800;
*     param x1_0 := 1;
*     param x2_0:= 0;
*     param alpha := 0.0;;
*
*     let {i in 0..nh}  u[i] := 0;
*     let {i in 0..nh} x1[i] := 1;
*     let {i in 0..nh} x2[i] := 0;
*
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
*     must be assigned by the user in the routine usrfun.
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
*         by usrfun.
*
*     2.  If the computed value of an element of G happens to be zero
*         at a given point, it must still be loaded in usrfun. (The
*         order of the coordinates is meaningful in SNOPT.)
*
*     On exit,
*     Errors      is 0 if there is enough storage, 1 otherwise.
*     nF        is the number of problem functions
*               (objective and constraints, linear and nonlinear).
*     n         is the number of variables.
*     neG       is the number of nonzeros in Jn.
*     neA       is the number of nonzeros in Jc.
*     neH       is the number of nonzeros in H(neH).
*     xlow      holds the lower bounds on x.
*     xupp      holds the upper bounds on x.
*     Flow      holds the lower bounds on F.
*     Fupp      holds the upper bounds on F.

*     xstate(1:n)   are the initial states for each x  (0,1,2,3,4,5).
*     Fstate(1:nF) are the initial states for each F  (0,1,2,3,4,5).
*     x (1:n)       are the initial values for x.
*     Fmul(1:nF)   are the initial values for the dual variables.
*
*     29 Jul 2001: First version of catmixa.
*     17 Jan 2003: Current version.
*     ==================================================================
      integer
     &     i, nOut, ju, jx1, jx2, nCon, nh, ode1, ode2, Obj
*     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)
      double precision   bplus,            bminus
      parameter         (bplus  = 1.0d+20, bminus = -bplus)
*     ------------------------------------------------------------------
      nOut = 6
*     ------------------------------------------------------------------
*     The following call fetches nh, which defines the number of
*     nonlinear constraints.  It is specified at runtime in the
*     SPECS file.
*     ------------------------------------------------------------------
      Errors = 0
      call sngeti
     &   ( 'Problem number', nh, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     Check if there is enough storage.

      if ( nh   .le.  1          .or.
     &     maxF .lt.  2* nh + 1  .or.
     &     maxn .lt.  3*(nh + 1) .or.
     &     lenA .lt.  6* nh      .or.
     &     lenG .lt. 13* nh          ) then
         write(nOut, *) 'Not enough storage to generate a problem ',
     &                  'with  nh =', nh
         Errors = 1
      end if

      if (Errors .ge. 1) go to 910

*     Write nh into the problem name.

      write(prob, '(i8)') nh
      if      (nh .lt.  100) then
         prob(1:6) = 'Catmix'
      else if (nh .lt. 1000) then
         prob(1:5) = 'Catmi'
      else
         prob(1:3) = 'Cat'
      end if

      write(nOut, *) 'Problem CATMIX.   nh =', nh

      n     = 3*(nh + 1)
      nCon  = 2* nh
      nF    = nCon + 1

*     The AMPL format variables are ordered as follows:
*     variables x1:      1: nh+1
*     variables x2:   nh+2:2nh+2
*     variables u:   2nh+2:3nh+2
*     jx1, jx2, ju are the pointers for the variables in x.

      neA  = 0
      neG  = 0

      jx1  = 1                  ! points to start of x1 in  x.
      jx2  = jx1 + nh + 1       ! points to start of x2 in  x.
      ju   = jx2 + nh + 1       ! points to start of u  in  x.

      ode1   = 1                ! points to ode1 constraints in F.
      ode2   = ode1 + nh        ! points to ode2 constraints in F.
      Obj    = ode2 + nh        ! Objective row of F.
      ObjRow = Obj

*     Linear terms first.

      neA = neA + 1
      iAfun(neA) =  ObjRow
      jAvar(neA) =  jx1  + nh
      A(neA)     =  one

      neA = neA + 1
      iAfun(neA) =  ObjRow
      jAvar(neA) =  jx2  + nh
      A(neA)     =  one

*     ObjAdd is a constant to be added to the objective.

      ObjAdd = - one

      do i = 0, nh-1
         neG = neG + 1
         iGfun(neG) =  ObjRow
         jGvar(neG) =  ju  + i
*        Gobj1      =  two*alpha*h*(x(ju+i+1) - x(ju+i))
*        G(neG)     =  Gobj0 - Gobj1
*        Gobj0      =  Gobj1

*        First ode constraint.

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  jx1  + i
*        G(neG)     =  - one + half*h*x(ju+i)

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  jx1  + i + 1
*        G(neG)     =    one + half*h*x(ju+i+1)

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  jx2  + i
*        G(neG)     =        - five*h*x(ju+i)

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  jx2  + i + 1
*        G(neG)     =        - five*h*x(ju+i+1)

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  ju   + i
*        G(neG)     =        - half*h*(ten*x(jx2+i)   - x(jx1+i))

         neG = neG + 1
         iGfun(neG) =  ode1 + i
         jGvar(neG) =  ju   + i + 1
*        G(neG)     =        - half*h*(ten*x(jx2+i+1) - x(jx1+i+1))

*        second ode constraint.

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  jx1  + i
*        G(neG)     =        - half*h*x(ju+i)

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  jx1  + i + 1
*        G(neG)     =        - half*h*x(ju+i+1)

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  jx2  + i
*        G(neG)     = - one  + half*h*(nine*x(ju+i)   + one)

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  jx2  + i + 1
*        G(neG)     =   one  + half*h*(nine*x(ju+i+1) + one)

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  ju   + i
*        G(neG)     = - half*h*(x(jx1+i)   - nine*x(jx2+i))

         neG = neG + 1
         iGfun(neG) =  ode2 + i
         jGvar(neG) =  ju   + i + 1
*        G(neG)     = - half*h*(x(jx1+i+1) - nine*x(jx2+i+1))

      end do

      neG = neG + 1
      iGfun(neG) =  ObjRow
      jGvar(neG) =  ju  + nh


*     ------------------------------------------------------------------
*     Initialize the bounds
*     ------------------------------------------------------------------
      do i = 0, nh

*        Initialize the bounds and values of the  x1 variables.

         xlow(jx1+i)   =  bminus
         xupp(jx1+i)   =  bplus
         x(jx1+i)      =  one
         xstate(jx1+i) =  0

*        Initialize the bounds and values of the  x2 variables.

         xlow(jx2+i)   =  bminus
         xupp(jx2+i)   =  bplus
         x(jx2+i)      =  zero
         xstate(jx2+i) =  0

*        Initialize the bounds and values of the   u variables.

         xlow(ju+i)    =  zero
         xupp(ju+i)    =  one
         x(ju+i)       =  zero
         xstate(ju+i)  =  0

      end do

*     Fix the boundary conditions.

      xlow(jx1)   = one
      xupp(jx1)   = one
      x(jx1)      = one

      xlow(jx2)   = zero
      xupp(jx2)   = zero
      x(jx2)      = zero

*     Bounds on F (all equalities).

      do i = 1, nCon
         Flow(i) = zero
         Fupp(i) = zero
      end do

      do i = 1, nCon
         Fmul(i) = zero
      end do

*     Set the objective and its bounds.

      Fmul(ObjRow) =  zero
      Flow(ObjRow) =  bminus
      Fupp(ObjRow) =  bplus

  910 return

      end ! subroutine catdat

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun
     &   ( Status, n, x,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     Status, needF, needG, n, nF, lenG,
     &     lencu, leniu, lenru, iu(leniu)
      double precision
     &     F(nF), G(lenG), x(n), ru(lenru)
      character*8
     &     cu(lencu)

*     ==================================================================
*     This is usrfun for problem catmix
*
*     param tf := 1;      # Final time
*     param nh := 800;
*     param x1_0 := 1;
*     param x2_0:= 0;
*     param alpha := 0.0;;
*
*     let {i in 0..nh}  u[i] := 0;
*     let {i in 0..nh} x1[i] := 1;
*     let {i in 0..nh} x2[i] := 0;
*
*     ==================================================================
      integer
     &     i, jx1, jx2, ju, neG, nh, ode1, ode2, ObjRow
      double precision
     &     alpha, FObj, Gobj0, Gobj1, h, rnh
*     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   two
      parameter         (two  = 2.0d+0)
      double precision   five,          nine,          ten
      parameter         (five = 5.0d+0, nine = 9.0d+0, ten  =10.0d+0)
      double precision   tf
      parameter         (tf   = one)
*     ------------------------------------------------------------------
      alpha  = zero
      nh     = n/3 - 1
      rnh    = nh
      h      = tf/rnh

*     The original catmix variables are ordered as follows:
*     variables x1:       1: nh+1
*     variables x2:    nh+2:2nh+2
*     variables u:    2nh+2:3nh+2
*     jx1, jx2, ju are the base indices for the corresponding variables

      neG  = 0

*     param x1_0;         # Initial condition for x1
*     param x2_0;         # Initial condition for x2
*     param nh;           # Number of subintervals
*     param alpha;        # smoothing parameter;

*     var u {0..nh} <= 1, >= 0;
*     var x1 {0..nh};
*     var x2 {0..nh};

      jx1    = 1                ! points to start of x1 in  x.
      jx2    = jx1 + nh + 1     ! points to start of x2 in  x.
      ju     = jx2 + nh + 1     ! points to start of u  in  x.

      ode1   = 1                ! points to ode1 constraints in F.
      ode2   = ode1 + nh        ! points to ode2 constraints in F.
      ObjRow = ode2 + nh        ! Objective row of F.

*     Fobj   =  - one + x(jx1+nh) + x(jx2+nh)
      Fobj   =    zero
      Gobj0  =    zero

      do i = 0, nh-1

         if (needF .gt. 0) then
            Fobj = Fobj + alpha*h*(x(ju+i+1) - x(ju+i))**2

*           subject to ode1 {i in 0..(nh-1)}
            F(ode1+i) = x(jx1+i+1) - x(jx1+i)
     &           - half*h*(  x(ju+i)  *(ten*x(jx2+i)   - x(jx1+i))
     &                     + x(ju+i+1)*(ten*x(jx2+i+1) - x(jx1+i+1)))
*           subject to ode2 {i in 0..(nh-1)}
            F(ode2+i) = x(jx2+i+1) - x(jx2+i)
     &           - half*h*( x(ju+i)  *(x(jx1+i)   - ten*x(jx2+i))
     &                                    - (one-x(ju+i))  *x(jx2+i)
     &                     +x(ju+i+1)*(x(jx1+i+1) - ten*x(jx2+i+1))
     &                                    - (one-x(ju+i+1))*x(jx2+i+1))
         end if

         if (needG .gt. 0) then

*           Objective gradient.

            neG = neG + 1
*           iGfun(neG) =  ObjRow
*           jGvar(neG) =  ju  + i
            Gobj1      =  two*alpha*h*(x(ju+i+1) - x(ju+i))
            G(neG)     =  Gobj0 - Gobj1
            Gobj0      =  Gobj1

*           First ode constraint.

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  jx1  + i
            G(neG)     =  - one + half*h*x(ju+i)

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  jx1  + i + 1
            G(neG)     =    one + half*h*x(ju+i+1)

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  jx2  + i
            G(neG)     =        - five*h*x(ju+i)

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  jx2  + i + 1
            G(neG)     =        - five*h*x(ju+i+1)

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  ju   + i
            G(neG)     =        - half*h*(ten*x(jx2+i)   - x(jx1+i))

            neG = neG + 1
*           iGfun(neG) =  ode1 + i
*           jGvar(neG) =  ju   + i + 1
            G(neG)     =        - half*h*(ten*x(jx2+i+1) - x(jx1+i+1))

*           second ode constraint.

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  jx1  + i
            G(neG)     =        - half*h*x(ju+i)

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  jx1  + i + 1
            G(neG)     =        - half*h*x(ju+i+1)

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  jx2  + i
            G(neG)     = - one  + half*h*(nine*x(ju+i)   + one)

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  jx2  + i + 1
            G(neG)     =   one  + half*h*(nine*x(ju+i+1) + one)

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  ju   + i
            G(neG)     = - half*h*(x(jx1+i)   - nine*x(jx2+i))

            neG = neG + 1
*           iGfun(neG) =  ode2 + i
*           jGvar(neG) =  ju   + i + 1
            G(neG)     = - half*h*(x(jx1+i+1) - nine*x(jx2+i+1))
         end if
      end do

      if (needG .gt. 0) then
         neG = neG + 1
*        iGfun(neG) =  ObjRow
*        jGvar(neG) =  ju  + nh
         G(neG)     = -Gobj0
      end if

      if (needF .gt. 0) then
         F(ObjRow) = FObj
      end if

      end ! subroutine usrfun

