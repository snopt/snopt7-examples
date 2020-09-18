*     ------------------------------------------------------------------
*     File MyProb.f
*     This is a main program to illustrate the use of  SNOPT with the
*     free format wrapper snOptA,  part of the SNOPT 7 package.
*
*     29 Dec 2002: First version for SNOPT 6
*     10 Jul 2004: Updated for SNOPT 7
*     02 Sep 2004: Current version.
*     ------------------------------------------------------------------
      program
     &     MyProb
      implicit
     &     none
      integer
     &     maxF, maxn, nxname, nFname, lenA, lenG
      parameter
     &   ( maxF   = 30,
     &     maxn   = 10,
     &     lenA   = 50, lenG   = 100,
     &     nxname =  1, nFname =   1 )
      integer
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     xstate(maxn), Fstate(maxF)
      character
     &     lfile*20, Prob*8, xnames(nxname)*8, Fnames(nFname)*8
      double precision
     &     ObjAdd, sInf, A(lenA), Flow(maxF), Fupp(maxF), F(maxF),
     &     xlow(maxn), xupp(maxn), x(maxn), Fmul(maxF), xmul(maxn)
      integer
     &     Errors, neA, neG, ObjRow, INFO, iPrt, iPrint, iSpecs, iSum,
     &     iSumm, Major, mincw, miniw, minrw, nF, n, nInf, nOut, nS
      external
     &     usrfun
*     ------------------------------------------------------------------
*     SNOPT workspace

      integer               lenrw
      parameter          (  lenrw = 20000)
      double precision   rw(lenrw)
      integer               leniw
      parameter          (  leniw = 10000)
      integer            iw(leniw)
      integer               lencw
      parameter          (  lencw =   500)
      character          cw(lencw)*8

      integer             Cold,       Basis,      Warm
      parameter          (Cold   = 0, Basis  = 1, Warm  = 2)
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used by the calling program.

      iSpecs =  4
      iSumm  =  6
      iPrint =  9
      nOut   =  6

      lfile = 'MyProb.spc'
      open( iSpecs, file=lfile, status='OLD',     err=800 )
      lfile = 'MyProb.out'
      open( iPrint, file=lfile, status='UNKNOWN', err=800 )

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
         go to 990
      end if

*     Set up the data structure for the sparse Jacobian.
*     Assign dummy values for the nonlinear elements.

      Errors = 0

      call Toy
     &   ( Errors, maxF, maxn,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     Prob, nF, n,
     &     ObjAdd, ObjRow, xlow, xupp, Flow, Fupp,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (Errors .gt. 0) go to 910

*     ------------------------------------------------------------------
*     Specify any options not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      Major    = 250
      iPrt     =   0
      iSum     =   0
      call snseti
     &   ( 'Major Iteration limit', Major, iPrt, iSum, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

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
         go to 920
      end if

      write(nOut, *) ' '
      write(nOut, *) 'snOptA finished.'
      write(nOut, *) 'Input errors  =', Errors
      write(nOut, *) 'INFO          =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (ObjRow .gt. 0)
     &write(nOut, *) 'Obj           =', ObjAdd + F(ObjRow)
      if (INFO .ge. 30) go to 920

      stop

*     ------------------------------------------------------------------
*     Error exit.
*     ------------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'Insufficient space to hold the problem'
      stop

  920 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )
 4010 format(/  a, 2x, i6 )

      end ! main program

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine Toy
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
     &     Errors, maxF, maxn, neA, lenA, neG, lenG, nF, n,
     &     ObjRow, lencw, leniw, lenrw, xstate(maxn), Fstate(maxF),
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iw(leniw)
      double precision
     &     ObjAdd,
     &     A(lenA), xlow(maxn), xupp(maxn), Flow(maxF), Fupp(maxF),
     &     x(maxn), F(maxF), xmul(maxn), Fmul(maxF), rw(lenrw)
      character
     &     Prob*8, cw(lencw)*8

*     ==================================================================
*     Toy defines input data for the toy problem discussed in the
*     snoptA Users Guide.
*
*     Minimize      3*x(1) + 5*x(2) + (x(1) + x(3) + x(4))**2
*
*     subject to      x(1)         +   x(3)**2 +   x(4)**2     = 2
*                                    4*x(3)    + 2*x(4)       >= 0
*                             x(2) +   x(3)**4 +   x(4)**4     = 4
*                     x(1) >= 0,                         x(4) >= 0.
*
*     On exit:
*        nF  is the number of objective and constraint functions
*               (including linear and nonlinear)
*        n    is the number of variables.
*
*        (iGfun(k),jGvar(k)), k = 1,2,...,neG, define the coordinates
*             of the nonzero problem derivatives.
*             If (iGfun(k),jGvar(k)) = (i,j), G(k) is the ijth element
*             of the problem vector F(i), i = 0,1,2,...,nF,  with
*             objective function in position 0 and constraint functions
*             in positions  1  through  m.
*
*        (iAfun(k),jAvar(k),a(k)), k = 1,2,...,neA, are the coordinates
*             of the nonzero constant problem derivatives.
*
*     ==================================================================
      integer
     &     i, Obj
*     ------------------------------------------------------------------
      double precision     zero,         one ,         two
      parameter           (zero =0.0d+0, one  =1.0d+0, two    =2.0d+0 )
      double precision     four,         five,         plInfy
      parameter           (four =4.0d+0, five =5.0d+0, plInfy =1.0d+20)
*     ------------------------------------------------------------------
*     Give the problem a name.

      Prob   = 'Toy prob'

*     Assign the dimensions of the constraint Jacobian.

      nF     = 4
      Obj    = 1                ! Toy problem objective row
      ObjRow = 1                ! Could be 0
      n      = 4

*     Check that there is enough storage.

      Errors = 0
      if (nF     .gt. maxF ) Errors = 1
      if (n      .gt. maxn ) Errors = 1
      if (Errors .gt.     0) return

      neG        =  0

      neG        =  neG + 1
*     G(neG)     =  two*sum + three
      iGfun(neG) =  Obj
      jGvar(neG) =  1

      neG        =  neG + 1
*     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  3

      neG        =  neG + 1
*     G(neG)     =  two*sum
      iGfun(neG) =  Obj
      jGvar(neG) =  4

*     Nonlinear constraints (derivatives by row)

      neG        =  neG + 1
*     G(neG)     =  two*x(3)
      iGfun(neG) =  2
      jGvar(neG) =  3

      neG        =  neG + 1
*     G(neG)     =  two*x(4)
      iGfun(neG) =  2
      jGvar(neG) =  4

      neG        =  neG + 1
*     G(neG)     =  four*x(3)**3
      iGfun(neG) =  4
      jGvar(neG) =  3

      neG        =  neG + 1
*     G(neG)     =  four*x(4)**3
      iGfun(neG) =  4
      jGvar(neG) =  4

*     neG        = 7 derivatives in all

*     -------------------------------------------------------
*     Next we assign the list of constant derivative entries.
*     -------------------------------------------------------
      neA        =  0

      neA        =  neA + 1
      iAfun(neA) =  Obj
      jAvar(neA) =  2
      A(neA)     =  five

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  4
      A(neA)     =  four

      neA        =  neA + 1
      iAfun(neA) =  3
      jAvar(neA) =  3
      A(neA)     =  two

      neA        =  neA + 1
      iAfun(neA) =  2
      jAvar(neA) =  1
      A(neA)     =  one

      neA        =  neA + 1
      iAfun(neA) =  4
      jAvar(neA) =  2
      A(neA)     =  one

*     neA        =  5  derivatives in all

*     ----------------
*     Initial x.
*     ----------------
      ObjAdd = zero

      x(1)   =  one
      x(2)   =  one
      x(3)   =  one
      x(4)   =  one

      do i = 1, n
         xlow(i)   = -plInfy
         xupp(i)   =  plInfy
         xstate(i) =  0
      end do

      xlow(1) = zero
      xlow(2) = zero

*     The objective row is a free row.

      Flow(Obj) = -plInfy
      Fupp(Obj) =  plInfy

      Flow(3)   = zero
      Fupp(3)   = plInfy

      Flow(2)   = two           ! Equality constraint
      Fupp(2)   = two

      Flow(4)   = four          ! Equality constraint
      Fupp(4)   = four

      do i = 1, nF
         Fmul(i) = zero
      end do

      end ! subroutine Toy

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine usrfun
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

*     ==================================================================
*     Computes the nonlinear objective and constraint terms for the toy
*     problem featured in the SNOPT users guide.
*     nF = 6, n = 5.
*
*     Minimize      3*x(1) + (x(1) + x(2) + x(3))**2 + 5*x(4)
*
*     subject to             4*x(2)    + 2*x(3)               >= 0
*                     x(1) +   x(2)**2 +   x(3)**2             = 2
*                              x(2)**4 +   x(3)**4   +   x(4)  = 4
*
*                     x(1) >= 0,                         x(4) >= 0.
*
*     The triples (g(k),iGfun(k),jGvar(k)), k = 1,2,...,neG, define
*     the sparsity pattern and values of the nonlinear elements
*     of the Jacobian.
*     ==================================================================
      integer
     &     neG, Out, Obj
      double precision
     &     sum
*     ------------------------------------------------------------------
      double precision     two,         three
      parameter           (two = 2.0d+0,three = 3.0d+0)
      double precision     four
      parameter           (four = 4.0d+0)
*     ------------------------------------------------------------------
      Out = 15

      Obj = 1                  ! Objective row of F

*     --------------------------------------------
*     Print something on the first and last entry.
*     -------------------------------------------
      if (Status .eq. 1) then       ! First
         if (Out .gt. 0) write(Out, '(/a)') ' This is problem  Toy'
      else  if (Status .ge. 2) then ! Last
         if (Out .gt. 0) write(Out, '(/a)') ' Finished problem Toy'
         return
      end if

      sum = x(1) + x(3) + x(4)

      if (needF .gt. 0) then
         F(Obj) =   three*x(1) + sum**2
         F(2) =         x(3)**2 +   x(4)**2
*        F(3) =       4*x(3)    + 2*x(4)     ! Linear constraint omitted
         F(4) =         x(3)**4 +   x(4)**4
      end if

      neG = 0

      if (needG .gt. 0) then

         neG        =  neG + 1
         G(neG)     =  two*sum + three
*        iGfun(neG) =  Obj
*        jGvar(neG) =  1

         neG        =  neG + 1
         G(neG)     =  two*sum
*        iGfun(neG) =  Obj
*        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  two*sum
*        iGfun(neG) =  Obj
*        jGvar(neG) =  4

*        Nonlinear constraints (derivatives by row)

         neG        =  neG + 1
         G(neG)     =  two*x(3)
*        iGfun(neG) =  2
*        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  two*x(4)
*        iGfun(neG) =  2
*        jGvar(neG) =  4

         neG        =  neG + 1
         G(neG)     =  four*x(3)**3
*        iGfun(neG) =  4
*        jGvar(neG) =  3

         neG        =  neG + 1
         G(neG)     =  four*x(4)**3
*        iGfun(neG) =  4
*        jGvar(neG) =  4

*        neG        = 7 derivatives in all
      end if
      end ! subroutine usrfun

