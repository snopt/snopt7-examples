*     ------------------------------------------------------------------
*     File t6wood.f  (Unix version)
*     Illustrates using SNOPT on a nonlinearly constrained problem
*     with MPS data read from a file.
*
*     16 May 1998: First   version.
*     28 Oct 2003: Current version.
*     ------------------------------------------------------------------
      program
     &     t6main

      implicit
     &     none
      integer
     &     maxm, maxn, maxne
      parameter
     &   ( maxm   = 1000,
     &     maxn   = 1000,
     &     maxne  = 3000 )
      character
     &     PrbNms(5)*8, Names(maxm+maxn)*8
      integer
     &     indA(maxne) , hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)

*     SNOPT workspace

      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lenrw = 20000)
      parameter
     &   ( leniw = 10000)
      parameter
     &   ( lencw =   500)
      double precision
     &     rw(lenrw)
      integer
     &     iw(leniw)
      character
     &     cw(lencw)*8
*     ------------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      integer
     &     Errors, iSpecs, iPrint, iSumm, i1, i2, INFO, iMPS, iObj,
     &     itnlim, m, mincw, miniw, minrw, n, ne, nInf, nName, nnCon,
     &     nnJac, nOut, nnObj, nS
      double precision
     &     Obj, ObjAdd, sInf
      external
     &     t6con, dummy
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by t6wood.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't6wood.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't6wood.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )
      end if

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

*     ------------------------------------------------------------------
*     Set up the data structure for the constraints.
*     MPSinp needs to know the number of nonlinear variables, etc.
*     The following calls fetch values set in the SPECS file.
*     Optionally, these values can be set in-line.
*     ------------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Nonlinear constraints',         nnCon, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'Nonlinear Jacobian  variables', nnJac, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'Nonlinear Objective variables', nnObj, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'MPS file',                       iMPS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     The problem name is not needed---it is set by MPSinp.
*     Specify the OBJECTIVE, RHS, RANGES and BOUNDS to be selected
*     from the MPS file.  Blank names mean "select the first one".

*     PrbNms(1) = '        '    ! PROBLEM   name
      PrbNms(2) = '        '    ! OBJECTIVE name
      PrbNms(3) = '        '    ! RHS       name
      PrbNms(4) = '        '    ! RANGES    name
      PrbNms(5) = '        '    ! BOUNDS    name

      if ( byname ) then

*        Unix and DOS systems.  Open the MPS file.

         lfile = 't6wood.mps'
         open( iMPS, file=lfile, status='OLD', err=800 )
      end if

      call MPSinp
     &   ( iMPS, maxm, maxn, maxne,
     &     nnCon, nnJac, nnObj,
     &     m, n, ne,
     &     iObj, ObjAdd, PrbNms,
     &     Acol, indA, locA, bl, bu, Names,
     &     hs, x, pi,
     &     INFO, mincw, miniw, minrw, nS,
     &     cw, lencw, iw, leniw, rw, lenrw )
      close( iMPS )
      if (INFO .ne. 103) go to 990

      nName = m + n

*     ------------------------------------------------------------------
*     Specify any options not set in the Specs file.
*     i1 and i2 may refer to the Print and Summary file respectively.
*     Setting them to 0 suppresses printing.
*     ------------------------------------------------------------------
      itnlim = 10000
      i1     =   0
      i2     =   0
      call snSeti
     &   ( 'Iterations', itnlim, i1, i2, Errors,
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
     &     iObj, ObjAdd, PrbNms(1),
     &     t6con, dummy,
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
      write(nOut, *) 't6wood finished.'
      write(nOut, *) 'Input  errors =', Errors
      write(nOut, *) 'snOptB INFO   =', INFO
      write(nOut, *) 'nInf          =', nInf
      write(nOut, *) 'sInf          =', sInf
      if (iObj .gt. 0) then
         write(nOut, *)
     &               'Obj           =', ObjAdd + x(n+iObj) + Obj
      else
         write(nOut, *)
     &               'Obj           =', ObjAdd + Obj
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

      end ! program  t6main

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t6con
     &   ( mode, m, n, neJac, x, f, g,
     &     nState, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none

      integer
     &     mode, m, n, neJac, nState, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     x(n), f(m),  g(m,n), rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     t6con  is funcon for test problem t6wood,
*     a chemical engineering design problem.
*     Originally called  woplant  (wood plant?).
*     m = 5,  n = 10.
*
*     For test purposes, we test  Derivative level
*     to decide whether or not to compute gradients.
*
*     Dec 1981: Original SNOPT version obtained via Bruce Murtagh.
*     Oct 1991: Converted to f77 for test problem t6wood.
*     ==================================================================
      integer
     &     i, iPrint, iSumm, lvlDer
      double precision
     &     ak1, ak2, ak3, b1t, b2t, b3t, fd, fg, fp, fr, fr2, fra,
     &     frb, frc, fre, frp, r1, r2, r3, recip, rr, temp,
     &     one, two, three, four, half, tenk, vrho
*     ------------------------------------------------------------------
      parameter (one   = 1.0d+0,  two  = 2.0d+0,
     &           three = 3.0d+0,  four = 4.0d+0,
     &           half  = 0.5d+0,  tenk = 10000.0d+0,  vrho = 3000.0d+0)
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level

      if (nState .eq. 1) then
         if (iPrint .gt. 0) write(iPrint, 5) lvlDer
         if (iSumm  .gt. 0) write(iSumm , 5) lvlDer
    5    format(/' Starting  problem t6wood.  Derivative level =', i3 /)
      else if  (nState .ge. 2) then
         if (iPrint .gt. 0) write(iPrint, 6)
         if (iSumm  .gt. 0) write(iSumm , 6)
    6    format(/' Finishing problem t6wood.')
      end if

*     Transform to original variables.

      fg     = tenk*(one + x(1))
      fp     = tenk*(one + x(2))
      fd     = tenk*(one + x(3))
      fra    = tenk*(one + x(4))
      frp    = tenk*(one + x(5))
      fre    = tenk*(one + x(6))
      frb    = tenk*(one + x(7))
      frc    = tenk*(one + x(8))
      fr     = tenk*(one + x(9))
      temp   = 630.0d+0 + 50.0d+0*x(10)

*     Rate constants.

      ak1    = 5.9755d+09 * dexp(-1.2d+4/temp)
      ak2    = 2.5962d+12 * dexp(-1.5d+4/temp)
      ak3    = 9.6283d+15 * dexp(-2.0d+4/temp)

*     Rate terms.

      fr2    = fr**2
      r1     = ak1*fra*frb*vrho/fr2
      r2     = ak2*frb*frc*vrho/fr2
      r3     = ak3*frc*frp*vrho/fr2

*     Nonlinear functions.

      recip  = one/(fr - fg - fp)
      f(1)   = two*r2       - fd*recip*fre
      f(2)   = r2 - half*r3 - fd*recip*(frp - fp) - fp
      f(3)   = - r1         - fd*recip*fra
      f(4)   = - r1 - r2    - fd*recip*frb
      f(5)   = 1.5d+0*r3    - fg

*     Scale them.

      do i = 1, m
         f(i) = f(i) / tenk
      end do

*     Compute the Jacobian (if SNOPT wants it).

      if (mode .eq. 0  .or.  lvlDer .lt. 2) return

      b1t    = 1.2d+4/temp**2
      b2t    = 1.5d+4/temp**2
      b3t    = 2.0d+4/temp**2
      rr     = recip**2

      g(1,1) = - fd*fre*rr
      g(1,2) =   g(1,1)
      g(1,3) = - fre*recip
      g(1,6) = - fd *recip
      g(1,7) =   two*r2/frb
      g(1,8) =   two*r2/frc
      g(1,9) = - four*r2/fr - g(1,1)
      g(1,10)=   two*r2*b2t

      g(2,1) = - fd*(frp - fp)*rr
      g(2,2) =   fd*(fr - frp - fg)*rr - one
      g(2,3) = - (frp - fp)*recip
      g(2,5) = - half*r3/frp - fd*recip
      g(2,7) =   r2/frb
      g(2,8) =   (r2 - half*r3)/frc
      g(2,9) = - two*(r2 - half*r3)/fr - g(2,1)
      g(2,10)=   r2*b2t - half*r3*b3t

      g(3,1) = - fd*fra*rr
      g(3,2) =   g(3,1)
      g(3,3) = - fra*recip
      g(3,4) = - r1/fra - fd*recip
      g(3,7) = - r1/frb
      g(3,9) =   two*r1/fr - g(3,1)
      g(3,10)= - r1*b1t

      g(4,1) = - fd*frb*rr
      g(4,2) =   g(4,1)
      g(4,3) = - frb*recip
      g(4,4) = - r1/fra
      g(4,7) = - (r1 + r2)/frb - fd*recip
      g(4,8) = - r2/frc
      g(4,9) =   two*(r1+r2)/fr - g(4,1)
      g(4,10)= - r1*b1t - r2*b2t

      g(5,1) = - 1.0d+0
      g(5,5) =   1.5d+0*r3/frp
      g(5,8) =   1.5d+0*r3/frc
      g(5,9) = - three *r3/fr
      g(5,10)=   1.5d+0*r3*b3t

*     Rescale the temperature derivatives.

      do i = 1, m
         g(i,10) = g(i,10) * 5.0d-3
      end do

      end ! subroutine t6con

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dummy
     &   ( mode, n, x, f, g, nState,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     mode, n, nState, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     f, x(n), g(n), rw(lenrw)
      character
     &     cw(lencw)*8

*     ==================================================================
*     Problem t6wood.
*     No nonlinear objective.
*     ==================================================================

*     Relax

      end ! subroutine dummy

