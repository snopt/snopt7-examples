*     ------------------------------------------------------------------
*     File t7etamacro.f  (Unix version)
*     Illustrates using SNOPT on a nonlinear problem with linear
*     constraints.
*
*     16 May 1998: First   version.
*     18 May 1998: Current version.
*     ------------------------------------------------------------------
      program            t7main
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

*     SQOPT workspace---------------------------------------------------
      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lencw = 500, leniw = 20000, lenrw = 50000 )
      character
     &     cw(lencw)*8
      integer
     &     iw(leniw)
      double precision
     &     rw(lenrw)
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
     &     obj, ObjAdd, sInf
      external
     &     dummy, t7obj
*     ------------------------------------------------------------------
*     Specify some of the SNOPT files.
*     iSpecs  is the Specs file   (0 if none).
*     iPrint  is the Print file   (0 if none).
*     iSumm   is the Summary file (0 if none).
*
*     nOut    is an output file used here by t7etamacro.

      iSpecs =  4
      iPrint =  9
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

*        Unix and DOS systems.  Open the Specs and print files.

         lfile = 't7etamacro.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't7etamacro.out'
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
*     Set up the data structure for the linear constraints.
*     MPSinp needs to know the number of nonlinear variables, etc.
*     The following calls fetch values set in the SPECS file.
*     Optionally, these values can be set in-line.
*     ------------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Nonlinear constraints',           nnCon, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &     ( 'Nonlinear Jacobian  variables', nnJac, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &     ( 'Nonlinear Objective variables', nnObj, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'MPS file',                         iMPS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     The problem name is not needed---it is set by MPSinp.
*     PrbNms(1) = '        '
*     Specify the OBJECTIVE, RHS, RANGES and BOUNDS to be selected
*     from the MPS file.  Blank names mean "select the first one".

      PrbNms(3) = '        '    ! RHS       name
      PrbNms(4) = '        '    ! RANGES    name
      PrbNms(5) = '        '    ! BOUNDS    name

*     Set the OBJECTIVE name to preempt the N row OPTIMALG.

      call snGetc
     &   ( 'Objective = ',                PrbNms(2), Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if ( byname ) then

*        Unix and DOS systems.  Open the MPS file.

         lfile = 't7etamacro.mps'
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
      call snOpt
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, PrbNms(1),
     &     dummy, t7obj,
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
      write(nOut, *) 't7etamacro finished.'
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
      if (INFO .ge. 30) go to 910

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

      end ! program t7main

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t7obj
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

*     ------------------------------------------------------------------
*     t7obj  is funobj for Alan Manne's energy model ETAMACRO.
*     ------------------------------------------------------------------
      integer
     &     it, iyr, j1, j2, j3, j4, j5, maxt, nOut, j
      double precision
     &     c, cdc10, cdc5, d, d5, dc, delta, dxe, dxk, dxn, e, gnp,
     &     one, pn, r, s5, s5j, spda, y, zel, zero, zkp, zne
*     ------------------------------------------------------------------
      intrinsic          log
      parameter        ( maxt = 16,  zero =  0.0d+0,  one = 1.0d+0 )

      double precision   a1, b1, alpha, beta, rho, rate, esub
      save               a1, b1, alpha, beta, rho, rate, esub

      double precision   deltat(maxt), zks(maxt), zls(maxt), zes(maxt),
     &                   zns   (maxt), zys(maxt), zlb(maxt)
      save               deltat      , zks      , zls      , zes      ,
     &                   zns         , zys      , zlb
      data               zlb   / 1.160d+0, 1.446d+0, 1.717d+0, 2.039d+0,
     &                           2.364d+0, 2.740d+0, 3.101d+0, 3.508d+0,
     &                           3.873d+0, 4.276d+0, 4.721d+0, 5.213d+0,
     &                           5.755d+0, 6.354d+0, 7.016d+0, 7.746d+0/

      data               rate  / 10.000d+0 /,  esub  / 0.2000d+0 /,
     &                   beta  / 0.4000d+0 /,  alpha / 0.3333d+0 /
*     ------------------------------------------------------------------
*     The PRINT file unit number is defined by the call to snInit
*     and saved in SNOPT workspace.

      nOut  = 6

      if (nState .eq. 1) then
*        ----------------------------------------------------
*        First entry.  Define some data.
*        ----------------------------------------------------
         spda   = 0.96d+0
         s5     = spda**5
         rho    = (esub - one) / esub
         delta  = one / (one + rate/100.0d+0)

*        Compute coefficients a1, b1 with given data at 1970.

         zel    = 1.650d+0
         zne    = 0.509d+0
         pn     = 0.080d+0
         gnp    = 1.360d+0
         zkp    = 2.5d+0 * gnp
         y      = gnp  +  zne * pn / (one - beta)
         b1     = pn * y**(rho - one) / (one - beta)
     &               * zel**(-rho*beta)
     &               * zne**(one - rho*(one - beta))
         a1     = b1 * zel**( rho*beta) * zne**(rho*(one - beta))
         a1     = (y**rho  -  a1) / zkp**(alpha*rho)
         if (nOut .gt. 0) then
            write(nOut, '(/a)') ' Starting  problem etamacro'
            write(nOut, '(/ 1p,   a   , e16.8, 4x,  a   , e16.8)')
     &                            ' a =', a1   ,    ' b =', b1
         end if

*        Compute deltat(j).

         d5        = delta**5
         deltat(1) = 1000
         do j = 2, maxt
            deltat(j) = deltat(j-1) * d5
         end do
         deltat(maxt) = deltat(maxt) / (one - d5)

*        Compute surviving quantities.

         do j = 1, maxt
            s5j    = s5**j
            zys(j) = y   * s5j
            zks(j) = zkp * s5j
            zls(j) = (zlb(j) - s5j)**(one - alpha)
            zes(j) = zel * s5j
            zns(j) = zne * s5j
         end do
      end if

*     -------------------------------------------
*     Normal entry.
*     Some output is produced on the final entry.
*     -------------------------------------------
      if (nState .ge. 2) then
         if (nOut .gt. 0) then
            write(nOut, '(/a)') ' Finishing problem etamacro'
            write(nOut, '(/ a / a)')
     &         ' Time series of Gross Output, Annual Consumption and',
     &         ' Cumulative Consumption discounted at 5% and 10% are'
         end if
      end if

      f      = zero
      cdc5   = zero
      cdc10  = zero

      do j     = 1, maxt
         j1    = j
         j2    = j1 + maxt
         j3    = j2 + maxt
         j4    = j3 + maxt
         j5    = j4 + maxt
         dxk   = x(j1) - zks(j)
         dxe   = x(j2) - zes(j)
         dxn   = x(j3) - zns(j)
         r     = (dxk**alpha * zls(j))**rho
         e     = ((dxe/dxn)**beta * dxn)**rho
         d     = (a1*r + b1*e)**(one/rho)
         y     = zys(j) + d
         c     = y - x(j4) - x(j5)
         f     = f + deltat(j) * log(c)

         if (nState .eq. 2) then
            it     = 5*(j - 1)
            cdc5   = cdc5  +  5*(one / 1.05d+0)**it * c
            cdc10  = cdc10 +  5*(one / 1.10d+0)**it * c
            iyr    = 1975  +  it
            if (nOut .gt. 0) then
               write(nOut, '(i5, 4f9.3)') iyr, y, c, cdc5, cdc10
            end if
         end if

*        Compute the gradient.

         dc     = deltat(j) / c
         d      = dc * d**(one - rho)
         g(j1)  = d  * a1 * r * alpha / dxk
         g(j2)  = d  * b1 * e * beta  / dxe
         g(j3)  = d  * b1 * e * (one - beta) / dxn
         g(j4)  = - dc
         g(j5)  = g(j4)
      end do

      end ! subroutine t7Obj

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dummy
     &   ( mode, nnCon, nnJac, neJac, x, fCon, gCon,
     &     nState, cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none

      integer
     &     mode, nnCon, nnJac, neJac, nState, lencu, leniu, lenru,
     &     iu(leniu)
      double precision
     &     x(nnJac), fCon(nnCon), gCon(neJac), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     Problem t7etamacro.
*     No nonlinear constraints.
*     ==================================================================

*     Relax

      end ! subroutine dummy

