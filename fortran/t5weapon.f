!     ------------------------------------------------------------------
!     File t5weapon.f
!     Illustrates using SNOPT on a linearly constrained problem.
!
!     15 May 1998: First   version.
!     27 Oct 2003: Current version.
!     24 Oct 2010: Comments updated
!     ------------------------------------------------------------------
      program
     &     t5main

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
     &     indA(maxne), hs(maxn+maxm)
      integer
     &     locA(maxn+1)
      double precision
     &     Acol(maxne), bl(maxn+maxm), bu(maxn+maxm),
     &     x(maxn+maxm), pi(maxm), rc(maxn+maxm)

      !SNOPT workspace--------------------------------------------------
      integer
     &     lencw, leniw, lenrw
      parameter
     &   ( lencw = 500, leniw = 10000, lenrw = 20000 )
      character
     &     cw(lencw)*8
      integer
     &     iw(leniw)
      double precision
     &     rw(lenrw)
      !-----------------------------------------------------------------
      logical
     &     byname
      character
     &     lfile*20
      integer
     &     Errors, iSpecs, iPrint, iSumm, i1, i2, INFO, iMPS, inewB,
     &     iObj, itnlim, m, mincw, miniw, minrw, n, ne, nInf, nName,
     &     nnCon, nnJac, nOut, nnObj, nS
      double precision
     &     Obj, ObjAdd, sInf
      external
     &     dummy, t5obj
      !-----------------------------------------------------------------
      ! Specify some of the SNOPT files.
      !  iSpecs  is the Specs file   (0 if none).
      !  iPrint  is the Print file   (0 if none).
      !  iSumm   is the Summary file (0 if none).
      !
      !  nOut    is an output file used here by t5weapon.

      iSpecs =  4
      iPrint =  9
      iNewB  = 11
      iSumm  =  6
      nOut   =  6

      byname = .true.

      if ( byname ) then

         ! Unix and DOS systems.  Open the Specs and print files.

         lfile = 't5weapon.spc'
         open( iSpecs, file=lfile, status='OLD',     err=800 )

         lfile = 't5weapon.out'
         open( iPrint, file=lfile, status='UNKNOWN', err=800 )

         lfile = 't5weapon.newbasis'
         open( iNewB, file=lfile, status='UNKNOWN', err=800 )
      end if

      !-----------------------------------------------------------------
      ! First,  snInit MUST be called to initialize optional parameters
      ! to their default values.
      !-----------------------------------------------------------------
      call snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Read a Specs file (Optional).
      !-----------------------------------------------------------------
      call snSpec
     &   ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ne. 101  .and.  INFO .ne. 107) then
         go to 990
      end if

      !-----------------------------------------------------------------
      ! Set up the data structure for the linear constraints.
      ! MPSinp needs to know the number of nonlinear variables, etc.
      ! The following calls fetch values set in the SPECS file.
      ! Optionally, these values can be set in-line.
      !-----------------------------------------------------------------
      Errors = 0

      call snGeti
     &   ( 'Nonlinear constraints        ', nnCon, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'Nonlinear Jacobian  variables', nnJac, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'Nonlinear Objective variables', nnObj, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snGeti
     &   ( 'MPS file                     ',  iMPS, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      ! The problem name is not needed---it is set by MPSinp.
      ! Specify the OBJECTIVE, RHS, RANGES and BOUNDS to be selected
      ! from the MPS file.  Blank names mean "select the first one".

*     PrbNms(1) = '        '    ! PROBLEM   name
      PrbNms(2) = '        '    ! OBJECTIVE name
      PrbNms(3) = '        '    ! RHS       name
      PrbNms(4) = '        '    ! RANGES    name
      PrbNms(5) = '        '    ! BOUNDS    name

      if ( byname ) then

         ! Unix and DOS systems.  Open the MPS file.

         lfile = 't5weapon.mps'
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

      !-----------------------------------------------------------------
      ! Specify any options not set in the Specs file.
      ! i1 and i2 may refer to the Print and Summary file respectively.
      ! Setting them to 0 suppresses printing.
      !-----------------------------------------------------------------
      itnlim = 1000
      i1     =    0
      i2     =    0
      call snSeti
     &   ( 'Iterations        ', itnlim, i1, i2, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Go for it, using a Cold start.
      ! hs    need not be set if a basis file is to be input.
      !       Otherwise, each hs(1:n) should be 0, 1, 2, 3, 4, or 5.
      !       The values are used by the Crash procedure s2crsh
      !       to choose an initial basis B.
      !       If hs(j) = 0 or 1, column j is eligible for B.
      !       If hs(j) = 2, column j is initially superbasic (not in B).
      !       If hs(j) = 3, column j is eligible for B and is given
      !                     preference over columns with hs(j) = 0 or 1.
      !       If hs(j) = 4 or 5, column j is initially nonbasic.
      !-----------------------------------------------------------------
      call snOptB
     &   ( 'Cold', m, n, ne, nName,
     &     nnCon, nnObj, nnJac,
     &     iObj, ObjAdd, PrbNms(1),
     &     dummy, t5obj,
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
      write(nOut, *) 't5weapon finished.'
      write(nOut, *) 'Input errors  =', Errors
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

      !-----------------------------------------------------------------
      ! Error exit.
      !-----------------------------------------------------------------
  800 write(nOut, 4000) 'Error while opening file', lfile
      stop

  910 write(nOut, *) ' '
      write(nOut, *) 'STOPPING because of error condition'

  990 stop

 4000 format(/  a, 2x, a  )

      end ! program t5main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine t5obj
     &   ( mode, n, x, f, g, nState,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     mode, n, nState, lencw, leniw, lenrw, iw(leniw)
      character
     &     cw(lencw)*8
      double precision
     &     f, x(n), g(n), rw(lenrw)

      !-----------------------------------------------------------------
      ! t5obj  is funobj for the Weapon Assignment problem t5weapon.
      ! It assumes the Specs file contains data after the End card.
      !-----------------------------------------------------------------
      intrinsic          log
      double precision   zero
      integer            nweapn,      ntargt
      parameter        ( nweapn = 5,  ntargt = 20,  zero   = 0.0d+0 )
      double precision   q(nweapn,ntargt), ql(nweapn,ntargt), u(ntargt)
      save               q               , ql               , u
      !-----------------------------------------------------------------
      integer
     &     i, iSpecs, j,k
      double precision
     &     xk, t
      !-----------------------------------------------------------------
      if (nState .eq. 1) then
         !----------------------------------------------------
         ! First entry.  Read some data defining the objective.
         !----------------------------------------------------
         ! Weapon data follows the SPECS data.
         ! The SPECS file unit number is defined by the call to snInit
         ! and saved in SNOPT workspace.

         iSpecs = iw( 11) ! Specs (options) file

         do i = 1, nweapn
            read (iSpecs, '(18f4.2)') (q(i,j), j = 1, ntargt)
            do j = 1, ntargt
               ql(i,j) = log( q(i,j) )
            end do
         end do
         read (iSpecs, '(18f4.0)') u
      end if

      !-------------
      ! Normal entry.
      !-------------
      k      = 0
      f      = zero

      do j = 1, ntargt
         t     = u(j)
         do i  = 1, nweapn
            xk = x(k+i)
            if (xk .gt. zero) t = t * q(i,j)**xk
         end do

         if (mode .eq. 2) then
            do i = 1, nweapn
               g(k+i) = t * ql(i,j)
            end do
         end if

         k     = k + nweapn
         f     = f + (t - u(j))
      end do

      end ! subroutine t5Obj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

      !=================================================================
      ! Problem t5weapon.
      !=================================================================

      ! Relax, no nonlinear constraints.

      end ! subroutine dummy

