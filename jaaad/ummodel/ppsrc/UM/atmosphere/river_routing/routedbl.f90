
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE ROUTEDBL(RUNOFF,RU,RATMED,NDEV,DT                      &
     &,     NX,NY,IGRCN,ISEQ,STO,JMAX,RMISS                             &
     &,     DOUT,DRUNIN)
!
!  Purpose: To route runoff
!  Method:
!  The runoff is routed along predetermined routes using an initial
!  water storage which is updated each timestep. The code and
! river direction and sequence files were produced by Taikan Oki
! (March 26, 1998 - see ) as
! standalone code and was adapted (as little as possible)
! by Cyndy Bunton to run in the UM GCM. The spinup was removed
! and initial water storage, river sequence and direction files
! input via the Atmosphere dump. Gridbox outflow, current water
! storage and total gridbox runoff are output via stash.
!        by Cyndy Bunton 28.02.03

!  Model            Modification history from model version 5.5:
! version  Date
!   6.0  12.08.03  Initialise new water storage in routine outflow1.
!                  Improve code efficiency for 'next point downstream'.
!                             C.Bunton
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!     In arrays dat(i, j),
!                         londitude = real(i)-0.5    [degree east]
!                         latitude = 90.5 - real(j)  [degree north]c
!
!     Initial condition = `rinit' [mm] equivalent water in river
!                         channel
!
!     calculate for Jan. 87 - Dec.88 and out put.
!
!     from 24.September.1996, by Taikan OKI
!     from 09.June.1997, by Taikan OKI
!     real and improve time integration scheme
!            from 11.June.1997, by Taikan Oki
!
!  at IIS, Univ. of Tokyo using g77 compiler, from 13.Feb.1998 by T.
!  using cindir*8 for input directory name, from 16.Feb.1998 by T. Ok
!
!--------------------------------------------------------------------
!
! For Export Version of STRIP (Scheme for Total Runoff Integrating
! Pathw
! from 26.March.1998 by Taikan Oki
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!
! IN ARGUMENTS:
!
      INTEGER                                                           &
     & NDEV                                                             &
                       ! IN number of river routing timesteps/day
     &,JMAX                                                             &
                       ! IN
     &,NX                                                               &
                       ! IN number of points in a row
     &,NY                                                               &
                       ! IN number of rows
     &,IGRCN(NX,NY)                                                     &
                       ! IN gridbox outflow direction
     &,ISEQ(NX,NY)     ! IN river flow sequence

      REAL                                                              &
     & RUNOFF(NX,NY)                                                    &
                       ! IN total runoff (Kg/m2/day)
     &,DT                                                               &
                       ! IN river routing timestep (s)
     &,RU                                                               &
                       ! IN effective river flow velocity (m/s)
     &,RATMED                                                           &
                       ! IN River meander ratio
     &,RMISS           ! IN missing data indicator

! IN/OUT Arguments
      REAL                                                              &
     & STO(NX,NY)      ! River water storage (Kg)
! OUT arguments
      REAL                                                              &
     & DOUT(NX,NY)                                                      &
                       ! OUT gridbox outflow (Kg/s)
     &,DRUNIN(NX,NY)   ! OUT gridbox inflow (Kg/s)
!
! LOCAL:

      real rin(nx, ny), rout(nx, ny), rsto(nx, ny), rinp(nx, ny)

      real rlen(nx,ny)  ! distance between grids [m]
      REAL rvel(nx, ny) ! flow velocity          [m/s]
      REAL rc(nx, ny)   ! transfer coefficient   [1/s]
      REAL area(nx, ny) ! grid area              [m^2]
      REAL din(nx, ny)  ! total input to river storage [kg/s]
      REAL sto2(nx, ny) ! river channel storage [kg] at end of step
      REAL stoall, dinall, sto2all, doutall, dtotal(nx,ny)              &
     &     , drunall, drivall, dinput(nx, ny)  ! totals

      integer inextx(nx, ny) ! next point downstream (x direction)
      integer inexty(nx, ny) ! next point downstream (x direction)
      integer iy, im, idec, ih ! year,month,day,hour
      integer nseqmax        ! maximum points in any river sequence
      integer i, j           ! loop counters

! Variables required for inland basin correction
      REAL basin_flux,flux2ocean

!     Preparation
!

! DEPENDS ON: initial0
      call initial0(rin, nx, ny, jmax, rmiss                            &
     &     , igrcn, iseq, nseqmax, inextx, inexty, rlen, sto, area)

!                                    fix. the velocity here
!      write(6, *) 'routetDBL: put uniform velovity.'
! DEPENDS ON: setrval
      call setrval(nx, ny, rvel, ru, jmax)
!      write(6, *) 'routeDBL: set coefficient.'
! DEPENDS ON: setcoef
      call setcoef(nx, ny, rlen, rvel, ratmed, rc, jmax)

! DEPENDS ON: mmd2kgs
            call mmd2kgs(runoff, igrcn, drunin, area, nx, ny, rmiss,    &
     &                   jmax)

! DEPENDS ON: setrval
            call setrval(nx, ny, dtotal, 0.0, ny)
! DEPENDS ON: setrval
            call setrval(nx, ny, dinput, 0.0, ny)


!**********************************************************
! Inland basin fix
! drunin is in kg/s

      basin_flux=0.
      basin_flux=sum(drunin,mask=igrcn .eq. 10)
      where(igrcn .eq. 10) drunin=0.

!**********************************************************

! DEPENDS ON: outflow1
              call outflow1(sto, drunin, rc, dt                         &
     &             , igrcn, iseq, inextx, inexty, nseqmax, nx, ny,jmax  &
     &             , sto2, din, dout                                    &
     &             , drunall, drivall, stoall, sto2all, doutall,dinall)
!
!
!***********************************************************
!Inland basin fix
!Add flux from inland basin to river outflow points
      flux2ocean=sum(dout,MASK=igrcn .eq. 9)
      where (igrcn .eq. 9) dout=dout+(basin_flux/flux2ocean)*dout


! Add in the runoff at TRIP seapoints produced by interpolation from
! atmos runoff to conserve due to mismatch between grids
! Note that this will not appear in RINP, ROUT written by

      do i = 1, nx
        do j = 1, jmax
          if(iseq(i,j) == rmiss)then
            if(drunin(i,j) /= 0.0)then
              dout(i,j) = drunin(i,j)
              din(i,j) = drunin(i,j)
! Add in to the totals for water balance check
              drunall = drunall + drunin(i,j)
              dinall = dinall + din(i,j)*dt
              doutall = doutall + dout(i,j)*dt
            else
              dout(i,j) = rmiss
              din(i,j) = rmiss
            endif
          endif
        enddo
      enddo


! DEPENDS ON: wrtwblog
              call wrtWBlog(6, iy, im, idec, 0, ndev                    &
     &           , stoall, sto2all, dinall, doutall, drunall, drivall   &
     &           , dt)
!
! DEPENDS ON: cp2
              call cp2(sto2, sto, nx, ny, jmax)
!
!
      END SUBROUTINE ROUTEDBL









!
!  returns latituds at (iy) in (nla)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     iy from SOUTH to NORTH.
!
!     from 23.Feb.1996, by Taikan OKI
!

!
!  returns longitude at (ix) in (nlo)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     ix from west to east.
!
!     from 25.Feb.1996, by Taikan OKI
!

!


!


!


!


!

!
!  returns latituds at (iy) in (nla)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     iy from SOUTH to NORTH.
!
!     from 23.Feb.1996, by Taikan OKI
!

!

!

! ******************************COPYRIGHT******************************
