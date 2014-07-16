#if defined(A26_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

       SUBROUTINE RIV_INTCTL(                                           &
     & XPA, XUA, XVA, YPA, YUA, YVA,                                    &
     & G_P_FIELD, G_R_FIELD, N_PROC, ME, RMDI,                          &
     & GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
     & INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
     & GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
     & RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
     & GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
     & FLANDG, RIV_STEP, RIV_VEL, RIV_MCOEF,                            &
     & TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
     &  DELTA_PHI,EARTH_RADIUS,FIRST,                                   &
     &  r_area, slope, flowobs1,r_inext,r_jnext,r_land,                 &
     &  substore,surfstore,flowin,bflowin,                              &
! IN/OUT accumulated runoff
     & TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
     & BOX_OUTFLOW, BOX_INFLOW, RIVEROUT_ATMOS                          &
! Add new arguments for inland basin outflow
! OUT INLAND BASINS
     &  ,INLANDOUT_ATMOS,INLANDOUT_RIV                                  &
     &  ,TOT_WBLAKE,TOT_SUBRUN,timestep,l_cable,wbratio)

USE mpl, ONLY: MPL_REAL


! Purpose:
! New Control routine for River routing for Global Model.
!
! Author:  C.B.Bunton 20/01/03
!
!  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.5:
! VERSION  DATE
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!  6.0    12/8/03   Initialise gather fields of gridbox outflow
!                   and inflow. C.Bunton
! 6.0  12.09.03  Extra arguments added for river routing 2A.V.A.Bell
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!-----------------------------------------------------------------

      IMPLICIT NONE

#include "parvars.h"
      INTEGER                                                           &
     & row_length                                                       &
                                 ! IN NO. OF COLUMNS IN ATMOSPHERE
     &,rows                                                             &
                                 ! IN NO. OF ROWS IN ATMOSPHERE
     &, global_row_length                                               &
                                 ! number of points on a row
     &, global_rows                                                     &
                                 ! NUMBER OF global rows
     &,land_points                                                      &
                                 ! IN number of landpoints
     &,RIVER_ROW_LENGTH                                                 &
                                 ! IN no. of columns in river grid
     &,RIVER_ROWS                                                       &
                                 ! IN no. of rows in river grid
     &,GLOBAL_RIVER_ROW_LENGTH                                          &
                                 ! IN global river row length
     &,GLOBAL_RIVER_ROWS                                                &
                                 ! IN GLOBAL river rows
     &,gather_pe_trip                                                   &
                                 ! IN pe River routing to be run on
     &, n_proc                                                          &
                                 ! IN Total number of processors
     &, me                       ! IN My processor number

      INTEGER                                                           &
     &  G_P_FIELD                                                       &
                                  ! IN size of global ATMOS field
     &, G_R_FIELD                                                       &
                                  ! IN Size of global river field
     &, land_index (land_points)  ! IN index of land to global points

      REAL                                                              &
     & TOT_WBLAKE(land_points)                                          &
     &,TOT_SUBRUN(land_points)                                         &
     &,WBLAKE_IN(row_length,rows)                                       &
     &,SUBROFF_IN(row_length,rows)                                      &
     &,GATHER_WBLAKE_IN(G_P_FIELD)                                      &
     &,GATHER_SUBROFF_IN(G_P_FIELD),                                    &
     & WBLAKE_P(G_P_FIELD),  & ! part of wblake subtracted from rnof
     & WBLAKE_T,             & ! total wb_lake
     & RNOFF2_T,             & ! total sub-runoff
     & wbratio                 ! WBLAKE_T/RNOFF2_T

      REAL                                                              &
     & TOT_SURF_RUNOFF(land_points)                                     &
                                   !IN Surf RUNOFF on land pts(KG/M2/S)
     &,TOT_SUB_RUNOFF(land_points)                                      &
                                   ! IN Subsurf.RUNOFF (KG/M2/S)
     &,RMDI                                                             &
                                  ! IN real missing data indicator
     &,XUA(0:row_length)                                                &
                                  ! IN Atmosphere UV longitude coords
     &,YUA(rows)                                                        &
                                  ! IN Atmosphere latitude coords
     &,XPA(row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     &,YPA(rows)                                                        &
                                  ! IN Atmosphere latitude coords
     &,XVA(row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     &,YVA(0:rows)                                                      &
                                  ! IN Atmosphere latitude coords
     &,A_BOXAREAS(row_length,rows)                                      &
                                  !IN ATMOS gridbox areas
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                  ! IN Land fraction on global field.
     &,DELTA_PHI                                                        &
                                  ! RCM gridsize (radians)
     &,EARTH_RADIUS               ! Mean earth radius (m)

      REAL                                                              &
     & TRIVDIR(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               ! IN river direction
     &,TRIVSEQ(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                               ! IN river sequence
     &,TWATSTOR(RIVER_ROW_LENGTH, RIVER_ROWS)                           &
                                               ! IN/OUT water store(Kg)
     &,RIV_VEL                                                          &
                                  ! IN river velocity
     &,RIV_MCOEF                                                        &
                                  ! IN meandering coefficient
     &,RIV_STEP                 & ! IN river timestep (secs)
     &,TIMESTEP                   ! IN timestep (secs)

      LOGICAL                                                           &
     & INVERT_ATMOS               ! IN True if ATMOS fields are S->N
!                                 ! for regridding runoff from ATMOS.

      REAL                                                              &
     & RIVEROUT_ATMOS(row_length,rows)                                  &
                                      ! OUT river flow out from each
!                           ! gridbox(KG/m2/S)
     &,BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                        &
                                                 ! OUT gridbox outflow
!                                ! river grid (Kg/s)
     &,BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)                         &
                                                  ! OUT gridbox runoff
!                                ! river grid(Kg/s)
! Declare new variables for inland basin outflow
     &,INLANDOUT_RIV(RIVER_ROW_LENGTH,RIVER_ROWS)                       &
! OUT TRIP OUTFLOW FROM INLAND BASINS ON TRIP GRID Kg/s
     &,INLANDOUT_ATMOS(ROW_LENGTH,ROWS)               !OUT
! TRIP OUTFLOW FROM  INLAND BASINS ON atmos GRID Kg/m2/s



! Ancillary variables for grid-to-grid model

       REAL                                                             &
     & R_AREA(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                           &
                       !ACCUMULATED AREAS FILE
     & R_INEXT(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                          &
                       ! X-COORDINATE OF DOWNSTREAM GRID PT
     & R_JNEXT(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                          &
                       ! Y-COORDINATE OF DOWNSTREAM GRID PT
     & SLOPE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                            &
                       ! SLOPES (NOT USED YET)
     & FLOWOBS1(GLOBAL_ROW_LENGTH,GLOBAL_ROWS),                         &
                       ! OPTIONAL INITIALISATION FOR FLOWS
     & R_LAND(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)
                       !LAND/RIVER DEPENDS ON VALUE OF A_THRESH

! PROGNOSTIC VARIABLES FOR GRID-TO-GRID MODEL

       REAL                                                             &
     &  SUBSTORE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                         &
                        ! ROUTING SUB_SURFACE STORE (MM)
     & ,SURFSTORE(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                        &
                        ! ROUTING SURFACE STORE (MM)
     & ,FLOWIN(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)                           &
                        !SURFACE LATERAL INFLOW (MM)
     & ,BFLOWIN(GLOBAL_ROW_LENGTH,GLOBAL_ROWS)
                        ! SUB-SURFACE LATERAL INFLOW (MM)

! Logical to detect first entry to routing code

       LOGICAL FIRST
       LOGICAL l_cable
! LOCAL

      INTEGER                                                           &
     & I,J,L,k
      REAL                                                              &
     & gather_TOT_RUNOFFIN(G_P_FIELD) ! TOTAL RATE OF RUNOFF (KG/M2/S)

      INTEGER                                                           &
     &  info                                                            &
     &, gc_info                                                         &
                                ! Return code from MPP
     &, icode                   ! Return code :0 Normal Exit :>0 Error

    real, dimension(1) :: sarray

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0
      Character(*) RoutineName
      Parameter ( RoutineName='RIV_INTCTL')

      LOGICAL                                                           &
     & INVERT_TRIP                                                      &
                               ! TRUE WHEN ROW INVERSION IS REQ
     &, REGRID                                                          &
                               ! TRUE if TRIP grid different to ATMOS
     &, CYCLIC_TRIP                                                     &
                               ! TRUE WHEN THE TRIP MODEL HAS CYCLIC
     &, GLOBAL_TRIP            ! TRUE WHEN TRIP GRID SURFACE IS SPHER
       PARAMETER(INVERT_TRIP=.FALSE.,CYCLIC_TRIP=.TRUE.,                &
     & GLOBAL_TRIP=.TRUE.,REGRID=.TRUE.)

!      REAL RMDI_TRIP
!      PARAMETER(RMDI_TRIP=-999)

      REAL                                                              &
     & SURF_RUNOFFIN(row_length,rows)                                   &
                                     !IN TOTAL RATE OF RUNOFF (KG/M2/S)
     &,SUB_RUNOFFIN(row_length,rows) !IN TOTAL RATE OF RUNOFF (KG/M2/S)

! Gathering and Scattering variables:

      REAL gather_riverout_ATMOS(G_P_FIELD)
                                     ! river outflow at seapoints on
!                                    ! the ATMOS grid
      REAL gather_trivdir(G_R_FIELD)
                                     ! global field of river direction
      REAL gather_trivseq (G_R_FIELD)
                                     ! global field of river sequence
      REAL gather_twatstor(G_R_FIELD)
                                     ! global field of water store (Kg)
      REAL gather_box_outflow (G_R_FIELD)
                                     ! global field of gridbox outflow
!                                    ! on river routing  grid  (Kg/s)
      REAL gather_box_inflow(G_R_FIELD)
                                     ! global field of gridbox runoff
!                                    ! on river routing  grid (Kg/s)

      REAL GATHER_SURF_RUNOFFIN(G_P_FIELD)
                                    ! field for gather runoffin to pe0
      REAL GATHER_SUB_RUNOFFIN(G_P_FIELD)
      REAL GATHER_SUB_RUNOFFIN1(G_P_FIELD)
                                    ! field for gather runoffin to pe0
      REAL gather_flandg(G_P_FIELD)
                                    ! field for gather to land/sea mask
      REAL gather_A_BOXAREAS(G_P_FIELD)
                                    ! field for gather to gridbox areas

! Declare new gather fields for inland basin outflow
       REAL gather_INLANDOUT_ATMOS(G_P_FIELD)
                                      ! TRIP outflow at INLAND points
!                                     on the ATMOS grid (Kg/m2/s)
       REAL gather_INLANDOUT_RIV(G_R_FIELD)
                                      !global field of gridbox runoff
!                                     on river routing  grid
!                                     INLAND BASINS ONLY (Kg/s)

       real  TOTWBLAKE,TOTORUNOFF,TOTNRUNOFF,TOTAREA

      EXTERNAL RIV_ROUT

! gather the TRIP variables to PE0 and call the TRIP river routing
! 1. Gather coupling fields from distributed processors onto a single
!    processor for input to river routing routines.

          info = 0

! Initialise gridbox outflow and inflow
        gather_box_outflow = 0.0
        gather_box_inflow = 0.0
! initialise new gather fields for inland basins

             gather_INLANDOUT_ATMOS=0.0
             gather_INLANDOUT_RIV=0.0

             GATHER_SUB_RUNOFFIN = 0.0
             GATHER_SUBROFF_IN = 0.0
             GATHER_WBLAKE_IN = 0.0
             gather_A_BOXAREAS = 0.0
             gather_flandg = 0.0


!********************************************************************
          DO J=1,rows
           DO I=1,row_length
             SURF_RUNOFFIN(i,j) = 0.0
           ENDDO
          ENDDO
          DO J=1,rows
           DO I=1,row_length
             SUB_RUNOFFIN(i,j) = 0.0
           ENDDO
          ENDDO

          if (l_cable) then
          DO J=1,rows
           DO I=1,row_length
             WBLAKE_IN(i,j)  = 0.0
             SUBROFF_IN(i,j) = 0.0
           ENDDO
          ENDDO
          endif

! Copy land points output back to full fields array.
          DO l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            SURF_RUNOFFIN(i,j) = TOT_SURF_RUNOFF(L)
            SUB_RUNOFFIN(i,j) = TOT_SUB_RUNOFF(L)
          if (l_cable) then
            WBLAKE_IN(i,j)  = TOT_WBLAKE(L)
            SUBROFF_IN(i,j) = TOT_SUBRUN(L)
          endif
          ENDDO

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SURF_RUNOFFIN,GATHER_SURF_RUNOFFIN,         &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SURF_RUNOFFIN'
            ICODE=1
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF


! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SUB_RUNOFFIN,GATHER_SUB_RUNOFFIN,           &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SUB_RUNOFFIN'
            ICODE=2
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF
          GATHER_SUB_RUNOFFIN1 = GATHER_SUB_RUNOFFIN

if (l_cable) then
! DEPENDS ON: gather_field
          CALL GATHER_FIELD(WBLAKE_IN,GATHER_WBLAKE_IN,                 &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of WBLAKE_IN'
            ICODE=8
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF
! DEPENDS ON: gather_field
          CALL GATHER_FIELD(SUBROFF_IN,GATHER_SUBROFF_IN,                 &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of SUBROFF_IN'
            ICODE=9
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF
endif

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(flandg,gather_flandg,                       &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of tile_frac_m1'
            ICODE=3
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF


! DEPENDS ON: gather_field
          CALL GATHER_FIELD(A_BOXAREAS,gather_A_BOXAREAS,               &
     &     lasize(1,fld_type_p,halo_type_no_halo),                      &
     &     lasize(2,fld_type_p,halo_type_no_halo),                      &
     &     glsize(1,fld_type_p),                                        &
     &     glsize(2,fld_type_p),                                        &
     &     fld_type_p,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of Boxareas'
            ICODE=4
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF


! DEPENDS ON: gather_field
          CALL GATHER_FIELD(TRIVDIR,gather_TRIVDIR,                     &
     &     lasize(1,fld_type_r,halo_type_no_halo),                      &
     &     lasize(2,fld_type_r,halo_type_no_halo),                      &
     &     glsize(1,fld_type_r),                                        &
     &     glsize(2,fld_type_r),                                        &
     &     fld_type_r,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of TRIVDIR'
            ICODE=5
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(TRIVSEQ,gather_TRIVSEQ,                     &
     &     lasize(1,fld_type_r,halo_type_no_halo),                      &
     &     lasize(2,fld_type_r,halo_type_no_halo),                      &
     &     glsize(1,fld_type_r),                                        &
     &     glsize(2,fld_type_r),                                        &
     &     fld_type_r,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of TRIVSEQ'
            ICODE=6
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! DEPENDS ON: gather_field
          CALL GATHER_FIELD(TWATSTOR,gather_TWATSTOR,                   &
     &     lasize(1,fld_type_r,halo_type_no_halo),                      &
     &     lasize(2,fld_type_r,halo_type_no_halo),                      &
     &     glsize(1,fld_type_r),                                        &
     &     glsize(2,fld_type_r),                                        &
     &     fld_type_r,halo_type_no_halo,                                &
     &     gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in gather of TWATSTOR'
            ICODE=7
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! Lestevens 31.10.13 - wblake fix
      WBLAKE_T = 0. ! real total
      RNOFF2_T = 0. ! real total
      WBLAKE_P = 0. ! real (g_p_field) rnof part
      wbratio  = 0. ! real ratio - Already init in atmos_phy2

! Set river routing to run on the last PE
          IF(ME == GATHER_PE_TRIP) THEN

!         k=7003
!         print 102,k,GATHER_SUB_RUNOFFIN(k),GATHER_SUBROFF_IN(k),GATHER_WBLAKE_IN(k),gather_A_BOXAREAS(k),gather_flandg(k)
!         k=7196
!         print 102,k,GATHER_SUB_RUNOFFIN(k),GATHER_SUBROFF_IN(k),GATHER_WBLAKE_IN(k),gather_A_BOXAREAS(k),gather_flandg(k),timestep
!      102   format('GSUBRUN0',i5,10e11.4)

    if (l_cable) then
      
      WBLAKE_T = sum( GATHER_WBLAKE_IN )
      RNOFF2_T = sum( GATHER_SUBROFF_IN )
      wbratio  = min(1.,WBLAKE_T/max(1.0,RNOFF2_T))

      !print *,'wblake 0', WBLAKE_T,RNOFF2_T,wbratio!,RIV_STEP,timestep

       WBLAKE_P          = wbratio*GATHER_SUBROFF_IN
       GATHER_SUBROFF_IN = (1.0 - wbratio)*GATHER_SUBROFF_IN
       RNOFF2_T= sum(GATHER_SUBROFF_IN) ! modified runoff
      !print *,'wblake 01', RNOFF2_T

       WBLAKE_T = WBLAKE_T - sum(WBLAKE_P)

       GATHER_SUB_RUNOFFIN =  0.0
        where (GATHER_FLANDG >= 1.e-3 )  
             GATHER_SUB_RUNOFFIN =  GATHER_SUBROFF_IN/    &
                 (GATHER_FLANDG*TIMESTEP*GATHER_A_BOXAREAS)
        endwhere 
    end if
       TOTAREA = sum( gather_A_BOXAREAS)
       TOTORUNOFF= sum(GATHER_SUB_RUNOFFIN1) ! modified runoff
       TOTNRUNOFF= sum(GATHER_SUB_RUNOFFIN) ! modified runoff
       TOTWBLAKE =  sum ( WBLAKE_IN )
       !print 121,TOTWBLAKE,TOTORUNOFF,TOTNRUNOFF,TOTORUNOFF,wbratio
!121    format(1x,'riverr,TOTWBLAKE', 10e11.6)
! Call the TRIP river routing scheme

!-------------------------------------------------------------------


! DEPENDS ON: riv_rout
      CALL RIV_ROUT(GATHER_SURF_RUNOFFIN,GATHER_SUB_RUNOFFIN,           &
     &      global_row_length, global_rows,gather_flandg,RMDI,          &
     &      INVERT_ATMOS,XPA,YPA,XUA,YUA,XVA,YVA,                       &
     &      gather_A_BOXAREAS,                                          &
     &      REGRID, CYCLIC_TRIP, GLOBAL_TRIP, INVERT_TRIP,              &
     &      global_river_row_length, global_river_rows,                 &
     &      RIV_VEL, RIV_MCOEF, RIV_STEP,                               &
     &      gather_TRIVDIR, gather_TRIVSEQ, gather_TWATSTOR,            &
     &      gather_BOX_OUTFLOW, gather_BOX_INFLOW,                      &
!  Add new gather fields for inland basin outflow
     &      gather_RIVEROUT_ATMOS,gather_INLANDOUT_RIV                  &
     &     ,gather_INLANDOUT_ATMOS)



! Set all total land pts to 0.0
          WHERE(gather_flandg == 1.0)gather_riverout_ATMOS=0.0


          ENDIF                         ! Single processor
!

if (l_cable) then
 !write(6,*) 'before bcast pe,wbratio=', me,wbratio
 sarray(1) = wbratio
 call gc_rbcast(1,1, GATHER_PE_TRIP,gc_all_proc_group,gc_info,sarray)
 wbratio = sarray(1)
 !write(6,*) 'after bcast pe,wbratio, ierr=', me,wbratio, gc_info
! CALL MPL_BCAST(wbratio,1,MPL_REAL,gather_pe_trip,gc_all_proc_group,gc_info)
!print 36, wbratio
!36 format('Lest ratio2',(1x,f14.9))
end if

! DEPENDS ON: scatter_field
         CALL SCATTER_FIELD(RIVEROUT_ATMOS,gather_riverout_ATMOS,       &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of RIVEROUT_ATMOS'
            ICODE=101
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! DEPENDS ON: scatter_field
          CALL SCATTER_FIELD(TWATSTOR,gather_TWATSTOR,                  &
     &    lasize(1,fld_type_r,halo_type_no_halo),                       &
     &    lasize(2,fld_type_r,halo_type_no_halo),                       &
     &    glsize(1,fld_type_r),                                         &
     &    glsize(2,fld_type_r),                                         &
     &    fld_type_r,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of TWATSTOR'
            ICODE=102
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF


! DEPENDS ON: scatter_field
          CALL SCATTER_FIELD(BOX_OUTFLOW,gather_BOX_OUTFLOW,            &
     &    lasize(1,fld_type_r,halo_type_no_halo),                       &
     &    lasize(2,fld_type_r,halo_type_no_halo),                       &
     &    glsize(1,fld_type_r),                                         &
     &    glsize(2,fld_type_r),                                         &
     &    fld_type_r,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of BOX_OUTFLOW'
            ICODE=103
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! DEPENDS ON: scatter_field
          CALL SCATTER_FIELD(BOX_INFLOW,gather_BOX_INFLOW,              &
     &    lasize(1,fld_type_r,halo_type_no_halo),                       &
     &    lasize(2,fld_type_r,halo_type_no_halo),                       &
     &    glsize(1,fld_type_r),                                         &
     &    glsize(2,fld_type_r),                                         &
     &    fld_type_r,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 : ERROR in scatter of BOX_INFLOW'
            ICODE=104
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF


! Scatter new inland basin fields

! DEPENDS ON: scatter_field
          CALL SCATTER_FIELD(INLANDOUT_ATMOS,                           &
     &    gather_INLANDOUT_ATMOS,                                       &
     &    lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    lasize(2,fld_type_p,halo_type_no_halo),                       &
     &    glsize(1,fld_type_p),                                         &
     &    glsize(2,fld_type_p),                                         &
     &    fld_type_p,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE=                                                   &
     &    'ATMPHB2 : ERROR in scatter of INLANDOUT_atmos'
            ICODE=105
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF

! DEPENDS ON: scatter_field
          CALL SCATTER_FIELD(INLANDOUT_RIV,                             &
     &    gather_INLANDOUT_RIV,                                         &
     &    lasize(1,fld_type_r,halo_type_no_halo),                       &
     &    lasize(2,fld_type_r,halo_type_no_halo),                       &
     &    glsize(1,fld_type_r),                                         &
     &    glsize(2,fld_type_r),                                         &
     &    fld_type_r,halo_type_no_halo,                                 &
     &    gather_pe_trip,GC_ALL_PROC_GROUP,info,cmessage)
          IF(info /= 0) THEN      ! Check return code
            CMESSAGE='ATMPHB2 :                                         &
     &      ERROR in scatter of INLANDOUT_TRIP'
            ICODE=106
! DEPENDS ON: ereport
            Call Ereport(RoutineName,icode,Cmessage)
          ENDIF



      RETURN
      END SUBROUTINE RIV_INTCTL
#endif
