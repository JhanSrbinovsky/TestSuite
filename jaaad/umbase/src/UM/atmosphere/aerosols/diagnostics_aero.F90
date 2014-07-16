#if defined(A17_2A) || defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      Subroutine diagnostics_aero(                                      &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      timestep                                   &
     &,                      at_extremity                               &
     &,                      L_SULPC_SO2, L_DMS_em                      &
     &,                      L_SULPC_DMS, L_SULPC_NEWDMS                &
     &,                      L_SULPC_OZONE, L_SULPC_NH3                 &
     &,                      L_SOOT                                     &
     &,                      MSA, NH3_DEP                               &
     &,                      DMS_emiss                                  &
     &,                      DELTAS_DMS                                 &
     &,                      F_DMS_TO_SO2                               &
     &,                      F_DMS_TO_SO4                               &
     &,                      F_DMS_TO_MSA                               &
     &,                      DELTAS_DRY                                 &
     &,                      DELTAS_WET                                 &
     &,                      DELTAS_WET_O3                              &
     &,                      DELTAS_EVAP                                &
     &,                      DELTAS_NUCL                                &
     &,                      DELTAS_DIFFUSE                             &
     &,                      DELTAS_COAG                                &
#if defined(A17_2B)
     &,                      DELTAS_MERGE                               &
#endif
     &,                      PSI                                        &
     &,                      PM10,      PM2p5                           &
     &,                      PM10_SO4,  PM2p5_SO4                       &
     &,                      PM10_BC,   PM2p5_BC                        &
     &,                      PM10_BB,   PM2p5_BB                        &
     &,                      PM10_OCFF, PM2p5_OCFF                      &
     &,                      PM10_SOA,  PM2p5_SOA                       &
     &,                      PM10_SS,   PM2p5_SS                        &
     &,                      PM10_DUST, PM2p5_DUST                      &
     &     ,                                                            &
#include "argsts.h"
     & STASHwork                                                        &
     &  )
!
!----------------------------------------------------------------------
! Purpose:  Calculates diagnostics from section 17 AERO_CTL2 routine
!           and outputs them.
!
! Current owner of code:      C E Johnson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20
!
!-----------------------------------------------------------------------
!
      Implicit None
!
! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
!
      Integer                                                           &
     &  row_length                                                      &
                            ! number of points on a row
     &, rows                                                            &
                            ! number of rows in a theta field
     &, n_rows                                                          &
                            ! number of rows in a v field
     &, model_levels                                                    &
                            ! number of model levels
     &, wet_model_levels    ! number of model levels where moisture
!
      Integer                                                           &
     &  global_row_length                                               &
                            ! NUMBER OF points on a global row
     &, global_rows                                                     &
                            ! NUMBER OF global rows
     &, me                                                              &
                            ! Processor number
     &, halo_i                                                          &
                            ! size of large halo in x direction
     &, halo_j                                                          &
                            ! size of large halo in y direction
     &, off_x                                                           &
                            ! size of small halo in x direction
     &, off_y                                                           &
                            ! size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

      Real                                                              &
     &  timestep
!
      Logical                                                           &
     &  L_SULPC_SO2                                                     &
                          ! T if S Cycle on
     &, L_SULPC_DMS                                                     &
                          ! T if DMS included
     &, L_DMS_em                                                        &
                          ! T if DMS emissions used
     &, L_SULPC_NEWDMS                                                  &
                          ! T if new DMS scheme used (requires OZONE)
     &, L_SULPC_OZONE                                                   &
                          ! T if OZONE field present
     &, L_SULPC_NH3                                                     &
                          ! T if NH3 field present
     &, L_SOOT            ! T if SOOT modelling on
!
! Arguments with intent IN/OUT (diagnostics):
      REAL                                                              &
     & MSA(row_length,rows,model_levels)                                &
                                                  !mmr S in MSA
     &,NH3_DEP(row_length,rows,model_levels)                            &
                                                  !NH3 depleted
     &,DMS_emiss(row_length,rows)                 !DMS emiss (kgSm-2s-1)
      REAL                                                              &
     & DELTAS_DMS(row_length,rows,model_levels)                         &
     &,F_DMS_TO_SO2(row_length,rows,model_levels)                       &
     &,F_DMS_TO_SO4(row_length,rows,model_levels)                       &
     &,F_DMS_TO_MSA(row_length,rows,model_levels)                       &
     &,DELTAS_DRY(row_length,rows,model_levels)                         &
     &,DELTAS_WET(row_length,rows,model_levels)                         &
     &,DELTAS_WET_O3(row_length,rows,model_levels)                      &
     &,DELTAS_EVAP(row_length,rows,model_levels)                        &
     &,DELTAS_NUCL(row_length,rows,model_levels)                        &
     &,DELTAS_DIFFUSE(row_length,rows,model_levels)                     &
     &,DELTAS_COAG(row_length,rows,model_levels)                        &
#if defined(A17_2B)
     &,DELTAS_MERGE(row_length,rows,model_levels)                       &
#endif
     &,PSI(row_length,rows,model_levels)                                &
     &,PM10(row_length,rows,model_levels)                               &
                                                  !PM10 (ug m-3)
     &,PM2p5(row_length,rows,model_levels)                              &  
                                                  !PM2.5 (ug m-3)
     &,PM10_SO4 (row_length, rows, model_levels)                        &
     &,PM2p5_SO4(row_length, rows, model_levels)                        &
                                    !Sulphate contributions to PM concs.
     &,PM10_BC (row_length, rows, model_levels)                         &
     &,PM2p5_BC(row_length, rows, model_levels)                         &
                                    !Black carbon contrib. to PM concs.
     &,PM10_BB (row_length, rows, model_levels)                         &
     &,PM2p5_BB(row_length, rows, model_levels)                         &
                                    !Biomass aerosol contrib to PM concs.
     &,PM10_OCFF (row_length, rows, model_levels)                       &
     &,PM2p5_OCFF(row_length, rows, model_levels)                       &
                                    !OCFF contributions to PM concs.
     &,PM10_SOA (row_length, rows, model_levels)                        &
     &,PM2p5_SOA(row_length, rows, model_levels)                        &
                                    !SOA contributions to PM concs.
     &,PM10_SS (row_length, rows, model_levels)                         &
     &,PM2p5_SS(row_length, rows, model_levels)                         &
                                    !Sea-salt contributions to PM concs.
     &,PM10_DUST (row_length, rows, model_levels)                       &
     &,PM2p5_DUST(row_length, rows, model_levels)
                                    !Dust contributions to PM concs.
!
#include "csubmodl.h"
#include "typsts.h"
#include "c_0_dg_c.h"

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)     ! STASH workspace for section 17

! Local variables
      Integer                                                           &
     & i, j, level, k                                                   &
                             ! loop counters
     &,    icode             ! Return code  =0 Normal exit  >1 Error

      Integer sect,item    ! STASH section, item no.s
      Parameter (sect = 17) !  for aero_ctl (S Cycle or soot)

      Real work_3d(row_length,rows,model_levels) ! work space

      Character*80  cmessage

      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_aero')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     &  copydiag_3d                                                     &
     &, copydiag                                                        &
     &, ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)
!
! Copy diagnostic information to STASHwork for STASH processing
!
! Write MSA to STASH if DMS included
!
      If (L_SULPC_DMS) Then
!
        item = 203                          !MSA
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! Convert to flux per sec
          Do level=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                work_3d(i,j,level) = MSA(i,j,level)/timestep
              End Do
            End Do
          End Do
!
! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(item,sect,im_index)),          &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 203)"//cmessage
          End if
!
        END IF
!
      End If
!
!
! Write NH3_DEP to STASH if NH3 included
!
      If (L_SULPC_NH3) Then
!
        item = 204                          !NH3_DEP
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! Convert to flux per sec
          Do level=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                work_3d(i,j,level) = NH3_DEP(i,j,level)/timestep
              End Do
            End Do
          End Do
!
! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(item,sect,im_index)),          &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 204)"//cmessage
          End if
!
        END IF
!
      End If         ! End L_SULPC_NH3 condn
!
!
!
! Diagnose DMS emissions if requested
!
      If (L_DMS_em) Then
!
        item = 205                          !DMS emissions
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! DEPENDS ON: copydiag
          Call copydiag (stashwork(si(item,sect,im_index)),             &
     &        DMS_emiss,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": Error in copydiag (item 205)"//cmessage
          End if
!
        END IF
!
      End If         ! End L_DMS_em condn

      IF (L_SULPC_SO2 .AND. L_SULPC_DMS) THEN

        item = 206
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 206)"//cmessage
          ENDIF

        ENDIF

        item = 207
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &              F_DMS_TO_SO2(i,j,level) / timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 207)"//cmessage
          ENDIF

        ENDIF

        item = 208
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &             F_DMS_TO_SO4(i,j,level) *                            &
     &             (1.0E00 - PSI(i,j,level))/ timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 208)"//cmessage
          ENDIF

        ENDIF

        item = 209
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &             F_DMS_TO_SO4(i,j,level) *                            &
     &             PSI(i,j,level) / timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 209)"//cmessage
          ENDIF

        ENDIF

        item = 210
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 210)"//cmessage
          ENDIF

        ENDIF

        item = 211
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_WET(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 211)"//cmessage
          ENDIF

        ENDIF

        item = 212
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_WET_O3(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 212)"//cmessage
          ENDIF

        ENDIF

        item = 213
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_EVAP(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 213)"//cmessage
          ENDIF

        ENDIF

        item = 214
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_NUCL(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 214)"//cmessage
          ENDIF

        ENDIF

        item = 215
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DIFFUSE(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 215)"//cmessage
          ENDIF

        ENDIF

        item = 216
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_COAG(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 216)"//cmessage
          ENDIF

        ENDIF

        item = 217
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)*             &
     &             (1.0E00-PSI(i,j,level))/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 217)"//cmessage
          ENDIF

        ENDIF

        item = 218
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)*             &
     &             PSI(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 218)"//cmessage
          ENDIF

        ENDIF

#if defined(A17_2B)
        item = 219
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_MERGE(i,j,level)            &
     &             /timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 219)"//cmessage
          ENDIF

        ENDIF

#endif

      ENDIF
      
! ---------------------------------------------------------------------------
!
!     Diagnose PM10, PM2.5 and the contributions of the different
!     aerosols species to them if requested. Note that PM10 & PM2.5
!     can be calculated as long as any of the aerosol species is used. 
     
! Diagnose PM10 if requested:
!
      item = 220
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10,                                                     &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 220)"//cmessage
        End if
      END IF
   
! Diagnose PM2.5 if requested:
!
      item = 221
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5,                                                    &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 221)"//cmessage
        End if
      END IF

! Diagnose PM10 & PM2.5 concs. due to different aerosol species
! if requested:
!
      item = 222
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SO4,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 222)"//cmessage
        End if
      END IF
!
      item = 223
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SO4,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 223)"//cmessage
        End if
      END IF
!
      item = 224
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_BC,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 224)"//cmessage
        End if
      END IF
!
      item = 225
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_BC,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 225)"//cmessage
        End if
      END IF
!
      item = 226
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_BB,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 226)"//cmessage
        End if
      END IF
!
      item = 227
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_BB,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 227)"//cmessage
        End if
      END IF
!
      item = 228
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_OCFF,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 228)"//cmessage
        End if
      END IF
!
      item = 229
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_OCFF,                                               &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 229)"//cmessage
        End if
      END IF
!
      item = 230
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SOA,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 230)"//cmessage
        End if
      END IF
!
      item = 231
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SOA,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 231)"//cmessage
        End if
      END IF
!
      item = 232
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SS,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 232)"//cmessage
        End if
      END IF
!
      item = 233
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SS,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 233)"//cmessage
        End if
      END IF
!
      item = 234
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_DUST,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 234)"//cmessage
        End if
      END IF
!
      item = 235
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_DUST,                                               &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 235)"//cmessage
        End if
      END IF

! ---------------------------------------------------------------------------

      IF (icode > 0) THEN
! DEPENDS ON: ereport
        CALL ereport ('DIAGNOSTICS_AERO', icode, cmessage)
      END IF
     
      RETURN
      END SUBROUTINE diagnostics_aero
!
#endif
