#if defined(CONTROL) || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Program: DERVSIZE -------------------------------------------------
!LL
!LL  Purpose: Calculate extra sizes required for dynamic allocation of
!LL           main memory in the model, derived from sizes passed by
!LL           READSIZE into the top level program UM_SHELL. These are
!LL           local sizes for each pe, except where explicitly stated,
!LL           having previously called decomposition routines.
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation: On-line UM document C1 - Dynamic allocation
!LL                          of primary fields
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE DERVSIZE(                                              &
     &             ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      IMPLICIT NONE
!
!  Argument list and comdecks
!
      INTEGER ICODE             ! OUT - Return code
      CHARACTER*80 CMESSAGE     ! OUT - Error message

#include "parvars.h"
#include "typsize.h"
! For STASH sizes [mos_mask_len]
#include "typstsz.h"
#include "cntlatm.h"
#include "cmaxsize.h"
#include "decomptp.h"
#include "decompdb.h"
#include "cprintst.h"
!  Local variables
       integer numside_rowso  ! # of rows in each ocean bdy field
       integer numside_colso  ! # of columns in each ocean bdy field
      INTEGER                                                           &
     &  iside                                                           &
                     ! loop counter for LBC sides
     &, ifld                                                            &
                     ! loop counter for field types
     &, ihalo                                                           &
                     ! loop counter for halo types
     &, iproc                                                           &
                     ! loop counter for processors
     &, irim                                                            &
                     ! loop counter for rim types
     &, info                                                            &
                     ! return code from GCOM
     &, lbc_row_len                                                     &
                     ! length of row of LBC
     &, lbc_nrows                                                       &
                     ! number of rows in LBC
     &, num_optional_lbcs_in                                            &
                              ! no. of optional lbc fields in input
     &, num_optional_lbcs_out ! no. of optional lbc fields output

!*----------------------------------------------------------------------

      INTEGER nohalo_IMT,nohalo_JMT,glob_IMT,glob_JMT

      ICODE=0

!     Initialise number of optional in/out lbcs to zero
      num_optional_lbcs_in  = 0
      num_optional_lbcs_out = 0
!L
!L   Atmosphere Boundary Datasets.
!L   2nd dimension of Level Dependent Constants.
      INTF_LEN2_LEVDEPC=4
!L   2nd dimension of Row/Col Dependent Constants.
      INTF_LEN2_ROWDEPC=2
      INTF_LEN2_COLDEPC=2      
!L
!L   Sizes applicable to all resolutions
#if defined(ATMOS)
      THETA_FIELD_SIZE=ROW_LENGTH*ROWS
! One less V row on the Northern most processors
      IF (at_extremity(PNorth) .and.                                    &
     &         model_domain  /=  mt_bi_cyclic_lam) THEN
        N_ROWS=ROWS-1
      ELSE
        N_ROWS=ROWS
      ENDIF
      U_FIELD_SIZE=THETA_FIELD_SIZE
      V_FIELD_SIZE=ROW_LENGTH*N_ROWS
      theta_off_size   = (row_length + 2*offx)   * (rows   + 2*offy)
      theta_halo_size  = (row_length + 2*halo_i) * (rows   + 2*halo_j)
      u_off_size       = (row_length + 2*offx)   * (rows   + 2*offy)
      u_halo_size      = (row_length + 2*halo_i) * (rows   + 2*halo_j)
      v_off_size       = (row_length + 2*offx)   * (n_rows + 2*offy)
      v_halo_size      = (row_length + 2*halo_i) * (n_rows + 2*halo_j)


!     No of levels for Convective Cloud Amount.
      IF (L_3D_CCA .OR. L_CCRAD) THEN
        N_CCA_LEV = WET_LEVELS
      ELSE
        N_CCA_LEV = 1
      ENDIF
      IF(PrintStatus >= PrStatus_Normal) THEN
        WRITE(6,*)                                                      &
     &  'DERVSIZE: Number of levels for convective clouds is ',         &
     &  N_CCA_LEV
      ENDIF

      ! ----------------------------------------------------------------
      ! Count number of optional lateral boundary fields expected in
      ! input dependent on whether the _lbc logicals are true or false
      ! ----------------------------------------------------------------

      ! Additional microphysics variables (ice crystals, rain, graupel)
      If (L_mcr_qcf2_lbc) Then  ! qcf2 lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If
      If (L_mcr_qrain_lbc) Then  ! qrain lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If
      If (L_mcr_qgraup_lbc) Then  ! qgraup lbcs active
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      End If

      ! Cloud fractions for PC2 (3 fields: bulk, liquid and frozen)
      If (L_pc2_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 3
      EndIf

      ! Murk aerosol
      If (L_murk_lbc) Then
        num_optional_lbcs_in = num_optional_lbcs_in + 1
      EndIf

      ! ----------------------------------------------------------------
      ! Count number of optional lateral boundary fields to write out
      ! dependent on whether the prognostics are active or not
      ! ----------------------------------------------------------------

      ! Additional microphysics variables (ice crystals, rain, graupel)
      If (L_mcr_qcf2) Then  ! qcf2 lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If
      If (L_mcr_qrain) Then  ! qrain lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If
      If (L_mcr_qgraup) Then  ! qgraup lbcs active
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      End If

      ! Cloud fractions for PC2 (3 fields: bulk, liquid and frozen)
      If (L_pc2) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 3
      EndIf

      ! Murk aerosol
      If (L_murk) Then
        num_optional_lbcs_out = num_optional_lbcs_out + 1
      EndIf

      ! ----------------------------------------------------------------


      A_LEN1_LEVDEPC=MODEL_LEVELS+1
! We use the global values here
      A_LEN1_ROWDEPC= global_rows
      A_LEN1_COLDEPC= global_row_length
      A_LEN1_FLDDEPC= global_rows * global_row_length
      MOS_MASK_LEN  = global_rows * global_row_length

!L Number of atmosphere model interface lookups
      INTF_LOOKUPSA = 13 + num_optional_lbcs_out + TR_VARS


!     No of LBC variables in 'old' LBC files
!     QCF LBC not catered for. Not used in 4.5 LBC files.

      OLD_INTF_LOOKUPSA = 5+TR_VARS

#else
      MOS_MASK_LEN=1
#endif
!
! Calculate sizes of Lateral Boundary Conditions Data:
! LENRIMA(fld_type,halo_type,rim_type)
!    size of a single level of LBC data for a given field type and
!    halo type and rim_type
! LBC_SIZEA(side,fld_type,halo_type,rim_type)
!    size of a single edge of LBC data for a given edge (North,
!    East,South,West) and a given field type and halo type and rim type
! LBC_STARTA(side,fld_type,halo_type,rim_type)
!    start address in LBC array of a given edge for a given
!    field type and halo type and rim_type
!

      DO irim=1,Nrima_max
      DO ifld=1,Nfld_max        ! loop over field types
        DO ihalo=1,NHalo_max    ! loop over halo types

          LENRIMA(ifld,ihalo,irim)=0
          global_LENRIMA(ifld,ihalo,irim)=0

          DO iside=1,4          ! loop over North,East,South,West

            DO iproc=0,nproc-1
              g_LBC_STARTA(iside,ifld,ihalo,irim,iproc)=0
            ENDDO
! First calculate the global_LENRIMA values - these are the sizes
! of the complete LBC as stored on disk before it is decomposed
! over processors.

            IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
! North or South boundaries

              IF (ifld  ==  fld_type_u) THEN
                lbc_row_len=glsize(1,ifld)-1
              ELSE
                lbc_row_len=glsize(1,ifld)
              ENDIF

              lbc_row_len=lbc_row_len+2*halosize(1,ihalo)

              lbc_nrows=halosize(2,ihalo)+RIMWIDTHA(irim)


            ELSE
! East or West boundaries
              lbc_row_len=halosize(1,ihalo)+RIMWIDTHA(irim)

              lbc_nrows=glsize(2,ifld)-2*RIMWIDTHA(irim)

            ENDIF ! North/South or East/West boundary

            IF (RIMWIDTHA(irim)  >   0) THEN

              global_LBC_SIZEA(iside,ifld,ihalo,irim)=                  &
     &          lbc_row_len*lbc_nrows

              IF (iside  ==  1) THEN
                global_LBC_STARTA(iside,ifld,ihalo,irim)=1
              ELSE
                global_LBC_STARTA(iside,ifld,ihalo,irim)=               &
     &            global_LBC_STARTA(iside-1,ifld,ihalo,irim)+           &
     &            global_LBC_SIZEA(iside-1,ifld,ihalo,irim)
              ENDIF

              global_LENRIMA(ifld,ihalo,irim)=                          &
     &          global_LENRIMA(ifld,ihalo,irim)+                        &
     &          global_LBC_SIZEA(iside,ifld,ihalo,irim)

            ELSE ! No LBCs if RIMWIDTH is  <=  0)

               global_LBC_SIZEA(iside,ifld,ihalo,irim)=0
               global_LBC_STARTA(iside,ifld,ihalo,irim)=1

            ENDIF ! IF (RIMWIDTHA  >   0)

! Now calculate local LENRIMA values, and the associated arrays
! LBC_SIZEA and LBC_STARTA

            IF (at_extremity(iside) .AND.                               &
     &          (RIMWIDTHA(irim)  >   0)) THEN
! This processor is at the edge of the grid

              IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
! North or South boundaries
! North/South boundaries can contain the corners of the LBCs

! Length of rows includes the halos.
! For U fields, there is one less point as the last point along
! each row is ignored. This only applies to the last processor along
! the row
                IF ((ifld  ==  fld_type_u) .AND.                        &
     &              (at_extremity(PEast))) THEN
                  lbc_row_len=lasize(1,ifld,ihalo)-1
                ELSE
                  lbc_row_len=lasize(1,ifld,ihalo)
                ENDIF

! And the number of rows is the size of the halo plus the rimwidth
                lbc_nrows=halosize(2,ihalo)+RIMWIDTHA(irim)

              ELSE
! East or West boundaries

! Length of row is the size of the halo plus the rimwidth
                lbc_row_len=halosize(1,ihalo)+RIMWIDTHA(irim)

! Number of rows is the number of rows on this PE minus any
! rows which are looked after by the North/South boundaries
! (ie. the corners).
                lbc_nrows=lasize(2,ifld,ihalo)
                IF (at_extremity(PNorth))                               &
     &            lbc_nrows=lbc_nrows-halosize(2,ihalo)-RIMWIDTHA(irim)
                IF (at_extremity(PSouth))                               &
     &            lbc_nrows=lbc_nrows-halosize(2,ihalo)-RIMWIDTHA(irim)

              ENDIF ! North/South or East/West boundary

              LBC_SIZEA(iside,ifld,ihalo,irim)=lbc_row_len*lbc_nrows

            ELSE
! This processor is not at the edge of the grid, or RIMWIDTH is
! zero (indicating no LBCs)
              LBC_SIZEA(iside,ifld,ihalo,irim)=0

            ENDIF

! LBC_STARTA contains the offset in the LBC array for each side
! (North,East,South,West) piece of data

            IF (iside  ==  1) THEN
              LBC_STARTA(iside,ifld,ihalo,irim)=1
            ELSE
              LBC_STARTA(iside,ifld,ihalo,irim)=                        &
     &          LBC_STARTA(iside-1,ifld,ihalo,irim)+                    &
     &          LBC_SIZEA(iside-1,ifld,ihalo,irim)
            ENDIF
            g_LBC_STARTA(iside,ifld,ihalo,irim,mype)=                   &
     &        LBC_STARTA(iside,ifld,ihalo,irim)

           LENRIMA(ifld,ihalo,irim)=LENRIMA(ifld,ihalo,irim)+           &
     &                              LBC_SIZEA(iside,ifld,ihalo,irim)

          ENDDO ! iside
        ENDDO ! ihalo
      ENDDO ! ifld
      ENDDO ! irim

! Now do some comms so that g_LBC_STARTA is filled with the
! value of LBC_STARTA on every processor

      CALL GC_IMAX(4*Nfld_max*NHalo_max*Nrima_max*nproc,nproc,info,     &
     &             g_LBC_STARTA)

! And set up a few other variables

      ! Includes one-off orography field at start
      RIM_LOOKUPSA = 13 + num_optional_lbcs_in
      BOUND_LOOKUPSA=RIM_LOOKUPSA+(NRIM_TIMESA-1)*(RIM_LOOKUPSA-1)
                        ! Total number of headers in the LBC file
!LL LAM DERIVED SIZES END
      RETURN
      END SUBROUTINE DERVSIZE
#endif
