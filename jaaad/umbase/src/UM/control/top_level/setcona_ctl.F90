#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SETCONA_CTL
!
! Purpose: Interface routine to SETCONA to allow SETCONA to be called
!          from INITAL. Also performs some intialisation of the
!          atmosphere STASH array that was previously carried out
!          in INITDUMP
!  Model            Modification history from model version 5.2:
! version  Date
! 5.2      13/09/00 New deck introduced so that SETCONA can be moved
!                   from INITDUMP and called from INITIAL after the
!                   intital LBCs have been read in.
!                                                        Paul Burton
!                   Replace eta_levels by zseak and Ck arrays for
!                   passing into STASH superarray. R Rawlins
! 5.3      31/10/01 Remove RimWeightsA_Orog. D.Robinson
! 5.4      02/09/02 Set up sea mask file in a_spsts(a_ixsts(12)).
!                   This is not the opposite of the land mask when
!                   coastal tiling is used.               K.Williams
! 6.1      24/06/04  Add tracer array sizes for super_array_sizes
!                    calculation                       A. Malcolm
!  6.1   04/08/04 Need global_rows to set diffusion coefficients
!                 in SETCONA. Add global_row_length for future.
!                                                         Terry Davies
!

      SUBROUTINE SETCONA_CTL(                                           &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argcona.h"
#include "arglndm.h"
     &  ICODE,CMESSAGE,isSTASH)     ! FLUME-STASH UM vn7.1

Use trignometric_mod, Only : cos_v_latitude, cos_theta_latitude

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "typcona.h"
#include "typlndm.h"
#include "cbound.h"
#include "cntlatm.h"
#include "atm_lsm.h"

      LOGICAL       isSTASH  ! FLUME-STASH UM vn7.1
      INTEGER                                                           &
     &  ICODE      ! Return code

      CHARACTER*80                                                      &
     &  CMESSAGE   ! Error message

! Local variables

      INTEGER                                                           &
     &  i,j,ij                                                          &
                   ! loop counters
     &, dummy      ! Dummy integer

      REAL                                                              &
     &  fland(theta_field_size)                                         &
                                 ! Land fraction un-compressed
     &, rsea                     ! Real equiv of logical lsea

      LOGICAL                                                           &
     &  lsea                     ! True if any of gridbox is sea

      EQUIVALENCE (rsea,lsea)

!----------------------------------------------------------------------

! 1.0 Call SETCONA

! DEPENDS ON: setcona
      CALL SETCONA( A_LEVDEPC(jetatheta), A_LEVDEPC(jetarho),           &
     &              D1(jvol_smc_sat), D1(jland), D1(jorog),             &
     &              D1(jorog_grad_x), D1(jorog_grad_y),                 &
     &              D1(jrho(1)), D1(jexner_rho_levels(1)),              &
     &              D1(jorog_lbc),d1(jexner_lbc),                       &
     &              LENRIMA,LBC_SIZEA,LBC_STARTA,                       &
     &              RIMWIDTHA, RIMWEIGHTSA,                             &
     &              global_row_length, global_rows,                     &
     &              MODEL_LEVELS, ROWS, N_ROWS, ROW_LENGTH,             &
     &              LAND_FIELD, WET_LEVELS, BL_LEVELS,                  &
     &              A_INTHD(ih_1_c_rho_level), CLOUD_LEVELS,            &
     &              A_REALHD(rh_z_top_theta),                           &
     &              tr_levels, tr_vars, tr_ukca,                        &
     &              A_INTHD(ih_height_gen),                             &
     &              A_REALHD(rh_deltaEW), A_REALHD(rh_deltaNS),         &
     &              A_REALHD(rh_baselat), A_REALHD(rh_baselong),        &
     &              A_REALHD(rh_rotlat), A_REALHD(rh_rotlong),          &
     &              A_COLDEPC(jlambda_input_p),                         &
     &              A_COLDEPC(jlambda_input_u),                         &
     &              A_ROWDEPC(jphi_input_p),                            &
     &              A_ROWDEPC(jphi_input_v),                            &
     &              A_FIXHD(fh_RowDepCStart),                           &
     &              D1(jexner_theta_levels(1)),                         &
     &              D1(jp_theta_levels(1)),                             &
     &              D1(jp(1)), D1(jpstar),                              &
#include "argcona.h"
#include "arglndm.h"
     &              ICODE, CMESSAGE, isSTASH)  ! FLUME-STASH UM vn7.1

! 2.0 Populate atmosphere STASH 'super' array - this is required to
!     pass atmosphere sub-model specific information into STASH lower
!     level routines in a generic way.

      DO i=1,MODEL_LEVELS
        a_spsts(a_ixsts(1) + i-1) = a_levdepc(Jzseak_rho+i-1)
        a_spsts(a_ixsts(2) + i-1) = a_levdepc(JCk_rho   +i-1)
      ENDDO ! i

      DO i=0,MODEL_LEVELS
        a_spsts(a_ixsts(3) + i  ) = a_levdepc(Jzseak_theta+i)
        a_spsts(a_ixsts(4) + i  ) = a_levdepc(JCk_theta   +i)
      ENDDO ! i

      ! Indicies 5-6 not currently used in atmosphere STASH array

      ! !!! NOTE AND BEWARE !!!
      ! Ingeger pointers passed as real variables
      a_spsts(a_ixsts(7)) = jp(1)
      a_spsts(a_ixsts(8)) = jpstar  ! not clear if this is still needed

      ! Note: No halos for cos_v_latitude passed into a_spsts
      DO j=1,n_rows
        DO i=1,row_length
          ij= (i-1) + (j-1)*row_length
          a_spsts(a_ixsts(9)+ij)=cos_v_latitude(i,j)
        ENDDO ! i
      ENDDO ! j

      ! Note: No halos for cos_theta_latitude passed into a_spsts
      DO j=1,rows
        DO i=1,row_length
          ij= (i-1) + (j-1)*row_length
          a_spsts(a_ixsts(10)+ij)=cos_theta_latitude(i,j)
        ENDDO ! i
      ENDDO ! j

      DO i=1,theta_field_size
        a_spsts(a_ixsts(11) + i-1 ) = d1(jland + i-1)
      ENDDO

!     Set up sea mask. Equivalent lsea and rsea required as a_spsts
!     is a real array and it needs to be filled with a logical.

      IF (l_ctile) then
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS(fland,d1(jfrac_land),                     &
     &   atmos_landmask_local, theta_field_size, dummy)
        DO i=1,theta_field_size
          lsea = (fland(i) /= 1.0)
          a_spsts(a_ixsts(12) + i-1 ) = rsea
        ENDDO
      ELSE
        DO i=1,theta_field_size
          lsea = (.NOT.ld1(jland + i-1))
          a_spsts(a_ixsts(12) + i-1 ) = rsea
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SETCONA_CTL
#endif
