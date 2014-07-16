#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_SUBDOMAIN: converts global subdomain boundaries
!                            to local subdomain boundaries
! GLOBAL_TO_LOCAL_RC: converts global row,column co-ordinates to
!                     processor co-ordinates plus local
!                     co-ordinates within the processor.
!
! Subroutine Interface:

! Subroutine Interface:

! Function Interface
      INTEGER FUNCTION GET_FLD_TYPE (grid_type_code)

      IMPLICIT NONE

!
! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/11/96 New deck created for MPP code.  P.Burton
!  5.0      22/6/99  Changed field type from p,uv to p,u,v
!  5.1      27/03/00 Deal with B-grid UV variables (for diagnostics)
!                    Add the new LBC variables types.       P.Burton
!  5.2      19/09/00 Added new ppx_atm_lbc_orog grid type   P.Burton
!  5.5      10/08/00 Modification for parallelisation of WAM.
!                 Author: Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  5.5      15/01/03 River routing support. P.Selwood.
!
! Subroutine arguments:

      INTEGER                                                           &
     &  grid_type_code     ! IN : STASH grid type code

! Parameters
#include "cppxref.h"
#include "parparm.h"

      IF ( (grid_type_code  ==  ppx_atm_tall) .OR.                      &
     &     (grid_type_code  ==  ppx_atm_tland) .OR.                     &
     &     (grid_type_code  ==  ppx_atm_tsea) .OR.                      &
     &     (grid_type_code  ==  ppx_atm_tzonal) .OR.                    &
     &     (grid_type_code  ==  ppx_atm_tmerid) .OR.                    &
     &     (grid_type_code  ==  ppx_atm_compressed) .OR.                &
     &     (grid_type_code  ==  ppx_atm_ozone) .OR.                     &
     &     (grid_type_code  ==  ppx_ocn_tall)    .OR.                   &
     &     (grid_type_code  ==  ppx_ocn_tfield)  .OR.                   &
     &     (grid_type_code  ==  ppx_ocn_tzonal)  .OR.                   &
     &     (grid_type_code  ==  ppx_atm_lbc_theta)  .OR.                &
     &     (grid_type_code  ==  ppx_atm_lbc_orog)  .OR.                 &
     &     (grid_type_code  ==  ppx_ocn_lbc_theta)  .OR.                &
     &     (grid_type_code  ==  ppx_wam_all) .OR.                       &
     &     (grid_type_code  ==  ppx_ocn_tmerid) ) THEN
        GET_FLD_TYPE=fld_type_p
      ELSEIF                                                            &
     &   ( (grid_type_code  ==  ppx_atm_cuall) .OR.                     &
     &     (grid_type_code  ==  ppx_ocn_uall) .OR.                      &
     &     (grid_type_code  ==  ppx_ocn_cuall) .OR.                     &
     &     (grid_type_code  ==  ppx_ocn_ufield) .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_uzonal) .OR.                    &
     &     (grid_type_code  ==  ppx_atm_lbc_u)  .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_lbc_u)  .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_umerid) ) THEN
        GET_FLD_TYPE=fld_type_u
      ELSEIF                                                            &
     &   ( (grid_type_code  ==  ppx_atm_cvall) .OR.                     &
     &     (grid_type_code  ==  ppx_atm_lbc_v)  .OR.                    &
     &     (grid_type_code  ==  ppx_atm_uall ) .OR.                     &
     &     (grid_type_code  ==  ppx_ocn_cvall) ) THEN
        GET_FLD_TYPE=fld_type_v
      ELSEIF                                                            &
     &   ( (grid_type_code  ==  ppx_atm_uall) .OR.                      &
     &     (grid_type_code  ==  ppx_atm_uland) .OR.                     &
     &     (grid_type_code  ==  ppx_atm_usea) .OR.                      &
     &     (grid_type_code  ==  ppx_atm_uzonal) .OR.                    &
     &     (grid_type_code  ==  ppx_atm_umerid) .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_uall) .OR.                      &
     &     (grid_type_code  ==  ppx_ocn_ufield) .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_uzonal) .OR.                    &
     &     (grid_type_code  ==  ppx_ocn_umerid) ) THEN
        GET_FLD_TYPE=fld_type_v  ! This is actually the B U/V (velocity)
                                 ! grid, but it has the same sizes as
                                 ! the C V grid.
      ELSEIF                                                            &
     &     (grid_type_code  ==  ppx_wam_sea) THEN
        GET_FLD_TYPE=fld_type_comp_wave
      ELSEIF                                                            &
     &     (grid_type_code  ==  ppx_wam_rim) THEN
        GET_FLD_TYPE=fld_type_rim_wave
      ELSEIF                                                            &
     &   ( (grid_type_code  ==  ppx_atm_river) ) THEN
        GET_FLD_TYPE=fld_type_r  ! River routing grid

      ELSE
        GET_FLD_TYPE=fld_type_unknown
      ENDIF

      RETURN

      END FUNCTION GET_FLD_TYPE

#endif
