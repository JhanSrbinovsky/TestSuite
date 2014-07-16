#if defined(CONTROL) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Determine STASH input length per vertical level for prog var
! Subroutine Interface:

      SUBROUTINE ADDRLN(IGPL,halo_type,LEN,size_type)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "version.h"
#include "parvars.h"
#include "typsize.h"
#include "cntlatm.h"
#include "cntlocn.h"
#include "model.h"
#include "grdtypes.h"
#include "cppxref.h"

! Subroutine arguments:
      INTEGER                                                           &
     &  IGPL                                                            &
                           ! IN : grid code
     &, halo_type                                                       &
                           ! IN : type of halo (none, single
                           !                    or extended)
     &, size_type                                                       &
                           ! IN : "local_data" or "global_data"
     &, LEN                ! OUT : length of field

! Local variables

      INTEGER                                                           &
     &  IX1                                                             &
                        ! Column number of start of area
     &, IX2                                                             &
                        ! Column number of end of area
     &, IY1                                                             &
                        ! Row number at start of area
     &, IY2                                                             &
                        ! Row number at end of area

     &, local_IX1                                                       &
                        ! local IX1 for this processor
     &, local_IX2                                                       &
                        ! local IX2 for this processor
     &, local_IY1                                                       &
                        ! local IY1 for this processor
     &, local_IY2                                                       &
                        ! local IY2 for this processor

     &, fld_type                                                        &
                        ! Which kind of variable (theta u or v?)

                        ! Variables from GT_DECODE:
     &, MODEL_TYPE                                                      &
                        ! model type of grid
     &, CONTENT                                                         &
                        ! content type of grid
     &, COVERAGE                                                        &
                        ! coverage of grid
     &, DOMAIN                                                          &
                        ! domain of grid
     &, CYCLIC                                                          &
                        ! does grid contain cyclic wrap columns

     &, rim_type                                                        &
                        ! Which type of rimwidth - normal or orography?
     &, ocean_extra_pts ! extra points added to ocean row length
!                       ! to allow for the wrap around pts at the
!                       ! start and end of each global row

! Functions
      INTEGER                                                           &
     &  GET_FLD_TYPE

! Externals

      EXTERNAL                                                          &
     &  GT_DECODE                                                       &
     &, LLTORC


!- End of Header ---------------------------------------------------

! Get information about grid type

! DEPENDS ON: gt_decode
      CALL GT_DECODE(IGPL,                                              &
     &               MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)

! Determine row/column nos. for global domain
! DEPENDS ON: lltorc
      CALL LLTORC(IGPL,90,-90,0,360,IY1,IY2,IX1,IX2)

      IF ((size_type  ==  local_data) .AND.                             &
     &    ((DOMAIN  /=  gt_compressed) .AND.                            &
     &     (DOMAIN  /=  gt_LBC))) THEN

! Convert the global subdomain limits to local subdomain limits
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE.,.TRUE.,                  &
     &                                  IGPL,halo_type,mype,            &
     &                                  IY1,IX2,IY2,IX1,                &
     &                                  local_IY1,local_IX2,            &
     &                                  local_IY2,local_IX1)


        IX1=local_IX1
        IX2=local_IX2
        IY1=local_IY1
        IY2=local_IY2

        ocean_extra_pts=0
        IF (at_extremity(PWest)) ocean_extra_pts=ocean_extra_pts+1
        IF (at_extremity(PEast)) ocean_extra_pts=ocean_extra_pts+1
      ELSE
        ocean_extra_pts=2 ! extra pt at start and end of row for
!                         ! wrap around
      ENDIF

      IF (DOMAIN  ==  gt_compressed) THEN

        IF (MODEL_TYPE  ==  gt_atmos) THEN

          IF (size_type  ==  local_data) THEN
            len=local_land_field
          ELSE
            len=global_land_field
          ENDIF

        ELSEIF (MODEL_TYPE  ==  gt_ocean) THEN

          len=-1 ! Set flag at this stage for a multi-level
                 ! compress

        ENDIF ! MODEL_TYPE

      ELSEIF (DOMAIN  ==  gt_ozone) THEN

        IF (ZonAvOzone) THEN ! zonal

          IF (size_type  ==  local_data) THEN
            LEN=IY2-IY1+1+2*halosize(2,halo_type)
          ELSE
            LEN=IY2-IY1+1
          ENDIF

        ELSE ! Full fields

          LEN=(IX2-IX1+1)*(IY2-IY1+1)

        ENDIF

      ELSEIF (DOMAIN  ==  gt_LBC) THEN

        IF (MODEL_TYPE  ==  gt_atmos) THEN

! DEPENDS ON: get_fld_type
          fld_type=GET_FLD_TYPE(IGPL)
          IF (IGPL  ==  ppx_atm_lbc_orog) THEN
            rim_type=rima_type_orog
          ELSE
            rim_type=rima_type_norm
          ENDIF

! NB - these sizes are just for a single level. The LBCs are stored
! with all vertical levels in one field. This number will need to
! be multiplied by the number of levels by the calling routine

          IF (size_type  ==  local_data) THEN
            LEN=LENRIMA(fld_type,halo_type,rim_type)
          ELSE
            LEN=global_LENRIMA(fld_type,halo_type,rim_type)
          ENDIF

        ELSEIF (MODEL_TYPE  ==  gt_ocean) THEN

         LEN=-1
         WRITE(6,*)'ADDRLN1: WARNING - ocean boundary data is not'
         WRITE(6,*)'stored in the dump from UM version 5.3'


        ELSEIF (MODEL_TYPE  ==  gt_wave) THEN

          IF (size_type  ==  local_data) THEN
            LEN=-1  ! Not yet supported
          ELSE
            LEN=-1
          ENDIF

        ENDIF

      ELSEIF (MODEL_TYPE  ==  gt_ocean) THEN

        IF (COX_O .AND. CYCLIC == gt_optcyclic) THEN

          IF (size_type  ==  local_data) THEN
            ! Primary fields are overdimensioned by one row
            ! to make them the same size as the mass grid

! DEPENDS ON: get_fld_type
            fld_type=GET_FLD_TYPE(IGPL)

            IF (at_extremity(PNorth) .AND.                              &
     &          (fld_type  ==  fld_type_u)) THEN
               LEN=(IX2-IX1+1+ocean_extra_pts)*                         &
     &             (IY2-IY1+2+2*halosize(2,halo_type))
            ELSE
               LEN=(IX2-IX1+1+ocean_extra_pts)*                         &
     &             (IY2-IY1+1+2*halosize(2,halo_type))
            ENDIF

          ELSE
            LEN=(IX2-IX1+1+ocean_extra_pts)*(IY2-IY1+1)
          ENDIF

        ELSEIF (CYCLIC  ==  gt_cyclic) THEN

          IF (size_type  ==  local_data) THEN

! DEPENDS ON: get_fld_type
            fld_type=GET_FLD_TYPE(IGPL)

            IF (at_extremity(PNorth) .AND.                              &
     &        (fld_type  ==  fld_type_u)) THEN

              LEN=(IX2-IX1+1+2*halosize(1,halo_type))*                  &
     &            (IY2-IY1+2+2*halosize(2,halo_type))

            ELSE

              LEN=(IX2-IX1+1+2*halosize(1,halo_type))*                  &
     &            (IY2-IY1+1+2*halosize(2,halo_type))

            ENDIF

          ELSE ! Full model size
            LEN=(IX2-IX1+1)*(IY2-IY1+1)
          ENDIF

        ELSE ! Ocean catch-all
          LEN=(IX2-IX1+1)*(IY2-IY1+1)
        ENDIF

      ELSE

        LEN =(IX2-IX1+1)*(IY2-IY1+1)

      ENDIF

 9999 CONTINUE

      RETURN
      END SUBROUTINE ADDRLN

!- End of subroutine code ------------------------------------------
#endif
