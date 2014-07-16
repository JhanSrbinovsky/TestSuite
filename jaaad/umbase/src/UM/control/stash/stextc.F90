#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STEXTC ---------------------------------------------------
!LL
!LL  Purpose: Extracts a weighted subfield within a region specified
!LL           by a lower left hand and upper right hand corner.
!LL           Single level at a time. (STASH service routine).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Allow level-by-level mass-weighting if mass-weights
!LL                  are so defined, otherwise use P*.
!LL   4.3  06/01/97  Moved weighting and masking calculations up to
!LL                  SPATIAL.                              P.Burton
!LL   5.0  13/07/99  Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0  17/11/99  Removed row_length,p_rows arguments
!LL                                             P.Burton
!LL   6.1  16/08/04  Increase first dimension in certain arrays
!LL                  - avoids bank conflicts. Jens-Olaf Beismann
!LL                  NEC.   Anthony A. Dickinson
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D711
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STEXTC(fieldin,vx,vy,fld_type,halo_type,               &
     &                  lwrap,lmasswt,                                  &
     &                  xstart,ystart,xend,yend,                        &
     &                  fieldout,                                       &
     &                  pstar_weight,                                   &
     &                  area_weight,mask,                               &
     &                  level_code,mask_code,weight_code,rmdi,          &
     &                  icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &    vx,vy,                                                        &
                                                ! IN  input field size
     &    fld_type,                                                     &
                                                ! IN  field type(u/v/p)
     &    halo_type,                                                    &
                                                ! IN  halo type
     &    xstart,ystart,                                                &
                                                ! IN  lower LH corner
     &    xend,yend,                                                    &
                                                ! IN  upper RH corner
     &    level_code,                                                   &
                                                ! IN  input level code
     &    mask_code,                                                    &
                                                ! IN  masking code
     &    weight_code,                                                  &
                                                ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER*(*)                                                     &
     &    cmessage                              ! OUT error return msg
      LOGICAL                                                           &
     &    lwrap,                                                        &
                                                ! IN  TRUE if wraparound
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    mask(vx+1,vy)                         ! IN  mask array
      REAL                                                              &
     &    fieldin(vx,vy),                                               &
                                                ! IN  input field
     &    in_aux(vx+1,vy),                                              &
                                                ! IN  input field
     &    fieldout(xstart:xend,ystart:yend),                            &
                                                ! OUT output field
     &    out_aux(xstart:xend+1,ystart:yend),                           &
     &    pstar_weight(vx+1,vy),                                        &
                                                ! IN  pstar mass weight
! (already interpolated to the correct grid and
!  set to 1.0 where no mass weighting is required)
     &    area_weight(vx+1,vy),                                         &
                                                ! IN  area weighting
! (already interpolated to the correct grid and
!  set to 1.0 where no area weighting is required)
     &    rmdi                                  ! IN  missing data indic
!*----------------------------------------------------------------------
!
! External subroutines called
!
!
#include "stparam.h"
#include "sterr.h"
#include "parvars.h"
!
! Local variables
!
      INTEGER i,j,ii   ! ARRAY INDICES FOR VARIABLE


!L----------------------------------------------------------------------

! Calculate the output field, by multiplying the input field by
! pstar_weight and area_weight. These arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0

      in_aux(1:vx,:) = fieldin(1:vx,:)

      DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
        ELSE
          ii=i
        ENDIF
        DO j=ystart,yend
          IF (mask(ii,j)) THEN
              out_aux(i,j) =                                            &
     &          in_aux(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
          ELSE
            out_aux(i,j)=rmdi
          ENDIF
        ENDDO
      ENDDO

      fieldout(:,:) = out_aux(xstart:xend,ystart:yend)


      RETURN
      END SUBROUTINE STEXTC
#endif
