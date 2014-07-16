#if defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine DPM_RPM ------------------------------------------------
!!!
!!! Purpose : Calculates the DPM_RPM ratio of litter input for use
!!!            in the RothC soil carbon sub-model
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  6.2     01/03/06     New Deck. Chris Jones.
!!!
!!!END ----------------------------------------------------------------

        SUBROUTINE DPM_RPM(LAND_PTS,TRIF_PTS,TRIF_INDEX,FRAC,           &
     &                     FRAC_AGRIC,LIT_C,DPM_RATIO)


      IMPLICIT NONE

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                ! IN Total number of land points.
     &,TRIF_PTS                                                         &
                                  ! IN Number of points on which
!                                 !    TRIFFID may operate
     &,TRIF_INDEX(LAND_PTS)                                             &
                                ! IN Indices of land points on
!                                 !    which TRIFFID may operate
     &,L,T                        ! WORK Loop counters

      REAL                                                              &
     & FRAC(LAND_PTS,NTYPE)                                             &
                                ! IN  Fractional cover of each
     &,FRAC_AGRIC(LAND_PTS)                                             &
                                ! IN  Agricultural (disturbed) frac
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                ! IN  PFT carbon litter
!                                 !    (kg C/m2/360days).
     &,DPM_RATIO(LAND_PTS)      ! OUT DPM ratio of litter input

      REAL                                                              &
     & grass(LAND_PTS)                                                  &
                                ! WORK  fractional grass cover
     &,crop(LAND_PTS)                                                   &
                                ! WORK  fractional crop cover
     &,Ltree(LAND_PTS)                                                  &
                                ! WORK  litter input from trees
     &,Lshrub(LAND_PTS)                                                 &
                                ! WORK  litter input from shrubs
     &,Lgrass_crop(LAND_PTS)                                            &
                                ! WORK  litter input from grass+crop
     &,Lgrass(LAND_PTS)                                                 &
                                ! WORK  litter input from grass
     &,Lcrop(LAND_PTS)                                                  &
                                ! WORK  litter input from crops
     &,DPM(LAND_PTS)            ! WORK  litter input to DPM

      REAL                                                              &
     & rt                                                               &
                                  ! DPM:RPM ratio for trees
     &,rs                                                               &
                                  ! DPM:RPM ratio for shrubs
     &,rg                                                               &
                                  ! DPM:RPM ratio for grass
     &,rc                         ! DPM:RPM ratio for crops

      PARAMETER(rt=0.25, rs=0.33, rg=0.67, rc=1.44)

! calculate grass/crop fractions from C3, C4 cover, and FRAC_AGRIC
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

        grass(L) = FRAC(L,3) + FRAC(L,4)
        if (grass(L)  >   FRAC_AGRIC(L)) then
          crop(L)  = FRAC_AGRIC(L)
          grass(L) = grass(L) - FRAC_AGRIC(L)
        else
          crop(L)  = grass(L)
          grass(L) = 0.0
        endif

      ENDDO

! calculate different litter-type amounts
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

        Ltree(L)       = LIT_C(L,1)*FRAC(L,1) + LIT_C(L,2)*FRAC(L,2)
        Lshrub(L)      = LIT_C(L,5)*FRAC(L,5)
        Lgrass_crop(L) = LIT_C(L,3)*FRAC(L,3) + LIT_C(L,4)*FRAC(L,4)
        Lgrass(L)      = Lgrass_crop(L)*grass(L) / (grass(L)+crop(L))
        Lcrop(L)       = Lgrass_crop(L)*crop(L) / (grass(L)+crop(L))

      ENDDO

! calculate DPM litter input and hence ratio
      DO T=1,TRIF_PTS
        L=TRIF_INDEX(T)

        DPM(L) = Ltree(L)  * rt/(1+rt) + Lshrub(L)*rs/(1+rs) +          &
     &           Lgrass(L) * rg/(1+rg) +  Lcrop(L)*rc/(1+rc)
        DPM_ratio(L) = DPM(L) /                                         &
     &               (Ltree(L) + Lshrub(L) + Lgrass(L) + Lcrop(L))

      ENDDO


      RETURN
      END SUBROUTINE DPM_RPM
#endif
