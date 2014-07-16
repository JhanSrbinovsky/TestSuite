#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE PROCESS_MET(p0,p_th,t0,tauu,tauv,hf,adr,ads,bl,        &
     &  land_fraction,cc,acr,acs,sice,sa,orog,gsf,lai_ft,canht,z0,gc,   &
     &  o3um,t0tile,sicetemp,soilmc,p,lnp,acp,adp,ra,rq,so4_vd)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Process met data
!-
!-   Inputs  : met data ...
!-   Outputs : met data ...
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    14/12/98  Created.  D.S. Stevenson
!  5.3    13/08/01  New dynamics version.  C.E. Johnson
!  5.3    17/11/01  Reinstated TROP. CB & CT removed. W.J. Collins
!  5.3    30/11/01  TROP now pressure not eta. W.J. Collins
!  6.1    21/10/04  Minor tidying of code. M.G. Sanderson
!  6.2   02/03/06  Extensive changes for new dry deposition scheme.
!                  Calls new subroutine AEROD to calculate ra, rq,
!                  so4_vd, and sets up global surface fraction array
!                  gsf. M.G. Sanderson
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!----------------------------------------------------------------------

      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: p0
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),   INTENT(IN)    :: p_th
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: t0
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: tauu
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: tauv
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: hf
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: adr
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: ads
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN)    :: bl
      REAL, DIMENSION(nlonpe,nlatpe),    INTENT(INOUT) :: land_fraction
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),   INTENT(INOUT) :: cc
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: acr
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: acs
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: sice
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: sa
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: orog
      REAL, DIMENSION(nlonpe,nlatpe,ntype),     INTENT(INOUT) :: gsf
      REAL, DIMENSION(nlonpe,nlatpe,ntype),     INTENT(INOUT) :: lai_ft
      REAL, DIMENSION(nlonpe,nlatpe,npft),      INTENT(INOUT) :: canht
      REAL, DIMENSION(nlonpe,nlatpe,ntype),     INTENT(INOUT) :: z0
      REAL, DIMENSION(nlonpe,nlatpe,npft),      INTENT(INOUT) :: gc
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),   INTENT(INOUT) :: o3um
      REAL, DIMENSION(nlonpe,nlatpe,ntype),     INTENT(INOUT) :: t0tile
      REAL, DIMENSION(nlonpe,nlatpe),         INTENT(INOUT) :: sicetemp
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(INOUT) :: soilmc
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(OUT)   :: p
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(OUT)   :: lnp
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(OUT)   :: acp
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(OUT)   :: adp
      REAL, DIMENSION(nlonpe,nlatpe,ntype),     INTENT(OUT)   :: ra
      REAL, DIMENSION(nlonpe,nlatpe,nc),        INTENT(OUT)   :: rq
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(OUT)   :: so4_vd

#include "soil_thick.h"

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: n
      INTEGER, DIMENSION(ntype-1) :: land_tiles =                       &
     &  (/ 1, 2, 3, 4, 5, 6, 8, 9 /)
      LOGICAL :: ldebug = .false.         ! Turns on/off debug o/p.
      REAL :: seafrac
      REAL :: lftotal
#include "c_mdi.h"

! Set sea ice, snow amount, lai and height of orography equal to zero
! where negative (mdi)
      WHERE(sice < 0.0) sice = 0.0
      WHERE(sa < 0.0) sa = 0.0
      WHERE(lai_ft < 0.0) lai_ft = 0.0
      WHERE(canht < 0.0) canht = 0.0
      WHERE(orog < 0.0) orog = 0.0
      WHERE(gc < 0.0) gc = 0.0
      WHERE(gsf < 0.0) gsf = 0.0
      WHERE(land_fraction < 0.0) land_fraction = 0.0
      WHERE(soilmc < 0.01) soilmc = 0.01 ! Don't want soilmc ~ 0

!     gsfcopy = gsf

! Set convective precipitation rates to zero if no cloud present
! Ensure cloud amount itself is zero or larger
      WHERE(MAXVAL(cc,DIM=3) <= 0.0)
        acr = 0.0
        acs = 0.0
      END WHERE
      cc = MAX(cc,0.0)

! Total precipitation = rain+snow
      acp = acr + acs
      adp = adr + ads

! Compose global pressure field, P
      p(:,:,0) = p0
      p(:,:,1:nmetlev) = p_th
      WHERE (p /= rmdi)
        lnp = LOG(p)
      ELSEWHERE
        lnp = rmdi
      END WHERE

! Convert ozone from mmr to vmr
      o3um = o3um * mair / mo3

! Tile fractions from STASH (in array gsf) are only defined over land.
! Set water fractions over the sea. Add sea ice and sea fractions to gsf
! Adjust tile fractions so they add up to 1. If gsf has the water tile
! fraction = 1.0, set the land fraction to 0.0

      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (gsf(i,k,nwater) == 1.0) land_fraction(i,k) = 0.0
          seafrac = 1.0 - land_fraction(i,k)
          IF (land_fraction(i,k) < 1.0 .AND. gsf(i,k,nwater) < 1.0) THEN
            gsf(i,k,nwater) = seafrac
            IF (land_fraction(i,k) > 0.0) THEN
              lftotal = 0.0
!CDIR UNROLL=8
              DO j = 1, ntype-1
                n = land_tiles(j)
                lftotal = lftotal + gsf(i,k,n)
              END DO
!CDIR UNROLL=8
              DO j = 1, ntype-1
                n = land_tiles(j)
                gsf(i,k,n) = gsf(i,k,n) * land_fraction(i,k) / lftotal
              END DO
            END IF
          END IF
          IF (sice(i,k) > 0.0) THEN
            gsf(i,k,nwater) = (1.0 - sice(i,k)) * seafrac
            gsf(i,k,ntype) = gsf(i,k,ntype) + sice(i,k) * seafrac
          END IF
        END DO
      END DO
!
! Set all tile temperatures to t0 where undefined.
!
      DO n = 1, ntype
        DO k = 1, rowspe
          DO i = 1, nlonpe
            IF (t0tile(i,k,n) < 100.0) t0tile(i,k,n) = t0(i,k)
          END DO
        END DO
      END DO
!
! Set up tile temperatures where a mixture of sea and sea ice is
! present. Set sea to freezing temperature (tfs) and ice to sea ice
! temperature.
!
      DO k = 1, rowspe
        DO i = 1, nlonpe
          IF (sice(i,k) > 0.0) THEN
            t0tile(i,k,ntype) = sicetemp(i,k)
            t0tile(i,k,nwater) = tfs
          END IF
        END DO
      END DO

! Calculate resistances for dry deposition
! DEPENDS ON: aerod
      CALL AEROD(t0,p0,hf,bl,tauu,tauv,canht,gsf,z0,ra,rq,so4_vd)

! Convert soil moisture values from kg m-2 to a volume ratio
! Note that kg (H2O) m-2 == depth in mm. Soil layer depth
! dzsoil is in metres - multiply by 1000 to convert to mm.

      DO k = 1, rowspe
        DO i = 1, nlonpe
          soilmc(i,k) = soilmc(i,k) / (dzsoil(1) * 1000.0)
        END DO
      END DO

! Print out some results for checking:
      IF (mype == 0 .AND. ldebug) THEN
        i=70
        j=13
        k=1
        WRITE(6,*) ' *** PROCMET: PE 0, Results for I,J,K=',i,j,k
        WRITE(6,*)  ' P0: ',p0(i,j)
        WRITE(6,*)  ' T0: ',t0(i,j)
        WRITE(6,*)  ' CC: ',cc(i,j,k)
        WRITE(6,*)  ' ACP: ',acp(i,j)
        WRITE(6,*)  ' ADP: ',adp(i,j)
        WRITE(6,*)  ' HF: ',hf(i,j)
        WRITE(6,*)  ' TAUU: ',tauu(i,j)
        WRITE(6,*)  ' TAUV: ',tauv(i,j)
        WRITE(6,*)  ' SICE: ',sice(i,j)
        WRITE(6,*)  ' SA: ',sa(i,j)
        WRITE(6,*)  ' '
      END IF

      END SUBROUTINE PROCESS_MET
#endif
