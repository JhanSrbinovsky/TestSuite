#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE STRATCALC2(xx,em,pos,ipos,nfill,trop,o3_ls,astep,lnp)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : CALCULATE O3 & HNO3 influx from stratosphere
!-
!-   Inputs  : POS,TROP,LNP,O3_LS
!-
!-   Outputs : XX
!
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.2    22/02/00  Created. W.J. Collins
!  5.5    13/02/04  Extensive modifications for vectorisation on SX6.
!                   K. Ketelsen
!  6.1    21/10/04  Refromatted code. M.G. Sanderson
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      USE HEIGHT_MOD                 !kk
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                          INTENT(IN) :: nfill
      INTEGER, DIMENSION(3,nclprc),     INTENT(IN) :: ipos

      REAL,                                     INTENT(IN) :: astep
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: trop
      REAL, DIMENSION(4,nclprc),                INTENT(IN) :: pos
      REAL, DIMENSION(nlnpe,nlpe,nmetlev),      INTENT(IN) :: o3_ls
      REAL, DIMENSION(nc,nclprc),               INTENT(INOUT) :: xx
      REAL, DIMENSION(nc+2,nclprc),             INTENT(INOUT) :: em

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l

      REAL :: o3strat
      REAL :: noystrat
      REAL :: deltao3
      REAL :: deltanoy
      REAL :: o3up
      REAL :: o3low
      REAL, PARAMETER :: noy_frac=1./1000.! Ratio of kg[N]/kg[O3]
                                          ! from Murphy and Fahey 1994
      REAL, PARAMETER :: stau=20.*daysec ! 20 day relaxation time
      REAL, DIMENSION(4,nfill) :: postop
      REAL, DIMENSION(4,nfill) :: posmid
      REAL, DIMENSION(4,nfill) :: posbot
      REAL, DIMENSION(4,nfill) :: posup
      REAL, DIMENSION(4,nfill) :: poslow
      REAL :: pup
      REAL :: plow
!kk
      LOGICAL,DIMENSION(nfill) :: flag
      INTEGER,DIMENSION(nfill) :: la
      REAL,DIMENSION(nfill)    :: val,valtop,valmid,valbot
      REAL,DIMENSION(nfill)    :: valup,vallow

      em(i_o3,:)   = 0.0
      em(i_hno3,:) = 0.0
      em(i_trop,:) = 0.0
      em(i_strat,:)= 0.0

! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,pos,lnp,val)
      DO j=1,nfill
! DEPENDS ON: getmetpoint
         flag(j) = (val(j) < GETMETPOINT(pos(:,j),trop,.TRUE.))
      END DO
      CALL HEIGHT_ETA_RHO(pos(:,1:nfill),la)

      DO i=1,4
        DO j=1,nfill
          IF (i == 3) THEN
            posbot(3,j) = 1.0
            posmid(3,j) = 1.0
            postop(3,j) = 1.0
            posup(3,j)  = 1.0
            poslow(3,j) = 1.0
          ELSE
            postop(i,j) = pos(i,j)
            posmid(i,j) = pos(i,j)
            posbot(i,j) = pos(i,j)
            posup(i,j)  = pos(i,j)
            poslow(i,j) = pos(i,j)
          END IF
        END DO
      END DO

! O3 columns are specified on Eta_theta levels.
! Taking differences between levels gives O3 mixing ratios on half-level
! Assume these half-levels are Eta_rho for simplicity.
      DO j=1,nfill
        IF (flag(j)) THEN
! above TPause
          l = la(j)
          i=ipos(1,j)-lndat+1
          k=ipos(2,j)-ltdat+1
! Eta_rho(L) lies between Eta_theta(L-1) and Eta_theta(L)
          posbot(3,j) = eta_theta(l-1)
          posmid(3,j) = eta_theta(l)
          poslow(3,j) = eta_rho(l)
          IF (l<nmetlev) THEN
            postop(3,j) = eta_theta(l+1)
            posup(3,j) = eta_rho(l+1)
          END IF
        END IF
      END DO

! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,posbot,lnp,valbot)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,posmid,lnp,valmid)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,postop,lnp,valtop)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,posup, lnp,valup)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,poslow,lnp,vallow)

! Calculate average ozone in upper and lower layers
! need to convert from molecules O3/m2 to moles O3/mole air
      DO j=1,nfill
        IF (flag(j)) THEN
! above TPause
          l = la(j)
          i = ipos(1,j)-lndat+1
          k = ipos(2,j)-ltdat+1

          IF (l < nmetlev) THEN
            postop(3,j) = eta_theta(L+1)
            o3up=(o3_ls(i,k,l)-o3_ls(i,k,l+1))*g*mair/                  &
     &       ((valmid(j)-valtop(j))*na)
            pup = valup(j)
          ELSE
            o3up = o3_ls(i,k,l)*g*mair/(valmid(j)*na)
            pup = valmid(j)/2.0
          END IF
          o3low=(o3_ls(i,k,l-1)-o3_ls(i,k,l))*g*mair/                   &
     &     ((valbot(j)-valmid(j))*na)
          plow = vallow(j)

! Interpolate between upper and lower layers to find ozone at
! parcel height
          o3strat=o3low+(o3up-o3low)*                                   &
     &      (val(j)-plow)/(pup-plow)
          noystrat=o3strat*noy_frac*(mo3/mnit)
          deltao3=o3strat-xx(i_o3,j)
          deltanoy=noystrat-xx(i_hno3,j)
          em(i_o3,j)=deltao3/stau
          em(i_hno3,j)=deltanoy/stau
          xx(i_strat,j)=1.0
          em(i_trop,j)=-xx(i_trop,j)/(2*daysec) ! 2 day lifetime
        ELSE                           ! in troposphere
          xx(i_trop,j)=1.0
          em(i_strat,j)=-xx(i_strat,j)/(2*daysec) ! 2 day lifetime
        END IF
      END DO

      END SUBROUTINE STRATCALC2
#endif
