#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DWCALC(dw,adp,acp,cc,pos,nfill,lnp)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Calculate wet deposition rates
!-
!-   Inputs  : ADP,ACP,CT,IPOS   ! dynamic & convective ppn, convective
!                                ! cloud top, eta level.
!-   Outputs : DW  (1/s)
!-   Controls: Scavenging coeffs.
!-
!
! Current Code Owner: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.1    06/02/96  Created.  C.E. Johnson
!  5.5    04/12/01  New dynamics version. W.J. Collins
!  6.0    12/09/03  Fix data initialisation problems encountered on
!                   some platforms. Introduce standard UM modification
!                   history. P.Dando
!  6.0    10/03/04  Vectorised code. K. Ketelsen
!  6.1    06/09/04  Calls ST_HEIGHT as function renamed. M.G. Sanderson
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE IN_STOCHEM_CHM
      USE HEIGHT_MOD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: nfill
      REAL, DIMENSION(4,nclprc),                INTENT(IN) :: pos
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: acp
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: adp
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),   INTENT(IN) :: cc
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
      REAL, DIMENSION(nc,nclprc),              INTENT(OUT) :: dw

      REAL, SAVE, DIMENSION(nc) :: dsc
      REAL, SAVE, DIMENSION(nc) :: csc
      REAL, DIMENSION(nc) :: rj
      REAL, DIMENSION(nfill) :: p1
      REAL, DIMENSION(nfill) :: p_top

      INTEGER :: k
      INTEGER :: kmax
      INTEGER :: j
      INTEGER, DIMENSION(nfill) :: i_w
      INTEGER, DIMENSION(nfill) :: j_w


      LOGICAL :: first_entry = .true.

      REAL            :: cprofile
      REAL, PARAMETER :: dts=3*hoursec ! 3 hour correlation period
      REAL, PARAMETER :: df=1.0 ! Fraction of cell affected by dyn ppn
      REAL, PARAMETER :: cf=0.3 ! Fraction of cell affected by conv ppn
      REAL, DIMENSION(nlev) :: dprofile=                                &
                                         ! Assumed dynamic rainfall prof
     &  (/1.0,1.0,1.0,1.0,0.78,0.78,0.56,0.56,0.33,0.33,0.11,0.11,      &
     &    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

!kk   Local variables defined by Klaus Ketelsen
      LOGICAL, DIMENSION(nfill) :: action
      INTEGER, DIMENSION(nfill) :: ka
      REAL, DIMENSION(4,nfill)  :: rpos

!      non-zero values of scavenging coefficient (1/cm)
!      Penner, Atherton & Graedel (1994).
!      species    large scale  convective
!      HCHO         2.0          4.0
!      HNO3         2.4          4.7
!      N2O5         1.0          2.0
!      H2O2         2.4          4.7
!      ROOH         2.0          4.0
!      NAER         5.0          1.0
!      SO2          0.8          1.5
!      SA,MSA       5.0          1.5
!      NAER         5.0          1.0
!      NH3, NH42SO4 as H2O2
!      BE7, BE10, PB210 as SA
!      CH3OH as ROOH


      IF (first_entry) THEN
! Initialise
         dsc = 0.0
         csc = 0.0
! Set components
         dsc(i_n2o5)    =  1.0
         csc(i_n2o5)    =  2.0
         dsc(i_hcho)    =  2.0
         csc(i_hcho)    =  4.0
         dsc(i_hno3)    =  2.4
         csc(i_hno3)    =  4.7
         dsc(i_h2o2)    =  2.4
         csc(i_h2o2)    =  4.7
         dsc(i_ch3ooh)  =  2.0
         csc(i_ch3ooh)  =  4.0
         dsc(i_so2)     =  0.8
         csc(i_so2)     =  1.5
         dsc(i_sa)      =  5.0
         csc(i_sa)      =  1.5
         dsc(i_c3h7ooh) =  2.0
         csc(i_c3h7ooh) =  4.0
         dsc(i_c2h5ooh) =  2.0
         csc(i_c2h5ooh) =  4.0
         dsc(i_c4h9ooh) =  2.0
         csc(i_c4h9ooh) =  4.0
      dsc(i_ch3oh)   = 2.0
      csc(i_ch3oh)   = 4.0
         dsc(i_ammsul)  =  2.4
         csc(i_ammsul)  =  4.7
         dsc(i_isopooh) =  2.0
         csc(i_isopooh) =  4.0
         dsc(i_mvkooh)  =  2.0
         csc(i_mvkooh)  =  4.0
         dsc(i_naer)    =  5.0
         csc(i_naer)    =  1.0
         dsc(i_msa)     =  5.0
         csc(i_msa)     =  1.5
         dsc(i_orgnit)  =  5.0
         csc(i_orgnit)  =  1.0
         dsc(i_nh3)     =  2.4
         csc(i_nh3)     =  4.7
         dsc(i_be7)     =  5.0
         csc(i_be7)     =  1.5
         dsc(i_be10)    =  5.0
         csc(i_be10)    =  1.5
         dsc(i_pb210)   =  5.0
         csc(i_pb210)   =  1.5

         first_entry = .false.
      ENDIF


! DEPENDS ON: st_height
      kmax = ST_HEIGHT(0.35,'Eta_rho') ! No Conv. depn. above this heigh

      DO j=1,nfill
! Indicies for met grids
! DEPENDS ON: imet
        i_w(j) = IMET(pos(1,j),.TRUE.)
! DEPENDS ON: jmet
        j_w(j) = JMET(pos(2,j),.TRUE.)
        dw(:,j) = 0.0
      END DO

      CALL HEIGHT_ETA_STOCHEM(pos,ka)
!kk      DO j=1,nfill
!kk         ka(j) = HEIGHT(pos(3,j),'Eta_stochem')
!kk      end do

      DO j=1,nfill
! Dynamic ppn.
        IF (adp(i_w(j),j_w(j)) > 1.0e-08 .AND. pos(3,j) < 0.15) THEN
! The factor of 10.0 converts ppn. from kg/(m^2.s) to cm/s.
          rj = dsc * (adp(i_w(j),j_w(j))/10.0) * dprofile(ka(j))/df
          dw(:,j) = rj
        END IF
      END DO

! Convective ppn: Ppn. profile 1.0 below eta=0.85,
! then declines linearly.
      DO j=1,nfill
        IF (acp(i_w(j),j_w(j)) > 1.0e-08) THEN
          action(j) = .true.
        ELSE
          action(j) = .false.
        END IF
      END DO

      DO j=1,nfill
        IF (action(j)) THEN
          DO k=nmetlev,1,-1
            IF(cc(i_w(j),j_w(j),k)>0.05 .AND. K<Kmax) EXIT
          END DO
          rpos(1,j)=longm_half(i_w(j)+lnbound-1)
          rpos(2,j)=latm_half(j_w(j)+lobound-1)
          rpos(3,j)=Eta_Theta(k)
        ELSE
          rpos(1,j)=longm_half(i_w(j)+lnbound-1)
          rpos(2,j)=latm_half(j_w(j)+lobound-1)
          rpos(3,j)=1.0
        END IF
      END DO

! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,rpos,lnp,p_top)
! DEPENDS ON: eta2p_v
      CALL ETA2P_V(nfill,pos,lnp,p1)

      DO j=1,nfill
        IF (action(j)) THEN
          IF (p1(j) > p_top(j)) THEN
            IF(p1(j) >= 85000.0) THEN
              cprofile = 1.0
            ELSE
              cprofile = 1.0 - ((85000.0-p1(j)) / (85000.0-p_top(j)))
            END IF
            rj = csc * acp(i_w(j),j_w(j))*cprofile / (cf*10.0)
            dw(:,j) = dw(:,j) - LOG(1.0-cf+cf*EXP(-dts*rj))/dts
          END IF
        END IF
      END DO

      END SUBROUTINE DWCALC
#endif
