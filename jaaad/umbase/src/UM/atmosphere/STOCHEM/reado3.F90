#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE reado3(o3_ls,month,z_top_of_model,                     &
     &                  first_constant_r_rho_level)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Read Li+Shine(+/- Piers trends) /Jeff Knights
!-                         O3 data
!-
!-   Inputs  : month,year
!-   Outputs : D
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.5    06/10/98  Created. D.S. Stevenson
!  5.5    05/10/01  Converted to ND eta grid and mods to reverse
!                   latitudes. C.E. Johnson
!  5.5    13/02/04  Trivial change to OPEN statement. K. Ketelsen
!  6.1    21/10/04  No change
!-
!VVV v5.2.2  5/X/01  Converted to ND eta grid
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: month
      INTEGER, INTENT(IN) :: first_constant_r_rho_level
      REAL,    INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlnpe,nlpe,nmetlev), INTENT(OUT) :: o3_ls

      INTEGER :: lun=30
      INTEGER :: i
      INTEGER :: ii
      INTEGER :: ii2
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: k
      INTEGER :: kk
      INTEGER, PARAMETER :: nlg_ugamp=144
      INTEGER, PARAMETER :: nlt_ugamp=73
      INTEGER, PARAMETER :: nlv_ugamp=47

      REAL :: fi
      REAL :: fj
      REAL :: fk
      REAL :: p
      REAL :: z
      REAL, DIMENSION(nlg_ugamp,nlt_ugamp,nlv_ugamp) :: ugamp
      REAL, DIMENSION(nlv_ugamp) :: p_ugamp
      REAL, DIMENSION(nlv_ugamp) :: eta_ugamp
      REAL, DIMENSION(nlg_ugamp), PARAMETER :: long_ugamp=              &
     & (/(i,i=0,nlg_ugamp-1)/)*(360./nlg_ugamp)
      REAL, DIMENSION(nlt_ugamp), PARAMETER :: lat_ugamp=               &
     & (/(i,i=0,nlt_ugamp-1)/)*(180./(nlt_ugamp-1))
      CHARACTER*72 :: cmessage
      CHARACTER(LEN=3), PARAMETER, DIMENSION(12) :: mth=                &
     & (/'jan','feb','mar','apr','may','jun',                           &
     &   'jul','aug','sep','oct','nov','dec'/)

      OPEN(LUN,FILE=trim(datdir)//mth(month)//'5ym',                    &
     &     STATUS='OLD')
! UGAMP vertical coordinate starts at the top. Need to reverse this.
      READ(lun,*)
      READ(lun,*)
      READ(lun,*) p_ugamp(nlv_ugamp:1:-1)
      READ(lun,*)
      READ(lun,*) ugamp(:,:,nlv_ugamp:1:-1)
      CLOSE(lun)
! Reverse latitudes for S => N grid. South to North anyway!
!      ugamp=ugamp(:,nlt_ugamp:1:-1,:)
!      lat_ugamp=lat0_ugamp(nlt_ugamp:1:-1)

! Find eta coordinates from P_UGAMP using standard atmosphere
      DO i=1,nlv_ugamp
        p=p_ugamp(i)*100.     ! mbar => Pa
! DEPENDS ON: rstd
        CALL RSTD(p,z)
        IF (z <= Z_Top_of_Model) THEN
! DEPENDS ON: rtoeta
          Eta_UGAMP(I)=RtoEta(z+Earth_Radius,0.0,z_top_of_model,        &
     &                        first_constant_r_rho_level)
        ELSE
          Eta_UGAMP(I)=z/Z_Top_of_Model
        ENDIF
      ENDDO
      WRITE(6,*) 'Eta_UGAMP:'
      WRITE(6,'(8F8.4)') Eta_UGAMP

      DO k=1,nmetlev
        kk=1
        DO
          IF (Eta_UGAMP(kk) > Eta_theta(k)) EXIT
          kk=kk+1
        ENDDO
        kk=kk-1
        fk=(ETA_theta(k)-Eta_UGAMP(kk))/(ETA_UGAMP(kk+1)-ETA_UGAMP(kk))
        DO j=1,nlpe
          jj=int((j+ltdat-1-.5)*REAL(nlt_ugamp-1)/REAL(mnlat))+1
          fj=(lat(j+ltdat-1-1)+dlat/2.-lat_ugamp(jj))/                  &
     &       (lat_ugamp(jj+1)-lat_ugamp(jj))
          DO i=1,nlnpe
            ii=int((i+lndat-1-.5)*REAL(nlg_ugamp)/REAL(nlong))+1
            ii2=MOD(ii,nlg_ugamp)+1
            fi=(long(i+lndat-1)+dlong/2.-long_ugamp(ii))/               &
     &         (360./nlg_ugamp)
            o3_ls(i,j,k)=                                               &
     &     fi *    fj *    fk *ugamp(ii2,jj+1,kk+1)+                    &
     & (1.-fi)*    fj *    fk *ugamp(ii ,jj+1,kk+1)+                    &
     &     fi *(1.-fj)*    fk *ugamp(ii2,jj  ,kk+1)+                    &
     & (1.-fi)*(1.-fj)*    fk *ugamp(ii ,jj  ,kk+1)+                    &
     &     fi *    fj *(1.-fk)*ugamp(ii2,jj+1,kk)+                      &
     & (1.-fi)*    fj *(1.-fk)*ugamp(ii ,jj+1,kk)+                      &
     &     fi *(1.-fj)*(1.-fk)*ugamp(ii2,jj  ,kk)+                      &
     & (1.-fi)*(1.-fj)*(1.-fk)*ugamp(ii ,jj  ,kk)
          IF(fi<0..or.fj<0..or.fk<0..or.fi>1..or.fj>1..or.fk>1.) THEN
            cmessage='FI etc out of bounds'
            WRITE(6,*) 'READO3: ',cmessage
            WRITE(6,*) 'FI,FJ,FK=',fi,fj,fk
! DEPENDS ON: ereport
            CALL EREPORT('READO3',1,cmessage)
          ENDIF
          END DO
        END DO
      END DO

      o3_ls=o3_ls*dob
      WRITE(6,*) 'O3_LS: (Dobsons)'
      WRITE(6,'(5e12.4)') o3_ls(1,1,:)/dob

      END SUBROUTINE READO3
#endif
