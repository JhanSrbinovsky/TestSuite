#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing all routines relating to 2D photolysis
!  used in UKCA sub-model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_PHOT2D

      USE ASAD_MOD,          ONLY: jpspj
      USE UKCA_CHEM1_DAT
      IMPLICIT NONE
      SAVE
      PRIVATE

#include "parvars.h"
#include "c_pi.h"

      INTEGER              :: myjppj  ! set = jppj below

      INTEGER,PARAMETER,PUBLIC        :: nolat=19   ! no of 2D lats
      INTEGER,PARAMETER,PUBLIC        :: nolev=17   ! no of 2D levs
      INTEGER,PARAMETER,PUBLIC        :: nlphot=51  ! no of 2D photol levs
      INTEGER,PARAMETER,PUBLIC        :: ntphot=3   ! no of times of day
      INTEGER,PARAMETER               :: jpphio=84  ! unit number
      INTEGER,PARAMETER               :: jpphin=58  ! unit number

!     Number of time intervals in data (74 corresponds to
!     a 5 day interval, with one extra set).

      INTEGER,PARAMETER               :: n_data_times=74  
      REAL,   PARAMETER               :: data_interval=5.0 ! interval in days

      LOGICAL, PARAMETER              :: L_optimised=.true.  ! T for optimised read        

!     Photolysis rates from 2-D model which have been interpolated onto 3-D

      REAL,ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC          :: pjin
 
      PUBLIC UKCA_PHOTIN, UKCA_CURVE, UKCA_INPR2D

      CONTAINS

!-----------------------------------------------------------------------

      SUBROUTINE UKCA_PHOTIN(i_day_number,row_lengthda,tot_p_rows,     &
                      p_levelsda,p_fieldda,first_row, glon, glat,      & 
                      sinlat,pl,jppj) 

! Purpose: Subroutine to read in 2D photolysis rates, and interpolate on
!          3-d latitudes, heights. reconstruct daily curves (every 5 day
!          plus code to account for hour angle (longitude)
!
!          3 values are stored for each level and each latitude
!          symmetrical distribution plus zero values at dawn and
!          dusk gives 7 points altogether.
!          Based on PHOTIN.F from Cambridge TOMCAT model.

! dataset: 51(levels)x19(latitudes)x3(points)x74(every 5days)
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: jppj
      INTEGER, INTENT(IN) :: i_day_number
      INTEGER, INTENT(IN) :: row_lengthda
      INTEGER, INTENT(IN) :: tot_p_rows
      INTEGER, INTENT(IN) :: p_levelsda
      INTEGER, INTENT(IN) :: p_fieldda
      INTEGER, INTENT(IN) :: first_row 
      INTEGER, INTENT(IN) :: glon 
      INTEGER, INTENT(IN) :: glat 

      REAL, INTENT(IN), DIMENSION(:)     :: SINLAT
      REAL, INTENT(IN), DIMENSION(:,:,:) :: PL

!     Local variables

      INTEGER :: i,j 
      INTEGER :: idofy                         ! Day number
      INTEGER :: ipos
      LOGICAL,SAVE  :: first_call=.true.

      REAL :: fpos
      REAL :: pr2d(nolev)                      ! 2D model level pressu
      REAL :: pr2dj(nlphot)                    ! 2D photolysis level p
      REAL :: zmean_pl(tot_p_rows,p_levelsda)  ! zonal mean of 3D press
      REAL,dimension(:,:,:,:),allocatable    :: pjin2d


      COMMON /STREQ/ PR2D

      IF(first_call)  THEN
        myjppj = jppj
        CALL phot2d_allocate_memory (tot_P_ROWS,P_LEVELSDA,ntphot)
        first_call = .false.
      END IF

      ALLOCATE(pjin2d(nolat,nlphot,ntphot,myjppj))

!     Find nearest day in photolysis dataset - day 1 is 31. December

      idofy = i_day_number
      fpos  = idofy/data_interval + 1.0
      ipos  = nint(fpos)
      IF ((fpos-ipos*1.0) < 0.0) ipos = ipos-1

!     Read in photolysis data

      IF (L_optimised) THEN
        CALL READ2D_OPT(ipos,fpos,pjin2d)
      ELSE
        CALL READ2D_ORIG(ipos,fpos,pjin2d)
      ENDIF

!     Set up 2-D pressure arrays

      CALL UKCA_INPR2D(pr2d,pr2dj)

!     Interpolate 2D photolysis rates onto 3-D levels and 
!     latitudes. Longitude comes later.

      CALL UKCA_CALC_ZMEAN_PRESS(row_lengthda, tot_p_rows,           & 
                                 p_levelsda, glon, pl, zmean_pl) 
         
      CALL UKCA_INTERPJ(pjin2d,pr2dj,zmean_pl,row_lengthda,          & 
                        tot_p_rows, p_levelsda,                      & 
                        p_fieldda,sinlat,Pi_Over_180 ) 
                   

      DEALLOCATE(pjin2d)

      RETURN
      END SUBROUTINE UKCA_PHOTIN

!---------------------------------------------------------------------

        SUBROUTINE UKCA_CURVE(                                         &
                       pjinda,tloc,dayl,p_field,p_levels,              &
                       tot_p_rows,row_length,wks)
!
! Purpose: Subroutine to interpolate tropospheric photolysis rates
!          in time. Based on curve.F from Cambridge TOMCAT model.
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: p_field                   ! no of points
        INTEGER, INTENT(IN) :: p_levels                  ! no of vert
        INTEGER, INTENT(IN) :: tot_p_rows                ! no of rows
        INTEGER, INTENT(IN) :: row_length                ! no of cols

        REAL, DIMENSION(:), INTENT(IN) :: dayl           ! day length
        REAL, DIMENSION(:), INTENT(IN) :: tloc           ! local time
        REAL, DIMENSION(:,:,:), INTENT(IN) :: pjinda     ! 2D photolys

        REAL, DIMENSION(:,:), INTENT(OUT) :: wks         ! interpolated

! Local variables

        INTEGER :: i                                       ! loop variab
        INTEGER :: j                                       ! loop variables
        INTEGER :: k                                       ! loop variables
        INTEGER :: jr                                      ! loop variables

        REAL, PARAMETER :: tfrac1 = 0.04691008             ! determines
        REAL, PARAMETER :: tfrac2 = 0.23076534             ! 2D photolys

        REAL :: dawn                ! time of dawn
        REAL :: dusk                ! time of dusk
        REAL :: timel               ! local time
        REAL :: slope               ! slope used in linear interpolation
        REAL :: const               ! intercept used in linear interpola
        REAL :: fgmt(7)             ! times at which photol rates are va

! Initialise wks

        wks = 0.0

! Calculate rates using a simple linear interpolation.

        DO j = 1,tot_p_rows
          DO i = 1,row_length
            k = i+(j-1)*row_length

! Non-Polar night

            IF (dayl(k) > 0.0) THEN
              dawn = 12.00 - (dayl(k)/2.0)
              dusk = 12.00 + (dayl(k)/2.0)

              fgmt(1) = dawn
              fgmt(2) = dawn + tfrac1*dayl(k)
              fgmt(3) = dawn + tfrac2*dayl(k)
              fgmt(4) = 12.00
              fgmt(5) = dawn + (1.0-tfrac2)*dayl(k)
              fgmt(6) = dawn + (1.0-tfrac1)*dayl(k)
              fgmt(7) = dusk

              timel = tloc(k)

              IF (timel > 24.0) timel = timel-24.0

              timel = min(timel,24.0)
              timel = max(timel, 0.0)

! Local Night-time

              IF (timel < dawn .OR. timel > dusk) THEN

                wks(k,1:MYjppj) = 0.0

! For the time between dawn and PJIN(1) or PJIN(5) and dusk

              ELSE IF ((timel >= dawn   .AND. timel < fgmt(2))         &
                   .OR.(timel > fgmt(6) .AND. timel <= dusk)) THEN

                IF (timel > fgmt(6)) timel = 24.00 - timel

                DO jr = 1, myjppj
                  slope = pjinda(j,1,jr)/(fgmt(2) - fgmt(1))
                  const = - slope* fgmt(1)
                  wks(k,jr)= slope*timel + const
                ENDDO

! For the time between PJIN(1) and PJIN(2) or PJIN(4) and PJIN(5)

              ELSE IF ((timel >= fgmt(2) .AND. timel < fgmt(3))        &
                   .OR.(timel >  fgmt(5) .AND. timel <= fgmt(6))) THEN

                IF (timel > fgmt(5)) timel = 24.00 - timel

                DO jr = 1, myjppj
                  slope = (pjinda(j,2,jr)-pjinda(j,1,jr))/             &
                       (fgmt(3) - fgmt(2))
                  const = pjinda(j,1,jr)- slope* fgmt(2)
                  wks(k,jr)= slope*timel + const
                ENDDO

! For the time between PJIN(2), PJIN(3) and PJIN(4)

              ELSE IF (timel >= fgmt(3) .AND. timel <= fgmt(5)) THEN

                IF (timel > fgmt(4)) timel = 24.00 - timel

                DO jr = 1, myjppj
                  slope = (pjinda(j,3,jr)-pjinda(j,2,jr))/             &
                      (fgmt(4) - fgmt(3))
                  const = pjinda(j,2,jr)- slope* fgmt(3)
                  wks(k,jr)= slope*timel + const
                ENDDO

              ENDIF    ! end of IF (timel < dawn .OR. timel > dusk)

! End of the condition on the non polar night

            ELSE

              wks(k,1:myjppj) = 0.0

            ENDIF     ! end of IF (dayl(k) > 0.0)

          ENDDO       ! end of looping over row length
        ENDDO         ! end of looping over rows

        RETURN
        END SUBROUTINE UKCA_CURVE

! ----------------------------------------------------------------------

        SUBROUTINE READ2D_OPT(ipos,fpos,pjin2d)

! ----------------------------------------------------------------------
! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

      INTEGER, INTENT(INOUT) :: ipos     ! integer position
      REAL,    INTENT(IN)    :: fpos     ! real position
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: pjin2d   ! Photol rate

! Local variables
      INTEGER :: ifinx(myjppj)
      INTEGER :: j                             ! Loop variable
      INTEGER :: jr                            ! Loop variable
      INTEGER :: k                             ! Loop variable
      INTEGER :: kk                            ! Loop variable
      INTEGER :: in                            ! Index
      INTEGER :: ierror                        ! Error flag
      INTEGER :: info                          ! Tag for communication

      REAL :: delpos
! Photol rates at all times
      REAL, allocatable, SAVE :: pjin2da(:,:,:,:,:)
      REAL :: pr(myjppj,3)
      REAL :: pfrac(1,3)

      LOGICAL :: L_exist                     ! Logical to check exist
      LOGICAL :: L_first=.true.              ! Logical for firstcall 

      CHARACTER(len=35) :: file2             ! Dir for photol files
      CHARACTER(len=60) :: fnme(myjppj)      ! inc length of string 
      CHARACTER(len=10) :: csp(myjppj,jpspj)
      CHARACTER(len=10) :: cmnt(myjppj)      ! Comment line
      CHARACTER(len=72) :: cmessage          ! Error message

! 1.  Determine filenames containing specified photolysis rates


      pjin2d = 0.0
! Reset ipos if nearest day=365 or zero
      delpos = fpos - ipos*1.0
      IF (ipos == 74) ipos=1

! 1.1  Read filenames from photolysis ratefile on first call

      IF (L_first) THEN
        ALLOCATE(pjin2da(myjppj,n_data_times,nlphot,ntphot,nolat))

!       use module to get cmnt

        IF (SIZE(ratj_defs) /= jppj) THEN
          cmessage='size of ratj_defs is not equal to jppj'
! DEPENDS ON: EREPORT
          CALL EREPORT('UKCA_PHOT2D.UKCA_READ2D',1,cmessage)
        ENDIF
        DO k=1,jppj
          cmnt(k)=ratj_defs(k)%fname
        ENDDO

!       1.2  Add '.bin' extension

        file2='/data/cr/cce/hadcj/tropdata/photol/'
        DO jr = 1, MYjppj
          in = INDEX(cmnt(jr),' ') - 1
          IF ( in .lt. 0 ) in = 10
          fnme(jr) = file2//cmnt(jr)(1:in)//'.bin'
          IF (mype == 0) PRINT*,'fnme =', fnme(jr), cmnt(jr)
        ENDDO

! 2.  Read specified photolysis rates

        IF (mype == 0) THEN
          DO jr = 1, MYjppj
            INQUIRE (FILE=FNME(JR), EXIST=L_exist)
            IF (.NOT. L_exist) THEN
              cmessage = 'File does not exist'
! DEPENDS ON: ereport
              CALL EREPORT('UKCA_PHOT2D.READ2D_OPT',jr,cmessage)
            ENDIF

          OPEN(jpphin, FILE=fnme(jr), FORM='UNFORMATTED', IOSTAT=ierror)
            IF (ierror /= 0) THEN
              cmessage = ' Error opening file '
! DEPENDS ON: ereport
              CALL EREPORT('UKCA_PHOT2D.READ2D_OPT',jr,cmessage)
            ENDIF

            READ(jpphin) pjin2da(jr,:,:,:,:)
            CLOSE(jpphin)

          ENDDO   ! jr
        ENDIF  ! mype
        l_first=.false.
      ENDIF  ! l_first

! Interpolate in time using the saved values
        IF (mype == 0) THEN
          DO jr = 1, MYjppj
            DO k = 1,nolat
              DO kk = 1,ntphot
                DO j = 1,nlphot
                  pjin2d(k,j,kk,jr) =                                  &
     &                  (pjin2da(jr,ipos+1,j,kk,k)-                    &
     &                   pjin2da(jr,ipos,j,kk,k))*delpos +             &
     &                   pjin2da(jr,ipos,j,kk,k)
              ENDDO
            ENDDO
          ENDDO
        ENDDO  ! jr
      ENDIF   ! mype

      CALL GC_RBCAST (1,nlphot*ntphot*nolat*MYjppj,0,nproc             &
                       ,info,pjin2d)
      CALL GC_SSYNC(nproc,info)

      RETURN
      END SUBROUTINE READ2D_OPT

! ---------------------------------------------------------------------

      SUBROUTINE READ2D_ORIG(ipos,fpos,pjin2d)

! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------

      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

      REAL, INTENT(IN) :: fpos

      INTEGER, INTENT(INOUT) :: ipos

      REAL, INTENT(OUT) :: pjin2d(nolat,nlphot,ntphot,myjppj) ! Photol rate

! Local variables

      INTEGER :: ifinx(myjppj)
      INTEGER :: ij                            ! Loop variable
      INTEGER :: j                             ! Loop variable
      INTEGER :: jr                            ! Loop variable
      INTEGER :: k                             ! Loop variable
      INTEGER :: kk                            ! Loop variable
      INTEGER :: in                            ! Index
      INTEGER :: ierror                        ! Error flag
      INTEGER :: info                          ! Tag for communication

      REAL :: delpos
      REAL :: pjin2d1(nlphot,ntphot,nolat)! Photol rates straddling time
      REAL :: pjin2d2(nlphot,ntphot,nolat)! Photol rates straddling time
      REAL :: pr(myjppj,3)
      REAL :: pfrac(1,3)

      LOGICAL :: L_exist                       ! Logical used to check ex

      CHARACTER(len=35) :: file2             ! Dir for photol files
      CHARACTER(len=60) :: fnme(myjppj)      ! inc length of string 
      CHARACTER(len=10) :: csp(myjppj,jpspj)
      CHARACTER(len=10) :: cmnt(myjppj)      ! Comment line
      CHARACTER(len=72) :: cmessage          ! Error message

! 1.  Determine filenames from module info

      IF (SIZE(ratj_defs) /= jppj) THEN
        cmessage='size of ratj_defs is not equal to jppj'
! DEPENDS ON: EREPORT
        CALL EREPORT('UKCA_PHOT2D.UKCA_READ2D',1,cmessage)
      ENDIF 

      DO k=1,jppj
        cmnt(k)=ratj_defs(k)%fname
      ENDDO 


! 1.2  Add '.d' extension

!        file2='/home/cr/cce/hadfo/um/tropdata/photol/'
      file2='/data/cr/cce/hadcj/tropdata/photol/'
      DO jr = 1, myjppj
        in = index(cmnt(jr),' ') - 1
        IF ( in .lt. 0 ) in = 10
        fnme(jr) = file2//cmnt(jr)(1:in)//'.dat'
        IF (mype == 0) PRINT*,'fnme =', fnme(jr), cmnt(jr)
      ENDDO

! Reset ipos if nearest day=365 or zero

      delpos = fpos - ipos*1.0
      IF (ipos == 74) ipos=1

! 2.  Read specified photolysis rates

      pjin2d = 0.0
      IF (mype == 0) THEN
        DO jr = 1, myjppj
          INQUIRE (FILE=FNME(JR), EXIST=L_exist)
          IF (.NOT. L_exist) THEN
            cmessage = 'File does not exist'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_PHOT2D.READ2D_ORIG',jr,cmessage)
          ENDIF

          OPEN(jpphin, FILE=fnme(jr), FORM='FORMATTED', IOSTAT=ierror)
          IF (ierror /= 0) THEN
            cmessage = 'Error reading file '
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_PHOT2D.READ2D_ORIG',jr,cmessage)
          ENDIF

          DO ij = 1,ipos
            read(jpphin,*) (((pjin2d1(j,kk,k),j=1,51),kk=1,3)          &
                           ,k=nolat,1,-1)
          ENDDO

          READ(jpphin,*) (((pjin2d2(j,kk,k),j=1,51),kk=1,3)            &
                        ,k=nolat,1,-1)

          DO k = 1,nolat
            DO kk = 1,ntphot
              DO j = 1,nlphot
                pjin2d(k,j,kk,jr) =                                    &
                         (pjin2d2(j,kk,k)-pjin2d1(j,kk,k))*            &
                          delpos + pjin2d1(j,kk,k)              
              ENDDO
            ENDDO
          ENDDO

          CLOSE(jpphin)

        ENDDO
      ENDIF      ! End of IF mype statement

      CALL GC_RBCAST (1,nlphot*ntphot*nolat*myjppj,0,nproc             &
                       ,info,pjin2d)
      CALL GC_SSYNC(nproc,info)

      RETURN
      END SUBROUTINE READ2D_ORIG

!------------------------------------------------------------

      SUBROUTINE UKCA_INPR2D(pr2d,pr2dj)
!
! Purpose: Subroutine to calculate the pressure levels for the
!          2D photolysis rates. Original version taken from the
!          Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        IMPLICIT NONE

!       Local variables

        INTEGER, PARAMETER :: maxlev = 30

        INTEGER :: j                        ! Loop variable
        INTEGER :: ij                       ! Loop variable

        REAL, PARAMETER :: fac  = 1.284025417  ! Used to ensure that pre
        REAL, PARAMETER :: psur = 1.0e5        ! Surface pressure in Pas
        REAL, PARAMETER :: ares = 3.0          ! Factor related to verti
        REAL, PARAMETER :: eps  = 1.0e-10      ! Factor used to calculat

        REAL :: ee                       ! Temporary store
        REAL :: fj                       ! Factor related to vertical re
        REAL :: za                       ! Factor related to vertical re
        REAL :: pr2d(nolev)              ! Pressures of 2D model levels
        REAL :: pr2dj(nlphot)            ! Pressures of 2D photolysis le
        REAL :: pp(nlphot+1)
        REAL :: pres(maxlev)

        DO j = 1, nolev
          ee = exp((j-1)/2.0)
          pr2d(j) = psur/(ee * fac)
        ENDDO

!       2D pressure levels - normal 2-D order - up to level 30
!       for photolysis

        DO j = 1,maxlev-1
          ee = exp((j-1) / 2.0)
          pres(j) = psur / (ee * fac)
        ENDDO
        pres(maxlev) = pres(maxlev-1)

        fj    = 2.0/ares
        pp(1) = (1.0-fj)*alog(psur)+fj*alog(pres(1))
        pp(1) = exp(pp(1))

        DO ij = 2,nlphot+1
          za     = ij/(2.0*ares)
          fj     = 2.0*za+0.5+eps
          j      = int(fj)
          fj     = fj-j-eps
          j      = j+1
          pp(ij) = (1.0-fj)*alog(pres(j-1))+ fj*alog(pres(j))
          pp(ij) = exp(pp(ij))
        ENDDO

        pr2dj(1)=(psur+pp(1))*0.5

        DO ij = 2,nlphot
          pr2dj(ij) = (pp(ij)+pp(ij-1))*0.5
        ENDDO

        RETURN
        END SUBROUTINE UKCA_INPR2D

!----------------------------------------------------------------

        SUBROUTINE UKCA_INTERPJ(pjin2d,pr2dj,zm_pl,lon,lat,lev,         & 
                                p_field,sinlat,degrad)
!
! Purpose: Subroutine to interpolate photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
        IMPLICIT NONE
        
        INTEGER :: lon                               ! No of longitudes
        INTEGER :: lat                               ! No of latitudes
        INTEGER :: lev                               ! No of levels
        INTEGER :: p_field                           ! No of spatial poi

        REAL, DIMENSION(:,:,:,:), INTENT(IN) :: pjin2d  ! 2D photo
        REAL, DIMENSION(:,:), INTENT(IN)     :: zm_pl   ! zmean 3D press 
        REAL, DIMENSION(:), INTENT(IN)       :: pr2dj   ! 2D photo
        REAL, DIMENSION(:), INTENT(IN)       :: sinlat  ! Sine (3D
        REAL, INTENT(IN)                     :: degrad  ! To conve
        
!       Local variables

        INTEGER :: i                     ! Loop variable
        INTEGER :: j                     ! Loop variable
        INTEGER :: jr                    ! Loop variable
        INTEGER :: k                     ! Loop variable
        INTEGER :: kk                    ! Loop variable
        INTEGER :: l                     ! Loop variable

        REAL, PARAMETER :: Npole = 90.0

        REAL :: p2d(nlphot)
        REAL :: lat2d(nolat)          ! 2D model latitudes
        REAL :: wks1(lat,nlphot)      ! Working array
        REAL :: wks2(nolat)           ! Working array
        REAL :: wks3(nlphot)          ! Working array
        REAL :: wks4(lat,lev)         ! Working array
        REAL :: lati
        REAL :: press
        REAL :: UKCA_FLUPJ            ! Function used for interpolation
        REAL :: delphi

!       Set up 2D latitudes. lat=1 pole nord
!       LAT2D() is the latitude in the centre of the 2D box in radians

        delphi = Npole/nolat
        DO i = 1,nolat
          lat2d(i)=sin((90.0 - (2*i-1)*delphi)*degrad)
        ENDDO

        DO jr = 1,MYjppj              ! Loop over photolysis reactions

!         Interpolate linearly in Sin(lat) (KK is the point through the

          DO kk = 1,ntphot          ! Loop over times of day

            DO j = 1,nlphot
              DO i = 1,nolat
                wks2(i) = pjin2d(i,j,kk,jr)
              ENDDO

              DO k = 1,lat
                lati = sinlat((k-1)*lon+1)
                lati = MAX(lati, lat2d(nolat))
                lati = MIN(lati, lat2d(    1))
!DEPENDS ON: ukca_flupj
                wks1(k,j) = UKCA_FLUPJ(lati,lat2d,wks2,nolat)
                wks1(k,j) = MAX(wks1(k,j), 0.0)
              ENDDO
            ENDDO

!           Interpolate linearly in log(P)

            DO k = 1,lat
              DO j = 1,nlphot
                wks3(j) = wks1(k,j)
                p2d(j)  = LOG(pr2dj(j))
              ENDDO

              DO l = 1,lev
                press = LOG(zm_pl(k,l)) 
                press = MAX(press, p2d(nlphot))
                press = MIN(press, p2d(1))
!DEPENDS ON: ukca_flupj
                wks4(k,l) = UKCA_FLUPJ(press,p2d,wks3,nlphot)
                wks4(k,l) = MAX(wks4(k,l), 0.0)
              ENDDO
            ENDDO

            DO l = 1,lev
              DO k = 1,lat
                pjin(k,l,kk,jr) = wks4(k,l)
              ENDDO
            ENDDO


          ENDDO     ! End of loop over times of day
        ENDDO       ! End of loop over photolysis reactions

        RETURN

      END SUBROUTINE UKCA_INTERPJ

!-----------------------------------------------------------------

      SUBROUTINE PHOT2D_ALLOCATE_MEMORY(lat,lev,ntphot)
!
! Purpose: Subroutine to allocate array to hold
!          2D photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,INTENT(IN)         :: lat,lev,ntphot
   
!       lat is the number of latitude on current PE, NOT nolat!
!       lev is numer of model level, NOT nolev

        ALLOCATE(pjin(lat,lev,ntphot,MYjppj))

        RETURN
      END SUBROUTINE PHOT2D_ALLOCATE_MEMORY

!---------------------------------------------------------------------- 
         
      SUBROUTINE UKCA_CALC_ZMEAN_PRESS(lon, lat, lev,              & 
     &                                  glon,pl, zmean_pl) 
! 
! Purpose: Subroutine to calculate zonal mean pressure profile 
! 
!          Called from UKCA_PHOTIN. 
! 
! Current code owner: Fiona O'Connor 
! 
! 
! Code description: 
!   Language: FORTRAN 90 
!   This code is written to UMDP3 v6 programming standards. 
! 
! --------------------------------------------------------------------- 
      IMPLICIT NONE 
         
#include "parvars.h" 
         
      INTEGER, INTENT(IN) :: lon        ! No of longitudes 
      INTEGER, INTENT(IN) :: lat        ! No of latitudes 
      INTEGER, INTENT(IN) :: lev        ! No of levels 
      INTEGER, INTENT(IN) :: glon       ! No of global lons 
 
      REAL, INTENT(IN)    :: pl(lon,lat,lev)   ! Pressure 
 
      REAL, INTENT(OUT)   :: zmean_pl (lat,lev) 
               
!     Local variables 
         
      INTEGER :: j,k                 ! Loop variables 
      INTEGER :: istat               ! Output status from GCOM routine 
         
      REAL    :: fac                 ! Factor used to calculate zmean 
      REAL    :: sumpl_row(lat,lev)  ! Global sum along rows 
         
      zmean_pl = 0.0 

      DO j=1,lev 
        DO k=1,lat 
          CALL GCG_RVECSUMR (lon, lon, 1, 1, pl(:,k,j),                     & 
                             gc_proc_row_group, istat, sumpl_row(k,j)) 
        ENDDO 
      ENDDO 
         
      fac      = 1.0/real(glon) 
      zmean_pl = sumpl_row*fac 
         
      RETURN 
      END SUBROUTINE UKCA_CALC_ZMEAN_PRESS 
         
!---------------------------------------------------------------------- 

      END MODULE UKCA_PHOT2D
#endif
