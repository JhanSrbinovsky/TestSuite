#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE JVALS(dj,dja,ipos,cloud,o3col,o3conc,land,snowfr,      &
     &  seaicefr,time,t0,t,p,nfill)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Evaluate photolysis rates.
!-
!-   Inputs  : CLOUD,OZONE,O3CONC,LAND,SNOWFR,SEAICEFR,P0,T
!-   Outputs : DJA
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.2    09/08/96  Created.  C.E. Johnson
!  4.5    27/04/98  Modifications to code structure.  W.J. Collins
!  5.3    29/11/01  Corrected LAT calculation. W.J. Collins
!  5.5    22/03/04  Extensive changes for vectorisation and load
!                   balancing across PEs. K. Ketelsen
!  6.2    06/04/05   Further changes for running on any PE
!                    configuration. R. Johanni

!VVV  V4.6  JVALS 18/1/01 Now returns DJ as well as DJA
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE EINTERP_MOD              !kk
      USE PHOT_MOD                 !mgs
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER,                      INTENT(IN) :: nfill
      INTEGER, DIMENSION(3,nclprc), INTENT(IN) :: ipos

      REAL, DIMENSION(nlonpe,nlatpe,nmetlev),   INTENT(IN) :: t
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: p
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: t0
      REAL, DIMENSION(nlnpe,nlpe,nlev),         INTENT(IN) :: cloud
      REAL, DIMENSION(nlnpe,nlpe,nlev),         INTENT(IN) :: o3conc
      REAL, DIMENSION(nlnpe,nlpe),              INTENT(IN) :: o3col
      REAL, DIMENSION(nlnpe,nlpe),              INTENT(IN) :: land
      REAL, DIMENSION(nlnpe,nlpe),              INTENT(IN) :: snowfr
      REAL, DIMENSION(nlnpe,nlpe),              INTENT(IN) :: seaicefr
      REAL,                                     INTENT(IN) :: time

      REAL, DIMENSION(ndj,nlnpe,nlpe,nlev),     INTENT(OUT) :: dja
      REAL, DIMENSION(ndj,nclprc),              INTENT(OUT) :: dj

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: l
      INTEGER :: m

      REAL, DIMENSION(0:nlev) :: c
      REAL, DIMENSION(0:nlev) :: pp
      REAL, DIMENSION(0:nlev) :: et
      REAL, DIMENSION(nlev)   :: o3

      REAL :: latit
      REAL :: longit
      REAL :: jtime

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev) :: x
!kk
      INTEGER,DIMENSION(nlnpe)     :: lon_a,lat_a
      REAL,DIMENSION(nlnpe,0:nlev) :: pp_a,et_a

      CHARACTER*72 :: cmessage

      integer         :: jlb       ! Number of day elements
      integer         :: jlb_phot  ! Numer of balanced elements

!kk   variables to collect data in 1-D lateral array

      REAL,DIMENSION(ndj,1,max_stochem_points,nlev) :: dja_lb
      REAL,DIMENSION(0:nlev,max_stochem_points)     :: c_lb
      REAL,DIMENSION(max_stochem_points)            :: o3col_lb
      REAL,DIMENSION(nlev,max_stochem_points)       :: o3_lb
      REAL,DIMENSION(max_stochem_points)            :: land_lb
      REAL,DIMENSION(max_stochem_points)            :: snowfr_lb
      REAL,DIMENSION(max_stochem_points)            :: seaicefr_lb
      REAL,DIMENSION(max_stochem_points)            :: zenith_lb
      REAL,DIMENSION(0:nlev,max_stochem_points)     :: pp_lb
      REAL,DIMENSION(0:nlev,max_stochem_points)     :: et_lb

      LOGICAL,DIMENSION(nlnpe,nlpe)         :: to_do
      REAL                                  :: zenith
      REAL                                  :: zen

!kk   variables for load balancing across PEs
!
      INTEGER, PARAMETER            :: move_limit=5
      INTEGER, SAVE                 :: data_move_direction
                               ! -1 get data from remote PE
                               !  0 no data move
                               !  1 put data to remote PE
      INTEGER, DIMENSION(0:nproc-1) :: all_jlb
      INTEGER                       :: my_partner
      INTEGER                       :: nr_ele_to_move

!kk   Move array from einterp to outside k and l loop

      x(:,:,0) = t0(:,:)
      x(:,:,1:nmetlev) = t(:,:,1:nmetlev)

! Calculate value for mid time_step
      jtime = time + stochem_advection_step / 2.0

! Establish cloud and ozone profiles.
      jlb = 0
      DO k=1,nlpe
        latit=(lat(k+ltdat-1-1)+lat(k+ltdat-1))/2.0 - 90.0
        DO l=1,nlnpe
           lon_a(l) = l+lndat-1
           lat_a(l) = k+ltdat-1
        END DO
!       Interpolate T and P profile to Eulerian grid centres.
        CALL EINTERP(nlnpe,lon_a,lat_a,p,et_a,pp_a,x)

        DO l=1,nlnpe
          pp = pp_a(l,:)
          et = et_a(l,:)
          IF (l+lndat-1 < nlong) THEN
            longit = (long(l+lndat-1) + long(l+lndat-1+1)) / 2.0
          ELSE
            longit = (long(l+lndat-1) + long(MOD(l+lndat-1,nlong)+1)    &
     &        + 360.0) / 2.0
          END IF
          c(1:nlev) = cloud(l,k,:)
          o3 = o3conc(l,k,:)
          c(0) = 0.0

! Check ozone remains positive
          IF (ANY(o3 < 0.0)) THEN
            cmessage='Negative o3 in jvals routine'
            WRITE(6,*) cmessage,'O3: ', O3
          END IF
          IF (o3col(l,k) <= 0.0) THEN
            cmessage='Negative or zero o3col in jvals routine'
            WRITE(6,*) cmessage,l,k,o3col(l,k)
! DEPENDS ON: ereport
            CALL EREPORT('Jvals',1,cmessage)
          END IF

! DEPENDS ON: zen
          zenith = ZEN(jtime,latit,longit)
          to_do(l,k) = (COS(ZENITH) > 0.0)      ! lateral day point

! Prepare data for column-photolysis model.
          IF (to_do(l,k)) THEN
             jlb = jlb+1
             c_lb(:,jlb)      = c
             o3col_lb(jlb)    = o3col(l,k)
             o3_lb(:,jlb)     = o3
             land_lb(jlb)     = land(l,k)
             snowfr_lb(jlb)   = snowfr(l,k)
             seaicefr_lb(jlb) = seaicefr(l,k)
             zenith_lb(jlb)   = zenith
             pp_lb(:,jlb)     = pp
             et_lb(:,jlb)     = et
          END IF
        END DO
      END DO

! Equalise number of photolysis calculations on each PE
      CALL BALANCE

! Call photolysis model.
      CALL PHOT(jlb_phot,dja_lb(:,1,1:jlb_phot,:),                      &
     &               c_lb(:,1:jlb_phot),                                &
     &               o3col_lb(1:jlb_phot),                              &
     &               o3_lb(:,1:jlb_phot),                               &
     &               land_lb(1:jlb_phot),                               &
     &               snowfr_lb(1:jlb_phot),                             &
     &               seaicefr_lb(1:jlb_phot),                           &
     &               zenith_lb(1:jlb_phot),                             &
     &               jtime,                                             &
     &               pp_lb(:,1:jlb_phot),                               &
     &               et_lb(:,1:jlb_phot))

! Get results from remote PE

! Move results back to original PE
      CALL MOVE_DATA_BACK

      j = 0
      DO k=1,nlpe
        DO l=1,nlnpe
          IF (to_do(l,k)) THEN
            j = j+1
            dja(:,l,k,:) = dja_lb(:,1,j,:)
          ELSE
            dja(:,l,k,:) = 0.0
          END IF
        END DO
      END DO

      DO j=1,nfill
        i=ipos(1,j)-lndat+1
        k=ipos(2,j)-ltdat+1
        l=ipos(3,j)
        dj(:,j)=dja(:,i,k,l)
      END DO

      CONTAINS

      SUBROUTINE BALANCE
!
! Move data to remote PE's to avoid load imbalance in phot.
! The load imbalance in phot is related to the different day-night
! distribution in the northern and southern hemispheres.
! Data exchange is only done between two PEs (mype-my_partner).
! The north-most PE exchanges data with the south-most PE...

! Data has been collected in arrays named ..._lb. These arrays have
! one dimension collecting all lateral points which are in the day time.

! The average nunber of elements on partner PE is computed and
! the difference between the actual and average number is moved.
!
! Current Owner of Code: M.G. sanderson
!
! History:
! Version   Date                    Comment
!  5.5    15/03/04  Created.  K. Ketelsen
!  6.1    21/10/04  Reformatted code. M.G. Sanderson
!
        IMPLICIT  NONE

        INTEGER :: tag, info
        INTEGER :: nr_data
        INTEGER :: recv_ind, send_ind

! Broadcast jbl to all PEs

        all_jlb       = 0
        all_jlb(mype) = jlb

        CALL GC_SSYNC(nproc,info)
        CALL GC_ISUM(nproc,nproc,info,all_jlb)
        CALL GC_SSYNC(nproc,info)

        my_partner = nproc-1-mype           ! partner to exchange Data

        nr_ele_to_move = (all_jlb(mype) - all_jlb(my_partner)) / 2

        IF (ABS(nr_ele_to_move) <= move_limit) THEN
          data_move_direction = 0           ! No action required
          jlb_phot = jlb
        ELSE IF (nr_ele_to_move < 0) THEN
          data_move_direction = -1
          jlb_phot = jlb-nr_ele_to_move
        ELSE
          data_move_direction = 1
          jlb_phot = jlb-nr_ele_to_move
        END IF

        nr_ele_to_move = ABS(nr_ele_to_move)

        tag = 12345
        IF (data_move_direction == -1) THEN ! Receive Data
          recv_ind = jlb+1
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RRECV(tag+1,nr_data,my_partner,info,                  &
     &      c_lb(0,recv_ind),c_lb(0,recv_ind))
          nr_data = nr_ele_to_move
          CALL GC_RRECV(tag+2,nr_data,my_partner,info,                  &
     &      o3col_lb(recv_ind),o3col_lb(recv_ind))
          nr_data = nlev*nr_ele_to_move
          CALL GC_RRECV(tag+3,nr_data,my_partner,info,                  &
     &                  o3_lb(1,recv_ind),o3_lb(1,recv_ind))
          nr_data = nr_ele_to_move
          CALL GC_RRECV(tag+4,nr_data,my_partner,info,                  &
     &                   land_lb(recv_ind),land_lb(recv_ind))
          nr_data = nr_ele_to_move
          CALL GC_RRECV(tag+5,nr_data,my_partner,info,                  &
     &                   snowfr_lb(recv_ind),snowfr_lb(recv_ind))
          nr_data = nr_ele_to_move
          CALL GC_RRECV(tag+6,nr_data,my_partner,info,                  &
     &                   seaicefr_lb(recv_ind),seaicefr_lb(recv_ind))
          nr_data = nr_ele_to_move
          CALL GC_RRECV(tag+8,nr_data,my_partner,info,                  &
     &                  zenith_lb(recv_ind),zenith_lb(recv_ind))
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RRECV(tag+9,nr_data,my_partner,info,                  &
     &                   pp_lb(0,recv_ind),pp_lb(0,recv_ind))
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RRECV(tag+0,nr_data,my_partner,info,                  &
     &                   et_lb(0,recv_ind),et_lb(0,recv_ind))
        ELSE IF (data_move_direction == 1) THEN   !send data
          send_ind = jlb-nr_ele_to_move+1
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RSEND(tag+1,nr_data,my_partner,info,                  &
     &                   c_lb(0,send_ind),c_lb(0,send_ind))
          nr_data = nr_ele_to_move
          CALL GC_RSEND(tag+2,nr_data,my_partner,info,                  &
     &                    o3col_lb(send_ind),o3col_lb(send_ind))
          nr_data = nlev*nr_ele_to_move
          CALL GC_RSEND(tag+3,nr_data,my_partner,info,                  &
     &                    o3_lb(1,send_ind),o3_lb(1,send_ind))
          nr_data = nr_ele_to_move
          CALL GC_RSEND(tag+4,nr_data,my_partner,info,                  &
     &                   land_lb(send_ind),land_lb(send_ind))
          nr_data = nr_ele_to_move
          CALL GC_RSEND(tag+5,nr_data,my_partner,info,                  &
     &                   snowfr_lb(send_ind),snowfr_lb(send_ind))
          nr_data = nr_ele_to_move
          CALL GC_RSEND(tag+6,nr_data,my_partner,info,                  &
     &                   seaicefr_lb(send_ind),seaicefr_lb(send_ind))
          nr_data = nr_ele_to_move
          CALL GC_RSEND(tag+8,nr_data,my_partner,info,                  &
     &                   zenith_lb(send_ind),zenith_lb(send_ind))
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RSEND(tag+9,nr_data,my_partner,info,                  &
     &                    pp_lb(0,send_ind),pp_lb(0,send_ind))
          nr_data = (nlev+1)*nr_ele_to_move
          CALL GC_RSEND(tag+0,nr_data,my_partner,info,                  &
     &                   et_lb(0,send_ind),et_lb(0,send_ind))
        END IF

        CALL GC_SSYNC(nproc,info)

      RETURN

      END SUBROUTINE BALANCE

      SUBROUTINE MOVE_DATA_BACK
!
! Current Owner of Code: M.G. sanderson
!
! History:
! Version   Date                    Comment
!  5.5    15/03/04  Created.  K. Ketelsen
!  6.1    21/10/04  Reformatted code. M.G. Sanderson
!
        IMPLICIT  NONE

        INTEGER :: tag, info
        INTEGER :: nr_data
        INTEGER :: recv_ind, send_ind
        INTEGER :: k

        tag = 23457
        CALL GC_SSYNC(NPROC,INFO)

        IF (data_move_direction == -1) THEN       ! Send Data
! send data
          send_ind = jlb+1
          nr_data  = ndj*nr_ele_to_move
          DO k=1,nlev
            CALL GC_RSEND(tag+k,nr_data,my_partner,info,                &
     &        dja_lb(1,1,send_ind,k),dja_lb(1,1,send_ind,k))
          END DO
        ELSE IF (data_move_direction == 1) THEN   ! Receive data
! receive data
           recv_ind = jlb-nr_ele_to_move+1
           nr_data  = ndj*nr_ele_to_move
           DO k=1,nlev
             CALL GC_RRECV(tag+k,nr_data,my_partner,info,               &
     &         dja_lb(1,1,recv_ind,k),dja_lb(1,1,recv_ind,k))
           END DO
        END IF

        CALL GC_SSYNC(NPROC,INFO)

        RETURN

      END SUBROUTINE MOVE_DATA_BACK

      END SUBROUTINE JVALS
#endif
