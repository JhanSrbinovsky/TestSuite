#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE RKO4_SINGLE_TIME(pos,u,v,w,bl,orog,seed2,              &
     &   nfill,cellswap,cellbase,z_top_of_model,                        &
     &   first_constant_r_rho_level)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods :  To advect cells in the windfield.
!-
!-   Inputs  : POS,U,V,W,BL,OROG,SEED2,NFILL
!-   Outputs : POS
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.3    21/09/01  Created. C.E. Johnson
!  5.5    21/01/04  Vectorised version. K. Ketelsen
!  6.1    23/09/04  Now uses integer random number seeds. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!-
!VVV  v5.2.1  RKO4_SINGLE_TIME 21/IX/01 -
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE VELOC_MOD                 !kk
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER,                INTENT(IN) :: nfill
      INTEGER,                INTENT(IN) :: cellbase
      INTEGER,                INTENT(IN) :: first_constant_r_rho_level

      REAL, INTENT(IN) :: z_top_of_model
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: u ! wind
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: v
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: w !
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe),           INTENT(IN) :: orog

      INTEGER, DIMENSION(2,nclprc*nproc),    INTENT(INOUT) :: cellswap

      INTEGER, DIMENSION(ransize),           INTENT(INOUT) :: seed2
      REAL, DIMENSION(4,nclprc),             INTENT(INOUT) :: pos

      INTEGER                :: k
      INTEGER                :: pe2
      INTEGER                :: im
      INTEGER                :: jm

      REAL                   :: half
      REAL                   :: sixth
      REAL                   :: vdiff
      REAL                   :: hdiff
      REAL                   :: Zorog
      REAL                   :: Zbl
      REAL, DIMENSION(4)     :: npos
      REAL, DIMENSION(3)     :: v1
      REAL, DIMENSION(3)     :: v2
      REAL, DIMENSION(3)     :: v3
      REAL, DIMENSION(3)     :: v4
      REAL, DIMENSION(nfill) :: su
      REAL, DIMENSION(nfill) :: sv
      REAL, DIMENSION(nfill) :: sw

!kk   local variables defined by Klaus Ketelsen

      REAL, DIMENSION(nfill)    :: vdiff_a, hdiff_a
      REAL, DIMENSION(3,nfill)  :: v1_a, v2_a, v3_a, v4_a
      REAL, DIMENSION(4,nfill)  :: kpos
      LOGICAL, DIMENSION(nfill) :: todo

      cellswap = 0

      half = Stochem_advection_step / 2.0
      sixth = Stochem_advection_step / 6.0

! S's are fractions of the standard deviations to add
! DEPENDS ON: gauss_ran
      su = GAUSS_RAN(seed2,ransize,nfill)
! DEPENDS ON: gauss_ran
      sv = GAUSS_RAN(seed2,ransize,nfill)
! DEPENDS ON: gauss_ran
      sw = GAUSS_RAN(seed2,ransize,nfill)

      DO k = 1, nfill
! DEPENDS ON: getmetpoint
        zorog = GETMETPOINT(pos(:,k),orog,.TRUE.)
! DEPENDS ON: getmetpoint
        zbl = GETMETPOINT(pos(:,k),bl,.TRUE.)
        IF (pos(4,k) - earth_radius - zorog <= zbl) THEN    ! In B.L.
          vdiff_a(k) = vdiff_bl
          hdiff_a(k) = hdiff_bl
        ELSE
          vdiff_a(k) = vdiff_trop
          hdiff_a(k) = hdiff_trop
        END IF
      END DO

! 1) Find winds at initial position and find new position at
!    time+astep/2.
      todo = .true.
      CALL VELOC(nfill,todo,pos(:,1:nfill),u,v,w,su,sv,sw,v1_a,         &
     &  hdiff_a,vdiff_a)

! Check parcel is not within 3 grid squares of the poles.
      DO k = 1, nfill
        todo(k)=(pos(2,k) > 3.0*dlatm .AND. pos(2,k) < 180.0-3.0*dlatm)
      END DO

      DO k=1,nfill
        v1 = v1_a(:,k)
        IF (todo(k)) THEN
          npos((/1,2,4/)) = pos((/1,2,4/),k) + v1*half
! DEPENDS ON: repos
          CALL REPOS(npos,bl,orog,z_top_of_model,                       &
     &      first_constant_r_rho_level)

          kpos(:,k) = npos
        END IF
      END DO

! 2) Estimate velocities at time+astep/2, using position from step 1,
!    and use to find new position at time+astep/2.
      CALL VELOC(nfill,todo,kpos,u,v,w,su,sv,sw,v2_a,hdiff_a,vdiff_a)

      DO k=1,nfill
        v2 = v2_a(:,k)

        IF (todo(k)) THEN
!kk          npos((/1,2,4/))=pos((/1,2,4/),k)+v2*half
          npos(1) = pos(1,k) + v2(1)*half
          npos(2) = pos(2,k) + v2(2)*half
          npos(4) = pos(4,k) + v2(3)*half
          npos(3) = kpos(3,k)
! DEPENDS ON: repos
          CALL REPOS(npos,bl,orog,z_top_of_model,                       &
     &      first_constant_r_rho_level)
          kpos(:,k) = npos
        END IF
      END DO

! 3) Estimate velocities at time+astep/2 using position from step 2,
!    and use to estimate position at time+astep.
      CALL VELOC(nfill,todo,kpos,u,v,w,su,sv,sw,v3_a,hdiff_a,vdiff_a)

      DO k=1,nfill
        v3 = v3_a(:,k)
        IF (todo(k)) THEN
          npos((/1,2,4/)) = pos((/1,2,4/),k) + v3*stochem_advection_step
          npos(3) = kpos(3,k)
! DEPENDS ON: repos
          CALL REPOS(npos,bl,orog,z_top_of_model,                       &
     &      first_constant_r_rho_level)
          kpos(:,k) = npos
        END IF
      END DO

! 4) Estimate velocities at time+astep using position from step 3.
      CALL VELOC(nfill,todo,kpos,u,v,w,su,sv,sw,v4_a,hdiff_a,vdiff_a)
      DO k = 1, nfill
        vdiff = vdiff_a(k)
        hdiff = hdiff_a(k)
        v1 = v1_a(:,k)
        v2 = v2_a(:,k)
        v3 = v3_a(:,k)
        v4 = v4_a(:,k)

        IF (todo(k)) THEN
! 5) Estimate final position.
          pos((/1,2,4/),k) = pos((/1,2,4/),k)+sixth*(v1+v4+2.0*(v2+v3))
        ELSE
! If near the poles, then just use simple forward timestep.
          pos((/1,2,4/),k) = pos((/1,2,4/),k)+v1*stochem_advection_step
        END IF

! Keep horizontal coordinates within globe.
        IF (pos(2,k) < 0.0) THEN
          pos(2,k) = -pos(2,k)
          pos(1,k) = pos(1,k) + 180.0
        END IF
        IF (pos(2,k) > 180.0) THEN
          pos(2,k) = 360.0 - pos(2,k)
          pos(1,k) = pos(1,k) + 180.0
        END IF
        pos(1,k) = AMOD(pos(1,k)+360.0,360.0)

        im = 1 + INT(pos(1,k)/dlong)
        jm = 1 + INT(pos(2,k)/dlat)

! If cell has moved out of the range of the processor,
! flag it to be swapped.
        IF (procmap(im,jm) /= mype) THEN
          pe2 = procmap(im,jm)

! swap cell K from mype to pe2
! Add 1 to CELLSWAP so that ME are greater than 0, then 0 can imply
! no swapping. Have to remember to subtract it off later.
          cellswap(1:2,k+cellbase) = (/ mype+1,pe2+1 /)
        END IF
      END DO

      END SUBROUTINE RKO4_SINGLE_TIME
#endif
