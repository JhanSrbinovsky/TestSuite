#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DRYDEP(bl,gsf,ra,rq,rc,o3_stom_frac,ddepo,sddat)
!
!  6.2   01/03/06   Rewritten to use resistances from new dry
!                   deposition scheme. M.G. Sanderson.
! Combine surface resistance rc with aerodynamic resistance ra and
! quasi-laminar resistance rq to get overall dry deposition velocity.
! Then calculate first-order loss rate, and interpolate to STOCHEM
! grid.
!
!  Model            Modification history from model version 6.0:
! version  Date
!
!  6.0   12/09/03   Fix data initialisation problems encountered on
!                   some platforms. Introduce standard UM modification
!                   history.                                   P.Dando
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
!
      IMPLICIT NONE
!
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: bl
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(IN) :: gsf
      REAL, DIMENSION(nlonpe,nlatpe,ntype), INTENT(IN) :: ra
      REAL, DIMENSION(nlonpe,nlatpe,nc), INTENT(IN) :: rq
      REAL, DIMENSION(nlonpe,nlatpe,ntype,nc), INTENT(IN) :: rc
      REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: o3_stom_frac
      REAL, DIMENSION(nc,nlnpe,nlpe), INTENT(OUT) :: ddepo
      REAL, DIMENSION(nc+1,nlnpe,nlpe), INTENT(OUT) :: sddat
!
      INTEGER :: i, j, k, n, nj
      REAL, DIMENSION(nlonpe,nlatpe,ntype,n_spec_deposit) :: vd
      REAL, DIMENSION(nlonpe,nlatpe,n_spec_deposit+1) :: sdd
      REAL, DIMENSION(nlonpe,nlatpe,n_spec_deposit+1) :: ddfolr
      REAL, DIMENSION(nlnpe,nlpe,n_spec_deposit+1) :: ddtemp

      DO j = 1, n_spec_deposit
        DO k = 1, rowspe
          DO i = 1, nlonpe
            sdd(i,k,j) = 0.0
            ddfolr(i,k,j) = 0.0
            vd(i,k,1,j) = 0.0
            vd(i,k,2,j) = 0.0
            vd(i,k,3,j) = 0.0
            vd(i,k,4,j) = 0.0
            vd(i,k,5,j) = 0.0
            vd(i,k,6,j) = 0.0
            vd(i,k,7,j) = 0.0
            vd(i,k,8,j) = 0.0
            vd(i,k,9,j) = 0.0
          END DO
        END DO
      END DO
      j = n_spec_deposit + 1
      DO k = 1, rowspe
        DO i = 1, nlonpe
          sdd(i,k,j) = 0.0
          ddfolr(i,k,j) = 0.0
        END DO
      END DO

      DO k = 1, rowspe
        DO j = 1, n_spec_deposit
          nj = dep_spec(j)
          DO n = 1, npft
            DO i = 1, nlonpe
!
! Do vegetated tiles first. Quasi-laminar resistance pre-multiplied by
! ln[z0m/z0] = 2.0 for vegetated areas, or 1.0 for smooth surfaces
! See Ganzeveld & Lelieveld, JGR 1995 Vol.100 No. D10 pp.20999-21012.
!
              IF (rc(i,k,n,nj) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
                vd(i,k,n,j) = 1.0 /                                     &
     &            (ra(i,k,n) + 2.0*rq(i,k,nj) + rc(i,k,n,nj))
                sdd(i,k,j) = sdd(i,k,j) + gsf(i,k,n) / rc(i,k,n,nj)
              END IF
            END DO
          END DO
!
! Now do calculation for non-vegetated tiles
!
          DO n = npft+1, ntype
            DO i = 1, nlonpe
              IF (rc(i,k,n,nj) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
                vd(i,k,n,j) = 1.0 /                                     &
     &            (ra(i,k,n) + rq(i,k,nj) + rc(i,k,n,nj))
                sdd(i,k,j) = sdd(i,k,j) + gsf(i,k,n) / rc(i,k,n,nj)
              END IF
            END DO
          END DO
        END DO
      END DO
!
! VD() now contains dry deposition velocities for each tile in each grid
! Calculate overall first-order loss rate over time ASTEP for each tile
! and sum over all tiles to obtain overall first-order loss rate DDFOLR
!
      DO n = 1, ntype
        DO j = 1, n_spec_deposit
          DO k = 1, rowspe
            DO i = 1, nlonpe
              IF (vd(i,k,n,j) > 0.0) THEN
                ddfolr(i,k,j) = ddfolr(i,k,j) + gsf(i,k,n) *            &
     &            (1.0-EXP(-vd(i,k,n,j) * stochem_advection_step /      &
     &            bl(i,k)))
              END IF
            END DO
          END DO
        END DO
      END DO

! DDFOLR() is loss rate constant over time stochem_advection_step.
! Divide by stochem_advection_step to get rate in s-1.

      DO k = 1, rowspe
        DO j = 1, n_spec_deposit
          DO i = 1, nlonpe
            ddfolr(i,k,j) = ddfolr(i,k,j) / stochem_advection_step
          END DO
        END DO
      END DO

! Interpolate the loss rate constants from the UM grid to the
! STOCHEM grid

! DEPENDS ON: met2data
      CALL MET2DATA(ddtemp,ddfolr,n_spec_deposit+1,n_spec_deposit+1)

      DO i = 1, nlnpe
        DO k = 1, nlpe
          DO j = 1, n_spec_deposit
            ddepo(dep_spec(j),i,k) = ddtemp(i,k,j)
          END DO
        END DO
      END DO

! Add fraction of ozone depositing via stomata to SDD

      j = n_spec_deposit+1
      DO k = 1, rowspe
        DO i = 1, nlonpe
          sdd(i,k,j) = o3_stom_frac(i,k)
        END DO
      END DO

! Interpolate the surface resistances from the UM grid to the
! STOCHEM grid

! DEPENDS ON: met2data
      CALL MET2DATA(ddtemp,sdd,n_spec_deposit+1,n_spec_deposit+1)

      DO i = 1, nlnpe
        DO k = 1, nlpe
          DO j = 1, n_spec_deposit
            sddat(dep_spec(j),i,k) = ddtemp(i,k,j)
          END DO
          sddat(nc+1,i,k) = ddtemp(i,k,n_spec_deposit+1) ! O3 stom frac
        END DO
      END DO

      END SUBROUTINE DRYDEP
#endif
