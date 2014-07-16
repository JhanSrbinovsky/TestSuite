#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE BTQ_INT
!!!  Purpose: To interpolate buoyancy parameters BT and BQ from full
!!!  levels to half levels
!!!
!!! HISTORY:
!!! DATE   VERSION   COMMENT
!!! ----   -------   -------
!!!
!!! 10/9/97  4.4     New Deck.  R.N.B.Smith
!!! 14/3/00  5.1     Addition of CFM calculation.  A.P.Lock
!!!
!!! CODE DESCRIPTION:
!!!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!!!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.
!!!
!!!---------------------------------------------------------------------
      SUBROUTINE BTQ_INT (                                              &
     & row_length, rows, BL_LEVELS                                      &
     &,z_full,z_half,BQ,BT,BQ_CLD,BT_CLD,A_QS,A_DQSDT                   &
     &,BQM,BTM,BQM_CLD,BTM_CLD,A_QSM,A_DQSDTM                           &
     &,CF,CFM                                                           &
     &,LTIMER                                                           &
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics

      INTEGER                                                           &
     & row_length                                                       &
     &,rows                                                             &
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed ! <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA
      REAL                                                              &
     & z_full(row_length,rows,BL_LEVELS)                                &
     &,z_half(row_length,rows,BL_LEVELS)                                &
     &,BQ(row_length,rows,BL_LEVELS)                                    &
                              ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BT(row_length,rows,BL_LEVELS)                                    &
                              ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BQ_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on p,T,q-levels (full levels).
     &,BT_CLD(row_length,rows,BL_LEVELS)                                &
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on p,T,q-levels (full levels).
     &,A_QS(row_length,rows,BL_LEVELS)                                  &
!                             ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).
     &,A_DQSDT(row_length,rows,BL_LEVELS)                               &
!                             ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).
     &,CF(row_length,rows,BL_LEVELS)
!                             ! IN Cloud fraction (decimal).

      REAL                                                              &
            ! OUT arrays,
     & BQM(row_length,rows,BL_LEVELS)                                   &
                              ! OUT A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM(row_length,rows,BL_LEVELS)                                   &
                              ! OUT A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BQM_CLD(row_length,rows,BL_LEVELS)                               &
!                             ! OUT A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM_CLD(row_length,rows,BL_LEVELS)                               &
!                             ! OUT A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_QSM(row_length,rows,BL_LEVELS)                                 &
!                             ! OUT Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_DQSDTM(row_length,rows,BL_LEVELS)                              &
!                             ! OUT Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,CFM(row_length,rows,BL_LEVELS)
!                             ! OUT Estimate of cloud fraction
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.


!-----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!    Local and other symbolic constants :-

!  Define local storage.

!  (b) Scalars.

      REAL                                                              &
     & WK                                                               &
               ! Temporary in weighting factor.
     &,WKM1                                                             &
               ! Temporary in weighting factor.
     &,weight1,weight2,weight3

      INTEGER                                                           &
     & I,J                                                              &
               ! Loop counter (horizontal field index).
     &,K       ! Loop counter (vertical level index).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BTQ_INT ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round levels.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        DO J=1, rows
        DO I=1, row_length

!-----------------------------------------------------------------------
!! 1.1 Calculate buoyancy parameters at half levels,
!!     i.e. at level K-1/2, if current level is level K.
!-----------------------------------------------------------------------

          WEIGHT1 = 1. / ( z_full(I,J,K) -                              &
     &                     z_full(I,J,K-1))
          WEIGHT2 = z_full(I,J,K) -                                     &
     &              z_half(I,J,K)
          WEIGHT3 = z_half(I,J,K) -                                     &
     &              z_full(I,J,K-1)
          WKM1 = weight3 * weight1         ! P243.C5 (2nd eqn)
          WK = weight2 * weight1             ! P243.C5 (1st eqn)

          BTM(I,j,K-1) = WKM1*BT(I,j,K) + WK*BT(I,j,K-1)
          BQM(I,j,K-1) = WKM1*BQ(I,j,K) + WK*BQ(I,j,K-1)
          BTM_CLD(I,j,K-1) = WKM1*BT_CLD(I,j,K) + WK*BT_CLD(I,j,K-1)
          BQM_CLD(I,j,K-1) = WKM1*BQ_CLD(I,j,K) + WK*BQ_CLD(I,j,K-1)
          A_QSM(I,j,K-1) = WKM1*A_QS(I,j,K) + WK*A_QS(I,j,K-1)
          A_DQSDTM(I,j,K-1) = WKM1*A_DQSDT(I,j,K) + WK*A_DQSDT(I,j,K-1)
!          CFM(I,j,K-1) = WKM1*CF(I,j,K) + WK*CF(I,j,K-1)
!         ! to compensate for lack of resolution of cloud layers
          CFM(I,j,K-1) = MAX( CF(I,j,K), CF(I,j,K-1) )


        ENDDO !
        ENDDO !
      ENDDO ! bl_levels

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('BTQ_INT ',4)
      ENDIF

      RETURN
      END SUBROUTINE BTQ_INT
#endif
