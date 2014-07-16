


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, N_LAYER, NCLDS    &
     &   , I_CLOUD, W_CLOUD_IN, TOTAL_CLOUD_COVER                       &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DECLARATION OF ARRAY SIZES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     COMDECKS INCLUDED
! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

! CLSCHM3A end
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             Number of layers seen in radiation
     &   , NCLDS                                                        &
!             NUMBER OF CLOUDY LAYERS
     &   , I_CLOUD
!             CLOUD SCHEME EMPLOYED
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD_IN(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
      REAL                                                              &
                !, INTENT(OUT)
     &     TOTAL_CLOUD_COVER(NPD_PROFILE)
!             TOTAL CLOUD COVER
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
      REAL                                                              &
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!
!
!
!     COPY W_CLOUD_IN TO W_CLOUD AND THEN CHECK THAT VALUES ARE
!     BETWEEN 0 AND 1.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              W_CLOUD(L,I) = MIN( MAX(W_CLOUD_IN(L,I),0.0) , 1.0 )
            ENDDO
         ENDDO
!
!     DIFFERENT OVERLAP ASSUMPTIONS ARE CODED INTO EACH SOLVER.
!
      IF (I_CLOUD == IP_CLOUD_MIX_MAX) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, N_LAYER+1-NCLDS)
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1) >  W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)             &
     &               *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)                &
     &            *(1.0E+00-W_CLOUD(L, I))
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               TOTAL_CLOUD_COVER(L)=MAX(TOTAL_CLOUD_COVER(L)            &
     &            , W_CLOUD(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_TRIPLE) THEN
!
!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-W_CLOUD(L, N_LAYER+1-NCLDS)
         ENDDO
         DO I=N_LAYER+1-NCLDS, N_LAYER-1
            DO L=1, N_PROFILE
               IF (W_CLOUD(L, I+1) >  W_CLOUD(L, I)) THEN
                  TOTAL_CLOUD_COVER(L)=TOTAL_CLOUD_COVER(L)             &
     &               *(1.0E+00-W_CLOUD(L, I+1))/(1.0E+00-W_CLOUD(L, I))
               ENDIF
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=1.0E+00-TOTAL_CLOUD_COVER(L)
         ENDDO
!
      ELSE IF (I_CLOUD == IP_CLOUD_CLEAR) THEN
!
         DO L=1, N_PROFILE
            TOTAL_CLOUD_COVER(L)=0.0E+00
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE R2_CALC_TOTAL_CLOUD_COVER
