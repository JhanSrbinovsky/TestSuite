#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to rescale optical depth and albedo.
!
! Method:
!       The standard rescaling formulae are applied.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.5             11-06-98                Optimised version
!                                               (P. Burton)
!  6.0  21/08/03  NEC SX-6 optimisation - add rewritten vectorisable
!                 loop under defined(NEC).  R Barnes & J-C Rioual.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_TAU_OMEGA(N_PROFILE                            &
     &   , I_LAYER_FIRST, I_LAYER_LAST                                  &
     &   , TAU, OMEGA, FORWARD_SCATTER                                  &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO RESCALE
     &   , I_LAYER_LAST
!             FIRST LAYER TO RESCALE
      REAL                                                              &
                !, INTENT(IN)
     &     FORWARD_SCATTER(NPD_PROFILE, NPD_LAYER)
!             FORWARD SCATTERING
      REAL                                                              &
                !, INTENT(INOUT)
     &     TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE

!     Temporary scalars for the rescheduled divide
      REAL                                                              &
     &  TEMP1,TEMP2,TEMP3,TEMP4,TEMP5



#if defined(NEC)
      IF(N_PROFILE >  0) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1,N_PROFILE
            TEMP1=1.0E+00 - OMEGA(L,I)*FORWARD_SCATTER(L,I)
            TEMP2=1.0E+00 - FORWARD_SCATTER(L,I)
            TEMP5=TEMP2/TEMP1
            TAU(L,I)=TAU(L,I)*TEMP1
            OMEGA(L,I) = OMEGA(L,I)*TEMP5
          ENDDO
        END DO
      END IF
#else
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         IF(N_PROFILE >  0) THEN
            TEMP1=1.0E+00 - OMEGA(1,I)*FORWARD_SCATTER(1,I)
            TEMP2=1.0E+00 - FORWARD_SCATTER(1,I)
            DO L=1,N_PROFILE-1
               TEMP5=TEMP2/TEMP1
               TAU(L,I)=TAU(L,I)*TEMP1
               TEMP3=1.0-OMEGA(L+1,I)*FORWARD_SCATTER(L+1,I)
               TEMP4=1.0-FORWARD_SCATTER(L+1,I)
               TEMP2=TEMP4
               TEMP1=TEMP3
               OMEGA(L,I) = OMEGA(L,I)*TEMP5
            ENDDO
            TEMP5=TEMP2/TEMP1
            TAU(N_PROFILE,I) = TAU(N_PROFILE,I)*TEMP1
            OMEGA(N_PROFILE,I)= OMEGA(N_PROFILE,I)*TEMP5
         END IF
      END DO
#endif

      RETURN
      END SUBROUTINE RESCALE_TAU_OMEGA
#endif
#endif
