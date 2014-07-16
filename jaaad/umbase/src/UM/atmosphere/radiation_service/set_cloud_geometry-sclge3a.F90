#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set geometry of clouds.
!
! Method:
!       For use in multi-column mode arrays are set for each layer
!       pointing to profiles which have non-negligible clear or
!       cloudy fractions. The topmost cloudy layers are also
!       detected.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT,WHENFLE
!                                  replaced by portable fortran code.
!                                                S.J.Swarbrick
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER                  &
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL                       &
     &   , W_CLOUD                                                      &
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE                             &
     &   , N_CLOUD_TOP                                                  &
     &   , N_FREE_PROFILE, I_FREE_PROFILE                               &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
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
     &   , N_LAYER
!             NUMBER OF LAYERS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CLOUD_TOP_GLOBAL
!             GLOBAL TOPMOST CLOUDY LAYER
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , N_FREE_PROFILE(NPD_LAYER)                                    &
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)                       &
!             CLOUD-FREE PROFILES
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             AMOUNTS OF CLOUD
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I,J,L
!             LOOP VARIABLES
      REAL                                                              &
     &     TOL_CLOUD
!             TOLERANCE FOR DECIDING THAT AMOUNTS OF CLOUD ARE
!             NEGLIGIBLE
!
!
!     SUBROUTINES CALLED:
!
!
!
!     SET A TOLERANCE FOR DECIDING THAT THE FRACTION OF CLOUD IN
!     A LAYER IS NOT NEGLIGIBLE.
      TOL_CLOUD=64.0E+00*EPSILON(TOL_CLOUD)
!
      DO I=1, N_LAYER
!
         J=1
         N_CLOUD_PROFILE(I)=0
         DO L=1,N_PROFILE
           IF (W_CLOUD(L,I) >  TOL_CLOUD) THEN
             I_CLOUD_PROFILE(J,I)=L
             J=J+1
             N_CLOUD_PROFILE(I)=N_CLOUD_PROFILE(I)+1
           END IF
         END DO
!
!
         J=1
         N_FREE_PROFILE(I)=0
         DO L=1,N_PROFILE
           IF (W_CLOUD(L,I) <= TOL_CLOUD) THEN
             I_FREE_PROFILE(J,I)=L
             J=J+1
             N_FREE_PROFILE(I)=N_FREE_PROFILE(I)+1
           END IF
         END DO
!
      ENDDO
!
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP=N_CLOUD_TOP_GLOBAL
      ELSE
         N_CLOUD_TOP=1
         DO WHILE ( (N_CLOUD_TOP <  N_LAYER).AND.                       &
     &              (N_CLOUD_PROFILE(N_CLOUD_TOP) == 0) )
            N_CLOUD_TOP=N_CLOUD_TOP+1
         ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SET_CLOUD_GEOMETRY
#endif
#endif
