#if defined(A70_1B) || defined(A70_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to rescale the asymmetry.
!
! Method:
!       The standard rescaling of the asymmetry is used.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_ASYMMETRY(N_PROFILE                            &
     &   , I_LAYER_FIRST, I_LAYER_LAST                                  &
     &   , ASYMMETRY, FORWARD_SCATTER                                   &
     &   , NPD_PROFILE, NPD_LAYER, l_pc2                                &
     &   )
!
!
      IMPLICIT NONE
!
!
      Logical                                                           &
                !, INTENT(IN)
     &     L_pc2  ! Use PC2 cloud scheme
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
!             LAST LAYER TO RESCALE
      REAL                                                              &
                !, INTENT(IN)
     &     FORWARD_SCATTER(NPD_PROFILE, NPD_LAYER)
!             FORWARD SCATTERING
      REAL                                                              &
                !, INTENT(INOUT)
     &     ASYMMETRY(NPD_PROFILE, NPD_LAYER)
!             ASYMMETRY
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
      IF (L_pc2) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
           DO L=1, N_PROFILE

             ! Protect against forward scatter being larger than 1
             IF (FORWARD_SCATTER(L, I)  >=  1.0) THEN
               ASYMMETRY(L, I) = 0.5
               FORWARD_SCATTER(L, I) = 1.0
             ELSE
              ASYMMETRY(L, I)=(ASYMMETRY(L, I)-FORWARD_SCATTER(L, I))   &
     &           /(1.0E+00-FORWARD_SCATTER(L, I))
             ENDIF
           ENDDO
        ENDDO

      ELSE  ! L_pc2
        DO I=I_LAYER_FIRST, I_LAYER_LAST
           DO L=1, N_PROFILE
            ASYMMETRY(L, I)=(ASYMMETRY(L, I)-FORWARD_SCATTER(L, I))     &
     &         /(1.0E+00-FORWARD_SCATTER(L, I))
           ENDDO
        ENDDO

      ENDIF  ! L_pc2
!
!
      RETURN
      END SUBROUTINE RESCALE_ASYMMETRY
#endif
