#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assemble lists of points with the same surfaces.
!
! Method:
!       The surfaces at the bottom of each profile are examined
!       and lists of points with the same type of surface are made.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENEQ replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE COLLECT_SURFACE(N_PROFILE                              &
     &   , I_SURFACE                                                    &
     &   , N_POINT_TYPE, INDEX_SURFACE                                  &
     &   , NPD_PROFILE, NPD_SURFACE                                     &
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
     &   , NPD_SURFACE
!             MAXIMUM NUMBER OF SURFACES
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_SURFACE(NPD_PROFILE)
!             SURFACE SPECIFICATIONS
      INTEGER                                                           &
                !, INTENT(OUT)
     &     N_POINT_TYPE(NPD_SURFACE)                                    &
!             NUMBER OF POINTS OF EEACH TYPE
     &   , INDEX_SURFACE(NPD_PROFILE, NPD_SURFACE)
!             LIST OF POINTS OF EACH TYPE
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , J
!
!
!     PASS THROUGH ALL THE COLUMNS COLLECTING LISTS OF POINTS WHICH
!     HAVE THE SAME SURFACE TYPE.
      DO K=1, NPD_SURFACE
         N_POINT_TYPE(K)=0
      ENDDO
!
!
      DO   K=1,NPD_SURFACE
           J=1
        DO L=1,N_PROFILE
          IF (I_SURFACE(L) == K) THEN
            INDEX_SURFACE(J,K)=L
            J=J+1
            N_POINT_TYPE(K)=N_POINT_TYPE(K)+1
          END IF
        END DO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE COLLECT_SURFACE
#endif
#endif
