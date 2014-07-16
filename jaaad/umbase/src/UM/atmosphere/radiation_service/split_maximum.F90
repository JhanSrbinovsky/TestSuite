#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to split the atmosphere into maximally overlapped columns.
!
! Method:
!       The layers are first ranked in order of increasing cloudiness.
!       This operation cannot be vectorized and is done for one profile
!       at a time. The areal extent of each column and the logical
!       cloud mask are then set.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
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
      SUBROUTINE SPLIT_MAXIMUM(N_PROFILE, N_LAYER                       &
     &   , W_CLOUD                                                      &
     &   , N_COLUMN, AREA_COLUMN, L_COLUMN                              &
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
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
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     DUMMY ARGUMENTS
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      INTEGER                                                           &
                !, INTENT(INOUT)
     &     N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             ARRAY OF TYPES
      REAL                                                              &
                !, INTENT(IN)
     &     AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)                         &
!             AREA OF EACH COLUMN
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
!     LOCAL ARGUMENTS
      INTEGER                                                           &
     &     IRANK(NPD_LAYER)                                             &
!             ARRAY TO RANK COLUMNS BY W
     &   , I                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIBLE
     &   , L
!             LOOP VARIBLE
      REAL                                                              &
     &     W_CLOUD_SINGLE(NPD_LAYER)                                    &
!             CLOUD AMOUNTS FOR SINGLE PROFILE
     &   , W                                                            &
!             SINGLE CLOUD AMOUNT
     &   , TOL_COLUMN
!             MINIMUM FRACTIONAL AREA NECESSARY TO DEFINE
!             A NEW COLUMN
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     RANK
!
!
!
!     A NEW COLUMN WILL BE DEFINED IF IT IS SUFFICIENTLY LARGE
!     IN RELATION TO THE PRECISION OF THE MACHINE.
      TOL_COLUMN=64.0E+00*EPSILON(TOL_COLUMN)
!
      DO L=1, N_PROFILE
!        GATHER THE CLOUD AMOUNTS FOR ONE PROFILE
         DO I=1, N_LAYER
            W_CLOUD_SINGLE(I)=W_CLOUD(L, I)
         ENDDO
!
!        FIRST FORM THE VECTOR IRANK, RANKING THE LAYERS IN ORDER OF
!        INCREASING CLOUD CONTENT.
! DEPENDS ON: rank
         CALL RANK(N_LAYER                                              &
     &      , W_CLOUD_SINGLE, IRANK                                     &
     &   , NPD_LAYER                                                    &
     &      )
!
!        PASS THROUGH ALL THE COLUMNS SETTING L_COLUMN EQUAL TO .FALSE.
!        IF THE COLUMN IS ACTUALLY CLEAR ON THAT LEVEL. THE ASSUMPTION
!        OF MAXIMUM OVERLAP IS USED HERE.
         N_COLUMN(L)=1
         W=0.0E+00
         I=1
30       IF (I <= N_LAYER) THEN
            IF ( W_CLOUD_SINGLE(IRANK(I)) <  (W+TOL_COLUMN) ) THEN
               I=I+1
               GOTO 30
            ELSE
               DO K=1, I-1
                  L_COLUMN(L, IRANK(K), N_COLUMN(L))=.FALSE.
               ENDDO
               DO K=I, N_LAYER
                  L_COLUMN(L, IRANK(K), N_COLUMN(L))=.TRUE.
               ENDDO
               AREA_COLUMN(L, N_COLUMN(L))                              &
     &            =W_CLOUD_SINGLE(IRANK(I))-W
               W=W_CLOUD_SINGLE(IRANK(I))
            ENDIF
            IF (W <  W_CLOUD_SINGLE(IRANK(N_LAYER))-TOL_COLUMN) THEN
               N_COLUMN(L)=N_COLUMN(L)+1
               GOTO 30
            ENDIF
         ENDIF
!
!        THERE IS A TOTALLY CLEAR COLUMN UNLESS AT LEAST ONE LAYER IS
!        TOTALLY CLOUDY.
         IF ((1.0E+00-W) >  TOL_COLUMN) THEN
!           INCREMENT THE NUMBER OF COLUMNS IF THE FIRST IS NOT BLANK.
            IF (W >= TOL_COLUMN) THEN
               N_COLUMN(L)=N_COLUMN(L)+1
            ENDIF
            DO K=1, N_LAYER
               L_COLUMN(L, IRANK(K), N_COLUMN(L))=.FALSE.
            ENDDO
            AREA_COLUMN(L, N_COLUMN(L))=1.0E+00-W
         ENDIF
!
      ENDDO
!
!
      RETURN
      END SUBROUTINE SPLIT_MAXIMUM
#endif
#endif
