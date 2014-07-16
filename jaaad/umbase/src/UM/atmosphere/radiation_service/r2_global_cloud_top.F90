#if defined(A70_1B) || defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to determine the global topmost cloudy layer.
!
! Purpose:
!   The routine determines the topmost cloudy layer over the whole
!   computational domain for use in configurations of the model where
!   results must be bit-reproducible irrespective of the number of
!   segments.
!
! Method:
!   The arrays LCA for layer cloud and CCT for convective cloud are
!   searched to determine the topmost layer occupied by cloud.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             17-12-96                Original Code
!                                               (J. M. Edwards)
!        4.3             17-04-97                Delete trailing blank l
!                                                (A. Brady)
!       4.4             23-09-97                Lower limit set for
!                                               searching.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_GLOBAL_CLOUD_TOP(N_POINTS, NLEVS, NCLDS             &
!                       Convective cloud Fields
     &   , CCA, CCT                                                     &
!                       Layer cloud Fields
     &   , LCA                                                          &
!                       Calculated top of cloud fields.
     &   , GLOBAL_CLOUD_TOP                                             &
!                       Size of arrays
     &   , ND_POINTS                                                    &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
!
!     SIZES OF ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     ND_POINTS
!             HORIZONTAL DIMENSION OF ARRAYS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_POINTS                                                     &
!             NUMBER OF POINTS TO CONSIDER
     &   , NLEVS                                                        &
!             NUMBER OF LEVELS CONSIDERED
     &   , NCLDS
!             NUMBER OF CLOUDY LAYERS
!
!     FIELDS FOR CONVECTIVE CLOUDS
      REAL                                                              &
                !, INTENT(IN)
     &     CCA(ND_POINTS)
!             AMOUNTS OF CONVECTIVE CLOUD
      INTEGER                                                           &
                !, INTENT(IN)
     &     CCT(ND_POINTS)
!             LEVEL OF TOP OF CONVECTIVE CLOUD
!             I.E. CONVECTIVE CLOUD REACHES THE BASE OF LAYER CCT,
!             BUT DOES NOT EXTEND INTO IT.
!
!     FIELDS FOR STRATIFORM CLOUDS
      REAL                                                              &
                !, INTENT(IN)
     &     LCA(ND_POINTS, NCLDS)
!             AMOUNTS OF LAYER CLOUD
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     GLOBAL_CLOUD_TOP
!             GLOBAL TOPMOST CLOUDY LEVEL
!
!
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     L
!             LOOP VARIABLE
      LOGICAL                                                           &
     &     L_LAYER_CLEAR
!             FLAG FOR LAYER FREE OF CLOUDS
!
!
!
!     INITIALIZE THE CLOUD-TOP TO THE LAYER ABOVE THAT IN WHICH CLOUDS
!     ARE PERMITTED AND REDUCE UNTIL NON-EMPTY LAYERS ARE FOUND.
      GLOBAL_CLOUD_TOP=NCLDS+1
      L_LAYER_CLEAR=.TRUE.
      DO WHILE ( (L_LAYER_CLEAR).AND.(GLOBAL_CLOUD_TOP >  1) )
         GLOBAL_CLOUD_TOP=GLOBAL_CLOUD_TOP-1
         DO L=1, N_POINTS
            L_LAYER_CLEAR=L_LAYER_CLEAR.AND.                            &
     &         (LCA(L, GLOBAL_CLOUD_TOP) <= 0.0E+00).AND.               &
     &         ( (CCA(L) <= 0.0E+00).OR.                                &
     &           (CCT(L) <  GLOBAL_CLOUD_TOP-1) )
         ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_GLOBAL_CLOUD_TOP
#endif
#endif
