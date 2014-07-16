#if defined(A05_4A) || defined(A05_5A)
      REAL C_DEEP,                                                      &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_DEEP     ! MASS FLUX FROM PARCEL BUOYANCY FOR DEEP
                      ! CONVECTION
!
      REAL C_SHALLOW,                                                   &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_SHALLOW  ! MASS FLUX FROM PARCEL BUOYANCY FOR SHALLOW
                      ! CONVECTION
!
      REAL C_MID,                                                       &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_MID      ! MASS FLUX FROM PARCEL BUOYANCY FOR MID-LEVEL
                      ! CONVECTION
!
      PARAMETER (C_DEEP = 5.17E-4, D_DEEP = 0.0)
      PARAMETER (C_SHALLOW = 5.17E-4, D_SHALLOW = 0.0)
      PARAMETER (C_MID = 5.17E-4, D_MID = 0.0)
#endif
