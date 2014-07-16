! TYPOCONE start
      REAL:: DXT  (IMT) ! Spacing between grid points along row (T)
      REAL:: DXTR (IMT) ! Reciprocal spacing of gridpoints along row
      REAL:: DXT2R(IMT) ! Half ------------"----------------
      REAL:: DXU  (IMT) ! Spacing of U points along row
      REAL:: DXUR (IMT) ! Reciprocal of spacing of U points along row
      REAL:: DXU2R(IMT) ! Half --------------------"----------------
      REAL:: DXU4R(IMT) ! Quarter -----------------"---------------
      REAL:: DXT4R(IMT) ! Quarter reciprocal T spacing along row

      REAL:: DYT  (JMT) ! Spacing of T gridpoints N/S
      REAL:: DYTR (JMT) ! Reciprocal of ---"-----
      REAL:: DYT2R(JMT) ! Half -------"---------
      REAL:: DYU  (JMT) ! Spacing of U gridpoints N/S
      REAL:: DYUR (JMT) ! Reciprocal ------"--------
      REAL:: DYU2R(JMT) ! Half -------------------"-----

      ! DYU2R for a particular row which is outside the scope of halos
      ! in mpp code.
      REAL:: DYU2RJ

      REAL:: DYU4R(JMT) ! Quarter -----------"--------
      REAL:: DYT4R(JMT) ! Quarter reciprocal T spacing N/S

      REAL:: CS   (JMT) ! Cosine at row (to calc gridlenth) U grid
      REAL:: CSR  (JMT) ! Reciprocal ------"-------------------

      ! CSR for a particular row which is outside the scope of halos in
      !  mpp code.
      REAL:: CSRJ

      REAL:: CST  (JMT) ! Cosine at T row (to calc grid length)
      REAL:: CSTR (JMT) ! Reciprocal -------"--------------
      REAL:: PHI  (JMT) ! Grid latitude (U grid) radians
      REAL:: PHIT (JMT) ! --"--         (T grid) --"--
      REAL:: SINE (JMT) ! Sine of grid latitude (U grid)
      REAL:: TNG  (JMT) ! Tangent of grid latitude (U grid)

      REAL:: C2DZ ( KM)  ! Twice vertical grid spacing (thickness *2)
      REAL:: DZ   ( KM)  ! Grid thickness
      REAL:: DZ2R ( KM)  ! Half reciprocal grid thickness

      REAL:: EEH  ( KM)  ! Upper vert mix coeff (T)
      REAL:: EEM  ( KM)  ! -----------"-------  (U)
      REAL:: FFH  (KM)   ! Lower vert mix coeff (T)
      REAL:: FFM  (KM)   ! ---------"---------- (U)

      REAL:: ZDZ  ( KM)  ! Position of layer bottom
      REAL:: DZZ (KMP1)  ! Spacing between layer centres
      REAL:: DZZ2R(KMP1) ! Half ----------"----------
      REAL:: ZDZZ(KMP1)  ! Depth of layer centres

      REAL:: SOL_PEN(0:KM) ! Penetration of solar heat to layer base

      REAL:: DTTSA(KM)   ! Array of tracer timestep lengths
      REAL:: RAT(99)     ! Array of ratios used in constructing DTTSA
      REAL:: RZ(KM)      ! Array of level thicknesses, scaled by RAT
      REAL:: RZZ(KMP1) !Array of separations of vert levels, based on RZ
      REAL:: RZZ2R(KMP1) ! The reciprocal of twice RZZ

      REAL:: DELPSL(0:KM)  ! Solar density change (ML model)
      REAL:: DECAY(KM)     ! Decay rate of WME with depth (ML model)

      REAL:: AHI(KM) ! Along-isopycnal diffusivity (AI in (1.4)) (cm2/s)
      REAL:: AMT(JMT)    ! Latitude-dependent viscosity on tracer grid
      REAL:: AMU(JMTM1)  ! Latitude-dependent viscosity on velocity grid

      ! Depth-dependent background vertical diffusivity
      REAL:: KAPPA_B_SI(KM)

      REAL:: COSINE(IMT) ! Cosine of real latitude at gridpoint
      REAL:: RLAMBDA(IMT)! Real longitude of point

      REAL:: EDDYDIFF(JMT)!Eddy diffusion coefficients for ice-ocean h.f
      REAL:: AMX(JMT)       ! Maximum ice concentrations
      REAL:: AIN_MIN(0:NICE) ! Min sea ice concentrations for categories
      REAL:: HIN_MAX(0:NICE) ! Sea ice thickness category limits

      REAL:: BBU(JMTM1)   !
      REAL:: CCU(JMTM1)   !
      REAL:: DDU(JMTM1)   ! used in calculation of horizontal diffusion
      REAL:: GGU(JMTM1)   !
      REAL:: HHU(JMTM1)   !

      INTEGER :: KKK(99)         ! Control for output of level printouts
      REAL    :: ATHKDF(KM)
      INTEGER :: kri(2)
      REAL    :: CSRJP            ! CSR for J+2 outside halo
      REAL    :: DYU2RJP          ! DYU2R for J+2 outside halo
      REAL    :: CSTJP            ! CST for J+2 outside halo
      REAL    :: DYTRJP           ! DYTR for J+2 outside halo
      REAL    :: CSJM             ! CS for J-1 outside halo
      REAL    :: DYURJM           ! DYUR for J-1 outside halo
      REAL    :: DYUJM   ! DYU for J-1 outside halo
      REAL    :: DYT2RJM   ! DYT2R for J-1 outside halo
      REAL    :: CSTRJP   ! CSTR for J+2 outside halo
      REAL    :: CSTRJM   ! CSTR for J-1 outside halo
      REAL    :: CSTJM    ! CST  for J-1 outside halo
      INTEGER :: MAX_LARGE_LEVELS ! max no of levs for large scheme

      ! no of levels within a layer considered to calculate mld_large
      INTEGER :: NO_LAYERS_IN_LEV

#include "comocone.h"
! TYPOCONE end
