
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ single line

Module Rcf_Calc_3D_CCA_Mod

!  Subroutine Name1  - 1 line
!  Subroutine Name2  - 1 line
!
! Description: Calculates a 3D convective cloud amount (i.e. on model
!             levels) from the 2D convective cloud amount array
!             according to parameters specified in the umui and the
!             position of cloud base, cloud top and freezing level.
!
! Method: The 2D convective cloud amount is expanded into the vertical
!         by applying it between the cloud base and top with the
!         additional constraints that
!         (i)   If the cloud base is in the boundary layer,
!         (ii)  cloud top is above the freezing level and
!         (iii) the cloud is more than 500mb deep
!         then the cloud below the freezing level will be multiplied
!         by TOWER_FACTOR, and the cloud above the freezing level
!         will be linearly (with model level) increased to cloud top
!         where it will be equal to the 2D fraction * ANVIL_FACTOR.
!
! Current Code Owner:
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.4   16/03/00   Remove comment lines from after #include
!                                                  S. Carroll
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_3D_CCA( cca_2d, ccb, cct, t, p, cca_3d )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Recon_Mod, Only : &
    Anvil_Factor,     &
    Tower_Factor

Use Rcf_Recon_Mod, Only : &
    L_Cloud_Deep

Implicit None

! Arguments
Type( field_type ), Intent(In)    :: cca_2d
Type( field_type ), Intent(In)    :: ccb
Type( field_type ), Intent(In)    :: cct
Type( field_type ), Intent(In)    :: t
Type( field_type ), Intent(In)    :: p
Type( field_type ), Intent(InOut) :: cca_3d

! Comdecks
! TM
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

! Local variables
Real, Parameter                   :: Deep_dp = 50000 ! Critical depth
                                                     !of clouds = 500hPa
Integer                           :: Freeze_Lev( t % level_size )
Integer                           :: anvil_lev
Integer                           :: i
Integer                           :: k
Real                              :: anvil_depth
Real                              :: p_cloud_base
Real                              :: p_cloud_top
Logical                           :: deep

!---------------------------------------------------------------------
! Calculate the freezing level field
!---------------------------------------------------------------------
Do k = 1, t % levels
  Do i = 1, t % level_size
    If ( t % data(i,k) < TM ) Then
      If ( k == 1 ) Then
        Freeze_Lev( i ) = k
      Else If ( t % data(i,k-1) > TM ) Then
        Freeze_Lev( i ) = k
      End If
    End If
  End Do
End Do

!======================================================================
!  ANVIL CLOUD CALCULATION:
!  If cloud base is in the PBL, and cloud top is above (or at)
!  the freezing level, then add an anvil cloud by increasing the
!  cloud fraction linearly from freezing lev to cloud top. Also
!  decrease the cloud fraction below this level to represent the
!  'tower'.
!======================================================================

! initial zero of data
cca_3d % data(:,:) = 0.0

If (L_CLOUD_DEEP) Then

  Do I = 1, cca_3d % level_size
    ANVIL_DEPTH = 0
    ANVIL_LEV = 0
    DEEP = .FALSE.
!----------------------------------------------------------------------
! Calculate cloud depth:
!----------------------------------------------------------------------
    If (CCA_2D % data(I,1) > 0.0) Then
      P_CLOUD_BASE = p % Data( i, ccb % data_int(i,1) )
      P_CLOUD_TOP  = p % Data( i, cct % data_int(i,1) )
      DEEP   = (P_CLOUD_BASE - P_CLOUD_TOP) .GE. DEEP_DP
!----------------------------------------------------------------------
! Check to see if cloud is deep and above freezing level:
!----------------------------------------------------------------------
      If ( ( ccb % data_int(i,1) .LT. Output_Grid % bl_levels)  .AND. &
           ( cct % data_int(i,1) .GT. FREEZE_LEV(I) ) .AND.           &
           ( DEEP ) ) Then
!----------------------------------------------------------------------
! Define anvil base level as freezing level or cloud base if above FL:
!----------------------------------------------------------------------
        ANVIL_DEPTH = ( cct % data_int(i,1) - FREEZE_LEV(I) )
        ANVIL_LEV   = FREEZE_LEV(I)
        If ( ANVIL_DEPTH .GT.                          &
              (cct % data_int(i,1) - ccb % data_int(i,1)) ) Then
          ANVIL_DEPTH = cct % data_int(i,1) - ccb % data_int(i,1)
          ANVIL_LEV   = ccb % data_int(i,1)
        End If
!----------------------------------------------------------------------
! Apply wedge-shaped anvil from anvil base to cloud top:
!----------------------------------------------------------------------
        Do K = ANVIL_LEV,(cct % data_int(i,1) - 1)
          CCA_3D % data(I,K) = (ANVIL_FACTOR - TOWER_FACTOR)       &
                             * CCA_2D % data(I,1)                  &
                             * (K - ANVIL_LEV + 1)/ANVIL_DEPTH     &
                             + (CCA_2D % data(I,1) * TOWER_FACTOR)

          If (CCA_3D % data(I,K) .GE. 1.0) Then
            CCA_3D % data(I,K) = 0.99
          End If
        End Do
!----------------------------------------------------------------------
! ...and tower below (i.e. from cloud base to anvil base):
!----------------------------------------------------------------------
        Do K = ccb % data_int(i,1), ANVIL_LEV-1
          CCA_3D % data(I,K) = TOWER_FACTOR * CCA_2D % data(I,1)
        End Do
      Else
!----------------------------------------------------------------------
! If cloud is not 'deep' keep old fraction, but put on model levels:
!----------------------------------------------------------------------
        Do K = ccb % data_int(i,1), (cct % data_int(i,1) - 1)
          CCA_3D % data(I,K) = CCA_2D % data(I,1)
        End Do
      End If
!----------------------------------------------------------------------
! Finally check there is no cloud below cloud base or above cloud top!
!----------------------------------------------------------------------
      Do K = 1, (ccb % data_int(i,1)-1)
        CCA_3D % data(I,K) = 0.0
      End Do

      Do K = cct % data_int(i,1), cca_3d % levels
        CCA_3D % data(I,K) = 0.0
      End Do
    Else
      CCA_3D % data(I,:) = 0.0
    End If
  End Do
Else
  Do I = 1, cca_3d % level_size
    ANVIL_DEPTH = 0
    ANVIL_LEV = 0
!----------------------------------------------------------------------
! Calculate cloud depth:
!----------------------------------------------------------------------
    If (CCA_2D % data(I,1) > 0.0) Then
!----------------------------------------------------------------------
! Check to see if cloud is deep and above freezing level:
!----------------------------------------------------------------------
      If ( ( ccb % data_int(i,1) .LT. Output_Grid % bl_levels )  .AND.&
                 ( cct % data_int(i,1) .GT. FREEZE_LEV(I) ) ) Then
!----------------------------------------------------------------------
! Define anvil base level as freezing level or cloud base if above FL:
!----------------------------------------------------------------------
        ANVIL_DEPTH = ( cct % data_int(i,1) - FREEZE_LEV(I) )
        ANVIL_LEV   = FREEZE_LEV(I)
        If ( ANVIL_DEPTH .GT.                                     &
               (cct % data_int(i,1)-ccb % data_int(i,1)) ) Then
          ANVIL_DEPTH = cct % data_int(i,1)-ccb % data_int(i,1)
          ANVIL_LEV   = ccb % data_int(i,1)
        End If
!----------------------------------------------------------------------
! Apply wedge-shaped anvil from anvil base to cloud top:
!----------------------------------------------------------------------
        Do K = ANVIL_LEV,(cct % data_int(i,1) - 1)
          CCA_3D % data(I,K) = (ANVIL_FACTOR - TOWER_FACTOR)        &
                             * CCA_2D % data(I,1)                   &
                             * (K - ANVIL_LEV + 1)/ANVIL_DEPTH      &
                             + (CCA_2D % data(I,1) * TOWER_FACTOR)

          If (CCA_3D % data(I,K) .GE. 1.0) Then
            CCA_3D % data(I,K) = 0.99
          End If
        End Do
!----------------------------------------------------------------------
! ...and tower below (i.e. from cloud base to anvil base):
!----------------------------------------------------------------------
        Do K = ccb % data_int(i,1), ANVIL_LEV-1
          CCA_3D % data(I,K) = TOWER_FACTOR * CCA_2D % data(I,1)
        End Do
      Else
!----------------------------------------------------------------------
! If cloud does not satisfy anvil criteria, keep old fraction, but put
! on model levels:
!----------------------------------------------------------------------
        Do K = ccb % data_int(i,1), (cct % data_int(i,1) - 1)
          CCA_3D % data(I,K) = CCA_2D % data(I,1)
        End Do
      End If
!----------------------------------------------------------------------
! Finally check there is no cloud below cloud base or above cloud top!
!----------------------------------------------------------------------
      Do K = 1, (ccb % data_int(i,1)-1)
        CCA_3D % data(I,K) = 0.0
      End Do
      Do K = cct % data_int(i,1), cca_3d % levels
        CCA_3D % data(I,K) = 0.0
      End Do
    Else
      CCA_3D % data(I,:) = 0.0
    End If
  End Do
End If

Return
End Subroutine Rcf_Calc_3D_CCA
End Module Rcf_Calc_3D_CCA_Mod
