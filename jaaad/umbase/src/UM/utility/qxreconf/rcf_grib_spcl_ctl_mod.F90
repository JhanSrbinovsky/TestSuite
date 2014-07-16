#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Establish what information is contained within the GRIB file

Module Rcf_Grib_Spcl_Ctl_Mod

! SUBROUTINE Rcf_Grib_Spcl_Ctl
!
! Description: This routine holds the 'logic' which results in calls
!              to further routines which perform 'minor' transformations
!              on the data read from the GRIB file.
!              e.g. to transform geopotential to Orography by dividing
!                   by the  gravitational constant g.
!
! Method: The logic is split into three blocks -
!         Block 1:
!            Simple transformations which affect only a single field,
!            not requiring any other fields to exist. Thus the 'switch'
!            _will_ be the STASH code.
!            e.g. The transformation of geopotential to Orography
!
!         Block 2:
!            More complex transformations which will rely on the
!            existance of multiple fields within the dump.
!            e.g. converting T to TH requires Pstar
!            **NOTE** Any fields 'created' here (i.e. not _replacing_ a
!            particular field read in from the GRIB data) Should also
!            have a corresponding entry in Rcf_Grib_Spcl_Hdr to allocate
!            header space for the field.
!
!         Block 3:
!            Transformations which may apply to _all_ fields based on a
!            given criteria rather than the STASH code.
!            e.g. The reversing of N->S grids to get S->N grids
!            -Note- This is done last specifically so all fields used
!            are still the same way up they were read in
!
!         REMEMBER - all 3 blocks are passed through _every_ time this
!                    routine is called. (i.e. for every field read)
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/07/02   Original code. Roddy Sharp (frtz)
!  6.0     16/09/03   Remove USE of RCF_GRIB_SPCL_THETA as deck
!                     has been purged.  R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

contains
Subroutine Rcf_Grib_Spcl_Ctl(FieldData,L_FldData,Lists,Current,       &
                             Hdr_Dmy,Hdr_Itm,List_no,Lkps_Posn)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    List_Marker,     Grib_Record,     &
    p_LatGrdPnt1,    p_LatExtrmPt,    &
    p_Param_ID,      LenArrayMax,     &
    Grb_Data_Log,    Grb_Data_Int,    &
    EID_Surf_Press,  EID_Log_Surf_Press

Use rcf_GRIB_lookups_Mod, Only : &
    grib_max_fields,             &
    GrbOrigECMWF

Use Rcf_HeadAddress_Mod, Only : &
    IC_YLen,IC_Xlen

Use Rcf_UMhead_Mod, Only :  &
    um_header_type            ! Derived containing UM header info

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Diag

Use Rcf_StashCodes_Mod, Only : &
    stashcode_orog,          stashcode_exner,     &
    stashcode_lsm,           stashcode_theta,     &
    stashcode_pstar,         stashcode_soil_temp, &
    stashcode_soil_moist

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Grib_Spcl_Orog_Mod, Only : &
    Rcf_Grib_Spcl_Orog

Use Rcf_Grib_Spcl_LSM_Mod, Only : &
    Rcf_Grib_Spcl_LSM

Use Rcf_Grib_Spcl_Exner_Mod, Only : &
    Rcf_Grib_Spcl_Exner

Use Rcf_Grib_Spcl_LPstar_Mod, Only : &
    Rcf_Grib_Spcl_LPstar

Use Rcf_Grib_Spcl_SoilM_Mod, Only : &
    Rcf_Grib_Spcl_SoilM

Use EReport_Mod, Only :     &
    EReport

Use Rcf_Reverse_Field_Mod, Only : &
    Rcf_Reverse_Field

Implicit None

! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable    !Description of variable
!
! Global variables (#include statements etc):

! Subroutine arguments

!< Scalar arguments with intent(In):>
Integer, Intent(IN)                 :: Lkps_Posn
Integer, Intent(IN)                 :: List_no

!< Array  arguments with intent(In):>
Type (List_Marker), Intent(In)      :: Lists(0:grib_max_fields)
Type (Um_Header_type),Intent(In)    :: Hdr_Dmy

!< Scalar arguments with intent(InOut):>

!< Array  arguments with intent(InOut):>
Real, Intent(InOut)                 :: FieldData(LenArrayMax)
Logical, Intent(InOut)              :: L_FldData(LenArrayMax)
Type (Grib_Record),Pointer          :: Current
Type (Um_Header_type),Intent(InOut) :: Hdr_Itm

!< Scalar arguments with intent(out):>

!< Array  arguments with intent(out):>

! Comdecks
#include "c_mdi.h"
#include "clookadd.h"
! contains LBLREC (amongst others)
#include "c_g.h"
! the gravitational constant 'G'

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_Ctl'

! Local variables
Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: I
Integer                          :: Count,Criteria

Real                             :: LandPacked(LenArrayMax)
Integer                          :: LandPoints

! Local variables for saved land sea mask
! Note: although L_GotLSM is initialised as false, it acts as a saved
! variable from then on.
Logical,Save  :: LandSeaMask(LenArrayMax)
Logical,Save  :: L_GotLSM=.false.

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Block 1 : Simple transformations
!=======================================================================

Select Case (Current % Stashcode)

  Case (stashcode_orog)
    If ( mype == 0 ) Then
      Write(6,*) "Special treatment for orography"
    End If
    Call Rcf_Grib_Spcl_Orog(Current,FieldData)

  Case (stashcode_lsm)
    If ( mype == 0 ) Then
      Write(6,*) "Special treatment for landsea mask"
    End If
    Call Rcf_Grib_Spcl_LSM(Current,FieldData,L_FldData,Hdr_Itm,       &
                           Lkps_Posn)

  Case (stashcode_exner)
    If ( mype == 0 ) Then
      Write(6,*) "Special treatment for exner pressure"
    End If
    Call Rcf_Grib_Spcl_Exner(Current,FieldData)

  Case (stashcode_pstar)
    If (Current % Block_1(p_Param_ID) == EID_Log_Surf_Press) Then
      If ( mype == 0 ) Then
        Write(6,*) "Special treatment for log surface pressure"
      End If
      Call Rcf_Grib_Spcl_LPstar(Current,FieldData)
    End If

  Case (stashcode_soil_moist)
    If ( mype == 0 ) Then
      Write(6,*) "Special treatment for soil moisture"
    End If
    Call Rcf_Grib_Spcl_SoilM(Current,FieldData)

End Select

!=======================================================================
!  Block 2 : Complex transformations
!=======================================================================

!=======================================================================
!  Block 3 : Global transformations
!=======================================================================

! Check data is stored South to North not North to South
If ( Current % Block_2(p_LatGrdPnt1) >                                &
                           Current % Block_2(p_LatExtrmPt)) Then
  If ( PrintStatus >= PrStatus_Diag   ) Then
    If ( mype == 0 ) Then
      Write (6,'(2A)') "About to reverse latitudes of ", Current % Desc
    End If
  End If

  Select Case(Current % Data_Type)

  Case (Grb_Data_Log)
    Call Rcf_Reverse_Field(L_FldData,     &     ! Data to reverse
                           Hdr_Itm % IntC (IC_XLen), & ! Row Length
                           Hdr_Itm % IntC (IC_YLen), & ! No. of rows
                           1,             &     ! No. of levels in Data
                           Lkps_Posn,     &     ! Position in lookups
                           Hdr_Itm)             ! UM Hdr for updating

  Case (Grb_Data_Int)
    ! At present there is no code to flip an Integer data field
    Cmessage(1) = 'Tried to Flip Integer Data Type without any code'
    ErrorStatus = 10
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )


  Case Default
    Call Rcf_Reverse_Field(FieldData,     &     ! Data to reverse
                           Hdr_Itm % IntC (IC_XLen), & ! Row Length
                           Hdr_Itm % IntC (IC_YLen), & ! No. of rows
                           1,             &     ! No. of levels in Data
                           Lkps_Posn,     &     ! Position in lookups
                           Hdr_Itm)             ! UM Hdr for updating

  End Select
End If

! If data needs to be land packed, then do so.
! Must be done after any other field manipulations.
! - Will be dependent on stash code
Select Case (Current % Stashcode)

  Case (stashcode_lsm)
    ! Store land/sea mask in case needed for land packing below
    If ( mype == 0 ) Then
      Write(6,*) "Store Land/Sea Mask for landpacking"
    End If
    LandSeaMask=L_FldData
    L_GotLSM=.true.

  Case (stashcode_soil_temp, &
        stashcode_soil_moist)

    If (.not. L_GotLSM) Then
      Cmessage(1) = 'LSM not available for land packing'
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage(1) )
    End If

    If ( mype == 0 ) Then
      Write(6,'(A,I5)') "Performing land-packing for stashcode ", &
          Current % Stashcode
    End If

    ! Use main UM routine To_Land_Points to do the land packing.
! DEPENDS ON: to_land_points
    Call To_Land_Points(FieldData,    & ! Data to landpack
                        LandPacked,   & ! Landpacked data
                        LandSeaMask,  & ! land/sea mask
                        LenArrayMax,  & ! size of input/output arrays
                        LandPoints)     ! number of landpoints

    ! Put the land packed values at the head of the FieldData
    ! array so that the reconfiguration sees the correct values.
    FieldData(1:LandPoints)=LandPacked(1:LandPoints)

End Select

!=======================================================================
!  Cleaning Up afterwards
!=======================================================================

Return

End Subroutine Rcf_Grib_Spcl_Ctl
End Module Rcf_Grib_Spcl_Ctl_Mod
#endif
