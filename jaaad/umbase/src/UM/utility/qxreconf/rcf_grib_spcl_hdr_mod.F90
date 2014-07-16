#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up lists and space to store Header info for data created not read

Module Rcf_Grib_Spcl_Hdr_Mod

! SUBROUTINE Rcf_Grib_Spcl_Hdr
!
! Description: Routind which creates lists and list members to hold
!              the header info for fields which are created rather than
!              read from the original GRIB data.
!              -NB this routine should not be used to alter the headers
!               of fields which are altered (e.g. geopotential to orog)
!
! Method:
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     21/08/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_Hdr(Lists)

Use Rcf_GRIB_Block_Params_Mod, Only :      &
  List_Marker,          Grib_Record,       &
  p_Lvl_Type,           p_Lvl_Desc_1,      &
  Tb3_Pressure,         p_Orig_cntr,       &
  p_Param_ID,           EID_Temperature,   &
  Grb_Data_Real

Use Rcf_GRIB_Lookups_Mod, Only :    &
  grib_max_fields,           grib_Exner_field,           &
  GrbOrigECMWF

Use EReport_Mod, Only :     &
    EReport

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_StashCodes_Mod, Only : &
   Stashcode_exner

Implicit None

! Global variables (#include statements etc):

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (List_Marker), Intent(InOut) :: Lists(0:grib_max_fields)

! Comdecks
#include "c_mdi.h"
#include "clookadd.h"
! contains LBLREC (amongst others)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_Hdr'

! Local variables

Type (Grib_Record),Pointer       :: New_Record,Current

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: I
Integer                          :: Count,Criteria

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Loop across all lists, finding those which will require extra fields
!=======================================================================
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) ) Then

  ! Lets only look at those ECMWF ones first
  If ( Lists(I) % Begin % Block_1 (p_Orig_cntr) == GrbOrigECMWF ) Then
    Select Case ( Lists(I) % Begin % Block_1 ( p_Param_ID ) )

      Case (EID_Temperature)
        ! Will need Exner for Height generation

        ! Set pointer to first entry in list
        Current => Lists(I) % Begin
        count =0

        ! Loop across members in list
        Do While (Associated(Current))
        count= count+1

        !Allocate a GRIB header to store info
        Write (6,*) 'About to allocate for Exner ', count
        Allocate(New_Record)

        New_Record % Block_1(:) = Current % Block_1(:)
        New_Record % Block_2(:) = Current % Block_2(:)
        New_Record % Block_3(:) = Current % Block_3(:)
        New_Record % Block_4(:) = Current % Block_4(:)
        New_Record % Block_R(:) = Current % Block_R(:)

        New_Record % Start_pos  = 0
        New_Record % Data_Type  = Grb_Data_Real
        New_Record % Num_Fp     = 0
        New_Record % Num_Vert   = 0
        New_Record % Num_Bitmap = 0
        New_Record % Num_Quasi  = 0

        ! Copy header info from Temp to Exner
        New_Record % StashCode = Stashcode_exner
        New_Record % Desc      = 'Exner field'

        New_Record % Block_1 (p_Orig_cntr) = (-1)
        New_Record % Block_1 (p_Param_ID)  = (1)

        ! Assign that record to it's list
        New_Record % Prev   => Lists(grib_Exner_field) % End
                               ! Point Prev pointer at end of
                               ! current list.(Null if first entry)
        If (Associated(Lists(grib_Exner_field) % End)) Then
                               ! If current end of list is a
                               ! valid record (Not first entry)
          Lists(grib_Exner_field) % End % Next  => New_Record
                               ! Point 'next' for previous entry
                               ! at current entry

        Else                         ! Else : must be 1st entry
          Lists(grib_Exner_field) % Begin  => New_Record
                               ! Point begining of List at New_Record
        End If

        Lists(grib_Exner_field) % End      => New_Record
                               ! Point End of List at (now complete)
                               ! New_Record Entry
        Nullify(New_Record % Next)   ! Ensure 'Next' is not associated

        Lists(grib_Exner_field) % LstCount =                          &
                                Lists(grib_Exner_field) % LstCount + 1
                               ! Add one to count of list size


      ! end do across list members
        Current => Current % Next
        End Do

    End Select

  End If
  End If

!=======================================================================
!  End Loop across all lists.
!=======================================================================
End Do  ! loop over all lists


Return

End Subroutine Rcf_Grib_Spcl_Hdr
End Module Rcf_Grib_Spcl_Hdr_Mod
#endif
