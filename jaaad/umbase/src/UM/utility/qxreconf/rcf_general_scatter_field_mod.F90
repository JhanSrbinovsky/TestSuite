#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Scatters any type of field from one processor to many processors

Module Rcf_General_Scatter_Field_Mod

!  Subroutine Rcf_General_Scatter_Field
!
! Description:
!   Chooses scattering method for a field based on grid_code.
!
! Method:
!   Land compressed fields are uncompressed, scattered and recompressed
!   LBCs and zonal fields have specialist subroutines
!   Other fields are "normal" and can be dealt with by Rcf_Scatter_Field
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   06/01/03   River routing support. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.



Contains

SUBROUTINE Rcf_General_Scatter_Field ( LOCAL_FIELD,  GLOBAL_FIELD , &
                                       LOCAL_SIZE,   GLOBAL_SIZE ,  &
                                       STASH_RECORD, SCATTER_PE )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only :   &
    Local_Land_Field,     &
    Local_Atmos_Landmask, &
    Glob_Atmos_Landmask

Use Rcf_Parvars_Mod, Only :                    &
    blsizep,                blsizeu,           &
    blsizev,                glsizeu,           &
    glsizev,                glsize,            &
    blsizer,                glsizer,           &
    lasize,                 gc_all_proc_group, &
    fld_type_p,             fld_type_u,        &
    fld_type_v,             fld_type_r,        &
    mype

Use Rcf_Scatter_Atmos_LBCs_Mod, Only : &
    Rcf_Scatter_Atmos_LBCs

Use Rcf_Scatter_Ocean_LBCs_Mod, Only : &
    Rcf_Scatter_Ocean_LBCs

Use Rcf_Scatter_Zonal_Field_Mod, Only : &
    Rcf_Scatter_Zonal_Field

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Get_Fld_Type

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: Global_Size  ! size of global_field
Integer, Intent(In)   :: Scatter_PE   ! PE on which global field lives
Integer, Intent(Out)  :: Local_Size   ! size of local_field

Real, Intent(In)      :: Global_Field( Global_Size ) ! field to scatter
Real, Intent(Out)     :: Local_Field( * )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Comdecks
#include "cppxref.h"

! Local variables
Integer               :: dummy       ! ignored argument
Integer               :: fld_type    ! P or U or V
Integer               :: local_x     ! local x size
Integer               :: local_y     ! local y size
Integer               :: global_x    ! global x size
Integer               :: global_y    ! global y size
Integer               :: ErrorStatus ! error code
Character (Len=*), Parameter   :: RoutineName='Rcf_General_Scatter_&
                                              &Field'
Character (Len=80)             :: Cmessage

Real                  :: buf_expand( glsize(1) * glsize(2) )
Real                  :: buf_expand_local( lasize(1) * lasize(2) )

!===================================================================

! Choose action depending on grid code
Select Case( Stash_Record % grid_type )

!-------------------------------------------------------------------
! Land compressed fields
!-------------------------------------------------------------------
  Case( ppx_atm_compressed )

    IF (mype .EQ. SCATTER_PE) THEN

! DEPENDS ON: from_land_points
      CALL FROM_LAND_POINTS(buf_expand, GLOBAL_FIELD, &
                            glob_atmos_landmask,      &
                            glsize(1)*glsize(2),dummy)

! SCATTER_PE now contains the expanded version of the full field

    ENDIF

! Now scatter this to all the other processors, putting the local
! part of the field into the array buf_expand_local

    CALL Rcf_Scatter_Field(buf_expand_local , buf_expand,       &
                           lasize(1) ,        lasize(2),        &
                           glsize(1) ,        glsize(2),        &
                           SCATTER_PE,        GC_ALL_PROC_GROUP )

! Pack the local field down to local land points and put
! the packed field into LOCAL_FIELD

! DEPENDS ON: to_land_points
    CALL TO_LAND_POINTS(buf_expand_local,LOCAL_FIELD, &
                        local_atmos_landmask,         &
                        lasize(1)*lasize(2), dummy)

    local_size = local_land_field

!-------------------------------------------------------------------
! Atmosphere Lateral boundary fields
!-------------------------------------------------------------------
  Case ( ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v )

    CALL Rcf_Scatter_atmos_LBCS( Global_Field, Global_Size, &
                                 Local_Field,  Local_Size,  &
                                 Stash_Record, Scatter_PE)


!-------------------------------------------------------------------
! Ocean Lateral boundary fields
!-------------------------------------------------------------------
  Case (ppx_ocn_rim)

    CALL Rcf_Scatter_Ocean_LBCS( Global_Field, Global_Size, &
                                 Local_Field,  Local_Size,  &
                                 Stash_Record, Scatter_PE)

!-------------------------------------------------------------------
! Zonal fields
!-------------------------------------------------------------------
  Case ( ppx_atm_tzonal, ppx_atm_uzonal, &
         ppx_ocn_tzonal, ppx_ocn_uzonal )

    If ( Stash_Record % grid_type == ppx_atm_tzonal .OR. &
         Stash_Record % grid_type == ppx_ocn_tzonal ) Then    ! P grid
      global_y = glsize(2)
      local_y  = blsizep(2)

    Else                ! U grid
      global_y = glsizeu(2)
      local_y  = blsizeu(2)
    End If

    local_size = local_y

    Call Rcf_Scatter_Zonal_Field( Local_Field, Global_Field,           &
                                  local_y,     global_y,               &
                                  1,          Stash_Record % grid_type,&
                                  Scatter_PE )

!-------------------------------------------------------------------
! Normal fields
!-------------------------------------------------------------------
  Case &
!     atmosphere grids
       ( ppx_atm_tall,   &! Atmos T points
         ppx_atm_tland,  &! Atmos T land points
         ppx_atm_tsea,   &! Atmos T sea points
         ppx_atm_uall,   &! Atmos U points
         ppx_atm_uland,  &! Atmos U land points
         ppx_atm_usea,   &! Atmos U sea points
         ppx_atm_cuall,  &! Atmos C grid U pts
         ppx_atm_cvall,  &! Atmos C grid V pts
         ppx_atm_ozone,  &! Atmos ozone field
         ppx_atm_river,  &! Atmos river routing field
!     ocean grids
         ppx_ocn_tcomp,  &! Ocean "Compressed" T
         ppx_ocn_ucomp,  &! Ocean "Compressed" u
         ppx_ocn_tall,   &! Ocean T points (cyc)
         ppx_ocn_uall,   &! Ocean U points (cyc)
         ppx_ocn_cuall,  &! Ocean C grid U pts
         ppx_ocn_cvall,  &! Ocean C grid V pts
         ppx_ocn_tfield, &! Ocean T points
         ppx_ocn_ufield)! Ocean U points


    fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

    IF (fld_type .EQ. fld_type_p) THEN
      global_x = glsize(1)
      global_y = glsize(2)
      local_x  = blsizep(1)
      local_y  = blsizep(2)
    ELSE IF (fld_type .EQ. fld_type_u) THEN
      global_x = glsizeu(1)
      global_y = glsizeu(2)
      local_x  = blsizeu(1)
      local_y  = blsizeu(2)
    ELSE IF (fld_type .EQ. fld_type_v) THEN
      global_x = glsizev(1)
      global_y = glsizev(2)
      local_x  = blsizev(1)
      local_y  = blsizev(2)
    ELSE IF (fld_type .EQ. fld_type_r) THEN
      global_x = glsizer(1)
      global_y = glsizer(2)
      local_x  = blsizer(1)
      local_y  = blsizer(2)
    ENDIF

    local_size = local_x * local_y

    Call Rcf_Scatter_Field( Local_Field, Global_Field,     &
                            local_x,     local_y,          &
                            global_x,    global_y,         &
                            Scatter_PE,  GC_ALL_PROC_GROUP )

!-------------------------------------------------------------------
! Any other type of field
!-------------------------------------------------------------------
  Case Default

    ErrorStatus = 10
    Cmessage='Field type not recognized for Scatter'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

RETURN
END Subroutine Rcf_General_Scatter_Field
End Module Rcf_General_Scatter_Field_Mod
#endif
