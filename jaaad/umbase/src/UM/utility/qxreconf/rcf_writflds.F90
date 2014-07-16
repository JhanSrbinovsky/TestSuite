#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This routine writes a number of fields out to dump

! Description:
!  This routine writes a number of fields out to dump.
!
! Method:
!  Setpos is used to set the position in the dump to write too.
!  The addressing is calculated and a number of variables
!  (grid code and lbc_levels) required by lower level routines are set.
!  Data is written by rcf_write_multi
!
! Documentation:
!    UMDP F3
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.4   13/06/02   Addition of optional 'single PE' flag  R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


!-------------------------------------------------------------------
! Note that this routine does not have a module as it is called
! will different type arguents for D1
!-------------------------------------------------------------------

Subroutine Rcf_WritFlds( NFTOUT,   NUMBER_OF_FIELDS,       &
                         POSITION, LOOKUP,LEN1_LOOKUP,     &
                         D1,   LEN_BUF,  FIXHD,            &
                         ErrorStatus,  CMESSAGE, Multi_PE)

Use Rcf_Parvars_Mod, Only :          &
    mype,                            &
    lasize,              blsizeu,    &
    blsizep,             blsizev

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus, &
    PrStatus_Diag

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Ppx_Info_Mod, Only : &
    STM_record_type


Use Rcf_Write_Multi_Mod, Only : &
    Rcf_Write_Multi

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Get_Fld_Type

IMPLICIT NONE

! Arguments
Integer, Intent(In)  :: nftout            ! Unit number for I/O
Integer, Intent(In)  :: number_of_fields ! No of fields to be
                                          ! written
Integer, Intent(In)  :: len_buf           ! Length of I/O buffer
Integer, Intent(In)  :: position          ! Field number from
                                          ! which to begin I/O
Integer, Intent(In)  :: fixhd(*)          ! fixed length header
Integer, Intent(In)  :: Len1_Lookup       ! 1st dim of lookup
Integer, Intent(In)  :: Lookup(Len1_Lookup,*) ! lookup table
Real,    Intent(In)  :: D1(*)             ! Data to write


Integer, Intent(Out)           :: ErrorStatus
Character(Len=80), Intent(Out) :: Cmessage

Logical, Intent(In)            :: Multi_PE

! Comdecks
#include "cppxref.h"
#include "clookadd.h"

! Local variables
Integer                :: k            ! index
Integer                :: D1_Off       ! Offset in D1 array
Integer                :: len_io
Integer                :: Field_Start  ! word address to begin I/O
Integer                :: Data_Write_Size ! data to write
Integer                :: Fld_Type
Integer                :: end_level
Integer                :: start_level
Integer                :: local_len    ! local size of field
Integer                :: grid_type
Integer                :: LBC_levels
Integer                :: field_model
Integer                :: field_sect
Integer                :: field_item
Integer                :: dummy_ppxi
Integer                :: dummy_ppxrecs
Character              :: dummy_ppxc
Character (Len=*), Parameter :: RoutineName='Writeflds'

Real                   :: A_IO
Type (STM_record_type) :: Stash_Record

! External subroutines called:---------------------------------
EXTERNAL PR_LOOK,IOERROR,SETPOS
! -------------------------------------------------------------

ErrorStatus = 0
Cmessage    = ' '

!L 2. Buffer out NUMBER_OF_FIELDS fields of data:
D1_Off = 0
DO K = POSITION, POSITION + NUMBER_OF_FIELDS-1

! Location on disk from which to begin I/O
  Field_Start=lookup(lbegin,k)

  Data_Write_Size = Lookup( lbnrec, K)

! Position file pointer
! DEPENDS ON: setpos
  Call setpos( nftout, Field_Start, ErrorStatus)


! Get some information about this field
  field_item  = MOD(LOOKUP(42,K),1000)
  field_sect  = (LOOKUP(42,K)-field_item)/1000
  field_model = LOOKUP(45,K)

  Stash_Record = Rcf_Exppx( field_model, field_sect, field_item )
  grid_type = Stash_Record % grid_type

! For atmosphere zonal ozone fields - set to zonal grid type
  If ( grid_type == ppx_atm_ozone .AND. LOOKUP(LBNPT,K) == 1) Then
       Stash_Record % grid_type = ppx_atm_tzonal
  End If

  Call Rcf_Write_Multi( nftout, d1( D1_Off + 1), Data_Write_Size,   &
                        len_io, local_len, a_io,                    &
                        LOOKUP(1,K), FIXHD(12),                     &
                        Stash_Record, Multi_PE )

! Reset the STASHmaster grid-type to what it really should be
  Stash_Record % grid_type = grid_type

! Check for I/O errors
 If (a_io /= -1.0 .OR. len_io /= Data_Write_Size ) Then
   Write(6,'('' *ERROR* Writing field no'',I5)')K
   If (FIXHD(5) < 6 .OR. FIXHD(5) > 10) Then ! Not AC/Cx/Cov/ObSt
! DEPENDS ON: pr_look
     CALL PR_LOOK( DUMMY_PPXI, DUMMY_PPXC, DUMMY_PPXRECS, &
                   LOOKUP,LOOKUP,LEN1_LOOKUP,K)
   End If

! DEPENDS ON: ioerror
   Call IOERROR('buffer out of real data', A_IO, LEN_IO,Data_Write_Size)
   ErrorStatus = NINT(A_IO)+1
   CMESSAGE    = 'Rcf_WRITFLDS:I/O error'
   Call Ereport( RoutineName, ErrorStatus, Cmessage )
 ENDIF

! Data summary used to be here - removed for time being

! Increment offset by size of level written out
  D1_Off = D1_Off + Local_Len

End Do

Return
End Subroutine Rcf_WritFlds


#endif
