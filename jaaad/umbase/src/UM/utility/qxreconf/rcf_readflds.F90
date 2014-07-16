#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This routine reads in a number of fields from the dump.

! Description:
!   reads a number of fields from a dump into the D1 array.
!
! Method:
!   Based on UM4.5 code.
!   UMDP F3 is the relevant documentation.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.1   15/11/00   Allow packed fieldsfiles for VAR. P.Selwood.
!   5.3   24/08/01   Improve error handling of ill-formed files.
!                    P.Selwood
!   6.1   01/08/04   Add default case for 'unkown gridtype' R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! No module for this routine as it is called with real/logical/integer
! arguments and this could cause problems. Also called from old
! F77 code and this would also cause problems here.

Subroutine Rcf_ReadFlds(NFTIN, NUMBER_OF_FIELDS,        &
                        POSITION, LOOKUP, LEN1_LOOKUP,  &
                        D1, LEN_BUF, FIXHD,             &
                        ICODE, CMESSAGE)

Use Rcf_Parvars_Mod, Only : &
    fld_type_p,             &
    fld_type_u,             &
    fld_type_v,             &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Diag

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Ppx_Info_Mod, Only : &
    STM_record_type

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Read_Multi_Mod, Only : &
    Rcf_Read_Multi

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Get_Fld_Type

Use Rcf_Recon_Mod, Only : &
    Var_Recon

IMPLICIT NONE

! Arguments
Integer, Intent(In)  :: nftin                 ! Unit no. for I/O
Integer, Intent(In)  :: number_of_fields      ! No of fields to read
Integer, Intent(In)  :: len_buf               ! Length of I/O buffer
Integer, Intent(In)  :: position              ! Field number from which
                                              ! to begin I/O
Integer, Intent(In)  :: fixhd(*)              ! Fixed lenght header
Integer, Intent(In)  :: len1_lookup           ! 1st dimension of lookup
Integer, Intent(In)  :: lookup(len1_lookup,*) ! PP Lookup table
Real,    Intent(Out) :: D1(*)                 ! data space to be filled

Integer, Intent(Out)            :: icode      ! error code
Character (Len=80), Intent(Out) :: Cmessage   ! Error Message

! Comdecks
#include "c_mdi.h"
#include "cppxref.h"
#include "clookadd.h"

! Local variables:---------------------------------------------
Integer, Parameter    :: unset = -1     ! flag
Integer   :: D1_Off         ! Offset in D1 array
Integer   :: k              ! index
Integer   :: len_io         ! Length of I/O returned by LENGTH
Integer   :: Data_Size      ! No of words of data on disk (inc WFIO pad)
Integer   :: Field_Start    ! word address to begin I/O
Integer   :: local_len      ! length of local part of field read in
Integer   :: Data_Read_Size ! No of words to read from disk
Integer   :: start_level    !  Level delimiter for  LBC size
Integer   :: end_level      !   "        "      "    "   "
Integer   :: lbc_levels     ! number of levels for LBC code.
Integer   :: field_model    !  Model code for field
Integer   :: field_sect     !  Section code for field
Integer   :: field_item     !  Item code for field
Integer   :: grid_type      !  Grid code from stash
Integer   :: ErrorStatus
Integer   :: Fld_Type
Integer   :: dummy_ppxi     ! Dummy values
Integer   :: dummy_ppxrecs  ! for old style subroutine call
Integer   :: expand_size    ! expanded data size

Real      :: A_IO

Character :: dummy_ppxc     ! dummy
Character (Len=*), Parameter  :: RoutineName = 'ReadFlds'

Type (STM_record_type)        :: Stash_Record

EXTERNAL PR_LOOK, IOERROR, SETPOS
! -------------------------------------------------------------

icode = 0
Cmessage = ' '

!  Buffer in NUMBER_OF_FIELDS fields of real data:
D1_Off = 0
DO  K=POSITION,POSITION+NUMBER_OF_FIELDS-1

! Location on disk from which to begin I/O
  Field_Start=lookup(lbegin,k)

! data_size contains the number of words of data used to store the
! field on disk
  If (mod(lookup(lbpack,k),10) == 2) Then
    Data_Size = (lookup(lblrec,k)+1)/2    ! 32 bit packed field
  Else
    Data_Size = lookup(lblrec,k)
  End If

! data_read_size contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added
  If (K /= (Position + Number_Of_Fields - 1) ) Then
    Data_Read_Size = LOOKUP(lbnrec,K)
  Else
    Data_Read_Size = data_size
  End If

! We will only deal with well-formed files. Thus an error is required
! if an old dump is encountered
  If ( LOOKUP(lblrec,K) /= 0 .AND.                             &
      (LOOKUP(lbnrec,K) == 0 .OR. LOOKUP(lbnrec,K) == IMDI) ) Then
    ErrorStatus = 10
    Cmessage = 'Invalid dump addressing: Possibly data file is &
               &ill-formed'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  ! Position file pointer
! DEPENDS ON: setpos
  Call Setpos(nftin, Field_Start, icode)

  ! Get some information about this field
  field_item  = MOD(LOOKUP(42,K),1000)
  field_sect  = (LOOKUP(42,K)-field_item)/1000
  field_model = LOOKUP(45,K)

  Stash_Record = Rcf_Exppx( field_model, field_sect, field_item )
  grid_type = Stash_Record % grid_type

  fld_type = Rcf_Get_Fld_Type( Stash_Record % grid_type)

  ! Special case for d1_grid_type for ancillary fields
  If ( FixHd(5) == 4) Then     ! Ancil data
    If ( (Mod( Lookup(lbpack,k)/10, 10 ) == 0 )  .AND. &
          grid_type == ppx_atm_compressed ) Then

      ! Compressed in stashmaster but uncompressed in header
      Select Case( fld_type )
        Case ( fld_type_p )
          Stash_Record % grid_type = ppx_atm_tall

        Case( fld_type_u )
          Stash_Record % grid_type = ppx_atm_cuall

        Case( fld_type_v )
          Stash_Record % grid_type = ppx_atm_cvall

        Case Default
          ErrorStatus = 20
          Cmessage = 'Unknown grid type for Ancil data field.'
          Call Ereport( RoutineName, ErrorStatus, Cmessage )

      End Select
    End If
  End If

  ! For atmosphere zonal ozone fields - set to zonal grid type
  IF (grid_type .EQ. ppx_atm_ozone .AND.  LOOKUP(LBNPT,K).EQ.1) THEN
    Stash_Record % grid_type = ppx_atm_tzonal
  ENDIF

! Set the size of the expanded data buffer
! Note that for compressed FieldsFiles the exact unpacked size isn't
! known and we need to guess a maximum for the memory needed.
  If ( Mod (Lookup(lbpack, k), 10) == 1) Then
    expand_size = Max( 2 * Data_Read_Size,                      &
                     2 * Lookup(lbrow, k) * Lookup(lbnpt, k) )
  Else
    expand_size = 2 * Data_Read_Size
  End If

  Call Rcf_Read_Multi(nftin, D1(D1_Off+1), Data_Read_Size,    &
                      expand_size, len_io, local_len, a_io,   &
                      LOOKUP(1,K), FIXHD(12), Stash_Record )

! Change the grid-type in the STASHmaster back to what it should
! be...
  Stash_Record % grid_type = grid_type

! Check for I/O errors
  If (a_io .ne. -1.0 .or. len_io .ne. Data_Read_Size) then
    Write(6,'('' *ERROR* Reading field no'',I5)')K
    If (FIXHD(5).LT.6 .OR. FIXHD(5).GT.10) THEN ! Not AC/Cx/Cov/ObSt
! DEPENDS ON: pr_look
      Call PR_Look( DUMMY_PPXI, DUMMY_PPXC, DUMMY_PPXRECS, &
                    LOOKUP,LOOKUP,LEN1_LOOKUP,K)
    End If

! DEPENDS ON: ioerror
    Call IOERROR('buffer in of real data',A_IO,LEN_IO, Data_Read_Size)

    ICODE=NINT(A_IO)+1
    Cmessage = 'Rcf_READFLDS:I/O error'
    Call Ereport( RoutineName, Icode, Cmessage )
  End If

  ! Data summary used to be here - removed for the time being

  ! increment size of level read in
  D1_Off = D1_Off + local_len

End Do

RETURN
End Subroutine Rcf_Readflds
#endif
