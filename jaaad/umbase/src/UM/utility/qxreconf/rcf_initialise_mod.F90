#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ initialisation for reconfiguration

Module Rcf_Initialise_Mod

!  Subroutine Rcf_Initialise - initialisation tasks
!
! Description:
! This module initialises the reconfiguration calculation,
! including Gcom, Namelists, Decompositions and File information
!
! Method:
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Initialise( hdr_in, hdr_out )

Use Rcf_Arch_Initialise_Mod, Only : &
    Rcf_Arch_Initialise

Use Rcf_Decompose_Mod, Only : &
    Rcf_Decompose

Use Rcf_Parvars_mod, Only : &
    mype,               nproc,              &
    nproc_x,            nproc_y,            &
    nproc_max,          Maxproc

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Read_Namelists_Mod, Only : &
    Rcf_Read_Namelists

Use Rcf_Stash_Init_Mod, Only : &
    Rcf_Stash_Init

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Files_Init_Mod, Only : &
    Rcf_Files_Init

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_DecompDB_Mod, Only : &
    decomp_DB

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input,    &
    decomp_rcf_output

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Min,           &
    LTimer,                 &
    L_IO_Timer

Use Rcf_HeadAddress_Mod, Only : &
    UM_Sector_Size

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_Set_Interp_Logicals_Mod, Only : &
    Rcf_Set_Interp_Logicals

Use Rcf_Recon_Mod, Only : &
    GRIB

Use Rcf_Grib_Control_Mod, Only : &
    Rcf_Grib_Control

Implicit None

! Arguments
Type (um_header_type), Intent(InOut)  :: hdr_in
Type (um_header_type), Intent(InOut)  :: hdr_out

! Local Vars/Params
Character (Len=*), Parameter    :: RoutineName = 'Rcf_Initialise'
Character (Len=80)              :: Cmessage
Character (Len=8)               :: c_nprocx
Character (Len=8)               :: c_nprocy
Character (Len=8)               :: c_printst
Character (Len=8)               :: c_timer
Character (Len=8)               :: c_um_sector_size
Integer                         :: io_timing
Integer                         :: ErrorStatus
Integer                         :: err

Integer                        :: ierr
Integer                        :: len_basename    ! length of file base
Integer                        :: len_dataw       ! length of $DATAW
Integer                        :: len_runid       ! length of $RUNID

Character (Len=180)            :: stdout_basename ! base of filename
Character (Len=170)            :: dataw_char      ! value of $DATAW
Character (Len=5)              :: runid_char      ! value of $RUNID
Character (Len=200)            :: stdout_filename ! file for stdout

External Fort_Get_Env

!----------------------------------------------------------------
! Initialise parallel variables
!----------------------------------------------------------------
! Get the x/y decomposition
Call Fort_Get_Env( 'RCF_NPROCX', 10, c_nprocx, 8, err )
If ( err /= 0 .OR. c_nprocx == 'UNSET' ) Then
  ErrorStatus = -30
  Cmessage = 'RCF_NPROCX not set: Defaults apply!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  If ( Mod( nproc, 2 ) == 0 ) Then
    nproc_x = 2
  Else
    nproc_x = 1
  End If
Else
  Read(c_nprocx,'(I4)') nproc_x
End If

Call Fort_Get_Env( 'RCF_NPROCY', 10, c_nprocy, 8, err )
If ( err /= 0 .OR. c_nprocy == 'UNSET' ) Then
  ErrorStatus = -40
  Cmessage = 'RCF_NPROCY not set: Defaults apply!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  nproc_y = nproc/nproc_x
Else
  Read(c_nprocy,'(I4)') nproc_y
End If


! Check the x/y decomp works ok
If ( nproc_x * nproc_y /= nproc ) Then
  ErrorStatus = 50
  Cmessage = 'Total number of processors does not fit EW/NS LPG'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Output the decomposition to STDOUT note that PrintStatus is not
! yet set
If (mype == 0) Then
  Write (6,*)
  Write (6,*) 'Parallel Reconfiguration using ', nproc,' processor(s)'
  Write (6,*) 'divided into a LPG with nproc_x=',nproc_x,           &
              'and nproc_y=',nproc_y
  Write (6,*)
End If

! Check on maximum size
nproc_max = nproc
If (nproc_max > Maxproc ) Then
  ErrorStatus = 45
  Cmessage = 'Using too many processors - need to edit MAXPROC in &
             & Rcf_Parvars'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!--------------------------------------------------------------------
! Set stdout unit (6) output to unique filename on every PE
!--------------------------------------------------------------------

CALL FORT_GET_ENV('RCF_STDOUT_FILE',15,stdout_basename,180, ierr)
IF (ierr .NE. 0) THEN
! Environment variable RCF_STDOUT_FILE has not been set, so we will
! construct a default stdout_basename of $DATAW/$RUNID.rcf.fort6.pe
  CALL FORT_GET_ENV('DATAW',5,dataw_char,170,ierr)
  IF (ierr .NE. 0) THEN
    Write (6,*) 'Default stdout filename not found - ',             &
                'set RCF_STDOUT_FILE?'
    WRITE(Cmessage,*) 'Failed to get value of $DATAW'
    ErrorStatus = 30
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  ENDIF

  CALL FORT_GET_ENV('RUNID',5,runid_char,5,ierr)
  IF (ierr .NE. 0) THEN
    Write (6,*) 'Default stdout filename not found - ',             &
                'set RCF_STDOUT_FILE?'
    Write (Cmessage,*) 'Failed to get value of $RUNID'
    ErrorStatus = 40
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  ENDIF

  len_dataw=Len_Trim(dataw_char)
  len_runid=Len_Trim(runid_char)
  stdout_basename = dataw_char(1:len_dataw)//'/'//               &
                    runid_char(1:len_runid)//'.rcf.fort6.pe'
ENDIF

! Now add PE number (mype) to stdout_basename to get the complete
! stdout_filename for this PE.

len_basename=Len_Trim(stdout_basename)
IF (mype .LT. 10) THEN
  WRITE(stdout_filename,'(A,I1)')  stdout_basename(1:len_basename), &
                                     mype
ELSEIF (mype .LT. 100) THEN
  WRITE(stdout_filename,'(A,I2)')  stdout_basename(1:len_basename), &
                                     mype
ELSEIF (mype .LT. 1000) THEN
  WRITE(stdout_filename,'(A,I3)')  stdout_basename(1:len_basename), &
                                     mype
ELSE
  WRITE(stdout_filename,'(A,I4)')  stdout_basename(1:len_basename), mype
ENDIF

! and close unit 6, then reopen to new filename
CLOSE(6)
OPEN(6,FILE=stdout_filename,STATUS='replace')

WRITE(6,*) nproc,' Processors initialised.'
WRITE(6,*) 'I am PE ',mype

!-----------------------------------------------------------------
! Architecture specific initialisation
!-----------------------------------------------------------------

Call Rcf_Arch_Initialise()

!----------------------------------------------------------------
! Initialise PrintStatus
!----------------------------------------------------------------
! Get level of diagnostics required
c_printst = 'UNSET'
Call Fort_Get_Env( 'RCF_PRINTSTATUS', 15, c_printst, 8, err )

If ( err /= 0 .OR. c_printst == 'UNSET' ) Then
  ErrorStatus = -60
  Cmessage = 'PrintStatus not set: Default applies!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  PrintStatus = PrStatus_Normal

Else   !  Env Var RCF_PRINTST has been set
  Read(c_printst,'(I4)') PrintStatus

End If

If (PrintStatus >= PrStatus_Min) Then
  write (6,*) 'PrintStatus is set to ',PrintStatus
End If

!----------------------------------------------------------------
! Initialise LTimer (logical for timer)
!----------------------------------------------------------------
c_timer = 'true'
Call Fort_Get_Env( 'RCF_TIMER', 9, c_timer, 8, err)

! Errors are ignored and default taken
If ( c_timer(1:4) == 'TRUE' .OR. c_timer(1:4) == 'true' ) Then
  LTimer = .TRUE.
  If (PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Timer is ON'
  End If
Else
  LTimer = .FALSE.
  If (PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Timer is OFF'
  End If
End If

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Reconfigure', 1)
! DEPENDS ON: timer
If (LTimer) Call Timer( 'Initialise', 3)

!------------------------------------------------------------------
! Initialise UM_Sector_Size from Env Var
!------------------------------------------------------------------
Call Fort_Get_Env('UM_SECTOR_SIZE',14,c_um_sector_size,8,err)
If ( err /=  0) THEN
  ErrorStatus = -70
  Write(Cmessage,*) 'UM_SECTOR_SIZE has not been set: Default = 2048'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
  um_sector_size = 2048
Else
  Read( c_um_sector_size, '(I4)') um_sector_size
End If

!------------------------------------------------------------------
! Set Namelist information
!------------------------------------------------------------------

Call Rcf_Read_Namelists( )

!-----------------------------------------------------------------
! Initialise IO Timer
!-----------------------------------------------------------------
If (L_IO_Timer) Then
  io_timing = 1
Else
  io_timing = 0
End If

Call io_timing_init( io_timing )

!-----------------------------------------------------------------
! Set STASH information
!-----------------------------------------------------------------

Call Rcf_Stash_Init()

!-----------------------------------------------------------------
! Handle GRIB data
!-----------------------------------------------------------------

If (GRIB) Then
  Call Rcf_Grib_Control( )
End If

!-----------------------------------------------------------------
! Set File Headers etc
!-----------------------------------------------------------------

Call Rcf_Files_Init( hdr_in, hdr_out )

!-----------------------------------------------------------------
! Set up Decompositions
!-----------------------------------------------------------------

Call Rcf_Decompose( Output_Grid, nproc_x, nproc_y, decomp_rcf_output )

Call Rcf_Decompose( Input_Grid, nproc_x, nproc_y, decomp_rcf_input )

Call Rcf_Change_Decomposition( decomp_rcf_input )

!-------------------------------------------------------------------
! Set interpolation logical flags
!-------------------------------------------------------------------
Call Rcf_Set_Interp_Logicals( Input_Grid, Output_Grid, Hdr_In )

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Initialise', 4)

Return
End Subroutine Rcf_Initialise
End Module Rcf_Initialise_Mod
#endif
