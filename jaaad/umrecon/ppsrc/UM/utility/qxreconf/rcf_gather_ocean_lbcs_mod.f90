
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers ocean LBCs from all PEs

Module Rcf_Gather_Ocean_Lbcs_Mod

!  Subroutine Rcf_Gather_Ocean_LBCs
!
! Description:
!   Gathers a distributed LBC field onto a single PE.
!
! Method:
!   Assumes that field is split into equal chunks on each PE.
!   Assumes input and output fields are same size.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Gather_Ocean_Lbcs(Full_LBC,     Full_LBC_Size,   &
                                 Part_LBC,     Part_LBC_Size,   &
                                 Stash_Record, Gather_pe)

Use Rcf_Parvars_Mod, Only : &
    mype,               &
    nproc,              &
    glsize

Use Rcf_Address_Length_Mod, Only : &
    Rcf_Address_Length

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

Implicit None

! Arguments
Integer, Intent(In)                 :: Full_LBC_Size
Integer, Intent(Out)                :: Part_LBC_Size
Integer, Intent(In)                 :: Gather_pe

Real, Intent(Out)                   :: Full_LBC( Full_LBC_Size )
Real, Intent(In)                    :: Part_LBC( * )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Comdecks
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! GC - General Communication primitives package. For use on
! multiprocessor shared memory and message passing systems.
!
!
! LICENSING TERMS
!
!  GC is provided free of charge. Unless otherwise agreed with SINTEF,
!  use and redistribution in source and binary forms are permitted
!  provided that
!
!      (1) source distributions retain all comments appearing within
!          this file header, and
!
!      (2) distributions including binaries display the following
!          acknowledgement:
!
!              "This product includes software developed by SINTEF.",
!
!          in the documentation or other materials provided with the
!          distribution and in all advertising materials mentioning
!          features or use of this software.
!
!  The name of SINTEF may not be used to endorse or promote products
!  derived from this software without specific prior written
!  permission.  SINTEF disclaims any warranty that this software will
!  be fit for any specific purposes. In no event shall SINTEF be liable
!  for any loss of performance or for indirect or consequential damage
!  or direct or indirect injury of any kind. In no case shall SINTEF
!  be liable for any representation or warranty make to any third party
!  by the users of this software.
!
!
! Fortran header file. PLEASE use the parameter variables in user
! routines calling GC and NOT the numeric values. The latter are
! subject to change without further notice.
!
!---------------------------------------------- ------------------------
! $Id: gps0h501,v 1.6 2000/04/17 10:05:47 t11ps Exp $
! (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

!    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
!                    value to be set from a shell variable.
!                      Author: Bob Carruthers  Cray Research.
!    5.1   17/04/00  Fixed/Free format. P.Selwood.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     GC general options
      INTEGER, PARAMETER :: GC_OK         =     0
      INTEGER, PARAMETER :: GC_FAIL       =    -1
      INTEGER, PARAMETER :: GC_NONE       =     0
      INTEGER, PARAMETER :: GC_ANY        =    -1
      INTEGER, PARAMETER :: GC_DONTCARE   =    -1
      INTEGER, PARAMETER :: GC_SHM_DIR    =     1
      INTEGER, PARAMETER :: GC_SHM_SAFE   =     2
      INTEGER, PARAMETER :: GC_NAM_TIMEOUT=     4
      INTEGER, PARAMETER :: GC_SHM_GET    = -9999
      INTEGER, PARAMETER :: GC_SHM_PUT    = -9998
      INTEGER, PARAMETER :: GC_USE_GET    = -9999
      INTEGER, PARAMETER :: GC_USE_PUT    = -9998

!     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

!     GC groups (GCG) support
      INTEGER, PARAMETER :: GC_ALLGROUP = 0
      INTEGER, PARAMETER :: GCG_ALL = GC_ALLGROUP

!     GC groups (GCG) functions
      INTEGER GCG_ME

!     GC reserved message tags
      INTEGER, PARAMETER :: GC_MTAG_LOW   = 999999901
      INTEGER, PARAMETER :: GC_MTAG_HIGH  = 999999999

!     GCG_RALLETOALLE index parameters
      INTEGER, PARAMETER :: S_DESTINATION_PE = 1
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_SEND_ARRAY = 2
      INTEGER, PARAMETER :: S_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: S_STRIDE_IN_SEND_ARRAY = 4
      INTEGER, PARAMETER :: S_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_RECV_ARRAY = 6
      INTEGER, PARAMETER :: S_STRIDE_IN_RECV_ARRAY = 7

      INTEGER, PARAMETER :: R_SOURCE_PE = 1
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_RECV_ARRAY = 2
      INTEGER, PARAMETER :: R_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: R_STRIDE_IN_RECV_ARRAY = 4
      INTEGER, PARAMETER :: R_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_SEND_ARRAY = 6
      INTEGER, PARAMETER :: R_STRIDE_IN_SEND_ARRAY = 7

! Local variables
Integer                      :: i
Integer                      :: global_size
Integer                      :: local_size
Integer                      :: local_extra
Integer                      :: part_size
Integer                      :: istat
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName='Rcf_Gather_Ocean_LBCs'
Character (Len=80)           :: Cmessage


! Global size of LBC - found from Rcf_Address_Length( input and
! output grids same size )

Call Rcf_Address_Length( Stash_Record % Grid_Type,             &
                         Stash_Record % Halo_Type, global_size )

local_size = global_size / nproc
local_extra = Mod( global_size, nproc )

! Use gc_rsend and gc_rrecv - could use gcg_ralltoall but this is
! simpler to work out

If (mype /= nproc - 1) Then
  Part_size = local_size
Else
  Part_size = local_size + local_extra
End If

! Now send to Processor gather_pe on all pes
Call Gc_rsend( 200, Part_size, gather_pe, istat, &
               Full_LBC( mype * local_size + 1), Part_LBC )
If (istat /= GC_OK) Then
  ErrorStatus = 30
  Cmessage = 'Failure in Gcom - gc_rsend'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Gc_Ssync( nproc, istat)

! Now recv on processor gather_pe
If (mype == gather_pe) Then
  Do i = 0, nproc - 2
    Call Gc_rrecv( 200, local_size, i, istat, &
                   Full_LBC( i * local_size + 1 ), Part_LBC )
    If (istat /= GC_OK) Then
      ErrorStatus = 10
      Cmessage = 'Failure in Gcom - gc_rrecv'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End Do

  ! Treat nproc - 1 seperately
  i = nproc - 1
  Call Gc_rrecv( 200, local_size + local_extra, i, istat, &
                 Full_LBC( i * local_size + 1 ), Part_LBC )
  If (istat /= GC_OK) Then
    ErrorStatus = 20
    Cmessage = 'Failure in Gcom - gc_rrecv'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If

Call Gc_Ssync( nproc, istat)

! Set the return size
part_LBC_size = local_size
If (mype == nproc - 1) part_LBC_size = part_LBC_size + local_extra

Return
End Subroutine Rcf_Gather_Ocean_Lbcs

End Module Rcf_Gather_Ocean_Lbcs_Mod


