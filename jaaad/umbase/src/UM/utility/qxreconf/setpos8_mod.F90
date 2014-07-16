#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
Module Setpos8_Mod
!+ Parallel UM version of SETPOS
!
contains
! Subroutine Interface:
Subroutine SetPos8(NFT,IPos,ICode)

!
! Description:
!  This routine provides an interface to SETPOS8 for the Parallel
!  Unified Model.
!
! Method:
!  The C SETPOS8 is renamed SETPOS8_SINGLE under #if defined(MPP).
!  This routine causes SETPOS8_SINGLE to be called by PE 0 only.
!
!
! History:
!  Model    Date     Comment and Author
!  version
!    5.4    22/05/02 Code copied and modified from SETPOS during
!                    development of GRIB reading code - R.Sharp
!
Use Rcf_Parvars_Mod

IMPLICIT NONE

! Subroutine Arguments Intent IN:

Integer, Intent(IN)  :: NFT    ! Fortran unit number
Integer, Intent(IN)  :: IPos   ! Position in file

! Subroutine Arguments Intent OUT:

Integer, Intent(OUT) :: ICode  ! Return code

! Parameters and Common blocks

!-------------------------------------------------------------------

ICode = 0

If ( mype == 0 ) Then    ! only PE 0 does any I/O
  Call SETPOS8_SINGLE(NFT,IPos,ICode)
End If

Return
End Subroutine SetPos8

End Module SetPos8_Mod
#endif
