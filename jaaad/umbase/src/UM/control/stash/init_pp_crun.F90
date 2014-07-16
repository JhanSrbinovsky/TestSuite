#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialises buffers for lookups/fixhd on a Crun

      Subroutine Init_PP_Crun( ftn_unit, env, nenv,                     &
     &                         len1_lookup, pp_len2_lookup,             &
     &                         len_fixhd, filetype)

      Use Field_Buff_mod, Only :                                        &
     &    Init_fxh,                                                     &
     &    Attach_fxh,                                                   &
     &    Init_ipplook,                                                 &
     &    Attach_ipplook

      Implicit None

!
! Description:
!   Initialises buffer space for lookups and headers for PP files
!   on a Crun.
!
! Method:
!    Simply allocates space, reads in the data and then attaches
!    data to the pointers in the PP buffer structures.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

! Subroutine arguments
      Integer :: ftn_unit
      Integer :: len1_lookup
      Integer :: pp_len2_lookup
      Integer :: len_fixhd
      Integer :: nenv

      Character*(*) :: env

      Character*1 :: filetype

! Local variables
      Integer :: icode
      Integer :: len_io
      Integer :: isize
      Integer :: address_ipplook
      Integer :: dummy
      Integer :: step
      Integer :: i

      Integer, Parameter :: current_io_pe=0

      Real :: ierr

! Pointer for the fixed length header
      Integer, Pointer :: pp_fixhd(:)

! Pointer for the lookup table
      Integer, Pointer :: ipplook(:,:)

      Character(Len=*), Parameter :: RoutineName='Init_PP_Crun'
      Character(Len=80)           :: Cmessage

#include "parvars.h"
#include "chsunits.h"
#include "cntlall.h"


! first we need to open the preattached file
! DEPENDS ON: file_open
      Call File_Open(ftn_unit, env, nenv, 1, 0, icode)
      If (icode /= 0) Then
        cmessage='Error opening file'
! DEPENDS ON: ereport
        Call Ereport(RoutineName, icode, cmessage)
      End If

! Is the file a pp file ?
      If ((filetype == 'p') .OR. (filetype == 'c')) Then

! allocate the arrays
        Allocate( pp_fixhd(len_fixhd) )
        Allocate( ipplook(len1_lookup,pp_len2_lookup) )

!  attach the fixed length header
        Call init_fxh(ftn_unit, len_fixhd)
        Call attach_fxh(pp_fixhd, ftn_unit)

! read the fixed length header
! DEPENDS ON: setpos
        Call setpos(ftn_unit, 0, icode )

! DEPENDS ON: buffin
        Call buffin(ftn_unit,pp_fixhd,len_fixhd,len_io,ierr)
        If (ierr /= -1.0) Then
          icode = -10
          cmessage = 'Error reading file'
! DEPENDS ON: ereport
          Call Ereport(RoutineName, icode, cmessage)
!     In a non-operational run, skip to the end of this routine.
          If (MODEL_STATUS /= 'Operational') THEN
             Go To 9999
          End If
        End If

! broadcast the fixed length header to all PEs
        Call gc_ibcast(77,len_fixhd,0,nproc,icode,pp_fixhd)
        If (icode /= 0) then
          Cmessage = 'Error in broadcast'
! DEPENDS ON: ereport
          Call Ereport(RoutineName, icode, cmessage)
        End If

! init the lookup table
        step=1
        Call init_ipplook(ipplook, ftn_unit, len1_lookup,               &
     &                    pp_len2_lookup,address_ipplook,step)

! the position of the lookup table is pp_fixhd(150)-1
        If (mype == 0) Then
          address_ipplook = pp_fixhd(150)-1
        End If

        step = 2
        Call init_ipplook(ipplook, ftn_unit, dummy, dummy,              &
     &                    address_ipplook, step)

! attach the lookup table
        If (mype == current_io_pe) Then
          Call attach_ipplook(ipplook,ftn_unit)
        End If

! read the lookup table
! DEPENDS ON: setpos
        Call setpos(ftn_unit, address_ipplook, icode )
        If (icode /= 0) Then
          Cmessage = 'Error in setpos'
! DEPENDS ON: ereport
          Call Ereport(RoutineName, icode, cmessage)
        End If

! Only current_io_pe should do the read as only this PE has 
! memory allocated for this operation!
        If (mype == current_io_pe) Then
          Call buffin_single(ftn_unit, ipplook,                         &
     &                len1_lookup*pp_len2_lookup, len_io, ierr )
          If (ierr /= -1.0) Then
            icode = -20
            Cmessage = 'Error reading lookup table'
! DEPENDS ON: ereport
            Call Ereport(RoutineName, icode, cmessage)
!     In a non-operational run, skip to the end of this routine.
            If (MODEL_STATUS /= 'Operational') THEN
               Go To 9999
            End If
          End If
        End If ! current_io_pe

! Nullify ipplook and pp_fixhd for safety
        Nullify(ipplook)
        Nullify(pp_fixhd)

      End If ! is_ppfile


9999  Continue
      End Subroutine init_pp_crun
#endif
