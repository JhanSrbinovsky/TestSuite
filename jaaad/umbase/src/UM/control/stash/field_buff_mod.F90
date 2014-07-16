#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A simple data module containing variables related to STASH buffering
!
! Current Code Owner: Paul Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


      Module field_buff_mod

      Integer, Parameter :: BUFF_MAX = 161     ! Max unit number
      Integer, Parameter :: BUFF_MIN = 0       ! Min unit number

      Type pp_buffer
!Fixed Length Header
        Integer, Pointer :: fixed_header(:)
!Lookup Table
        Integer, Pointer :: pp_lookup_table(:,:)
        Integer          :: len1_pp_lookup
        Integer          :: len2_pp_lookup
        Integer          :: pp_lookup_address
      End Type pp_buffer

      ! Static buffers for PP data - 1 per permissible unit
      Type(pp_buffer) :: pp_buff( BUFF_MIN : BUFF_MAX )

      ! Dummy lookup for target on non I/O PEs
      Integer, Target :: dummy_lookup(1,1)

      Contains

!----------------------------------------------------------------------
! Subroutine Init_FXH
!
! Initialises fixed header structures
!----------------------------------------------------------------------
      Subroutine Init_fxh(ftn_unit, len_fixhd)

      Implicit None

      Integer, Intent(In) :: ftn_unit
      Integer, Intent(In) :: len_fixhd

      Allocate(pp_buff(ftn_unit) % fixed_header(len_fixhd))

      End Subroutine init_fxh



!----------------------------------------------------------------------
! Subroutine Attach_FXH
!
! Attaches a fixed header to a given unit for PP output
!----------------------------------------------------------------------
      Subroutine Attach_fxh(pp_fixhd, ftn_unit)

      Implicit None

      Integer, Intent(In) :: ftn_unit
      Integer, Pointer :: pp_fixhd(:)

      pp_fixhd => pp_buff(ftn_unit) % fixed_header

      End Subroutine attach_fxh



!----------------------------------------------------------------------
! Subroutine Init_Ipplook
!
! Initialises lookup table structures
!----------------------------------------------------------------------
      Subroutine Init_ipplook(ipplook, unit, len1_lookup,               &
     &                        pp_len2_lookup, address, step)

      Implicit None

      ! Subroutine Arguments
      Integer, Pointer    :: ipplook(:,:)
      Integer, Intent(In) :: unit
      Integer, Intent(In) :: step
      Integer, Intent(In) :: len1_lookup
      Integer, Intent(In) :: pp_len2_lookup
      Integer, Intent(In) :: address

      ! Local Variables
      Integer :: ii, jj
      Integer :: icode              ! error code
      Integer, parameter :: current_io_pe = 0

      Character (len=*), parameter :: RoutineName='Init_Ipplook'
      Character (len=80)           :: Cmessage


#include "parvars.h"

      If (mype == current_io_pe) Then
        If (step == 1) Then

!
! Check the size parameters for a size match
!
          If (pp_buff(unit) % len1_pp_lookup == len1_lookup .AND.       &
     &        pp_buff(unit) % len2_pp_lookup == pp_len2_lookup) Then

!
! The size is okay - check if the table is allocated
!
            If (.NOT. Associated(pp_buff(unit) % pp_lookup_table)) Then

! Table is not allocated - allocate it
!
              Allocate(pp_buff(unit) % pp_lookup_table(len1_lookup,     &
     &                                                pp_len2_lookup))
              pp_buff(unit) % len1_pp_lookup = len1_lookup
              pp_buff(unit) % len2_pp_lookup = pp_len2_lookup
            End If
!
! Size does not match
          Else
!
! Deallocate a current one if it exists
!
            If( Associated(pp_buff(unit) % pp_lookup_table)) Then
              Deallocate(pp_buff(unit) % pp_lookup_table)
            End If
!
! Now allocate the new buffer
!
            Allocate(pp_buff(unit) % pp_lookup_table( len1_lookup,      &
     &                                                pp_len2_lookup))
            pp_buff(unit) % len1_pp_lookup = len1_lookup
            pp_buff(unit) % len2_pp_lookup = pp_len2_lookup

          End If ! Check if the stored sizes match

          ipplook => pp_buff(unit) % pp_lookup_table

          Do ii=1, pp_len2_lookup
            Do jj=1, len1_lookup
              ipplook(jj,ii) = -99
            End Do
          End Do

        Else If (step == 2) Then
          pp_buff(unit) % pp_lookup_address = address

        Else   ! step /= 1 or 2
          Cmessage = 'init_ipplook step should be 1 or 2'
          Icode = 10
! DEPENDS ON: ereport
          Call Ereport(RoutineName, Icode, Cmessage)
        End If ! Test on step

      Else ! mype isn't the current_io_pe

        dummy_lookup(:,:) = -99
        ipplook => dummy_lookup

      End If ! mype == current_io_pe

      End Subroutine init_ipplook


!----------------------------------------------------------------------
! Subroutine Attach_Ipplook
!
! This routine attaches a buffered lookup table to a pointer
!----------------------------------------------------------------------
      Subroutine attach_ipplook(ipplook,unit)

      Implicit None

      !Subroutine Arguments
      Integer, Pointer    :: ipplook(:,:)
      Integer, Intent(In) :: unit

      ipplook => pp_buff(unit) % pp_lookup_table

      End Subroutine attach_ipplook

      End Module field_buff_mod
#endif
