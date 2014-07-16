#if defined(C70_1A) || defined(MAKEBC)
#if !defined(SCMA) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INTF_NEW_FILES -----------------------------------------
!
!LL  Purpose: To test whether a new boundary data file needs to be
!LL           opened
!LL
!LL  Model            Modification history from model version 4.5
!LL version  Date
!LL  4.5   3/09/98    New deck added M.J.Bell
!LL  5.3  04/08/01    Allow sub-hourly output of LBCs. Peter Clark.
!    6.0  05/09/03    Add defs to allow use in makebc. R.Sempers
!    6.2  15/08/05    Free format fixes. P.Selwood
!LL
!LLEND ---------------------------------------------------------------
! ---------------------------------------------------------------------
      subroutine intf_new_files(first_unit, last_unit, max_n_intf, im,  &
     &    TYPE_LETTER_1, FT_OUTPUT, INTF_FREQ_HR, INTF_FREQ_MN,         &
     &    INTF_FREQ_SC, FT_STEPS, STEP, FT_FIRSTSTEP, INTERFACE_STEPS,  &
     &    LNEWBND )

!L  Purpose: determines new output interface files for a submodel.
!L
      implicit none
      integer first_unit ! IN first unit to test
      integer last_unit  ! IN last unit to test
      integer max_n_intf ! IN number of interface files
      integer im         ! IN  sub-model identifier
      character*1 TYPE_LETTER_1(20:last_unit) ! IN
      character*1 FT_OUTPUT(20:last_unit)     ! IN
      integer INTF_FREQ_HR(max_n_intf)     ! IN
      integer INTF_FREQ_MN(max_n_intf)     ! IN
      integer INTF_FREQ_SC(max_n_intf)     ! IN
      integer FT_STEPS(20:last_unit)          ! IN
      integer STEP                         ! IN model step no.
      integer FT_FIRSTSTEP(20:last_unit)      ! IN
      integer INTERFACE_STEPS(max_n_intf)  ! IN
      logical LNEWBND(max_n_intf)          ! OUT
!-----------------------------------------------------------------------
!L Declaration of local variables
      integer iunit
!     logical ll_intf_type   ! OUT T => file is an output interface file
      integer jintf          ! boundary file area number
      INTEGER INTF_FREQ_SECS  ! Interface frequency in seconds
!-----------------------------------------------------------------------
      do iunit = first_unit, last_unit

       if (type_letter_1(iunit) == 'b') then  !  Boundary file

         IF (FT_OUTPUT(IUNIT) == 'Y') THEN ! Intf. data output?
! DEPENDS ON: intf_area
          call intf_area ( im, iunit, JINTF)

      INTF_FREQ_SECS=3600*INTF_FREQ_HR(JINTF) +                         &
     &  60*INTF_FREQ_MN(JINTF)+INTF_FREQ_SC(JINTF)
      IF (INTF_FREQ_SECS >  0) THEN

            IF (STEP == 0) THEN

              LNEWBND(JINTF) = .TRUE. ! New intf data file required at
                                      ! first entry to ININTF1
            ELSE

              IF (FT_STEPS(IUNIT) == 0) LNEWBND(JINTF) = .FALSE. !False
                                         ! if incomplete single file

              IF (FT_STEPS(IUNIT) >  0) LNEWBND(JINTF) = .NOT.(         &
!               step = first timestep to get boundary data
     &          (STEP-FT_FIRSTSTEP(IUNIT) == 0 .OR.                     &
!               step = timestep to start new file
     &          MOD(STEP-FT_FIRSTSTEP(IUNIT),FT_STEPS(IUNIT)) /= 0)     &
     &          .AND.                                                   &
     &     STEP >  FT_FIRSTSTEP(IUNIT)-INTERFACE_STEPS(JINTF))
!                 ! False if incomplete file in sequence
            END IF  ! STEP

           ENDIF  !  INTF_FREQ_HR
!
          ELSE    !  FT_OUTPUT(IUNIT)
!
!  Possible place for setting switches for reading in interface data
!
          ENDIF   !  FT_OUTPUT
        ENDIF  !  TYPE_LETTER_1

      ENDDO   ! IUNIT

      return
      END SUBROUTINE intf_new_files
#endif
#endif
