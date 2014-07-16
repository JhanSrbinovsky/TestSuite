#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up variables before the call to the physics routines.
!
! Subroutine Interface:
      SUBROUTINE PRE_PHYSICS(                                           &
! ARGUMENTS IN
     &    row_length, rows, model_levels, wet_model_levels, qcl, qcf,   &
     &    ug, vg, f_coriolis, timestep, obs, geoforce, lcal360,         &
     &    co2start, co2rate, co2end, daycount, stepcount, ntrad1,       &
     &    ntrad, a_sw_radstep_diag, a_sw_radstep_prog,                  &
     &    l_triffid, npft,                                              &
! ARGUMENTS IN/OUT
     &    u, v, npft_trif,                                              &
! ARGUMENTS WITH INTENT OUT
     &    co2_mmr, L_rad_step, L_rad_step_prog,                         &
     &    l_rad_step_diag)


      IMPLICIT NONE
!
! Description: This routine sets several variables before the call to
!              the two physics routines.
! Method: It takes the relevant parts from the two physics routines at
!         4.5
!
!
! Current Code Owner: Zoe Gardner
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.3  15/05/01   Original Code       Z. Gardner
!  6.1  06/08/04   Correction of the geostrophic forcing
!                  option.
!  6.2  09/02/06 Replace ntrad by a_sw_radstep_prog and
!                a_sw_radstep_diag in order to implement
!                versions 3C and 3Z of the radiation code.
!                                          (J.-C. Thelen)
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!  ARGUMENTS WITH INTENT IN

! Parameters

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, npft

! Arrays

      Real                                                              &
     &  qcf(row_length, rows, wet_model_levels)                         &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, ug(row_length, rows)                                            &
     &, vg(row_length, rows)                                            &
     &, f_coriolis(row_length, rows)

! constants

      Real                                                              &
     &  co2start                                                        &
     &, co2end                                                          &
     &, co2rate                                                         &
     &, timestep

! time information

      Integer                                                           &
     &  daycount                                                        &
     &, stepcount                                                       &
     &, ntrad                                                           &
                  ! No. of timesteps between calls to radiation
     &, a_sw_radstep_diag                                               &
     &, a_sw_radstep_prog                                               &
     &, ntrad1    ! 1st timestep on which radiation called

! Logicals

      Logical                                                           &
     &  obs                                                             &
     &, geoforce                                                        &
     &, lcal360                                                         &
     &, local_time                                                      &
                    ! IN T if diagnostics required
                    !    for local time rather than GMT
     &, l_triffid

! ARGUMENTS WITH INTENT IN/OUT

      Real                                                              &
     &  u(row_length, rows, model_levels)                               &
     &, v(row_length, rows, model_levels)

      Integer                                                           &
     & npft_trif

! ARGUMENTS WITH INTENT OUT

      Real                                                              &
     &  co2_mmr

      Logical                                                           &
     &  L_rad_step                                                      &
     & ,L_rad_step_prog                                                 &
     & ,L_rad_step_diag

! local variables.

! loop counters
      Integer                                                           &
     &  i, j, k
      Real ::  utmp      ! Temporary u-wind

! External routines:

! --------------------------------------------------------------------
      If (l_triffid) then
        NPFT_TRIF = NPFT
      Else
        NPFT_TRIF = 1
      End If

!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     If Geostrophic forcing is chosen
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
      If (geoforce) then
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
!             Store current u-wind to ensure temporally consistent
!             updating of v-wind.
              utmp=u(i,j,k)
              u(i,j,k) = u(i,j,k) - f_coriolis(i,j) *                   &
     &                  (vg(i,j)-v(i,j,k))*timestep
              v(i,j,k) = v(i,j,k) + f_coriolis(i,j) *                   &
     &                  (ug(i,j)-utmp)*timestep
            enddo
          enddo
        enddo
      endif

!
!---------------------------------------------------------------------
!     Calculate the CO2 mass mixing ratio using the rate of change
!     (per year)
!---------------------------------------------------------------------
!
      If (lcal360) then
        co2_mmr = co2start + co2rate                                    &
     &    * ((daycount-1)*86400 + (stepcount-1)*timestep)               &
     &    / 360*86400
      else
        co2_mmr = co2start + co2rate                                    &
     &    * ((daycount-1)*86400 + (stepcount-1)*timestep)               &
     &    / 365*86400
      endif
      If (co2_mmr  >   co2end) co2_mmr = co2end

!---------------------------------------------------------------------
!     Is this a radiation timestep?
!---------------------------------------------------------------------

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      l_rad_step_diag=.false.
      l_rad_step_prog=.false.

      If (stepcount==ntrad1) Then
        l_rad_step_diag=.true.
        l_rad_step_prog=.true.
      Endif

      If (mod(stepcount-ntrad1, a_sw_radstep_diag)) Then
        l_rad_step_diag=.true.
      Endif
      If (mod(stepcount-ntrad1, a_sw_radstep_prog)) Then
        l_rad_step_prog=.true.
      Endif
#else
      If ((stepcount==ntrad1) .or. (mod(stepcount-ntrad1, ntrad)        &
     &                         == 0)) Then
        L_rad_step = .true.
      Else
        L_rad_step = .false.
      Endif
#endif

      Return
      END SUBROUTINE PRE_PHYSICS
#endif
