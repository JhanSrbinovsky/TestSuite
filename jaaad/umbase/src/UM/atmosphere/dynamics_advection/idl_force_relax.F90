#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Force_Relax
!

      Subroutine IDL_Force_Relax(                                       &
     &             row_length, rows, model_levels                       &
     &,            halo_x, halo_y                                       &
     &,            off_x, off_y, global_row_length, global_rows         &
     &,            nproc, timestep, timestep_number                     &
     &,            max_model_levels, max_num_force_times                &
     &,            newtonian_timescale                                  &
     &,            num_force_times, force_time_interval                 &
     &,            force_data_modlev                                    &
     &,            fld, fld_incr)


! Purpose: To apply Newtonian relaxation to a field
!
! Method:  Interpolates idealised forcing data to the current time
!          Calculates the domain average of the field for each level
!          Calculates the difference between the domain avg profile and
!           the 1D forcing profile.
!          Applies Newtonian relaxation of this increment using a
!           defined relaxation timescale
!
! Original Programmer:   Richard M. Forbes
! Current code owner:    Andrew J. Malcolm
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.1    01/08/04   Original code. R. Forbes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None


! Variables with Intent (In)

      Integer                                                           &
     &  row_length                                                      &
                            ! Number of points on a row
     &, rows                                                            &
                            ! Number of rows
     &, model_levels                                                    &
                            ! Number of model levels
     &, off_x                                                           &
                            ! Size of halo in i
     &, off_y                                                           &
                            ! Size of halo in j.
     &, halo_x                                                          &
                            ! Size of halo in i for fld
     &, halo_y                                                          &
                            ! Size of halo in j for fld
     &, global_row_length                                               &
                            ! Points per global row
     &, global_rows                                                     &
                            ! No. of global (theta) rows
     &, nproc                                                           &
                            ! Number of processors
     &, num_force_times                                                 &
                            ! Number of times in forcing array
     &, timestep_number                                                 &
                            ! Timestep number in run
     &, max_model_levels                                                &
                            ! Maximum number of model levels
     &, max_num_force_times ! Maximum number of forcing times

      Real                                                              &
     &  timestep                                                        &
                               ! Timestep interval (secs)
     &, newtonian_timescale                                             &
                               ! Relaxation timescale
     &, force_time_interval    ! Interval between forcing data

      ! Forcing data on model levels
      Real                                                              &
     &  force_data_modlev(max_model_levels, max_num_force_times)

      Real                                                              &
     &  fld(1-halo_x:row_length+halo_x, 1-halo_y:rows+halo_y,           &
     &         model_levels)     ! Full field


! Variables with Intent (InOut)

      Real                                                              &
     &  fld_incr(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &         model_levels)     ! Field increment


! Local variables

      Integer                                                           &
     &  i,j,k                                                           &
                   ! Loop counters
     &, istat                                                           &
                   ! Return status from GC routine
     &, t_before                                                        &
                   ! Time index in forcing array before current time
     &, t_after    ! Time index in forcing array after current time

      Real                                                              &
     &  fld_avg                                                         &
                              ! Global average of field on a level
     &, relaxts_recip                                                   &
                              ! Reciprocal of relaxation coefficient
     &, current_time                                                    &
                              ! Time since start of run in seconds
     &, weight                                                          &
                              ! Weight for linear time interpolation
     &, force_data_modlev_int ! Time interpolated forcing data

!---------------------------------------------------------------------
! Setup time interpolation variables
!---------------------------------------------------------------------

      ! Setup time interpolation variables
      current_time = timestep*timestep_number
      t_before     = INT(current_time/force_time_interval) + 1
      t_after      = t_before + 1

      ! Set t_after appropriately if end of forcing dataset
      ! coincides with end of model run
      If (t_after  ==  num_force_times+1) t_after = t_before

      ! Set linear interpolation weights
      weight   = 1.0 - MOD(current_time,force_time_interval)            &
     &                        / force_time_interval

      ! Set reciprocal of relaxation timescale*timestep
      relaxts_recip = timestep/newtonian_timescale

!---------------------------------------------------------------------
! Loop over levels, calculate domain average and apply relaxation
!---------------------------------------------------------------------


      Do k = 1, model_levels

        !-------------------------------------------------------------
        ! Calculate horizontal average for whole domain
        ! Assumes trivial trigs. For lat/lon grid functionality in
        ! the future need to include cos_v_latitude and cos_u_latitude
        ! for area weights
        !-------------------------------------------------------------

        fld_avg = 0.0

        ! Sum field on this processor for this level
        Do j = 1, rows
          Do i = 1, row_length
            fld_avg = fld_avg + fld(i,j,k)
          End Do
        End Do


        ! Calculate sum over all processors
        !   Note: this summation is not reproducible.
        !   Extra code is required within #if defined(REPROD)

        CALL GC_RSUM(1,nproc,istat,fld_avg)


        ! Calculate average over global field
        !   Note: No spatial weighting done yet

        fld_avg = fld_avg/(global_rows*global_row_length)

        !-------------------------------------------------------------
        ! Interpolate forcing data to current time
        !-------------------------------------------------------------

        force_data_modlev_int = force_data_modlev(k,t_before)*weight    &
     &      + force_data_modlev(k,t_after)*(1-weight)

        !-------------------------------------------------------------
        ! Add relaxation term to tendency
        !-------------------------------------------------------------

        Do j = 1,rows
          Do i = 1,row_length
            fld_incr(i,j,k) = fld_incr(i,j,k) -                         &
     &           ( fld_avg - force_data_modlev_int )*relaxts_recip
          End Do
        End Do

      End Do  ! loop over levels


! End of routine.
      Return
      END SUBROUTINE IDL_Force_Relax

#endif
