
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_Exner_at_theta

!
! Subroutine Calc_P_from_Exner
!

      Subroutine Calc_P_from_Exner(                                     &
     &                      p, kappa, p_zero,                           &
     &                      row_length, rows, model_levels,             &
     &                      off_x, off_y,                               &
     &                      exner,L_include_halos)

! Purpose:
!          Calculates pressure from exner pressure.
!
! Method:
!          p=p_zero*exner**(1/kappa)
!
!
!
! Original Programmer: Clive Wilson
! Current code owner:
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 19/06/01  5.3    New Routine        C. Wilson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y      ! Size of small halo in j.
      LOGICAL                                                           &
     &  L_include_halos  ! If .TRUE. then include halo regions
                         ! when performing the calculations

! Physical constants

      Real                                                              &
     &  kappa                                                           &
     &, p_zero


      Real                                                              &
     &  Exner(1-off_x:row_length+off_x,                                 &
     &        1-off_y:rows+off_y, model_levels)

! Arguments with Intent OUT. ie: Output variables.
       Real                                                             &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &,  i_start,i_end                                                  &
                         ! Loop bounds
     &,  j_start,j_end                                                  &
                         ! Loop bounds
     &,  vector_length   !

      Real                                                              &
     &  recip_kappa

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate p from exner.
! ----------------------------------------------------------------------

      IF (L_include_halos) THEN
        i_start=1-off_x
        i_end=row_length+off_x
        j_start=1-Off_y
        j_end=rows+Off_y
      ELSE
        i_start=1
        i_end=row_length
        j_start=1
        j_end=rows
      ENDIF

      vector_length = i_end -  i_start +1
      recip_kappa = 1./ kappa

      Do k = 1, model_levels
        Do j = j_start, j_end







          Do i = i_start, i_end
            p(i,j,k) =   p_zero * exner(i,j,k) ** recip_kappa
          End Do

        End Do
      End Do

! End of routine
      return
      END SUBROUTINE Calc_P_from_Exner

