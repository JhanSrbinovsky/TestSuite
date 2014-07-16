#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine COPY_FIELD

      SUBROUTINE COPY_FIELD(                                            &
     &         field_in, field_out                                      &
     &,        row_length_in, row_length_out, rows_in, rows_out         &
     &,        levels_in, levels_out, level_start, level_end            &
     &,        haloi_in, haloj_in, haloi_out, haloj_out                 &
     &,        fld_type, L_haloes, L_swap, L_vector)

! Purpose:
!     This routine copies one field into another, allowing for a
!     different halo size in the two fields.
!    If L_haloes is true then haloes are copied explicitly BUT NB
!    NB input haloes must be greater or equal to output haloes
!    If L_haloes is falsee then only non-halo regions copied
!    If L_SWAP is true the halos on the destination field will be
!    updated using the standard swap_bounds
!
! Method:
!          Is described in ;
!
! Original Progammer: T. Davies
!
! History:
! Version    Date     Comment
! ----     -------     -------
! 5.3       10/10/01    This deck created.      Terry Davies
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &   row_length_in                                                  &
     &,  row_length_out                                                 &
     &,  rows_in                                                        &
     &,  rows_out                                                       &
     &,  levels_in                                                      &
     &,  levels_out                                                     &
     &,  haloi_in                                                       &
     &,  haloj_in                                                       &
     &,  haloi_out                                                      &
     &,  haloj_out                                                      &
     &,  level_start                                                    &
     &,  level_end                                                      &
     &,  fld_type     ! Data type for swap bounds

      Logical                                                           &
     &  L_haloes                                                        &
                   ! If true fill haloes explicitly
     &, L_swap                                                          &
                   ! If true use swap bounds to fill haloes
     &, L_vector   ! Vector switch for swap bounds

! Primary Arrays
      Real                                                              &
     &  field_in(1-haloi_in: row_length_in+haloi_in                     &
     &,          1-haloj_in: rows_in+haloj_in, levels_in)               &
     &, field_out(1-haloi_out: row_length_out+haloi_out                 &
     &,           1-haloj_out: rows_out+haloj_out, levels_out)

! Local Variables.
! scalars
      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, levels

! No External Routines:

! Functions: None
! ----------------------------------------------------------------------

      if (L_haloes) then
! ----------------------------------------------------------------------
! Section 1.  If L_haloe = .true. then explicitly copy halo region
!             In this case halo_in must be >= halo_out
! ----------------------------------------------------------------------

        if (       haloi_out  >   haloi_in                              &
     &       .or.   haloj_out  >   haloj_in) then
          print*,' Copying fields using SUBROUTINE COPY_FIELD '
          print*,' You must have output halo <= input halo '
          STOP
        endif

        Do k = level_start, level_end
          Do j = 1-haloj_out, rows_out+haloj_out
            Do i = 1-haloi_out, row_length_out+haloi_out
              field_out(i,j,k) = field_in(i,j,k)
            End Do
          End Do
        End Do

      else
! ----------------------------------------------------------------------
! Section 2.  If L_ haloes= .false. then only copy non-halo region
!              To fill haloes set L_swap = .true to use swap bounds
! ----------------------------------------------------------------------
        Do k = level_start, level_end
          Do j = 1, rows_out
            Do i = 1, row_length_out
              field_out(i,j,k) = field_in(i,j,k)
            End Do
          End Do
        End Do

        if (L_swap) then

          levels = level_end - level_start + 1

! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   field_out(1-haloi_out,1-haloj_out,level_start) &
     &,                  row_length_out, rows_out, levels               &
     &,                  haloi_out, haloj_out, fld_type, L_vector)
        endif !  L_swap

      endif ! L_haloes

      return
      END SUBROUTINE COPY_FIELD

#endif
#endif
