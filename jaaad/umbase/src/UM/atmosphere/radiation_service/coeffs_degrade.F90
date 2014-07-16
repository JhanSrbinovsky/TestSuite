#if defined(A70_1B)||defined(A70_1C)||defined(A70_1Z)
!+ Calculate the interpolation coefficients for spatial degradation
!-
!  Purpose: To calculate the coefficients used to interpolate the
!           radiation quantities from neighbouring grid boxes to
!           a grid box where the radiation calculations were not
!           done on a radiation timestep.
!
!  Method:  In a chequer-board pattern, a box which has not done a
!           radiation calculation is surrounded on the north, south,
!           east and west by boxes which have calculated the radiation
!           quantities therefore a coefficient of 0.25 for these four
!           points would seem to be sufficient. However, the land-sea
!           contrast can sometimes be so sharp that we would not wish
!           to interpolate from one grid box to another if their
!           surface types do not match. Furthermore, at the boundary
!           of the whole model domain, halo points do not contain
!           meaningful values, so we do not use them.
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   16/08/01   Original code.  S. Cusack.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
!
!
!
      SUBROUTINE COEFFS_DEGRADE(es_space_interp, land, first_row,       &
     &                          last_row, model_domain, mt_lam,         &
     &                          at_top_of_LPG, at_base_of_LPG,          &
     &                          at_left_of_LPG,at_right_of_LPG,         &
     &                          p_field, row_length, rows, offx, offy)
!
      IMPLICIT NONE
!
      Integer, Intent(IN) :: p_field
!                              Full length of data on a PE
      Integer, Intent(IN) :: row_length
!                              No of points per row
      Integer, Intent(IN) :: first_row
!                              first row of non-polar data
      Integer, Intent(IN) :: last_row
!                              last row of non-polar data
      Integer, Intent(IN) :: model_domain
!
      Integer, Intent(IN) :: mt_lam
!
      Integer, Intent(IN) :: offx
!
      Integer, Intent(IN) :: offy
!
      Integer, Intent(IN) :: rows
!
      Logical, Intent(IN) :: land(1-offx:row_length+offx,               &
     &                                             1-offy:rows+offy)
!            If true, then grid box is a land point
      Logical, Intent(IN) ::  at_top_of_LPG
!              Logical variable indicating if this PE at edge of LPG
      Logical, Intent(IN) ::  at_right_of_LPG
!              Logical variable indicating if this PE at edge of LPG
      Logical, Intent(IN) ::  at_base_of_LPG
!              Logical variable indicating if this PE at edge of LPG
      Logical, Intent(IN) ::  at_left_of_LPG
!              Logical variable indicating if this PE at edge of LPG

      Real, Intent(OUT) :: es_space_interp(4, row_length, rows)
!                The coefficients for radiation interpolation
!
!     Local variables
!
      Integer :: i              ! loop variable
      Integer :: j              ! loop variable
      Integer :: k              ! loop variable
      Integer :: row_length_tot ! Total row_length including halos

      Logical :: Surf_match(4, row_length, rows)
!                   True if central grid box and its neighbour
!                   have same surface type

      Real :: fac1              ! Factors used for
      Real :: fac2              ! computation of final
      Real :: fac3              ! coefficients
      Real :: fac4              !
      Real :: num1              ! Temporary storage
      Real :: tempnum           ! Temporary storage

! externals: none

      Do j=first_row,last_row
        Do i=1,row_length
          Do k=1,4
            es_space_interp(k,i,j)=0.
          End Do
        End Do
      End Do

!  Firstly, set all points assuming there are no boundary problems on
!  the whole domain.

      Do j=first_row,last_row
        Do i=1,row_length
          Surf_match(1,i,j)=(land(i,j).EQV.land(I,J+1))
          Surf_match(2,i,j)=(land(i,j).EQV.land(I+1,J))
          Surf_match(3,i,j)=(land(i,j).EQV.land(I,J-1))
          Surf_match(4,i,j)=(land(i,j).EQV.land(I-1,J))
        End Do
      End Do
!
!  Now the problem at the edges of the whole model domain are sorted
!  out. Note that global and non-global runs must be treated
!  differently - the latter type have a problem at the left and right
!  edges of the whole model domain, whereas global models only have
!  problems at the top and bottom of the whole domain.
!
      If (at_top_of_LPG) Then
        Do i=1, row_length
          Surf_match(1,i,last_row)=.FALSE.
        End Do
      End If
!
      If (at_base_of_LPG) Then
        Do i=1, row_length
          Surf_match(3,i,first_row)=.FALSE.
        End Do
      End If
!
      If (model_domain  ==  mt_lam) Then
        If (at_right_of_LPG) Then
          Do j=first_row,last_row
            Surf_match(2,row_length,j)=.FALSE.
          End Do
        End If
!
        If (at_left_of_LPG) Then
          Do j=first_row,last_row
            Surf_match(4,1,j)=.FALSE.
          End Do
        End If
      End If
!
! Now calculate the coefficients for interpolation.
!
      Do j = first_row,last_row
        Do i = 1,row_length
          fac1=0.
          fac2=0.
          fac3=0.
          fac4=0.
          If (Surf_match(1,i,j)) Then
            fac1=1.0
          End If
          If (Surf_match(2,i,j)) Then
            fac2=1.0
          End If
          If (Surf_match(3,i,j)) Then
            fac3=1.0
          End If
          If (Surf_match(4,i,j)) Then
            fac4=1.0
          End If
          tempnum=fac1+fac2+fac3+fac4
          If (tempnum <  0.1) Then
            fac1=1.
            fac2=1.
            fac3=1.
            fac4=1.
          End If
          num1=1./(fac1+fac2+fac3+fac4)
          es_space_interp(1,i,j)=fac1*num1
          es_space_interp(2,i,j)=fac2*num1
          es_space_interp(3,i,j)=fac3*num1
          es_space_interp(4,i,j)=fac4*num1
        End Do
      End Do
!
      RETURN
      END SUBROUTINE COEFFS_DEGRADE
#endif
