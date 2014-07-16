
!+ Interpolation routine for radiation calculations
!-
!  Purpose: To interpolate the fields from grid boxes where radiation
!           calculations were performed to grid boxes where the
!           calculations should have been done, but weren't because
!           of spatial degradation of calculations.
!
!  Method:  Quite straightforward. FIRST_DATA_INTERP contains the
!           information on where the first grid box, on the first row,
!           needs interpolation, and from the knowledge of a chequer-
!           board pattern we can locate in a simple arithmetic fashion
!           all points which need interpolation. Note that calculations
!           performed inside data part of the PE only (i.e. not halo).
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
      SUBROUTINE INTERP_FIELD(es_space_interp,field,levels,             &
     &                        first_data_interp,row_length,num_rows,    &
     &                        first_row, last_row, offx, offy )
!
      IMPLICIT NONE
!
      Integer  ,Intent(IN) :: levels
!            The no. of levels in the input field
      Integer  ,Intent(IN) :: row_length
!            The length of a row on the PE (control level)
      Integer  ,Intent(IN) :: num_rows
!            The number of rows excluding halo and polar rows
      Integer  ,Intent(IN) :: first_data_interp
!              The first data point, in data co-ords (ie excluding
!              the halo), which needs to be interpolated on a PE
      Integer  ,Intent(IN) :: first_row
!            The first row that needs to be interpolated
      Integer  ,Intent(IN) :: last_row
!            The last row that needs to be interpolated
      Integer  ,Intent(IN) :: offx, offy
!            The halo widths of the input field
!
      Real  ,Intent(IN) :: es_space_interp(4, row_length, num_rows)
!            The coefficients for radiation interpolation

      Real  ,Intent(INOUT) :: field(1-offx:row_length+offx,             &
     &                                  1-offy:num_rows+offy,levels)
!            The field to be interpolated to neighbouring grid boxes
!
!     Local variables
!
      Integer :: i              ! loop variable
      Integer :: j              ! loop variable
      Integer :: k              ! loop variable
      Integer :: delta          ! temporary storage variable

      If (first_row == 1) Then
        delta=first_data_interp
      Else
        delta=MOD(first_data_interp+1,2)
      End If
!
      Do k=1,levels
        Do j=first_row,last_row,2
          Do i=1+delta,row_length,2
            field(i,j,k)=es_space_interp(1,i,j)*field(i,j+1,k)          &
     &            + es_space_interp(2,i,j)*field(i+1,j,k)               &
     &            + es_space_interp(3,i,j)*field(i,j-1,k)               &
     &            + es_space_interp(4,i,j)*field(i-1,j,k)
          End Do
        End Do
        Do j=first_row+1,last_row,2
          Do i=2-delta,row_length,2
            field(i,j,k)=es_space_interp(1,i,j)*field(i,j+1,k)          &
     &            + es_space_interp(2,i,j)*field(i+1,j,k)               &
     &            + es_space_interp(3,i,j)*field(i,j-1,k)               &
     &            + es_space_interp(4,i,j)*field(i-1,j,k)
          End Do
        End Do
      End Do
!
!
      RETURN
      END SUBROUTINE INTERP_FIELD
