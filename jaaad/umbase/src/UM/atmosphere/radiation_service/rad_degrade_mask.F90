#if defined(A70_1B)||defined(A70_1C)||defined(A70_1Z)
!+ Create a chequerboard array of logicals over the whole domain
!-
!  Purpose: To create a logical array which has a true chequer-board
!           pattern over the whole model DOMAIN (not just the PE).
!
!  Method:  The datastart variable stores the co-ordinate of the top-
!           left grid box of a PE in terms of the DOMAIN co-ordinates.
!           It can be used in conjunction with the co-ordinates of the
!           grid box in terms of the PE to ensure that the whole domain,
!           and not just a single PE, has a chequer-board pattern of
!           radiation calculations.
!
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
      SUBROUTINE RAD_DEGRADE_MASK(rad_mask, datastart, Ndim_max,        &
     &                  first_row, last_row, offx, row_length, rows)
!
      IMPLICIT NONE
!
      Integer, Intent(IN) :: row_length
!            The data length of a row on the PE (control level)
      Integer, Intent(IN) :: rows
!            The number of rows on the PE (control level)
      Integer, Intent(IN) :: first_row
!            first non-polar row of field (misses halo for MPP code)
      Integer, Intent(IN) :: last_row
!            last non-polar row of field (misses halo for MPP code)
      Integer, Intent(IN) :: offx
!            halo size in EW direction
      Integer, Intent(IN) :: Ndim_max
!            maximum number of spatial dimensions
      Integer, Intent(IN) :: datastart(Ndim_max)
!            position of personal data in global data
!            (set at very high level in control code)
!
      Logical, Intent(OUT) :: rad_mask(row_length, rows)
!            A chequer-board array of values such that the whole DOMAIN
!            is a consistent chequerboard


!     Local variables

      Integer :: i              ! loop variable
      Integer :: j              ! loop variable
      Integer :: tempint        ! temporary storage variable

! Note that all calculations must be done in the data part of the array,
! and that the co-ordinates of such data points are relative to the
! first data point of the whole domain.
! The formula below is a fast method of producing a chequer-board
! pattern over the whole domain.

      Do J=1,last_row
        Do I=1,row_length
          tempint=MOD((datastart(1)+I+datastart(2)+J),2)
          rad_mask(i,j)=(tempint == 0)
        End Do
      End Do
!
      RETURN
      END SUBROUTINE RAD_DEGRADE_MASK
#endif
