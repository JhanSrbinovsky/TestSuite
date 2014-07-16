
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Change in total cloud fraction
! Subroutine Interface:
      SUBROUTINE PC2_TOTAL_CF(                                          &
!      Number of points
     & points                                                           &
!      Input fields
     &,CFL, CFF, DELTACL, DELTACF                                       &
!      Input Output fields
     &,CF)
!
      IMPLICIT NONE
!
! Purpose:
!   Update the total cloud fraction due to changes in
!   liquid and ice cloud fractions.
!
! Method:
!   Assumes a random overlap of the forced cloud with already existing
!   conditions in the gridbox. See Annex D of the PC2 cloud scheme
!   documentation.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!  6.2    04-10-04   Change to minimum overlap (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & POINTS
!       No. of points being processed.
!
      REAL                                                              &
                        !, INTENT(IN)
     & CFL(points)                                                      &
!       Liquid cloud fraction
     &,CFF(points)                                                      &
!       Ice cloud fraction
     &,DELTACL(points)                                                  &
!       Change in liquid cloud fraction
     &,DELTACF(points)
!       Change in ice cloud fraction
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & CF(points)
!       Total cloud fraction
!
!  External functions:
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & DELTAC_1                                                         &
                   ! Change in total cloud fraction due to liquid cloud
     &,DELTAC_2    ! Change in total cloud fraction due to ice cloud
!
!  (b) Others.
      INTEGER I           ! Loop counter
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
! Points_do1:
      DO i=1,points
!
! ----------------------------------------------------------------------
! 1. Update total cloud fraction.
! ----------------------------------------------------------------------
!
! Calculate change in total cloud fraction due to a change in liquid
! cloud fraction. This depends upon the sign of the change of liquid
! cloud fraction.
!
        IF (DELTACL(i)  >   0.0 .AND. CFL(i)  <   1.0) THEN
! Random overlap
!         DELTAC_1 = DELTACL(i) *(1.0 - CF(i))/(1.0 - CFL(i))
! Minimum overlap
          DELTAC_1 = MIN( DELTACL(i) , (1.0-CF(i)) )
        ELSE IF (DELTACL(i)  <   0.0 .AND. CFL(i)  >   0.0) THEN
! Random overlap
!         DELTAC_1 = DELTACL(i) * (CF(i)-CFF(i)) / CFL(i)
! Minimum overlap
          DELTAC_1 = MAX( DELTACL(i), (CFF(i)-CF(i)) )
        ELSE
          DELTAC_1 = 0.0
        ENDIF
!
! Calculate change in total cloud fraction due to a change in ice
! cloud fraction. This depends upon the sign of the change of ice
! cloud fraction.
!
        IF (DELTACF(i)  >   0.0 .AND. CFF(i)  <   1.0) THEN
! Random overlap
!         DELTAC_2 = DELTACF(i) *(1.0 - CF(i))/(1.0 - CFF(i))
! Minimum overlap
          DELTAC_2 = MIN( DELTACF(i) , (1.0-CF(i)) )
        ELSE IF (DELTACF(i)  <   0.0 .AND. CFF(i)  >   0.0) THEN
! Random overlap
!         DELTAC_2 = DELTACF(i) * (CF(i)-CFL(i)) / CFF(i)
! Minimum overlap
          DELTAC_2 = MAX( DELTACF(i), (CFL(i)-CF(i)) )
        ELSE
          DELTAC_2 = 0.0
        ENDIF
!
! Sum the two changes
!
        CF(i) = CF(i) + DELTAC_1 + DELTAC_2
!
! For minimum overlap we need to check that the total cloud
! fraction is constrained within 0 and 1
!
          CF(i) = MAX(  MIN( CF(i),1.0 )  , 0.0)
!
! Points_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_TOTAL_CF
