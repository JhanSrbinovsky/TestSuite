
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate QT from Q, QCF and QCL

      Subroutine Q_to_QT (                                              &
     &           Q, QCF, QCL, QT, Q_Field_Size, Wet_Levels)

! Description:
!   This routine derives QT from the Q prognostics. QT is a prognostic
!   for the UM prior to the New Dynamics.
!
! Method:
!   Formula used QT = Q + QCL + QCF
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 6.1     18/08/04 Original code. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

      Implicit None

      Integer :: wet_levels
      Integer :: q_field_size

      Real, Intent(IN)  :: Q   (q_field_size, wet_levels)
      Real, Intent(IN)  :: QCF (q_field_size, wet_levels)
      Real, Intent(IN)  :: QCL (q_field_size, wet_levels)
      Real, Intent(OUT) :: QT  (q_field_size, wet_levels)

      Integer :: i, level

      Do level = 1, wet_levels
        Do i = 1, q_field_size
          QT(i,level) = Q(i,level) + QCF(i,level) + QCL(i,level)
        End Do
      End Do

      Return
      END SUBROUTINE Q_to_QT
