#if defined(A05_4A) || defined(A05_5A) || defined(A05_0A)
! Description:
!  Parameters for the calculation of T at the lifting condensation level.
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
! See Bolton 1980 Mon Wea Rev 108, P1046-1053 for details of empirical 
! relationship

      REAL,PARAMETER::                                                  &
     &  a_bolton = 55.0                                                 &
     &, b_bolton = 2840.0                                               &
     &, c_bolton = 3.5                                                  &
     &, d_bolton = 4.805
      
#endif
