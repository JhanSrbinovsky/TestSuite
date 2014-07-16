#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE CON_SCAV-----------------------------------------------
!LL
!LL  Purpose: Scavenge aerosol by convective precipitation.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.4  15/08/94  New routine. Pete Clark.
!LL  4.5  01/05/98  Restrict murk aerosol calculations to aerosol
!LL                 levels=boundary levels. P.Clark
!LL  5.2  03/11/00  Changed addressing for new dynamics. P.Selwood.
!    5.5  17/04/03 Removal of references to obsolete sections
!                  A05_2A,2C,3B. T. White
!    6.2  03/02/05 Added section 5A
!    6.4  12/12/06 Removed def 3C. R A Stratton
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL  Logical component covered: Part of P26.
!LL
!LL  System task:
!LL
!LL  Documentation: In preparation
!*
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE CON_SCAV(                                              &
     & TIMESTEP,row_length,rows,LEVELS,                                 &
     & CBLEVEL,CTLEVEL,                                                 &
     & RAIN,SNOW,AEROSOL                                                &
     &)
      IMPLICIT NONE

      Integer, Intent(In) ::                                            &
     & row_length                                                       &
                                ! Number of points in row
     &,rows                                                             &
                                ! Number of rows
     &,LEVELS                                                           &
                                ! Number of levels.
     &,CBLEVEL(row_length,rows)                                         &
                                ! Convective cloud base level.
     &,CTLEVEL(row_length,rows) ! Convective cloud top level.

      Real, Intent(In) ::                                               &
     & TIMESTEP                                                         &
                               ! Timestep (s).
     &,RAIN(row_length,rows)                                            &
                               ! Rate of rainfall in this layer from
                               !       above
                               !       (kg per sq m per s).
     &,SNOW(row_length,rows)   ! Rate of snowfall in this layer from
                               !       above
                               !       (kg per sq m per s).

      Real, Intent(InOut) ::                                            &
     & AEROSOL(row_length,rows,levels) ! Aerosol mixing ratio

!*L   External subprogram called :-
!     EXTERNAL None
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
      Real, Parameter ::                                                &
     & krain = 1.0E-4                                                   &
     &,ksnow = 1.0E-4

      Real                                                              &
                      ! Real workspace.
     & RRAIN,RSNOW

!  (b) Others.
      Integer I,J,K    ! Loop counters

! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
      RRAIN=KRAIN*TIMESTEP*3600.0
      RSNOW=KSNOW*TIMESTEP*3600.0

      Do j = 1, rows
        Do i = 1,row_length
          If (CTLEVEL(i,j) > 0) Then
            Do K=1,MIN(CTLEVEL(i,j),LEVELS)
              AEROSOL(i,j,k)=AEROSOL(i,j,k)/                            &
     &          (1.0+RRAIN*RAIN(i,j)+RSNOW*SNOW(i,j))
            End Do
          End If
        End Do
      End Do

      Return
      END SUBROUTINE CON_SCAV
#endif
