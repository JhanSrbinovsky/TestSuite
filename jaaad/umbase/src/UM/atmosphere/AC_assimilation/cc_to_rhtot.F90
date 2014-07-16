#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES VANMOPS_MIXED_PHASE and CC_TO_RHTOT-------------------
!LL
!LL  Purpose : Performs vertical analysis of MOPS cloud data
!LL            for (2A) mixed phase cloud microphysics scheme
!LL     When IPASS=1,
!LL     model field is interpolated to ob locations and increments
!LL     calculated. Preliminary weight normalisation factors are also
!LL     calculated.
!LL
!LL     When IPASS=2,
!LL     data density is interpolated to ob locations and
!LL     final weight normalisation factors are calculated.
!LL
!LL     Documentation in a Working Paper by B Macpherson & D Wilson is
!LL     on Metweb at URL
!LL     http://fr2010/~frff/papers/MOPS_for_new_micro.ps.gz
!LL
!LL  For use on Cray
!LL
!LL
!LL B Macpherson<- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history :
!LL version  Date
!LL   4.5   18/11/97  New routine developed from VANMOPS B Macpherson
!LL   5.2   30/11/00  amend for new dynamics      B Macpherson
!LL                   NB oper 4.5 mods not yet included
!LL   5.3   08/06/01  Remove duplicate declarations.  A van der Wal
!     5.3   17/09/01  fix incorrect changes to
!                     theta and q.  B Macpherson
!     6.2   07/11/05  Include empirically adj. cloud frac. D. Wilson
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND-------------------------------------------------------------
!
!*L  Arguments:---------------------------------------------------


      SUBROUTINE CC_TO_RHTOT  (CC,NPTS,RHC,RHTOT,K,L_eacf,BL_LEVELS)

      IMPLICIT NONE

      INTEGER                                                           &
     & NPTS                                                             &
                           ! IN No. of points on level.
     &,K                                                                &
                           ! IN: Level no.
     &,BL_LEVELS           ! IN: Number of BL levels

      REAL                                                              &
     & RHC                                                              &
                           ! IN Critical relative humidity (fraction).
     &,CC(NPTS)                                                         &
                           ! IN Cloud cover (fraction).
     &,RHTOT(NPTS)                                                      &
                           ! OUT Total relative humidity (fraction).
     &,TEMPNUM

      LOGICAL                                                           &
                           !, INTENT(IN)
     & L_eacf              ! Use empirically adjusted cloud fraction

!-------------------------------------------------------------------
!     External subroutine calls NONE
!-------------------------------------------------------------------

! Local variables---------------------------------------------------

      INTEGER I     ! Do loop index

! Code in calling routine restricts CC to range 0-1
! but check to be safe

       DO I=1,NPTS
         IF (CC(I) >  1.0) CC(I) = 1.0
         IF (CC(I) <  0.0) CC(I) = 0.0

         If (L_eacf) Then
           ! Use empirically adjusted cloud fraction

           IF (K <= BL_LEVELS) THEN
             IF (CC(I) <= 0.5) THEN
               TEMPNUM = (-1.+SQRT(2.*CC(I)))*0.816 - 0.184
               RHTOT(I) = TEMPNUM*(1.-RHC) + 1.
             ELSE
               TEMPNUM = (1.-SQRT(2.*(1.-CC(I))))*0.816 - 0.184
               RHTOT(I) = TEMPNUM*(1.-RHC) + 1.
             ENDIF
           ELSE
             IF (CC(I) <= 0.5) THEN
               TEMPNUM = (-1.+SQRT(2.*CC(I)))*0.9045 - 0.0955
               RHTOT(I) = TEMPNUM*(1.-RHC) + 1.
             ELSE
               TEMPNUM = (1.-SQRT(2.*(1.-CC(I))))*0.9045 - 0.0955
               RHTOT(I) = TEMPNUM*(1.-RHC) + 1.
             ENDIF
           ENDIF

         Else  ! L_eacf
           ! Use normal method to calculate RHtot

         IF (CC(I) >  0.5) THEN
!  working paper eqn 9b
           RHTOT(I) = 1. + (1.-RHC) * (1. - SQRT(2.*(1.-CC(I))  )  )
         ELSE
!  working paper eqn 9a
           RHTOT(I) = 1. + (1.-RHC) * ( SQRT(2.*CC(I))-1. )
         ENDIF

         End If  ! L_eacf
       END DO

      RETURN
      END SUBROUTINE CC_TO_RHTOT

#endif
