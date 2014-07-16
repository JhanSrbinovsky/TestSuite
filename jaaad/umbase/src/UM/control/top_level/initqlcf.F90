#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INITQLCF------------------------------------------------
!LL
!LL  Purpose: Diagnose (initialise) cloud water (frozen + liquid
!LL           separately) and cloud amount, from standard prognostic
!LL           model variables (Q, T, P*).
!LL
!LL      For initialisation of RHC see the documentation
!LL      of sub-component P292 (large-scale cloud). SRHC1 and SRHC2
!LL      contain the critical relative humidity discussed in the
!LL      paragraph incorporating equations P292.11 - P292.14.  The
!LL      values below are based on those used in the old GCM.
!LL
!LL C.Wilson    <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 4.5      01/05/1998 modify prototype in order to introduce SCM.
!LL                     (JC Thil)
!LL 5.3      02/05/2001 Update to 5.3.  Take out calculation of
!LL                     pressure and pass it in through the argument
!LL                     list instead.                   Z. Gardner
!   5.5      11/02/2003 Fix to allow RHcrit = 1.   R.M.Forbes
!LL
!LL Programming standard :
!LL
!LL Logical components covered : P292
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!LL
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE INITQLCF                                               &
     &  (P,RHCRIT,Q,T,P_LEVELS,ROW_LENGTH, ROWS,OCF,OQCF,OQCL,NBL,JLEV)

      IMPLICIT NONE

      INTEGER                                                           &
     & NBL                                                              &
                           ! IN No of levels in b.l.
     &,JLEV                                                             &
                           ! IN Model level to be processed
     &,ROW_LENGTH                                                       &
                           ! IN No. of points in x direction.
     &,ROWS                                                             &
                           ! IN No. of points in y direction.
     &,P_LEVELS            ! IN No. of hybrid levels in the model.

      REAL                                                              &
     &  RHCRIT(P_LEVELS)                                                &
                           ! Critical relative humidities
     &,P(ROW_LENGTH, ROWS)                                              &
                           ! IN Pressure (Pa).
     &,Q(ROW_LENGTH, ROWS)                                              &
                            ! IN Sp humidity (kg water per kg air).
     &,T(ROW_LENGTH, ROWS)                                              &
                            ! IN Temperature (K).
     &,OCF(ROW_LENGTH, ROWS)                                            &
                            ! OUT Cloud fraction (decimal fraction).
     &,OQCF(ROW_LENGTH, ROWS)                                           &
                              ! OUT Cloud ice content (kg per kg air).
     &,OQCL(ROW_LENGTH, ROWS) ! OUT Cloud liquid water (kg per kg air).

!  External subroutine called-------------------------------------------
      EXTERNAL QSAT
! Workspace usage:-----------------------------------------------------
      REAL                                                              &
                            ! Workspace (see later comments for usage):-
     & WQSAT(ROW_LENGTH, ROWS)  ! WORK
!*---------------------------------------------------------------------
! Local varables:------------------------------------------------------
      REAL                                                              &
                            ! Local (see later comments for usage):-
     & WAC                                                              &
                            ! LOCAL
     &,WRH                                                              &
                            ! LOCAL
     &,ALPHALRCP                                                        &
                            ! LOCAL
     &,QC                   ! LOCAL Cloud water/ice content

#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_0_dg_c.h"
#include "c_pi.h"
      REAL PFIDDLE,PALCON2E,PALDEP2E,PRCP,                              &
                                           ! Derived + local constants.
     &     PC1,PC2,PC3
      PARAMETER (                                                       &
     & PFIDDLE=0.25                                                     &
                                           ! Used in 8/8 cloud case.
     &,PALCON2E=LC*LC*EPSILON                                           &
                                           !
     &,PALDEP2E=(LC+LF)*(LC+LF)*EPSILON                                 &
                                           !
     &,PRCP=R*CP                                                        &
                                           !
     &,PC1=1.060660172                                                  &
                                           ! 3/sqrt(8).
     &,PC2=2.0*PC1                                                      &
                                           !
     &,PC3=PI/3.0                                                       &
                                           ! pi/3
     &)

      REAL SRHC1,SRHC2     !critical relative humidity consts
      DATA SRHC1/0.925/                                                 &
                           !in b.l.
     &,SRHC2/0.85/         !above b.l.
      REAL RHC             !critical relative humidity

      INTEGER I, J     ! Do loop index

!-----------------------------------------------------------------------
!LL 0. Initialise critical relative humidity according to level
!-----------------------------------------------------------------------

        IF(JLEV <  NBL)THEN
           RHC=SRHC1
        ELSE
           RHC=SRHC2
        ENDIF
      RHC = RHCRIT(JLEV)

!-----------------------------------------------------------------------
!LL 1. Calculate pressure (in array W), hence QSAT, hence relative
!LL    humidity in WRH.
!-----------------------------------------------------------------------

! DEPENDS ON: qsat
        CALL QSAT(WQSAT,T(1,1),P,ROW_LENGTH*ROWS)

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            WRH=Q(I,J)/WQSAT(I,J)
          IF(WRH >  1.0)WRH=1.0
!-----------------------------------------------------------------------
!LL 2. Calculate cloud fraction OCF.  (Known as C in formulae).
!-----------------------------------------------------------------------
            OCF(I,J)=0.0

            IF (RHC  <   1.) THEN

              IF (WRH  >   RHC .AND. WRH  <   (5.+RHC)/6.) THEN
                OCF(I,J) = 2.*COS(PC3+ACOS(PC1*(WRH-RHC)/(1.-RHC) )/3.)
                OCF(I,J) = OCF(I,J)*OCF(I,J)
              ENDIF
              IF (WRH  >=  (5.+RHC)/6.) THEN
                OCF(I,J) = PC2*(1.-WRH)/(1.-RHC)
                OCF(I,J) = 1.-OCF(I,J)**(2./3.)
              ENDIF
              IF (OCF(I,J) <  0.0) OCF(I,J)=0.0
              IF (OCF(I,J) >  1.0) OCF(I,J)=1.0

            ELSE ! RHC = 1. so set cloud fraction to 0 or 1

              IF (WRH  <   1.0) THEN
                OCF(I,J) = 0.0
              ELSE
                OCF(I,J) = 1.0
              ENDIF

            ENDIF

!-----------------------------------------------------------------------
!LL 3. Calculate F(C) - store in WAC.
!-----------------------------------------------------------------------
          WAC=0.0
            IF(OCF(I,J) <= 0.5 .AND. OCF(I,J) >  0.0)THEN
              WAC=2.*OCF(I,J)
            WAC=SQRT(WAC*WAC*WAC)/6.
          ENDIF
            IF(OCF(I,J) >  0.5)THEN
              WAC=2.*(1.-OCF(I,J))
            WAC=1.+                                                     &
     &               SQRT(WAC*WAC*WAC)/6.-                              &
     &               SQRT(WAC)
          ENDIF
!-----------------------------------------------------------------------
!LL 4. Calculate A(C) - store in WAC.
!-----------------------------------------------------------------------
          WAC=WAC*(1.-RHC)
!-----------------------------------------------------------------------
!LL 5. Calculate total cloud water - store in QC
!-----------------------------------------------------------------------
            IF(T(I,J) >  TM)THEN
              ALPHALRCP=PALCON2E*WQSAT(I,J)/(PRCP*T(I,J)*T(I,J))
          ELSE
              ALPHALRCP=PALDEP2E*WQSAT(I,J)/(PRCP*T(I,J)*T(I,J))
          ENDIF
!      Special treatment of full liquid cloud cover case.
            IF(OCF(I,J) >= 1. .AND. T(I,J) >= TM)THEN
            QC=PFIDDLE/(1.+ALPHALRCP*(1.+PFIDDLE))
          ELSE
            QC=WAC/(1.+ALPHALRCP*(WRH+WAC))
          ENDIF
            QC=QC*WQSAT(I,J)
!-----------------------------------------------------------------------
!LL 6. Partition cloud water into liquid and ice, and store in output
!LL    arrays OQCL, OQCF.
!-----------------------------------------------------------------------
            OQCL(I,J)=0.0
            OQCF(I,J)=0.0
            IF(T(I,J) >  TM)THEN
              OQCL(I,J)=QC
          ELSE
              OQCF(I,J)=QC
          ENDIF
          ENDDO ! I
        ENDDO ! J


      RETURN
      END SUBROUTINE INITQLCF
#endif
