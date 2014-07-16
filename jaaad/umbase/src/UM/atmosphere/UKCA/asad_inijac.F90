#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
! Purpose: To set up the indexing arrays used in the routine asad_jac.f 
!          to compute the Jacobian matrix elements for the IMPACT time i
!          integration routine
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Interface
!     ---------
!     All variables used and set are in common blocks.
!
!     Method
!     ------
!     The basic procedure is to scan the reactions counting those that
!     will contribute to the Jacobian and storing an entry in a list to
!     that reaction. For species of type 'FM' and 'TR' their
!     contribution will not change during the model run. However, for
!     species that go into and out of a family (type 'FT'), we cannot
!     precompute their contribution, since this is dependent on their
!     lifetime which will vary with time and space. It is more efficient
!     therefore in the Jacobian matrix calculation not to be using
!     species of type 'FT'. Another complication is introduced by the
!     use of chemical families. These can make a positive contribution
!     to the main diagonal of the Jacobian matrix and we take account of
!     these separately (see Carver & Stott, 1997, Annales Geophysicae
!     for more details).
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INIJAC

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

!       Local variables

        INTEGER :: ifam(jpmsp)
        INTEGER :: itr(jpmsp)
        INTEGER :: ifneg(2*jpnr,jpctr)
        INTEGER :: ifpos(jppjac,jpctr)
        INTEGER :: ifn(jpctr)
        INTEGER :: ifp(jpctr)
        INTEGER :: ir1                 ! Index
        INTEGER :: ir2                 ! Index
        INTEGER :: ip1                 ! Index
        INTEGER :: ip2                 ! Index
        INTEGER :: ip3                 ! Index
        INTEGER :: is                  ! Index
        INTEGER :: njx                 ! Index
        INTEGER :: j                   ! Loop variable
        INTEGER :: jc                  ! Loop variable
        INTEGER :: je                  ! Loop variable
        INTEGER :: jg                  ! Loop variable
        INTEGER :: jr                  ! Loop variable
        INTEGER :: i
        INTEGER :: irem
        INTEGER :: jpnjcx3

        CHARACTER (LEN=2)  :: itype(jpmsp)
        CHARACTER (LEN=72) :: cmessage     ! Error message

!       1. Loop over reactions and set index arrays.
!          ---- ---- --------- --- --- ----- -------

        jpnjcx3 = (jpnr/(3*3))+3*3

        DO jc = 1, jpctr
          ifn(jc) = 0
          ifp(jc) = 0
        ENDDO

!       1.1  Compute contribution from each reaction.

        DO jr = 1, jpnr
          ir1 = nspi(jr,1)
          ir2 = nspi(jr,2)
          ip1 = nspi(jr,3)
          ip2 = nspi(jr,4)
          ip3 = nspi(jr,5)

!         1.1.1  Set indices to determine species type.

          DO j = 1, jpmsp
            ifam(j)  = 0
            itype(j) = '  '
            itr(j)   = 0
          ENDDO
          IF ( ir1 /= 0 ) then
            ifam(1)  = moffam(ir1)
            itype(1) = ctype(ir1)
            itr(1)   = madvtr(ir1)
          ENDIF
          IF ( ir2 /= 0 ) then
            ifam(2)  = moffam(ir2)
            itype(2) = ctype(ir2)
            itr(2)   = madvtr(ir2)
          ENDIF
          IF ( ip1 /= 0 ) then
            ifam(3)  = moffam(ip1)
            itype(3) = ctype(ip1)
            itr(3)   = madvtr(ip1)
          ENDIF
          IF ( ip2 /= 0 ) then
            ifam(4)  = moffam(ip2)
            itype(4) = ctype(ip2)
            itr(4)   = madvtr(ip2)
          ENDIF
          IF ( ip3 /= 0 ) then
            ifam(5)  = moffam(ip3)
            itype(5) = ctype(ip3)
            itr(5)   = madvtr(ip3)
          ENDIF

!         1.2  If first reactant is a family member add the
!              reaction to the negative and positive lists.

          IF ( ifam(1) /= 0 ) then
            is = -nodd(ir1)
            IF ( ifam(2) == ifam(1) .AND. itype(2) == jpfm ) THEN
              is = is - nodd(ir2)
            END IF
            IF ( ifam(3) == ifam(1) .AND. itype(3) == jpfm ) THEN
              is = is + nodd(ip1)
            END IF
            IF ( ifam(4) == ifam(1) .AND. itype(4) == jpfm ) THEN
              is = is + nodd(ip2)
            END IF
            IF ( ifam(5) == ifam(1) .AND. itype(5) == jpfm ) THEN
              is = is + nodd(ip3)
            END IF
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(ifam(1)) = ifn(ifam(1)) + 1
                IF ( ifn(ifam(1)) >  2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifneg dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(ifam(1)),ifam(1)) = jr
              END DO
            ELSE IF ( is > 0 ) then
              DO je = 1, is
                ifp(ifam(1)) = ifp(ifam(1)) + 1
                IF ( ifp(ifam(1)) >  jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(ifam(1)),ifam(1)) = jr
              END DO
            END IF
          END IF

!         1.3  If second reactant is a family member add to the lists

          IF ( ifam(2) /= 0 ) THEN
            is = -nodd(ir2)
            IF ( ifam(1) == ifam(2) .AND. itype(2) == jpfm ) THEN
              is = is - nodd(ir1)
            END IF
            IF ( ifam(3) == ifam(2) .AND. itype(3) == jpfm ) THEN 
              is = is + nodd(ip1)
            END IF
            IF ( ifam(4) == ifam(2) .AND. itype(4) == jpfm ) THEN
              is = is + nodd(ip2)
            END IF
            IF ( ifam(5) == ifam(2) .AND. itype(5) == jpfm ) THEN
              is = is + nodd(ip3)
            END IF
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(ifam(2)) = ifn(ifam(2)) + 1
                IF ( ifn(ifam(2)) > 2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage = 'ARRAY TOO SMALL: Increase ifneg'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(ifam(2)),ifam(2)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(ifam(2)) = ifp(ifam(2)) + 1
                IF ( ifp(ifam(2)) > jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(ifam(2)),ifam(2)) = jr
              END DO
            END IF
          END IF

!         1.4  Calc. contribution from non-family species.

          IF ( itr(1) /= 0 ) then
            is = -1
            IF ( itr(2) == itr(1) .AND. ir2 /= 0 ) THEN
              is = is - 1
            END IF
            IF ( itr(3) == itr(1) .AND. ip1 /= 0 ) THEN
              is = is + 1
            END IF
            IF ( itr(4) == itr(1) .AND. ip2 /= 0 ) THEN 
              is = is + 1
            END IF
            IF ( itr(5) == itr(1) .AND. ip3 /= 0 ) THEN
              is = is + 1
            END IF

            IF ( is < 0 ) then
              DO je = 1, abs(is)
                ifn(itr(1)) = ifn(itr(1)) + 1
                IF ( ifn(itr(1)) > 2*jpnr ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                  ' in inijac.f. Increase ifneg dimension.'
                  cmessage = 'ARRAY TOO SMALL: Increase ifneg'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(itr(1)),itr(1)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(itr(1)) = ifp(itr(1)) + 1
                IF ( ifp(itr(1)) >  jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(itr(1)),itr(1)) = jr
              END DO
            END IF
          END IF

          IF ( itr(2) /= 0 ) then
            is = -1
            IF ( itr(2) == itr(1) .AND. ir1 /= 0 ) THEN 
              is = is - 1
            END IF
            IF ( itr(3) == itr(2) .AND. ip1 /= 0 ) THEN
              is = is + 1
            END IF
            IF ( itr(4) == itr(2) .AND. ip2 /= 0 ) THEN 
              is = is + 1
            END IF
            IF ( itr(5) == itr(2) .AND. ip3 /= 0 ) THEN 
              is = is + 1
            END IF
            
            IF ( is < 0 ) THEN
              DO je = 1, abs(is)
                ifn(itr(2)) = ifn(itr(2)) + 1
                IF ( ifn(itr(2)) >  2*jpnr ) then
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                 ' in inijac.f. Increase ifneg dimension.'
                  cmessage='ARRAY TOO SMALL: Increase igneg dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifneg(ifn(itr(2)),itr(2)) = jr
              END DO
            ELSE IF ( is > 0 ) THEN
              DO je = 1, is
                ifp(itr(2)) = ifp(itr(2)) + 1
                IF ( ifp(itr(2)) > jppjac ) THEN
                  WRITE (6,*) 'INTERNAL ASAD ERROR: array too small ', &
                  ' in inijac.f. Increase ifpos dimension.'
                  cmessage='ARRAY TOO SMALL: Increase ifpos dimension'
! DEPENDS ON: ereport
                  CALL EREPORT('ASAD_INIJAC',je,cmessage)
                END IF
                ifpos(ifp(itr(2)),itr(2)) = jr
              END DO
            END IF
          END IF
        END DO


!       2.  Partition lists into 3, 2, and 1 sum loops taking
!           positive and negative contributions into account.
!           -------- -------- ------------- ---- --------

        DO jc = 1, jpctr
          njx = ifn(jc) / 3
          IF ( njx > jpnjcx3 ) THEN
            WRITE (6,*) 'FATAL ASAD ERROR: array njacx3 too ',         &
                     'small in inijac.f. Increase jpnjcx3 to at ',     &
                     'least ', njx
            cmessage = 'ARRAY NJACX3 TOO SMALL: Increase jpnjcx3'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_INIJAC',jc,cmessage)
          END IF

!         2.1. Negative 3 term sums.

          i = 1
          njcgrp(jc,3) = njx
          DO jg = 1, njcgrp(jc,3)
            njacx3(1,jg,jc) = ifneg(i,jc)
            njacx3(2,jg,jc) = ifneg(i+1,jc)
            njacx3(3,jg,jc) = ifneg(i+2,jc)
            i               = i + 3
          END DO

!         2.2  Negative 2 and single term sums.

          irem = mod(ifn(jc),3)
          IF ( irem == 2 ) THEN
            njcgrp(jc,2) = 1
            njacx2(1,jc) = ifneg(i,jc)
            njacx2(2,jc) = ifneg(i+1,jc)
            i = i + 2
          ELSE IF ( irem  ==  1 ) THEN
            njcgrp(jc,1) = 1
            njacx1(jc) = ifneg(i,jc)
          END IF
        END DO

!       2.3  Terms which make a positive contribution to the
!            Jacobian. Don't bother to unroll these as there
!            are unlikely to be many of them.

        DO jc = 1, jpctr
          nmpjac(jc) = ifp(jc)
          IF ( ifp(jc) > jppjac ) THEN
            WRITE (6,*) 'FATAL ASAD ERROR: array npjac1 too ',         &
           ' small in injac.f. Increase parameter jppjac to at ',      &
           ' least ',ifp(jc)
            cmessage='ARRAY NPJAC1 TOO SMALL: Increase jppjac'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_INIJAC',jc,cmessage)
          END IF
          DO jr = 1, ifp(jc)
            npjac1(jr,jc) = ifpos(jr,jc)
          END DO
        END DO

        RETURN
        END SUBROUTINE ASAD_INIJAC
#endif
