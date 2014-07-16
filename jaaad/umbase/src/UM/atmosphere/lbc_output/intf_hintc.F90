#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INTF_HINTC
!LL
!LL Purpose : Calculate horizontal interpolation coefficients to get
!LL           required interface data for a limited area grid.
!LL
!LL Control routine for CRAY YMP
!LL
!LL  Model            Modification history :
!LL version  Date
!LL   3.1   15/12/92  New routine written by D. Robinson
!LL   3.2   14/05/93  Dynamic allocation changes. D Robinson
!LL   3.4   30/03/94  DEF LBOUTA replaced by LOGICAL LLBOUTA
!LL                                          S.J.Swarbrick
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL  3.4  29/11/94  Add lengths to arg.list for portable dyn.allocn.
!    4.1  19/01/96  Replaced references of model dimensions with the
!                   DA versions (for MPP code)    P.Burton
!LL  4.3  19/02/97  Correct coeff3/4 for MPP code. RTHBarnes.
!LL  4.5  29/07/98  Rename CINTF to CINTFA. D. Robinson.
!LL  5.0  19/05/99  Remove DA variabels
!LL                 Change variable names
!LL                                                        P.Burton
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!    6.0  05/09/03  Add new def to allow use in makebc. R.Sempers
!    6.1  18/08/04  Cater for interpolation between C v-grid and B
!                   v grid. D Robinson.
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL System components covered : D810
!LL
!LL System task : D81
!LL
!LL Documentation : Unified Model Documentation Paper No D8
!LLEND

!*L   Arguments:

      SUBROUTINE INTF_HINTC (                                           &
#include "argduma.h"
#include "arginfa.h"
     & JINTF,LEN_INTF_P,LEN_INTF_U,                                     &
     & ICODE,CMESSAGE,LLBOUTA)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typinfa.h"

      INTEGER                                                           &
     &  JINTF                                                           &
                       ! Index to interface area
     & ,LEN_INTF_P                                                      &
                       ! Length of interface p* field
     & ,LEN_INTF_U     ! Length of interface u  field

      INTEGER                                                           &
     &       ICODE     ! Return code

      CHARACTER*(80)                                                    &
     &       CMESSAGE  ! Error message

      LOGICAL LLBOUTA  ! Replaces DEF LBOUTA

!    Local variables
      INTEGER                                                           &
     &        I,J,IJ,                                                   &
                         ! DO loop indices
     &        IP_P,IP_U,                                                &
                         ! Pointers to int coeffs for area JINTF
     &        ROW,                                                      &
                         ! Loop index for rows in rimwidth
     &        IRIM       ! Rim point number

      LOGICAL                                                           &
     &       CYCLIC,                                                    &
                       ! =T, if input grid cyclic
     &       ROT_IN    ! =T, if input grid rotated

      REAL                                                              &
     &  AP_LAMBDA_TARG(LEN_INTF_P),                                     &
                                    ! Longitude coordinates of target
!                                   ! P grid in degrees using same
!                                   ! rotation as source grid
     &  AU_LAMBDA_TARG(LEN_INTF_U),                                     &
                                    ! Longitude coordinates of target
!                                   ! V grid in degrees using same
!                                   ! rotation as source grid
     &  AP_PHI_TARG(LEN_INTF_P),                                        &
                                    ! latitude coordinates of target
!                                   ! P grid in degrees using same
!                                   ! rotation as source grid
     &  AU_PHI_TARG(LEN_INTF_U)     ! latitude coordinates of target
!                                   ! V grid in degrees using same
!                                   ! rotation as source grid

#if defined(ATMOS)
!     Use of U_ROWS and U_FIELD gives compile errors if
!     *DEF OCEAN is activated. Needs attention when routine
!     adapted for ocean use.
      REAL                                                              &
     & LAMBDA(LEN_INTF_P)                                               &
                                    !
     &,PHI   (LEN_INTF_P)           !

!     Latitudes and longitudes for model p, u and v grids.
      Real, dimension (:), allocatable :: ap_phi_srce
      Real, dimension (:), allocatable :: au_phi_srce
      Real, dimension (:), allocatable :: av_phi_srce
      Real, dimension (:), allocatable :: ap_lambda_srce
      Real, dimension (:), allocatable :: au_lambda_srce
      Real, dimension (:), allocatable :: av_lambda_srce

!     Arrays required if model grid is rotated.
      Real, dimension (:), allocatable :: lambda_rot
      Real, dimension (:), allocatable :: lambda_inn
      Real, dimension (:), allocatable :: phi_inn

      Integer :: istat
      Integer :: info

#endif

!*----------------------------------------------------------------------
!L    Subroutines called:
      EXTERNAL H_INT_CO,EQTOLL,LLTOEQ,W_COEFF

!L    Internal Structure

#if defined(ATMOS)
      IF (LLBOUTA) THEN

      ICODE=0
      CMESSAGE=' '

!L 1.0 Get positions in interpolation coeff. arrays for this area

      IP_P = 1
      IP_U = 1
      IF (JINTF >  1) THEN
        DO J=1,JINTF-1
          If (LBC_ND(J) == 0) Then
          IP_P = IP_P + LEN_INTFA_P(J)
          IP_U = IP_U + LEN_INTFA_U(J)
          End If
        ENDDO
      ENDIF

!L 2.0 Set up interpolation constants

! Logical to indicate input grid is rotated
      ROT_IN=A_REALHD(5) /= 90..OR.A_REALHD(6) /= 0.

! Logical to indicate if input data cyclic
      CYCLIC=A_FIXHD(4) <  3

!L 2.1 Calculate coordinates of limited area P grid

!L     Northern points

      IRIM=1

      DO ROW=1,INTFWIDTHA(JINTF)
        DO I=1,INTF_ROW_LENGTH(JINTF)
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-1)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-1)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Eastern points

      DO ROW=INTFWIDTHA(JINTF)+1,INTF_P_ROWS(JINTF)-INTFWIDTHA(JINTF)
        DO I=INTF_ROW_LENGTH(JINTF)-INTFWIDTHA(JINTF)+1,                &
     &       INTF_ROW_LENGTH(JINTF)
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-1)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-1)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Southern points

      DO ROW=INTF_P_ROWS(JINTF)+1-INTFWIDTHA(JINTF),INTF_P_ROWS(JINTF)
        DO I=1,INTF_ROW_LENGTH(JINTF)
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-1)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-1)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Western points

      DO ROW=INTFWIDTHA(JINTF)+1,INTF_P_ROWS(JINTF)-INTFWIDTHA(JINTF)
        DO I=1,INTFWIDTHA(JINTF)
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-1)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-1)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

! DEPENDS ON: eqtoll
      CALL EQTOLL(PHI,LAMBDA,AP_PHI_TARG,AP_LAMBDA_TARG,                &
     &     INTF_POLELAT(JINTF),INTF_POLELONG(JINTF),LEN_INTF_P)

      IF (ROT_IN) THEN

! DEPENDS ON: lltoeq
        CALL LLTOEQ (AP_PHI_TARG,AP_LAMBDA_TARG,AP_PHI_TARG,            &
     &               AP_LAMBDA_TARG,A_REALHD(5),A_REALHD(6),LEN_INTF_P)

        DO I=1,LEN_INTF_P
          IF (AP_LAMBDA_TARG(I) >  180.) THEN
            AP_LAMBDA_TARG(I) = AP_LAMBDA_TARG(I)-360.
          ENDIF
        ENDDO

      ENDIF

!L 2.2 Calculate coordinates of limited area U grid

!L     Northern points

      IRIM=1

      DO ROW=1,INTFWIDTHA(JINTF)
        DO I=1,INTF_ROW_LENGTH(JINTF)-1
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-.5)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-.5)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Eastern points

      DO ROW=INTFWIDTHA(JINTF)+1,INTF_P_ROWS(JINTF)-1-INTFWIDTHA(JINTF)
        DO I=INTF_ROW_LENGTH(JINTF)-INTFWIDTHA(JINTF),                  &
     &       INTF_ROW_LENGTH(JINTF)-1
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-.5)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-.5)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Southern points

      DO ROW=INTF_P_ROWS(JINTF)-INTFWIDTHA(JINTF),INTF_P_ROWS(JINTF)-1
        DO I=1,INTF_ROW_LENGTH(JINTF)-1
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-.5)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-.5)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

!L     Western points

      DO ROW=INTFWIDTHA(JINTF)+1,INTF_P_ROWS(JINTF)-1-INTFWIDTHA(JINTF)
        DO I=1,INTFWIDTHA(JINTF)
          LAMBDA(IRIM)=INTF_FIRSTLONG(JINTF)+(I-.5)*INTF_EWSPACE(JINTF)
          PHI(IRIM)=INTF_FIRSTLAT(JINTF)-(ROW-.5)*INTF_NSSPACE(JINTF)
          IRIM=IRIM+1
        ENDDO
      ENDDO

! DEPENDS ON: eqtoll
      CALL EQTOLL(PHI,LAMBDA,AU_PHI_TARG,AU_LAMBDA_TARG,                &
     &     INTF_POLELAT(JINTF),INTF_POLELONG(JINTF),LEN_INTF_U)

!L 2.3 Calculate coefficients for rotation of winds

! DEPENDS ON: w_coeff
      CALL W_COEFF(COEFF1(IP_U),COEFF2(IP_U),AU_LAMBDA_TARG,            &
     &     LAMBDA,INTF_POLELAT(JINTF),INTF_POLELONG(JINTF),LEN_INTF_U)

! Block of code deleted only required if input model grid is rotated.
! Skipped to speed up code development.
! If required, then code needs upgrading.

!--------------------------------------------------------------------
      allocate ( ap_phi_srce(glsize(2,fld_type_p)) )
      allocate ( au_phi_srce(glsize(2,fld_type_u)) )
      allocate ( av_phi_srce(glsize(2,fld_type_v)) )
      allocate ( ap_lambda_srce(glsize(1,fld_type_p)) )
      allocate ( au_lambda_srce(glsize(1,fld_type_u)) )
      allocate ( av_lambda_srce(glsize(1,fld_type_v)) )

! Source grid P-points latitude and longitude
      Do I=1,glsize(2,fld_type_p)
        AP_PHI_SRCE(I)=A_REALHD(3)+A_REALHD(2)*(I-1)
      ENDDO
      Do I=1,glsize(1,fld_type_p)
       AP_LAMBDA_SRCE(I)=A_REALHD(4)+A_REALHD(1)*(I-1)
      ENDDO

      IF (ROT_IN) THEN
        DO I=1,ROW_LENGTH
          IF (AP_LAMBDA_SRCE(I) >  180.) THEN
            AP_LAMBDA_SRCE(I) = AP_LAMBDA_SRCE(I)-360.
          ENDIF
        ENDDO
      ENDIF

!L Interpolation constants for P_grid interpolation

! DEPENDS ON: h_int_co
      CALL H_INT_CO(AP_INDEX_B_L(IP_P), AP_INDEX_B_R(IP_P),             &
     &              AP_WEIGHT_T_R(IP_P),AP_WEIGHT_B_R(IP_P),            &
     &              AP_WEIGHT_T_L(IP_P),AP_WEIGHT_B_L(IP_P),            &
     &              AP_LAMBDA_SRCE,AP_PHI_SRCE,AP_LAMBDA_TARG,          &
     &              AP_PHI_TARG,                                        &
     &              glsize(1,fld_type_p),glsize(2,fld_type_p),          &
     &              LEN_INTF_P,CYCLIC)

!     Source grid U-points latitude and longitude

      Do I=1,glsize(2,fld_type_u)
        AU_PHI_SRCE(I)=A_REALHD(3)+A_REALHD(2)*(I-1)
      ENDDO
      Do I=1,glsize(1,fld_type_u)
       AU_LAMBDA_SRCE(I)=A_REALHD(4)+A_REALHD(1)*(I-1+0.5)
      ENDDO

      IF (ROT_IN) THEN
        DO I=1,ROW_LENGTH
          IF (AU_LAMBDA_SRCE(I) >  180.) THEN
            AU_LAMBDA_SRCE(I) = AU_LAMBDA_SRCE(I)-360.
          ENDIF
        ENDDO
      ENDIF

!L Interpolation constants for U_grid interpolation

! DEPENDS ON: h_int_co
      CALL H_INT_CO(AU_INDEX_B_L(IP_U), AU_INDEX_B_R(IP_U),             &
     &              AU_WEIGHT_T_R(IP_U),AU_WEIGHT_B_R(IP_U),            &
     &              AU_WEIGHT_T_L(IP_U),AU_WEIGHT_B_L(IP_U),            &
     &              AU_LAMBDA_SRCE,AU_PHI_SRCE,AU_LAMBDA_TARG,          &
     &              AU_PHI_TARG,                                        &
     &              glsize(1,fld_type_u),glsize(2,fld_type_u),          &
     &              LEN_INTF_U,CYCLIC)

!     Source grid V-points latitude and longitude

      DO I=1,glsize(2,fld_type_v)
        AV_PHI_SRCE(I)=A_REALHD(3)+A_REALHD(2)*(I-1+0.5)
      ENDDO
      DO I=1,glsize(1,fld_type_v)
        AV_LAMBDA_SRCE(I)=A_REALHD(4)+A_REALHD(1)*(I-1)
      ENDDO

      IF (ROT_IN) THEN
        DO I=1,glsize(1,fld_type_v)
          IF (AV_LAMBDA_SRCE(I) >  180.) THEN
            AV_LAMBDA_SRCE(I) = AV_LAMBDA_SRCE(I)-360.
          ENDIF
        ENDDO
      ENDIF

!L Interpolation constants for V_grid interpolation

! DEPENDS ON: h_int_co
      CALL H_INT_CO(AV_INDEX_B_L(IP_U), AV_INDEX_B_R(IP_U),             &
     &              AV_WEIGHT_T_R(IP_U),AV_WEIGHT_B_R(IP_U),            &
     &              AV_WEIGHT_T_L(IP_U),AV_WEIGHT_B_L(IP_U),            &
     &              AV_LAMBDA_SRCE,AV_PHI_SRCE,                         &
     &              AU_LAMBDA_TARG,AU_PHI_TARG,                         &
     &              glsize(1,fld_type_v),glsize(2,fld_type_v),          &
     &              LEN_INTF_U,CYCLIC)

      deallocate ( ap_phi_srce )
      deallocate ( au_phi_srce )
      deallocate ( av_phi_srce )
      deallocate ( ap_lambda_srce )
      deallocate ( au_lambda_srce )
      deallocate ( av_lambda_srce )

      END IF     !   LLBOUTA
#endif


!L 4   End of routine

      RETURN
      END SUBROUTINE INTF_HINTC

!-----------------------------------------------------------------------


#endif
