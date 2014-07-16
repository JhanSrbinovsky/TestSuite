#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Set the STASH addresses for D1
! Subroutine Interface:

!- End of subroutine code -------------------------------------------


!+Compute data lengths and addresses for primary fields
! Subroutine Interface:
      SUBROUTINE PRIMARY(ISEC,IITM,Im_index,Im_ident,Sm_ident,          &
     &                  RLEVS,RADDRESS,PIrow,ErrorStatus,CMESSAGE)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
!   4.0     Oct. 95                    S.J.Swarbrick
!   4.1     Apr. 96    Generalisation
!                      of routine      S.J.Swarbrick
!   4.2     28/11/96   MPP code : Added calculation of global (dump)
!                                 lengths
! Generalise code for dual-time level prognostics
!                             S.J.Swarbrick
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!   5.2     25/08/00   Changes for implementation of LBCs. The changes
!                      are required because the LBCs now have a section
!                      of their own (ie. not section zero), but still
!                      behave as section zero (ie. primary) fields.
!                      PRIMARY now takes the section number as an
!                      argument.
!                      P.Burton
!
! Code description:
!   FORTRAN 77 + common Fortran 90 extensions.
!   Written to UM programming standards version 7.
!
! System component covered:
! System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "version.h"
#include "parparm.h"
#include "typsize.h"
#include "model.h"
#include "cstash.h"
#include "stextend.h"

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER ISEC      ! Section number
      INTEGER IITM      ! Item number
      INTEGER Im_ident  ! Current internal model number
      INTEGER Im_index  ! Current position in internal model list
      INTEGER Sm_ident  ! Submodel identifier (absolute)
!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      LOGICAL MODEL_LEV
      LOGICAL LADDR
      LOGICAL LMASK
      INTEGER RLEVS      ! No. of levels for reconfiguration
      INTEGER DLEVS      ! No of levels inc pseudo levels
      INTEGER RPLEVS     ! & of pseudo-levels
      INTEGER RADDRESS   ! Address for reconfiguration
      INTEGER I
      INTEGER IL1,IL2
      INTEGER IPL1,IPL2
      INTEGER LEN        ! Data length for primary field
#if defined(MPP)
      INTEGER global_LEN ! Global data length for primary field
#endif
      INTEGER PIrow      ! Counter for ProgItems array

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      EXTERNAL TSTMSK,ADDRLN,LEVCOD,OCNVOL

!- End of Header ---------------------------------------------------

! Find out whether the primary is included for this version
! DEPENDS ON: tstmsk
      CALL TSTMSK(Im_ident,ISEC,LMASK,LADDR,ErrorStatus,CMESSAGE)
      IF (LADDR) THEN
       IF (ISPACE == 10) THEN
! Space code 10 means: no space is required for this item in D1 or
!  the dump, but stashmaster data is required, so an "address" of
!  -1 is set to ensure that the corresponding record will be read
!  into PPXI in routine GET_PPX_PART (called by U_MODEL).
         IN_S(1,Im_ident,ISEC,IITM)=-1
       ELSE

        IF (ISEC  ==  0) THEN
! Start address for model levels in PP array
        PPIND_S(Im_ident,IITM) = NHEAD(Im_ident)+1
        ENDIF ! IF (ISEC  ==  0)

! Find address length per level
#if !defined(MPP)
! DEPENDS ON: addrln
        CALL ADDRLN(IGP,LEN,ErrorStatus)
#else
! DEPENDS ON: addrln
        CALL ADDRLN(IGP,halo_type,LEN,local_data)
! DEPENDS ON: addrln
        CALL ADDRLN(IGP,halo_type,global_LEN,                           &
     &              global_dump_data)
#endif

! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV .OR.(ILEV == 5 .AND. IPSEUDO /= 0)) THEN
! Field has model levels - decode level codes
         IF (ILEV  /=  5) THEN
! DEPENDS ON: levcod
          CALL LEVCOD(IBOT,IL1,ErrorStatus,CMESSAGE)
! DEPENDS ON: levcod
          CALL LEVCOD(ITOP,IL2,ErrorStatus,CMESSAGE)
         ELSE
          IL1=1
          IL2=1
         END IF
! No. of model levels for D1 addressing
          DLEVS=IL2-IL1+1
! Initialise first & last pseudo level indices
          IPL1 =0
          IPL2 =0
          IF (IFLAG == 0.AND.IPSEUDO /= 0) THEN
! Primary with input on all available pseudo levels -
!   decode pseudo level codes
! DEPENDS ON: pslevcod
            CALL PSLEVCOD(IPFIRST,IPL1,'F',ErrorStatus,CMESSAGE)
! DEPENDS ON: pslevcod
            CALL PSLEVCOD(IPLAST ,IPL2,'L',ErrorStatus,CMESSAGE)
            DLEVS=DLEVS*(IPL2-IPL1+1)
          END IF
          RPLEVS=IPL2-IPL1+1
! Multiply length per level by no. of levels
          IF(LEN == -1) THEN         !Grid codes 31,32
! DEPENDS ON: ocnvol
            CALL OCNVOL(LEN,IL1,IL2)
#if defined(MPP)
! DEPENDS ON: ocnvol
            CALL OCNVOL(global_LEN,IL1,IL2)
#endif
          ELSE
            LEN=LEN*(IL2-IL1+1)*(IPL2-IPL1+1)
#if defined(MPP)
            global_LEN=global_LEN*(IL2-IL1+1)*(IPL2-IPL1+1)
#endif
          END IF
          IF (ISPACE /= 4.AND.ISPACE /= 9) THEN
! Increment no. of headers
            NHEAD   (Im_ident)=   NHEAD(Im_ident)                       &
     &                                +(IL2-IL1+1)*(IPL2-IPL1+1)
            NHeadSub(Sm_ident)=NHeadSub(Sm_ident)                       &
     &                                +(IL2-IL1+1)*(IPL2-IPL1+1)
          END IF
        ELSE
! Not model levels
          RLEVS=1
          DLEVS=1
          RPLEVS=1
          IF (ISPACE /= 4.AND.ISPACE /= 9) THEN
            NHEAD   (Im_ident)=NHEAD   (Im_ident)+1
            NHeadSub(Sm_ident)=NHeadSub(Sm_ident)+1
          END IF
        END IF

! The input start address for primary (m,0,i) is assigned
!  to IN_S(1,m,0,i).
! Addresses are set relative to the beginning of the primary data,
!  since the primary data starts at the beginning of D1.
        IF(ISPACE /= 5) THEN
          IF(ISPACE /= 9) THEN
! Start address for this primary field
          IN_S(1,Im_ident,ISEC,IITM)=LPRIM(Sm_ident)+1
! Increment LPRIM by LEN (=data length for this primary field)
          LPRIM  (Sm_ident)      =LPRIM  (Sm_ident)+LEN
          LPrimIM(Im_ident)      =LPrimIM(Im_ident)+LEN
! Information for preliminary D1 addressing array
          N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
          IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
            D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=prog
            D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
            D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IITM
            D1_PADDR(d1_levs,N_OBJ_D1(Sm_ident),Sm_ident)=DLEVS
            D1_PADDR(d1_sect,N_OBJ_D1(Sm_ident),Sm_ident)=ISEC
          ENDIF
#if defined(MPP)
          global_LPRIM  (Sm_ident) =global_LPRIM  (Sm_ident)+global_LEN
          global_LPrimIM(Im_ident) =global_LPrimIM(Im_ident)+global_LEN
#endif
! Dual addresses for ocean fields with dual time level
          IF(ISPACE == 8) THEN
            LPRIM_O2             =LPRIM_O2+LEN
! Information for preliminary D1 addressing array
            N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
            IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
              D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=extra_d1
              D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
              D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IITM
              D1_PADDR(d1_levs,N_OBJ_D1(Sm_ident),Sm_ident)=DLEVS
              D1_PADDR(d1_sect,N_OBJ_D1(Sm_ident),Sm_ident)=ISEC
            ENDIF
          END IF
          ELSE ! Space = 9
!           These are EXNER etc items. Record the address relative
!           to start of LEXTRA space in D1. A loop in ADDRES
!           will then add on LPRIM and LDUMP
            IN_S(1,Im_ident,ISEC,IITM)=LEXTRA(Sm_ident)+1
            LEXTRA(Sm_ident) = LEXTRA(Sm_ident)+LEN
! Information for preliminary D1 addressing array
            N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
            IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
              D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=extra_d1
              D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
              D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IITM
              D1_PADDR(d1_levs,N_OBJ_D1(Sm_ident),Sm_ident)=DLEVS
              D1_PADDR(d1_sect,N_OBJ_D1(Sm_ident),Sm_ident)=ISEC
            ENDIF
          ENDIF
        ELSE
! ISP=5 means: set address of prim var in dump only.
! D1 address is then set to same address as previous item
          IF (IITM  ==  1) THEN
            IN_S(1,Im_ident,ISEC,IITM)=1
          ELSE
            IN_S(1,Im_ident,ISEC,IITM)=IN_S(1,Im_ident,ISEC,IITM-1)
          ENDIF
        END IF
! The input length for primary (m,0,i) is assigned to IN_S(2,m,0,i).
        IN_S(2,Im_ident,ISEC,IITM)=LEN

       END IF  ! ISPACE  /=  10
      END IF   ! LADDR

      RETURN
      END SUBROUTINE PRIMARY

!+Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
!- End of Function code --------------------------------------------

!+Decode the STASH pseudo level code
! Subroutine Interface:
#endif
