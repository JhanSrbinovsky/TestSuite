#if defined(C70_1A) || defined(FLDOP) || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE READSTM
!LL     PURPOSE  TO READ A RECORD FROM THE PRE STASH-MASTERS FILE
!LL   AND RETURN THE PPXREF CODES AND NAME OF A GIVEN DIAGNOSTIC
!LL   TESTED UNDER CFT77 ON OS 5.1
!LL
!LL   AUTHOR            M.J.CARTER      DATE 23/01/92
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                   Add read of new record for sub-model type.
!LL   3.4    01/11/94 M.Carter. Change to stash-master file format
!LL                   to allow 4 digits for PP field code.
!LL   3.5    13/03/95 Sub-Models Project:
!LL                   Mods which facilitate expansion of PPXREF file to
!LL                   incorporate the entire preSTASHmaster record, and
!LL                   to have "Internal Model no." as an extra dimension
!LL                   in the file (in addition to section,item).
!LL                   S.J.Swarbrick
!     4.0   15/12/95 Changed interface : Added Int_Model_no as an
!                    extra argument, rather than using CODES(1)
!                    which can overwrite DNAM due to equivalencing
!                    in MKPPXRF1.  P.Burton
!     4.1   Apr. 96  Modified to read new format of STASH master files
!                     S.J.Swarbrick
! vn4.4  10/4/97  Changes to FORMAT statements to allow
!                           NAG f90 compiled code to run. Ian Edmond
!     5.0   29/06/99 New halo type parameter included.
!                    D.M. Goddard
!     5.5   07/04/03 For SX, avoid reading in STASHmaster option code
!                                                          E.Leung
!     5.5    3/02/03  STASHmaster option codes increased from 20 digits
!                     to 30. Code adjusted accordingly. W Roseblade.
!     6.0   18/06/03  Removing UTILIO, FLDOP defs    E.Leung
!LL
!LL   LOGICAL COMPONENT R913
!LL
!LL   PROJECT TASK: C4
!LL
!LL   PROGRAMMING STANDARD  UMDP 4
!LL
!LL   EXTERNAL DOCUMENT C4
!LL
!LLEND
      SUBROUTINE READSTM                                                &
     &          (IMASK,DNAM,CODES,NFT,ICODE,CMESSAGE)
      IMPLICIT NONE

#include "csubmodl.h"
#include "cppxref.h"
#include "version.h"

!     ARGUMENTS
      INTEGER NFT             !IN:    UNIT NUMBER STMSTS

      INTEGER ICODE           !OUT: RETURN CODE
      CHARACTER*(80) CMESSAGE !OUT: RETURN MESSAGE IF THERE IS A FAILURE

      CHARACTER*1 DNAM(PPXREF_CHARLEN) !OUT: VARIABLE NAME FROM RECORD
      INTEGER CODES(PPXREF_CODELEN)    !OUT: PPXREF CODES FROM RECORD
      INTEGER IMASK(20)                !OUT: VERSION MASK

      INTEGER IMSK          ! Decimal equivalent of binary IMASK
      INTEGER II                       !LOCAL: loop mark.
      integer opcod(30)
!
      ICODE=0
      CMESSAGE=' '

      READ(NFT,2010,END=3100,ERR=3200)                                  &
     & CODES(ppx_model_number)  ,                                       &
     & CODES(ppx_section_number),                                       &
     & CODES(ppx_item_number)   , DNAM
 2010 FORMAT(2X,3(I5,2X),36A1)

      IF(CODES(ppx_model_number) == -1) GO TO 9999

      READ(NFT,2110,END=3100,ERR=3200)                                  &
     & CODES(ppx_space_code),                                           &
     & CODES(ppx_ptr_code),                                             &
     & CODES(ppx_timavail_code),                                        &
     & CODES(ppx_grid_type),                                            &
     & CODES(ppx_lv_code),                                              &
     & CODES(ppx_lb_code),                                              &
     & CODES(ppx_lt_code),                                              &
     & CODES(ppx_pt_code),                                              &
     & CODES(ppx_pf_code),                                              &
     & CODES(ppx_pl_code),                                              &
     & CODES(ppx_lev_flag)
 2110 FORMAT(2X,11(I5,2X))

! 30-digit option code read in as 6x5 digit groups instead of 4x5
      READ(NFT,2120,END=3100,ERR=3200)                                  &
     &(opcod(ii),ii=1,30),                                              &
     &(IMASK(II),II=1,20),                                              &
     &CODES(ppx_halo_type)
 2120 FORMAT(3X,30(I1),3X,20(I1),2X,I5)

      CODES(ppx_opt_code  )=                                            &
     & opcod(30)+opcod(29)*10  +opcod(28)*100  +                        &
     &           opcod(27)*1000+opcod(26)*10000
      CODES(ppx_opt_code+1)=                                            &
     & opcod(25)+opcod(24)*10  +opcod(23)*100  +                        &
     &           opcod(22)*1000+opcod(21)*10000
      CODES(ppx_opt_code+2)=                                            &
     & opcod(20)+opcod(19)*10  +opcod(18)*100  +                        &
     &           opcod(17)*1000+opcod(16)*10000
      CODES(ppx_opt_code+3)=                                            &
     & opcod(15)+opcod(14)*10  +opcod(13)*100  +                        &
     &           opcod(12)*1000+opcod(11)*10000
      CODES(ppx_opt_code+4)=                                            &
     & opcod(10)+opcod( 9)*10  +opcod( 8)*100  +                        &
     &           opcod( 7)*1000+opcod( 6)*10000
      CODES(ppx_opt_code+5)=                                            &
     & opcod( 5)+opcod( 4)*10  +opcod( 3)*100  +                        &
     &           opcod( 2)*1000+opcod( 1)*10000

!   Binary version mask was read into array IMASK
!   Convert version mask to decimal form IMSK
          IMSK = 0
          DO II=20,1,-1
            IF((IMASK(II) /= 0).AND.(IMASK(II) /= 1)) THEN
              WRITE(6,*) 'READSTM: improper IMASK in user diag'
              WRITE(6,*) 'Model, Section, Item ',                       &
     &        CODES(ppx_model_number)  ,                                &
     &        CODES(ppx_section_number),                                &
     &        CODES(ppx_item_number)
            ELSE
              IF(IMASK(II) == 1) THEN
                IMSK=IMSK+2**(20-II)
              END IF
            END IF
          END DO
!     Insert decimal value of version mask
          CODES(ppx_version_mask)=IMSK

      READ(NFT,2130,END=3100,ERR=3200)                                  &
     & CODES(ppx_data_type),                                            &
     & CODES(ppx_dump_packing),                                         &
     &(CODES(II),                                                       &
     & II= ppx_pack_acc, ppx_pack_acc+PPXREF_PACK_PROFS-1)
 2130 FORMAT(2X,I5,2X,I5,3X,I3,9(2X,I3))

      READ(NFT,2140,END=3100,ERR=3200)                                  &
     & CODES(ppx_rotate_code),                                          &
     & CODES(ppx_field_code),                                           &
     & CODES(ppx_user_code),                                            &
     & CODES(ppx_lbvc_code),                                            &
     & CODES(ppx_base_level),                                           &
     & CODES(ppx_top_level),                                            &
     & CODES(ppx_ref_lbvc_code),                                        &
     & CODES(ppx_cf_levelcode),                                         &
     & CODES(ppx_cf_fieldcode)
 2140 FORMAT(2X,9(I5,2X))
 3100 GO TO 9999 ! Normal completion
 3200 WRITE(6,*)' MESSAGE FROM ROUTINE READSTM: '
      WRITE(6,*)' ERROR OCCURRED WHILE READING STASHmaster FILE '
      CMESSAGE=' READSTM: ERROR READING STASHMASTERS FILE'
      ICODE=2

 9999 CONTINUE
      RETURN

      END SUBROUTINE READSTM

#endif
