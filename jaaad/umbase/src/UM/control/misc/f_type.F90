#if defined(C70_1A) || defined(PUMF)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE F_TYPE-----------------------------------------------
!LL
!LL  Purpose:  Returns each field code and associated field length from
!LL            the PP header and a count of the number of fields
!LL            of each type.
!LL
!LL  Written by A. Dickinson
!LL
!LL
!LL Programming standard :
!LL
!LL Logical components covered :
!LL
!LL Project task :
!LL
!LL  Documentation: None
!LL
!LLEND----------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE F_TYPE(LOOKUP,LEN2_LOOKUP,PP_NUM,N_TYPES               &
     &,PP_LEN,PP_STASH,PP_TYPE,PP_POS,PP_LS,FIXHD,                      &
#include "argppx.h"
     &TITLE)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN2_LOOKUP                                                      &
                               !IN 2nd dimension of LOOKUP
     &,N_TYPES                                                          &
                               !IN No of separate field types in file
     &,LOOKUP(64,LEN2_LOOKUP)                                           &
                               !IN LOOKUP record
     &,PP_NUM(LEN2_LOOKUP)                                              &
                               !OUT No of successive fields with same co
     &,PP_LEN(LEN2_LOOKUP)                                              &
                               !OUT Length of field
     &,PP_STASH(LEN2_LOOKUP)                                            &
                               !OUT PP code of field
     &,PP_TYPE(LEN2_LOOKUP)                                             &
                               !OUT Integer/real/timeseries
     &,PP_POS(LEN2_LOOKUP)                                              &
                               !OUT Pointer to number of PP field
     &,PP_LS(LEN2_LOOKUP)                                               &
                                !OUT Data stored on land or sea pts
     &,FIXHD(*)

      CHARACTER*(80)TITLE


! Comdecks: ------------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
! Local variables: -----------------------------------------------------
      INTEGER MODEL             !Internal model number from LOOKUP

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & PP_XREF(PPXREF_CODELEN)  !PPXREF codes for a given section/item

! External subroutines called:------------------------------------------
      CHARACTER*36 EXPPXC
      EXTERNAL EXPPXC
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      INTEGER                                                           &
     & ICODE                                                            &
                  ! Error code
     &,ITEM_CODE                                                        &
                  ! STASH item code
     &,SECTION    ! STASH section number

      CHARACTER                                                         &
     & CMESSAGE*80                                                      &
                   ! Error message
     &,PHRASE*(PPXREF_CHARLEN) ! Name of field

      INTEGER I,K
      character(len=30) :: valid(len2_lookup)
      integer :: fc_time(len2_lookup)
      integer :: lbc_levs(len2_lookup)
!*----------------------------------------------------------------------

! Initialise arrays
      DO K=1,LEN2_LOOKUP
        PP_NUM(K)=1
        PP_LEN(K)=0
        PP_STASH(K)=0
        PP_TYPE(K)=0
        PP_POS(K)=0
        PP_LS(K)=0
      ENDDO

      write(valid(1),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i4)')              &
     &                 lookup(1,1),':',lookup(2,1),':',                 &
     &                 lookup(3,1),':',lookup(4,1),':',                 &
     &                 lookup(5,1),':',lookup(6,1)
      fc_time(1)=lookup(14,1)

      N_TYPES=1
      PP_LEN(1)=LOOKUP(15,1)
      PP_STASH(1)=LOOKUP(42,1)
      PP_TYPE(1)=LOOKUP(39,1)
      PP_POS(1)=1
        IF(MOD(INT(LOOKUP(21,1)/10),10) == 2)THEN
          PP_LS(1)=MOD(INT(LOOKUP(21,1)/100),10)
        ENDIF

      DO K=2,LEN2_LOOKUP
        IF(LOOKUP(42,K) == LOOKUP(42,K-1).AND.                          &
     &     (LOOKUP(18,K)*lookup(19,k)) ==                               &
     &     (LOOKUP(18,K-1)*lookup(19,k-1))                              &
     &            .and.                                                 &
     &       lookup(1,k) == lookup(1,k-1) .and.                         &
     &       lookup(2,k) == lookup(2,k-1) .and.                         &
     &       lookup(3,k) == lookup(3,k-1) .and.                         &
     &       lookup(4,k) == lookup(4,k-1) .and.                         &
     &       lookup(5,k) == lookup(5,k-1) .and.                         &
     &       lookup(6,k) == lookup(6,k-1) .and.                         &
     &       lookup(14,k) == lookup(14,k-1) )THEN
          PP_NUM(N_TYPES)=PP_NUM(N_TYPES)+1
        ELSE
          N_TYPES=N_TYPES+1
          PP_LEN(N_TYPES)=LOOKUP(15,K)
          PP_STASH(N_TYPES)=LOOKUP(42,K)
          PP_TYPE(N_TYPES)=LOOKUP(39,K)
          PP_POS(N_TYPES)=PP_POS(N_TYPES-1)+PP_NUM(N_TYPES-1)
          write(valid(n_types),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i4)')    &
     &lookup(1,k),':',lookup(2,k),':',                                  &
     &lookup(3,k),':',lookup(4,k),':',                                  &
     &lookup(5,k),':',lookup(6,k)
          fc_time(n_types)=lookup(14,k)
          if(fixhd(5)==5)then
            lbc_levs(n_types)=lookup(17,k)-100
          endif
        IF(MOD(INT(LOOKUP(21,K)/10),10) == 2)THEN
          PP_LS(N_TYPES)=MOD(INT(LOOKUP(21,K)/100),10)
        ENDIF
        ENDIF
      ENDDO

! Print out details of fields
      WRITE(6,'(''  '',/,'' '',A/)')TITLE

       select case(fixhd(5))
         case(3)
           write(6,*)'Key to output- fields file'
           write(6,*)'Data stored on land/sea points'
           write(6,*)'No of fields'
           write(6,*)'Length of field before unpacking (==LBLREC)'
           write(6,*)'Data type'
           write(6,*)'STASH Code'
           write(6,*)'Start of first field'
           write(6,*)'Description'
           write(6,*)'Validity time (yyyy:mm:dd:hh:mn:day_no)'
           write(6,*)'Forecast period'
         case(5)
           write(6,*)'Key to output - lbc file'
           write(6,*)'Data stored on land/sea points'
           write(6,*)'No of fields'
           write(6,*)'Length of field before unpacking (==LBLREC)'
           write(6,*)'Data type'
           write(6,*)'STASH Code'
           write(6,*)'Start of first field'
           write(6,*)'Description'
           write(6,*)'Validity time (yyyy:mm:dd:hh:mn:day_no)'
         case default
           write(6,*)'Key to output - dump or other'
           write(6,*)'Data stored on land/sea points'
           write(6,*)'No of fields'
           write(6,*)'Length of field before unpacking (==LBLREC)'
           write(6,*)'Data type'
           write(6,*)'STASH Code'
           write(6,*)'Start of first field'
           write(6,*)'Description'
           write(6,*)'Validity time (yyyy:mm:dd:hh:mn:day_no)'
           write(6,*)'Forecast period'
       end select

      I=1
      DO K=1,N_TYPES
        if(lookup(42,i) == -99)then
          exit
        endif
        PHRASE=' '
        ITEM_CODE=MOD(LOOKUP(42,I),1000)
        SECTION=(LOOKUP(42,I)-ITEM_CODE)/1000
        MODEL=LOOKUP(45,I)
        ICODE = 0

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        if(model == 10)then
          model = 1
        endif

! DEPENDS ON: exppxc
        PHRASE=EXPPXC(MODEL,SECTION,ITEM_CODE,                          &
#include "argppx.h"
     &              ICODE,CMESSAGE)
        IF(ICODE /= 0)THEN
          PHRASE='NON-STANDARD FIELD'
        ENDIF
        I=I+PP_NUM(K)
        select case(fixhd(5))
          case(3)
            WRITE(6,'('' '',I4,I5,I8,I4,2I6,2x,A36,a30,i3)')            &
     &      PP_LS(K),PP_NUM(K),PP_LEN(K),                               &
     &      PP_TYPE(K),PP_STASH(K),PP_POS(K),                           &
     &      PHRASE,valid(k),fc_time(k)
          case(5)
            WRITE(6,'('' '',I4,I5,I8,I4,2I6,2x,A36,a30)')               &
     &      PP_LS(K),lbc_levs(K),PP_LEN(K),                             &
     &      PP_TYPE(K),PP_STASH(K),PP_POS(K),                           &
     &      PHRASE,valid(k)
          case default
            WRITE(6,'('' '',I4,I5,I8,I4,2I6,2x,A36,a30,i3)')            &
     &      PP_LS(K),PP_NUM(K),PP_LEN(K),                               &
     &      PP_TYPE(K),PP_STASH(K),PP_POS(K),                           &
     &      PHRASE,valid(k),fc_time(k)
        end select

      ENDDO

      RETURN
      END SUBROUTINE F_TYPE
#endif
