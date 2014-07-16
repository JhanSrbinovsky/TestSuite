#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Add inactive records to STASH list, when space is required
!
! Subroutine Interface:



!+Find whether ST_list entry (Im_ident,ISEC,ITEM) is an implied diag
! Subroutine Interface:



!+Add diagnostic to the STASH list (LIST_S)
! Subroutine Interface:

      SUBROUTINE ADDIN                                                  &
     &(NRECS,ITEM,ISEC,Im_ident,LBVC,ErrorStatus,CMESSAGE)
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
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.1     May. 96    Various improvements. S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER NRECS
      INTEGER ITEM
      INTEGER ISEC
      INTEGER Im_ident
      INTEGER LBVC

!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! Local scalars:
      LOGICAL MODEL_LEV
      INTEGER IBOT1
      INTEGER ITOP1

! ErrorStatus
      INTEGER ErrorStatus

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      EXTERNAL LEVCOD

!- End of Header ----------------------------------------------------


      LIST_S(st_item_code   ,NRECS)=ITEM
      LIST_S(st_sect_no_code,NRECS)=ISEC
      LIST_S(st_model_code  ,NRECS)=Im_ident
      LIST_S(st_proc_no_code,NRECS)=0

        IF(IFLAG == 1)THEN
! Attempt to do vertical compression - not allowed
          WRITE(6,*)                                                    &
     & 'INACTR: SPACECODE ',ISPACE,' INDICATES IMPLIED DIAGNOSTIC.'
          WRITE(6,*)                                                    &
     & 'NOT ALLOWED WITH LEVEL COMPRESSION FLAG 1 - CHECK STASHMASTER'
          WRITE(6,*)                                                    &
     &   'MODEL ',Im_ident,' SECTION ',ISEC,' ITEM ',ITEM
        ELSE
! DEPENDS ON: disct_lev
          MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
          IF (MODEL_LEV) THEN
! Model levels
! Set bottom level
! DEPENDS ON: levcod
            CALL LEVCOD(IBOT,IBOT1,ErrorStatus,CMESSAGE)
! Set top level
! DEPENDS ON: levcod
            CALL LEVCOD(ITOP,ITOP1,ErrorStatus,CMESSAGE)
            LIST_S(st_input_bottom,NRECS)=IBOT1
            LIST_S(st_input_top   ,NRECS)=ITOP1
          ELSE IF(ILEV == 5) THEN
            LIST_S(st_input_bottom,NRECS)=100
            LIST_S(st_input_top   ,NRECS)=LBVC
          ELSE
            WRITE(6,*)                                                  &
     &     'INACTR: LEVEL TYPE ERROR ON IMPLIED DIAGNOSTIC ',           &
     &     ' - ONLY MODEL LEVELS OR SINGLE LEVEL ALLOWED '
            WRITE(6,*)  'MODEL ',Im_ident,                              &
     &     ' SECT ',ISEC,' ITEM ',ITEM,' LEV CODE ',ILEV
          END IF
        END IF

      LIST_S(st_freq_code      ,NRECS)=1
      LIST_S(st_start_time_code,NRECS)=0
      LIST_S(st_end_time_code  ,NRECS)=0
      LIST_S(st_period_code    ,NRECS)=0
      LIST_S(st_gridpoint_code ,NRECS)=1
      LIST_S(st_weight_code    ,NRECS)=0
      LIST_S(st_north_code     ,NRECS)=0
      LIST_S(st_south_code     ,NRECS)=0
      LIST_S(st_west_code      ,NRECS)=0
      LIST_S(st_east_code      ,NRECS)=0
      LIST_S(st_input_code     ,NRECS)=1
      LIST_S(st_input_length   ,NRECS)=0
      LIST_S(st_output_code    ,NRECS)=0
      LIST_S(st_output_length  ,NRECS)=0
      LIST_S(st_output_addr    ,NRECS)=0
      LIST_S(st_output_bottom  ,NRECS)=0
      LIST_S(st_output_top     ,NRECS)=0
      LIST_S(st_lookup_ptr     ,NRECS)=-1
      LIST_S(st_series_ptr     ,NRECS)=0
      LIST_S(st_macrotag       ,NRECS)=0
      LIST_S(st_pseudo_in      ,NRECS)=0
      LIST_S(st_pseudo_out     ,NRECS)=0
      LIST_S(NELEMP+1          ,NRECS)=NRECS

      RETURN
      END SUBROUTINE ADDIN
#endif
