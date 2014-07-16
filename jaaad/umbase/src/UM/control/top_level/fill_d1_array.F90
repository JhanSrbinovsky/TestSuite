#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE INITCTL-------------------------------------------------
!LL
!LL   programmers of some or all of previous code & changes include:
!LL    M J CARTER  S TETT   P.TREVELYAN  C WILSON  T JOHNS
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  16/02/93  Pass pseudo-level info to DIAGDESC for printout.
!LL  3.1   03/02/93 : added comdeck CHSUNITS to define NUNITS for i/o.
!LL 3.2    27/03/93 Dynamic allocation of main data arrays. R. Rawlins
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Remove A_MAX_VARIABLES, Add read of PPINDEX.
!LL   3.4    07/12/94 M.Carter. Change to SI_LEN  to calculate
!LL                   STASH_MAXLEN because of STOCGT and different
!LL                   lengths needed in SI_LEN and STASH_LIST
!LL   3.5    02/02/95 M.Carter. Correction to the calculation of
!LL                   MAX_STASH_LEVELS to properly account for levels
!LL                   and pseudo-levels. Bug Fix.
!LL   3.5    Apr. 95    Submodels project.
!LL                   Routine substantially modified. No longer reads
!LL                   from STASH control file; instead, the STASH list,
!LL                   STASH index, and STASH addresses and lengths are
!LL                   passed in via arrays set up by the STASH_PROC code
!LL                   PP_LEN2_LOOK, FT_OUTPUT values also passed in from
!LL                   STASH_PROC.
!LL                     S.J.Swarbrick
!LL  4.0  18/10/95  Remove GET_FILE from EXTERNAL statement. RTHBarnes
!LL  4.1     Apr. 96  Rationalise *CALLs & SI addressing for
!LL                    atmos items 4&5 and 10&11         S.J.Swarbrick
!LL   4.4    05/09/97 Step over space code 10 items S.D.Mullerworth
!LL 4.3-4.4   16/09/97 Added subroutine FILL_D1_ARRAY at 4.3. Plus
!LL                    minor correction at 4.4 S.D.Mullerworth
!LL   5.0    08/06/99 Set up D1_ADDR(halo_type) entry      P.Burton
!    5.0   29/04/99  Introduce conditional printing of messages
!                    dependent on PrintStatus variable. R Rawlins
!LL  5.0  21/05/99 Remove refererences to ..DA (dynamically allocated)
!LL                variables previously needed for portability.
!LL                Also include switch for to avoid calling INITMOS
!LL                if no MOS output requested. R.Rawlins
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL  5.1  13/07/00  Change PrintStatus test for formatted STASH
!LL                 descriptions. R Rawlins
!LL  5.2  25/08/00  Change to FILL_D1_ARRAY to take account of
!LL                 section number information contained in D1_PADDR
!LL                                                         P.Burton
!LL  5.2  18/09/00  Remove redundant code re thetal,qt. R Rawlins
!    5.3  20/08/01  Add call to initialise peer output files
!                   S.D.Mullerworth
!LL  5.3  25/09/00  Add halotype to printout. D Robinson
!    6.1  24/06/04  Extend Fill_D1_array for section 33 tracers. RBarnes
!    6.2  11/04/05  Allow users to override the number of
!                   fields in a fieldsfile. P.Selwood.
!LL
!LL
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DP NO. 3, VERSION 3
!LL
!LL  SYSTEM TASK: C4
!LL
!LL  SYSTEM COMPONENTS: C30, C40
!LL
!LL  PURPOSE:   Initialises STASH control arrays from STASH control
!LL            file.
!LL
!LL  EXTERNAL DOCUMENTATION: UMDP NO. C4 VERSION NO. 4
!LL
!LLEND-------------------------------------------------------------


!LL  SUBROUTINE FILL_D1_ARRAY------------------------------------------
!LL
!LL  PURPOSE: Fill D1 addressing array with useful information.
!LL           S.D.Mullerworth

      SUBROUTINE FILL_D1_ARRAY(                                         &
#include "argsts.h"
#include "argppx.h"
#include "argd1.h"
     &                  ICODE,CMESSAGE)

      IMPLICIT NONE

#include "csubmodl.h"
#include "parparm.h"
#include "typsize.h"
! Contains *CALL CPPXREF
#include "typsts.h"
! Contains *CALL VERSION
#include "ppxlook.h"
#include "chsunits.h"
#include "chistory.h"
#include "stparam.h"
#include "c_mdi.h"
#include "cstash.h"
#include "model.h"
! Declares arrays used in STASH_PROC code (LIST_S etc.);
#include "stextend.h"
                !   also contains common block STEXTEND
! For accessing D1 addressing array
#include "typd1.h"
! Print status information in CPRINTST:
#include "cprintst.h"

      INTEGER                                                           &
     &  II,                                                             &
                 ! Addresses preliminary array
     &  SM,                                                             &
                 ! Addresses final array=1 for 1st submod =2 for 2nd
                 ! submodel etc
     &  TYPE,                                                           &
                 ! Code for prognostic, diagnostic, secondary or other
     &  IOBJ,                                                           &
                 ! Addresses final array
     &  ISEC,                                                           &
                 ! Section number
     &  ITM,                                                            &
                 ! Item number
     &  LEVS,                                                           &
                 ! No of levels
     &  INF,                                                            &
              ! Diagnostic STASHlist number or prognosic item number
     &  Im_ident,                                                       &
     &  Sm_ident,                                                       &
     &  LOOKUP_PTR,                                                     &
                    ! Pointer to lookup table
     &  EXT_ADDR,                                                       &
                  ! Temporary pointer
     &  ICODE                   ! OUT: Error return code
!
      CHARACTER*256                                                     &
     &    CMESSAGE               ! OUT: Error return message

      INTEGER EXPPXI
      EXTERNAL EXPPXI

! Initialise array
      DO Sm_ident=1,N_SUBMODEL_PARTITION
        DO II=1,N_OBJ_D1_MAX
          DO INF=1,D1_LIST_LEN
            D1_ADDR(INF,II,Sm_ident)=-1
            NO_OBJ_D1(Sm_ident)=0
          ENDDO
        ENDDO
      ENDDO

      IF(PrintStatus >= PrStatus_Oper) THEN
! Set up addressing of D1
      WRITE(6,*)'Addressing of D1 array'
      WRITE(6,*)'Key to Type:'
      WRITE(6,*)'Type=0: Prognostic'
      WRITE(6,*)'Type=1: Diagnostics in dump'
      WRITE(6,*)'Type=2: Secondary diagnostics'
      WRITE(6,*)'Type=3: Others (eg P_EXNER in atmos or 2nd of dual '
      WRITE(6,*)'        time level ocean fields)'
      ENDIF  ! PrintStatus test
      SM=0
      DO Sm_ident=1,N_SUBMODEL_PARTITION_MAX
        IOBJ=0
        SM=SUBMODEL_FOR_SM(Sm_ident)
        IF (SM /= 0) THEN
         IF (NO_OBJ_D1(SM) == 0) THEN
          NO_OBJ_D1(SM)=N_OBJ_D1(Sm_ident)
      WRITE(6,*)'Submodel id ',Sm_ident
      WRITE(6,*)'Submodel Number ',SM
          WRITE(6,*)'No of objects in this submodel: ',NO_OBJ_D1(SM)
! Address if submodel not empty and not already addressed
          DO II=1,NO_OBJ_D1(SM)
!           Preliminary array held in D1_PADDR - full array in D1_ADDR
!           Index II in D1_PADDR goes into index IOBJ of D1_ADDR
!           First add prognostics followed by diagnostics...
            Im_ident=D1_PADDR(d1_im,II,Sm_ident)
            INF=D1_PADDR(d1_extra_info,II,Sm_ident)
            ISEC=D1_PADDR(d1_sect,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == prog) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_section,IOBJ,SM)=ISEC
              D1_ADDR(d1_no_levels,IOBJ,SM)=                            &
     &          D1_PADDR(d1_levs,II,Sm_ident)
              D1_ADDR(d1_object_type,IOBJ,SM)=prognostic
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= IN_S(1,Im_ident,ISEC,INF)
            ELSEIF (TYPE == diag) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_object_type,IOBJ,SM)=diagnostic
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= STLIST(st_output_addr,INF)
            ENDIF
          ENDDO
!         Calculate end position of progs and diags for ocean
          IF(SM_IDENT == O_SM)THEN
            EXT_ADDR=LPrimIM(O_IM)+LDumpIM(O_IM)+1
          ENDIF

!         Extra data between primary and secondary diagnostics
          DO II=1,NO_OBJ_D1(SM)
            ISEC=D1_PADDR(d1_sect,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == extra_d1) THEN
              Im_ident=D1_PADDR(d1_im,II,Sm_ident)
              INF=D1_PADDR(d1_extra_info,II,Sm_ident)
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
             D1_ADDR(d1_section,IOBJ,SM)=ISEC
              D1_ADDR(d1_no_levels,IOBJ,SM)=                            &
     &          D1_PADDR(d1_levs,II,Sm_ident)
              D1_ADDR(d1_object_type,IOBJ,SM)=other
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              IF(SM_IDENT /= O_IM)THEN
!               NOT OCEAN: Address was calculated in ADDRES
                D1_ADDR(d1_address,IOBJ,SM)=IN_S(1,Im_ident,ISEC,INF)
              ELSE
!               OCEAN: This is first time 2nd timestep prognostics
!               have been addressed so calculate
                D1_ADDR(d1_address,IOBJ,SM)=EXT_ADDR
                EXT_ADDR=EXT_ADDR+IN_S(2,Im_ident,ISEC,INF)
              ENDIF
            ENDIF
          ENDDO
!         Finally add secondary diagnostics
          DO II=1,NO_OBJ_D1(SM)
!           Preliminary array held in D1_PADDR - full array in D1_ADDR
            Im_ident=D1_PADDR(d1_im,II,Sm_ident)
            INF=D1_PADDR(d1_extra_info,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == seco) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_object_type,IOBJ,SM)=secondary
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= STLIST(st_output_addr,INF)
            ENDIF
          ENDDO

          LOOKUP_PTR=0
          DO II=1,NO_OBJ_D1(SM)
            TYPE= D1_ADDR(d1_object_type,II,SM)
            ISEC= D1_ADDR(d1_section,II,SM)
            INF = D1_ADDR(d1_stlist_no,II,SM)
            Im_ident = D1_ADDR(d1_imodl,II,SM)
            IF((TYPE == prognostic).OR.(TYPE == other))THEN
! Prognostics don't have STASHlist numbers
              D1_ADDR(d1_stlist_no,II,SM)= -1
              D1_ADDR(d1_item,II,SM)   = INF
              D1_ADDR(d1_length,II,SM) = IN_S(2,Im_ident,ISEC,INF)
              ISEC = D1_ADDR(d1_section,II,SM)
              ITM  = INF
!-------------------------------------------------------------------
! Prognostic items:
! Additional items can be added to the array here. Its code (eg
! d1_item, d1_levels) should be added to the TYPD1 comdeck and
! set as a parameter. The D1_LIST_LEN parameter should be changed
! as required
!-------------------------------------------------------------------
            ELSE
              D1_ADDR(d1_section,II,SM)= STLIST(st_sect_code,INF)
              D1_ADDR(d1_item,II,SM)   = STLIST(st_item_code,INF)
              D1_ADDR(d1_length,II,SM) = STLIST(st_output_length,INF)
              ISEC=D1_ADDR(d1_section,II,SM)
              ITM=D1_ADDR(d1_item,II,SM)
! STASH list pointer to D1 address information
              STLIST(st_position_in_d1,INF) = II
!-------------------------------------------------------------------
! Diagnostic items
! Add items as per prognostics
!-------------------------------------------------------------------
              D1_ADDR(d1_north_code,II,SM)    =STLIST(st_north_code,INF)
              D1_ADDR(d1_south_code,II,SM)    =STLIST(st_south_code,INF)
              D1_ADDR(d1_east_code,II,SM)     =STLIST(st_east_code,INF)
              D1_ADDR(d1_west_code,II,SM)     =STLIST(st_west_code,INF)
              D1_ADDR(d1_gridpoint_code,II,SM)=STLIST(s_grid,INF)
              D1_ADDR(d1_proc_no_code,II,SM)  =STLIST(s_proc,INF)
! 1. Number of levels
              IF(STLIST(st_output_bottom,INF) == 100) THEN
! Special levels
                LEVS=1
              ELSE IF(STLIST(st_series_ptr,INF) /= 0) THEN
! Time series domain
                LEVS=1
              ELSE IF(STLIST(st_gridpoint_code,INF) >= 10               &
     &            .AND.STLIST(st_gridpoint_code,INF) <  20) THEN
! Vertical ave.
                LEVS=1
              ELSE  IF(STLIST(st_output_bottom,INF) <  0) THEN
! Levels list
                LEVS=LEVLST_S(1,-STLIST(st_output_bottom,INF))
              ELSE
! Range of model levels
                LEVS=STLIST(st_output_top   ,INF)                       &
     &            -STLIST(st_output_bottom,INF)+1
              END IF

              IF (STLIST(st_pseudo_out,INF) >  0) THEN
! Pseudo levels
                LEVS=LEVS*LENPLST(STLIST(st_pseudo_out,INF))
              END IF
              D1_ADDR(d1_no_levels,II,SM) = LEVS
            ENDIF
!-------------------------------------------------------------------
! Items whose settings are common to progs and diags (eg from PPXREF)
! Add items as per prognostics
! ISEC and ITM set above
!-------------------------------------------------------------------
            D1_ADDR(d1_grid_type,II,SM) =                               &
! DEPENDS ON: exppxi
     &        EXPPXI(Im_ident,ISEC,ITM,ppx_grid_type,                   &
#include "argppx.h"
     &        ICODE, CMESSAGE)
            D1_ADDR(d1_halo_type,II,SM) =                               &
! DEPENDS ON: exppxi
     &        EXPPXI(Im_ident,ISEC,ITM,ppx_halo_type,                   &
#include "argppx.h"
     &        ICODE, CMESSAGE)
            LOOKUP_PTR=LOOKUP_PTR+D1_ADDR(d1_no_levels,II,SM)
            D1_ADDR(d1_lookup_ptr,II,SM)=LOOKUP_PTR
          ENDDO
          IF(PrintStatus >= PrStatus_Normal) THEN
          WRITE(6,*)                                                    &
     &'      Type Modl Sect Item   Address   Length Levels Gridtype',   &
     &' Halotype'
          DO II=1,NO_OBJ_D1(SM)
            WRITE(6,'(5I5,I11,I9,I6,I7,I8)')                            &
     &     II,D1_ADDR(d1_object_type,II,SM),D1_ADDR(d1_imodl,IOBJ,SM),  &
     &        D1_ADDR(d1_section,II,SM),D1_ADDR(d1_item,II,SM),         &
     &        D1_ADDR(d1_address,II,SM),D1_ADDR(d1_length,II,SM),       &
     &        D1_ADDR(d1_no_levels,II,SM),D1_ADDR(d1_grid_type,II,SM),  &
     &        D1_ADDR(d1_halo_type,II,SM)

          ENDDO
          ENDIF  ! PrintStatus test
         ENDIF ! IF (NO_OBJ_D1(SM) == 0) THEN
        ENDIF

      ENDDO ! DO Sm_ident=1,N_SUBMODEL_PARTITION_MAX

  999 CONTINUE
      RETURN
      END SUBROUTINE FILL_D1_ARRAY

#endif
