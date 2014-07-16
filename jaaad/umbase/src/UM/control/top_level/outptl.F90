#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Calc stash list output lens; reset boundary spec for full area output.

! Subroutine Interface:

      SUBROUTINE OUTPTL(                                                &
#include "argppx.h"
     &                  NRECS,ErrorStatus,CMESSAGE)
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
!   4.1     Apr. 96    Inclusion of wave model.  S.J.Swarbrick
!   4.2     03/09/96   MPP code : Made addressing local - a processor
!                      only holds data in its geographical area.
!                                                            P. Burton
!   4.3     13/03/97   MPP fixes                             P.Burton
!LL 4.4     22/11/96   Altered for new case of daily mean timeseries
!                      R A Stratton.
!   4.4     12/06/97   Corrected code for global wrap around P.Burton
!   4.5     03/09/98   Don't set the decomposition for the
!                      slab model. Prior to this there would be an
!                      error message if the slab model was selected
!                      for mpp runs. Slab model can now be run with
!                      mpp selected.        C. D. Hewitt
!   4.5     23/01/98   Set up st_dump_level_output_length.  P.Burton
!LL 5.0     23/06/99   Added halo_type argument to GLOBAL_TO_LOCAL
!LL                    Update names of PARVARS variables
!LL                                                       P.Burton
!   5.0     25/06/99   Remove MPP precompile switches        P.Burton
!   5.1     15/05/00   S-N ordering consistency correction. R Rawlins
!   5.3     04/10.01   Removed indication that slab was running non-MPP
!                                                           K.Williams
!   5.5     08/08/00   Modification for parallelisation of WAM.
!                      Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!   6.2     15/08/05   Free format fixes. P.Selwood
!
!   6.2     23/11/05   Removed all references to the wavemodel.
!                      T.Edwards
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"
#include "parvars.h"
#include "decomptp.h"

! Subroutine arguments
!   Array arguments with intent(in):
      INTEGER NRECS

!   Array arguments with intent(out):
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local variables
      INTEGER output_length
      INTEGER IE
      INTEGER IN
      INTEGER IP_DIM
      INTEGER IREC
      INTEGER IS
      INTEGER MODL
      INTEGER ISEC
      INTEGER ITEM
      INTEGER IT_DIM
      INTEGER IW
      INTEGER IX_DIM
      INTEGER IY_DIM
      INTEGER IZ_DIM
      INTEGER                                                           &
! local versions of the global subdomain boundaries
     &  local_north,local_east,local_south,local_west                   &
     &, local_IN,local_IE,local_IS,local_IW                             &
! global versions of the X and Y horizontal dimensions, and
! total output size
     &, global_IX_DIM,global_IY_DIM,global_output_length                &
! variables indicating the decomposition type at various stages
     &, orig_decomp,decomp_type

! Function and subroutine calls:
      INTEGER  EXPPXI
      EXTERNAL EXPPXI,LLTORC

!- End of Header --------------------------------------------------

      orig_decomp=current_decomp_type

! Loop over STASH records
      DO IREC=1,NRECS

! Obtain model, section, item for this record
        MODL = LIST_S(st_model_code  ,IREC)
        ISEC = LIST_S(st_sect_no_code,IREC)
        ITEM = LIST_S(st_item_code   ,IREC)

! Set the correct decomposition type for this model
      IF (MODL  ==  ATMOS_IM) THEN
        decomp_type=decomp_standard_atmos
      ELSEIF (MODL  ==  OCEAN_IM) THEN
        decomp_type=decomp_nowrap_ocean
      ELSEIF (MODL  ==  SLAB_IM) THEN
        decomp_type=decomp_standard_atmos
      ELSE
!       Shouldn't get to this
        decomp_type=decomp_unset
        WRITE(6,*) 'OUTPTL : Error'
        WRITE(6,*) 'Unsupported Model ',MODL,' for MPP code'
        CMESSAGE='Unsupported Model for MPP code'
        ErrorStatus=-1
        GOTO 999
      ENDIF

      IF (current_decomp_type  /=  decomp_type) THEN
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_type,ErrorStatus)
        IF (ErrorStatus  /=  0) THEN
          WRITE(6,*) 'OUTPUTL : Error'
          WRITE(6,*) 'Call to CHANGE_DECOMPOSITION failed with ',       &
     &               'decomposition type ',decomp_type
          CMESSAGE='Unsupported decomposition for MPP code'
          GOTO 999
        ENDIF
      ENDIF


! Extract level code, grid type code from ppx lookup array
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL,ISEC,ITEM,ppx_lv_code,                    &
#include "argppx.h"
     &                                ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL,ISEC,ITEM,ppx_grid_type,                  &
#include "argppx.h"
     &                                ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL ,ISEC ,ITEM,ppx_pt_code      ,            &
#include "argppx.h"
     &                                ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        HALO_TYPE = EXPPXI(MODL,ISEC,ITEM,ppx_halo_type,                &
#include "argppx.h"
     &                     ErrorStatus,CMESSAGE)
        IF (LIST_S(st_proc_no_code,IREC) == 0) THEN
! Dummy record - output length zero
          LIST_S(st_output_length,IREC)=0
        ELSE IF(LIST_S(st_input_code,IREC) <  0.and.                    &
     &         LIST_S(st_proc_no_code,IREC) /= 8.and.                   &
     &           (LIST_S(st_output_code,IREC)  /=  -89)) THEN
! Only copy parent's length for non MOS data. MOS data is
! stored as a global field, so is longer than the parent data
! which is distributed across fields
! Child record - get output length from parent
          LIST_S(st_output_length,IREC)=                                &
     &    LIST_S(st_output_length,-LIST_S(st_input_code,IREC))

        ELSE
! Neither dummy nor child - calculate output length
!   T dimension (equals 1 except for the time series case)
          IF((LIST_S(st_proc_no_code,IREC) == 1).OR.                    &
     &       (LIST_S(st_proc_no_code,IREC) == 2).OR.                    &
     &       (LIST_S(st_proc_no_code,IREC) == 3).OR.                    &
     &       (LIST_S(st_proc_no_code,IREC) == 5).OR.                    &
     &       (LIST_S(st_proc_no_code,IREC) == 6))     THEN
            IT_DIM=1
          ELSE IF (LIST_S(st_proc_no_code,IREC) == 4.or.                &
     &            LIST_S(st_proc_no_code,IREC) == 8) THEN
! Time series case
            IT_DIM=                                                     &
     &      LIST_S(st_period_code,IREC)/LIST_S(st_freq_code,IREC)

          ELSE
            WRITE(6,*)'OUTPTL: ERROR UNEXPECTED PROCESSING CODE',       &
     &      LIST_S(st_proc_no_code,IREC),'   FOR RECORD ',IREC
          END IF

          IF( LIST_S(st_series_ptr,IREC) == 0) THEN
! Set up local versions of the boundaries of the subdomain

! DEPENDS ON: global_to_local_subdomain
            CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE., .TRUE.,             &
     &                               IGP,HALO_TYPE,mype,                &
     &                               LIST_S(st_south_code,IREC),        &
     &                               LIST_S(st_east_code,IREC),         &
     &                               LIST_S(st_north_code,IREC),        &
     &                               LIST_S(st_west_code,IREC),         &
     &                               local_south,local_east,            &
     &                               local_north,local_west)

! Not a time series profile
!   X dimension
            IF(LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  20) THEN
              IX_DIM=local_east-local_west+1
              global_IX_DIM=LIST_S(st_east_code,IREC)-                  &
     &                      LIST_S(st_west_code,IREC)+1
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  30) THEN
              IX_DIM=1
              global_IX_DIM=1
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  40) THEN
              IX_DIM=local_east-local_west+1
              global_IX_DIM=LIST_S(st_east_code,IREC)-                  &
     &                      LIST_S(st_west_code,IREC)+1
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IX_DIM=1
              global_IX_DIM=1
            ELSE
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! X dim

            IF(IX_DIM <  1) THEN
! Area cut by global model
! DEPENDS ON: lltorc
              CALL LLTORC(IGP,90,-90,0,360,IN,IS,IW,IE)
! DEPENDS ON: global_to_local_subdomain
              CALL GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
     &          .TRUE. , .TRUE. , IGP ,halo_type, mype ,                &
     &          IN,IE,IS,IW,                                            &
     &          local_IN,local_IE,local_IS,local_IW)
              IX_DIM=IX_DIM+local_IE-2*halosize(1,halo_type)
! Subtract two halos, because we don't want wrap around to include
! the halo at the end, and the beginning of field

            END IF
            IF (global_IX_DIM <  1) THEN
! DEPENDS ON: lltorc
              CALL LLTORC(IGP,90,-90,0,360,IN,IS,IW,IE)
              global_IX_DIM=global_IX_DIM+IE
            ENDIF

!   Y dimension
            IF(LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  30) THEN
              IF (IGP >= 60.AND.IGP <  70) THEN
! Wave model grid - first lat is southern most
                IY_DIM=local_north-local_south+1
                global_IY_DIM=LIST_S(st_north_code,IREC)-               &
     &                        LIST_S(st_south_code,IREC)+1
              ELSE
! Atmos grid - first lat is southern most
                IY_DIM=local_north-local_south+1
                global_IY_DIM=LIST_S(st_north_code,IREC)-               &
     &                        LIST_S(st_south_code,IREC)+1
              END IF
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <= 40) THEN
              IY_DIM=1
              global_IY_DIM=1
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IY_DIM=1
              global_IY_DIM=1
            ELSE
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! Y dim

!   Z dimension
            IF(LIST_S(st_gridpoint_code,IREC) <  0) THEN
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  10) THEN
              IF(ILEV == 5) THEN
                IZ_DIM=1
              ELSE IF(LIST_S(st_output_bottom,IREC) <  0) THEN
                IZ_DIM=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
              ELSE
                IZ_DIM=LIST_S(st_output_top,IREC)-                      &
     &                 LIST_S(st_output_bottom,IREC)+1
              END IF
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <  20) THEN
              IZ_DIM=1
            ELSE IF(LIST_S(st_gridpoint_code,IREC) <= 43) THEN
              IF(ILEV == 5) THEN
                IZ_DIM=1
              ELSE IF(LIST_S(st_output_bottom,IREC) <  0) THEN
                IZ_DIM=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
              ELSE
                IZ_DIM=LIST_S(st_output_top,IREC)-                      &
     &                 LIST_S(st_output_bottom,IREC)+1
              END IF
            ELSE
              WRITE(6,*)'OUTPTL: ERROR UNEXPECTED GRIDPOINT CODE',      &
     &        LIST_S(st_gridpoint_code,IREC),'   FOR RECORD ',IREC
            END IF  ! Z dim

!   P dimension - pseudo levels
            IF(IPSEUDO >  0) THEN
              IP_DIM=LENPLST(LIST_S(st_pseudo_out,IREC))
            ELSE
              IP_DIM=1
            END IF

! Output length - total number of points
            IF (LIST_S(st_output_code,IREC)  ==  -89) THEN
! this is MOS data
!              output_length = IT_DIM*global_IX_DIM*global_IY_DIM*
!     &                        IZ_DIM*IP_DIM
! The previous two lines replaced by the following - as this
! size is just used to dimension work arrays in stwork which
! deal with the field level by level. This makes the work array
! considerably smaller.

              output_length = IT_DIM*global_IX_DIM*global_IY_DIM*       &
     &                        IP_DIM

            ELSE
              output_length = IT_DIM*IX_DIM*IY_DIM*IZ_DIM*IP_DIM
            ENDIF

            global_output_length =                                      &
     &        IT_DIM*global_IX_DIM*global_IY_DIM*IZ_DIM*IP_DIM
            LIST_S(st_output_length,IREC) = output_length
            LIST_S(st_dump_output_length,IREC) =                        &
     &        global_output_length
            LIST_S(st_dump_level_output_length,IREC) =                  &
     &        global_IX_DIM*global_IY_DIM  ! size of horizontal field
          ELSE    ! Time series profile
            LIST_S(st_output_length,IREC)=                              &
     &       NRECS_TS(LIST_S(st_series_ptr,IREC))*IT_DIM+               &
     &      (NRECS_TS(LIST_S(st_series_ptr,IREC))+1)*6

            LIST_S(st_dump_output_length,IREC)=                         &
     &        LIST_S(st_output_length,IREC)
          END IF
        END IF    ! Neither dummy nor child
      END DO      ! Loop over STASH records
      IF ((orig_decomp  /=  current_decomp_type) .AND.                  &
     &    (orig_decomp  /=  decomp_unset)) THEN
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(orig_decomp,ErrorStatus)
      ENDIF

 999  RETURN
      END SUBROUTINE OUTPTL

!- End of subroutine code -------------------------------------------
#endif
