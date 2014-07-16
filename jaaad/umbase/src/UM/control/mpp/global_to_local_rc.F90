#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_SUBDOMAIN: converts global subdomain boundaries
!                            to local subdomain boundaries
! GLOBAL_TO_LOCAL_RC: converts global row,column co-ordinates to
!                     processor co-ordinates plus local
!                     co-ordinates within the processor.
!
! Subroutine Interface:

! Subroutine Interface:
      SUBROUTINE GLOBAL_TO_LOCAL_RC(grid_code,halo_type,                &
     &                              global_column_in , global_row,      &
     &                              processor_x , processor_y,          &
     &                              local_column, local_row)

      IMPLICIT NONE
!
! Description:
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.
!
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      17 /09/96 New deck created for MPP code.  P.Burton
!  4.3      13/03/97  Various bug fixes               P.Burton
!  4.4      18/06/97  Check that row number is valid      P.Burton
!           06/10/97  Set correct row length and n_rows
!                     in dowhile loop.                    P.Burton
!  5.0      22/6/99  Added halo_type argument
!                    Change references to PARVARS variables to include
!                      the new dimensions
!                                                          P.Burton
!  5.1      28/01/00 Simplified calculations using g_pe_index arrays
!                                                           P.Burton
!
! Subroutine arguments:

      INTEGER                                                           &
     &  grid_code                                                       &
                           ! IN : STASH grid type code
     &, halo_type                                                       &
                           ! IN : which type of halo has the field got
     &, global_column_in                                                &
                           ! IN : global column number
     &, global_row                                                      &
                           ! IN : global row number
     &, processor_x                                                     &
                           ! OUT : processor X (EW) co-ordinate
!                          !       (0->nproc_x)
     &, processor_y                                                     &
                           ! OUT : processor Y (NS) co-ordinate
!                               (0->nproc_y)
     &, local_column                                                    &
                           ! OUT : local column number on processor
     &, local_row          ! OUT : local row number on processor

! Parameters and COMMON blocks
#include "parvars.h"
#include "sterr.h"

! Local variables

      INTEGER                                                           &
     &  global_column                                                   &
                      ! modified version of global_column_in which
!                     ! takes account of domains wrapping over
!                     ! 0 degree longitude
     &, fld_type                                                        &
                      ! field stored on P grid or U grid?
     &, row_len_nh,nrows_nh                                             &
                             ! row_len and n_rows when halos removed
     &, proc  ! loop counter for loop over processors

      INTEGER GET_FLD_TYPE  ! function

! ------------------------------------------------------------------

! Find out if the data is on a mass or velocity grid

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(grid_code)

      IF (fld_type  ==  fld_type_unknown) THEN
        WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',                   &
     &    'field with gridtype code ',grid_code
        WRITE(6,*) 'Unable to process this field.'
        processor_x=st_no_data
        processor_y=st_no_data
        local_column=st_no_data
        local_row=st_no_data
        GOTO 9999
      ENDIF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

      IF (global_column_in  >   glsize(1,fld_type)) THEN
        global_column=MOD(global_column_in+1,glsize(1,fld_type))-1
      ELSE
        global_column=global_column_in
      ENDIF

      IF ((global_column  <   1) .OR.                                   &
     &    (global_row  <   1) .OR.                                      &
     &    (global_row  >   glsize(2,fld_type))) THEN

        WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',                   &
     &  'impossible global row/column co-ordinates ',                   &
     &  'row: ',global_row,' column: ',global_column

        processor_x=st_no_data
        processor_y=st_no_data
        local_column=st_no_data
        local_row=st_no_data

      ENDIF

      processor_x=g_pe_index_EW(global_column)
      processor_y=g_pe_index_NS(global_row)

      proc=processor_x+processor_y*gridsize(1)

      local_column=halosize(1,halo_type)+global_column-                 &
     &             g_datastart(1,proc)+1
      local_row=halosize(2,halo_type)+global_row-                       &
     &          g_datastart(2,proc)+1
 9999 CONTINUE

      RETURN
      END SUBROUTINE GLOBAL_TO_LOCAL_RC

! Function Interface

#endif
