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
      SUBROUTINE GLOBAL_TO_LOCAL_SUBDOMAIN(                             &
     &                              L_include_halosEW,                  &
     &                              L_include_halosNS,                  &
     &                              grid_code,halo_type,procid,         &
     &                              global_start_row_in,                &
     &                              global_end_col_in,                  &
     &                              global_end_row_in,                  &
     &                              global_start_col_in,                &
     &                              local_start_row,local_end_col,      &
     &                              local_end_row,local_start_col)
      IMPLICIT NONE
!
! Description:
! Takes a global definition of a subdomain region (in terms of
! model gridpoints) and translates it into local numbers.
! This effectively means local co-ordinates of the region of the
! subdomain which intersects with this processor's area.
!
! Method:
! Use the datastart variable in PARVARS to see if the requested
! subdomain intersects with this processor's area, if it does
! then use datastart to convert to local co-ordinate and do a bit
! of logic using MAX and MIN to ensure the local co-ordinates
! actually lie within the local area  Then make any corrections
! necessary to account for a subdomain which crosses over the
! 0 longitude line. Finally, if L_include_halos is set to
! .TRUE. - include any relevant halo regions.
!
! Current code owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      03/09/96 New deck created for MPP code.  P.Burton
!  4.3      13/03/97 Various bug fixes               P.Burton
!  4.4      12/06/97 Another bug fix                 P.Burton
!  5.0      9/6/99   Added halo_type argument
!                    Change references to PARVARS variables to include
!                      the new dimensions
!                    Changed logic to allow for indexing from South
!                                                          P.Burton
!  5.1      28/01/00 Changed variable names to deal with start+end
!                    rather than North,South,East+West which get
!                    confusing with grids which have different
!                    orderings.                           P.Burton
!   5.1     25/01/00 Corrections for uv points on B grid. R.Rawlins
!   5.3     22/11/01 Enable MPP as the only option for
!                    small executables         E.Leung
!  5.3      21/08/01 Added ocean boundary data grid types M. J. Bell
!  5.5      10/08/00 Modification for parallelisation of WAM.
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine arguments:

      LOGICAL                                                           &
     &  L_include_halosEW                                               &
                           ! IN : include East-West halos in local
!                          !      region if set to .TRUE.
     &, L_include_halosNS  ! IN : include North-South halos in local
!                          !      region if set to .TRUE.
      INTEGER                                                           &

     &  grid_code                                                       &
                         ! IN : STASH grid type of field
     &, halo_type                                                       &
                         ! IN : which type of halo has the field got
     &, procid                                                          &
                         ! IN : processor to produce result for
     &, global_start_row_in                                             &
                            ! IN : first row of global subdomain
     &, global_end_col_in                                               &
                            ! IN : last column of global subdomain
     &, global_end_row_in                                               &
                            ! IN : last row of global subdomain
     &, global_start_col_in                                             &
                            ! IN : first column of global subdomain

     &, local_start_row                                                 &
                            ! OUT : first row of local subdomain
     &, local_end_col                                                   &
                            ! OUT : last column of local subdomain
     &, local_end_row                                                   &
                            ! OUT : last row of local subdomain
     &, local_start_col     ! OUT : first column of local subdomain

! Parameters and Common blocks

#include "parvars.h"
#include "sterr.h"

! Local variables
      INTEGER                                                           &
! Copies of the input arguments, that can be modified for
! wrap-around calculations
     &  global_start_row,global_end_col                                 &
     &, global_end_row,global_start_col                                 &
     &, fld_type                                                        &
                  ! is field on P or U grid?
     &, row_len_nh                                                      &
                      ! row length when halos are removed
     &, nrows_nh                                                        &
                      ! number of rows when halos are removed
     &, first_global_pt_EW                                              &
                           ! global point number of first and last
     &, last_global_pt_EW                                               &
                           ! local points in local area
     &, first_global_pt_NS                                              &
                           ! in the East-West and
     &, last_global_pt_NS  ! North-South directions

      LOGICAL                                                           &
! Logicals indicating if this processor contains part of a
! subdomain
     &  NS_intersect,EW_intersect                                       &
     &, wrap                                                            &
             ! set to .TRUE. if the subdomain passes over the
!            ! the 0 degree longitude line
     &, fullfield ! if the field is NOT a subdomain

      INTEGER GET_FLD_TYPE  ! function

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

      global_start_row=global_start_row_in
      global_end_col=global_end_col_in
      global_end_row=global_end_row_in
      global_start_col=global_start_col_in

! Find out if the data is on a mass or velocity grid

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(grid_code)

      IF (fld_type  ==  fld_type_unknown) THEN
        WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',            &
     &    'field with gridtype code ',grid_code
        WRITE(6,*) 'Unable to process this field.'
        local_start_row=st_no_data
        local_end_col=st_no_data
        local_end_row=st_no_data
        local_start_col=st_no_data
        GOTO 9999
      ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

      fullfield= ((global_start_col  ==  1) .AND.                       &
     &            (global_end_col  ==  glsize(1,fld_type)) .AND.        &
     &            (global_start_row  ==  1) .AND.                       &
     &            (global_end_row  ==  glsize(2,fld_type)))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

      IF (fullfield) THEN

        IF (L_include_halosNS) THEN
          local_start_row=1
          local_end_row=g_lasize(2,fld_type,halo_type,procid)
        ELSE
          local_start_row=1+halosize(2,halo_type)
          local_end_row=g_lasize(2,fld_type,halo_type,procid)-          &
     &                  halosize(2,halo_type)
        ENDIF
        IF (L_include_halosEW) THEN
          local_start_col=1
          local_end_col=g_lasize(1,fld_type,halo_type,procid)
        ELSE
          local_start_col=1+halosize(1,halo_type)
          local_end_col=g_lasize(1,fld_type,halo_type,procid)-          &
     &                 halosize(1,halo_type)
        ENDIF

      ELSE ! a subdomain requires some careful analysis:

        row_len_nh=g_blsize(1,fld_type,procid)
        nrows_nh=g_blsize(2,fld_type,procid)

! Set up variables giving the global point numbers of the
! start and end of this processor's subdomain

        first_global_pt_EW=g_datastart(1,procid)
        last_global_pt_EW=first_global_pt_EW+row_len_nh-1

        first_global_pt_NS=g_datastart(2,procid)
        last_global_pt_NS=first_global_pt_NS+nrows_nh-1

! If global_east is greater than the global row length, this
! indicates a wrap around - but not in the format this code
! wants - where it expects a wrap around to be indicated by
! the east column being less than the west column.

        IF (global_end_col  <   global_start_col) THEN
          wrap=.TRUE.
        ELSEIF (global_end_col  >   glsize(1,fld_type)) THEN
          wrap=.TRUE.
          global_end_col=global_end_col-glsize(1,fld_type)
        ELSE
          wrap=.FALSE.
        ENDIF

        EW_intersect =                                                  &
     &    (( .NOT. wrap) .AND.                                          &
     &     ((global_end_col  >=  first_global_pt_EW) .AND.              &
     &      (global_start_col  <=  last_global_pt_EW)))                 &
     &    .OR.                                                          &
     &    ((wrap) .AND.                                                 &
     &     ((global_start_col  <=  last_global_pt_EW) .OR.              &
     &      (global_end_col  >=  first_global_pt_EW)))

        NS_intersect =                                                  &
     &    ((global_start_row  <=  last_global_pt_NS) .AND.              &
     &     (global_end_row  >=  first_global_pt_NS))

        IF (NS_intersect) THEN

          IF ((global_start_row  >=  first_global_pt_NS) .AND.          &
     &        (global_start_row  <=  last_global_pt_NS)) THEN
! This processor contains the NS start of the subarea
            local_start_row=global_start_row-first_global_pt_NS+        &
     &                    halosize(2,halo_type)+1
          ELSE
! This processor is to the North of the start of the subarea
            local_start_row=1+halosize(2,halo_type)
          ENDIF

          IF ((global_end_row  >=  first_global_pt_NS) .AND.            &
     &        (global_end_row  <=  last_global_pt_NS)) THEN
! This processor contains the NS end of the subarea
            local_end_row=global_end_row-first_global_pt_NS+            &
     &                    halosize(2,halo_type)+1
          ELSE
! This processor is to the South of the subarea
            local_end_row=halosize(2,halo_type)+nrows_nh
          ENDIF

        ELSE

          local_start_row=st_no_data
          local_end_row=st_no_data

        ENDIF

        IF (EW_intersect) THEN

          IF ((global_start_col  >=  first_global_pt_EW) .AND.          &
     &        (global_start_col  <=  last_global_pt_EW)) THEN
! This processor contains the EW start of the subarea
            local_start_col=global_start_col-first_global_pt_EW+        &
     &                   halosize(1,halo_type)+1
          ELSE
! This processor is to the right of the start of the subarea
            local_start_col=1+halosize(1,halo_type)
          ENDIF

          IF ((global_end_col  >=  first_global_pt_EW) .AND.            &
     &        (global_end_col  <=  last_global_pt_EW)) THEN
! This processor contains the EW end of the subarea
            local_end_col=global_end_col-first_global_pt_EW+            &
     &                   halosize(1,halo_type)+1
          ELSE
! This processor is to the left of the end of the subarea
            local_end_col=halosize(1,halo_type)+row_len_nh
          ENDIF

        ELSE

          local_start_col=st_no_data
          local_end_col=st_no_data

        ENDIF

      ENDIF ! is this a fullfield?

 9999 CONTINUE

      RETURN
      END SUBROUTINE GLOBAL_TO_LOCAL_SUBDOMAIN

! Subroutine Interface:

! Function Interface

#endif
