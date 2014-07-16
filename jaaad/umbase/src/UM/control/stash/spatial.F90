#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SPATIAL --------------------------------------------------
!LL
!LL  Purpose: Performs general spatial processing on an input field to
!LL           produce an output field or scalar within STASH.  Lower-
!LL           level routines are called to perform the various spatial
!LL           processing options.
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.3  30/03/94  Explicitly declare (sub-addressed) output field
!LL                  fieldout using 'lenout' dimension.  Tim Johns.
!LL   3.3  16/09/93  Pass LOGICAL lmasswt to processing routines to
!LL                  denote that level-by-level mass-weights exist.
!LL   4.3   9/12/96  Added MPP code.
!LL                  Moved calculation of weighting and masking terms
!LL                  up from processing routines.            P.Burton
!LL   4.4   13/06/97 MPP : Where reduction spatial meaning takes place
!LL                  processors not getting results should set
!LL                  their diagnostic space to zeros.          P.Burton
!LL   4.4   22/10/97 MPP : Prevent uninitialised points when
!LL                  pstar_weight on U or C grid S.D.Mullerworth
!LL   5.0   13/07/99 Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0   22/06/99 Added halo_type argument            P.Burton
!LL   5.0   01/11/99 Added halo size arguments
!LL                  Add halos to weighting/mask arrays  P.Burton
!LL   5.0   17/11/99 Added no_rows as an argument, and use to
!LL                  dimension local weighting/masking arrays
!LL                                                      P.Burton
!LL   5.0    9/11/99 Correct to south-north order for ystart,yend
!LL                                                       P.Burton
!LL   5.1   15/05/0 S-N ordering consistency correction. R Rawlins
!     5.3   05/10/01 Correction to mes diagnostics for fields with halos
!                    ie prognostics u,v,etc. R Rawlins
!     5.4   02/09/02 Use sea mask for grid type=3. This is not the
!                    reverse of the land mask when coastal tiling
!                    is used.                            K.Williams
!     5.5   17/02/03 Amendment for Wave model. D.Holmes-Bell
!     6.0   06/10/03 Cater for river grid 23. C.Bunton
!     6.1   16/08/04 Increase first dimension in certain arrays
!                    that are passed to STEXTC. Jens-Olaf Beismann
!                    NEC.   Anthony A. Dickinson
!     6.2   18/10/05 Fix bugs introduced by gan1f601      Andy Malcolm
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D71
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!LL-----------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE SPATIAL(fieldin,vx,vy,vz,GR,st_grid,                   &
     &                   fld_type,halo_type,                            &
     &                   halo_x,halo_y,                                 &
     &                   lcyclic,lmasswt,                               &
     &      n_cols_out,n_rows_out,                                      &
     &      this_index_lev,level_list,index_lev,no_of_levels,           &
     &      no_of_levels_masswt,                                        &
     &      p,pstar,                                                    &
     &      cos_v_latitude,cos_theta_latitude,land,sea,                 &
     &      row_length,rows,n_rows,no_rows,model_levels,                &
     &      fieldout,lenout,                                            &
     &      control,control_size,rmdi,                                  &
     &      icode,cmessage)
!
      IMPLICIT NONE
!
#include "parvars.h"

      INTEGER                                                           &
     &    vx,vy,vz,                                                     &
                                            ! IN size of fieldin
     &    lenout,                                                       &
                                            ! IN size of fieldout
     &    GR,                                                           &
                                            ! IN ppxref gridtype code
     &    st_grid,                                                      &
                                            ! IN STASH gridtype code
     &    fld_type,                                                     &
                                            ! IN field type (u/v/p)
     &    halo_type,                                                    &
                                            ! IN halo type
     &    halo_x,                                                       &
                                            ! IN EW halo of input
     &    halo_y,                                                       &
                                            ! IN NS halo of input
     &    n_rows_out,                                                   &
                                            ! OUT no. of output rows
     &    n_cols_out,                                                   &
                                            ! OUT no. of output cols
     &    this_index_lev,                                               &
                                            ! IN level index, this field
     &    row_length,rows,n_rows,                                       &
                                            ! IN horiz. sizes (C grid)
     &    no_rows,                                                      &
                                            ! IN number of rows used
     &    model_levels,                                                 &
                                            ! IN vertical size
     &    control_size,                                                 &
                                            ! IN size of control record
     &    control(control_size),                                        &
                                            ! IN control record
     &    icode,                                                        &
                                            ! OUT error code 0 if ok
     &    no_of_levels,                                                 &
                                            ! IN no of levels
     &    no_of_levels_masswt,                                          &
                                ! IN levels for mass weighting array
                                ! lmasswt F: =1; T: =no_of_levels
     &    index_lev(no_of_levels),                                      &
                                            ! IN index to levels
     &    level_list(no_of_levels)          ! IN model level list
      REAL                                                              &
     &    fieldin(vx,vy,vz),                                            &
                             ! IN fieldin which is to be acted on
     &    p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
                                            ! IN pressure (rho levels)
     &    pstar(row_length,rows),                                       &
                                            ! IN surface pressure
     &    cos_v_latitude(row_length,n_rows),                            &
                                            ! IN v-grid area fn
     &    cos_theta_latitude(row_length,rows),                          &
                                                ! IN T-grid area fn
     &    fieldout(lenout),                                             &
                                                ! OUT output field
     &    rmdi                                  ! IN  missing data indic
      LOGICAL                                                           &
     &    lcyclic,                                                      &
                                                ! IN .true. if cyclic EW
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    land(row_length,rows),                                        &
                                                ! IN land mask
     &    sea(row_length,rows)                  ! IN sea mask
      CHARACTER*(*) cmessage                    ! OUT error message

!*----------------------------------------------------------------------
!
#include "stparam.h"
#include "sterr.h"
#include "cppxref.h"
!L
!L external routines
!L
      EXTERNAL stextc ! extracts the field
      EXTERNAL stcolm ! computes the column mean
      EXTERNAL stzonm ! computes the zonal mean
      EXTERNAL stmerm ! computes the meridional mean
      EXTERNAL stglom ! computes the global mean
      EXTERNAL stfieldm ! computes the field mean
!L
!L local variables
!L
      LOGICAL lwrap                ! TRUE if output field wraparound EW
      LOGICAL lmasswt_strict       ! copy of lmasswt - but set to false
!                                  ! if mass weighting is not requested
      INTEGER xstart,ystart        ! lower left hand corner coords
      INTEGER xend,yend            ! upper right hand corner coords
      INTEGER processing_code      ! what kind of mean  will be done
      INTEGER what_level           ! what type of input level
      INTEGER what_mask            ! what mask is used
      INTEGER what_weight          ! what weighting is used

      INTEGER i,j,k                                                     &
                                   ! loop counters
     &,model_level                                                      &
                                   ! model level
     &,this_index_lev_masswt       ! level index for mass weighting
                                   ! (=1 if no mass weighting or no
                                   !  model level weights available)

      INTEGER                                                           &
! global versions of the extracted area domain limits
     &  global_xstart,global_xend,global_ystart,global_yend

! workspace arrays containining weighting factors and masks.
      REAL                                                              &
     &  area_weight(1-halo_x:row_length+halo_x+1,                       &
     &              1-halo_y:no_rows+halo_y)                            &
     &, pstar_weight(1-halo_x:row_length+halo_x+1,                      &
     &               1-halo_y:no_rows+halo_y,                           &
     &               no_of_levels_masswt)                               &
     &, pstar_interp(row_length,no_rows)
      LOGICAL                                                           &
     &  mask(1-halo_x:row_length+halo_x+1,                              &
     &       1-halo_y:no_rows+halo_y)


!L----------------------------------------------------------------------
!L 1. Set up local variables
!L
      xstart=control(st_west_code)
      xend=control(st_east_code)
      ystart=control(st_south_code)  ! NOTE: Grid is assumed to be
      yend=control(st_north_code)    !       oriented south-to-north

      global_xstart=xstart
      global_ystart=ystart
      global_xend=xend
      global_yend=yend

! and calculate what the local subdomain limits are:
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE.,.TRUE.,                  &
     &                                  GR,halo_type,mype,              &
     &                                  global_ystart,global_xend,      &
     &                                  global_yend,global_xstart,      &
     &                                  ystart,xend,yend,xstart)

! Check if wraparound field
      IF (xstart >  xend) THEN

        IF (lcyclic) THEN
          xend=xend+blsize(1,fld_type)
! subtract two halos as we don't wish to include halos at the end
! and start of the row within the wrap around domain
          lwrap=.TRUE.
        ELSE
          icode=st_bad_wraparound    ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF

      ELSE
        lwrap=.FALSE.
      ENDIF

      IF (global_xstart  >   global_xend) THEN
        IF (lcyclic) THEN
          global_xend=global_xend+glsize(1,fld_type)
        ELSE
          icode=st_bad_wraparound  ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF
      ENDIF

      processing_code=control(st_gridpoint_code)
      what_level=control(st_input_bottom)
      what_mask=mod(processing_code,block_size)
      what_weight=control(st_weight_code)
!L
!L 1.1 Prevent masking or weighting if input field is not 2D in extent
!L     - weighting and masking is assumed to have been done outside.
!L
      IF ( (.NOT.(st_grid == st_tp_grid.OR.st_grid == st_uv_grid.OR.    &
     &            st_grid == st_cu_grid.OR.st_grid == st_cv_grid.OR.    &
     &             st_grid == st_riv_grid))                             &
     &      .AND.(what_mask   /= stash_null_mask_code   .OR.            &
     &            what_weight /= stash_weight_null_code) ) THEN
       icode=st_not_supported
       cmessage='SPATIAL : Masking/weighting unsupported - non 2D field'
       GOTO 999
      ENDIF

! Check for supported weighting and masking options

      IF (.NOT. ((what_weight  ==  stash_weight_null_code) .OR.         &
     &           (what_weight  ==  stash_weight_area_code) .OR.         &
     &           (what_weight  ==  stash_weight_volume_code) .OR.       &
     &           (what_weight  ==  stash_weight_mass_code) ) ) THEN
        cmessage='SPATIAL : Unrecognized weighting option'
        icode=unknown_weight
        GOTO 999
      ENDIF

      IF (.NOT. ((what_mask  ==  stash_null_mask_code) .OR.             &
     &           (what_mask  ==  stash_land_mask_code) .OR.             &
     &           (what_mask  ==  stash_sea_mask_code ) ) ) THEN
        cmessage='SPATIAL : Unrecognized masking option'
        icode=unknown_mask
        GOTO 999
      ENDIF

      IF (what_weight  ==  stash_weight_volume_code) THEN
        cmessage='SPATIAL : Volume-weighting not supported'
        icode=st_illegal_weight
        GOTO 999
      ENDIF

! Set lmasswt_strict - copy of lmasswt, but set to false is mass
! weighting not requested

      lmasswt_strict=                                                   &
     &  (lmasswt .AND. (what_weight  ==  stash_weight_mass_code))

! Precalculate weighting and masking arrays
! I've used IF tests inside the loops, but since the logical
! expressions are invariant wrt i and j, the compiler will
! move them outside the DO loops. It makes the code a lot shorter!

! NOTE that neither area-weights or mass-weights are normalised so
! that the interpretation of weighted diagnostics is non-obvious. Also
! the PP lookup header has no switch for indicating whether or not such
! processing has taken place. More careful treatment of horizontal
! interpolation is required for stricter accuracy.


! area weighting
      IF (what_weight  ==  stash_weight_null_code) THEN
          ! Ensure initialisation of weight array including halos
          area_weight (:,:)   = 1.0
      ELSE ! some form of area weighting will be required
          DO j=1,no_rows
          DO i=1,row_length
             IF (st_grid  ==  st_cv_grid) THEN
                 area_weight(i,j)=cos_v_latitude(i,j)
             ELSE
! NOTE that this is will not be accurate for C-u grid variables for
! LAMs, since cos_theta_latitude will differ between theta,u positions
                 area_weight(i,j)=cos_theta_latitude(i,j)
             ENDIF
         ENDDO ! i
         ENDDO ! j

      ENDIF    ! what_weight


! mass weighting
      IF ((what_weight  ==  stash_weight_null_code) .OR.                &
     &    (what_weight  ==  stash_weight_area_code)) THEN
! No mass weighting is required
          this_index_lev_masswt = 1

          ! Ensure initialisation of weight array including halos
          pstar_weight(:,:,this_index_lev_masswt) = 1.0
      ELSE

! Mass weighting requested

! Ensure that halos are initialised
        DO j=no_rows-1,no_rows
        DO i=1,row_length
           pstar_interp(i,j) =1.0
        ENDDO !i
        ENDDO !j

! Interpolate pstar to correct horizontal staggering
        IF     (st_grid  ==  st_cu_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_U(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSEIF (st_grid  ==  st_cv_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_V(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSE
          DO j=1,no_rows
            DO i=1,row_length
              pstar_interp(i,j)=pstar(i,j)
            ENDDO
          ENDDO
        ENDIF

        IF(lmasswt) THEN  ! model level mass weights available

          this_index_lev_masswt = this_index_lev

          DO k=1,no_of_levels_masswt
            model_level = level_list(k)
            IF(model_level == model_levels) THEN  ! top level

              DO j=1,no_rows
              DO i=1,row_length
                pstar_weight(i,j,k) = p(i,j,model_level)
              ENDDO !i
              ENDDO !j

            ELSE      ! not top level

              DO j=1,no_rows
              DO i=1,row_length
! Only accurate for variables on theta levels
                pstar_weight(i,j,k) = p(i,j,model_level) -              &
     &                                  p(i,j,model_level+1)
              ENDDO !i
              ENDDO !j
            ENDIF

          ENDDO ! k

        ELSE              ! no model level mass weights available:
                          ! weight by pstar

            this_index_lev_masswt = 1

            DO j=1,no_rows
            DO i=1,row_length
               pstar_weight(i,j,this_index_lev_masswt) =                &
     &                                               pstar_interp(i,j)
            ENDDO !i
            ENDDO !j
! Horizontal interpolation may be required - this should be done here:

        ENDIF   ! lmasswt
      ENDIF

! masking

      DO j=1,no_rows
        DO i=1,row_length
          IF (what_mask  ==  stash_land_mask_code) THEN
            mask(i,j)=land(i,j)
          ELSEIF (what_mask  ==  stash_sea_mask_code) THEN
            mask(i,j)=sea(i,j)
          ELSE
            mask(i,j)=.TRUE.
          ENDIF
        ENDDO
      ENDDO

! Update the halo areas of the weighting/masking arrays

!  [Note that for lams at UM5.2 and before, valid values for rim
!   boundaries would be needed using FILL_EXTERNAL_HALOS calls for
!   area_weight and pstar_weight arrays, but this should now be
!   superseded by initialising full arrays]

      ! Update halos only if halos present (standard diagnostics have no
      ! halos, and if some weighting required (probably defunct
      ! functionality)
      IF( (halo_x  /=  0 .OR. halo_y  /=  0) .AND.                      &
     &    (what_weight  /=  stash_weight_null_code .OR.                 &
     &     what_mask    /=  stash_null_mask_code )     ) THEN

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( AREA_WEIGHT(1-halo_x:row_length+halo_x,:),      &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( PSTAR_WEIGHT(1-halo_x:row_length+halo_x,:,:),   &
     &     ROW_LENGTH,NO_ROWS,                                          &
     &                no_of_levels_masswt,halo_x,halo_y,                &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( MASK(1-halo_x:row_length+halo_x,:),             &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
      ENDIF
!L----------------------------------------------------------------------
!L 2. Call service routine to perform required processing
!L
!L 2.1 Extract sub field (single level at a time)
!L
      IF (processing_code <  extract_top.and.                           &
     &    processing_code >= extract_base) THEN
        n_rows_out=(yend+1)-ystart
        n_cols_out=(xend+1)-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stextc
        CALL STEXTC(fieldin,vx,vy,fld_type,halo_type,                   &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.2 Calculate column mean (over multiple levels indexed by index_lev)
!L
      ELSEIF (processing_code <  vert_mean_top.and.                     &
     &        processing_code >  vert_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=xend+1-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stcolm
        CALL STCOLM(fieldin,vx,vy,vz,fld_type,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.3 Calculate zonal mean (single level at a time)
!L
      ELSEIF (processing_code <  zonal_mean_top.and.                    &
     &        processing_code >  zonal_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=1
! DEPENDS ON: stzonm
        CALL STZONM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.4 Calculate meridional mean (single level at a time)
!L
      ELSEIF (processing_code <  merid_mean_top.and.                    &
     &        processing_code >  merid_mean_base) THEN
        n_rows_out=1
        n_cols_out=xend+1-xstart
! DEPENDS ON: stmerm
        CALL STMERM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.5 Calculate field mean (single level at a time)
!L
      ELSEIF (processing_code <  field_mean_top.and.                    &
     &        processing_code >  field_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stfieldm
        CALL STFIELDM(fieldin,vx,vy,fld_type,gr,halo_type,              &
     &                lwrap,lmasswt_strict,                             &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.6 Calculate global mean (over multiple levels)
!L
      ELSEIF (processing_code <  global_mean_top.and.                   &
     &        processing_code >  global_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stglom
        CALL STGLOM(fieldin,vx,vy,vz,fld_type,gr,halo_type,             &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.7 Invalid processing option
!L
      ELSE
        icode=unknown_processing
        write(cmessage,111)'unknown processing option',                 &
     &    processing_code
      ENDIF
!L
  999 CONTINUE
111   format('SPATIAL : >>FATAL ERROR <<',a40,i5)
!
      RETURN
      END SUBROUTINE SPATIAL
#endif
