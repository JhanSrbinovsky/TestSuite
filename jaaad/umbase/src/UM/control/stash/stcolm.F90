#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STCOLM ---------------------------------------------------
!LL
!LL  Purpose: Calculate weighted column mean within a region specified
!LL           by a lower left hand and upper right hand corner.
!LL           (STASH service routine).
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.3  16/09/93  Allow mass-weighting only if level type has a well-
!LL                  defined mass-weight (ie. model levels/half-levels).
!LL   4.3  06/01/97  Moved calculation of  weighting and masking arrays
!LL                   up to SPATIAL.                          P.Burton
!LL   5.0  13/07/99  Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0  17/11/99  Removed row_length,p_rows arguments
!LL                                             P.Burton
!     6.2   18/10/05 Fix bugs introduced by gan1f601      Andy Malcolm
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D712
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LLEND ---------------------------------------------------------------
!
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STCOLM(fieldin,vx,vy,vz,fld_type,halo_type,            &
     &                  lwrap,lmasswt,                                  &
     &                  xstart,ystart,xend,yend,                        &
     &                  fieldout,index_lev,zsize,                       &
     &                  pstar_weight,                                   &
     &                  area_weight,mask,                               &
     &                  level_code,mask_code,weight_code,rmdi,          &
     &                  icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &    vx,vy,vz,                                                     &
                                                ! IN  input field size
     &    fld_type,                                                     &
                                                ! IN  field type(u/v/p)
     &    halo_type,                                                    &
                                                ! IN  halo type
     &    xstart,ystart,                                                &
                                                ! IN  lower LH corner
     &    xend,yend,                                                    &
                                                ! IN  upper RH corner
     &    zsize,                                                        &
                                      ! IN no of horiz levels to process
     &    index_lev(zsize),                                             &
                                      ! IN offset for each horiz level
     &    level_code,                                                   &
                                                ! IN  input level code
     &    mask_code,                                                    &
                                                ! IN  masking code
     &    weight_code,                                                  &
                                                ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER*(*)                                                     &
     &    cmessage                              ! OUT error return msg
      LOGICAL                                                           &
     &    lwrap,                                                        &
                                                ! IN  TRUE if wraparound
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    mask(vx+1,vy)                         ! IN  mask array
      REAL                                                              &
     &    fieldin(vx,vy,vz),                                            &
                                                ! IN  input field
     &    fieldout(xstart:xend,ystart:yend),                            &
                                                ! OUT output field
     &    pstar_weight(vx+1,vy,zsize),                                  &
                                                ! IN  mass weight factor
     &    area_weight(vx+1,vy),                                         &
                                                ! IN  area weight factor
     &    rmdi                                  ! IN  missing data indic
!*----------------------------------------------------------------------
!
! External subroutines called
!
!
#include "stparam.h"
#include "sterr.h"
#include "parvars.h"
!
! Local variables
!
        INTEGER i,j,k,ii,kk       ! ARRAY INDICES FOR VARIABLE

        REAL SUMCTOP(xstart:xend,ystart:yend)
        REAL SUMCBOT(xstart:xend,ystart:yend)

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
      DO 1 i=xstart,xend
      DO j=ystart,yend
        SUMCTOP(i,j)=0.0
        SUMCBOT(i,j)=0.0
      ENDDO
 1    CONTINUE
!L----------------------------------------------------------------------
!L 1. Form column sums
!L
!L 1.1 NULL weighting or area weighting
!L
      IF (weight_code == stash_weight_null_code.or.                     &
     &    weight_code == stash_weight_area_code) THEN
        DO 11 kk=1,zsize
          k=index_lev(kk)
          DO i=xstart,xend

            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
            ELSE
              ii=i
            ENDIF

            DO j=ystart,yend
              SUMCBOT(i,j)=SUMCBOT(i,j)+1.0
              SUMCTOP(i,j)=SUMCTOP(i,j)+fieldin(ii,j,k)
            ENDDO
          ENDDO
 11     CONTINUE
!L
!L 1.2 mass weighting
!L
      ELSEIF (weight_code == stash_weight_mass_code) THEN
        IF (.NOT.lmasswt) THEN
! Mass-weighting on level types with no mass-weight defined is not
! supported - should be prevented by UI
          cmessage='STCOLM  : mass-weights not defined for this diag'
          icode=st_illegal_weight
          goto 999
        ELSE
          DO 121 kk=1,zsize
            k=index_lev(kk)
            DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
              ELSE
                ii=i
              ENDIF
              DO j=ystart,yend
                SUMCBOT(i,j)=SUMCBOT(i,j) +  pstar_weight(ii,j,kk)
                SUMCTOP(i,j)=SUMCTOP(i,j) +  fieldin(ii,j,k)*           &
     &                                       pstar_weight(ii,j,kk)
              ENDDO
            ENDDO
 121      CONTINUE
        ENDIF
      ELSE
        cmessage='STCOLM  : Invalid weighting code detected'
        icode=unknown_weight
        goto 999
      ENDIF
!L----------------------------------------------------------------------
!L 2. Perform masking (set missing data at masked points) - compute mean
!L
      DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
        ELSE
          ii=i
        ENDIF
        DO j=ystart,yend
          IF (mask(ii,j)) THEN
            fieldout(i,j)=SUMCTOP(i,j)/SUMCBOT(i,j)
          ELSE
            fieldout(i,j)=rmdi
          ENDIF
        ENDDO
      ENDDO
!L
  999 CONTINUE
      RETURN
      END SUBROUTINE STCOLM
#endif
