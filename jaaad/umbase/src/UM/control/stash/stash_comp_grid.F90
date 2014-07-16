#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE STASH_COMP_GRID----------------------
!LL
!LL  Compute grid descriptors
!LL  out_bdx,out_bdy,out_bzx and out_bzy from inputs
!LL  1)  samples
!LL  2)  grid type code
!LL  3)  ocean (true if an ocean grid)
!LL  4)  real and integer headers
!LL  5)  w_mostcol and n_mostrow
!LL  6) processing code
!LL
!LL  Tested under compiler CFT77
!LL  Tested under OS version 6.1
!LL
!LL  Author Simon Tett (Based on  code in pp_head by Tim Johns)
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!    5.4    11/04/02 Set up details for variable horizontal
!                    grids.                         R. Hill
!    6.0  01/10/03 Correct header values for river fields. C.Bunton
!    6.2  12/07/05 Correct header values for timeseries. R A Stratton
!LL
!LL  Logical components covered: D40
!LL
!LL  Project TASK: C4
!LL
!LL  Programming standard: U M DOC  Paper NO. 4,
!LL
!LL  External documentation  C4
!LL
!LLEND-------------------------------------------------------------

!
!*L  INTERFACE and ARGUMENTS:------------------------------------------
      SUBROUTINE STASH_COMP_GRID(                                       &
     &  out_bzx,out_bzy,out_bdx,out_bdy,                                &
     &  samples,st_grid,ocean,                                          &
     &  w_mostcol,n_mostrow,                                            &
     &  realhd,len_realhd,inthd,len_inthd,gr,                           &
     &  ICODE,CMESSAGE)
!
      IMPLICIT NONE


!
      LOGICAL  OCEAN          !IN TRUE if processing an ocean diagnostic
!
      INTEGER samples                ! IN no of samples in period (times
!
      INTEGER                                                           &
     &  ICODE                                                           &
                          !OUT   Return code from the routine
     &, LEN_REALHD                                                      &
                          !IN    Length of the Real Constants
     &, LEN_INTHD                                                       &
                          !IN    Length of integer constants
     &, INTHD(LEN_INTHD)  !IN    Integer constants
!
      INTEGER                                                           &
     &  st_grid                                                         &
                          !IN    STASH horizontal grid type
     &, N_MOSTROW                                                       &
                          !IN    The most nrthrly row.
     &, W_MOSTCOL                                                       &
                          !IN    The most westerly column
     &, gr                !IN    The type of processing done
!
      REAL                                                              &
     &  REALHD(LEN_REALHD)  !IN  Real header
      REAL                                                              &
     &  out_bzx,out_bdx,out_bzy,out_bdy  ! OUT grid descriptors
      CHARACTER*(*)                                                     &
     &  cmessage           ! OUT error messages
!*
!
!L Local Variables
!
      INTEGER mean_code
!
!*---------------------------------------------------------------------
#include "clookadd.h"
#include "stparam.h"
#include "cppxref.h"
#include "parvars.h"
#include "comvgrid.h"
!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: None
!
!*---------------------------------------------------------------------
!
!LL   Construct PP header
      IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
        OUT_BZX=0.0
        OUT_BDX=0.0
        OUT_BZY=0.0
        OUT_BDY=0.0
      ELSE
        IF (OCEAN) THEN       !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX fo
          IF (st_grid == st_uv_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)/2.0
            OUT_BZX=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid == st_tp_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)
            OUT_BZX=REALHD(4)-REALHD(1)
          ELSEIF (st_grid == st_cu_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)
            OUT_BZX=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid == st_cv_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)/2.0
            OUT_BZX=REALHD(4)-REALHD(1)
          ENDIF
          IF (REALHD(32) >  REALHD(29)) THEN !   greater than RMDI
            OUT_BDY=0.0
            OUT_BDX=REALHD(32)
          ELSE
            OUT_BDY=REALHD(2)
            OUT_BDX=REALHD(1)
          ENDIF
      ! If this is variable horizontal grid information
      ! then the horizontal spacings in the headers must be
      ! set to zero in the appropriate direction.
      IF (VAR_GRID_TYPE >  0) THEN
         IF (X_VAR_GRID) OUT_BDX = 0.0
         IF (Y_VAR_GRID) OUT_BDY = 0.0
      ENDIF
        ELSE                 !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX for
          IF(st_grid == st_riv_grid)THEN
            OUT_BDY = -180.0/180
            OUT_BZY = -REALHD(3) - OUT_BDY*0.5
            OUT_BDX = 360.0/360
            OUT_BZX = REALHD(4) - OUT_BDX*0.5
          ELSE
           IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)/2.0 ! UV pts
           ELSE
            OUT_BZY=REALHD(3)-REALHD(2) ! Zeroth Lat OUT_BZY
           ENDIF
!
           IF(st_grid == st_uv_grid.OR.st_grid == st_cu_grid) THEN
             OUT_BZX=REALHD(4)-REALHD(1)/2.0 !UV points
           ELSE
             OUT_BZX=REALHD(4)-REALHD(1) ! Zeroth Long OUT_BZX
           ENDIF
           OUT_BDX=REALHD(1) ! Long intvl OUT_BDX
          OUT_BDY=REALHD(2) ! Lat intvl OUT_BDY

          ENDIF
        ENDIF
!
! Add on offset for fields not starting from the origin
!
      ! If this is variable horizontal grid information
      ! then the horizontal start points are NOT equally spaced
      ! we must get our actual start point some other way!
      ! This is all theoretical since this routine currently
      ! gets called with hard coded values of N_MOSTROW
      ! and W_MOSTCOL. Furthermore, the use of N_MOSTROW
      ! is entirely wrong in this routine!
      ! It should be the southern-most row for this to do
      ! anything useful!
      IF (VAR_GRID_TYPE >  0) THEN
         IF (X_VAR_GRID) OUT_BZX = X_BOUNDARY(W_MOSTCOL,VAR_GRID_TYPE)
         IF (Y_VAR_GRID) OUT_BZY = Y_BOUNDARY(N_MOSTROW,VAR_GRID_TYPE)
      ELSE
        OUT_BZY=OUT_BZY                                                 &
     &       +(N_MOSTROW-1)*OUT_BDY
        OUT_BZX=OUT_BZX                                                 &
     &     +(W_MOSTCOL-1)*OUT_BDX

      ENDIF
        IF(OUT_BZX >= 360.0)                                            &
     &     OUT_BZX=OUT_BZX-360.0
!
! If horizontal averaging has been applied to the output field,
! set OUT_BDX and/or OUT_BDY to the full domain extent
!
        mean_code=(GR/block_size)*block_size
        IF (mean_code == zonal_mean_base .OR.                           &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          OUT_BDX=REAL(INTHD(6))*REALHD(1)
        ENDIF
        IF (mean_code == merid_mean_base .OR.                           &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          OUT_BDY=REAL(INTHD(7))*REALHD(2)
        ENDIF
      ENDIF
!
  999 CONTINUE
      RETURN
      END SUBROUTINE STASH_COMP_GRID
#endif
