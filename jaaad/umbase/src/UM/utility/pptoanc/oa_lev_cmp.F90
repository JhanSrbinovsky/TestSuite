#if defined(PPTOANC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routines: OA_PACK, OA_UNPACK and OA_LEV_CMP
!LL
!LL  Author: M. J. Bell   Date: 7 January 1992
!LL
!LL  Logical components covered:
!LL
!LL  Tested under compiler: VAX: Vax Fortran 5.4   Cray: cf77
!LL  Tested under OS version: VAX: VAX/VMS V5.5   Cray: Unicos 5.0
!LL
!LL  Programming standard : FOAM Doc Paper 3/2/1
!LL
      SUBROUTINE OA_LEV_CMP(NO_ROWS_M, NO_COLS_M, NO_LEVS_M             &
     &      ,NO_SEG, INDX_CMP, INDX_TO_ROWS, NO_CMP, KLEV               &
     &      ,KLEV_OFF_CMP, NO_CMP_KLEV, KLEV_OFF_UV, NO_UV_KLEV)
!*
!LL  Purpose:  This subroutine calculates the index of the first value
!LL            and the number of points in the 3D compressed array
!LL            or 3D uncompressed array for a given level
!
      IMPLICIT NONE
!
!*L  ARGUMENT LIST
!
! Logicals
! Dimensions
      INTEGER NO_ROWS_M  ! IN number of rows on model grid
      INTEGER NO_COLS_M  ! IN number of columns on model grid
      INTEGER NO_LEVS_M  ! IN number of levels on model grid
      INTEGER NO_SEG     ! IN number of sea segments in 3D compressed fi
! Compression and expansion indices
      INTEGER INDX_CMP(NO_SEG) ! IN index in compressed array of start o
!                                   each sea segment
      INTEGER INDX_TO_ROWS(NO_ROWS_M*NO_LEVS_M) ! IN index of first/next
!                                   sea segment for each row and level
      INTEGER NO_CMP           ! IN total no of points in a 3D
! Level
      INTEGER KLEV             ! IN model level
! Output
      INTEGER KLEV_OFF_CMP ! OUT offset of 1st compressed point on this
!                                level from 1st compressed point
      INTEGER NO_CMP_KLEV  ! OUT number of compressed points on this lev
      INTEGER KLEV_OFF_UV  ! OUT offset to u-v grid data on this level
      INTEGER NO_UV_KLEV   ! OUT no. of u-v grid points on this level
!*
!    LIST OF OTHER VARIABLES
      INTEGER ISEG_ST     ! index for first segment on this level
      INTEGER ISEG_NXT    ! index for first segment on next level
!*
!-----------------------------------------------------------------------
!

        NO_CMP_KLEV = NO_ROWS_M * NO_COLS_M
        KLEV_OFF_CMP = (KLEV-1) * NO_CMP_KLEV
#if defined(MPP)
        ! u and v fields have same no. of rows as t field in MPP mode
        NO_UV_KLEV = NO_ROWS_M * NO_COLS_M
#else
        ! u and v fields have 1 row less than t field in non-MPP mode
        NO_UV_KLEV = (NO_ROWS_M-1) * NO_COLS_M
#endif
        KLEV_OFF_UV = (KLEV-1) * NO_UV_KLEV


      RETURN
      END SUBROUTINE OA_LEV_CMP
!
!-----------------------------------------------------------------------
#endif
