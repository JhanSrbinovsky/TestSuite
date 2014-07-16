#if defined(OASIS3) || defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE MEAN_POLAR_ROW(field_local)

  IMPLICIT NONE
  !
  ! Descrition : This routine calculates the mean value of all grid points
  !              on the North polar row of an incoming field
  !              which may be decomposed over all PEs. The mean value
  !              is then written to all points on the polar row in order
  !              to ensure only one value is present at the singularity.
  !              This routine will only work for fields on p points of
  !              the C grid (i.e. for points defined actually ON the pole).
  !
  ! Author:      R. Hill
  !=======================================================================
#include "parvars.h"

  ! Incoming 2D field
  REAL :: field_local(lasize(1,fld_type_p,halo_type_no_halo), &
                      lasize(2,fld_type_p,halo_type_no_halo))

  ! Local fields
  REAL :: mean

  INTEGER :: info

  ! 1st: Employ RVECSUMR in working out a full polar row sum for
  ! polar averaging purposes - it should be faster than a full
  ! gather/scatter.

  CALL GC_SSYNC ( nproc, info )

  ! We only need to involve the nortern-most PEs
  IF (at_extremity(PNorth)) THEN
    mean = 0.0 ! initialise mean value

    CALL GCG_RVECSUMR (lasize(1,fld_type_p,halo_type_no_halo),    &
       lasize(1,fld_type_p,halo_type_no_halo),                    &
       1,                                                         &
       1,                                                         &
       field_local(:,lasize(2,fld_type_p,halo_type_no_halo)),     &
       gc_proc_row_group,                                         &
       info,                                                      &
       mean)

    ! 2nd: Do the meaning on the polar row
    mean=mean/glsize(1,fld_type_p)

    ! Copy our meaned field into the polar row
    field_local(1:lasize(1,fld_type_p,halo_type_no_halo),     &
       lasize(2,fld_type_p,halo_type_no_halo)) = mean

  END IF ! at_extremity(PNorth)=.true.


END SUBROUTINE MEAN_POLAR_ROW
#endif
