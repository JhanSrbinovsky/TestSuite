#if defined(A19_1A) || defined(A19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Initialises accumulated carbon fluxes to zero if new calling period
!
! Subroutine Interface:
      SUBROUTINE INIT_ACC(LAND_PTS,                                     &
     &       NPP_PFT_ACC,G_LEAF_PHEN_PFT_ACC,                           &
     &       RESP_W_PFT_ACC,RESP_S_ACC,ICODE,CMESSAGE)

      IMPLICIT NONE
!
! Description:
!   Resets accumulation prognostics to zero if a new TRIFFID calling
!   period is starting.  This routine is needed when starting an NRUN
!   from an initial dump created in either of the following situations:
!
!   i)  Initial dump created from a non-TRIFFID run
!
!   ii) Initial dump created in a TRIFFID run mid-way through a TRIFFID
!       calling period.  The NRUN may re-start at the same point within
!       this calling period and continue with the accumulation already
!       part-completed in this dump; in this case this routine will not
!       be used.  Alternatively, the NRUN may start a new calling
!       period, in which case the accumulation must begin; this routine
!       allows this by re-setting the relevant prognostics to zero.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4   10/10/97   Original code.  Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics   M. Best
!  6.2  01/03/06  initialise all 4 respiration pools.       C.D. Jones
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Arguments

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_PTS                            ! IN number of land points


      REAL                                                              &
     & NPP_PFT_ACC(LAND_PTS,NPFT)                                       &
                                         !INOUT Accumulated NPP on
!                                        !      Plant Functional Types
     &,G_LEAF_PHEN_PFT_ACC(LAND_PTS,NPFT)                               &
                                         !INOUT Accum. phenological
!                                        !      leaf turnover rate PFTs
     &,RESP_W_PFT_ACC(LAND_PTS,NPFT)                                    &
                                         !INOUT Accumulated wood
!                                        !      respiration on PFTs
     &,RESP_S_ACC(LAND_PTS,4)         !INOUT


      INTEGER                                                           &
     & L                                                                &
                               ! Loop counter for land points
     &,N                       ! Loop counter for plant functional types

      INTEGER ICODE            ! Work - Internal return code
      CHARACTER*80 CMESSAGE    ! Work - Internal error message


      WRITE (6,*)                                                       &
     & 'INIT_ACC: setting accumulation prognostics to zero'

      DO L=1,LAND_PTS
        DO N=1,NPFT
          NPP_PFT_ACC(L,N) = 0.0
          G_LEAF_PHEN_PFT_ACC(L,N) = 0.0
          RESP_W_PFT_ACC(L,N) = 0.0
        ENDDO
      ENDDO
      DO N=1,4
        DO L=1,LAND_PTS
          RESP_S_ACC(L,N) = 0.0
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE INIT_ACC
#endif
