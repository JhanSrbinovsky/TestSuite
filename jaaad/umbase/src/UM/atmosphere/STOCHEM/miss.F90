#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE MISS(missings,missingbe7,missingbe10,                  &
     &                missinglnox,missingacnox,so2em,be7em,             &
     &                be10em,lnoxem,acnoxem,totso2em,totbe7em,          &
     &                totbe10em,totlnoxem,totacnoxem,nnn)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : ACCOUNT FOR GRID VOLUMES WITH NO CELLS
!-
!-   Inputs  : NO2EM,SO2EM,BE7EM,BE10EM,LNOxEM,
!-             TOTNO2EM,TOTSO2EM,TOTBE7EM,TOTBE10EM,TotLNOxEM,
!-   Outputs : MISSING,MISSINGS,MISSINGBE7,MISSINGBE10,MissingLNOx
!-   Controls:
!-
!-   Created   23-Jan-2001 W.J. Collins
!VVV  V4.6  MISS 23/1/01
!    6.2   21/10/05  Replace GSYNC with SSYNC. P.Selwood
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: nnn

      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: so2em
      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: be7em
      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: be10em
      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: lnoxem
      REAL, DIMENSION(nlnpe,nlpe,nlev), INTENT(IN) :: acnoxem
      REAL, INTENT(IN) :: totso2em
      REAL, INTENT(IN) :: totbe7em
      REAL, INTENT(IN) :: totbe10em
      REAL, INTENT(IN) :: totlnoxem
      REAL, INTENT(IN) :: totacnoxem
      REAL, INTENT(OUT) :: missings
      REAL, INTENT(OUT) :: missingbe7
      REAL, INTENT(OUT) :: missingbe10
      REAL, INTENT(OUT) :: missinglnox
      REAL, INTENT(OUT) :: missingacnox

      INTEGER :: info
      LOGICAL, DIMENSION(nlnpe,nlpe,nlev) :: mask

! Find grid boxes with no cells, increase emissions in all
! others to make up only for boxes on this PE.
! Spread just copies the 2D array to give a 3D array with
! NLEV identical levels
        mask=(nnn==0)
        missings=   SUM(so2em(:,:,:) ,MASK=mask)
        MissingLNOx=SUM(LNOxEM(:,:,:) ,MASK=MASK)
        MissingAcNOx=SUM(AcNOxEM(:,:,:) ,MASK=MASK)
        missingbe7= SUM(be7em(:,:,:) ,MASK=mask)
        missingbe10=SUM(be10em(:,:,:),MASK=mask)
!       Sum all the `MISSING's over all processors
        info=GC_SHM_PUT
        CALL GC_SSYNC(nproc,info)
        CALL GC_RSUM(1,nproc,info,missings)
        CALL GC_RSUM(1,nproc,info,missingbe7)
        CALL GC_RSUM(1,nproc,info,missingbe10)
        CALL GC_RSUM(1,NPROC,INFO,MissingLNOx)
        CALL GC_RSUM(1,NPROC,INFO,MissingAcNOx)
! TotLNOxEM must be added up over all PEs here
        CALL GC_RSUM(1,NPROC,INFO,TotLNOxEM)
        CALL GC_SSYNC(nproc,info)

        IF(TotLNOxEM>0) MissingLNOx=TotLNOxEM/( TotLNOxEM- MissingLNOx)
        IF(TotAcNOxEM>0) MissingAcNOx=TotAcNOxEM/                       &
     &                   ( TotAcNOxEM- MissingAcNOx)
        IF(totso2em>0)  missings=   totso2em/( totso2em- missings)
        IF(totbe7em>0)  missingbe7= totbe7em/( totbe7em- missingbe7)
        IF(totbe10em>0) missingbe10=totbe10em/(totbe10em-missingbe10)


      END SUBROUTINE MISS
#endif
