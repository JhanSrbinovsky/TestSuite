! CONANC start
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  04/10/94  Increase no.of datasets for murk & user ancil. RTHB
!  4.1  02/04/96  Introduce wave sub-model.  RTHBarnes.
!  4.1  26/02/96  Increase NANCIL_DATASETSA. D. Robinson.
!  4.3   18/3/97  Increase NANCIL_DATASETSA to 19.  William Ingram.
!  4.4   12/9/97  Increase NANCIL_DATASETSA to 22.  R.A.Betts
!  5.3   19/06/01 Increase NANCIL_DATASETSA to 25 (from 24) for
!                 tropopause-based ozone               Dave Tan
!  5.5   30/12/02 Increase NANCIL_DATASETSA to 28. Dave Robinson.
!  6.1   07/04/04 Increase NANCIL_DATASETSA to 29. Andy Jones
!  6.1   08/11/04 Increase NANCIL_DATASETSA to 32 for RR files. R.Sharp
!  7.3   31/05/10 Increase NANCIL_DATASETSA to 67 from 47 for
!                 tracer fluxes - kdcorbin
!
! To be used in conjunction with file TYPANC, defining dimensions
! of headers from ancillary files. A separate file is needed to
! allow calculation of super array sizes in UM_INDEX.
      INTEGER, PARAMETER :: NANCIL_DATASETSA=67
      INTEGER, PARAMETER :: NANCIL_DATASETSO=10
      INTEGER, PARAMETER :: NANCIL_DATASETSW=1
      INTEGER, PARAMETER :: NANCIL_DATASETS=NANCIL_DATASETSA+           &
     &  NANCIL_DATASETSO+NANCIL_DATASETSW
! CONANC end
