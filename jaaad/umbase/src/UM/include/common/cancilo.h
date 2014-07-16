!*L--------------------COMDECK  CANCILO ---------------------------
!
!Contains PARAMETERS, headers, and index blocks for control of update
!of ancillary fields.
!System component F0171 Parameters determined by model version and
!size parameters.
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  04/10/94  Increase NANCIL_FIELDS from 23 to 37. RTHBarnes
!  3.4  04/08/94  LEVELS array included. Mike Bell.
!
      INTEGER                                                           &
     &  NANCIL_FIELDS,                                                  &
                         ! Maximum total number of ancillary fields
     &  FILEANCIL,                                                      &
                         ! File number associated with ancillary fields
     &  NLOOKUP,                                                        &
                         ! Position of given ancillary field in lookup
!                        ! tables ( Set by INANCCTL from parameters,
!                        ! and CSIZEATM AND CSIZEOCN).
     &  LOOKUP_STEP,                                                    &
                         ! Interval between PP Headers refering to
!                        ! to the same ancillary fields at diferent time
     &  LEVELS,                                                         &
                         ! Number of levels of data in each ancillary
!                        ! field (Set by INANCILO )
     &  STASHANCIL,                                                     &
                         ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block

!PARAMETER statement, fixed by model version and information from
!CSIZEATM,CSIZEOCN

      PARAMETER                                                         &
     &          (NANCIL_FIELDS=37) ! Hard-wired, add user ancillaries.

!*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL                                                           &
     &         UPDATE
      INTEGER  FIELDCODE,                                               &
     &         STEPS
!*----------------------------------------------------------------------
      COMMON/CTANCILO/                                                  &
     &         FIELDCODE(2,NANCIL_FIELDS),                              &
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS)
      COMMON/IXANCILO/ FILEANCIL(NANCIL_FIELDS),                        &
     &           NLOOKUP(NANCIL_FIELDS),                                &
     &           LOOKUP_STEP(NANCIL_FIELDS),                            &
     &           LEVELS(NANCIL_FIELDS),                                 &
     &           STASHANCIL(NANCIL_FIELDS),                             &
     &           D1_ANCILADD(NANCIL_FIELDS)
