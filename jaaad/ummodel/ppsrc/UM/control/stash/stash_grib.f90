
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PP2GRIB------------------------------------------------
!LL
!LL  Purpose:
!LL   to code pp_header and un-packed data into grib
!LL
!LL  Written by G.Ross/ P.Smith
!LL
!LL  Model            Modification history from model version 3.3:
!LL version  Date
!LL   3.4   6/10/94 : Correct so that able to encode data other
!LL                   than just CF (m08) fields ie climate data.
!LL                   Also return error code and message.
!LL   3.4   2/12/94 : Extra argument introduced in subroutine
!LL                   CODER. MSG_LVL set to 2 ie Errors only
!     4.0  20/01/95 : Further changes to improve encoding of climate
!                     fields in grib & correct errors. (R. A. Stratton)
!     4.0  23/03/95 : Allow alternative packing method to be used for
!                     ppxref profile 6. (R.A.Stratton)
!     4.5  20/03/98   Correction for year 2K.
!                     Author D.M. Goddard
!     5.5  25/04/03   Grib data format not supported on non-CRAY
!                     platform                         E.Leung
!     6.1  12/09/03   Re-enable for FLDOP and FLDMOD on non-CRAY
!                     platforms.                       P.Dando
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL
!LL  System component:
!LL
!LL  System task:
!LL
!LL  Documentation:
!LL
!LLEND---------------------------------------------------------
!*L Arguments:-------------------------------------------------
!LL  SUBROUTINE GRIB_STASH---------------------------------------------
!LL
!LL  Purpose:
!LL   GRIB_STASH is a subroutine to indentify the stash parameter
!LL   value and section number from the grib header codes
!LL
!LL   octet 4 of the grib product definition section is the version
!LL   number of the table 2 (parameter code description)
!LL   values from 128 to 254 are available for local use, and we
!LL   use them to describe the stash section number of the field. for
!LL   each stash section number there are two octet 4 values. the first
!LL   is for stash parameter values from 0 to 255, the second for values
!LL   256 to 511.
!LL   octet 9 is the code value in table 2, ie stash parameter value, or
!LL   stash parameter value -256 if it is more than 255.
!LL
!LL  Written by G.Ross/ P.Smith
!LL
!LL  Model            Modification history from model version 3.3:
!LL version  Date
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL
!LL  System component:
!LL
!LL  System task:
!LL
!LL  Documentation:
!LL
!LLEND---------------------------------------------------------
!*L Arguments:-------------------------------------------------
!LL  SUBROUTINE STASH_GRIB---------------------------------------------
!LL
!LL  Purpose:
!LL   STASH_GRIB is a subroutine to code the stash section number and
!LL   parameter value in elements of the grib header.
!LL
!LL   octet 4 of the grib product definition section is the version
!LL   number of the table 2 (parameter code description)
!LL   values from 128 to 254 areavailable for local use, and we
!LL   use them to describe the stash section number of the field. for
!LL   each stash section number there are two octet 4 values. the first
!LL   is for stash parameter values from 0 to 255, the second for values
!LL   256 to 511.
!LL   octet 9 is the code value in table 2, ie stash parameter value, or
!LL   stash parameter value -256 if it is more than 255.
!LL
!LL
!LL  Written by G.Ross/ P.Smith
!LL
!LL  Model            Modification history from model version 3.3:
!LL version  Date
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL
!LL  System component:
!LL
!LL  System task:
!LL
!LL  Documentation:
!LL
!LLEND---------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE STASH_GRIB(STASH_SECTION_NUMBER,STASH_ITEM_NUMBER,     &
     &                      GRIB_BLOCK1_OCTET4,GRIB_BLOCK1_OCTET9,      &
     &                      ERROR)
      INTEGER                                                           &
     &   STASH_SECTION_NUMBER                                           &
                              ! STASH SECTION NUMBER         INPUT
     &  ,STASH_ITEM_NUMBER                                              &
                              ! STASH PARAMETER VALUE        INPUT
     &  ,GRIB_BLOCK1_OCTET4                                             &
                              ! OCTET 4 FROM GRIB PDB        OUTPUT
     &  ,GRIB_BLOCK1_OCTET9                                             &
                              ! OCTET 9 FROM GRIB PDB        OUTPUT
     &  ,ERROR                ! ERROR OUTPUT CODE            OUTPUT
!     LOCAL VARIABLES
      INTEGER                                                           &
     &   CARRY   ! CARRY VALUE FROM ODD VALUES OF GRIB_BLOCK1_OCTET4
!****
      IF(STASH_ITEM_NUMBER >  511.OR.STASH_ITEM_NUMBER <  0) THEN
        ERROR = 999
        RETURN
      ELSE IF(STASH_ITEM_NUMBER >  255) THEN
        CARRY = 1
        GRIB_BLOCK1_OCTET9 = STASH_ITEM_NUMBER - 256
      ELSE
        CARRY = 0
        GRIB_BLOCK1_OCTET9 = STASH_ITEM_NUMBER
      ENDIF
      IF((STASH_SECTION_NUMBER >= 0).AND.                               &
     &   (STASH_SECTION_NUMBER <= 62)) THEN
        GRIB_BLOCK1_OCTET4 = STASH_SECTION_NUMBER*2 + 128 + CARRY
      ELSE
        ERROR = 999
      ENDIF
      RETURN
!****
      END SUBROUTINE STASH_GRIB
