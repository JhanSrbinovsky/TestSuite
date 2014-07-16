

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Function returns full UM version id
!
! Function interface:
      INTEGER Function GET_UM_VERSION_ID(model_id)

      IMPLICIT NONE
!
! Description:
!   Function generates the full UM version ID
!
! Method:
!   The function uses function GET_UM_VERSION to extract UM version
!   number from environment variable $VN, and then generates the
!   UM version ID as an eight digit number xxxxyyyy, where xxxx is
!   the UM version number e.g. 0503 and yyyy is the model identifier
!   e.g. 1111 for unified model.
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version    Date     Comment
! -------    ----     -------
! 5.3        29/01/02 Original code. D.M. Goddard
! 6.0        05/09/03 Added new def for use with makebc. R. Sempers
! 6.0        17/10/03 Optimisation - only calculate once. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + some CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Function arguments
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN) ::                                            &
     &  model_id                        ! Model ID identifier

! Local parameters:
      INTEGER, PARAMETER ::                                             &
     &  multiplication_factor = 10000                                   &
     &, no_version = -9999                                              &
     &, unset_version = -5555

      CHARACTER (LEN=*), PARAMETER ::                                   &
     &  routinename ='get_um_version_id'

! Local scalars:
      INTEGER      ::   um_version
      INTEGER,SAVE ::   um_version_id=unset_version
      INTEGER      ::   get_um_version

! Function & Subroutine calls:
      External get_um_version

!- End of header

      IF (um_version_id == unset_version) THEN
! DEPENDS ON: get_um_version
        um_version = get_um_version()

        IF (um_version  /=  no_version ) THEN
          um_version_id = um_version*multiplication_factor              &
     &                        +model_id
        ELSE
          um_version_id = model_id
        END IF
      ENDIF

      get_um_version_id = um_version_id

      RETURN
      END FUNCTION GET_UM_VERSION_ID
