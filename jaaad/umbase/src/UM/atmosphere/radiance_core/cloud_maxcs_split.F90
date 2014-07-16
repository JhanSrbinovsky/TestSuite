#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to split cloud into maximally overlapped C/S.
!
! Method:
!
!   Convective cloud is left-justified in the grid-box while
!   stratiform cloud is right-justified.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CLOUD_MAXCS_SPLIT(IERR, N_PROFILE, N_LAYER, N_CLOUD_TOP&
     &  , W_CLOUD, FRAC_CLOUD                                           &
     &  , N_CLOUD_TYPE                                                  &
     &  , N_COLUMN_CLD, N_COLUMN_SLV, LIST_COLUMN_SLV                   &
     &  , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                     &
     &  , ND_PROFILE, ND_LAYER, ID_CT, ND_COLUMN, ND_CLOUD_TYPE)
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_COLUMN                                                     &
!           Size allocated for columns at a point
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for columns at a point
     &  , ID_CT
!           Topmost allocated cloudy layer
!
!     Include header files.
#include "c_kinds.h"
#include "error_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "cloud_type_pcf3z.h"
!
!
!     Dummy variables
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_TYPE
!           Number of types of cloud
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Amount of cloud
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)
!           Fractions of different types of cloud
!
      INTEGER, INTENT(OUT) ::                                           &
     &    N_COLUMN_CLD(ND_PROFILE)                                      &
!           Number of columns in each profile (including those of
!           zero width)
     &  , N_COLUMN_SLV(ND_PROFILE)                                      &
!           Number of columns to be solved in each profile
     &  , LIST_COLUMN_SLV(ND_PROFILE, ND_COLUMN)                        &
!           List of columns requiring an actual solution
     &  , I_CLM_LYR_CHN(ND_PROFILE, ND_COLUMN)                          &
!           Layer in the current column to change
     &  , I_CLM_CLD_TYP(ND_PROFILE, ND_COLUMN)
!           Type of cloud to introduce in the changed layer
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    AREA_COLUMN(ND_PROFILE, ND_COLUMN)
!           Area of each column
!
!
!     Local variables
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , II                                                            &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , N_CLD_LAYER                                                   &
!           Number of cloudy layers
     &  , PTR_ST                                                        &
!           Pointer to stratiform cloud in arrays
     &  , PTR_CNV                                                       &
!           Pointer to convective cloud in arrays
     &  , KEY_ST(ND_LAYER+1-ID_CT)                                      &
!           Pointers to layers listing left edge of stratiform cloud
!           in increasing order
     &  , KEY_CNV(ND_LAYER+1-ID_CT)                                     &
!           Pointers to layers listing right edge of convective cloud
!           in increasing order
     &  , I_KEY_CNV                                                     &
!           Current pointer to convective cloud
     &  , I_KEY_ST
!           Current pointer to stratiform cloud
      REAL  (Real64) ::                                                 &
     &    CNV_RIGHT(ND_LAYER+1-ID_CT)                                   &
!           Right edges of convective cloud
     &  , STRAT_LEFT(ND_LAYER+1-ID_CT)                                  &
!           Left edges of stratiform cloud
     &  , X_CNV                                                         &
!           Right edge of current convective cloud
     &  , X_ST                                                          &
!           Left edge of current stratiform cloud
     &  , X_DONE                                                        &
!           Fraction of the column treated
     &  , X_NEW_DONE                                                    &
!           Fraction of the column treated after adding new column
     &  , DX_COL
!           Width of the current column
      REAL  (Real64) ::                                                 &
     &    TOL_CLOUD
!           Tolerance used in neglecting cloudy columns
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SHELL_SORT
!
!
!
!     Set the tolerance used for clouds.
      TOL_CLOUD=1.0E+04_Real64*EPSILON(TOL_CLOUD)
!
      IF (N_CLOUD_TYPE == 2) THEN
        PTR_ST=0
        PTR_CNV=0
        DO K=1, N_CLOUD_TYPE
          IF (K == IP_CLOUD_TYPE_STRAT) PTR_ST=K
          IF (K == IP_CLOUD_TYPE_CONV) PTR_CNV=K
        ENDDO
        IF ( (PTR_ST == 0).OR.(PTR_CNV == 0) ) THEN
          WRITE(IU_ERR, '(/A)')                                         &
     &      '*** Error: A type of cloud is missing.'
          IERR=I_ERR_FATAL
          RETURN
        ENDIF
      ELSE IF (N_CLOUD_TYPE == 1) THEN
!       Only stratiform cloud is present.
        PTR_CNV=0
        PTR_ST=1
      ELSE
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: There are too many types of cloud for '           &
     &    //'the type of overlap.'
        IERR=I_ERR_FATAL
        RETURN
      ENDIF
!
!
      N_CLD_LAYER=N_LAYER+1-N_CLOUD_TOP
!     We decompose a column at a time, as this is algorithmically
!     easier, even if not compatible with vectorization.
      DO L=1, N_PROFILE
!
!       Cloud is aligned with convective cloud against the left-hand
!       edge of the grid-box and stratiform cloud to the right. We
!       therefore need to find the right-hand edge of convective
!       cloud and the left-hand edge of stratiform cloud to partition.
        DO I=N_CLOUD_TOP, N_LAYER
          II=I+1-N_CLOUD_TOP
!
          IF (N_CLOUD_TYPE == 2) THEN
!           Calculate an explicit position for convective cloud if
!           included.
            CNV_RIGHT(II)=W_CLOUD(L, I)*FRAC_CLOUD(L, I, PTR_CNV)
          ELSE
!           In the absence of convective cloud set its width to 0.
            CNV_RIGHT(II)=0.0E+00_Real64
          ENDIF
!
          STRAT_LEFT(II)=1.0E+00_Real64-W_CLOUD(L, I)*FRAC_CLOUD(L, I   &
     &      , PTR_ST)
!         Initialize the sorting key.
          KEY_ST(II)=II
          KEY_CNV(II)=II
        ENDDO
!
!       Find the key ranking these edges in increasing order.
! DEPENDS ON: shell_sort
        CALL SHELL_SORT(N_CLD_LAYER, KEY_CNV, CNV_RIGHT)
! DEPENDS ON: shell_sort
        CALL SHELL_SORT(N_CLD_LAYER, KEY_ST, STRAT_LEFT)
!
!
!       Build up the list of notional columns and the list of those
!       where a solution is required.
        N_COLUMN_CLD(L)=0
        N_COLUMN_SLV(L)=0
!       Set the changes from a totally clear column to the first
!       actually used.
        DO I=N_CLOUD_TOP, N_LAYER
          II=I+1-N_CLOUD_TOP
          IF (CNV_RIGHT(II) >  TOL_CLOUD) THEN
            IF (N_COLUMN_CLD(L) <  ND_COLUMN) THEN
              N_COLUMN_CLD(L)=N_COLUMN_CLD(L)+1
            ELSE
              WRITE(IU_ERR, '(/A)')                                     &
     &          '*** Error: ND_COLUMN is too small for the cloud '      &
     &          //'geometry selected.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            I_CLM_LYR_CHN(L, N_COLUMN_CLD(L))=I
            I_CLM_CLD_TYP(L, N_COLUMN_CLD(L))=PTR_CNV
            AREA_COLUMN(L, N_COLUMN_CLD(L))=0.0E+00_Real64
          ENDIF
          IF (STRAT_LEFT(II) <= TOL_CLOUD) THEN
            IF (N_COLUMN_CLD(L) <  ND_COLUMN) THEN
              N_COLUMN_CLD(L)=N_COLUMN_CLD(L)+1
            ELSE
              WRITE(IU_ERR, '(/A)')                                     &
     &          '*** Error: ND_COLUMN is too small for the cloud '      &
     &          //'geometry selected.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            I_CLM_LYR_CHN(L, N_COLUMN_CLD(L))=I
            I_CLM_CLD_TYP(L, N_COLUMN_CLD(L))=PTR_ST
            AREA_COLUMN(L, N_COLUMN_CLD(L))=0.0E+00_Real64
          ENDIF
        ENDDO
!
!       Now set up the mapping changing the contents of each layer
!       proceeding to the right.
        X_DONE=0.0E+00_Real64
!       Set the positions of the next convective and stratiform
!       changes, together with their corresponding indices.
        I_KEY_CNV=1
        X_CNV=CNV_RIGHT(KEY_CNV(1))
        DO WHILE ( (I_KEY_CNV <  N_CLD_LAYER).AND.(X_CNV <  TOL_CLOUD) )
          I_KEY_CNV=I_KEY_CNV+1
          X_CNV=CNV_RIGHT(KEY_CNV(I_KEY_CNV))
        ENDDO
        I_KEY_ST=1
        X_ST=STRAT_LEFT(KEY_ST(1))
        DO WHILE ( (I_KEY_ST <  N_CLD_LAYER).AND.(X_ST <  TOL_CLOUD) )
          I_KEY_ST=I_KEY_ST+1
          X_ST=STRAT_LEFT(KEY_ST(I_KEY_ST))
        ENDDO
!
!       Proceed throught the grid-box making the changes.
        DO WHILE (X_DONE <  1.0E+00_Real64-TOL_CLOUD)
!
          IF (X_CNV <= X_ST) THEN
!           The next change involves clearing convective cloud.
            IF (N_COLUMN_CLD(L) <  ND_COLUMN) THEN
              N_COLUMN_CLD(L)=N_COLUMN_CLD(L)+1
            ELSE
              WRITE(IU_ERR, '(/A)')                                     &
     &          '*** Error: ND_COLUMN is too small for the cloud '      &
     &          //'geometry selected.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            I_CLM_LYR_CHN(L, N_COLUMN_CLD(L))=KEY_CNV(I_KEY_CNV)        &
     &        +N_CLOUD_TOP-1
            I_CLM_CLD_TYP(L, N_COLUMN_CLD(L))=0
            X_NEW_DONE=X_CNV
            I_KEY_CNV=I_KEY_CNV+1
            IF (I_KEY_CNV <= N_CLD_LAYER) THEN
              X_CNV=CNV_RIGHT(KEY_CNV(I_KEY_CNV))
            ELSE
!             There are no further changes to convective cloud
!             right of this.
              X_CNV=1.0E+00_Real64
            ENDIF
          ELSE IF (X_CNV >  X_ST) THEN
!           The next change involves introducing stratiform cloud.
            IF (N_COLUMN_CLD(L) <  ND_COLUMN) THEN
              N_COLUMN_CLD(L)=N_COLUMN_CLD(L)+1
            ELSE
              WRITE(IU_ERR, '(/A)')                                     &
     &          '*** Error: ND_COLUMN is too small for the cloud '      &
     &          //'geometry selected.'
              IERR=I_ERR_FATAL
              RETURN
            ENDIF
            I_CLM_LYR_CHN(L, N_COLUMN_CLD(L))=KEY_ST(I_KEY_ST)          &
     &        +N_CLOUD_TOP-1
            I_CLM_CLD_TYP(L, N_COLUMN_CLD(L))=PTR_ST
            X_NEW_DONE=X_ST
            I_KEY_ST=I_KEY_ST+1
            IF (I_KEY_ST <= N_CLD_LAYER) THEN
              X_ST=STRAT_LEFT(KEY_ST(I_KEY_ST))
            ELSE
!             There are no further changes to stratiform cloud
!             right of this.
              X_ST=1.0E+00_Real64
            ENDIF
          ENDIF
!
!         If both convective and stratiform right markers have
!         reached 1 we have a closing column.
          IF ( (X_ST >  1.0E+00_Real64-TOL_CLOUD).AND.                  &
     &         (X_CNV >  1.0E+00_Real64-TOL_CLOUD) )                    &
     &      X_NEW_DONE=1.0E+00

!
!         If this new column is wide enough we solve within it.
          DX_COL=X_NEW_DONE-X_DONE
          IF (DX_COL >  TOL_CLOUD) THEN
            N_COLUMN_SLV(L)=N_COLUMN_SLV(L)+1
            LIST_COLUMN_SLV(L, N_COLUMN_SLV(L))=N_COLUMN_CLD(L)
            AREA_COLUMN(L, N_COLUMN_CLD(L))=DX_COL
            X_DONE=X_NEW_DONE
          ELSE
            AREA_COLUMN(L, N_COLUMN_CLD(L))=0.0E+00_Real64
          ENDIF
!
        ENDDO
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CLOUD_MAXCS_SPLIT
#endif
