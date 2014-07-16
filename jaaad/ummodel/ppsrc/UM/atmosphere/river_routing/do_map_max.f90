
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


!*********************************************************************
      SUBROUTINE DO_MAP_MAX(GAPS_LAMBDA_SRCE,GAPS_PHI_SRCE,LROW_SRCE    &
     &,INVERT_SRCE,DATA_SRCE,GAPS_LAMBDA_TARG,GAPS_PHI_TARG,COUNT_TARG  &
     &,BASE_TARG,LROW_TARG,WANT,MASK_TARG,INDEX_SRCE,WEIGHT,ADJUST      &
!   Add new variable to do_map_max arguments for inland
!   basin calculations
     &,DATA_TARG,ADJUST_TARG,ICODE,CMESSAGE,                            &
     & TARGET_COUNT)

!   Subroutine DO_AREAVER -------------------------------------------
!
! Purpose:
!
!   Uses the count and weights obtained by prearav but puts the
!   variable into the box mostly mapped onto.
!   Perform area-averaging to transform data from the source grid to
!   the target grid, or adjust the values on the source grid to have
!   the area-averages supplied on the target grid. The latter mode
!   is intended for adjusting values obtained by interpolating from
!   "target" to "source" in order to conserve the area-averages.
!   This mode should be used ONLY if each source box belongs in
!   exactly one target box. ADJUST=0 selects normal area-averaging,
!   ADJUST=1 selects adjustment by addition (use this mode for fields
!   which may have either sign), ADJUST=2 selects adjustment by
!   multiplication (for fields which are positive-definite or
!   negative-definite).
!
!   The shape of the source and target grids are specified by their
!   dimensions GAPS_aa_bb, which give the number of gaps in the
!   aa=LAMBDA,PHI coordinate in the bb=SRCE,TARG grid. (The product
!   of GAPS_LAMBDA_bb and GAPS_PHI_bb is the number of boxes in the
!   bb grid.)
!
!   The input and output data are supplied as 2D arrays DATA_SRCE and
!   DATA_TARG, whose first dimensions should also be supplied. Speci-
!   fying these sizes separately from the actual dimensions of the
!   grids allows for columns and rows in the arrays to be ignored.
!   A target land/sea mask should be supplied in MASK_TARG, with the
!   value indicating wanted points specified in WANT. Points which
!   are unwanted or which lie outside the source grid are not altered
!   in DATA_TARG. DATA_SRCE can optionally be supplied with its rows
!   in reverse order (i.e. with the first row corresponding to
!   minimum LAMBDA).
!
!   The arrays COUNT_TARG, BASE_TARG, INDEX_SRCE and WEIGHT should be
!   supplied as returned by PRE_AREAVER q.v.
!
!   Programming Standard, paper 4 version 4 (14.12.90)
!
! Modification history:
!   5.2    05.09.01   ADJUST_TARG calculated for use with CSRV_TWO_WAY
!                     according to ADJUST value.   R.Thorpe
!   5.2    05.09.01   And initialised for use with regridding runoff.
!                                                  C.Bunton
!
!
! Logical components covered :
!
! Project task :
!
! External documentation: Unified Model documentation paper No:
!                         Version:
!
!END -----------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     & GAPS_LAMBDA_SRCE                                                 &
                               !IN number lambda gaps in source grid
     &,GAPS_PHI_SRCE                                                    &
                               !IN number phi gaps in source grid
     &,LROW_SRCE                                                        &
                               !IN first dimension of source arrays
     &,GAPS_LAMBDA_TARG                                                 &
                               !IN number lambda gaps in target grid
     &,GAPS_PHI_TARG                                                    &
                               !IN number phi gaps in target grid
     &,LROW_TARG                                                        &
                               !IN first dimension of target arrays
     &,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                       &
!                              !IN no. of source boxes in target box
!   Define  new variable in do_map_max arguments for inland basin
!   calculations
     &,TARGET_COUNT(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                     &
!                    COUNT_TARG FOR OUTPUT  TO RIVER1A

     &,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                        &
!                              !IN first index in list for target box
     &,INDEX_SRCE(*)                                                    &
                               !IN list of source box indices
     &,ADJUST                                                           &
                               !IN selects normal or adjust mode
     &,ICODE                   !OUT return code
      LOGICAL                                                           &
     & INVERT_SRCE                                                      &
                               !IN DATA_SRCE rows in reverse order
     &,WANT                                                             &
                               !IN indicator of wanted points in mask
     &,MASK_TARG(LROW_TARG,*)  !IN land/sea mask for target grid

!     NB alternative intents below apply for normal/adjust mode
      REAL                                                              &
     & DATA_SRCE(LROW_SRCE,*)                                           &
                               !IN/INOUT data on source grid
     &,WEIGHT(*)                                                        &
                               !IN list of weights for source boxes
     &,DATA_TARG(LROW_TARG,*)                                           &
                               !INOUT/IN data on target grid
     &,ADJUST_TARG(LROW_TARG,*)

      CHARACTER                                                         &
     & CMESSAGE*(*)            !OUT error message
!
      INTEGER                                                           &
     & IP                                                               &
                               ! pointer into lists
     &,I,J,K                                                            &
                                   ! loop index
     &,IX1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)                              &
!                              ! working SRCE LAMBDA indices
     &,IY1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)                              &
!                              ! working SRCE PHI indices
     &,IX2,IY2                 ! working TARG LAMBDA/PHI indices
      REAL                                                              &
     & TEMP_TARG                                                        &
                               ! workspace for area-average
     &,DELTA                                                            &
                               ! additive adjustment
     &,RATIO                                                            &
                               ! multiplicative adjustment
     &,WEIGHT_MAX              ! max weighting
!
!    Loop over all target boxes and calculate values as required.
!
!
!     The weights and source box indices are recorded in continuous
!     lists. COUNT_TARG indicates how many consecutive entries in these
!     lists apply to each target box.
!
      DO IY2=1,GAPS_PHI_TARG
        DO IX2=1,GAPS_LAMBDA_TARG
!   Calculate new counter for inland basin outflow
         TARGET_COUNT(IX2,IY2)=COUNT_TARG(IX2,IY2)

          IF (MASK_TARG(IX2,IY2).EQV.WANT) THEN
          IF (COUNT_TARG(IX2,IY2) /= 0) THEN
            IP=BASE_TARG(IX2,IY2)+1
            WEIGHT_MAX=WEIGHT(IP)
            IX1(1)=MOD(INDEX_SRCE(IP)-1,GAPS_LAMBDA_SRCE)+1
            IY1(1)=(INDEX_SRCE(IP)-1)/GAPS_LAMBDA_SRCE+1
            IF (INVERT_SRCE) IY1(1)=GAPS_PHI_SRCE-IY1(1)+1

            DO I=1,COUNT_TARG(IX2,IY2)
              IP=BASE_TARG(IX2,IY2)+I
              IX1(I)=MOD(INDEX_SRCE(IP)-1,GAPS_LAMBDA_SRCE)+1
              IY1(I)=(INDEX_SRCE(IP)-1)/GAPS_LAMBDA_SRCE+1
              IF (INVERT_SRCE) IY1(I)=GAPS_PHI_SRCE-IY1(I)+1
              IF(WEIGHT(IP) >  WEIGHT_MAX)THEN
                WEIGHT_MAX = WEIGHT(IP)
                IX1(1)=IX1(I)
                IY1(1)=IY1(I)
              ENDIF
            ENDDO

            DATA_SRCE(IX1(1),IY1(1)) = DATA_SRCE(IX1(1),IY1(1)) +       &
     &                                     DATA_TARG(IX2,IY2)

          ENDIF
         ENDIF
        ENDDO
      ENDDO

!
      ICODE=0
      CMESSAGE=' '

      RETURN
      END SUBROUTINE DO_MAP_MAX

