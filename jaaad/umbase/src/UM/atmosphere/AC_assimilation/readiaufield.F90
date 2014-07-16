#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Main routine for IAU scheme.


!+ Read in and check a single IAU increment field.

      SUBROUTINE ReadIAUField (                                         &
#include "argppx.h"
#include "argcona.h"
#include "arglndm.h"
     &                           L_FirstUpdate,                         &
                                                ! in
     &                           IAU_lookup,                            &
                                                ! in
     &                           FieldNum,                              &
                                                ! in
     &                           LocFldLen,                             &
                                                ! in
     &                           Field )        ! out

! Description:
!
!   1. Check field dimensions for compatibility with the model.
!   2. Read in field from IAU increment file.
!   3. If required, reset polar rows to their mean values.
!   4. If required, write basic field stats to standard output.
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   17/10/01   Original code. Adam Clayton
!   6.2   21/10/05   Replace GSYNC with SSYNC. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

#include "cppxref.h"
#include "csubmodl.h"
#include "ppxlook.h"
#include "parvars.h"
#include "typsize.h"
#include "cmaxsize.h"
#include "typcona.h"
#include "typlndm.h"
#include "clookadd.h"
#include "ctfilt.h"
#include "cntlatm.h"

! Subroutine arguments:

      LOGICAL, INTENT(IN)  ::                                           &
      
     &  L_FirstUpdate         ! .TRUE. if called during 1st IAU update

      INTEGER, INTENT(IN)  ::                                           &
                              ! Lookup tables of IAU increment file
      
     &  IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      INTEGER, INTENT(IN)  ::                                           &
      
     &  FieldNum,                                                       &
                              ! Field number in IAU increment file
     &  LocFldLen             ! Local length of field

      REAL,    INTENT(OUT) ::                                           &
      
     &  Field(LocFldLen)      ! Local part of field

! Local variables:

      INTEGER :: i,                                                     &
     &           isec,                                                  &
     &           item,                                                  &
     &           Code,                                                  &
     &           Level,                                                 &
     &           BufLen,                                                &
     &           IAU_rows,                                              &
     &           IAU_cols,                                              &
     &           Model_rows,                                            &
     &           Model_cols,                                            &
     &           ICode,                                                 &
     &           IM_index,                                              &
     &           ppx_row,                                               &
     &           halo_type_orig,                                        &
     &           s_addr,                                                &
     &           e_addr

      REAL    :: polar_sum(1),                                          &
     &           polar_row(row_length, 1),                              &
     &           polar_mean,                                            &
     &           Global_max,                                            &
     &           Global_min,                                            &
     &           Global_mean,                                           &
     &           Global_RMS

      CHARACTER(80) :: CMessage
      CHARACTER(7)  :: FldDesc

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='ReadIAUField' )

      LOGICAL, SAVE :: L_FirstCall = .TRUE.

! External subroutines called:

      EXTERNAL EReport, READFLDS, GCG_RVECSUMR, FieldStats

!- End of header ------------------------------------------------------

      Code   = IAU_lookup(ITEM_CODE, FieldNum)
      Level  = IAU_lookup(LBLEV,     FieldNum)
      BufLen = IAU_lookup(LBLREC,    FieldNum)

!----------------------------------------------------------------------
! [1]: Check field dimensions for compatibility with the model.
!----------------------------------------------------------------------

      IF (L_FirstUpdate) THEN

        IAU_rows = IAU_lookup(LBROW, FieldNum)
        IAU_cols = IAU_lookup(LBNPT, FieldNum)

        Model_rows = global_rows
        Model_cols = global_row_length

        IF (Code == 3) Model_rows = global_rows-1 ! v

        IF (IAU_rows /= Model_rows .OR.                                 &
     &      IAU_cols /= Model_cols) THEN

          ICode    = 1
          CMessage = 'Field dimension mis-match.'

          WRITE (6,*) ''
          WRITE (6,*) 'ReadIAUField: ', CMessage
          WRITE (6,*) ''
          WRITE (6,*) '  STASH code:    ', Code
          WRITE (6,*) '  Field level:   ', Level
          WRITE (6,*) '  IAU field num: ', FieldNum
          WRITE (6,*) '  IAU   rows:    ', IAU_rows
          WRITE (6,*) '  IAU   cols:    ', IAU_cols
          WRITE (6,*) '  Model rows:    ', Model_rows
          WRITE (6,*) '  Model cols:    ', Model_cols
          WRITE (6,*) ''

! DEPENDS ON: ereport
         CALL EReport (RoutineName, ICode, CMessage)

        END IF

      END IF

!----------------------------------------------------------------------
! [2]: Read in field from IAU increment file.
!----------------------------------------------------------------------

      IM_index  = INTERNAL_MODEL_INDEX(A_IM)

      ! Temporarily change PPXI halo entry so that haloes are not read:
      IF (Code  <=  1000) THEN
        ppx_row        = PPXPTR(IM_index, 0, Code)
      ELSE
        isec = Code/1000
        item = Code-1000*isec
        ppx_row        = PPXPTR(IM_index, isec, item)
      END IF

      halo_type_orig = PPXI(ppx_row, ppx_halo_type)
      PPXI(ppx_row, ppx_halo_type) = halo_type_no_halo


! DEPENDS ON: readflds
      CALL READFLDS (IAU_unit,                                          &
                                     ! in
     &               1,                                                 &
                                     ! in
     &               FieldNum,                                          &
                                     ! in
     &               IAU_lookup,                                        &
                                     ! in
     &               IAU_Len1Lookup,                                    &
                                     ! in
     &               Field,                                             &
                                     ! out
     &               BufLen,                                            &
                                     ! in
     &               IAU_FixHd,                                         &
                                     ! in
#include "argppx.h"
     &               ICode,                                             &
                                     ! out
     &               CMessage)       ! out

      IF (ICode > 0) THEN
        WRITE (6,*) 'ReadIAUField: Error reading IAU field no. ',       &
     &               FieldNum
! DEPENDS ON: ereport
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      ! Restore PPXI halo entry:
      PPXI(ppx_row, ppx_halo_type) = halo_type_orig

!----------------------------------------------------------------------
! [3]: If required, reset polar rows to their mean values.
!----------------------------------------------------------------------

      IF (L_IAU_ResetPoles          .AND.                               &
     &    Model_domain == mt_global .AND.                               &
     &    PPXI(ppx_row, ppx_grid_type) == ppx_atm_tall) THEN

        CALL GC_SSYNC ( nproc, ICode )
        IF (ICode > 0) THEN
          CMessage = 'Error from GC_SSYNC'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        IF (at_extremity(PNorth)) THEN

          s_addr = 1 + row_length * (rows-1)
          e_addr = s_addr + row_length - 1

          polar_row(:,1) = Field(s_addr:e_addr)

#if defined(REPROD)
          CALL GCG_RVECSUMR (row_length,                                &
                                                ! in
     &                       row_length,                                &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       polar_row,                                 &
                                                ! in
     &                       gc_proc_row_group,                         &
                                                ! in
     &                       ICode,                                     &
                                                ! out
     &                       polar_sum)         ! out

          IF (ICode > 0) THEN
            CMessage = 'GCG_RVECSUMR error for North Pole'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF
#else
          CALL GCG_RVECSUMF (row_length,                                &
                                                ! in
     &                       row_length,                                &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       polar_row,                                 &
                                                ! in
     &                       gc_proc_row_group,                         &
                                                ! in
     &                       ICode,                                     &
                                                ! out
     &                       polar_sum)         ! out

          IF (ICode > 0) THEN
            CMessage = 'GCG_RVECSUMF error for North Pole'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF
#endif

          polar_mean = polar_sum(1) / REAL(global_row_length)

          Field(s_addr:e_addr) = polar_mean

        END IF ! (at_extremity(PNorth))

        IF (at_extremity(PSouth)) THEN

          s_addr = 1
          e_addr = s_addr + row_length - 1

          polar_row(:,1) = Field(s_addr:e_addr)

#if defined(REPROD)
          CALL GCG_RVECSUMR (row_length,                                &
                                                ! in
     &                       row_length,                                &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       polar_row,                                 &
                                                ! in
     &                       gc_proc_row_group,                         &
                                                ! in
     &                       ICode,                                     &
                                                ! out
     &                       polar_sum)         ! out

          IF (ICode > 0) THEN
            CMessage = 'GCG_RVECSUMR error for South Pole'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF
#else
          CALL GCG_RVECSUMF (row_length,                                &
                                                ! in
     &                       row_length,                                &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       1,                                         &
                                                ! in
     &                       polar_row,                                 &
                                                ! in
     &                       gc_proc_row_group,                         &
                                                ! in
     &                       ICode,                                     &
                                                ! out
     &                       polar_sum)         ! out

          IF (ICode > 0) THEN
            CMessage = 'GCG_RVECSUMF error for South Pole'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)
          END IF
#endif

          polar_mean = polar_sum(1) / REAL(global_row_length)

          Field(s_addr:e_addr) = polar_mean

        END IF ! (at_extremity(PSouth))

      END IF ! Reset poles?

!----------------------------------------------------------------------
! [4]: If required, write basic field stats to standard output.
!----------------------------------------------------------------------

      IF (L_FirstUpdate .AND. L_IAU_Diags) THEN

        IF (L_FirstCall) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'ReadIAUField: Summary of IAU increment fields:'
          WRITE (6,*) ''
          WRITE (6,*) '  Field   Level  Max          Min         '      &
     &                              //' Mean         RMS'
          WRITE (6,*) '  -----   -----  ---          ---         '      &
     &                              //' ----         ---'
        END IF

! DEPENDS ON: fieldstats
        CALL FieldStats (                                               &
#include "argcona.h"
#include "arglndm.h"
     &                    LocFldLen,                                    &
                                                        ! in
     &                    Field,                                        &
                                                        ! in
     &                    PPXI(ppx_row, ppx_grid_type),                 &
                                                        ! in
     &                    halo_type_no_halo,                            &
                                                        ! in
     &                    Global_max,                                   &
                                                        ! out
     &                    Global_min,                                   &
                                                        ! out
     &                    Global_mean,                                  &
                                                        ! out
     &                    Global_RMS )                  ! out

        ! Get field description:
        DO i = 1, IAU_NumFldCodes
          IF (IAU_FldCodes(i) == Code)                                  &
     &      FldDesc = IAU_FldDescs(i)
        END DO

        WRITE (6,'(3A,I4,A,4(A,E12.5))')                                &
     &    '   ', FldDesc, ' ', Level, ' ',                              &
     &    ' ', Global_max,  ' ', Global_min,                            &
     &    ' ', Global_mean, ' ', Global_RMS

      END IF

      L_FirstCall = .FALSE.


      RETURN
      END SUBROUTINE ReadIAUField

!+ Calculate updated cloud variables with modified Smith (1990) scheme.


#endif
