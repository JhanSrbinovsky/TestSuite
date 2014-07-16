#if defined(C80_1A) || defined(UTILIO) || defined(FLDOP)               \
 || defined(VAROPSVER)
#if !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Write out fields to a UM format file.

      SUBROUTINE WRITFLDS ( Unit,                                       &
                                          ! in
     &                      NumFields,                                  &
                                          ! in
     &                      Position,                                   &
                                          ! in
     &                      Lookup,                                     &
                                          ! in
     &                      Len1Lookup,                                 &
                                          ! in
     &                      D1,                                         &
                                          ! in
     &                      BufLen,                                     &
                                          ! in
     &                      FixHd,                                      &
                                          ! in
#include "argppx.h"
     &                      ICode,                                      &
                                          ! out
     &                      CMessage )    ! out

! Description:
!
!   Buffers out NumFields fields from D1 to a UM format file on unit
!   Unit, starting at field number Position.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   3.1   19/02/93   Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!   3.3   25/11/93   Use PR_LFLD to print logical fields. Skip DIAG81
!                    diagnostics for observation files. Skip field
!                    summaries for boundary data. D. Robinson.
!   3.3   08/12/93   Extra argument - first dimension of lookup table.
!                    Remove hard-wired value of 64. D. Robinson
!   4.1   11/05/96   Allowed for Var and OPS files. Author Colin Parrett
!   4.1   03/01/96   Relace Char*100 with Char*80 (N Farnon)
!   4.1   18/06/96   Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!   4.4   25/04/97   Changes to write well-formed records if the
!                    input dumpfile is in that format (almost PP file
!                    format)
!                    Author: Bob Carruthers, Cray Research
!   4.5   08/07/98   Corrected error, when writing last
!                    field could cause data to be written from past
!                    the end of the input array.        Paul Burton
!   4.5   28/10/98   Introduce Single Column Model. J-C Thil.
!   5.1   13/04/00   Rewritten. New dynamics MPP capability added.
!                    Adam Clayton
!   5.3   22/11/01   Enable MPP as the only option for
!                    small executables         E.Leung
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "c_mdi.h"
#include "cntl_io.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"

#include "parvars.h"
#include "d1_addr.h"

! Subroutine arguments:

      INTEGER,       INTENT(IN) ::                                      &
     &  Len1Lookup            ! Lookup dimension.

      INTEGER,       INTENT(IN) ::                                      &
      
     &  Unit,                                                           &
                              ! Unit number of UM file.
     &  NumFields,                                                      &
                              ! Number of fields to write out.
     &  Position,                                                       &
                              ! Field number from which to begin I/O.
     &  Lookup(Len1Lookup,*)  ! Lookup table for UM file.

      REAL,          INTENT(IN) ::                                      &
      
     &  D1(*)                 ! Array containing fields to be written.

      INTEGER,       INTENT(IN) ::                                      &
      
     &  BufLen,                                                         &
                              ! Length of I/O buffer.
     &  FixHd(*)              ! Fixed-length header for UM file.

      INTEGER,       INTENT(OUT) ::                                     &
      
     &  ICode                 ! Return code. >0 => error.

      CHARACTER(80), INTENT(OUT) ::                                     &
      
     &  CMessage              ! Error message.

! Local variables:

      INTEGER :: i, j, k,                                               &
     &           LenIO,                                                 &
     &           ipts,                                                  &
                                 ! No. of values to be written to disk.
     &           WordAddress,                                           &
                                 ! Address from which to begin I/O.
     &           l_ipts,                                                &
                                 ! Record length during index search.
     &           um_sector_ipts,                                        &
                                 ! No. of words to write, rounded up.
     &           ipts_write      ! No. of words written to disk.

      REAL :: Buf( ( (BufLen+um_sector_size-1)/um_sector_size ) *       &
     &             um_sector_size )

!dir$ cache_align Buf

      INTEGER :: Address,                                               &
                                 ! Start address of D1 field.
     &           LocalLen,                                              &
     &           item,                                                  &
     &           section,                                               &
     &           model,                                                 &
     &           JCode,                                                 &
                                 ! Return code.
     &           fld_type,                                              &
     &           grid_type,                                             &
     &           halo_type,                                             &
     &           fake_D1_addr(d1_list_len) ! Fake D1_addr record to
                                           ! be fed to write_multi.

! External functions called:

      INTEGER  EXPPXI, GET_FLD_TYPE

      EXTERNAL EXPPXI, GET_FLD_TYPE

! External subroutines called:

      EXTERNAL SETPOS,  write_multi

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Initialize.
!----------------------------------------------------------------------

      Address = 1

      ICode    = 0
      J        = 0
      CMessage = ' '

!----------------------------------------------------------------------
! [2]: Write fields.
!----------------------------------------------------------------------

      DO k = Position, Position + NumFields - 1

       ! See whether data is stored in 32-bit format on disk.
        IF (MOD(Lookup(LBPACK, k), 10) == 2) THEN
          ipts = ( Lookup(LBLREC, k) + 1 ) / 2
        ELSE
          ipts = ( Lookup(LBLREC, k) + 1 )
        END IF

        IF ( Lookup(LBNREC, k) == 0    .OR.                             &
                                            ! Old format dumps.
     &       Lookup(LBNREC, k) == IMDI .OR.                             &
                                            ! \ Ocean ACOBS
     &       Lookup(LBEGIN, k) == IMDI .OR.                             &
                                            ! / files?
     &       ( Lookup(LBNREC, k) == IMDI                                &
                                            ! \ Prognostic lookup
     &         .AND. FixHd(12) < 301 )                                  &
                                            ! / in pre-vn3.2 dump.
     &     ) THEN

          WordAddress = 1
          DO i = 2, k
            IF (MOD(Lookup(LBPACK, i-1), 10) == 2) THEN
              l_ipts = ( Lookup(LBLREC, i-1) + 1 ) / 2
            ELSE
              l_ipts =   Lookup(LBLREC, i-1)
            END IF
            WordAddress = WordAddress + l_ipts
          END DO
          WordAddress = WordAddress + FixHd(160) - 2
          um_sector_ipts = ipts

        ELSE ! PP type files and post vn4.3 dumps.

          WordAddress = Lookup(LBEGIN, k)

         ! Use the stored rounded-up value:
          um_sector_ipts = Lookup(LBNREC, k)

        END IF

        ipts_write = um_sector_ipts

       ! Position file pointer:
! DEPENDS ON: setpos
        CALL SETPOS (Unit, WordAddress, ICode)

!----------------------------------------------------------------------
! [2.2]: MPP write.
!----------------------------------------------------------------------

! Set up fake D1_addr record:

        DO i = 1, d1_list_len
          fake_D1_addr(i) = 0
        END DO

        item    = MOD(Lookup(ITEM_CODE, k), 1000)
        section = ( Lookup(ITEM_CODE, k) - item ) / 1000
        model   = Lookup(MODEL_CODE, k)

! DEPENDS ON: exppxi
        halo_type = EXPPXI ( model, section, item, ppx_halo_type,       &
#include "argppx.h"
     &                       JCode, CMessage )

! DEPENDS ON: exppxi
        grid_type = EXPPXI ( model, section, item, ppx_grid_type,       &
#include "argppx.h"
     &                       JCode, CMessage )

#if defined(FLDOP) || defined(MERGE) || defined(PPTOANC)
        grid_type=1
        halo_type=3
#endif
        IF (JCode /= 0) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'WRITFLDS: Failed to get PPXREF info.'
          WRITE (6,*) ''
          WRITE (6,*) '  Model ID:      ', model
          WRITE (6,*) '  Section:       ', section
          WRITE (6,*) '  Item:          ', item
          WRITE (6,*) '  Error code:    ', JCode
          WRITE (6,*) '  Error message: ', CMessage
          WRITE (6,*) ''
          ICode    = 1
          CMessage = 'WRITFLDS: Failed to get PPXREF info.'
          RETURN
        END IF

        fake_D1_addr(d1_object_type) = prognostic
        fake_D1_addr(d1_imodl)       = model
        fake_D1_addr(d1_section)     = section
        fake_D1_addr(d1_item)        = item
        fake_D1_addr(d1_halo_type)   = halo_type

       ! Grid type: for LBCs we need some special logic...
        IF (Lookup(LBHEM, k) == 99) THEN

          IF ( Lookup(MODEL_CODE, k) == ATMOS_IM ) THEN
            fake_D1_addr(d1_grid_type) = ppx_atm_rim
          ELSE IF ( Lookup(MODEL_CODE, k) == OCEAN_IM ) THEN
            fake_D1_addr(d1_grid_type) = ppx_ocn_rim
          ELSE IF ( Lookup(MODEL_CODE, k) == WAVE_IM ) THEN
            fake_D1_addr(d1_grid_type) = ppx_wam_rim
          ELSE
            ICode = 2
            WRITE (6,*) ''
            WRITE (6,*) 'WRITFLDS: Cannot process LBC for model type',  &
     &                             Lookup(MODEL_CODE, k)
            WRITE (6,*) ''
            CMessage = 'Cannot write LBCs for this model type.'
            RETURN
          END IF

        ELSE ! Not an LBC.

          fake_D1_addr(d1_grid_type) = grid_type

        END IF

! DEPENDS ON: get_fld_type
        fld_type = GET_FLD_TYPE(grid_type)

        fake_D1_addr(d1_length)    = lasize(1, fld_type, halo_type) *   &
     &                               lasize(2, fld_type, halo_type)
        fake_D1_addr(d1_no_levels) = 1

! Write field:

! DEPENDS ON: write_multi
        CALL write_multi (                                              &
     &    Unit,         D1(Address), um_sector_ipts,                    &
                                                       ! in
     &    LenIO,        LocalLen,                                       &
                                                       ! out
     &    Lookup(1,k),  FixHd(12),   Buf,                               &
                                                       ! in
     &    fake_D1_addr,                                                 &
                                                       ! in
     &    JCode,        CMessage )                     ! out

        Address = Address + LocalLen
        if(locallen == 0)then
          address=address+LOOKUP(LBLREC,K)
        endif

! Check for errors:

        IF (LenIO /= um_sector_ipts) THEN

          WRITE (6,*) ''
          WRITE (6,*) 'WRITFLDS: Error writing field number ', k,       &
     &                         ' on unit ', Unit
          WRITE (6,*) ''
          WRITE (6,*) '  write_multi error code:    ', JCode
          WRITE (6,*) '  write_multi error message: ', CMessage
          WRITE (6,*) ''
          ICode    = JCode + 1
          CMessage = 'WRITFLDS: Error from write_multi.'
          RETURN

        END IF

      END DO ! k


      RETURN
      END SUBROUTINE WRITFLDS
#endif
#endif
