#if defined(CONVIEEE) || defined(CONVPP) || defined(CUMF)              \
 || defined(PUMF)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine FIND_MAX_FIELD_SIZE ---------------------------------
!  
!    Purpose:  Reads in and searches LOOKUP header for maximum field
!              size
!  
!    Written by A. Dickinson
!  
!  
!    Logical component number: E5
!  
!    External Documentation: None
!  
!     
!    Arguments:--------------------------------------------------------

      SUBROUTINE FIND_MAX_FIELD_SIZE(                                   &
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)
     &  NFTIN,LEN1_LOOKUP,LEN2_LOOKUP,FIXHD,MAX_FIELD_SIZE,             &
     &  wgdos_expand                                                    &
#else
     &  NFTIN,LEN1_LOOKUP,LEN2_LOOKUP,FIXHD,MAX_FIELD_SIZE              &
#endif
     &  )

      IMPLICIT NONE

      INTEGER                                                           &
     & NFTIN                                                            &
                      !IN Unit number of file
     &,LEN1_LOOKUP                                                      &
                      !IN 1st dim of LOOKUP array
     &,LEN2_LOOKUP                                                      &
                      !IN 2nd dim of LOOKUP array
     &,FIXHD(*)                                                         &
                      !IN Fixed length header
     &,MAX_FIELD_SIZE !OUT Maximum size of field held on file
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)
      integer wgdos_expand
#endif
! -------------------------------------------------------------
! Workspace usage:---------------------------------------------

      INTEGER                                                           &
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Lookup header

! -------------------------------------------------------------
! Local variables:---------------------------------------------

      INTEGER                                                           &
     & LEN_IO                                                           &
                 ! No of words transferred by BUFFIN
     &,K                                                                &
                 ! Loop index

     &,ICODE     !Return code from setpos
      REAL                                                              &
     & A         ! BUFFIN error code
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)
      integer pack_type

#include "clookadd.h"
#endif
#include "cntl_io.h"

! -------------------------------------------------------------
!   External subroutines called:-------------------------------
      EXTERNAL SETPOS,BUFFIN,IOERROR,ABORT
!*-------------------------------------------------------------

!  Internal structure: none

! Move to start of Look Up Table
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN,FIXHD(150)-1,ICODE)

! Read in fields from LOOKUP table
! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,LOOKUP(1,1),FIXHD(151)*FIXHD(152),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD(151)*FIXHD(152))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of lookup table',A,LEN_IO,              &
     &               FIXHD(151)*FIXHD(152))
! DEPENDS ON: ereport
        CALL EREPORT('FIND_MAX_FIELD_SIZE', 1000,                       &
     &   'buffer in of lookup table wrong size')

      ENDIF

! Find maximum field size
      MAX_FIELD_SIZE=0
      DO K=1,LEN2_LOOKUP
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)
        pack_type=mod(lookup(lbpack, k),10)
        if(pack_type == 1 .or. pack_type == 3                           &
     &                    .or. pack_type == 4) then
          max_field_size=max(max_field_size,                            &
     &     lookup(lbrow, k)*lookup(lbnpt, k), lookup(lblrec,k))
        else if(pack_type == 2) then
          max_field_size=max(max_field_size,                            &
     &     ((lookup(lblrec,k)+1)/2)*2)
        else
          max_field_size=max(max_field_size,lookup(lblrec,k))
        endif
#else
        MAX_FIELD_SIZE=MAX0(MAX_FIELD_SIZE,LOOKUP(15,K))
#endif
      ENDDO
#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)

      max_field_size=((max_field_size+um_sector_size-1)/                &
     & um_sector_size)*um_sector_size
      write(6,'(//''Maximum Field Size = '',i7/)') max_field_size
#endif

      RETURN
      END SUBROUTINE FIND_MAX_FIELD_SIZE
#endif
