#if defined(C70_1A) && !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INTF_OUT -----------------------------------------------
!
!LL  Purpose: To open boundary files
!LL           Ocean has 4 on Fortran unit numbers 100-103
!LL           Atmos has 8 on Fortran unit numbers 140-147
!LL
!LL  Model            Modification history from model version 4.5
!LL version  Date
!LL  4.5   3/09/98    New deck added M.J.Bell
!LL  5.3  17/08/01    Output data in WFIO format. Also output grid code
!LL                   and number of levels for each field. M.J.Bell
!LL  5.4  31/08/02    Changes to allow for flexibility over which fields
!LL                   output in lbcs, also to allow for extra choices:
!LL                   SSH/rigid-lid pressure and barotropic currents.
!LL                   D.Storkey
!LL  5.5  17/02/03    Prevent Wave model data from being packed.
!LL                   D.Holmes-Bell
!    5.5  05/02/03    Portability changes allowing for big_endian
!                     I/O on little_endian platforms.        P.Dando
!    6.2  23/11/05    Removed all references to the wavemodel.
!                     T.Edwards
!LL
!LLEND ---------------------------------------------------------------
       subroutine intf_out (                                            &
#include "adumlen.h"
#include "ainflen.h"
#include "argdum.h"
#include "arginf.h"
#include "argppx.h"
     & NFTOUT, JINTF, im, mype, INTF_LOOKUPS_OUT,                       &
     & INTF_PACK, INTFWIDTH, LEN_INTF_P, LEN_INTF_U,                    &
     & len_intf_data, item_intf, len_bdy_flds,                          &
     & no_levs_bdy_flds, grid_code_bdy_flds, len_buffer,                &
     & dump_lookup_intf, intf_data, icode, cmessage )
!---------------------------------------

      implicit none

#include "cdumlen.h"
#include "cinflen.h"
#include "typdum.h"
#include "typinf.h"

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

      integer                                                           &
     &    NFTOUT,                                                       &
                          ! unit to write to
     &    JINTF,                                                        &
                          ! number of this boundary file
     &    im,                                                           &
                          ! internal model identifier
     &    mype,                                                         &
                          ! number of "my" processor
     &    INTF_LOOKUPS_OUT(N_INTF),                                     &
                                   ! Number of fields output
     &    INTF_PACK(N_INTF),                                            &
                                   ! Packing Indicator for data
     &    INTFWIDTH(N_INTF),                                            &
                                   ! Width of interface zone
     &    LEN_INTF_P(N_INTF),                                           &
                                   ! Length of interface p field
     &    LEN_INTF_U(N_INTF),                                           &
                                   ! Length of interface u field
     &    len_intf_data,                                                &
                           ! length of field of data to output
     &    item_intf(INTF_LOOKUPS),                                      &
                                      ! stash item numbers of fields
     &    len_bdy_flds(INTF_LOOKUPS),                                   &
                                       ! length of interface fields
     &    no_levs_bdy_flds (intf_lookups),                              &
                                           ! number of levels in fields
     &    grid_code_bdy_flds (intf_lookups),                            &
                                             ! grid codes of fields
     &    len_buffer,                                                   &
                           ! length of last WFIO field written out
     &    dump_lookup_intf(INTF_LOOKUPS) ! numbers of corresponding
                                          ! dump lookup tables

      real intf_data( len_intf_data )  ! boundary data for output
      integer icode
       character*80 cmessage

!-----------------------------------------------------------
#include "cmaxsize.h"
#include "cmaxsizo.h"
#include "chsunits.h"
! for FT_LASTFIELD
#include "ihisto.h"
! for  MODEL_FT_UNIT
#include "clfhist.h"
#include "clookadd.h"
! for ft_steps and lcal360
#include "cntlall.h"
! for STEPS_PER_SEC etc.
#include "cntlgen.h"
! for basis_time_days, STEP_im etc.
#include "ctime.h"
! stores um_sector_size
#include "cntl_io.h"

!L Local variables

      real   buffer (len_buffer) !  temporary buffer for WFIO data

      integer                                                           &
     &  LEN_PPNAME,                                                     &
                     ! length of pp file name
     &  NTIME,                                                          &
                     ! number of field in output file
     &  lookup_start,                                                   &
                     ! start location to write lookup table
     &  LEN_IO,                                                         &
     &  disk_address,                                                   &
     &  start_addr,                                                     &
     &  len_data,                                                       &
                        !
     &  var,                                                            &
                        ! loop index for variable
     &  i,                                                              &
                        ! loop index
     &  N1,                                                             &
                        ! local packing index
     &  disk_length,                                                    &
     &  iaddr,                                                          &
     &  j,                                                              &
                        ! loop index
     &  data_start,                                                     &
                        ! start location for writing interface data

     &  EXPPXI          ! function
       EXTERNAL EXPPXI

      integer SEC,YY,MM,DD,HR,MN,SS,DAY_NO

      real a_io

      LOGICAL LPACK_32B,                                                &
                            ! pack as 32 bit numbers
     &        LPACK_PPXREF  !

      CHARACTER*80 STRING         ! work array
      CHARACTER*14 PPNAME         ! boundary output filename
!-------------------------------------------------------------

!L 0. Miscellaneous Preliminaries

      LPACK_32B = INTF_PACK(JINTF) == 1
      LPACK_PPXREF = INTF_PACK(JINTF) == 2

!L 1.0   Open file; determine where to write new data

!L     Open boundary output file if reinitialised during run

      IF (FT_STEPS(NFTOUT) >  0) THEN
        STRING = MODEL_FT_UNIT(NFTOUT)
        PPNAME = STRING(18:31)
        LEN_PPNAME = LEN(PPNAME)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,PPNAME,LEN_PPNAME,1,1,ICODE)
        IF (ICODE /= 0) THEN
          CMESSAGE="INTF_OUT: Error opening preassigned boundary file"
          GO TO 999   !  Return
        ENDIF
      ENDIF

!      Determine position where to Buffer out data to
      NTIME=FT_LASTFIELD(NFTOUT)+1

!L 2.  Set up headers

!L 2.1 Fixed length header
      FIXHD_INTF(152,JINTF) = INTF_LOOKUPS_OUT(JINTF)*NTIME
      FIXHD_INTF(161,JINTF) = LEN_INTF_DATA*NTIME

!L 2.2 Integer Constants
      INTHD_INTF(3,JINTF) = NTIME

!L 2.3 LOOKUP Table

!  2.3.1   Determine position in LOOKUP table
      LOOKUP_START = FIXHD_INTF(150,JINTF)                              &
     &   + FIXHD_INTF(151,JINTF)*INTF_LOOKUPS_OUT(JINTF)*(NTIME-1) - 1

! 2.3.2  For well-formed I/O re-read the last lookup
!       table on disk to find disk_address
!       also set initial start address

      if(ntime /= 1) then
! DEPENDS ON: setpos
        call setpos(nftout, lookup_start-len1_lookup, icode)
! DEPENDS ON: buffin
        call buffin(nftout, lookup_intf(1, 1, jintf), len1_lookup,      &
     &   len_io, a_io)

!--check for errors
        if(a_io /= -1.0 .or. len_io /= len1_lookup) then
! DEPENDS ON: ioerror
          call ioerror('intf_out: Buffer in of Last Lookup Header',     &
     &     a_io, len_io, len1_lookup)
          cmessage=' intf_out: I/O Error on reading last lookup'
          icode=5
          goto 999
        endif

!--compute the new disk address from the last address and length
        disk_address=lookup_intf(lbegin, 1, jintf)+                     &
     &               lookup_intf(lbnrec, 1, jintf)

      else     ! ntime

        disk_address=fixhd_intf(160, jintf)-1
      endif     ! ntime

!--round this disk address to ensure we start on a sector boundary
      disk_address=((disk_address+um_sector_size-1)/                    &
     & um_sector_size)*um_sector_size

! - start address (not used by well formed I/O ?)
      START_ADDR = FIXHD_INTF(161,JINTF)-LEN_INTF_DATA+1

! 2.3.3  Check that there is enough space for this entry in LOOKUP table

      IF (FIXHD_INTF(150,JINTF)+                                        &
     &    FIXHD_INTF(151,JINTF)*FIXHD_INTF(152,JINTF) >                 &
     &   FIXHD_INTF(160,JINTF)) THEN
        CMESSAGE=' INTF_OUT: Insufficient space for headers in boundary &
     &                       dataset.'
        ICODE=1
        GO TO 999   !  Return
      ENDIF

! 2.3.5 Set validity times

      SEC = STEPim(im) * SECS_PER_PERIODim(im) /                        &
     &      STEPS_PER_PERIODim(im)

! DEPENDS ON: sec2time
      CALL SEC2TIME(0,SEC,BASIS_TIME_DAYS,BASIS_TIME_SECS,              &
     &                  YY,MM,DD,HR,MN,SS,DAY_NO,LCAL360)

      DO VAR = 1, INTF_LOOKUPS_OUT(JINTF)

! 2.3.6 Initialise lookup tables (with values from dump lookup tables)

        DO I=1,LEN1_LOOKUP
          LOOKUP_INTF(I,VAR,JINTF)=LOOKUP(I,dump_lookup_intf(var))
        ENDDO

! 2.3.7 Set times in lookup tables

        LOOKUP_INTF(LBYR ,VAR,JINTF) = YY
        LOOKUP_INTF(LBMON,VAR,JINTF) = MM
        LOOKUP_INTF(LBDAT,VAR,JINTF) = DD
        LOOKUP_INTF(LBHR ,VAR,JINTF) = HR
        LOOKUP_INTF(LBMIN,VAR,JINTF) = MN
        LOOKUP_INTF(LBDAY,VAR,JINTF) = DAY_NO

        LOOKUP_INTF(LBYRD ,VAR,JINTF) = FIXHD(21)
        LOOKUP_INTF(LBMOND,VAR,JINTF) = FIXHD(22)
        LOOKUP_INTF(LBDATD,VAR,JINTF) = FIXHD(23)
        LOOKUP_INTF(LBHRD ,VAR,JINTF) = FIXHD(24)
        LOOKUP_INTF(LBMIND,VAR,JINTF) = FIXHD(25)
        LOOKUP_INTF(LBDAYD,VAR,JINTF) = FIXHD(27)

!  2.3.8 Set the length of the field in LOOKUP table
! (simpler than in original atmosphere code)  !! CHECK THIS !!

        LOOKUP_INTF(LBLREC,VAR,JINTF) = len_bdy_flds(var)

! 2.3.9 Set packing info
        N1 = 0   !  Data not packed
        IF (LPACK_32B) N1 = 2  ! Data packed as 32 bits
        IF (LPACK_PPXREF) THEN
! DEPENDS ON: exppxi
          N1 = EXPPXI(im,0,item_intf(jintf),ppx_dump_packing,           &
#include "argppx.h"
     &                 icode,cmessage)
          if (icode  >   0) then
             write(6,*) 'exppxi failed in intf_out'
             go to 999
          end if
        END IF
        LOOKUP_INTF(LBPACK,VAR,JINTF)= N1

! 2.3.10 Store the disk address; and calculate for next field
        lookup_intf(lbegin, var, jintf)=disk_address

!--fetch the data field length, allowing for packing
        if(mod(lookup_intf(lbpack, var, jintf), 10) == 2) then
          disk_length=(lookup_intf(lblrec, var, jintf)+1)/2
        else
          disk_length=lookup_intf(lblrec, var, jintf)
        endif

        disk_length=((disk_length+um_sector_size-1)/um_sector_size)     &
     &                                            * um_sector_size

!--store the rounded-up length
        lookup_intf(lbnrec, var, jintf)=disk_length

!--update the disk address
        disk_address=disk_address+disk_length

! 2.3.11 Set other elements in the lookup table

! grid code (set using cppxref definitions at 5.3)
        LOOKUP_INTF(LBCODE,VAR,JINTF) = grid_code_bdy_flds(var)

        LOOKUP_INTF(LBHEM,VAR,JINTF) = no_levs_bdy_flds(var) + 100
        LOOKUP_INTF(ITEM_CODE,VAR,JINTF) = item_intf(var)
        LOOKUP_INTF(LBROW,VAR,JINTF)=INTFWIDTH(JINTF)


!L The lookup table for rigid-lid pressure and barotropic currents are
!L initialized using those for the temperature field and baroclinic
!L currents. So need to reset field and STASH codes here. Set STASH
!L codes equal to ITEM_INTF(var) for all fields to save using IF tests.

         if (item_intf(var)  ==  31285) then
            LOOKUP_INTF(LBFC,var,JINTF)=617
         endif
         if (item_intf(var)  ==  135) then
            LOOKUP_INTF(LBFC,var,JINTF)=711
         endif
         if (item_intf(var)  ==  136) then
            LOOKUP_INTF(LBFC,var,JINTF)=712
         endif

         LOOKUP_INTF(ITEM_CODE,VAR,JINTF) = ITEM_INTF(VAR)

        LOOKUP_INTF(LBNPT,VAR,JINTF) =                                  &
     &           LEN_INTF_P(JINTF)/INTFWIDTH(JINTF)
        IF (  ( IM  ==  1. .AND. (VAR == 2.OR.VAR == 3) )               &
     &  .OR. (IM == 2. .AND. LOOKUP_INTF(LBFC,VAR,JINTF) >  700)) THEN

          LOOKUP_INTF(LBNPT,VAR,JINTF) =                                &
     &           LEN_INTF_U(JINTF)/INTFWIDTH(JINTF)
        END IF

        LOOKUP_INTF(LBLEV,VAR,JINTF)=-1
        LOOKUP_INTF(NADDR,VAR,JINTF) = START_ADDR
        START_ADDR = START_ADDR + LOOKUP_INTF(LBLREC,VAR,JINTF)

      END DO  ! VAR


!L 3.  Write out headers/data

!L 3.1 Fixed length header

        IADDR = 0
! DEPENDS ON: setpos
        CALL SETPOS (NFTOUT,IADDR,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FIXHD_INTF(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of fixed length header',A_IO,LEN_IO, &
     &                  LEN_FIXHD)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=2
          GO TO 999   !  Return
        END IF

!L 3.2 Integer constants

! DEPENDS ON: buffout
        CALL BUFFOUT (NFTOUT,INTHD_INTF(1,JINTF),                       &
     &                PP_LEN_INTHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= PP_LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of integer header',A_IO,LEN_IO,      &
     &                  PP_LEN_INTHD)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=3
          GO TO 999   !  Return
        END IF

!L 3.3 PP headers in LOOKUP table
! DEPENDS ON: setpos
        CALL SETPOS(NFTOUT,LOOKUP_START,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LOOKUP_INTF(1,1,JINTF),                     &
     &               LEN1_LOOKUP*INTF_LOOKUPS_OUT(JINTF),LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.                                             &
     &     LEN_IO /= LEN1_LOOKUP*INTF_LOOKUPS_OUT(JINTF)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of PP header',A_IO,LEN_IO,           &
     &                  LEN1_LOOKUP*INTF_LOOKUPS)
          CMESSAGE=' intf_out: I/O ERROR '
          ICODE=4
          GO TO 999   !  Return
        END IF

!L 4. Write out data, packing it if necessary

        IADDR = 1

        DO var = 1,INTF_LOOKUPS_OUT(JINTF)

         DATA_START =   lookup_intf(lbegin, VAR, jintf)
         len_data =     LOOKUP_INTF(LBLREC,VAR,JINTF)
         disk_length =  lookup_intf(lbnrec, var, jintf)

         IF (mype  ==  0) THEN

          do i = 1, disk_length
            buffer(i) = 0.0
          end do

          IF (MOD(LOOKUP_INTF(LBPACK,VAR,JINTF),10) == 2) THEN

!L 4.1 Write data to a buffer

!L 4.1.1 if data is to be packed, pack it in buffer array

! DEPENDS ON: pack21
           CALL PACK21(len_data, INTF_DATA(IADDR),                      &
                       BUFFER )

!--check that we are not packing an odd nuber of words
            if((lookup_intf(lblrec,var,jintf)/2)*2  /=                  &
     &       lookup_intf(lblrec,var,jintf)) then
              if(mype == 0) then
                write(6,7734) lookup_intf(lblrec,var,jintf)
              end if
7734          format(/'LBC Data contains ',i10,' Words, which is',      &
     &         ' an Odd Number which is not allowed for 32-bit',        &
     &         ' Packing')
            endif

!L 4.1.2 else just write it to the buffer array

          ELSE        !    LOOKUP_INTF(LBPACK..

            do i=1,len_data
              buffer(i) = INTF_DATA(IADDR-1+i)
            end do

          END IF   ! LOOKUP_INTF(LBPACK

         END IF  ! mype  ==  0

!L 4.2 Write out the buffer array

! DEPENDS ON: setpos
         CALL SETPOS(NFTOUT,DATA_START,ICODE)

         IF (MOD(LOOKUP_INTF(LBPACK,VAR,JINTF),10)  ==  2) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*disk_length 32 bit words using BUFFO32
! DEPENDS ON: buffo32
           CALL BUFFO32(NFTOUT,buffer,2*disk_length,LEN_IO,A_IO)
! And then halve LEN_IO to satisfy tests against disk_length
           LEN_IO = LEN_IO/2
         ELSE
! For non-packed data
! DEPENDS ON: buffout
           CALL BUFFOUT(NFTOUT,buffer,disk_length,LEN_IO,A_IO)
         ENDIF


! Check for I/O Errors

         IF(A_IO /= -1.0.OR.LEN_IO /= disk_length) THEN
! DEPENDS ON: ioerror
            CALL IOERROR('buffer out of boundary data',A_IO,LEN_IO,     &
     &                  disk_length)
            CMESSAGE=' intf_out: I/O ERROR '
            ICODE=51
            GO TO 999   !  Return
         END IF

         IADDR = IADDR+ len_data

        ENDDO          ! VAR



!L 5.    Close boundary output file if reinitialised during run
      IF (FT_STEPS(NFTOUT) >  0) THEN
        LEN_PPNAME=LEN(PPNAME)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,PPNAME,LEN_PPNAME,1,0,ICODE)
      END IF

!L 6.  Update FT_LASTFIELD
      FT_LASTFIELD(NFTOUT) = FT_LASTFIELD(NFTOUT) + 1

 999  RETURN
      END SUBROUTINE intf_out
#endif
