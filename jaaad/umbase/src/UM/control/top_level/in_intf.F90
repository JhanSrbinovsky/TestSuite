#if defined(CONTROL) && defined(ATMOS) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine IN_INTF
!LL
!LL Purpose : Takes as input, codes set by the user interface defining
!LL           the start time, end time, and interval for creating
!LL           interface data for a limited area model, and data
!LL           defining the limited area grid. The source model may also
!LL           be limited area. Sets up fixed length, integer & real
!LL           headers and level dependent constants for the interface
!LL           data set. All prognostic variables for which horizontal
!LL           differencing is performed are included, ie all tracers
!LL           but no surface fields.
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL Logical components covered : D810
!LL
!LL System task : D81
!LL
!LL Documentation : Unified Model Documentation Paper No D8
!LL
!LLEND -------------------------------------------------------------

      SUBROUTINE IN_INTF (                                              &

!*L   Arguments:
#include "argduma.h"
#include "arginfa.h"
     &           NFTOUT,ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typinfa.h"

      INTEGER                                                           &
     &       NFTOUT,                                                    &
                       ! FTN Number to write interface data
     &       ICODE     ! Return code

      CHARACTER*80                                                      &
     &       CMESSAGE  ! Error message
!*

#include "chsunits.h"
#include "csubmodl.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "c_mdi.h"
#include "cntl_io.h"

!    Local variables
      INTEGER                                                           &
     &        I,J,IJ,                                                   &
                       ! DO loop indices
     &        ROW,                                                      &
     &        IRIM,                                                     &
     &        DAYS,SECS,                                                &
                         ! Time difference relative to reference point
     &        LEVEL,                                                    &
     &        LEN_PP,                                                   &
                       ! Total length of PP headers
     &        LEN_IO,                                                   &
                       ! Total length of data read in
     &        NTIMES,                                                   &
                       ! Number of times for which interface data
!                      ! is required.
     &        ZERO,                                                     &
     &        JINTF          ! Interface area index
!L----------------------------------------------------------------------

      REAL                                                              &
     &       A_IO
      REAL, ALLOCATABLE :: Dumm1(:)   !dummy array to fill file

      INTEGER YY,MM,DD,HH,MN,SEC,DAY_NO
      INTEGER A_SECS_PER_STEP ! secs per step for atmos. sub_model
      INTEGER A_INTF_FREQ_SECS  ! Interface frequency in seconds
      Integer :: ntime
      Integer :: lookup_start

      Integer lbc_len_inthd, lbc_len_realhd, lbc_lookupsa
      CHARACTER*80 STRING           ! MODEL_FT_UNIT value
      CHARACTER*14 PPNAME           ! Boundary data name
      Character (Len=*), Parameter :: RoutineName = 'IN_INTF'

!L    Internal Structure

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER ('IN_INTF',3)

#if defined(ATMOS)
      IF (LLBOUTim(A_IM)) THEN

! DEPENDS ON: intf_area
      CALL INTF_AREA (A_IM,NFTOUT,JINTF)

      A_INTF_FREQ_SECS=3600*A_INTF_FREQ_HR(JINTF) +                     &
     &  60*A_INTF_FREQ_MN(JINTF)+A_INTF_FREQ_SC(JINTF)
      IF (A_INTF_FREQ_SECS >  0) THEN

        IF (LNEWBND(JINTF)) THEN

!L      LNEWBND = true (New dataset to be set up)

! 1.1   Compute secs per step for atmosphere sub_model
        A_SECS_PER_STEP = SECS_PER_PERIODim(a_im)/                      &
     &                    STEPS_PER_PERIODim(a_im)

!L  1.0 Set up headers

        IF (FT_STEPS(NFTOUT) == 0) THEN
          NTIMES = ((A_INTF_END_HR(JINTF)-A_INTF_START_HR(JINTF))*3600)/&
     &              A_INTF_FREQ_SECS+1
        ELSE ! reinitialisation
          NTIMES = FT_STEPS(NFTOUT)/                                    &
     &               (A_INTF_FREQ_SECS / A_SECS_PER_STEP)
          IF (STEPim(a_im)-1 <= FT_FIRSTSTEP(NFTOUT)) THEN ! 1st file
            NTIMES = NTIMES + 1
          ENDIF
        ENDIF

      If (lbc_nd(jintf) == 1) Then
        lbc_len_inthd  = pp_len_inthd
        lbc_len_realhd = pp_len_realhd
        lbc_lookupsa   = intf_lookupsa
      Else
        lbc_len_inthd  = 15
        lbc_len_realhd = 6
        lbc_lookupsa   = old_intf_lookupsa
      Endif
!L  1.1 Fixed length header
!       Copy main dump header as first step

      DO I=1,LEN_FIXHD
        FIXHD_INTFA(I,JINTF)=A_FIXHD(I)
      ENDDO

!        Set boundary data identifiers

      FIXHD_INTFA(1,JINTF)  = IMDI
      FIXHD_INTFA(5,JINTF)  = 5
      FIXHD_INTFA(10,JINTF) = 1
      If (lbc_nd(jintf) == 0) Then
        Fixhd_Intfa(3,jintf)  = 1
        Fixhd_Intfa(9,jintf)  = 2
        Fixhd_Intfa(12,jintf) = 405
      End If

!        Modify individual items

!L      Calculate start validity time

      IF (FT_STEPS(NFTOUT) == 0) THEN
        SECS = A_INTF_START_HR(JINTF)*3600
        DAYS = 0
      ELSE
       IF (STEPim(a_im)-1 <= FT_FIRSTSTEP(NFTOUT)) THEN   ! First file
! DEPENDS ON: stp2time
        CALL STP2TIME(FT_FIRSTSTEP(NFTOUT),STEPS_PER_PERIODim(a_im),    &
     &                SECS_PER_PERIODim(a_im),DAYS,SECS)
       ELSE ! not first file
! DEPENDS ON: stp2time
        CALL STP2TIME(STEPim(a_im)-1,STEPS_PER_PERIODim(a_im),          &
     &                SECS_PER_PERIODim(a_im),DAYS,SECS)
        SECS=SECS+ A_INTF_FREQ_SECS
       ENDIF
      ENDIF

! DEPENDS ON: sec2time
      CALL SEC2TIME(DAYS,SECS,BASIS_TIME_DAYS,BASIS_TIME_SECS,          &
     &              YY,MM,DD,HH,MN,SEC,DAY_NO,LCAL360)

      FIXHD_INTFA(21,JINTF) = YY
      FIXHD_INTFA(22,JINTF) = MM
      FIXHD_INTFA(23,JINTF) = DD
      FIXHD_INTFA(24,JINTF) = HH
      FIXHD_INTFA(25,JINTF) = MN
      FIXHD_INTFA(26,JINTF) = SEC
      FIXHD_INTFA(27,JINTF) = DAY_NO

!     Data interval
      FIXHD_INTFA(35,JINTF) = 0
      FIXHD_INTFA(36,JINTF) = 0
      FIXHD_INTFA(37,JINTF) = A_INTF_FREQ_HR(JINTF)/24
      FIXHD_INTFA(38,JINTF) = MOD(A_INTF_FREQ_HR(JINTF),24)
      FIXHD_INTFA(39,JINTF) = A_INTF_FREQ_MN(JINTF)
      FIXHD_INTFA(40,JINTF) = A_INTF_FREQ_SC(JINTF)
      FIXHD_INTFA(41,JINTF) = A_INTF_FREQ_HR(JINTF)/24

!L      Calculate last validity time

      IF (FT_STEPS(NFTOUT) == 0) THEN
        SECS=A_INTF_END_HR(JINTF)*3600
        DAYS = 0
      ELSE
       IF (STEPim(a_im)-1 <= FT_FIRSTSTEP(NFTOUT)) THEN   ! First file
! DEPENDS ON: stp2time
        CALL STP2TIME(FT_FIRSTSTEP(NFTOUT),STEPS_PER_PERIODim(a_im),    &
     &                SECS_PER_PERIODim(a_im),DAYS,SECS)
        SECS=SECS + FT_STEPS(NFTOUT)*A_SECS_PER_STEP
       ELSE ! not first file
! DEPENDS ON: stp2time
        CALL STP2TIME(STEPim(a_im)-1+FT_STEPS(NFTOUT),                  &
     &      STEPS_PER_PERIODim(a_im),SECS_PER_PERIODim(a_im),DAYS,SECS)
       ENDIF
      ENDIF

! DEPENDS ON: sec2time
      CALL SEC2TIME(DAYS,SECS,BASIS_TIME_DAYS,BASIS_TIME_SECS,          &
     &              YY,MM,DD,HH,MN,SEC,DAY_NO,LCAL360)

      FIXHD_INTFA(28,JINTF) = YY
      FIXHD_INTFA(29,JINTF) = MM
      FIXHD_INTFA(30,JINTF) = DD
      FIXHD_INTFA(31,JINTF) = HH
      FIXHD_INTFA(32,JINTF) = MN
      FIXHD_INTFA(33,JINTF) = SEC
      FIXHD_INTFA(34,JINTF) = DAY_NO

!L      Modify header lengths

      FIXHD_INTFA(101,JINTF)=lbc_len_inthd
      FIXHD_INTFA(105,JINTF) =                                          &
     &FIXHD_INTFA(100,JINTF)+FIXHD_INTFA(101,JINTF)
      FIXHD_INTFA(106,JINTF)=lbc_len_realhd
      FIXHD_INTFA(110,JINTF)=                                           &
     &FIXHD_INTFA(105,JINTF)+FIXHD_INTFA(106,JINTF)

!L      Set length of level dependent constants

      FIXHD_INTFA(112,JINTF)=4
      IF (INTF_VERT_INTERP(JINTF)) THEN
        If (lbc_nd(jintf) == 1) Then
          FIXHD_INTFA(111,JINTF)=INTF_P_LEVELS(JINTF)+1
        Else
          FIXHD_INTFA(111,JINTF)=INTF_P_LEVELS(JINTF)
        End If
      ELSE
        If (lbc_nd(jintf) == 1) Then
          FIXHD_INTFA(111,JINTF)=MODEL_LEVELS+1
        Else
          FIXHD_INTFA(111,JINTF)=MODEL_LEVELS
        End If
      END IF

!  NO row and column dependent constants
      FIXHD_INTFA(115,JINTF)=IMDI
      FIXHD_INTFA(116,JINTF)=IMDI
      FIXHD_INTFA(117,JINTF)=IMDI
      FIXHD_INTFA(120,JINTF)=IMDI
      FIXHD_INTFA(121,JINTF)=IMDI
      FIXHD_INTFA(122,JINTF)=IMDI
      
      IF (INTF_L_VAR_LBC(JINTF)) THEN
        FIXHD_INTFA(115,JINTF)=                                         &
       &        FIXHD_INTFA(110,JINTF)+                                 &
       &        FIXHD_INTFA(111,JINTF)*FIXHD_INTFA(112,JINTF)
        FIXHD_INTFA(116,JINTF)= INTF_P_ROWS(JINTF)
        FIXHD_INTFA(117,JINTF)=2
        FIXHD_INTFA(120,JINTF)=                                         &
       &        FIXHD_INTFA(115,JINTF)+                                 &
       &        FIXHD_INTFA(116,JINTF)*FIXHD_INTFA(117,JINTF)
        FIXHD_INTFA(121,JINTF)=INTF_ROW_LENGTH(JINTF)
        FIXHD_INTFA(122,JINTF)=2
      END IF
      
!  NO field_constants,extra constants,temp_historyfile or compressed
!  indexes
      DO I=125,145
        FIXHD_INTFA(I,JINTF)=IMDI
      ENDDO

!     Start address and length of look up table
      FIXHD_INTFA(150,JINTF) = FIXHD_INTFA(110,JINTF) +                 &
     &                  FIXHD_INTFA(111,JINTF) * FIXHD_INTFA(112,JINTF)

      IF (INTF_L_VAR_LBC(JINTF)) THEN
        FIXHD_INTFA(150,JINTF) = FIXHD_INTFA(120,JINTF) +               &
     &                  FIXHD_INTFA(121,JINTF) * FIXHD_INTFA(122,JINTF)
      END IF
       
      FIXHD_INTFA(152,JINTF) = NTIMES * LBC_LOOKUPSA
      FIXHD_INTFA(153,JINTF) = IMDI

!     Start address and length of data section
!--make sure the data starts on a sector bndry
      fixhd_intfa(160, jintf)=((fixhd_intfa(150, jintf)+                &
     & fixhd_intfa(151, jintf)*fixhd_intfa(152, jintf)+                 &
     & um_sector_size-1)/um_sector_size)*um_sector_size+1
      FIXHD_INTFA(161,JINTF) = 0

!L  1.2 Integer header

      DO I=1,PP_LEN_INTHD
        INTHD_INTFA(I,JINTF) = A_INTHD(I)
      ENDDO

      INTHD_INTFA(1,JINTF) = INTERFACE_FSTEPim(JINTF,a_im)
      INTHD_INTFA(2,JINTF) = A_INTF_FREQ_HR(JINTF)
      INTHD_INTFA(3,JINTF) = NTIMES
      INTHD_INTFA(6,JINTF) = INTF_ROW_LENGTH(JINTF)
      INTHD_INTFA(7,JINTF) = INTF_P_ROWS(JINTF)
      INTHD_INTFA(8,JINTF) = INTF_P_LEVELS(JINTF)
      INTHD_INTFA(9,JINTF) = INTF_Q_LEVELS(JINTF)
      INTHD_INTFA(12,JINTF)= INTF_TR_LEVELS(JINTF)
      INTHD_INTFA(15,JINTF)= LBC_LOOKUPSA

      If (lbc_nd(jintf) == 1) Then

      Do I=16,46
        INTHD_INTFA(I,JINTF) = IMDI
      End Do

!     Algorithm used for generating heights.
      INTHD_INTFA(17,JINTF) = 2

!     First rho level at which height is constant
      If (Intf_Vert_Interp(jintf)) Then
        INTHD_INTFA(24,JINTF) = LBC_FIRST_R_RHO(JINTF)
      Else
        INTHD_INTFA(24,JINTF) = A_INTHD(24)
      End If

      End If

!L  1.3 Real header

      DO I=1,PP_LEN_REALHD
       REALHD_INTFA(I,JINTF) = A_REALHD(I)
      ENDDO

      REALHD_INTFA(1,JINTF) = INTF_EWSPACE(JINTF)
      REALHD_INTFA(2,JINTF) = INTF_NSSPACE(JINTF)
      REALHD_INTFA(3,JINTF) = INTF_FIRSTLAT(JINTF)
      REALHD_INTFA(4,JINTF) = INTF_FIRSTLONG(JINTF)
      REALHD_INTFA(5,JINTF) = INTF_POLELAT(JINTF)
      REALHD_INTFA(6,JINTF) = INTF_POLELONG(JINTF)
      If (lbc_nd(jintf) == 1) Then

      Do I=7,38
        REALHD_INTFA(I,JINTF) = RMDI
      End Do

!     Height at top of model
      If (Intf_Vert_Interp(jintf)) Then
        REALHD_INTFA(16,JINTF)= LBC_Z_TOP_MODEL(JINTF)
      Else
        REALHD_INTFA(16,JINTF)= A_REALHD(16)
      End If

      End If


!L  1.4 Level dependent constants

        IF (INTF_VERT_INTERP(JINTF)) THEN
!L      Set level dependent constants from namelist INTFCNST  if
!L      vertical interpolation required

          If (lbc_nd(jintf) == 1) Then

          levdepc_intfa(:,:,jintf) = rmdi

          DO LEVEL=1,INTF_P_LEVELS(JINTF)+1
            LEVDEPC_INTFA(LEVEL,1,JINTF) = LBC_ETA_THETA(LEVEL,JINTF)
          ENDDO
          DO LEVEL=1,INTF_P_LEVELS(JINTF)
            LEVDEPC_INTFA(LEVEL,2,JINTF) = LBC_ETA_RHO(LEVEL,JINTF)
          ENDDO

          Else

          DO LEVEL=1,INTF_P_LEVELS(JINTF)
            LEVDEPC_INTFA(LEVEL,1,JINTF) = INTF_AK(LEVEL,JINTF)
            LEVDEPC_INTFA(LEVEL,2,JINTF) = INTF_BK(LEVEL,JINTF)
            LEVDEPC_INTFA(LEVEL,3,JINTF) =                              &
     &              INTF_AKH(LEVEL+1,JINTF)-INTF_AKH(LEVEL,JINTF)
            LEVDEPC_INTFA(LEVEL,4,JINTF) =                              &
     &      INTF_BKH(LEVEL+1,JINTF)-INTF_BKH(LEVEL,JINTF)
          ENDDO

          Endif

        ELSE

!       Copy level dependent constants from source model

          LevDepC_Intfa(:,:,jintf) = rmdi

!         Eta - theta levels
          Do Level = 1, Model_Levels+1
            LevDepC_Intfa(Level,1,jintf) = A_LevDepc(Level)
          End Do

!         Eta - rho levels
          Do Level = 1, Model_Levels
            LevDepC_Intfa(Level,2,jintf) =                              &
     &      A_LevDepc(Model_Levels+1+Level)
          End Do

!         RH_Crit and Soil Moisture levels not copied.


        END IF

!L  1.5 Row/Col dependent constants
        
        rowdepc_intfa(:,:,jintf) = rmdi
        coldepc_intfa(:,:,jintf) = rmdi
        
        IF (INTF_L_VAR_LBC(JINTF)) THEN
 
          DO I=1, INTF_P_ROWS(JINTF)
            ROWDEPC_INTFA(I,1,JINTF) =  PHI_INTF_P(I, JINTF)
            ROWDEPC_INTFA(I,2,JINTF) =  PHI_INTF_V(I, JINTF)
          ENDDO
          DO J=1, INTF_ROW_LENGTH(JINTF)
            COLDEPC_INTFA(J,1,JINTF) =  LAMBDA_INTF_P(J, JINTF)
            COLDEPC_INTFA(J,2,JINTF) =  LAMBDA_INTF_U(J, JINTF)
          ENDDO 
        END IF
        
!L  2.0 Write out headers

!L  2.1 Fixed length header

! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FIXHD_INTFA(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of fixed length header',A_IO,LEN_IO, &
     &                  LEN_FIXHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=1
          RETURN
        END IF

!L  2.2 Integer header

! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,INTHD_INTFA(1,JINTF),                       &
     &               LBC_LEN_INTHD,LEN_IO,A_IO)

!       Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LBC_LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of integer header',A_IO,LEN_IO,      &
     &                  LBC_LEN_INTHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=2
          RETURN
        END IF

!L  2.3 Real header

! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,REALHD_INTFA(1,JINTF),                      &
     &               LBC_LEN_REALHD,LEN_IO,A_IO)

!       Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LBC_LEN_REALHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of real header',A_IO,LEN_IO,         &
     &                  LBC_LEN_REALHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=3
          RETURN
        END IF

!L  2.4 Level dependent constants

!       Write out each variable separately as second dimension
!       of LEVDEPC_INTFA is now a maximum dimension for all
!       interface areas being generated

        DO I=1,4

! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LEVDEPC_INTFA(1,I,JINTF),                   &
     &               FIXHD_INTFA(111,JINTF),LEN_IO,A_IO)

!       Check for I/O Errors

        IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(111,JINTF)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of level dependent constants',A_IO,  &
     &                  LEN_IO,FIXHD_INTFA(111,JINTF))
          CMESSAGE='IN_INTF : I/O ERROR'
          ICODE=4
          RETURN
        END IF

        ENDDO

!L  2.5 Row/Col dependent constants

        IF (INTF_L_VAR_LBC(JINTF)) THEN
! ROW first        
          DO I=1,2
        
! DEPENDS ON: buffout
            CALL BUFFOUT(NFTOUT,ROWDEPC_INTFA(1,I,JINTF),               &
     &                   FIXHD_INTFA(116,JINTF),LEN_IO,A_IO)

!           Check for I/O Errors

            IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(116,JINTF)) THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of row dependent constants',     &
     &                      A_IO,LEN_IO,FIXHD_INTFA(116,JINTF))
              CMESSAGE='IN_INTF : I/O ERROR'
              ICODE=5
              RETURN
            END IF

          ENDDO

! ColDepC        

          DO I=1,2
        
! DEPENDS ON: buffout
            CALL BUFFOUT(NFTOUT,COLDEPC_INTFA(1,I,JINTF),               &
     &                   FIXHD_INTFA(121,JINTF),LEN_IO,A_IO)

!           Check for I/O Errors

            IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(121,JINTF)) THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of col dependent constants',     &
     &                      A_IO,LEN_IO,FIXHD_INTFA(121,JINTF))
              CMESSAGE='IN_INTF : I/O ERROR'
              ICODE=6
              RETURN
            END IF

          ENDDO 
        END IF     ! INTF_L_VAR_LBC(JINTF) 

!L  2.6 Write dummy record to reserve space for PP headers

        If (lbc_nd(jintf) == 1) Then
          LEN_PP =  (INTHD_INTFA(3,JINTF)-1)  * ( LBC_LOOKUPSA-1 )      &
     &            + LBC_LOOKUPSA                                        &
     &            + 1
        Else
          LEN_PP=INTHD_INTFA(3,JINTF)*INTF_LOOKUPSA+1
        End If

        IF(LEN_PP >  LEN_TOT) THEN
          CMESSAGE='IN_INTF:Insufficient space for PP headers'
          ICODE=5
          RETURN
        END IF

        allocate (Dumm1(LEN_PP))
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,Dumm1(1),LEN_PP,LEN_IO,A_IO)
        deallocate (Dumm1)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_PP) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of dummy PP headers',A_IO,LEN_IO,    &
     &                  LEN_PP)
          CMESSAGE='IN_INTF=I/O ERROR'
          ICODE=6
          RETURN
        END IF

!    Remainder of headers not used

       ELSE
!L      LNEWBND = False (dataset exists)

!L  3.0 Read in headers
!L      If reinitialised boundary output file to be processed
!L      then open it.
!L
        IF (FT_STEPS(NFTOUT) >  0) THEN
        STRING=MODEL_FT_UNIT(NFTOUT)
        PPNAME=STRING(18:31)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,PPNAME,14,1,1,ICODE)
        IF (ICODE /= 0) THEN
          CMESSAGE="IN_INTF: Error opening preassigned boundary file"
          GO TO 999   !  Return
        ENDIF
!
        ENDIF

!L  3.1  Fixed length header

! DEPENDS ON: buffin
        CALL BUFFIN(NFTOUT,FIXHD_INTFA(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)
! Check for I/O Errors
         IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of fixed length header',A_IO,LEN_IO, &
     &                   LEN_FIXHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=7
           RETURN
         ENDIF

!L  3.2  Integer header

! DEPENDS ON: buffin
         CALL BUFFIN(NFTOUT,INTHD_INTFA(1,JINTF),PP_LEN_INTHD,          &
     &               LEN_IO,A_IO)

! Check for I/O Errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_LEN_INTHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of integer header',A_IO,LEN_IO,      &
     &                   PP_LEN_INTHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=8
           RETURN
         END IF

!L  3.3  Real header

! DEPENDS ON: buffin
         CALL BUFFIN(NFTOUT,REALHD_INTFA(1,JINTF),PP_LEN_REALHD,        &
     &               LEN_IO,A_IO)

! Check for I/O Errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_LEN_REALHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of real header',A_IO,LEN_IO,         &
     &                   PP_LEN_REALHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=9
           RETURN
         END IF

!   3.4  Look up Headers

!        Read in just one header ; this will be the last header that
!        corresponds to the last LBC variable before the start of a
!        CRUN. Required to continue disk and start addresses in LOOKUP.

         ntime = ft_lastfield(nftout)

         If (ntime > 0) Then

           lookup_start = fixhd_intfa(150,jintf) +                      &
     &                    fixhd_intfa(151,jintf) *                      &
     &                   (intf_lookupsa + (ntime-1)*(intf_lookupsa-1))

!          Point to last header
! DEPENDS ON: setpos
           Call Setpos (nftout,lookup_start-len1_lookup-1,icode)

           If (ICode /= 0) Then
             CMESSAGE = 'Error with SETPOS for LBC Lookup Table'
             ICode    = 10
! DEPENDS ON: ereport
             Call EReport (RoutineName, ICode, CMessage)
           End If

!          Read in one lookup header
! DEPENDS ON: buffin
           Call Buffin (nftout,lookup_intfa(1,intf_lookupsa,jintf),     &
     &                  len1_lookup,len_io,a_io)

           If (A_IO /= -1.0 .or. LEN_IO /= len1_lookup) Then
! DEPENDS ON: ioerror
             CALL IOERROR('buffer in of lookup header',A_IO,LEN_IO,     &
     &                   len1_lookup)
             CMESSAGE = 'I/O ERROR with buffin of LBC Lookup Header'
             ICode    = 11
! DEPENDS ON: ereport
             Call EReport (RoutineName, ICode, CMessage)
           End If

         Endif   !  If ntime > 0

!L Reset LNEWBND to true after reading in header information to allow
!L writing of new headers for subsequent reinitialised files.
!
         LNEWBND(JINTF) = .TRUE.
       END IF

       If (LBC_ND(JINTF) == 0) Then   !  For 4.5 LBCs only

!L  4.0 Calculate the interpolation cofficients for interface area JINTF
! DEPENDS ON: timer
         IF (LTIMER) CALL TIMER ('INTF_HIC',3)

! DEPENDS ON: intf_hintc
         CALL INTF_HINTC (                                              &
#include "argduma.h"
#include "arginfa.h"
     &   JINTF,LEN_INTFA_P(JINTF),LEN_INTFA_U(JINTF),                   &
     &   ICODE,CMESSAGE,LLBOUTim(A_IM))

! DEPENDS ON: timer
         IF (LTIMER) CALL TIMER ('INTF_HIC',4)

         IF (ICODE /= 0) THEN
           CMESSAGE = 'IN_INTF : Error in routine INTF_HINTC'
           GO TO 999   !  Return
         ENDIF

       End If
      END IF

      END IF         !  LLBOUTim(A_IM)
#endif


!L  5.0  End of routine

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER ('IN_INTF',4)
 999  RETURN
      END SUBROUTINE IN_INTF

!-----------------------------------------------------------------------


#endif
