#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Setup for LBC coupling
Subroutine lbc_coup_setup(ltimer_lbc)

  Implicit None
  

! Description
!  Does the setup in preparation for LBC coupling.
!
! Method
!  The main function of this coupling code is to allow the driving model
! (generating the LBCs) and the LAM model (receiving the LBCs) to run in
! parallel. In a nutshell, this coupling ensures that the LBCs are
! available every time it reaches the stage when it needs to read in more
! LBCs. Communication between the global and LAM was done through a
! communication file - which the global wrote to and the LAM read from it.
!
! This functionality was used for a number of years in the operational
! suite to enable the global and the old UK mesoscale run alongside rather 
! than the UK Mes wait for the global to complete its run. It stopped being 
! used in the oper suite during 2006 when the NAE took over the old
! UK Mes, and the global/NAE main runs are not run in parallel. The
! global/NAE update runs are run in parallel but the NAE uses a LBC file
! from an earlier global run so no need for coupling.
!
! Current Code Owner:
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


  Logical,Intent(in) :: ltimer_lbc  ! Timer logical

#include "parvars.h"
#include "typsize.h"
  !L       Control for boundary updating.
#include "cmaxsize.h"
#include "cintfa.h"
#include "cbound.h"

  Integer           :: um_lbc_coup ! LBC Coupling Switch : 1/0 is on/off
  Character(len=80) :: filename    ! Filename of communication file.
  Character(len=8)  :: c_lbc_coup  ! Character variable to read env var
  Character(len=8)  :: ch_date2    ! Date returned from date_and_time
  Character(len=10) :: ch_time2    ! Time returned from date_and_time
  Character(len=5)  :: runid_char  ! RUNID for job
  Character(len=4)  :: runtype_char! Run Type (ie. NRUN, CRUN)
  Integer           :: lbc_ntimes  ! No of BCs in communication file.
  Integer           :: len_wait_tot! Total wait for availability of BCs
  Integer*8         :: isleep      ! Return value from SLEEP
  Integer*8         :: sleep       ! SLEEP function to make job wait
  Integer           :: j           ! Loop counter

  Logical           :: l_exist     !  T : Communication File exists
  Logical           :: l_active    !  T : Output stream active for LBCs.

#include "cprintst.h"
  ! Error reporting
  Integer           :: ICODE       ! =0 normal exit; >0 error exit
  Character(len=256):: Cmessage    ! Error message
  Character*(*)        RoutineName
  Parameter (RoutineName='LBC_SETUP_FILES')

#include "lbc_coup.h"

  ! Find out if LBC Coupling has been switched on this run
  ! from the env. variable UM_LBC_COUP.

  Call fort_get_env('UM_LBC_COUP',11,c_lbc_coup,8,icode)
  If (icode /= 0) Then
     um_lbc_coup = 0 ! No coupling
     If (PrintStatus  >=  PrStatus_Normal) Then
        Write (6,*) ' Env Var UM_LBC_COUP not set.'
        Write (6,*) ' Setting UM_LBC_COUP to ',um_lbc_coup
     End If
  Else
     Read(c_lbc_coup,'(i8)') um_lbc_coup
     If (PrintStatus  >=  PrStatus_Normal) Then
        Write (6,*) ' UM_LBC_COUP is set to ',UM_LBC_COUP
     End If
  End If
  If (um_lbc_coup == 0 .Or. um_lbc_coup == 1) Then
     l_lbc_coup = um_lbc_coup == 1
  Else
     Write (6,*) ' Invalid value given to UM_LBC_COUP ',             &
          UM_LBC_COUP
     Write (6,*) ' Valid values are 0 or 1'
     Write (6,*) ' L_LBC_COUP set to F. No LBC Coupling ',           &
          'in this run.'
     cmessage = 'U_MODEL : Invalid value given to UM_LBC_COUP'
     icode = 100
! DEPENDS ON: ereport
     Call Ereport(RoutineName,ICODE,Cmessage)
  End If

  If (PrintStatus  >=  PrStatus_Normal) Then
     If (l_lbc_coup) Then
        Write (6,*) ' LBC COUPLING switched on in this run.'
     Else
        Write (6,*) ' LBC COUPLING switched off in this run.'
     End If
  End If

  If (l_lbc_coup) Then

#if defined(ATMOS) && defined(GLOBAL)

     ! Find out which LBC output stream is providing the data
     ! from the env. variable UM_LBC_STREAM.

     Call fort_get_env('UM_LBC_STREAM',13,c_lbc_coup,8,icode)
     If (icode /= 0) Then
        um_lbc_stream = 0 ! No coupling
        If (PrintStatus  >=  PrStatus_Min) Then
           Write (6,*) ' gl : Env Var UM_LBC_STREAM not set.'
           Write (6,*) ' gl : Setting UM_LBC_STREAM to ',um_lbc_stream
        End If
     Else
        Read(c_lbc_coup,'(i8)') um_lbc_stream
        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : UM_LBC_STREAM is set to ',UM_LBC_STREAM
        End If
     End If

     ! Check validity of UM_LBC_STREAM
     
     If (um_lbc_stream <  1.Or.um_lbc_stream >  max_n_intf_a) Then
        Write (6,*) ' gl : UM_LBC_STREAM = ',UM_LBC_STREAM,           &
             ' is an invalid value.'
        Write (6,*) ' gl : Valid values are 1-',MAX_N_INTF_A
        cmessage = 'U_MODEL : Invalid value given to UM_LBC_STREAM'
        icode = 101
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ICODE,Cmessage)
     End If

     ! Check if this output stream is active.
     l_active = .False.
     Do j=1,n_intf_a
        l_active = l_active .Or. um_lbc_stream == lbc_stream_a(j)
     End Do
     If (.Not.l_active) Then
        Write (6,*) ' gl : Output LBC stream ',UM_LBC_STREAM,         &
             ' is inactive. Check UM_LBC_STREAM.'
        Write (6,*) ' gl : Active LBC streams are ',                  &
             (LBC_STREAM_A(j),j=1,n_intf_a)
        cmessage = 'U_MODEL : Output LBC stream is inactive.'
        icode = 101
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ICODE,Cmessage)
     End If


#endif
#if defined(ATMOS) && !defined(GLOBAL)

     ! Find out how long the mesoscale is to wait if there
     ! are insufficient boundary conditions to proceed.
     
     Call fort_get_env('UM_LBC_WAIT',11,c_lbc_coup,8,icode)
     If (icode /= 0) Then
        um_lbc_wait = 0 ! No waiting
        If (PrintStatus  >=  PrStatus_Min) Then
           Write (6,*) ' ms : Env Var UM_LBC_WAIT not set.'
           Write (6,*) ' ms : Setting UM_LBC_WAIT to ',um_lbc_wait
        End If
     Else
        Read(c_lbc_coup,'(i8)') um_lbc_wait
        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' ms : UM_LBC_WAIT is set to ',um_lbc_wait
        End If
     End If
     
     ! Find out maximum wait if there are insufficient
     ! boundary conditions to proceed.
     
     Call fort_get_env('UM_LBC_WAIT_MAX',15,c_lbc_coup,8,icode)
     If (icode /= 0) Then
        um_lbc_wait_max = 0 ! No waiting
        If (PrintStatus  >=  PrStatus_Min) Then
           Write (6,*) ' ms : Env Var UM_LBC_WAIT_MAX not set.'
           Write (6,*) ' ms : Setting UM_LBC_WAIT_MAX to ',       &
                um_lbc_wait_max
        End If
     Else
        Read(c_lbc_coup,'(i8)') um_lbc_wait_max
        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' ms : UM_LBC_WAIT_MAX is set to ',UM_LBC_WAIT_MAX
        End If
     End If
#endif

  End If  !  if l_lbc_coup

  ICODE=0  ! (incase FORT_GET_ENV has left with non-zero value)
#if defined(ATMOS) && defined(GLOBAL)

  If (l_lbc_coup) Then

     If (mype == 0) Then

        ! Get filename attached to Unit 190
        Call GET_FILE(190,filename,80,ICODE)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with GET_FILE',                  &
                ' for Unit No 190.'
           Write (6,*) ' gl : Return code from GET_FILE ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Error in GET_FILE for Unit No 190.'
           icode = 102
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : Filename for unit no 190 ',FILENAME
        End If

        ! Open the file with WRITE permission.
        Open(unit=190,FILE=filename,action="write",iostat=icode)
        
        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with OPEN for Unit 190.'
           Write (6,*) ' gl : Return code from OPEN ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with OPEN for Unit No 190.'
           icode = 103
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : File(unit 190) has been opened.'
        End If

        ! Send info to meso that L_LBC_COUP=T in Global.
        
        lbc_ntimes = 1000 + um_lbc_stream
        Write (190,*) lbc_ntimes

! DEPENDS ON: um_fort_flush
        Call um_fort_flush (190,icode)

        If (icode /= 0) Then
           Write (6,*) 'Return Code from FLUSH ',icode
           icode = 104
           Write (cmessage,*) 'U_MODEL : Error flushing out '//        &
                'contents for Unit 190.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : lbc_ntimes ', lbc_ntimes,                  &
                ' sent to LBC_FILE.'
        End If

        ! Unit 191 : File with information for operators to
        ! monitor progress with Boundary Conditions generated.
           
        ! Get filename attached to Unit 191
        Call GET_FILE(191,filename,80,ICODE)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with GET_FILE for Unit 191.'
           Write (6,*) ' gl : Return code from GET_FILE ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Error in GET_FILE for Unit No 191.'
           icode = 105
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : Filename for unit no 191 ',FILENAME
        End If
           
        ! Open the file with WRITE permission only.
        Open(unit=191,FILE=filename,action="write",iostat=icode)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with OPEN for Unit 191.'
           Write (6,*) ' gl : Return code from OPEN ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with OPEN for Unit No 191.'
           icode = 106
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Write (6,*) ' gl : File opened on unit 191.'

        ! Send RUNID and date to file on unit 191
        Call fort_get_env('RUNID',5,runid_char,5,icode)
        If (icode /= 0) Then
           Write (6,*) ' Problem with FORT_GET_ENV for RUNID.'
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with FORT_GET_ENV for RUNID.'
           icode = 107
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Call fort_get_env('TYPE',4,runtype_char,5,icode)
        If (icode /= 0) Then
           Write (6,*) ' Problem with FORT_GET_ENV for TYPE.'
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with FORT_GET_ENV for TYPE.'
           icode = 108
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Call Date_and_time(ch_date2, ch_time2)
        Write (191,*) ' RUNID : ',Trim(runid_char),                   &
             '  RUN TYPE : ',Trim(runtype_char),                      &
             '  on ',ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4)
        
! DEPENDS ON: um_fort_flush
        Call um_fort_flush (191,icode)
        
        If (icode /= 0) Then
           Write (6,*) 'Return Code from FLUSH ',icode
           icode = 109
           Write (cmessage,*) 'U_MODEL : Error flushing out '//        &
                'contents for Unit 191.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

     End If  !  if mype == 0
        
  End If  !   if l_lbc_coup

#endif
#if defined(ATMOS) && !defined(GLOBAL)

  If (l_lbc_coup) Then

! DEPENDS ON: timer
     If (LTIMER_LBC) Call Timer ('LBC_WAIT',3)

     If (mype == 0) Then

        ! Get filename attached to Unit 190
        Call GET_FILE(190,lbc_filename,80,ICODE)
        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' ms : Filename from GET_FILE ',lbc_filename
        End If

        If (icode /= 0) Then
           Write (6,*) ' Return code from GET_FILE ',icode
           icode = 600
           Write (cmessage,*) 'U_MODEL : Problem with GET_FILE '//     &
                'for Unit No 190.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Call Date_and_time(ch_date2, ch_time2)

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write(6,*) 'LBC_COUP: ',                                      &
                ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',&
                ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),       &
                ' Wait to see if file exists.'
        End If

        len_wait_tot = 0
149     Continue

        ! Check that the file exists.
        Inquire (file=lbc_filename,exist=l_exist,iostat=icode)

        If (l_exist) Then  !  file exists

           Call Date_and_time(ch_date2, ch_time2)
           If (PrintStatus  >=  PrStatus_Normal) Then
              Write(6,*) 'LBC_COUP: ',                                      &
                   ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',&
                   ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),       &
                   ' File exists - Proceed to open file.'
           End If
      
           ! Open the file with READ ONLY permission.
           Open (unit=190,file=lbc_filename,action="read",             &
                iostat=icode)

           ! Check return code from OPEN.
           If (icode /= 0) Then
              Write (6,*) ' Return code from OPEN ',icode
              icode = 601
              Write (cmessage,*) 'U_MODEL : Problem with OPEN '//       &
                   'for Unit No 190.'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ICODE,Cmessage)
           End If

        Else  !  file does not exist

           If (len_wait_tot >= um_lbc_wait_max) Then

              ! Maximum wait time has been reached/exceeded.

              Write (6,*) ' ms : lbc_file does not exist.'
              Write (6,*) ' ms : Maximum wait time reached',            &
                   ' after ',um_lbc_wait_max,' seconds.'
              icode = 602
              cmessage = 'U_MODEL : LBC_FILE does not exist.'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ICODE,Cmessage)

           Else

              ! Wait for um_lbc_wait seconds before another attempt
              ! to see if file exists.

              If (PrintStatus  >=  PrStatus_Normal) Then
                 Write (6,*) ' ms : lbc_file does not exist yet.'
                 Write (6,*) ' ms : wait for ',um_lbc_wait,              &
                      ' seconds and retry.'
              End If
              isleep = sleep(um_lbc_wait)
              len_wait_tot = len_wait_tot+um_lbc_wait
              If (PrintStatus  >=  PrStatus_Normal) Then
                 Write (6,*) ' ms : Total Wait so far ',len_wait_tot,    &
                      ' seconds.'
              End If
              Goto 149   !  Retry to see if LBC_FILE exists
              
           End If

        End If  ! if l_exist

     End If  !  if mype=0


     ! Check that LBC_COUPLING has been switched on in Global.
     If (mype == 0) Then

        len_wait_tot = 0
150     Continue

        ! Close the communication file and re-open
        Close(190)
        Open (190,file=lbc_filename,action="read",iostat=icode)

        ! Check retrun code from OPEN
        If (icode /= 0) Then
           Write (6,*) ' Return code from OPEN ',icode
           icode = 603
           Write (cmessage,*) 'U_MODEL : Problem with OPEN '//         &
                'for Unit No 190.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        ! Read in the first value
        Read (190,*,iostat=icode) lbc_ntimes

        ! Check return code from READ
        If (icode /= 0) Then

           If (PrintStatus  >=  PrStatus_Normal) Then
              Write (6,*) ' ms : Return code from READ ',icode
           End If

           If (len_wait_tot >= um_lbc_wait_max) Then

              ! Maximum wait time has been reached or exceeded.
              ! Give up waiting and abort.

              Write (6,*) ' ms : Required LBC_NTIMES not read in',      &
                   ' after ',um_lbc_wait_max,' seconds.'
              icode = 604
              cmessage = 'U_MODEL : Required LBC_NTIMES '//             &
                      'not found in LBC_FILE.'
! DEPENDS ON: ereport
              Call Ereport(RoutineName,ICODE,Cmessage)

           Else

              ! Wait for um_lbc_wait seconds abefore another attempt
              ! to read a value.
              
              If (PrintStatus  >=  PrStatus_Normal) Then
                 Write (6,*) ' ms : wait for ',um_lbc_wait,              &
                      ' seconds and retry.'
              End If
              isleep = sleep(um_lbc_wait)
              len_wait_tot = len_wait_tot+um_lbc_wait
              If (PrintStatus  >=  PrStatus_Normal) Then
                 Write (6,*) ' ms : Total Wait so far ',len_wait_tot,    &
                      ' seconds.'
              End If
              goto 150   !  Retry to see if required LBC_NTIMES exists

           End If

        End If  !  if icode /= 0

        ! The first value in the file is >1000.
        If (lbc_ntimes >  1000) Then
           um_lbc_stream = lbc_ntimes - 1000
           If (PrintStatus  >=  PrStatus_Normal) Then
              Write (6,*) ' ms : l_lbc_coup = T in Global'
              Write (6,*) ' ms : global output stream is ',um_lbc_stream
           End If
        End If  !  if l_lbc_ntimes

     End If  !  if mype == 0

! DEPENDS ON: timer
     If (LTIMER_LBC) Call Timer ('LBC_WAIT',4)

  End If  !  if l_lbc_coup

#endif

  Return
End Subroutine lbc_coup_setup
#endif
