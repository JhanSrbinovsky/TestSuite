#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Update LBC coupling
Subroutine lbc_coup_update( &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
     submodel,ICODE,CMESSAGE)

  Implicit None

! Description
!  Does an LBC coupling update
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


  Integer,           Intent(in   ) :: submodel! Submodel id 
  Integer,           Intent(  out) :: ICODE   ! =0 normal exit; >0 error exit
  Character(len=256),Intent(  out) :: Cmessage! Error message

#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"

#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typptra.h"
#include "typbnd.h"
#include "typsts.h"

#include "chsunits.h"
#include "cintfa.h"
#include "ccontrol.h"
#include "ctime.h"
#include "ppxlook.h"

  Integer           :: um_lbc_coup ! LBC Coupling Switch : 1/0 is on/off
  Character(len=8)  :: ch_date2    ! Date returned from date_and_time
  Character(len=10) :: ch_time2    ! Time returned from date_and_time
  Integer           :: lbc_ntimes  ! No of BCs in communication file.
  Integer           :: ms_ntimes   ! No of BCs required in mesoscale.
  Integer           :: len_wait_tot!  Total wait for availability of BCs
  Integer*8         :: isleep      !  Return value from SLEEP
  Integer*8         :: sleep       !  SLEEP function to make job wait

  Logical           :: l_active    !  T : Output stream active for LBCs.

#include "cprintst.h"
  ! Error reporting
  Character*(*)        RoutineName
  Parameter (RoutineName='LBC_UPDATE')

#include "lbc_coup.h"
  
#if defined(ATMOS) && !defined(GLOBAL)
  
  ! Swap to second atmos boundary file?
  If (Num_ALBCs == 2 .And. STEPim(a_im) == ALBC_SwapStep) Then
     ALBC_num = 2
     If (PrintStatus >= PrStatus_Normal) Then
        Write (6,*) ''
        Write (6,*) 'U_MODEL: Swapping to 2nd atmos boundary '//&
             'file'
        Write (6,*) '         Step = ', STEPim(a_im)
        Write (6,*) ''
     End If
  End If

  ! NOTE: If using two atmos boundary files, coupling is only
  !       activated for the second.

  If (l_lbc_coup .And..Not.(Num_ALBCs == 2 .And. ALBC_num == 1)) Then

     ! The boundary conditions (BCs) are updated every N hours.
     ! The BCs required to proceed N hours are read in. If the
     ! model has run M hours, then BCs must be available at
     ! least M+N hours.

     Call Date_and_time(ch_date2, ch_time2)

     If (PrintStatus  >=  PrStatus_Normal) Then
        Write(6,*)  'LBC_COUP: ',                                      &
             ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ', &
             ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),        &
             ' Wait to call INBOUND/UPBOUND in U_MODEL.'
     End If

     ! Determine which boundary data is required to proceed
     ! the next period.
     If (ALBC_num == 2) Then
        ms_ntimes = 1 + ( (stepim(a_im)-ALBC_SwapStep) &
             / boundary_stepsim(a_im) )
     Else
        ms_ntimes = 2 + (stepim(a_im)/boundary_stepsim(a_im))
     End If

     If (mype == 0) Then

        len_wait_tot = 0
160     Continue

        ! Close the communication file and re-open.
        Close(190)
        Open (190,file=lbc_filename,action="read",iostat=icode)
        
        ! Check return code from OPEN.
        If (icode /= 0) Then
           Write (6,*) ' Return code from OPEN ',icode
           icode = 701
           Write (cmessage,*) 'U_MODEL : Problem with OPEN '//   &
                'for Unit No 190.'
           ! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

161     Continue

        ! Read next value.
        Read (190,*,iostat=icode) lbc_ntimes

        ! Check return code from READ.
        If (icode /= 0) Then

           Write (6,*) ' ms : Return code from READ ',icode

           If (len_wait_tot >= um_lbc_wait_max) Then
              ! Maximum wait time has been reached or exceeded.
              Write (6,*) ' ms : Maximum wait time reached'//  &
                   ' after ',um_lbc_wait_max,' seconds.'
              icode = 702
              Write (cmessage,*) 'U_MODEL : Maximum wait time '// &
                   'reached while reading from LBC_FILE.'
              ! DEPENDS ON: ereport
              Call Ereport(RoutineName,ICODE,Cmessage)
           End If

           ! Wait for um_lbc_wait seconds before re-trying.
           
           If (PrintStatus  >=  PrStatus_Diag) Then
              Write (6,*) ' ms : Wait for ',um_lbc_wait,' seconds and retry.'
           End If
           isleep = sleep(um_lbc_wait)
           len_wait_tot = len_wait_tot+um_lbc_wait
           If (PrintStatus  >=  PrStatus_Diag) Then
              Write (6,*) ' ms : Total Wait so far ',len_wait_tot,' seconds.'
           End If

           go to 160  ! Retry finding required lbc_ntimes

        End If  !  if icode /= 0

        ! See if required lbc_ntimes has been read in.
        If (lbc_ntimes >= 1000) Then

           ! First value in file is always >1000. Read next value.
           go to 161

        Elseif (lbc_ntimes <  ms_ntimes) Then

           If (PrintStatus  >=  PrStatus_Diag) Then
              Write (6,*) ' ms : gl_ntimes = ',lbc_ntimes,  &
                   ' read in. gl_ntimes >= ',ms_ntimes,     &
                   ' is required. Read next value.'
           End If
           go to 161

        Elseif (lbc_ntimes >= ms_ntimes) Then

           If (PrintStatus  >=  PrStatus_Diag) Then
              Write (6,*) ' ms : gl_ntimes = ',lbc_ntimes,               &
                   ' read in. gl_ntimes >= ',ms_ntimes,' is required.',  &
                   ' Proceed.'
           End If

           Call Date_and_time (ch_date2, ch_time2)
           If (PrintStatus  >=  PrStatus_Diag) Then
              Write(6,*)  'LBC_COUP: ',                                  &
                   ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),    &
                   ' on ',                                               &
                   ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),    &
                   ' Proceed to call INBOUND in U_MODEL.'
           End If

        End If  !  if lbc_ntimes

     End If !  if mype=0


     ! Call IN_BOUND to update the headers/lookup-table


     !  call inbound for the atmos only
     ! DEPENDS ON: timer
     If (LTIMER) Call TIMER('IN_BOUNDA',3)

     ! DEPENDS ON: inbounda
     Call INBOUNDA(                      &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
          A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,  &
          A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,  &
          A_LEN1_COLDEPC,A_LEN2_COLDEPC )

     ! DEPENDS ON: timer
     If (LTIMER) Call TIMER('IN_BOUNDA',4)

     Call Date_and_time(ch_date2, ch_time2)

     If (PrintStatus  >=  PrStatus_Diag) Then
        Write(6,*)  'LBC_COUP: ',                                      &
             ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ', &
             ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),        &
             ' Proceed to call UPBOUND in U_MODEL.'
     End If

  End If ! (l_lbc_coup .AND.
  !  .NOT.(Num_ALBCs == 2 .AND. ALBC_num == 1))

  If (.Not.l_lbc_coup .And. Num_ALBCs == 2) Then
     If (STEPim(a_im) == ALBC_SwapStep) Then ! Swap to 2nd bndy file

        ! DEPENDS ON: inbounda
        Call INBOUNDA(                      &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
             A_LEN1_LEVDEPC,A_LEN2_LEVDEPC,  &
             A_LEN1_ROWDEPC,A_LEN2_ROWDEPC,  &
             A_LEN1_COLDEPC,A_LEN2_COLDEPC )

     End If
  End If
#endif


  ! DEPENDS ON: timer
  If (LTIMER) Call TIMER('UP_BOUND',3)
  ! DEPENDS ON: up_bound
  Call UP_BOUND(submodel, &
#include "argd1.h"
#include "argduma.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
       ICODE,CMESSAGE)
  ! DEPENDS ON: timer
  If (LTIMER) Call TIMER('UP_BOUND',4)
  ! DEPENDS ON: ereport
  If (ICODE  /=  0) Call Ereport(RoutineName,ICODE,Cmessage)

  Return
End Subroutine lbc_coup_update
#endif
