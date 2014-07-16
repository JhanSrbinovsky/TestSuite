#if defined(CONTROL) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL ----------- SUBROUTINE GEN_INTF -----------------------------------
!LL
!LL Purpose: To generate a PP header and boundary data from a global
!LL          or regional model field at a particular time. Creates an
!LL          interface file for use by a limited area model.
!LL
!LL          GEN_INTF determines whether interface data is required
!LL          for each area on this timestep and calls GEN_INTF_A or
!LL          GEN_INTF_O generate the interface data.
!LL
!LL Control routine for Cray YMP
!LL
!LL Programing standard: UM Documentation paper No. 3,
!LL                      Version No 1, dated 15/01/90
!LL
!LL System components covered: D81
!LL
!LL System task: D81
!LL
!LL Documentation: UM Documentation paper No D8,
!LL
!LLEND ----------------------------------------------------------------

!*L Argument list for GEN_INTF
      SUBROUTINE GEN_INTF (                                             &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "arginfa.h"
#include "argppx.h"
     &                     internal_model,ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "cmaxsizo.h"
#include "csubmodl.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typsts.h"
#include "typptra.h"
#include "typinfa.h"

      INTEGER                                                           &
     &     internal_model,                                              &
                           ! Sub-model indicator
     &       ICODE         ! Return code : =0 Normal exit

      CHARACTER*(80) CMESSAGE   ! Error message if ICODE > 0

!*
!  For call to FINDPTR
#include "c_mdi.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "ctime.h"
#include "ppxlook.h"
#include "chistory.h"
#include "cenvir.h"
#include "cprintst.h"
#include "lbc_coup.h"
      integer lbc_ntimes     ! No of times BC's have been generated.
#if defined(MPP)
      integer ierr           ! Error code
      integer info           ! Return code from GCOM routine.
#endif
      character*8 ch_date2   ! Date returned from date_and_time.
      character*10 ch_time2  ! Time returned from date_and_time.

#include "typocdpt.h"

      INTEGER NFTOUT

      INTEGER      JINTF            ! Interface area index
      INTEGER      LEN_INTF_DATA_DA ! Length of workspace
      INTEGER      IM_IDENT   ! internal model identifier
      INTEGER      IM_INDEX   ! internal model index for STASH arrays
      INTEGER      JSH   ! pointer for ocean surface height field.
      INTEGER      STASHMACRO_TAG ! tag for ocean surface height in D1.

      Integer :: lasize_p, lasize_u, lasize_v
!*L External subroutines called :
      EXTERNAL :: INTF_UNIT
#if defined(ATMOS)
      EXTERNAL :: GEN_INTF_A
#endif
!*
!L Internal structure:

!  Set up internal model identifier and STASH index
      im_ident = internal_model
      im_index = internal_model_index(im_ident)

       ICODE=0
       CMESSAGE=' '

#if defined(ATMOS)
      IF (LLBOUTim(A_IM)) THEN

! Atmosphere Interface

! 1.0 Loop over all areas
        DO JINTF=1,N_INTF_A

!     Determine if interface data required this timestep
          IF (INTERFACE_STEPSim(JINTF,a_im) >  0) THEN

            IF ( MOD(STEPim(a_im)-INTERFACE_FSTEPim(JINTF,a_im),        &
     &      INTERFACE_STEPSim(JINTF,a_im)) == 0                         &
     &      .AND. STEPim(a_im) >= INTERFACE_FSTEPim(JINTF,a_im)         &
     &      .AND. STEPim(a_im) <= INTERFACE_LSTEPim(JINTF,a_im) ) THEN

               write (6,*) 'Gen_Intf: Timestep ',STEPim(a_im),          &
     &         ' : Generate Atmos LBCs for Area ',lbc_stream_a(jintf)

! DEPENDS ON: intf_unit
               call intf_unit (1,jintf,nftout)

               If (STEPim(a_im) == 0 .and.                              &
     &             model_domain  ==  mt_LAM) Then

!                U and V are not set in the halos. Fill with sensible
!                data. Required to un-rotate U & V in Gen_Intf if
!                generating LBCs from a LAM run.

! DEPENDS ON: fill_external_halos
                 Call Fill_External_Halos (d1(ju(1)),                   &
     &                row_length, rows, model_levels, offx, offy)
! DEPENDS ON: fill_external_halos
                 Call Fill_External_Halos (d1(jv(1)),                   &
     &                row_length, n_rows, model_levels, offx, offy)
! DEPENDS ON: fill_external_halos
                 Call Fill_External_Halos (d1(ju_adv(1)),               &
     &                row_length, rows, model_levels, halo_i, halo_j)
! DEPENDS ON: fill_external_halos
                 Call Fill_External_Halos (d1(jv_adv(1)),               &
     &                row_length, n_rows, model_levels, halo_i, halo_j)

               End If

               If (STEPim(a_im) > 0 .and.                               &
     &             model_domain == mt_LAM) Then

! DEPENDS ON: swap_bounds
                 Call Swap_Bounds (d1(ju(1)),                           &
     &                row_length, rows, model_levels, offx, offy,       &
     &                fld_type_u,.true.)

! DEPENDS ON: swap_bounds
                 Call Swap_Bounds (d1(jv(1)),                           &
     &                row_length, n_rows, model_levels, offx, offy,     &
     &                fld_type_v,.true.)

! DEPENDS ON: swap_bounds
                 Call Swap_Bounds (d1(ju_adv(1)),                       &
     &                row_length, rows, model_levels, halo_i, halo_j,   &
     &                fld_type_u,.true.)

! DEPENDS ON: swap_bounds
                 Call Swap_Bounds (d1(jv_adv(1)),                       &
     &                row_length, n_rows, model_levels, halo_i, halo_j, &
     &                fld_type_v,.true.)

               End If

!     Call GEN_INTF_A to generate Atmosphere LBCs for this area.

               If ( lbc_nd(jintf) == 1 ) then       !  LBCs for ND

! DEPENDS ON: gen_intf_a
               call gen_intf_a (                                        &
#include "argduma.h"
#include "arginfa.h"
#include "argptra.h"
#include "argsts.h"
#include "argppx.h"
     &              jintf                                               &
     &,             nftout                                              &
     &,             d1                                                  &
     &,             atmos_im                                            &
     & )

               Elseif ( lbc_nd(jintf) == 0 ) then   !  LBCs for old UM

!     Length of workspace to be dynamic allocated
                LEN_INTF_DATA_DA = LEN_INTFA_DATA(JINTF)

!     Generate Atmos LBCs for pre-NDs UM (ie. pre 5.0)
! DEPENDS ON: gen_intf_a_old_lbcs
                CALL GEN_INTF_A_OLD_LBCS (                              &
#include "argduma.h"
#include "argptra.h"
#include "arginfa.h"
#include "argsts.h"
#include "argppx.h"
#if defined(MPP)
     &  glsize(1,1)*glsize(2,1),                                        &
#endif
     &  JINTF,NFTOUT,                                                   &
     &  D1(jpstar),                                                     &
     &  D1(ju(1)),                                                      &
     &  D1(jv(1)),                                                      &
     &  D1(jtheta(1)),                                                  &
     &  D1(jq(1)),                                                      &
     &  D1(jqcf(1)),                                                    &
     &  D1(jqcl(1)),                                                    &
     &  D1(jexner_rho_levels(1)),                                       &
     &  D1(jexner_theta_levels(1)),                                     &
     &  D1(jtracer(1,1)),                                               &
     &  len_intfa_p(jintf),len_intfa_u(jintf),len_intfa_data(jintf),    &
     &  intf_p_levels(jintf),                                           &
     &  atmos_im,                                                       &
     &  ICODE,CMESSAGE)

                 End If

                IF (ICODE /= 0) THEN
                  GO TO 9999   !  Return
                ENDIF

          if (l_lbc_coup .and.                                          &
     &        lbc_stream_a(jintf) == um_lbc_stream) then

!           Flush buffer to send latest BC's to file.

#if defined(MPP)
            if (mype == 0) then
#endif

!             Flush out all boundary data from buffer.
              call flush_buffer(nftout,icode)

              if (icode /= 0) then
                write (6,*) 'Return Code from FLUSH_BUFFER ',icode,     &
     &          ' for unit number ',nftout
                icode = 501
                write (cmessage,*) 'GENINTF : Error flushing out '//    &
     &          'Boundary Data.'
              endif
#if defined(MPP)

            endif  !  if mype=0

!           Broadcast ICODE to all PE's
            ierr=icode
            call gc_ibcast(450,1,0,nproc,info,ierr)
            icode = ierr
#endif

!           Check ICODE before proceeding.
            if (icode /= 0) then
              write (6,*) ' GENINTF - Error detected'
              write (6,*) ' CMESSAGE : ',CMESSAGE
              write (6,*) ' ICODE : ',ICODE
              go to 9999  !  Return
            endif

!           Get the number of times BC's have been generated.
            lbc_ntimes = ft_lastfield(nftout)

            write (6,*) ' gl : after gen_intf - lbc_ntimes ',           &
     &      lbc_ntimes

#if defined(MPP)
            if (mype == 0) then
#endif

!             Send message to communication file that next lot of
!             BC's have been generated.
              write (190,*) lbc_ntimes

!             Flush message out.
! DEPENDS ON: um_fort_flush
              call um_fort_flush (190,icode)
              if (icode /= 0) then
                write (6,*) 'Return Code from FLUSH ',icode
                icode = 503
                write (cmessage,*) 'GENINTF : Error flushing out '//    &
     &          'contents for Unit 190.'
                go to 150
              endif

!             Write a text message that next lot of BC's have
!             been generated.
              call date_and_time(ch_date2, ch_time2)
              if (lbc_fc_hrs >= 0) then
                write (191,*)                                           &
     &          ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),      &
     &          ' Boundary data has been generated for T+',lbc_fc_hrs
              else
                write (191,*)                                           &
     &          ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),      &
     &          ' Boundary data has been generated for T',lbc_fc_hrs
              endif

!             Flush message out.
! DEPENDS ON: um_fort_flush
              call um_fort_flush (191,icode)
              if (icode /= 0) then
                write (6,*) 'Return Code from FLUSH ',icode
                icode = 504
                write (cmessage,*) 'GENINTF : Error flushing out '//    &
     &          'contents for Unit 191.'
                go to 150
              endif

#if defined(MPP)
            endif  !  if mype=0
#endif

            call date_and_time(ch_date2, ch_time2)

            write(6,*) 'LBC_COUP: ',                                    &
     &      ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',   &
     &      ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),          &
     &      ' Timestep ',stepim(a_im),                                  &
     &      ' Boundary conditions generated.'

!           If all boundary conditions have been generated, add
!           value 7777 to end of file.

            if (stepim(a_im)  ==  interface_LSTEPim(jintf,a_im)) then

              lbc_ntimes = 7777

              write (6,*) ' gl : after gen_intf - lbc_ntimes ',         &
     &        lbc_ntimes

#if defined(MPP)
              if (mype == 0) then
#endif

!               Write to communication file and flush.
                write (190,*) lbc_ntimes
! DEPENDS ON: um_fort_flush
                call um_fort_flush (190,icode)
                if (icode /= 0) then
                  write (6,*) 'Return Code from FLUSH ',icode
                  icode = 506
                  write (cmessage,*) 'GENINTF : Error flushing out '//  &
     &            'contents for Unit 190.'
                  go to 150
                endif

!               Write text message and flush.
                write (191,*) ' All Boundary data has been generated.'
! DEPENDS ON: um_fort_flush
                call um_fort_flush (191,icode)
                if (icode /= 0) then
                  write (6,*) 'Return Code from FLUSH ',icode
                  icode = 507
                  write (cmessage,*) 'GENINTF : Error flushing out '//  &
     &            'contents for Unit 191.'
                  go to 150
                endif

#if defined(MPP)
              endif  !  if mype=0
#endif

              write (6,*)                                               &
     &        'LBC_COUP: GEN_INTF - All Boundary Conditions generated ',&
     &        'for stream ',jintf

            endif   !  if stepim(a_im)

 150        continue

#if defined(MPP)
!           Broadcast ICODE to all PEs.
            ierr=icode
            call gc_ibcast(450,1,0,nproc,info,ierr)
            icode = ierr
#endif

!           Check ICODE before proceeding.
            if (icode /= 0) then
              write (6,*) ' GENINTF - Error detected'
              write (6,*) ' ICODE : ',ICODE
              write (6,*) ' CMESSAGE : ',CMESSAGE
              go to 9999  !  Return
            endif

          endif  !  if l_lbc_coup

!         Close file if no more LBCs
!         Reinitialised files are closed in Gen_Intf
          if (stepim(a_im) == interface_LSTEPim(jintf,a_im) .and.       &
     &        ft_steps(nftout) == 0) then

! DEPENDS ON: file_close
            call file_close (nftout,ft_environ(nftout),                 &
     &                       len_ft_envir(nftout),1,0,icode)

            If (PrintStatus >= PrStatus_Normal) Then
              write(6,*) 'Gen_Intf: Timestep ',stepim(a_im),            &
     &        ' : LBC file closed for Area ',lbc_stream_a(jintf)
            End If

          End If

          ENDIF  !  if (mod(step...
        ENDIF    !  if (interface_steps...
      ENDDO      !  Loop over JINTF

      END IF ! LLBOUTim(A_IM)
#endif

 9999 RETURN
      END SUBROUTINE GEN_INTF

#endif
