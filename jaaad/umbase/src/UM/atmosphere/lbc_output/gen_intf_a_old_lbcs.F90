#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL ----------- SUBROUTINE GEN_INTF_A_OLD_LBCS -------------------------
!LL
!LL Purpose: To generate a PP header and boundary data from a global
!LL          or regional model field at a particular time. Creates an
!LL          interface file for use by a limited area model.
!LL
!LL Control routine for Cray YMP
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!
!     History prior to 4.5 removed.
!     See 5.1 version for full history.
!
!LL   4.5  15/04/98   Added start-end arguments to V_INT routines
!LL                   S.D.Mullerworth
!    4.5   17/10/97 Parallelise the horizontal and vertical
!                   interpolations, and use a multi-level
!                   gather.
!
!                   *** T3E Specific Code ***
!
!                     Authors: Bob Carruthers, Cray Research
!                              Paul Burton
!LL  4.5  29/07/98  Rename CINTF to CINTFA. Call INTF_UNIT.
!LL                 Remove DEF,RECONF. D. Robinson.
!    4.5   17/08/98 Global/Mes parallel running. Send messages to
!                   communication & information files. D. Robinson.
!LL  5.0   20/05/99 Change variable names
!LL                 Remove DA arguments from GEN_INTF_A
!LL                                                        P.Burton
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.2  13/11/00  GEN_INTF deleted. GEN_INTF_A renamed to
!LL                 GEN_INTF_A_OLD_UM. D.Robinson
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!    6.1  18/08/04  Adapt GENINTF1 to produce old LBCs from NDs
!                   prognostics. D Robinson
!    6.2  21/10/05  Replace GSYNC with SSYNC. P.Selwood.
!
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


!*L Argument list for GEN_INTF_A_OLD_UM
! =================================================================
! Old GEN_INTF_A renamed to GEN_INTF_A_OLD_LBCS
! Retain until decision made on generating LBCs for UM 4.5
! New GEN_INTF_A (in deck GENINTA1) generates LBCs for New Dynamics
! =================================================================
      SUBROUTINE GEN_INTF_A_OLD_LBCS (                                  &
#include "argduma.h"
#include "argptra.h"
#include "arginfa.h"
#include "argsts.h"
#include "argppx.h"
     &  WORK_FLD_SIZE,                                                  &
     &  JINTF,NFTOUT,                                                   &
     &  PSTAR,U,V,THETA,Q,QCF,QCL,                                      &
     &  EXNER_P_RHO,EXNER_P_THETA,TRACER,                               &
     &  LEN_INTF_P,LEN_INTF_U,LEN_INTF_DATA,                            &
     &  INTF_P_LEVS,                                                    &
     &  INTERNAL_MODEL,                                                 &
     &  ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typptra.h"
#include "typinfa.h"
#include "csubmodl.h"
#include "typsts.h"

      INTEGER                                                           &
     &  JINTF                                                           &
                        !  Index to interface area
     & ,NFTOUT                                                          &
                        !  Unit number for interface data
     & ,LEN_INTF_P                                                      &
                        !  Length of interface p* grid
     & ,LEN_INTF_U                                                      &
                        !  Length of interface u  grid
     & ,LEN_INTF_DATA                                                   &
                        !  Length of interface data
     & ,INTF_P_LEVS                                                     &
                        !  No of model levels in interface data
     & ,WORK_FLD_SIZE                                                   &
                        ! Size of full global field
     & ,INTERNAL_MODEL

      INTEGER                                                           &
     &       ICODE         ! Return code : =0 Normal exit
!                          !               >0 Error condition
      CHARACTER*(80) CMESSAGE ! Error message if ICODE>0

      REAL                                                              &
     &  PSTAR(theta_field_size)                                         &
                                             !  Model P* data
     & ,U(u_off_size,MODEL_LEVELS)                                      &
                                             !  Model u components
     & ,V(v_off_size,MODEL_LEVELS)                                      &
                                             !  Model v components
     & ,THETA(theta_field_size,MODEL_LEVELS)                            &
                                             !  Model theta data
     & ,Q(theta_halo_size,WET_LEVELS)                                   &
                                             !  Model Q data
     & ,QCF(theta_halo_size,WET_LEVELS)                                 &
                                             !  Model QCF data
     & ,QCL(theta_halo_size,WET_LEVELS)                                 &
                                             !  Model QCL data
     & ,EXNER_P_RHO(theta_off_size,MODEL_LEVELS)                        &
                                                 !
     & ,EXNER_P_THETA(theta_off_size,MODEL_LEVELS)                      &
     & ,TRACER(TR_VARS*theta_off_size,TR_LEVELS)
                                             !  Model tracer data

      LOGICAL                                                           &
     &  LPACK_32B                                                       &
                                       !  Packing Indicator
     & ,LPACK_PPXREF
!*
      INTEGER                                                           &
     &       I,                                                         &
     &       J,                                                         &
     &       IADDR,IADDR_V,                                             &
     &       LEVEL,                                                     &
     &       VAR,                                                       &
     &       LOOKUP_START,                                              &
     &       LEN_IO,                                                    &
     &       SEC,                                                       &
     &       DATA_START,                                                &
     &       CODE                                                       &
     &       ,NTIME                                                     &
                            ! postion number of interface data
     &,LEN_PPNAME                                                       &
     &      ,im_ident                                                   &
                            !  Internal model identifier
     &      ,im_index       !  Internal model index in STASH arrays

      INTEGER YY,MM,DD,HR,MN,SS,DAY_NO

      REAL TEMP

      LOGICAL                                                           &
     & ROT_IN   ! =T, if input model grid rotated
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "clookadd.h"
#include "ctracera.h"
#include "ppxlook.h"
#include "c_r_cp.h"
#include "c_mdi.h"
#include "gccom.h"
      INTEGER info
!
      integer len_intf_uv_data            ! Length of U or V data to be
                                          ! collected after processing
                                          ! or copying
      integer len_intf_theta_data         ! Length of THETA data to be
                                          ! collected after processing
                                          ! or copying
      integer len_intf_qt_data            ! Length of QT data to be
                                          ! collected after processing
                                          ! or copying
      integer lbc_address_work(maxproc)   ! Address of the work area
                                          ! on remote PE's
      integer lbc_address_data(maxproc)   ! Address of the data area
                                          ! on remote PE's
      integer iaddr_good                  ! The real address in
                                          ! intf_data, maintained
                                          ! correctly across PE's
      integer iaddr_u                     ! Address of U in 'intf_data'
      integer iaddr_theta                 ! Address of THETA in
                                          ! 'intf_data'
      integer iaddr_qt                    ! Address of QT in 'intf_data'
!
      integer pe_for_level_uv(MODEL_LEVELS)
                                    ! PE to work on a given level
                                          ! for U/V
      integer pe_for_level_theta(MODEL_LEVELS)
                                    ! PE to work on a given level
                                          ! for THETA
      integer pe_for_level_qt(MODEL_LEVELS)
                                    ! PE to work on a given level
                                          ! for QT
      integer local_level(MODEL_LEVELS)
                                    ! Index of a global level on a
                                          ! local PE when the levels are
                                          ! distributed over several PEs
      integer pe_for_var(intf_lookupsa)   ! PE numbers for each variable
      real t1, t2, t3                     ! local Timers
!
!
#if defined(TIME_LBC)
#include "t3eclktk.h"
#endif

!
      INTEGER EXPPXI
      EXTERNAL EXPPXI
      INTEGER N1,N2,N3,N4,N5,NPACK
      INTEGER LEN_DATA
#include "cntl_io.h"
!
      integer                                                           &
     & disk_address                                                     &
                            ! Current rounded disk address
     &,disk_length                                                      &
                            ! Current data record length on disk
     &,len_buf              ! Maximum Record length to O/P
!
!*L Workspace used

!     Dynamic allocated workspace
      REAL                                                              &
     & intf_data(((len_intf_data+um_sector_size-1)/                     &
     & um_sector_size)*um_sector_size),                                 &
     & INTF_WORK(LEN_INTF_P,2),                                         &
     & WORK_U     (LEN_INTF_P,MODEL_LEVELS),                            &
     & WORK_V     (LEN_INTF_P,MODEL_LEVELS),                            &
     & WORK_THETAL(LEN_INTF_P,MODEL_LEVELS),                            &
     & WORK_QT    (LEN_INTF_P,WET_LEVELS),                              &
     & WORK1      (LEN_INTF_P),                                         &
     & WORK2      (LEN_INTF_P),                                         &
     & INTF_PSTAR(LEN_INTF_P),                                          &
     & P_OUT(LEN_INTF_P),                                               &
     & P_HALF_TMP(LEN_INTF_P,INTF_P_LEVS+1),                            &
     & P_EXNER_HALF_TMP(LEN_INTF_P,INTF_P_LEVS+1),                      &
     & U_TEMP(u_off_size, model_levels),                                &
     & V_TEMP(v_off_size, model_levels),                                &
     & WORK_GLOBAL1(WORK_FLD_SIZE),                                     &
     & WORK_GLOBAL2(WORK_FLD_SIZE),                                     &
     &       A_IO
!     Workspace for T3E vector function rtor_v
      integer n_input                  ! No. of inputs for rtor_v
      REAL    P_HALF_TMP_wk(LEN_INTF_P,INTF_P_LEVS+1)
      REAL    KAPPA_HALF_wk(LEN_INTF_P,INTF_P_LEVS+1)
      REAL    P_TMP_wk(LEN_INTF_P,MODEL_LEVELS)
      REAL    KAPPA_wk(LEN_INTF_P,MODEL_LEVELS)
!dir$ cache_align intf_data

!*---------------------------------------------------------------------
      INTEGER IP_P,IP_U
      CHARACTER*80 STRING         ! work array
      CHARACTER*14 PPNAME         ! boundary output filename
      INTEGER START_ADDR          ! Start Address in LOOKUP table

!*---------------------------------------------------------------------
!     Stash item numbers for interface fields
!     Any change to code generating and testing ITEM_INTFA should also
!     consider the corresponding use of ITEM_BOUNDA in INBOUND1/CHKLKBA1
      INTEGER ITEM_INTFA (INTF_LOOKUPSA)
!*---------------------------------------------------------------------
!*
      Integer            :: k
      Integer            :: pstar_fld_type
      Integer            :: pstar_halo_type
      Integer            :: u_halo_type, v_halo_type
      Integer            :: theta_halo_type, q_halo_type
      Integer, Parameter :: lbc_len_inthd = 15

      Real, allocatable :: ThetaL (:,:)
      Real, allocatable :: QT (:,:)
      Real, allocatable :: p_lbc (:,:)
      Real, allocatable :: p_src (:,:)

      Integer :: ErrorStatus  !  Return Code
      Integer :: gather_pe

      Integer, Parameter :: pe0 = 0
!     Integer, Parameter :: pe1 = 1
!     Integer, Parameter :: pe2 = 2
!     Integer, Parameter :: pe3 = 3

      Integer :: level_count
      Integer :: istat

      Character (Len=*), Parameter :: RoutineName = 'Gen_Intf_old_lbcs'

#include "p_exnerc.h"

!L Internal structure:

       ICODE=0
       CMESSAGE=' '

       im_ident = internal_model
       im_index = internal_model_index(im_ident)
       LPACK_32B = INTF_PACK(JINTF) == 1
       LPACK_PPXREF = INTF_PACK(JINTF) == 2

!L     Atmosphere interface

! Logical to indicate if model grid is rotated
      ROT_IN=A_REALHD(5) /= 90..OR.A_REALHD(6) /= 0.

        item_intfa(:)=0

! Set up list of variables expected to be boundary updated.
          ITEM_INTFA(1) = 1   ! Pstar
          ITEM_INTFA(2) = 2   ! u-compt wind
          ITEM_INTFA(3) = 3   ! v-compt wind
          ITEM_INTFA(4) = 5   ! thetal
          ITEM_INTFA(5) = 11  ! qt
          IF (TR_VARS  >   0) THEN
! Find STASH item no. for each tracer in use.
            I=0          ! count tracers in use
            DO J = A_TRACER_FIRST,A_TRACER_LAST
              IF (SI(J,0,im_index) /= 1) THEN  ! tracer is in use
                I = I+1
                ITEM_INTFA(5+I) = J
              END IF
            END DO
! Number of tracers found should correspond to TR_VARS
            IF (I /= TR_VARS) THEN
              WRITE(6,*)' GEN_INTF_A: no.of tracers found, ',I,         &
     & ', differs from TR_VARS, ',TR_VARS
            CMESSAGE=' GEN_INTF_A: inconsistency in number of tracers'
              ICODE = 100
              GO TO 999
            END IF
          END IF
!       QCF not catered for ; not used in 'old' LBC files.


!L 1.0 Generate data on the boundary zone of limited area grid

      IADDR=1

!     Set up pointers to interpolation coefficients for this area
      IP_P = 1
      IP_U = 1
      IF (JINTF >  1) THEN
        DO J=1,JINTF-1
          If (LBC_ND(J) == 0) Then
          IP_P = IP_P + LEN_INTFA_P(J)
          IP_U = IP_U + LEN_INTFA_U(J)
          End If
        ENDDO
      ENDIF
!
!--check if there are enough processors to use a parallel algorithm
!  (remember that PE's 0 and 1 hold only U/V data, PE 2 holds
!   only THETA data, and PE 3 holds only QT data, etc, so that
!   the vertical interpolations can be done in parallel after
!   the horizontal gathers and interpolations)

!
!--set up the pe data for the multi-level gathers and horizontal
!  interpolations
      do i=1, MODEL_LEVELS
        pe_for_level_uv(i)    = mod(i-1,nproc)
        pe_for_level_theta(i) = mod(i-1,nproc)
        pe_for_level_qt(i)    = mod(i-1,nproc)
        local_level(i)=1
      end do
!
!--initialise the 'pe_for_var' array (Vertical Interpolation)
      do i=1, intf_lookupsa
        pe_for_var(i) = mod(i-1,nproc)
      end do

!L 1.0 Generate data on the boundary zone of limited area grid

      iaddr_good=1
      iaddr=iaddr_good

!L 1.1 P Star

#if defined(TIME_LBC)
      t1=rtc()
      t2=0.
      t3=rtc()
#endif

      pstar_fld_type  = 1
      pstar_halo_type = 3
      gather_pe       = 0

! DEPENDS ON: gather_field
      Call Gather_Field (PStar, Work_Global1,                           &
     &                   lasize(1,pstar_fld_type,pstar_halo_type),      &
     &                   lasize(2,pstar_fld_type,pstar_halo_type),      &
     &                   glsize(1,pstar_fld_type),                      &
     &                   glsize(2,pstar_fld_type),                      &
     &                   pstar_fld_type,                                &
     &                   pstar_halo_type,                               &
     &                   gather_pe,                                     &
     &                   gc_all_proc_group,                             &
     &                   info, cmessage )

#if defined(TIME_LBC)
      if(mype == 0)                                                     &
     & write(0,*)'Time for Gather on PSTAR was ',                       &
     & (rtc()-t3)/ticks_per_second
#endif

#if defined(TIME_LBC)
      t2=t2-rtc()
      t3=rtc()
#endif
      IF (mype  ==  0) THEN
!L 1.1.1 Horizontal interpolation
! DEPENDS ON: h_int_bl
        CALL H_INT_BL(glsize(2,pstar_fld_type),                         &
     &                glsize(1,pstar_fld_type),                         &
     &                LEN_INTF_P,                                       &
     &                AP_INDEX_B_L(IP_P),AP_INDEX_B_R(IP_P),            &
     &                WORK_GLOBAL1,                                     &
     &                AP_WEIGHT_B_L(IP_P),AP_WEIGHT_B_R(IP_P),          &
     &                AP_WEIGHT_T_L(IP_P),AP_WEIGHT_T_R(IP_P),          &
     &                INTF_DATA(IADDR) )

!L 1.1.2 Save pstar for vertical interpolation
        IF (INTF_VERT_INTERP(JINTF)) THEN

          DO I=1,LEN_INTF_P
            INTF_PSTAR(I)=INTF_DATA(IADDR+I-1)
          ENDDO
        ENDIF

      ENDIF  ! IF (mype  ==  0)
#if defined(TIME_LBC)
      if(mype == 0)                                                     &
     & write(0,*)'Time for Horizontal on PSTAR was ',                   &
     & (rtc()-t3)/ticks_per_second
      t2=t2+rtc()
#endif
!
!--broadcast this result to everyone
      call gc_rbcast(8491, len_intf_p, 0, nproc, info, intf_pstar)
!--update the address in 'intf_data'
      iaddr_good=iaddr_good+len_intf_p

      iaddr=iaddr_good

#if defined(TIME_LBC)
      t3=rtc()
#endif
!--loop over levels
      do level=1, MODEL_LEVELS

!L 1.2 U and V components

!L 1.2.1 Rotate winds to standard lat-lon if input grid is rotated.

        IF (ROT_IN) THEN

! DEPENDS ON: w_eqtoll
          call w_eqtoll(coeff3, coeff4, u(1, level), v(1, level),       &
     &                  u_temp(1, level),v_temp(1, level),              &
     &                  U_FIELD_SIZE, U_FIELD_SIZE)

        ELSE

          Do i=1,u_off_size
            u_temp(i,level) = u(i,level)
          End Do
          Do i=1,v_off_size
            v_temp(i,level) = v(i,level)
          End Do

        ENDIF

      end do
!L 1.2.2 Horizontal interpolation - winds

!     Gather u and v onto PEs for horizontal interpolation

      Do level = 1, model_levels

        u_halo_type = 1

! DEPENDS ON: gather_field
        call gather_field (u_temp(1,level),                             &
     &                     work_global1,                                &
     &                     lasize(1,fld_type_u,u_halo_type),            &
     &                     lasize(2,fld_type_u,u_halo_type),            &
     &                     glsize(1,fld_type_u),                        &
     &                     glsize(2,fld_type_u),                        &
     &                     fld_type_u,                                  &
     &                     u_halo_type,                                 &
     &                     pe_for_level_uv(level),                      &
     &                     gc_all_proc_group,                           &
     &                     icode,                                       &
     &                     cmessage )

        v_halo_type = 1

! DEPENDS ON: gather_field
        call gather_field (v_temp(1,level),                             &
     &                     work_global2,                                &
     &                     lasize(1,fld_type_v,v_halo_type),            &
     &                     lasize(2,fld_type_v,v_halo_type),            &
     &                     glsize(1,fld_type_v),                        &
     &                     glsize(2,fld_type_v),                        &
     &                     fld_type_v,                                  &
     &                     v_halo_type,                                 &
     &                     pe_for_level_uv(level),                      &
     &                     gc_all_proc_group,                           &
     &                     icode,                                       &
     &                     cmessage )


#if defined(TIME_LBC)
      if(mype == pe_for_var(1))                                         &
     & write(0,*)'Time for Gather on U was ',                           &
     & (rtc()-t3)/ticks_per_second
      if(mype == pe_for_var(2))                                         &
     & write(0,*)'Time for Gather on V was ',                           &
     & (rtc()-t3)/ticks_per_second
#endif

#if defined(TIME_LBC)
      t2=t2-rtc()
      t3=rtc()
#endif

        If (mype == pe_for_level_uv(level) ) Then

! DEPENDS ON: h_int_bl
          CALL H_INT_BL(glsize(2,fld_type_u),                           &
     &                  glsize(1,fld_type_u),                           &
     &                  LEN_INTF_U,                                     &
     &                  AU_INDEX_B_L(IP_U),AU_INDEX_B_R(IP_U),          &
     &                  WORK_GLOBAL1,                                   &
     &                  AU_WEIGHT_B_L(IP_U),AU_WEIGHT_B_R(IP_U),        &
     &                  AU_WEIGHT_T_L(IP_U),AU_WEIGHT_T_R(IP_U),        &
     &                  WORK1)

! DEPENDS ON: h_int_bl
          CALL H_INT_BL(glsize(2,fld_type_v),                           &
     &                  glsize(1,fld_type_v),                           &
     &                  LEN_INTF_U,                                     &
     &                  AV_INDEX_B_L(IP_U),AV_INDEX_B_R(IP_U),          &
     &                  WORK_GLOBAL2,                                   &
     &                  AV_WEIGHT_B_L(IP_U),AV_WEIGHT_B_R(IP_U),        &
     &                  AV_WEIGHT_T_L(IP_U),AV_WEIGHT_T_R(IP_U),        &
     &                  WORK2)

!         work1/work2 contain u/v after horizontal interpolation

! DEPENDS ON: w_lltoeq
          CALL W_LLTOEQ(COEFF1(IP_U),COEFF2(IP_U),                      &
     &                  WORK1,                                          &
                                                   ! u_in
     &                  WORK2,                                          &
                                                   ! v_in
     &                  WORK_U(1,LEVEL),                                &
                                                   ! u_out
     &                  WORK_V(1,LEVEL),                                &
                                                   ! v_out
     &                  LEN_INTF_U,LEN_INTF_U)

!         work_u/work_v contain rotated u/v

        End If

      End Do

!--make sure everyone has finished computing their data
      Call GC_SSync (nproc, istat)

!  Collect U data on PE Zero and V data on PE 1

      Do level = 1, model_levels

!--U data

        If (mype == pe_for_level_uv(level) ) Then
          If (pe_for_level_uv(level) /= pe_for_var(1) ) Then

!           Send a copy to PE 0

            Call gc_rsend ( 91, len_intf_p,                             &
     &                      pe_for_var(1), info,                        &
     &                      work_u(1,level), work_u(1,level) )

          End If
        End If

        If (mype == pe_for_var(1) ) Then
          If (pe_for_level_uv(level) /= pe_for_var(1) ) Then

            ! Receive copy on PE 0

            Call gc_rrecv ( 91, len_intf_p,                             &
     &                      pe_for_level_uv(level), info,               &
     &                      work_u(1,level), work_u(1,level) )

          End If
        End If

!--V data

       If (mype == pe_for_level_uv(level) ) Then
         If ( pe_for_level_uv(level) /= pe_for_var(2) ) Then

!           Send a copy to PE 1

            Call gc_rsend ( 92, len_intf_p,                             &
     &                      pe_for_var(2), info,                        &
     &                      work_v(1,level), work_v(1,level) )

          End If
        End If

        If (mype == pe_for_var(2) ) Then
          If ( pe_for_level_uv(level) /= pe_for_var(2) ) Then

!           Receive copy on PE 1

            Call gc_rrecv ( 92, len_intf_p,                             &
     &                      pe_for_level_uv(level), info,               &
     &                      work_v(1,level), work_v(1,level) )

          End If
        End If

      End Do

      Call GC_SSync (nproc, istat)

#if defined(TIME_LBC)
      t2=t2+rtc()
      if(mype == pe_for_var(1))                                         &
     & write(0,*)'Time for Horizontal on U was ',                       &
     & (rtc()-t3)/ticks_per_second
      if(mype == pe_for_var(2))                                         &
     & write(0,*)'Time for Horizontal on V was ',                       &
     & (rtc()-t3)/ticks_per_second
#endif

!--compute transfer addresses and lengths
      if (intf_vert_interp(jintf)) then
        len_intf_uv_data=intf_p_levels(jintf)*len_intf_u
      else
        len_intf_uv_data=MODEL_LEVELS*len_intf_u
      endif
!--set up the addresses
      iaddr_u=iaddr_good
      iaddr_v=iaddr_good+len_intf_uv_data
!--now update 'iaddr_good'
      iaddr_good=iaddr_good+len_intf_uv_data*2

      iaddr=iaddr_good

!L 1.3 THETAL

! Derive ThetaL from Model Variables.

      Allocate ( ThetaL (theta_off_size, Model_Levels) )

! DEPENDS ON: theta_to_thetal
      Call Theta_to_ThetaL (Theta, QCF, QCL, Exner_P_Theta, ThetaL,     &
     &     row_length, rows, Offx, Offy, halo_i, halo_j,                &
     &     wet_levels, model_levels)

!L 1.3.1 Horizontal interpolation - thetal

#if defined(TIME_LBC)
      t3=rtc()
#endif


      theta_halo_type = 1

      Do level = 1, model_levels

! DEPENDS ON: gather_field
        Call gather_field (ThetaL(1,level),                             &
     &                     work_global1,                                &
     &                     lasize(1,fld_type_p,theta_halo_type),        &
     &                     lasize(2,fld_type_p,theta_halo_type),        &
     &                     glsize(1,fld_type_p),                        &
     &                     glsize(2,fld_type_p),                        &
     &                     fld_type_p,                                  &
     &                     theta_halo_type,                             &
     &                     pe_for_level_theta(level),                   &
     &                     gc_all_proc_group,                           &
     &                     icode,                                       &
     &                     cmessage )


#if defined(TIME_LBC)
      if(mype == pe_for_var(3))                                         &
     & write(0,*)'Time for Gather on THETA was ',                       &
     & (rtc()-t3)/ticks_per_second
#endif

#if defined(TIME_LBC)
      t2=t2-rtc()
      t3=rtc()
#endif
!--now do the horizontal interpolation in parallel

        If (mype == pe_for_level_theta(level) ) Then

! DEPENDS ON: h_int_bl
          Call H_INT_BL ( glsize(2,fld_type_p),                         &
     &                    glsize(1,fld_type_p),                         &
     &                    len_intf_p,                                   &
     &                    AP_INDEX_B_L(IP_P),AP_INDEX_B_R(IP_P),        &
     &                    WORK_GLOBAL1,                                 &
     &                    AP_WEIGHT_B_L(IP_P),AP_WEIGHT_B_R(IP_P),      &
     &                    AP_WEIGHT_T_L(IP_P),AP_WEIGHT_T_R(IP_P),      &
     &                    WORK_THETAL(1,LEVEL) )

        End If

      End Do

!--make sure everyone has finished computing their data
      Call GC_SSync (nproc, istat)

      Deallocate (ThetaL)

!  Collect ThetaL data on PE 2

      Do level = 1, model_levels

        If (mype == pe_for_level_theta(level) ) Then
          If (pe_for_level_theta(level) /= pe_for_var(3) ) Then

!           Send a copy to PE 2

            Call gc_rsend ( 93, len_intf_p,                             &
     &                      pe_for_var(3), info,                        &
     &                      work_thetal(1,level), work_thetal(1,level) )

          End If
        End If

        If (mype == pe_for_var(3) ) Then
          If (pe_for_level_theta(level) /= pe_for_var(3) ) Then

!           Receive copy on PE 2

            Call gc_rrecv ( 93, len_intf_p,                             &
     &                      pe_for_level_theta(level), info,            &
     &                      work_thetal(1,level), work_thetal(1,level) )

          End If
        End If

      End Do

      Call GC_SSync (nproc, istat)

#if defined(TIME_LBC)
      t2=t2+rtc()
      if(mype == pe_for_var(3))                                         &
     & write(0,*)'Time for Horizontal on THETA was ',                   &
     & (rtc()-t3)/ticks_per_second
#endif

!--compute transfer addresses and lengths
      if (intf_vert_interp(jintf)) then
        len_intf_theta_data=intf_p_levels(jintf)*len_intf_p
      else
        len_intf_theta_data=MODEL_LEVELS*len_intf_p
      endif
      iaddr_theta=iaddr_good
!--now update 'iaddr_good'
      iaddr_good=iaddr_good+len_intf_theta_data

      iaddr=iaddr_good

!L 1.4 QT

! Set up workspace for QT

       Allocate ( QT ( theta_halo_size, Wet_Levels) )

! Derive QT from Q, QCF and QCL.

! DEPENDS ON: q_to_qt
       Call Q_To_QT (Q, QCF, QCL, QT, theta_halo_size, Wet_Levels)

!L 1.4.1 Horizontal interpolation - QT

#if defined(TIME_LBC)
      t3=rtc()
#endif


      q_halo_type = 2

      Do level = 1, wet_levels

! DEPENDS ON: gather_field
        Call gather_field (QT(1,level),                                 &
     &                     work_global1,                                &
     &                     lasize(1,fld_type_p,q_halo_type),            &
     &                     lasize(2,fld_type_p,q_halo_type),            &
     &                     glsize(1,fld_type_p),                        &
     &                     glsize(2,fld_type_p),                        &
     &                     fld_type_p,                                  &
     &                     q_halo_type,                                 &
     &                     pe_for_level_qt(level),                      &
     &                     gc_all_proc_group,                           &
     &                     icode,                                       &
     &                     cmessage )


#if defined(TIME_LBC)
      if(mype == pe_for_var(4))                                         &
     & write(0,*)'Time for Gather on QT was ',                          &
     & (rtc()-t3)/ticks_per_second
#endif

#if defined(TIME_LBC)
      t2=t2-rtc()
      t3=rtc()
#endif
!--now do the horizontal interpolation in parallel

        If (mype == pe_for_level_qt(level) ) Then

! DEPENDS ON: h_int_bl
          CALL H_INT_BL (glsize(2,fld_type_p),                          &
     &                   glsize(1,fld_type_p),                          &
     &                   LEN_INTF_P,                                    &
     &                   AP_INDEX_B_L(IP_P),AP_INDEX_B_R(IP_P),         &
     &                   WORK_GLOBAL1,                                  &
     &                   AP_WEIGHT_B_L(IP_P),AP_WEIGHT_B_R(IP_P),       &
     &                   AP_WEIGHT_T_L(IP_P),AP_WEIGHT_T_R(IP_P),       &
     &                   WORK_QT(1,LEVEL) )

        End If

      End Do

      Call GC_SSync (nproc, istat)

      Deallocate (QT)

!  Collect QT data on PE 3

      Do level = 1, wet_levels

        If ( mype == pe_for_level_qt(level) ) Then
          If ( pe_for_level_qt(level) /= pe_for_var(4) ) Then

!           Send a copy to PE 3

            Call gc_rsend ( 94, len_intf_p,                             &
     &                      pe_for_var(4), info,                        &
     &                      work_qt(1,level), work_qt(1,level) )

          End If
        End If

        If ( mype == pe_for_var(4) ) Then
          If ( pe_for_level_qt(level) /= pe_for_var(4) ) Then

!           Receive copy on PE 3

            Call gc_rrecv ( 94, len_intf_p,                             &
     &                      pe_for_level_qt(level), info,               &
     &                      work_qt(1,level), work_qt(1,level) )

          End If
        End If

      End Do

      Call GC_SSync (nproc, istat)

#if defined(TIME_LBC)
      t2=t2+rtc()
      if(mype == pe_for_var(4))                                         &
     & write(0,*)'Time for Horizontal on QT was ',                      &
     & (rtc()-t3)/ticks_per_second

      if(mype == 0) write(0,*)'Time for Horizontal Interpolation was ', &
     & t2/ticks_per_second
#endif

!--compute transfer addresses and lengths
      if (intf_vert_interp(jintf)) then
        len_intf_qt_data=intf_q_LEVELS(jintf)*len_intf_p
      else
        len_intf_qt_data=WET_LEVELS*len_intf_p
      endif
      iaddr_qt=iaddr_good
!--now update 'iaddr_good'
      iaddr_good=iaddr_good+len_intf_qt_data

      IF (INTF_VERT_INTERP(JINTF)) THEN

! --------------------------------------
! Allocate space for the pressure fields
! --------------------------------------

        allocate ( p_src (len_intf_p,model_levels)          )
        allocate ( p_lbc (len_intf_p,intf_p_levels(jintf))  )

! ----------------------------------
! Calculate pressures for lbc levels
! ----------------------------------

        Do k = 1, intf_p_levels(jintf)
          Do i = 1, len_intf_p
            p_lbc(i,k) =                                                &
     &      intf_ak(k,jintf) + intf_pstar(i) * intf_bk(k,jintf)
          End Do
        End Do

! ----------------------------------------
! Calculate pressures for model rho levels
! ----------------------------------------

        Gather_PE = pe_for_var(1)

! DEPENDS ON: calc_p_src
        Call Calc_P_Src (                                               &
#include "arginfa.h"
     &       ip_p, theta_off_size, len_intf_p,                          &
     &       gather_PE, pe_for_level_uv, local_level,                   &
     &       exner_p_rho, p_src )

! --------------------------------------
! Vert Int for u is done on PE 0
! Vert Int for v is done on PE 1
! P_SRC is set up on PE 0 ; copy to PE 1
! --------------------------------------

        If (mype==pe_for_var(1)) Then ! Send a copy to PE 1
          info=GC_NONE
          Call gc_rsend (95,len_intf_p*model_levels,pe_for_var(2),      &
     &                   info,p_src,p_src)
          End If

        If (mype==pe_for_var(2)) Then ! Receive a copy from PE 0
          info=GC_NONE
          Call gc_rrecv (95,len_intf_p*model_levels,pe_for_var(1),      &
     &                   info,p_src,p_src)
        End If

      End If ! Vert Interp

!L
!L Now do the vertical interpolation in parallel
!L
#if defined(TIME_LBC)
      t2=rtc()
#endif

      iaddr=iaddr_u
      if (mype == pe_for_var(1) .or. mype == pe_for_var(2)) then
#if defined(TIME_LBC)
        t3=rtc()
#endif

!L 1.2.3 Vertical interpolation - winds

        IF (INTF_VERT_INTERP(JINTF)) THEN

          DO LEVEL=1,MODEL_LEVELS
            DO I=LEN_INTF_U+1,LEN_INTF_P
              INTF_WORK(I,1) = 0.0
              INTF_WORK(I,2) = 0.0
            ENDDO
          ENDDO

          DO LEVEL=1,INTF_P_LEVELS(JINTF)

            If (mype == pe_for_var(1) ) Then     !  u
! DEPENDS ON: v_int
              Call V_INT (p_src,p_lbc(1,level),                         &
     &                    work_u(1,1),intf_work(1,1),                   &
     &                    LEN_INTF_P,MODEL_LEVELS,TEMP,TEMP,.FALSE.,    &
     &                    1,LEN_INTF_P)
            End If

            If (mype == pe_for_var(2) ) Then     !  v
! DEPENDS ON: v_int
              Call V_INT (p_src,p_lbc(1,level),                         &
     &                    work_v(1,1),intf_work(1,2),                   &
     &                    LEN_INTF_P,MODEL_LEVELS,TEMP,TEMP,.FALSE.,    &
     &                    1,LEN_INTF_P)
            End If

            IADDR_V = IADDR + INTF_P_LEVELS(JINTF)*LEN_INTF_U
!
            if(mype == pe_for_var(1)) then
              do i=1,len_intf_u
                intf_data(iaddr+i-1)   = intf_work(i,1)
              end do
            endif

            if(mype == pe_for_var(2)) then
              do i=1,len_intf_u
                intf_data(iaddr_v+i-1) = intf_work(i,2)
              enddo
            endif

            IADDR = IADDR + LEN_INTF_U

          ENDDO

          IADDR = IADDR + INTF_P_LEVELS(JINTF)*LEN_INTF_U

        ELSE

          DO LEVEL=1,MODEL_LEVELS

            DO I=1,LEN_INTF_U
              INTF_DATA(IADDR+I-1) = WORK_U(1,LEVEL)
            ENDDO
            IADDR_V = IADDR + MODEL_LEVELS*LEN_INTF_U
            DO I=1,LEN_INTF_U
              INTF_DATA(IADDR_V+I-1)=WORK_V(1,LEVEL)
            ENDDO
            IADDR = IADDR + LEN_INTF_U

          ENDDO

          IADDR = IADDR + MODEL_LEVELS*LEN_INTF_U

        ENDIF

!     Collect the interpolated V data on PE 0

      iaddr_v = iaddr_u + len_intf_uv_data

      If (mype==pe_for_var(2)) Then ! Send a copy to PE 0
        info=GC_NONE
        call gc_rsend (96, len_intf_u*model_levels,                     &
     &                 pe_for_var(1),info,                              &
     &                 intf_data(iaddr_v), intf_data(iaddr_v) )
      End If

      If (mype==pe_for_var(1)) Then ! Receive a copy from PE 1
        info=GC_NONE
        call gc_rrecv (96, len_intf_u*model_levels,                     &
     &                 pe_for_var(2), info,                             &
     &                 intf_data(iaddr_v), intf_data(iaddr_v) )
      End If

#if defined(TIME_LBC)
        if(mype == pe_for_var(1))                                       &
     &   write(0,*)'Time for Vertical on U was ',                       &
     &   (rtc()-t3)/ticks_per_second
        if(mype == pe_for_var(2))                                       &
     &   write(0,*)'Time for Vertical on V was ',                       &
     &   (rtc()-t3)/ticks_per_second
#endif
      ENDIF  ! IF (mype == pe_for_var(1) .or. mype == pe_for_var(2))


      If (Intf_Vert_Interp(jintf)) Then

! -----------------------------------------
! ThetaL and QT are on theta levels
! Update P_SRC to contain P on theta levels
! -----------------------------------------

        Gather_PE = pe_for_var(3)

! DEPENDS ON: calc_p_src
        Call Calc_P_Src (                                               &
#include "arginfa.h"
     &       ip_p, theta_off_size, len_intf_p,                          &
     &       gather_PE, pe_for_level_uv, local_level,                   &
     &       exner_p_theta, p_src )

! --------------------------------------
! Vert Int for Theta is done on PE 2
! Vert Int for QT    is done on PE 3
! P_SRC is set up on PE 2 ; copy to PE 3
! --------------------------------------

        If (mype==pe_for_var(3)) Then ! Send a copy to PE 3
          info=GC_NONE
          call gc_rsend (97,len_intf_p*model_levels,                    &
     &                   pe_for_var(4),info,                            &
     &                   p_src,p_src)
        End If

        If (mype==pe_for_var(4)) Then ! Receive a copy from PE 2
          info=GC_NONE
          call gc_rrecv (97,len_intf_p*model_levels,                    &
     &                   pe_for_var(3),info,                            &
     &                   p_src,p_src)
        End If

      End If !  If Vert Interp

      iaddr=iaddr_theta
      if (mype  ==  pe_for_var(3)) then
#if defined(TIME_LBC)
        t3=rtc()
#endif

!L 1.3.2 Vertical interpolation - thetal

        IF(INTF_VERT_INTERP(JINTF)) THEN

! input level pressures already set up for winds

! Calculate pressure and exner pressure at output half levels
          DO LEVEL=1,INTF_P_LEVELS(JINTF)+1
            DO I=1,LEN_INTF_P
              P_HALF_TMP(I,LEVEL)=                                      &
     &          INTF_AKH(LEVEL,JINTF)+INTF_BKH(LEVEL,JINTF)*            &
     &          INTF_PSTAR(I)
              P_HALF_TMP_wk(I,LEVEL)=(P_HALF_TMP(I,LEVEL)/PREF)
            ENDDO
          ENDDO
          DO LEVEL=1,INTF_P_LEVELS(JINTF)+1
            DO I=1,LEN_INTF_P
              KAPPA_HALF_wk(I,LEVEL)=KAPPA
            ENDDO
          ENDDO
          n_input=LEN_INTF_P*(INTF_P_LEVELS(JINTF)+1)
#if defined(VECTLIB)
          call rtor_v                                                   &
     &     (n_input,P_HALF_TMP_wk,KAPPA_HALF_wk,P_EXNER_HALF_TMP)
#else
          DO LEVEL=1,INTF_P_LEVELS(JINTF)+1
            DO I=1,LEN_INTF_P
              P_EXNER_HALF_TMP(I,LEVEL)=P_HALF_TMP_wk(I,LEVEL)**KAPPA
            ENDDO
          ENDDO
#endif

! Convert input theta to temperature
          DO LEVEL=1,MODEL_LEVELS
            DO I=1,LEN_INTF_P
              P_TMP_wk(I,LEVEL)=P_SRC(I,LEVEL)/PREF
            ENDDO
          ENDDO
          DO LEVEL=1,MODEL_LEVELS
            DO I=1,LEN_INTF_P
              KAPPA_wk(I,LEVEL)=KAPPA
            ENDDO
          ENDDO
          n_input=LEN_INTF_P*MODEL_LEVELS
#if defined(VECTLIB)
          call rtor_v(n_input,P_TMP_wk,KAPPA_wk,P_TMP_wk)
#else
          DO LEVEL=1,MODEL_LEVELS
            DO I=1,LEN_INTF_P
              P_TMP_wk(I,LEVEL)=P_TMP_wk(I,LEVEL)**KAPPA
            ENDDO
          ENDDO
#endif
          DO LEVEL=1,MODEL_LEVELS
            DO I=1,LEN_INTF_P
              Work_ThetaL(I,LEVEL) = Work_ThetaL(I,LEVEL) *             &
     &                               P_TMP_wk(I,LEVEL)
            ENDDO
          ENDDO

          DO LEVEL=1,INTF_P_LEVELS(JINTF)

! DEPENDS ON: v_int
            CALL V_INT(p_src,p_lbc(1,level),                            &
     &                 work_thetal(1,1),intf_data(IADDR),               &
     &                 LEN_INTF_P,MODEL_LEVELS,TEMP,TEMP,.FALSE.,       &
     &                 1,LEN_INTF_P)

! Convert output temperature to theta

            DO I=1,LEN_INTF_P
              INTF_DATA(IADDR+I-1) = INTF_DATA(IADDR+I-1) /             &
     &        P_EXNER_C(P_EXNER_HALF_TMP(I,LEVEL+1),                    &
     &        P_EXNER_HALF_TMP(I,LEVEL),P_HALF_TMP(I,LEVEL+1),          &
     &        P_HALF_TMP(I,LEVEL),KAPPA)
            ENDDO

            IADDR=IADDR+LEN_INTF_P

          ENDDO

        ELSE

          DO LEVEL=1,MODEL_LEVELS
            DO I=1,LEN_INTF_P
              INTF_DATA(IADDR+I-1)=WORK_ThetaL(I,LEVEL)
            ENDDO
            IADDR=IADDR+LEN_INTF_P
          ENDDO

        ENDIF


#if defined(TIME_LBC)
        write(0,*)'Time for Vertical on THETA was ',                    &
     &   (rtc()-t3)/ticks_per_second
#endif
      ENDIF  ! IF (mype  ==  pe_for_var(3))

! Send the interpolated ThetaL data to PE 0

      If (mype==pe_for_var(3)) Then ! Send a copy to PE 0
        info=GC_NONE
        call gc_rsend (98, len_intf_p*model_levels,                     &
     &                 pe_for_var(1), info,                             &
     &                 intf_data(iaddr_theta), intf_data(iaddr_theta) )
      End If

      If (mype==pe_for_var(1)) Then ! Receive a copy from PE 2
        info=GC_NONE
        call gc_rrecv (98, len_intf_p*model_levels,                     &
     &                 pe_for_var(3), info,                             &
     &                 intf_data(iaddr_theta), intf_data(iaddr_theta) )
      End If

      iaddr=iaddr_qt
      IF (mype  ==  pe_for_var(4)) THEN
#if defined(TIME_LBC)
        t3=rtc()
#endif

!L 1.4.2 Vertical interpolation - QT

        IF (INTF_VERT_INTERP(JINTF)) THEN

!  input level pressures already calculated for thetal

          DO LEVEL=1,INTF_Q_LEVELS(JINTF)
!  set up output level pressure
            DO I=1,LEN_INTF_P
              P_OUT(I) =                                                &
     &        INTF_AK(LEVEL,JINTF)+INTF_PSTAR(I)*INTF_BK(LEVEL,JINTF)
            ENDDO

! DEPENDS ON: v_int
            CALL V_INT(p_src,p_lbc(1,level),                            &
     &                 work_qt(1,1),intf_data(IADDR),                   &
     &                 LEN_INTF_P,WET_LEVELS,TEMP,TEMP,.FALSE.,         &
     &                 1,LEN_INTF_P)


            IADDR=IADDR+LEN_INTF_P

          ENDDO

        ELSE

          DO LEVEL=1,WET_LEVELS
            DO I=1,LEN_INTF_P
              INTF_DATA(IADDR+I-1)=WORK_QT(I,LEVEL)
            ENDDO
            IADDR=IADDR+LEN_INTF_P
          ENDDO

        ENDIF

#if defined(TIME_LBC)
        write(0,*)'Time for Vertical on QT was ',                       &
     &   (rtc()-t3)/ticks_per_second
#endif
      ENDIF  ! IF (mype  ==  pe_for_var(4))

! Send the interpolated QT data to PE 0

      If (mype==pe_for_var(4)) Then ! Send a copy to PE 0
        info=GC_NONE
        call gc_rsend (99, len_intf_p*wet_levels,                       &
     &                 pe_for_var(1), info,                             &
     &                 intf_data(iaddr_qt), intf_data(iaddr_qt) )
      End If

      If (mype==pe_for_var(1)) Then ! Receive a copy from PE 3
        info=GC_NONE
        call gc_rrecv (99, len_intf_p*wet_levels,                       &
     &                 pe_for_var(4), info,                             &
     &                 intf_data(iaddr_qt), intf_data(iaddr_qt) )
      End If

!--ensure that everyone has sent acress their data
      Call GC_SSync (nproc, istat)

#if defined(TIME_LBC)
      if(mype == 0) write(0,*)'Time for Vertical Interpolation was ',   &
     & (rtc()-t2)/ticks_per_second

      if(mype == 0) write(0,*)'Time to Collect LBC Variables was ',     &
     & (rtc()-t1)/ticks_per_second
#endif


!L 1.5 TRACERS

!--resume collective processing
6000  continue




!L 2.0 Update Information in Headers

!L     Open boundary output file if reinitialised during run

      IF (FT_STEPS(NFTOUT) >  0) THEN
        STRING=MODEL_FT_UNIT(NFTOUT)
        PPNAME=STRING(18:31)
        LEN_PPNAME=LEN(PPNAME)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,PPNAME,LEN_PPNAME,1,1,ICODE)
        IF (ICODE /= 0) THEN
          CMESSAGE="GEN_INTF: Error opening preassigned boundary file"
          GO TO 999   !  Return
        ENDIF

!      Determine position where to Buffer out data to

       NTIME=FT_LASTFIELD(NFTOUT)+1
      ELSE
       NTIME=FT_LASTFIELD(NFTOUT)+1
!      A_IO=UNIT(NFTOUT)  ! Only valid with CRAY-specific BUFFER IN/OUT

      ENDIF

!L 2.1 Fixed length header
      FIXHD_INTFA(152,JINTF) = OLD_INTF_LOOKUPSA*NTIME
      FIXHD_INTFA(161,JINTF) = LEN_INTF_DATA*NTIME

!L 2.2 Integer Constants
      INTHD_INTFA(3,JINTF) = NTIME

!L 2.3 LOOKUP Table

!     Determine position in LOOKUP table
      LOOKUP_START=FIXHD_INTFA(150,JINTF) +                             &
     &             FIXHD_INTFA(151,JINTF)*OLD_INTF_LOOKUPSA*(NTIME-1)-1

!
!--for well-formed I/O, we must read back the
!  last lookup table entry on disk
!
      if(ntime /= 1) then
! DEPENDS ON: setpos
        call setpos(nftout, lookup_start-len1_lookup, icode)
! DEPENDS ON: buffin
        call buffin(nftout, lookup_intfa(1, 1, jintf), len1_lookup,     &
     &   len_io, a_io)
!--check for errors
        if(a_io /= -1.0 .or. len_io /= len1_lookup) then
! DEPENDS ON: ioerror
          call ioerror('GEN_INTF_A: Buffer in of Last Lookup Header',   &
     &     a_io, len_io, len1_lookup)
          cmessage=' GEN_INTF_A: I/O Error on Read'
          icode=5
          goto 999
        endif
!--compute the new disk address from the last address and length
        disk_address=lookup_intfa(lbegin, 1, jintf)+                    &
     &               lookup_intfa(lbnrec, 1, jintf)
      else
        disk_address=fixhd_intfa(160, jintf)-1
      endif
!--round this disk to ensure we start on a sector boundary
      disk_address=((disk_address+um_sector_size-1)/                    &
     & um_sector_size)*um_sector_size
!--zero the maximum output record size
      len_buf=0
!     Check that there is enough space for this entry in LOOKUP table
      IF (FIXHD_INTFA(150,JINTF)+                                       &
     &    FIXHD_INTFA(151,JINTF)*FIXHD_INTFA(152,JINTF) >               &
     &   FIXHD_INTFA(160,JINTF)) THEN
        CMESSAGE=' GEN_INTF: Insufficient space for headers in boundary &
     &                       dataset.'
        ICODE=1
        GO TO 999   !  Return
      ENDIF

        START_ADDR = FIXHD_INTFA(161,JINTF)-LEN_INTF_DATA+1

        do var=1,OLD_INTF_LOOKUPSA

! Set STASHCODE for variable required for interfacing

          CODE=ITEM_INTFA(VAR)

          Do I=1,45
           LOOKUP_INTFA(I,VAR,JINTF)=IMDI
          ENDDO
          Do I=46,LEN1_LOOKUP
           LOOKUP_INTFA(I,VAR,JINTF)=0
          ENDDO

          SEC = STEPim(a_im) * SECS_PER_PERIODim(a_im) /                &
     &          STEPS_PER_PERIODim(a_im)

! DEPENDS ON: sec2time
          CALL SEC2TIME(0,SEC,BASIS_TIME_DAYS,BASIS_TIME_SECS,          &
     &                  YY,MM,DD,HR,MN,SS,DAY_NO,LCAL360)

!         Validity time of this field
          LOOKUP_INTFA(LBYR ,VAR,JINTF) = YY
          LOOKUP_INTFA(LBMON,VAR,JINTF) = MM
          LOOKUP_INTFA(LBDAT,VAR,JINTF) = DD
          LOOKUP_INTFA(LBHR ,VAR,JINTF) = HR
          LOOKUP_INTFA(LBMIN,VAR,JINTF) = MN
          LOOKUP_INTFA(LBDAY,VAR,JINTF) = DAY_NO

          LOOKUP_INTFA(LBYRD ,VAR,JINTF) = A_FIXHD(21)
          LOOKUP_INTFA(LBMOND,VAR,JINTF) = A_FIXHD(22)
          LOOKUP_INTFA(LBDATD,VAR,JINTF) = A_FIXHD(23)
          LOOKUP_INTFA(LBHRD ,VAR,JINTF) = A_FIXHD(24)
          LOOKUP_INTFA(LBMIND,VAR,JINTF) = A_FIXHD(25)
          LOOKUP_INTFA(LBDAYD,VAR,JINTF) = A_FIXHD(27)

!! GENERALISED FOR MIXED PHASE PRECIP SCHEME
          IF(VAR == 1) THEN
            LEN_DATA = LEN_INTF_P
          ELSE IF(VAR == 2.OR.VAR == 3) THEN
            LEN_DATA = LEN_INTF_U*INTF_P_LEVELS(JINTF)
          ELSE IF(VAR == 4) THEN
            LEN_DATA = LEN_INTF_P*INTF_P_LEVELS(JINTF)
          ELSE IF(VAR == 5) THEN
            LEN_DATA = LEN_INTF_P*INTF_Q_LEVELS(JINTF)
          ELSE IF(VAR >  5 .AND. CODE >  60) THEN
            LEN_DATA = LEN_INTF_P*INTF_TR_LEVELS(JINTF)
          ELSE IF(VAR >  5 .AND. CODE == 12) THEN
            LEN_DATA = LEN_INTF_P*INTF_Q_LEVELS(JINTF)
          END IF
          LOOKUP_INTFA(LBLREC,VAR,JINTF) = LEN_DATA
!         New packing information from UM Version 2.8
          N1 = 0   !  Data not packed
          IF (LPACK_32B) N1 = 2  ! Data packed as 32 bits
          IF (LPACK_PPXREF) THEN
! DEPENDS ON: exppxi
            N1 = EXPPXI(atmos_im,0,item_intfa(var),ppx_dump_packing,    &
#include "argppx.h"
     &                  icode,cmessage)
          ENDIF
          N2 = 0   !  Data not compressed
          N3 = 0   !  Compression definition
          N4 = 0   !  Number format
          N5 = 0   !  Not used
          NPACK = N5*10000 + N4*1000 +N3*100 + N2*10 + N1
          LOOKUP_INTFA(LBPACK,VAR,JINTF)= NPACK
          LOOKUP_INTFA(LBREL,VAR,JINTF)=2

          If (var == 1) Then
            lookup_intfa(lbfc,var,jintf) = 8
          Else If (var == 2) Then
            lookup_intfa(lbfc,var,jintf) = 56
          Else If (var == 3) Then
            lookup_intfa(lbfc,var,jintf) = 57
          Else If (var == 4) Then
            lookup_intfa(lbfc,var,jintf) = 19
          Else If (var == 5) Then
            lookup_intfa(lbfc,var,jintf) = 95
          End If
!
!--make sure that the LBC complete record is well formed
!
!--set the disk address
          lookup_intfa(lbegin, var, jintf)=disk_address
!--fetch the data field length, allowing for packing
          if(mod(lookup_intfa(lbpack, var, jintf), 10) == 2) then
            disk_length=(lookup_intfa(lblrec, var, jintf)+1)/2
          else
            disk_length=lookup_intfa(lblrec, var, jintf)
          endif
!--update the maximum record size
          len_buf=len_buf+disk_length
!--store the rounded-up length
          lookup_intfa(lbnrec, var, jintf)=disk_length
!--update the disk address
          disk_address=disk_address+disk_length
          LOOKUP_INTFA(LBCODE,VAR,JINTF)=1
            IF(VAR == 2.OR.VAR == 3) THEN
              LOOKUP_INTFA(LBCODE,VAR,JINTF)=2
            END IF
          LOOKUP_INTFA(LBHEM,VAR,JINTF)=99
          LOOKUP_INTFA(LBROW,VAR,JINTF)=INTFWIDTHA(JINTF)
          LOOKUP_INTFA(LBNPT,VAR,JINTF)=LEN_INTF_P/INTFWIDTHA(JINTF)
          IF (VAR == 2.OR.VAR == 3) THEN
            LOOKUP_INTFA(LBNPT,VAR,JINTF)=LEN_INTF_U/INTFWIDTHA(JINTF)
          END IF
          LOOKUP_INTFA(LBLEV,VAR,JINTF)=-1
          LOOKUP_INTFA(DATA_TYPE,VAR,JINTF)=1
          LOOKUP_INTFA(NADDR,VAR,JINTF) = START_ADDR
          LOOKUP_INTFA(ITEM_CODE,VAR,JINTF) = CODE
          START_ADDR = START_ADDR + LOOKUP_INTFA(LBLREC,VAR,JINTF)

        ENDDO


!L 3.0 Pack data as required

        IADDR = 1
        LEN_DATA = 0
        do var=1,OLD_INTF_LOOKUPSA
          IF (MOD(LOOKUP_INTFA(LBPACK,VAR,JINTF),10) == 2) THEN
            IF (mype  ==  0) THEN
! DEPENDS ON: pack21
            CALL PACK21(LOOKUP_INTFA(LBLREC,VAR,JINTF),                 &
                        INTF_DATA(IADDR),INTF_DATA(LEN_DATA+1))
            ENDIF
!--the (+1) in the expression below is unnecessary, since
!  LBC data is composed of two rows NS and two rows EW, and
!  thus always has an even number of data points.  If this
!  is not true, then READFLDS will either get the data one
!  out downwards if the (+1) is omitted, or one word upwards
!  if the (+1) is added.  In other words, the packing will
!  cause either one word to be omitted or one word added in
!  the data after the read.  This is because READFLDS reads
!  and converts the whole LBC record at one go, rather than
!  as a series of separate records.
            len_data = len_data+(lookup_intfa(lblrec,var,jintf)+1)/2
!--check that we are not packing an odd nuber of words
            if((lookup_intfa(lblrec,var,jintf)/2)*2  /=                 &
     &       lookup_intfa(lblrec,var,jintf)) then
              write(6,7734) lookup_intfa(lblrec,var,jintf)
7734          format(/'LBC Data contains ',i10,' Words, which is',      &
     &         ' an Odd Number which is not allowed for 32-bit',        &
     &         ' Packing')
            endif
          ELSE
            IF (LEN_DATA+1 <  IADDR) THEN
            IF (mype  ==  0) THEN
              DO J = 1,LOOKUP_INTFA(LBLREC,VAR,JINTF)
                INTF_DATA(LEN_DATA+J) = INTF_DATA(IADDR+J-1)
              ENDDO
            ENDIF
            ENDIF
            LEN_DATA = LEN_DATA+LOOKUP_INTFA(LBLREC,VAR,JINTF)
          ENDIF
          IADDR = IADDR+LOOKUP_INTFA(LBLREC,VAR,JINTF)
        ENDDO
!L 4.0 Write out headers/data

!L 4.1 Fixed length header

        IADDR = 0
! DEPENDS ON: setpos
        CALL SETPOS (NFTOUT,IADDR,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FIXHD_INTFA(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of fixed length header',A_IO,LEN_IO, &
     &                  LEN_FIXHD)
          CMESSAGE=' GEN_INTF: I/O ERROR '
          ICODE=2
          GO TO 999   !  Return
        END IF

!L 4.2 Integer constants

! DEPENDS ON: buffout
        CALL BUFFOUT (NFTOUT,INTHD_INTFA(1,JINTF),                      &
     &                LBC_LEN_INTHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LBC_LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of integer header',A_IO,LEN_IO,      &
     &                  LBC_LEN_INTHD)
          CMESSAGE=' GEN_INTF: I/O ERROR '
          ICODE=3
          GO TO 999   !  Return
        END IF

!L 4.3 PP headers in LOOKUP table
! DEPENDS ON: setpos
        CALL SETPOS(NFTOUT,LOOKUP_START,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LOOKUP_INTFA(1,1,JINTF),                    &
     &               LEN1_LOOKUP*OLD_INTF_LOOKUPSA,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN1_LOOKUP*OLD_INTF_LOOKUPSA)THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of PP header',A_IO,LEN_IO,           &
     &                  LEN1_LOOKUP*OLD_INTF_LOOKUPSA)
          CMESSAGE=' GEN_INTF: I/O ERROR '
          ICODE=4
          GO TO 999   !  Return
        END IF

!L 4.4 Interface data
!       Determine position in data section

        DATA_START =                                                    &
     &   lookup_intfa(lbegin, 1, jintf)
!--round this disk length to a multiple of the sector size
        len_data=((len_data+um_sector_size-1)/                          &
     &    um_sector_size)*um_sector_size
! DEPENDS ON: setpos
        CALL SETPOS(NFTOUT,DATA_START,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,INTF_DATA(1),LEN_DATA,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_DATA) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of boundary data',A_IO,LEN_IO,       &
     &                  LEN_DATA)
          CMESSAGE=' GEN_INTF: I/O ERROR '
          ICODE=51
          GO TO 999   !  Return
        END IF

!L     Close boundary output file if reinitialised during run
      IF (FT_STEPS(NFTOUT) >  0) THEN
        LEN_PPNAME=LEN(PPNAME)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,PPNAME,LEN_PPNAME,1,0,ICODE)
      END IF

!L     Update FT_LASTFIELD
      FT_LASTFIELD(NFTOUT) = FT_LASTFIELD(NFTOUT) + 1

 999  RETURN
      END SUBROUTINE GEN_INTF_A_OLD_LBCS
#endif
