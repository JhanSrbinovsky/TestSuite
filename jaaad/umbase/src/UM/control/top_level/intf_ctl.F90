#if defined(CONTROL) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INTF_CTL ----------------------------------------------
!LL
!LL  Purpose: Initialises namelist and sets switch for the writing of
!LL           boundary data.
!LL
!LL  Tested under complier:   cft77 5.0.2
!LL  Tested under OS version: Unicos 6.5.1.a
!LL
!LL  Author:  R.G. Jones         Date: 7 July 1992
!LL
!LL  Programming standard: UM Doc Paper 3, version 4 (05/02/92)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:                 D8
!LL
!LL  External documentation:       UM Doc Paper D8
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments--------------------------------------------

      SUBROUTINE INTF_CTL (                                             &
#include "arginfa.h"
#include "argduma.h"
#include "argptra.h"
     &                     ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parparm.h"
#include "typsize.h"
#include "typinfa.h"
#include "typduma.h"
#include "typptra.h"

      INTEGER                                                           &
     &  ICODE                     ! Out - Return Code

      CHARACTER*80 CMESSAGE     ! Out - Error message on failure
!*----------------------------------------------------------------------
!  Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "cmaxsizo.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#include "cprintst.h"
#if defined(ATMOS)
#include "cintfa.h"
#endif
!
!-----------------------------------------------------------------------
#include "c_mdi.h"
!-----------------------------------------------------------------------
!
!  Local variables
!
      INTEGER JINTF               ! Interface area index
      INTEGER INTF_AREA_NO        ! Interface area number
      INTEGER INTF_UNIT_NO        ! Interface unit number
      INTEGER J                   ! For setting LBC_UNIT_NO
      INTEGER J_LBC               ! For setting lambda, phi

#if defined(ATMOS)
      INTEGER STEP                ! Current atmos. step
      INTEGER A_STEPS_PER_HR ! steps per hour for atmosphere
      INTEGER A_SECS_PER_STEP ! secs per step for atmos. sub_model

      INTEGER                                                           &
     & IERR                                                             &
                               !Return code from ABCALC
     &,NQ               ! number of water variables in interface files

!  Namelist for interface constants
!
#include "cnaminfa.h"
#endif

        Character (Len=*), Parameter :: RoutineName = 'IntfCtl'
        Integer, Parameter  :: nft_vertlevs = 90
        Integer, Parameter  :: nft_HorzGrid = 129
        Integer, Parameter  :: idummy = 1
        Logical             :: l_exist
        Character (Len=80)  :: Filename
        Integer             :: ErrorStatus
        Integer             :: status
        Integer             :: jlev

        ICODE=0
        CMESSAGE=' '

!L 1. Atmosphere data

#if defined(ATMOS)
      IF (LLBOUTim(a_im)) THEN

!     Initialise namelist parameters to cater for no boundary datasets
      DO JINTF = 1,MAX_N_INTF_A
        A_INTF_START_HR(JINTF) = 0
        A_INTF_END_HR  (JINTF) = 0
        A_INTF_FREQ_HR (JINTF) = 0
        A_INTF_FREQ_MN (JINTF) = 0
        A_INTF_FREQ_SC (JINTF) = 0
        INTF_METH_LEV_CALC(JINTF) = 5
        INTF_MAX_SIG_HLEV(JINTF)  = 0
        INTF_MIN_PRS_HLEV(JINTF)  = 0
        INTF_PACK(JINTF) = 1
        LBC_ND(JINTF)    = 1
        INTF_PACK(JINTF) = 1            !  32 bit packing
        INTF_L_VAR_LBC(JINTF)=.false.
        INTF_HORZGRID(JINTF)=' '
        INTF_VERTLEVS(JINTF)=' '
        DO J=1,MAX_INTF_LEVELS+1
          INTF_ETAH(J,JINTF)=RMDI
        ENDDO
        DO J=1,MAX_INTF_LBCROW_LENGTH
          LAMBDA_INTF_P(J,JINTF)=RMDI
          LAMBDA_INTF_U(J,JINTF)=RMDI
        ENDDO
        DO J=1,MAX_INTF_LBCROWS
          PHI_INTF_P(J,JINTF)=RMDI
          PHI_INTF_V(J,JINTF)=RMDI
        ENDDO
      ENDDO
      
      DO J=1,MAX_INTF_LBCROW_LENGTH
        LAMBDA_LBC_P(J)=RMDI
        LAMBDA_LBC_U(J)=RMDI
      ENDDO
      DO J=1,MAX_INTF_LBCROWS
        PHI_LBC_P(J)=RMDI
        PHI_LBC_V(J)=RMDI 
      ENDDO
      LBC_Q_MIN = 0.0

!  Read namelist for interface

      REWIND (5)
      READ(5,INTFCNSTA)

      If (PrintStatus >= PrStatus_Oper ) Then
        Write (6,INTFCNSTA)
      End If
      
      Do jintf = 1, n_intf_a 
         
        If ( Intf_L_Var_Lbc(jintf) ) Then 
          FileName = Trim( INTF_HorzGrid(jintf) ) 
! Check file exists  
          Inquire( file=FileName, exist=l_exist )   
          If ( .Not. l_exist ) Then 
            Write (6,*) 'LBC VarRes Grids: ',FileName   
            ErrorStatus = 20 
            Cmessage = 'LBC VarRes Grids Namelist file does not exist!'
            Call Ereport( RoutineName, ErrorStatus, Cmessage )  
          End If       
! Open the file containing horizontal grids 
          Open( Unit=nft_HorzGrid, File=FileName, IOstat=status )

! Read LBC VarGrid namelist 
          Read ( Unit=nft_HorzGrid, Nml=Lbcgrids) 

! Write out namelist for diagnostic  
          If (PrintStatus >= PrStatus_Normal) Then  
            Write ( Unit = 6, Nml = Lbcgrids) 
          End If
          
          do J_lbc = 1, INTF_ROW_LENGTH(jintf)
            Lambda_intf_p(J_lbc,jintf)= Lambda_lbc_p(J_lbc)
            Lambda_intf_u(J_lbc,jintf)= Lambda_lbc_u(J_lbc)
          end do
          do J_lbc = 1, INTF_P_ROWS(jintf)
            Phi_intf_p(J_lbc,jintf)= Phi_lbc_p(J_lbc)
            Phi_intf_v(J_lbc,jintf)= Phi_lbc_v(J_lbc)
          end do
          
          Intf_FirstLong(jintf) =RMDI
          Intf_ewspace(jintf)   =RMDI
          Intf_FirstLat(jintf)  =RMDI
          Intf_nsspace(jintf)   =RMDI
            
        Else     ! Regular Grid LBC as UNUI input 
          
          If (PrintStatus >= PrStatus_Normal) Then    
            Write (6,*) ' Regular Grid  LBC output'  
          End If
                 
          do J_lbc = 1, INTF_ROW_LENGTH(jintf)
            Lambda_intf_p(J_lbc,jintf) = Intf_FirstLong(jintf) +        &
                         (J_lbc - 1.0)*Intf_ewspace(jintf)  
            Lambda_intf_u(J_lbc,jintf) = Intf_FirstLong(jintf) +        & 
                         (J_lbc - 0.5)*Intf_ewspace(jintf) 
          end do
          do J_lbc = 1, INTF_P_ROWS(jintf)
            Phi_intf_p(J_lbc,jintf) = Intf_FirstLat(jintf) +            &
                         (J_lbc - 1.0)*Intf_nsspace(jintf) 
            Phi_intf_v(J_lbc,jintf) = Intf_FirstLat(jintf) +            &
                         (J_lbc - 0.5)*Intf_nsspace(jintf) 
          end do 
        End if
      End do

      If (PrintStatus >= PrStatus_Normal ) Then
        write (6,*) ' '
        Do jintf =1, n_intf_a
          write (6,*) 'LBC Area ',lbc_stream_a(jintf),' : ',            &
     &    intf_vertlevs(jintf)
        End Do
      End If

      Intf_Vert_Interp(:) = .false.
      do jintf = 1, n_intf_a


          FileName = Trim ( Intf_VertLevs(jintf) )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together

          Inquire( file=FileName, exist=l_exist )

          If ( .Not. l_exist ) Then
            Write (6,*) ' LBC Vertical Levels file: ',FileName
            ErrorStatus = 20
            Cmessage = 'Vertical Levels Namelist file does not exist!'
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

! Open the file containing vertical levels
          Open( Unit=nft_vertlevs, File=FileName, IOstat=status )

       If (lbc_nd(jintf) == 1) Then
         eta_theta(:)=0.0
         eta_rho(:)  =0.0
! Read VERTLEVS namelist
          Read ( Unit=nft_vertlevs, nml=vertlevs )

!         store contents
          lbc_z_top_model(jintf) = z_top_of_model
          lbc_first_r_rho(jintf) =                                      &
     &        first_constant_r_rho_level
          do jlev = 1, intf_p_levels(jintf) + 1
            lbc_eta_theta(jlev,jintf) = eta_theta(jlev)
          enddo
          do jlev = 1, intf_p_levels(jintf)
            lbc_eta_rho(jlev,jintf) = eta_rho(jlev)
          enddo

       Else  !  Old LBCs

         etah(:) = 0.0
         Read (Unit=nft_vertlevs, nml = OLDVERT)

         intf_vert_interp(jintf)   = vert_interp
         intf_meth_lev_calc(jintf) = meth_lev_calc
         intf_max_sig_hlev(jintf)  = max_sig_hlev
         intf_min_prs_hlev(jintf)  = min_prs_hlev
         Do jlev = 1, intf_p_levels(jintf)+1
           intf_etah(jlev,jintf) = etah(jlev)
         End Do

      End If  !  lbc_nd
          Close ( Unit=nft_vertlevs)

      enddo

! Check if vertical interpolation required.
! When running makebc, this is called by get_bc
#if !defined(MAKEBC)
! DEPENDS ON: lbc_chk_vert_interp
      Call LBC_Chk_Vert_Interp (                                        &
#include "arginfa.h"
#include "argduma.h"
#include "argptra.h"
     & idummy )
#endif

!L Initialise LBC_UNIT_NO_A and FT_OUTPUT.

      DO JINTF=1,MAX_N_INTF_A
        LBC_UNIT_NO_A(JINTF) = 0
      ENDDO
      DO JINTF=1,N_INTF_A
        J = 139+LBC_STREAM_A(JINTF)
        LBC_UNIT_NO_A(JINTF) = J
        FT_Output(J) = 'Y'
      ENDDO

      write (6,*) ' '
      write (6,*) 'LBC_STREAM_A and LBC_UNIT_NO_A'
      do jintf=1,n_intf_a
      write (6,*) jintf,' stream no ',lbc_stream_a(jintf),              &
     &                  ' unit_no ',lbc_unit_no_a(jintf)
      enddo
      write (6,*) ' '
      write (6,*) 'FT_Output for Boundary Files.'
      do j=140,147
      write (6,*) ' Unit No ',J,' FT_Output ',ft_output(j)
      enddo

!L Calculate INTF_AK,INTF_BK,INTF_AKH,INTF_BKH from INTF_ETAH
!L
! Initialise arrays
      DO JINTF = 1,N_INTF_A
        DO J=1,MAX_INTF_MODEL_LEVELS
          INTF_AK(J,JINTF)=RMDI
          INTF_BK(J,JINTF)=RMDI
        END DO
        DO J=1,MAX_INTF_MODEL_LEVELS+1

          INTF_AKH(J,JINTF)=RMDI
          INTF_BKH(J,JINTF)=RMDI
        ENDDO
      ENDDO

      DO JINTF = 1,N_INTF_A

        if ( lbc_nd(jintf) /= 1 ) then  !  Not ND LBCs

        IF(INTF_VERT_INTERP(JINTF)) THEN ! Vertical interpolation needed

! DEPENDS ON: abcalc
          CALL ABCALC(INTF_METH_LEV_CALC(JINTF),1,1                     &
     &,               INTF_P_LEVELS(JINTF)                              &
     &,               INTF_ETAH(INTF_MIN_PRS_HLEV(JINTF),JINTF)         &
     &,               INTF_ETAH(INTF_MAX_SIG_HLEV(JINTF),JINTF)         &
     &,               INTF_ETAH(1,JINTF)                                &
     &,               INTF_AK(1,JINTF),INTF_BK(1,JINTF)                 &
     &,               INTF_AKH(1,JINTF),INTF_BKH(1,JINTF),IERR)

          IF(IERR /= 0) THEN
            WRITE(6,*) ' *ERROR*  IN ABCALC FROM INTF_CTL. IERR  '      &
     &        ,IERR
            WRITE(6,*) '   CHECK YOUR ATMOS LEVEL SPEC FOR MODEL'
            ICODE = 2
            WRITE (cmessage,*) 'INTF_CTL : Error IN ABCALC.'
            GO TO 9999   !  Return
          END IF

          WRITE(6,*) 'INTF_AK='
          WRITE(6,'(3(E22.15,'',''))')                                  &
     &       (INTF_AK(J,JINTF),J=1,INTF_P_LEVELS(JINTF))
          WRITE(6,*) 'INTF_BK='
          WRITE(6,'(3(E22.15,'',''))')                                  &
     &       (INTF_BK(J,JINTF),J=1,INTF_P_LEVELS(JINTF))
          WRITE(6,*) 'INTF_AKH='
          WRITE(6,'(3(E22.15,'',''))')                                  &
     &     (INTF_AKH(J,JINTF),J=1,INTF_P_LEVELS(JINTF)+1)
          WRITE(6,*) 'INTF_BKH='
          WRITE(6,'(3(E22.15,'',''))')                                  &
     &     (INTF_BKH(J,JINTF),J=1,INTF_P_LEVELS(JINTF)+1)

          endif   !  lbc_nd

        END IF

       IF (INTF_PACK(JINTF) <  0 .OR. INTF_PACK(JINTF) >  2) THEN
          CMESSAGE = 'INTFCTL : Invalid value for INTF_PACK.'
          ICODE = 1
          WRITE (6,*) ' INTF_PACK ',INTF_PACK
          WRITE (6,*) ' Valid values for INTF_PACK are 0 and 1'
          GO TO 9999   !  Return
        ENDIF

      END  DO

!     Initialise variables in CTIME comdeck.
      A_STEPS_PER_HR = 3600*STEPS_PER_PERIODim(a_im)/                   &
     &                      SECS_PER_PERIODim(a_im)
! Compute secs per step for atmosphere sub_model
      A_SECS_PER_STEP = SECS_PER_PERIODim(a_im)/                        &
     &                   STEPS_PER_PERIODim(a_im)
      DO JINTF = 1,N_INTF_A
       INTERFACE_FSTEPim(JINTF,a_im) =                                  &
     &   A_INTF_START_HR(JINTF)*A_STEPS_PER_HR
       INTERFACE_LSTEPim(JINTF,a_im) =                                  &
     &   A_INTF_END_HR(JINTF)  *A_STEPS_PER_HR
       INTERFACE_STEPSim(JINTF,a_im) =                                  &
     &  ( 3600*A_INTF_FREQ_HR(JINTF) +                                  &
     &    60*A_INTF_FREQ_MN(JINTF)+A_INTF_FREQ_SC(JINTF) ) /            &
     &    A_SECS_PER_STEP
      ENDDO
!
!     Initialise variables in CINTF comdeck.
      DO JINTF = 1,N_INTF_A

        If (lbc_nd(jintf)==0) Then  !  4.5 LBCs

!       Length of interface field on p* grid
        LEN_INTFA_P(JINTF) =                                            &
     &  (INTF_ROW_LENGTH(JINTF)+INTF_P_ROWS(JINTF)-2*INTFWIDTHA(JINTF)) &
     &  *2*INTFWIDTHA(JINTF)

!       Length of interface field on wind grid
        LEN_INTFA_U(JINTF) = LEN_INTFA_P(JINTF) - 4*INTFWIDTHA(JINTF)

!       Length of interface data
!       Only qT in 4.5 LBCs
        NQ = 1
        LEN_INTFA_DATA(JINTF) = LEN_INTFA_P(JINTF) *                    &
     &   (INTF_P_LEVELS(JINTF) + INTF_Q_LEVELS(JINTF)*NQ + 1 +          &
     &    INTF_TR_LEVELS(JINTF)*TR_VARS) +                              &
     &    LEN_INTFA_U(JINTF)*INTF_P_LEVELS(JINTF)*2

        End If  !   lbc_nd

      ENDDO
!
      STEP = STEPim(a_im)
!
!  Loop over all possible fortran unit numbers to find
!  new output boundary files (i.e. set LNEWBND)

! DEPENDS ON: intf_new_files
      call intf_new_files (20, NUNITS, MAX_N_INTF_A, A_IM,              &
     &    TYPE_LETTER_1, FT_OUTPUT, A_INTF_FREQ_HR, A_INTF_FREQ_MN,     &
     &    A_INTF_FREQ_SC, FT_STEPS, STEP,                               &
     &    FT_FIRSTSTEP, INTERFACE_STEPSim(1,A_IM),                      &
     &    LNEWBND )

      END IF ! LLBOUTim(a_im)
#endif

!L 2. Ocean data (redundant)

!
!L End of routine
!
 9999 CONTINUE
      RETURN
      END SUBROUTINE INTF_CTL
#endif
