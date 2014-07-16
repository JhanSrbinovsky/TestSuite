
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_BL_PT2 ----------------------------------------------
!!!
!!!  Purpose: Calculate increments for
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!  6.1  17/05/04  Pass switches for diagnostic calculations.
!                                                       M. Diamantakis
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_BL_PT2 (                                            &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS                     &
     &,l_ftl,l_fqw,l_taux,l_tauy                                        &
     &,RHOKH,RHOKM_U,RHOKM_V                                            &
     &,RDZ_CHARNEY_GRID,RDZ_U,RDZ_V,GAMMA                               &
     &,CT_CTQ,DQW,DTL,CQ_CM_U,CQ_CM_V,DU,DV                             &
     &,FQW,FTL,TAU_X,TAU_Y,QW,TL                                        &
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                                 ! Local number of points on a row
     &, ROWS                                                            &
                                 ! Local number of rows in a theta field
     &, N_ROWS                                                          &
                                 ! Local number of rows in a v field
     &, OFF_X                                                           &
                                 ! Size of small halo in i
     &, OFF_Y                                                           &
                                 ! Size of small halo in j.
     &,BL_LEVELS                 ! IN No. of atmospheric levels for
!                                     which boundary layer fluxes are
!                                     calculated.

      Logical                                                           &
     & l_ftl                                                            &
     &,l_fqw                                                            &
     &,l_taux                                                           &
     &,l_tauy


      REAL                                                              &
     & RHOKH(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                ! IN Exchange coeff for FTL above
!                                !    surface.
     &,RHOKM_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                             &
!                                ! IN Exchange coefficients for
!                                !    momentum, on U-grid with
!                                !    first and last rows ignored.
!                                !    for K>=2 (from KMKH).
     &,RHOKM_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                           &
!                                ! IN Exchange coefficients for
!                                !    momentum, on V-grid with
!                                !    first and last rows ignored.
!                                !    for K>=2 (from KMKH).
     &,RDZ_CHARNEY_GRID(ROW_LENGTH,ROWS,BL_LEVELS)                      &
!                                ! IN RDZ(,1) is the reciprocal of the
!                                ! height of level 1, i.e. of the
!                                ! middle of layer 1.  For K > 1,
!                                ! RDZ(,K) is the reciprocal
!                                ! of the vertical distance
!                                ! from level K-1 to level K.
     &,RDZ_U(ROW_LENGTH,ROWS,2:BL_LEVELS)                               &
!                                ! IN Reciprocal of the vertical
!                                !    distance from level K-1 to
!                                !    level K. (K > 1) on wind levels
     &,RDZ_V(ROW_LENGTH,N_ROWS,2:BL_LEVELS)                             &
!                                ! IN Reciprocal of the vertical
!                                !    distance from level K-1 to
!                                !    level K. (K > 1) on wind levels
     &,GAMMA(BL_LEVELS)          ! IN Implicit weighting.


      REAL                                                              &
     & CT_CTQ(ROW_LENGTH,ROWS,BL_LEVELS)                                &
!                                ! INOUT Coefficient in T and q
!                                !       tri-diagonal implicit matrix
     &,CQ_CM_U(ROW_LENGTH,ROWS,BL_LEVELS)                               &
!                                ! INOUT Coefficient in U tri-diagonal
!                                !       implicit matrix
     &,CQ_CM_V(ROW_LENGTH,N_ROWS,BL_LEVELS)                             &
!                                ! INOUT Coefficient in V tri-diagonal
!                                !       implicit matrix
     &,DQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                ! INOUT BL increment to q field
     &,DTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                ! INOUT BL increment to T field
     &,DU(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:ROWS+OFF_Y,BL_LEVELS)                                 &
!                                ! INOUT BL increment to u wind field
     &,DV(1-OFF_X:ROW_LENGTH+OFF_X,                                     &
     &    1-OFF_Y:N_ROWS+OFF_Y,BL_LEVELS)                               &
!                                ! INOUT BL increment to v wind field
     &,FQW(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                ! INOUT Flux of QW (ie., for surface,
!                                !       total evaporation). Kg/sq m/s
     &,FTL(ROW_LENGTH,ROWS,BL_LEVELS)                                   &
!                                ! INOUT Flux of TL (ie., for surface,
!                                !       H/Cp where H is sensible heat
!                                !       in W per sq m).
     &,TAU_X(ROW_LENGTH,ROWS,BL_LEVELS)                                 &
!                                ! INOUT x-component of turbulent
!                                !       stress at levels k-1/2;
!                                !       eg. TAUX(,1) is surface stress
!                                !       on UV-grid, 1st and last rows
!                                !       set to "missing data". (N/sq m)
!                                !       IN as "explicit" fluxes from
!                                !       ex_flux_uv, OUT as "implicit
     &,TAU_Y(ROW_LENGTH,N_ROWS,BL_LEVELS)                               &
!                                ! INOUT y-component of turbulent
!                                !       stress at levels k-1/2;
!                                !       eg. TAUX(,1) is surface stress
!                                !       on UV-grid, 1st and last rows
!                                !       set to "missing data". (N/sq m)
!                                !       IN as "explicit" fluxes from
!                                !       ex_flux_uv, OUT as "implicit
     &,QW(ROW_LENGTH,ROWS,BL_LEVELS)                                    &
!                                ! INOUT Total water content (kg per
!                                !       kg air).  From P243.
     &,TL(ROW_LENGTH,ROWS,BL_LEVELS)
!                                ! INOUT Liquid/frozen water
!                                !       temperature (K).  From P243.


! eak

       Real                                                             &
     & xx1(ROW_LENGTH,ROWS),                                            &
     & xx2(ROW_LENGTH,ROWS),                                            &
     & xx3(ROW_LENGTH,ROWS),                                            &
     & xx4(ROW_LENGTH,ROWS),                                            &
     & xx5(ROW_LENGTH,ROWS),                                            &
     & xx6(ROW_LENGTH,ROWS),                                            &
     & xx7(ROW_LENGTH,ROWS),                                            &
     & xx8(ROW_LENGTH,ROWS)

!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------

!  External references :-
      EXTERNAL TIMER

!  Local scalars :-
      INTEGER                                                           &
     & I,J                                                              &
                ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).




!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_BL_PT2 ',3)
      ENDIF
!      print *,'imblpt2 1',dqw(:,:,1),dtl(:,:,1),du(:,:,1),dv(:,:,1)
      xx1 = dqw(:,:,1)
      xx2 = dtl(:,:,1)
      xx3 = du(:,:,1)
      xx4 = dv(:,:,1)
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        DU(I,J,1) = DU(I,J,1) - CQ_CM_U(I,J,1)*TAU_X(I,J,1)
       ENDDO
      ENDDO

      DO J=1,N_ROWS
       DO I=1,ROW_LENGTH
        DV(I,J,1) = DV(I,J,1) - CQ_CM_V(I,J,1)*TAU_Y(I,J,1)
       ENDDO
      ENDDO

      DO K=2,BL_LEVELS

        DO J=1,ROWS
         DO I=1,ROW_LENGTH
          DU(I,J,K) = DU(I,J,K) - CQ_CM_U(I,J,K)*DU(I,J,K-1)
         ENDDO
        ENDDO

        DO J=1,N_ROWS
         DO I=1,ROW_LENGTH
          DV(I,J,K) = DV(I,J,K) - CQ_CM_V(I,J,K)*DV(I,J,K-1)
         ENDDO
        ENDDO

      ENDDO
!      print *,'imblpt2 2',dqw(:,:,1),dtl(:,:,1),du(:,:,1),dv(:,:,1)
!      print *,'imblpt2 3',tl(:,:,1),qw(:,:,1)
      xx5 = dqw(:,:,1)
      xx6 = dtl(:,:,1)
      xx7 = tl(:,:,1)
      xx8 = qw(:,:,1)
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        DTL(I,J,1) = DTL(I,J,1) - CT_CTQ(I,J,1)*FTL(I,J,1)/CP
        TL(I,J,1) = TL(I,J,1) + DTL(I,J,1)
        DQW(I,J,1) = DQW(I,J,1) - CT_CTQ(I,J,1)*FQW(I,J,1)
        QW(I,J,1) = QW(I,J,1) + DQW(I,J,1)
       ENDDO
      ENDDO
!      print 99,dqw(4,2,1),dtl(4,2,1),CT_CTQ(4,2,1),FTL(4,2,1), &
!               TL(4,2,1),QW(4,2,1),DU(4,2,1),DV(4,2,1)
!      print 99,dqw(28,57,1),dtl(28,57,1),CT_CTQ(28,57,1),FTL(28,57,1), &
!               TL(28,57,1),QW(28,57,1),DU(28,57,1),DV(28,57,1)
!      print 99,dqw(35,56,1),dtl(35,56,1),CT_CTQ(35,56,1),FTL(35,56,1), &
!               TL(35,56,1),QW(35,56,1),DU(35,56,1),DV(35,56,1)
!!      print 99,dqw(13,29,1),dtl(13,29,1),CT_CTQ(13,29,1),FTL(13,29,1), &
!!               TL(13,29,1),QW(13,29,1),DU(13,29,1),DV(13,29,1)
!99    format(1x,'imblpt2 4 ',f9.6,f7.2,2x,f8.2,f7.2,2x,f7.1,f9.6,2x,2f6.2) 
!!      print 100,tl(:,:,1),xx7,xx2,dtl(:,:,1),qw(:,:,1),xx8,xx1,         &
!!     & dqw(:,:,1),du(:,:,1),xx3,dv(:,:,1),xx4
!!	100   format('imblpt2',2f7.2,1x,2f6.3,2x,2f9.6,1x,2f10.7,2x,2f6.3,     &
!!	     & 2x,2f6.3)


	      DO K=2,BL_LEVELS
		DO J=1,ROWS
		 DO I=1,ROW_LENGTH

		  DTL(I,J,K) = DTL(I,J,K) - CT_CTQ(I,J,K)*DTL(I,J,K-1)
		  TL(I,J,K) = TL(I,J,K) + DTL(I,J,K)
		  DQW(I,J,K) = DQW(I,J,K) - CT_CTQ(I,J,K)*DQW(I,J,K-1)
		  QW(I,J,K) = QW(I,J,K) + DQW(I,J,K)

		 ENDDO
		ENDDO
	      ENDDO !bl_levels


	! MD
	! Calculate stress and flux diagnostics. The fluxes are
	! calculated only when requested.
	!
	      If ( l_taux ) Then
		DO K=2,BL_LEVELS
		  DO J=1,ROWS
		    DO I=1,ROW_LENGTH
		      TAU_X(I,J,K) = TAU_X(I,J,K)                               &
	     &                     + GAMMA(K)*RHOKM_U(I,J,K)*RDZ_U(I,J,K)       &
	     &                     * ( DU(I,J,K)-DU(I,J,K-1) )
		    ENDDO
		  ENDDO
		ENDDO ! bl_levels
	      Endif

	      If ( l_tauy ) Then
		DO K=2,BL_LEVELS
		  DO J=1,N_ROWS
		    DO I=1,ROW_LENGTH
		      TAU_Y(I,J,K) = TAU_Y(I,J,K)                               &
	     &                     + GAMMA(K)*RHOKM_V(I,J,K)*RDZ_V(I,J,K)       &
	     &                     * (DV(I,J,K)-DV(I,J,K-1))
		    ENDDO
		  ENDDO
		ENDDO ! bl_levels
	      Endif

	      If ( l_ftl ) Then
		DO K=2,BL_LEVELS
		  DO J=1,ROWS
		    DO I=1,ROW_LENGTH
		      FTL(I,J,K) = FTL(I,J,K) - GAMMA(K)*RHOKH(I,J,K)           &
	     &                   * RDZ_CHARNEY_GRID(I,J,K)                      &
	     &                   * (DTL(I,J,K)-DTL(I,J,K-1))
		    ENDDO
		  ENDDO
		ENDDO
	      Endif

	      If ( l_fqw ) Then
		DO K=2,BL_LEVELS
		  DO J=1,ROWS
		    DO I=1,ROW_LENGTH
		      FQW(I,J,K) = FQW(I,J,K) - GAMMA(K)*RHOKH(I,J,K)           &
	     &                   * RDZ_CHARNEY_GRID(I,J,K)                      &
	     &                   * (DQW(I,J,K)-DQW(I,J,K-1))
		    ENDDO
		  ENDDO
		ENDDO ! bl_levels
	      Endif



	      IF (LTIMER) THEN
	! DEPENDS ON: timer
		CALL TIMER('IM_BL_PT2 ',4)
	      ENDIF

	      RETURN
	      END SUBROUTINE IM_BL_PT2
