#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper routine for ISCCP_CLOUDTYPES
!
! Purpose:
!   To loop over lit columns, calling ISCCP_CLOUD_TYPES in each
!
! Method:
!   Trivial
!
! Current Owner of Code: M. J. Webb
!
! History:
! Version       Date             Comment
!  4.4(as mod)  20-09-99         Original Code based on swrad3a.f
!                                (M. J. Webb)
!  5.5          21-06-02         Upgraded to work with new dynamics.
!                                (A.Keen / K.Williams)
!  6.1          30-09-03         Replaced with vector version for SX6
!                                (K.Williams / M.Webb)
!  6.1          05-05-04         Make mean cloud optical depth, top
!                                pressure and total cloud area
!                                diagnostics available.
!                                (A. Jones)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE ISCCP_CLOUDTYPES_FLD(I_SEGMENT                         &
!          Input
     &   , NPD_FIELD                                                    &
                        ! Number of points in field
     &   , N_POINTS                                                     &
                       ! Number of points in segment
     &   , NLIT                                                         &
                        ! Number of lit points
     &   , LIST                                                         &
                        ! Indexes of lit points
     &   , NLEV                                                         &
                        ! Number of levels
     &   , NCOL                                                         &
                        ! Number of columns for gridbox decomposition
     &   , PFULL_FLD                                                    &
                        ! Pressure on full levels
     &   , PHALF_FLD                                                    &
                        ! Pressure on full levels
     &   , QV_FLD                                                       &
                        ! Water vapour
     &   , CC_FLD                                                       &
                        ! Total cloud amount (cc+(1-cc)*ls)
     &   , CONV_FLD                                                     &
                        ! Convective cloud amount
     &   , DTAU_S_FLD                                                   &
                        ! Large-scale optical thickness
     &   , DTAU_C_FLD                                                   &
                        ! Convective optical thickness
     &   , TOP_HEIGHT                                                   &
                        ! Flag to specify cloud top method
     &   , OVERLAP                                                      &
                        ! Flag to specify overlap method
     &   , SKT_FLD                                                      &
                        ! Surface temperature
     &   , EMSFC_LW_FLD                                                 &
                        ! Surface emissivity
     &   , AT_FLD                                                       &
                        ! Atmospheric temperature
     &   , DEM_S_FLD                                                    &
                        ! Large-scale cloud emissivities
     &   , DEM_C_FLD                                                    &
                        ! Convective cloud emissivities
     &   , TRINDX_FLD                                                   &
                        ! Tropopause index
!          Output
     &   , FQ_ISCCP_FLD                                                 &
                        ! Cloud amounts in various ISCCP classes
     &   , MEANALBEDOCLD                                                &
                        ! Weighted mean cloud albedo
     &   , MEANTAUCLD                                                   &
                        ! Weighted mean cloud optical depth
     &   , MEANPTOP                                                     &
                        ! Weighted mean cloud top pressure
     &   , TOTALCLDAREA                                                 &
                        ! Total cloud area
     &   )
!
      IMPLICIT NONE
!
!     COMDECKS INCLUDED

#include "isccpdata.h"

!
!     ARGUMENTS

!     Input
      INTEGER DEBUG
      INTEGER DEBUGCOL
      INTEGER I_SEGMENT
      INTEGER NPD_FIELD
      INTEGER N_POINTS
      INTEGER NLIT
      INTEGER LIST(NPD_FIELD)
      INTEGER NLEV
      INTEGER NCOL
      REAL PFULL_FLD(NPD_FIELD, NLEV)
      REAL PHALF_FLD(NPD_FIELD, NLEV+1)
      REAL QV_FLD(NPD_FIELD, NLEV)
      REAL CC_FLD(NPD_FIELD, NLEV)
      REAL CONV_FLD(NPD_FIELD, NLEV)
      REAL DTAU_S_FLD(NPD_FIELD, NLEV)
      REAL DTAU_C_FLD(NPD_FIELD, NLEV)
      INTEGER TOP_HEIGHT
      !  1 = adjust top height, that is compute infrared
      !  brightness temperature and adjust cloud top
      !  pressure accordingly
      !  2 = do not adjust top height, that is cloud top
      !  pressure is the actual cloud top pressure
      !  in the model
      INTEGER OVERLAP
      !  overlap type
      !  1=max
      !  2=rand
      !  3=max/rand
      REAL SKT_FLD(NPD_FIELD)
      REAL EMSFC_LW_FLD(NPD_FIELD)
      REAL AT_FLD(NPD_FIELD, NLEV)
      REAL DEM_S_FLD(NPD_FIELD, NLEV)
      REAL DEM_C_FLD(NPD_FIELD, NLEV)
      INTEGER TRINDX_FLD(NPD_FIELD)
!     Output
      REAL FQ_ISCCP_FLD(NPD_FIELD, 7, 7)
      REAL MEANALBEDOCLD(NPD_FIELD)
      REAL MEANTAUCLD(NPD_FIELD)
      REAL MEANPTOP(NPD_FIELD)
      REAL TOTALCLDAREA(NPD_FIELD)

!     LOCAL VARIABLES.
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE

      REAL PFULL_FLD_LIT(NLIT, NLEV)
      REAL PHALF_FLD_LIT(NLIT, NLEV+1)
      REAL QV_FLD_LIT(NLIT, NLEV)
      REAL CC_FLD_LIT(NLIT, NLEV)
      REAL CONV_FLD_LIT(NLIT, NLEV)
      REAL DTAU_S_FLD_LIT(NLIT, NLEV)
      REAL DTAU_C_FLD_LIT(NLIT, NLEV)
      REAL AT_FLD_LIT(NLIT, NLEV)
      REAL DEM_S_FLD_LIT(NLIT, NLEV)
      REAL DEM_C_FLD_LIT(NLIT, NLEV)
      INTEGER SEED_FLD_LIT(NLIT)
      REAL SKT_FLD_LIT(NLIT)
      REAL EMSFC_LW_FLD_LIT(NLIT)
      INTEGER TRINDX_FLD_LIT(NLIT)
      INTEGER SUNLIT_FLD(NLIT)
      INTEGER SEED_FLD(NLIT)
      REAL FQ_ISCCP_FLD_LIT(NLIT, 7, 7)
      REAL MEANALBEDOCLD_LIT(NLIT)
      REAL MEANTAUCLD_LIT(NLIT)
      REAL MEANPTOP_LIT(NLIT)
      REAL TOTALCLDAREA_LIT(NLIT)

!     Initialise output arrays to zero
      DO I=1,7
        DO J=1,7
          DO L=1,N_POINTS
            FQ_ISCCP_FLD(L,I,J)=0
          ENDDO
        ENDDO
      ENDDO
      DO L=1,N_POINTS
        MEANALBEDOCLD(L)=0.0
        MEANTAUCLD(L)=0.0
        MEANPTOP(L)=0.0
        TOTALCLDAREA(L)=0.0
      ENDDO

      DO L=1,NLIT
        DO I=1,NLEV+1
          PHALF_FLD_LIT(L,NLEV+1+1-I)=PHALF_FLD(LIST(L),I)
        ENDDO
        DO I=1,NLEV
          PFULL_FLD_LIT(L,NLEV+1-I)=PFULL_FLD(LIST(L),I)
          QV_FLD_LIT(L,NLEV+1-I)=QV_FLD(LIST(L),I)
          CC_FLD_LIT(L,NLEV+1-I)=CC_FLD(LIST(L),I)
          CONV_FLD_LIT(L,NLEV+1-I)=CONV_FLD(LIST(L),I)
          DTAU_S_FLD_LIT(L,NLEV+1-I)=DTAU_S_FLD(LIST(L),I)
          DTAU_C_FLD_LIT(L,NLEV+1-I)=DTAU_C_FLD(LIST(L),I)
          AT_FLD_LIT(L,NLEV+1-I)=AT_FLD(LIST(L),I)
          DEM_S_FLD_LIT(L,NLEV+1-I)=DEM_S_FLD(LIST(L),I)
          DEM_C_FLD_LIT(L,NLEV+1-I)=DEM_C_FLD(LIST(L),I)
        ENDDO
      ENDDO

!     change tropopause index to count the other way up
      DO L=1,NLIT
        IF (TRINDX_FLD(L) >  0) THEN
          TRINDX_FLD_LIT(L)=NLEV+1-TRINDX_FLD(LIST(L))
        ELSE
          TRINDX_FLD_LIT(L)=-1
        ENDIF
        SKT_FLD_LIT(L)=SKT_FLD(LIST(L))
        EMSFC_LW_FLD_LIT(L)=EMSFC_LW_FLD(LIST(L))
        SEED_FLD_LIT(L)=(PFULL_FLD_LIT(L,nlev)                          &
     &  -int(PFULL_FLD_LIT(L,nlev)))*1000000
        SUNLIT_FLD(L)=1
      ENDDO

            write (6,*)                                                 &
     & 'isccp_cloud_types_fld-3.4: Calling isccp_cloudtypes: '

      debug=0
      debugcol=0
! DEPENDS ON: isccp_cloud_types
      CALL ISCCP_CLOUD_TYPES(                                           &
     &           debug                                                  &
     &          ,debugcol                                               &
     &          ,nlit                                                   &
     &          ,sunlit_fld                                             &
     &          ,nlev                                                   &
     &          ,ncol                                                   &
     &          ,seed_fld_lit                                           &
     &          ,pfull_fld_lit                                          &
     &          ,phalf_fld_lit                                          &
     &          ,qv_fld_lit                                             &
     &          ,cc_fld_lit                                             &
     &          ,conv_fld_lit                                           &
     &          ,dtau_s_fld_lit                                         &
     &          ,dtau_c_fld_lit                                         &
     &          ,top_height                                             &
     &          ,overlap                                                &
     &          ,tautab                                                 &
     &          ,invtau                                                 &
     &          ,skt_fld_lit                                            &
     &          ,emsfc_lw_fld_lit                                       &
     &          ,at_fld_lit                                             &
     &          ,dem_s_fld_lit                                          &
     &          ,dem_c_fld_lit                                          &
     &          ,fq_isccp_fld_lit                                       &
     &          ,trindx_fld_lit                                         &
     &          ,totalcldarea_lit                                       &
     &          ,meanptop_lit                                           &
     &          ,meantaucld_lit                                         &
     &          ,meanalbedocld_lit                                      &
!     &          ,boxtau
!     &          ,boxptop
     & )


!      invert output arrays
      DO L=1,NLIT
       MEANALBEDOCLD(LIST(L))=MEANALBEDOCLD_LIT(L)
       MEANTAUCLD(LIST(L))=MEANTAUCLD_LIT(L)
       MEANPTOP(LIST(L))=MEANPTOP_LIT(L)
       TOTALCLDAREA(LIST(L))=TOTALCLDAREA_LIT(L)
       DO I=1,7
          DO J=1,7
            FQ_ISCCP_FLD(LIST(L),J,7+1-I)=FQ_ISCCP_FLD_LIT(L,J,I)
          ENDDO
        ENDDO
      ENDDO

!     sychronise PEs
!     call barrier()
!     Force core dump
!      if (node == -1) nlit=sqrt(real((-nlit-1)))
!     Force other PE's stop here
!      call barrier()

      RETURN
      END SUBROUTINE ISCCP_CLOUDTYPES_FLD
#endif
