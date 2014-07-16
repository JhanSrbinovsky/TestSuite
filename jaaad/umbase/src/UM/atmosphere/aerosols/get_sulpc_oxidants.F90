#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE GET_SULPC_OXIDANTS(                                    &

! Arguments IN
     &  L_SULPC_ONLINE_OXIDANTS, L_UKCA, L_UKCA_TROP, L_UKCA_TROPISOP,  &
     &  L_UKCA_STRATTROP, first_atmstep_call,                           &
#include "arg_atm_fields.h"
#include "argd1.h"
! Arguments OUT
     &  O3_MMR_out, H2O2_MMR_out, OH_conc_out, HO2_conc_out)
!
!---------------------------------------------------------------------
! Purpose: To extract from the D1 array the appropriate concentrations 
!          of oxidants for the sulphur cycle, either prescribed 
!          oxidants from ancillary file or on-line oxidants from UKCA, 
!          depending on the value of L_SULPC_ONLINE_OXIDANTS.
!          Called by Atm_Step.
!
! Current owner of code:   Jamie Rae
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation:
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE

! Global variables:
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "csubmodl.h"

! Variables with intent IN:
#include "ctracera.h"
#include "typ_atm_fields.h"
      LOGICAL, INTENT(IN) :: first_atmstep_call
      LOGICAL, INTENT(IN) :: L_SULPC_ONLINE_OXIDANTS
      LOGICAL, INTENT(IN) :: L_UKCA
      LOGICAL, INTENT(IN) :: L_UKCA_TROP
      LOGICAL, INTENT(IN) :: L_UKCA_TROPISOP
      LOGICAL, INTENT(IN) :: L_UKCA_STRATTROP

! Variables with intent OUT:
! Oxidant arrays to be passed back to ATM_STEP have no haloes:
      REAL, INTENT(OUT) :: O3_MMR_out(row_length,rows,model_levels)
      REAL, INTENT(OUT) :: H2O2_MMR_out(row_length,rows,model_levels)
      REAL, INTENT(OUT) :: OH_conc_out(row_length,rows,model_levels)
      REAL, INTENT(OUT) :: HO2_conc_out(row_length,rows,model_levels)

! Local scalars:
      INTEGER :: m_atm_modl ! Atmosphere model identifier
      INTEGER :: section    ! Section number of diagnostic
      INTEGER :: item       ! Item number of diagnostic
      INTEGER :: halo_typ   ! Halo type of diagnostic
      INTEGER :: error_code ! For passing to subroutine EREPORT
      INTEGER :: array_size 
      INTEGER, SAVE :: x_halo_OH
      INTEGER, SAVE :: x_halo_O3
      INTEGER, SAVE :: x_halo_H2O2
      INTEGER, SAVE :: x_halo_HO2
      INTEGER, SAVE :: y_halo_OH
      INTEGER, SAVE :: y_halo_O3
      INTEGER, SAVE :: y_halo_H2O2
      INTEGER, SAVE :: y_halo_HO2
      INTEGER :: i          ! Loop variable
      INTEGER :: j          ! Loop variable
      INTEGER :: k          ! Loop variable
      INTEGER, PARAMETER :: ukca_sect = 33  ! UKCA STASH section.

#include "c_v_m.h"
#include "c_sulchm.h"
#include "c_rmol.h"

! Local arrays:
      INTEGER, PARAMETER :: ukca_item(4) = (/1,8,251,252/)
!                                 STASH item numbers for UKCA oxidants
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: O3_MMR
!                     O3 mass mixing ratio (kg[O3]/kg[air])
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: H2O2_MMR
!                     H2O2 mass mixing ratio (kg[H2O2]/kg[air])
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: OH_MMR
!                     OH mass mixing ratio (kg[OH]/kg[air])
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: HO2_MMR
!                     HO2 mass mixing ratio (kg[HO2]/kg[air])
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: OH_conc
!                     OH concentration (molecules/cm3)
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: HO2_conc
!                     HO2 concentration (molecules/cm3)

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_theta_levels_local
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_theta_levels_local
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: theta_local 
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: T
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: air_density  
!                                                     molecules / cm3

      CHARACTER*72 :: cmessage        ! Error message

      LOGICAL, DIMENSION(4) :: L_UKCA_FOUND ! Denotes whether or not each 
                                            ! UKCA oxidant has been found 
                                            ! in D1 array
!----------------------------------------------------------------------

! Check that values of logical variables are consistent:
      IF((L_UKCA_TROP .OR. L_UKCA_TROPISOP .OR. L_UKCA_STRATTROP)     &
     &       .AND. (.NOT. L_UKCA)) THEN
        cmessage = 'Inconsistency in UKCA logicals.'
        CALL EREPORT('GET_SULPC_OXIDANTS',100,cmessage)
      END IF

      IF(L_SULPC_ONLINE_OXIDANTS .AND. (.NOT. (L_UKCA_TROP .OR.       &
     &         L_UKCA_TROPISOP .OR. L_UKCA_STRATTROP))) THEN
        cmessage = 'UKCA oxidants are selected for sulphur cycle, '   &
     &              //'but UKCA is turned off.'
        CALL EREPORT('GET_SULPC_OXIDANTS',101,cmessage)
      END IF

! Determine halo sizes for UKCA oxidant arrays:
      IF(L_SULPC_ONLINE_OXIDANTS) THEN
        IF(first_atmstep_call) THEN
          L_UKCA_FOUND(:) = .FALSE.             ! Initalise to F
          m_atm_modl = SUBMODEL_FOR_SM(a_im)
          DO i=1,no_obj_D1(m_atm_modl)
            section = D1_ADDR(D1_section,i,m_atm_modl)
            IF(section == UKCA_sect) THEN
              item = D1_ADDR(D1_item,i,m_atm_modl)
              IF(item == ukca_item(1)) THEN
                halo_typ = D1_ADDR(D1_halo_type,i,m_atm_modl) 
                x_halo_O3 = HALOSIZE(1,halo_typ)
                y_halo_O3 = HALOSIZE(2,halo_typ)
                L_UKCA_FOUND(1) = .TRUE.
              ELSE IF(item == ukca_item(2)) THEN
                halo_typ = D1_ADDR(D1_halo_type,i,m_atm_modl) 
                x_halo_H2O2 = HALOSIZE(1,halo_typ)
                y_halo_H2O2 = HALOSIZE(2,halo_typ)
                L_UKCA_FOUND(2) = .TRUE.
              ELSE IF(item ==  ukca_item(3)) THEN
                halo_typ = D1_ADDR(D1_halo_type,i,m_atm_modl) 
                x_halo_OH = HALOSIZE(1,halo_typ)
                y_halo_OH = HALOSIZE(2,halo_typ)
                L_UKCA_FOUND(3) = .TRUE.
              ELSE IF(item ==  ukca_item(4)) THEN
                halo_typ = D1_ADDR(D1_halo_type,i,m_atm_modl)
                x_halo_HO2 = HALOSIZE(1,halo_typ)
                y_halo_HO2 = HALOSIZE(2,halo_typ)
                L_UKCA_FOUND(4) = .TRUE.
              END IF ! IF condition on item
            END IF   ! IF condition on section
          END DO

          DO i=1,4
            IF(.NOT. L_UKCA_FOUND(i)) THEN 
! If any UKCA oxidants cannot be found in D1 array, return an error 
! message and stop.
              cmessage = 'At least one UKCA oxidant not found in D1 '   &
     &                     //'array.'
              error_code = (ukca_sect * 1000) + ukca_item(i)
              CALL EREPORT('GET_SULPC_OXIDANTS',error_code,cmessage)
            END IF
          END DO
        END IF ! IF(first_atmstep_call)

! Extract UKCA oxidants from D1 array:


        ALLOCATE(O3_MMR(1-x_halo_O3:row_length+x_halo_O3,               &
     &                     1-y_halo_O3:rows+y_halo_O3, 1:model_levels))
        array_size = SIZE(O3_MMR)
        O3_MMR = RESHAPE(O3_UKCA,                                       &
     &          (/ row_length+(2*x_halo_O3), rows+(2*y_halo_O3),        &
     &            model_levels /))
        ALLOCATE(H2O2_MMR(1-x_halo_H2O2:row_length+x_halo_H2O2,         &
     &              1-y_halo_H2O2:rows+y_halo_H2O2, 1:model_levels))
        array_size = SIZE(H2O2_MMR)
        H2O2_MMR = RESHAPE(                                             &
     &                 H2O2_UKCA,                                       &
     &                 (/ row_length+(2*x_halo_H2O2),                   &
     &                  rows+(2*y_halo_H2O2),model_levels /))     
        ALLOCATE(OH_MMR(1-x_halo_OH:row_length+x_halo_OH,               &
     &                    1-y_halo_OH:rows+y_halo_OH, 1:model_levels))
        array_size = SIZE(OH_MMR)
        OH_MMR = RESHAPE(OH_UKCA,                                       &
     &            (/ row_length+(2*x_halo_OH), rows+(2*y_halo_OH),      &
     &            model_levels /) )      
        ALLOCATE(HO2_MMR(1-x_halo_HO2:row_length+x_halo_HO2,            &
     &                 1-y_halo_HO2:rows+y_halo_HO2, 1:model_levels))
        array_size = SIZE(HO2_MMR)
        HO2_MMR = RESHAPE(HO2_UKCA,                                     &
     &            (/ row_length+(2*x_halo_HO2), rows+(2*y_halo_HO2),    &
     &             model_levels /))

! Extract from D1 the arrays necessary for calculation of air density:
        ALLOCATE(p_theta_levels_local(1-offx:row_length+offx,                 &
     &            1-offy:rows+offy, model_levels))
        array_size = SIZE(p_theta_levels_local)
        p_theta_levels_local = RESHAPE(p_theta_levels,                  &
     &                      (/ row_length+(2*offx),rows+(2*offy),       &
     &                       model_levels /) )
        ALLOCATE(theta_local(1-offx:row_length+offx,1-offy:rows+offy,         &
     &              model_levels))
        array_size = SIZE(theta_local)
        theta_local = RESHAPE(theta,                                    &
     &      (/row_length+(2*offx), rows+(2*offy), model_levels/))     
        ALLOCATE(exner_theta_levels_local(1-offx:row_length+offx,             &
     &          1-offy:rows+offy, model_levels))
        array_size = SIZE(exner_theta_levels_local)
        exner_theta_levels_local = RESHAPE(exner_theta_levels,          &
     &         (/ row_length+(2*offx), rows+(2*offy), model_levels /) )
        ALLOCATE(OH_conc(row_length,rows,model_levels))
        ALLOCATE(HO2_conc(row_length,rows,model_levels))  
        ALLOCATE(T(row_length,rows,model_levels)) 
        ALLOCATE(air_density(row_length,rows,model_levels)) 

! Convert OH and HO2 from mass-mixing ratios to molecules/cm3, as 
! required for sulphur-cycle routine:
        DO k=1,model_levels                                  
          DO j=1,rows                                      
            DO i=1,row_length        
              T(i,j,k)=theta_local(i,j,k)*exner_theta_levels_local(i,j,k) 
              air_density(i,j,k) = p_theta_levels_local(i,j,k) * avogadro     &
     &             * 1.0E-06 / (rmol * T(i,j,k))               
              OH_conc(i,j,k) = OH_MMR(i,j,k) * air_density(i,j,k)       &
     &                             / c_OH               
              HO2_conc(i,j,k) = HO2_MMR(i,j,k) * air_density(i,j,k)     &
     &                           / c_HO2            
            END DO                                        
          END DO                                           
        END DO 
        DEALLOCATE(p_theta_levels_local)
        DEALLOCATE(theta_local)
        DEALLOCATE(exner_theta_levels_local)
        DEALLOCATE(T)
        DEALLOCATE(air_density)             
      ELSE

! Prescribed oxidants from ancillary selected. Extract prescribed 
! oxidants from D1 array.
        ALLOCATE(O3_MMR(row_length, rows, model_levels))
        ALLOCATE(H2O2_MMR(row_length, rows, model_levels))
        ALLOCATE(OH_conc(row_length, rows, model_levels))
        ALLOCATE(HO2_conc(row_length, rows, model_levels))
        array_size = SIZE(O3_MMR)
        O3_MMR = RESHAPE(O3_chem,                                     &
     &               (/ row_length, rows, model_levels /) )
        H2O2_MMR = RESHAPE(H2O2_limit,(/ row_length, rows, model_levels/) )
        OH_conc = RESHAPE(OH, (/ row_length, rows, model_levels /) )
        HO2_conc = RESHAPE(HO2, (/ row_length, rows, model_levels /) )
      END IF

      O3_MMR_out(:,:,:) = O3_MMR(1:row_length,1:rows,1:model_levels)
      H2O2_MMR_out(:,:,:) = H2O2_MMR(1:row_length,1:rows,1:model_levels)
      OH_conc_out(:,:,:) = OH_conc(1:row_length,1:rows,1:model_levels) 
      HO2_conc_out(:,:,:) = HO2_conc(1:row_length,1:rows,1:model_levels)

      IF(ALLOCATED(O3_MMR)) DEALLOCATE (O3_MMR)
      IF(ALLOCATED(H2O2_MMR)) DEALLOCATE (H2O2_MMR)
      IF(ALLOCATED(OH_MMR)) DEALLOCATE (OH_MMR)
      IF(ALLOCATED(HO2_MMR)) DEALLOCATE (HO2_MMR)
      IF(ALLOCATED(OH_conc)) DEALLOCATE (OH_conc)
      IF(ALLOCATED(HO2_conc)) DEALLOCATE (HO2_conc)

      END SUBROUTINE GET_SULPC_OXIDANTS

#endif
