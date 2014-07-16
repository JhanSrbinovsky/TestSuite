#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!     Budget calculation for UKCA gas-phase chemistry.
!     Contains FUNCTION elm.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern
!                            Fiona O'Connor
!
! Method:
!   This update contains code necessary for budget calculations. We
!   calculate throughput through every reaction channel, wet and dry
!   deposition and emission.
!
!   For emission and dry deposition, the full surface fields of of the
!   amount of tracer deposited or emitted per gridcell and per timestep
!   is calculated in units of mol. There is one 3-D field for dry
!   deposition and one for emission. Integrated totals are achieved by
!   selecting "Accumulation every 1 timestep" in the STASH profile.

!   For bi- and termolecular reactions, photolysis and wet deposition
!   zonal totals are calculated. (Full 3-D fields would take too much
!   space but the user is free to add 3-D budget code for specific
!   reactions. Two 3-D fields are assigned for bimolecular reactions
!   and 1 each for termolecular, photolysis and wet deposition reactions.
!
!   The budgets are calculated as (tendency at end of timestep due to
!   a specific reaction) * volume * timestep, summed up over each row
!   and expressed as moles. They are approximate as they do not take
!   into account nonlinearities and coupling between reactions which
!   ASAD would consider. ASAD variables are evaluated after chemistry.
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE ukca_budget2(  rows, row_length, n_pnts, k, p_levelsda,&
                      secs_per_step,                                    &
                      zdryrt, zwetrt, zprt,                             &
                      nabund, nbimol, nphotf,                           &
                      abundance,                                        &
                      bimol_budget,                                     &
                      termol_budget,                                    &
                      hetero_budget, dd_budget,                         &
                      wd_budget,                                        &
                      phot_budget, volume, pres,                        &
                      n_cpl, n_o3,                                      &
                      cpl_o3,                                           &
                      global_rows, global_row_length,                   &
                      first_row, first_column,                          &
                      within_pe_domain )

      USE ASAD_MOD
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "c_v_m.h"

! parallel variables
#include "parvars.h"
! sulphur chemistry constants
#include "c_sulchm.h"

! Subroutine interface

      INTEGER, INTENT(IN) :: rows            ! number of rows
      INTEGER, INTENT(IN) :: row_length      ! length of row
      INTEGER, INTENT(IN) :: n_pnts          ! length of chem. vector
      INTEGER, INTENT(IN) :: k               ! level of budgeting
      INTEGER, INTENT(IN) :: p_levelsda      ! number of levels
      INTEGER, INTENT(IN) :: nabund          ! number of 3-D fields for
                                             ! abundace
      INTEGER, INTENT(IN) :: nbimol          ! number of 3-D fields for
      INTEGER, INTENT(IN) :: nphotf          ! bimol and phot. reactions

      INTEGER, INTENT(IN) :: n_o3            ! number of O3 variable in
                                             ! species field
      INTEGER, INTENT(IN) :: n_cpl           ! number of chemical production
                                             ! and loss fields for O3

! variables needed for performing the zonal total
      INTEGER, INTENT(IN) :: global_rows
      INTEGER, INTENT(IN) :: global_row_length
      INTEGER, INTENT(IN) :: first_row
      INTEGER, INTENT(IN) :: first_column

      LOGICAL, INTENT(IN) :: within_pe_domain

      REAL, INTENT(IN) :: secs_per_step      ! timestep

! Input fields
      REAL, INTENT(IN) :: zdryrt(row_length, rows, jpdd)
      REAL, INTENT(IN) :: zwetrt(row_length, rows, p_levelsda, jpdw)
      REAL, INTENT(IN) :: zprt(row_length, rows, jppj)
      REAL, INTENT(IN) :: volume(row_length, rows, p_levelsda)
      REAL, INTENT(IN) :: pres(row_length, rows, p_levelsda)

! Output budget fields
      REAL, INTENT(INOUT) ::                                            &
                      bimol_budget(row_length, rows, p_levelsda, nbimol)
      REAL, INTENT(INOUT) :: termol_budget(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: hetero_budget(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: wd_budget(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) :: dd_budget(row_length, rows, p_levelsda)
      REAL, INTENT(INOUT) ::                                            &
                       phot_budget(row_length, rows, p_levelsda, nphotf)
      REAL, INTENT(INOUT) :: cpl_o3(row_length, rows, p_levelsda, n_cpl)
      REAL, INTENT(INOUT) ::                                            &
                         abundance(row_length, rows, p_levelsda, nabund)

! Local variables
      INTEGER :: budget_size

! Conversion factor, to revert from molecules/cm^3 to moles for throughput
      REAL, PARAMETER :: convfac = 1.0E6 / avogadro

      INTEGER :: icode

      INTEGER :: i
      INTEGER :: j
      INTEGER :: l
      INTEGER :: l1
      INTEGER :: l2
      INTEGER :: m
      INTEGER :: m1
      INTEGER :: m2
      INTEGER :: m3
      INTEGER :: index
      INTEGER :: p1
      INTEGER :: p2

      LOGICAL, SAVE :: firstcall = .TRUE.

      REAL :: reaction_throughput(row_length)

! These transferred to ASAD_MOD because of save attribute
!      INTEGER, SAVE :: rpartner_bimol(jpbk, 2)
!      INTEGER, SAVE :: rpartner_termol(jptk, 2)
! set to zero to cope with jphk=0 condition
!      INTEGER, SAVE :: rpartner_hetero(0:jphk, 2)
!      INTEGER, SAVE :: rpartner_wd(jpdw)
!      INTEGER, SAVE :: rpartner_dd(jpdd)
!      INTEGER, SAVE :: rpartner_phot(jppj)
!      REAL, ALLOCATABLE, SAVE :: zonal_total(:,:,:)

      budget_size = jpnr + jpspec + jpdw
! main block
      IF (firstcall) THEN ! initializing block

        WRITE(6,*) 'FIRST_ROW  =',first_row, first_column

        rpartner_bimol = -1
        rpartner_termol = -1
        IF (jphk > 0) THEN
          rpartner_hetero = -1
          rpartner_hetero(0,:) = 0
        ENDIF
        IF (jpdw > 0) rpartner_wd = -1
        IF (jpdd > 0) rpartner_dd = -1
        rpartner_phot = -1

! calculate reaction pairing information
        DO m=1,2      ! count over reaction partners for reactions

          DO j=1,jpspec
            DO i=1,jpbk ! count over reactions
! find species number corresponding to m-th reaction partner
              IF (speci(j) == spb(i,m)) rpartner_bimol(i,m) = j
            END DO

            DO i=1,jptk
              IF (speci(j) == spt(i,m)) rpartner_termol(i,m) = j
            END DO

            IF (jphk > 0) THEN
              DO i=1,jphk
                IF (speci(j) == sph(i,m)) rpartner_hetero(i,m) = j
              END DO
            END IF
          END DO

        END DO

! find species number corresponding to photolysis reactions
        DO i=1,jppj ! count over reactions
          DO j=1,jpspec
            IF (speci(j) == spj(i,1)) rpartner_phot(i) = j
          END DO
        END DO

! find species number corresponding to dry deposition reaction
        i=1
        DO j=1,jpspec
          IF (ldepd(j)) THEN
            rpartner_dd(i) = j
            i = i + 1
          END IF
        END DO

! find species number corresponding to wet deposition reaction
        IF (jpdw > 0) THEN
          i=1
          DO j=1,jpspec
            IF (ldepw(j)) THEN
              rpartner_wd(i) = j
              i = i + 1
            END IF
          END DO
        END IF

        IF (minval(rpartner_bimol) < 0) THEN
          DO i=1,jpbk
            WRITE (6,*) i,(rpartner_bimol(i,m),m=1,2)
          END DO
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',1,                                     &
           'Undefined bimol. reaction partner')
        END IF

! Don't check 2nd reaction partner of termol. reactions; could be "m"
! which is not in the list of species.
        IF (minval(rpartner_termol(:,1)) < 0)                           &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',2,                                     &
            'Undefined termol. reaction partner')

        IF (minval(rpartner_phot) < 0) THEN
          IF (mype == 0)                                                &
            WRITE (6,*) rpartner_phot
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',3,                                     &
            'Undefined photol. reaction partner')
        END IF

        IF (jphk > 0) THEN
          IF (minval(rpartner_hetero) < 0)                              &
! DEPENDS ON: ereport
            CALL ereport('BUDGET2',4,                                   &
              'Undefined termol. reaction partner')
        END IF

        IF (minval(rpartner_dd) < 0) THEN
          IF (mype == 0) THEN
            WRITE (6,*) (rpartner_dd(j),j=1,jpdd)
          END IF
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',5,                                     &
            'jpdd does not match chch.d')
        END IF

        IF (jpdw > 0) THEN
          IF (minval(rpartner_wd) < 0)                                  &
! DEPENDS ON: ereport
            CALL ereport('BUDGET2',6,                                   &
              'jpdw does not match chch.d')
        END IF

        IF (jpbk > nbimol*global_row_length)                            &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',7,                                     &
            'Increase size of bimol_budget')

        IF (jppj > nphotf*global_row_length)                            &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',8,                                     &
            'Increase size of phot_budget')

        IF (jptk > global_row_length)                                   &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',9,                                     &
            'Increase size of termol_budget')

        IF (jphk > global_row_length)                                   &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',10,                                    &
            'Increase size of hetero_budget')

        IF (jpspec > global_row_length*nabund)                          &
! DEPENDS ON: ereport
          CALL ereport('BUDGET2',11,                                    &
            'Increase size of abundance')

        firstcall = .FALSE.

      END IF !   end  initialization

! calculate abundances for all species and produce zonal total

      IF (within_pe_domain) THEN
        IF (.NOT. ALLOCATED(zonal_total)) THEN
          ALLOCATE(zonal_total(budget_size, global_rows, p_levelsda))
          zonal_total = 0.
        END IF

        index = 0

        DO m=1,jpspec
          index = index + 1
          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! calculate species abundance. y is in molecules per cm^3. volume is in m^3
! Hence factor 10^6 necessary. Divide by Avogadro's number to get moles.

            reaction_throughput = y(l1:l2,m) * volume(:,i,k)*convfac

! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=SUM(reaction_throughput)
          END DO
        END DO


! calculate 3-D reaction rate for each bimol. reaction and produce
! column integral

        DO m=1,jpbk  ! loop for bimol. reactions
          index = index + 1
          p1 = rpartner_bimol(m,1)
          p2 = rpartner_bimol(m,2)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! reaction turnaround in mol per gridcell per timestep;
! Factor 1E6 necessary because volume is in m^3 but concentration in
! molecules per cm^3. Divide by Avogadro's number to revert to mol.
! The summation is over latitude circles (zonal total).

            reaction_throughput =                                       &
              rk(l1:l2,nbrkx(m))*y(l1:l2,p1)*y(l1:l2,p2)                &
              * secs_per_step * volume(:,i,k) * convfac

! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=SUM(reaction_throughput)
          END DO
        END DO

! calculate throughput for termolecular reactions

        DO m=1,jptk
          index = index + 1
          p1 = rpartner_termol(m,1)
          p2 = rpartner_termol(m,2)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length
            IF (p2 > 0) THEN
! regular termolecular reaction
              reaction_throughput=                                      &
                rk(l1:l2,ntrkx(m))*y(l1:l2,p1)*y(l1:l2,p2)              &
                * secs_per_step * volume(:,i,k) * convfac
            ELSE
! thermal decomposition
              reaction_throughput=rk(l1:l2,ntrkx(m))*y(l1:l2,p1)        &
                  * secs_per_step * volume(:,i,k) * convfac
            END IF
 ! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=SUM(reaction_throughput)
          END DO
        END DO

! calculate throughput for photolysis reactions

        DO m=1,jppj
          index = index + 1
          p1 = rpartner_phot(m)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! reaction rate in molecules per gridcell per timestep;
! Factor 1E6 necessary because volume is in m^3 but concentration in
! moelcules per cm^3.

            reaction_throughput=zprt(:,i,m) * y(l1:l2,p1)               &
              * secs_per_step * volume(:,i,k) * convfac

! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=SUM(reaction_throughput)
          END DO
        END DO

! calculate throughput for wet deposition reactions

        DO m=1,jpdw
          index = index + 1
          p1 = rpartner_wd(m)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! reaction rate in molecules per gridcell per timestep;
! Factor 1E6 necessary because volume is in m^3 but concentration in
! moelcules per cm^3.

            reaction_throughput=zwetrt(:,i,k,m) * y(l1:l2,p1)           &
                * secs_per_step * volume(:,i,k) * convfac

! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=                         &
              SUM(reaction_throughput)
          END DO
        END DO
! calculate throughput for dry deposition.
! Here: Stack different species with dry deposition
! into different levels of budget variable. DO NOT TAKE ZONAL MEAN.
        DO m=1,jpdd
          p1 = rpartner_dd(m)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! reaction rate in molecules per gridcell per timestep;
! Factor 1E6 necessary because volume is in m^3 but concentration in
! moelcules per cm^3.

            dd_budget(:,i,m)=zdryrt(:,i,m) * y(l1:l2,p1)                &
              * secs_per_step * volume(:,i,k) * convfac

          END DO
        END DO

! calculate throughput for heterogeneous reactions

        DO m=1,jphk
          index = index + 1
          p1 = rpartner_hetero(m,1)
          p2 = rpartner_hetero(m,2)

          DO i=1,rows
            l1 = (i-1) * row_length + 1
            l2 = i * row_length

! reaction rate in molecules per gridcell per timestep;
! Factor 1E6 necessary because volume is in m^3 but concentration in
! moelcules per cm^3.
            reaction_throughput=                                        &
              rk(l1:l2,nhrkx(m))*y(l1:l2,p1)*y(l1:l2,p2)                &
                * secs_per_step * volume(:,i,k) * convfac

! do sum across latitude circles within processor domain
            zonal_total(index,i+first_row-1,k)=                         &
              SUM(reaction_throughput)
          END DO
        END DO

      ELSE

! End of zonal total taking for each call to budget.
! Next is taking of totals across processor boundaries, done once per
! chemical timestep
        CALL GC_RSUM(budget_size*global_rows*p_levelsda,nproc,          &
                     icode,zonal_total)
        index = 0
        DO m=1,jpspec
          index = index + 1
          m1 = (m-1)/global_row_length + 1 ! field 1 or 2?
          m2 = m - (m1 - 1) * global_row_length ! index within global field
          m3 = m2 - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            abundance(m3,:,:,m1) = zonal_total(index,                   &
                                    first_row:first_row+rows-1,:)
        END DO
        DO m=1,jpbk
          index = index + 1
          m1 = (m-1)/global_row_length + 1 ! field 1 or 2?
          m2 = m - (m1 - 1) * global_row_length ! index within global field
          m3 = m2 - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            bimol_budget(m3,:,:,m1) = zonal_total(index,                &
                                    first_row:first_row+rows-1,:)
        END DO
        DO m=1,jptk
          index = index + 1
          m3 = m - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            termol_budget(m3,:,:) = zonal_total(index,                  &
                                    first_row:first_row+rows-1,:)
        END DO
        DO m=1,jppj
          index = index + 1
          m1 = (m-1)/global_row_length + 1 ! field 1 or 2?
          m2 = m - (m1 - 1) * global_row_length ! index within global field
          m3 = m2 - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            phot_budget(m3,:,:,m1) = zonal_total(index,                 &
                                    first_row:first_row+rows-1,:)
        END DO
        DO m=1,jpdw
          index = index + 1
          m3 = m - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            wd_budget(m3,:,:) = zonal_total(index,                      &
                                    first_row:first_row+rows-1,:)
        END DO
        DO m=1,jphk
          index = index + 1
          m3 = m - first_column + 1       ! index locally
          IF ((m3 >= 1) .AND. (m3 <= row_length))                       &
            hetero_budget(m3,:,:) = zonal_total(index,                  &
                                    first_row:first_row+rows-1,:)
        END DO
        DEALLOCATE(zonal_total)
      END IF

! end of 2-D diagnostics of all chemical channels.
! Next is 3-D diagnostics (full fields) for both tropospheric and strato-
! spheric reactions. This is to allow construction of an ozone budget.

! Make sure the number of requested 3-D diagnostics is not 0
      IF (n_cpl > 0) THEN

! calculate 3-D diagnostic fields for ozone production and loss
        cpl_o3(:,:,k,:) = 0.

        DO m=1,jpbk
          p1 = rpartner_bimol(m,1)
          p2 = rpartner_bimol(m,2)

          DO i=1,rows
            DO j=1,rows
              l=(i-1) * row_length + j
! Make sure we're in the troposphere defined by pressure > 600 hPa or
! pressure > 50 hPa and O3 < 150 ppbv.

! Check if pressure > 50 hPa and
              IF (p(l) > 5000.) THEN

! Check if ozone < 150 ppbv or pres > 600 hPa
                IF ((p(l) > 60000.) .OR. (y(l,n_o3)/tnd(l) < 1.5e-7))   &
                THEN

! In the troposphere all production of O3 involves NO as reactant
                  IF (elem('NO        ',spb(1,m),spb(2,m))) THEN

! 1. calculate net O3 production in troposphere by HO2 + NO
                    IF (elem('HO2       ',spb(1,m),spb(2,m)))           &
                      cpl_o3(j,i,k,1) = rk(l,nbrkx(m))*y(l,p1)*y(l,p2)  &
                      * secs_per_step * volume(j,i,k) * convfac

! 2. calculate net O3 production in troposphere by MeOO + NO.
! Do not use channel that produces MeONO2.
                    IF (n_cpl > 1) THEN
                      IF ((spb(3,m) /= 'MeONO2    ') .AND.              &
                        elem('MeOO      ',spb(1,m),spb(2,m)))           &
                        cpl_o3(j,i,k,2) =                               &
                          rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                &
                          * secs_per_step * volume(j,i,k)               &
                          * convfac
                    END IF

! 3. calculate net O3 production in troposphere by RO2 + NO, where
!    RO2 = MeCOCH2OO, n-PrOO, i-PrOO, or EtOO
                    IF (n_cpl > 2) THEN
                      IF (elem('MeCOCH2OO ',spb(1,m),spb(2,m)) .OR.     &
                          elem('n-PrOO    ',spb(1,m),spb(2,m)) .OR.     &
                          elem('i-PrOO    ',spb(1,m),spb(2,m)) .OR.     &
                          elem('EtOO      ',spb(1,m),spb(2,m)))         &
                        cpl_o3(j,i,k,3) = cpl_o3(j,i,k,3) +             &
                          rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                &
                          * secs_per_step * volume(j,i,k)               &
                          * convfac
                    END IF
                  END IF  ! of code for tropospheric chem. production

! 4. calculate net O3 loss in troposphere by H2O + O(1D)
                  IF (n_cpl > 3) THEN
                    IF (elem('H2O       ',spb(1,m),spb(2,m)) .AND.      &
                        elem('O(1D)     ',spb(1,m),spb(2,m)))           &
                      cpl_o3(j,i,k,4) = 2.* rk(l,nbrkx(m))*y(l,p1)      &
                        * y(l,p2)                                       &
                        * secs_per_step * volume(j,i,k)                 &
                        * convfac
                  END IF

! 5. calculate net O3 loss in troposphere by HO2 + O3
                  IF (n_cpl > 4) THEN
                    IF (elem('HO2       ',spb(1,m),spb(2,m)) .AND.      &
                        elem('O3        ',spb(1,m),spb(2,m)))           &
                      cpl_o3(j,i,k,5) =                                 &
                        rk(l,nbrkx(m)) * y(l,p1) * y(l,p2)              &
                        * secs_per_step * volume(j,i,k)                 &
                        * convfac
                  END IF

! 6. calculate net O3 loss in troposphere by OH + O3

                  IF (n_cpl > 5) THEN
                    IF (elem('OH        ',spb(1,m),spb(2,m)) .AND.      &
                        elem('O3        ',spb(1,m),spb(2,m)))           &
                      cpl_o3(j,i,k,6) =                                 &
                        rk(l,nbrkx(m)) * y(l,p1) * y(l,p2)              &
                        * secs_per_step * volume(j,i,k)                 &
                        * convfac
                  END IF

                END IF
              END IF   ! end of troposphere-only diagnostics

! Diagnose stratospheric loss cycles following Lee et al., JGR, 2002.
! 7. Calculate net O3 loss from Cl2O2 cycle. Uses photolysis not
!    bimolecular rates. Done outside this loop
!
! 8. Calculate net O3 loss from BrO-ClO cycle.
              IF (n_cpl > 7) THEN
                IF (elem('BrO       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('ClO       ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,8) = cpl_o3(j,i,k,8) +                   &
                     2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                &
                     * secs_per_step * volume(j,i,k)                    &
                     * convfac
              END IF

! 9. Calculate net O3 loss from HO2 cycle.
              IF (n_cpl > 8) THEN
                IF (elem('HO2       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O3        ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,9) =                                     &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 10. Calculate net O3 loss from ClO-HO2 cycle.
              IF (n_cpl > 9) THEN
                IF (elem('HO2       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('ClO       ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,10) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 11. Calculate net O3 loss from BrO-HO2 cycle.
              IF (n_cpl > 10) THEN
                IF (elem('HO2       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('BrO       ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,11) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 12. Calculate net O3 loss from ClO-O(3P) cycle.
              IF (n_cpl > 11) THEN
                IF (elem('O(3P)     ',spb(1,m),spb(2,m)) .AND.          &
                    elem('ClO       ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,12) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 13. Calculate net O3 loss from NO2-O(3P) cycle.
              IF (n_cpl > 12) THEN
                IF (elem('NO2       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O(3P)     ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,13) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 14. Calculate net O3 loss from BrO-O(3P) cycle.
              IF (n_cpl > 13) THEN
                IF (elem('BrO       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O(3P)     ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,14) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 15. Calculate net O3 loss from HO2-O(3P) cycle.
              IF (n_cpl > 14) THEN
                IF (elem('HO2       ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O(3P)     ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,15) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 16. Calculate net O3 loss from H-O3 cycle.
              IF (n_cpl > 15) THEN
                IF (elem('H         ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O3        ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,16) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 17. Calculate net O3 loss from NO3 photolysis. Done outside this loop.

! 18. Calculate net O3 loss from Chapman loss O3 + O -> 2 O2
              IF (n_cpl > 17) THEN
                IF (elem('O(3P)     ',spb(1,m),spb(2,m)) .AND.          &
                    elem('O3        ',spb(1,m),spb(2,m)))               &
                  cpl_o3(j,i,k,18) =                                    &
                    2. * rk(l,nbrkx(m))*y(l,p1)*y(l,p2)                 &
                    * secs_per_step * volume(j,i,k)                     &
                    * convfac
              END IF

! 19. Calculate net O3 production from Chapman. Done outside the loop.

            END DO   ! close j loop
          END DO     ! close i loop
        END DO       ! close loop over bimolecular reactions

        DO m=1,jppj
          p1 = rpartner_phot(m)

          DO i=1,rows
            DO j=1,row_length
              l = (i-1) * row_length + j

! Consider Cl2O2 loop: Cl2O2 + h nu -> Cl + ClOO

              IF ((n_cpl > 6) .AND. (spj(1,m) == 'Cl2O2     '))         &
                cpl_o3(j,i,k,7) = 2. * zprt(j,i,m) * y(l,p1)            &
                * secs_per_step * volume(j,i,k) * convfac

! Consider NO3 loop: NO3 + h nu -> NO + O2. Make sure correct branch is used.
              IF ((n_cpl > 17)                  .AND.                   &
                  (spj(1,m) == 'NO3       ')    .AND.                   &
                  elem('NO        ',spj(3,m),spj(4,m)))                 &
                cpl_o3(j,i,k,18) = 2. * zprt(j,i,m) * y(l,p1)           &
                * secs_per_step * volume(j,i,k) * convfac

! Consider ozone production: O2 + h nu -> {O(3P), O(1D)} + O(3P)
              IF ((n_cpl > 19) .AND. (spj(1,m) == 'O2        '))        &
                cpl_o3(j,i,k,20) = cpl_o3(j,i,k,20) +                   &
                2. * zprt(j,i,m) * y(l,p1)                              &
              * secs_per_step * volume(j,i,k) * convfac

            END DO
          END DO
        END DO
      END IF

      CONTAINS
!============================================================================
! This function checks whether string a is in {b,c}.
      LOGICAL FUNCTION elem(a, b, c)

      IMPLICIT NONE

      CHARACTER*10 :: a,b,c

      elem = ((a == b) .OR. (a == c))
      END FUNCTION elem

      END SUBROUTINE ukca_budget2
#endif
