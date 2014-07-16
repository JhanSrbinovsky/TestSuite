#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module to contain tracer and species numbers
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Current code owner: Colin Johnson/Olaf Morgenstern
!                     Fiona O'Connor
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      MODULE UKCA_CSPECIES
      IMPLICIT NONE
      PRIVATE

      INTEGER :: m           ! Loop counter

      INTEGER, SAVE, PUBLIC :: n_ox    ! tracer numbers
      INTEGER, SAVE, PUBLIC :: n_o3
      INTEGER, SAVE, PUBLIC :: n_o3s
      INTEGER, SAVE, PUBLIC :: n_nox
      INTEGER, SAVE, PUBLIC :: n_no
      INTEGER, SAVE, PUBLIC :: n_no2
      INTEGER, SAVE, PUBLIC :: n_no3
      INTEGER, SAVE, PUBLIC :: n_n2o5
      INTEGER, SAVE, PUBLIC :: n_ho2no2
      INTEGER, SAVE, PUBLIC :: n_hono2
      INTEGER, SAVE, PUBLIC :: n_ch4
      INTEGER, SAVE, PUBLIC :: n_sx
      INTEGER, SAVE, PUBLIC :: n_h2
      INTEGER, SAVE, PUBLIC :: n_h2o
      INTEGER, SAVE, PUBLIC :: n_ox1
      INTEGER, SAVE, PUBLIC :: n_sx1
      INTEGER, SAVE, PUBLIC :: n_ro3 ! ozone tracer for radiative purposes

      INTEGER, SAVE, PUBLIC :: nn_o3    ! Species numbers
      INTEGER, SAVE, PUBLIC :: nn_o3s
      INTEGER, SAVE, PUBLIC :: nn_oh
      INTEGER, SAVE, PUBLIC :: nn_ho2
      INTEGER, SAVE, PUBLIC :: nn_no
      INTEGER, SAVE, PUBLIC :: nn_no2
      INTEGER, SAVE, PUBLIC :: nn_o1d
      INTEGER, SAVE, PUBLIC :: nn_o3p
      INTEGER, SAVE, PUBLIC :: nn_meoo
      INTEGER, SAVE, PUBLIC :: nn_meco3
      INTEGER, SAVE, PUBLIC :: nn_etoo
      INTEGER, SAVE, PUBLIC :: nn_etco3
      INTEGER, SAVE, PUBLIC :: nn_nproo
      INTEGER, SAVE, PUBLIC :: nn_nprooh
      INTEGER, SAVE, PUBLIC :: nn_iProo
      INTEGER, SAVE, PUBLIC :: nn_iProoh
      INTEGER, SAVE, PUBLIC :: nn_mecoch2oo
      INTEGER, SAVE, PUBLIC :: nn_no3
      INTEGER, SAVE, PUBLIC :: nn_n2o5
      INTEGER, SAVE, PUBLIC :: nn_ho2no2
      INTEGER, SAVE, PUBLIC :: nn_hono2
      INTEGER, SAVE, PUBLIC :: nn_ch4
      INTEGER, SAVE, PUBLIC :: nn_so2
      INTEGER, SAVE, PUBLIC :: nn_n
      INTEGER, SAVE, PUBLIC :: nn_cl
      INTEGER, SAVE, PUBLIC :: nn_clo
      INTEGER, SAVE, PUBLIC :: nn_hcl
      INTEGER, SAVE, PUBLIC :: nn_cl2o2
      INTEGER, SAVE, PUBLIC :: nn_meo
      INTEGER, SAVE, PUBLIC :: nn_bro
      INTEGER, SAVE, PUBLIC :: nn_br


      REAL, SAVE, ALLOCATABLE, PUBLIC :: c_species(:)

      PUBLIC UKCA_CALC_CSPECIES

      CONTAINS

      SUBROUTINE UKCA_CALC_CSPECIES

      USE ASAD_MOD,     ONLY: advt, speci
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "c_v_m.h"

      INTEGER :: i
      CHARACTER(LEN=72) :: cmessage

!     Compute list of conversion factors from mmr to vmr.
!     Change this if new tracers are introduced !!!!!!!!!!!!!!!!!!

      c_species = 0.
      WHERE (advt.eq.'Ox        ') c_species = c_o3
      WHERE (advt.eq.'O3        ') c_species = c_o3
      WHERE (advt.eq.'NO        ') c_species = c_no
      WHERE (advt.eq.'NO2       ') c_species = c_no2
      WHERE (advt.eq.'NO3       ') c_species = c_no3
      WHERE (advt.eq.'NOx       ') c_species = c_no2
      WHERE (advt.eq.'N2O5      ') c_species = c_n2o5
      WHERE (advt.eq.'HO2NO2    ') c_species = c_ho2no2
      WHERE (advt.eq.'HONO2     ') c_species = c_hono2
      WHERE (advt.eq.'H2O2      ') c_species = c_h2o2
      WHERE (advt.eq.'CH4       ') c_species = c_ch4
      WHERE (advt.eq.'CO        ') c_species = c_co
      WHERE (advt.eq.'HCHO      ') c_species = c_hcho
      WHERE (advt.eq.'MeOOH     ') c_species = c_meooh
      WHERE (advt.eq.'HONO      ') c_species = c_hono
      WHERE (advt.eq.'C2H6      ') c_species = c_c2h6
      WHERE (advt.eq.'EtOOH     ') c_species = c_etooh
      WHERE (advt.eq.'MeCHO     ') c_species = c_mecho
      WHERE (advt.eq.'PAN       ') c_species = c_pan
      WHERE (advt.eq.'C3H8      ') c_species = c_c3h8
      WHERE (advt.eq.'i-PrOOH   ') c_species = c_prooh
      WHERE (advt.eq.'n-PrOOH   ') c_species = c_prooh
      WHERE (advt.eq.'EtCHO     ') c_species = c_etcho
      WHERE (advt.eq.'Me2CO     ') c_species = c_me2co
      WHERE (advt.eq.'MeCOCH2OOH') c_species = c_mecoch2ooh
      WHERE (advt.eq.'PPAN      ') c_species = c_ppan
      WHERE (advt.eq.'MeONO2    ') c_species = c_meono2
      WHERE (advt.eq.'Sx        ') c_species = c_o3
      WHERE (advt.eq.'O3S       ') c_species = c_o3
      WHERE (advt.eq.'HOx       ') c_species = c_ho2
      WHERE (advt.eq.'N2O       ') c_species = c_n2o
      WHERE (advt.eq.'CFCl3     ') c_species = c_cfcl3
      WHERE (advt.eq.'CF2Cl2    ') c_species = c_cf2cl2
      WHERE (advt.eq.'H2O       ') c_species = c_h2o
      WHERE (advt.eq.'ClONO2    ') c_species = c_clono2
      WHERE (advt.eq.'Clx       ') c_species = c_clo
      WHERE (advt.eq.'HCl       ') c_species = c_hcl
      WHERE (advt.eq.'HOCl      ') c_species = c_hocl
      WHERE (advt.eq.'OClO      ') c_species = c_oclo
      WHERE (advt.eq.'Brx       ') c_species = c_bro
      WHERE (advt.eq.'HOBr      ') c_species = c_hobr
      WHERE (advt.eq.'BrONO2    ') c_species = c_brono2
      WHERE (advt.eq.'BrCl      ') c_species = c_brcl
      WHERE (advt.eq.'MeBr      ') c_species = c_mebr
      WHERE (advt.eq.'HBr       ') c_species = c_hbr
      WHERE (advt.eq.'CF2ClCFCl2') c_species = c_cf2clcfcl2
      WHERE (advt.eq.'CHF2Cl    ') c_species = c_chf2cl
      WHERE (advt.eq.'MeCCl3    ') c_species = c_meccl3
      WHERE (advt.eq.'CCl4      ') c_species = c_ccl4
      WHERE (advt.eq.'MeCl      ') c_species = c_mecl
      WHERE (advt.eq.'CF2ClBr   ') c_species = c_cf2clbr
      WHERE (advt.eq.'CF3Br     ') c_species = c_cf3br
      WHERE (advt.eq.'C5H8      ') c_species = c_c5h8
      WHERE (advt.eq.'ISO2      ') c_species = c_iso2
      WHERE (advt.eq.'ISOOH     ') c_species = c_isooh
      WHERE (advt.eq.'ISON      ') c_species = c_ison
      WHERE (advt.eq.'MACR      ') c_species = c_macr
      WHERE (advt.eq.'MACRO2    ') c_species = c_macro2
      WHERE (advt.eq.'MACROOH   ') c_species = c_macrooh
      WHERE (advt.eq.'MPAN      ') c_species = c_mpan
      WHERE (advt.eq.'HACET     ') c_species = c_hacet
      WHERE (advt.eq.'MGLY      ') c_species = c_mgly
      WHERE (advt.eq.'NALD      ') c_species = c_nald
      WHERE (advt.eq.'HCOOH     ') c_species = c_hcooh
      WHERE (advt.eq.'MeCO3H    ') c_species = c_meco3h
      WHERE (advt.eq.'MeCO2H    ') c_species = c_meco2h
      WHERE (advt.eq.'SO2       ') c_species = c_so2
      WHERE (advt.eq.'DMS       ') c_species = c_dms
      WHERE (advt.eq.'Me2S      ') c_species = c_me2s
      WHERE (advt.eq.'MSA       ') c_species = c_msa
      WHERE (advt.eq.'H2SO4     ') c_species = c_h2so4
      WHERE (advt.eq.'NH3       ') c_species = c_nh3

!     Initialise tracer numbers

      n_h2o    = 0
      n_sx     = 0
      n_n2o5   = 0
      n_ox     = 0
      n_o3     = 0
      n_o3s    = 0
      n_nox    = 0
      n_no     = 0
      n_no2    = 0
      n_no3    = 0
      n_ho2no2 = 0
      n_hono2  = 0
      n_ch4    = 0
      n_h2     = 0

!     Find tracer numbers

      DO m=1,jpctr
        SELECT CASE (advt(m))
          CASE ('Ox        ')
            n_ox     = m
            n_ox1    = m
          CASE ('O3        ')
            n_o3     = m
          CASE ('O3S       ')
            n_o3s    = m
          CASE ('Sx        ')   ! stratospheric ozone tracer
            n_sx     = m
            n_sx1    = m
          CASE ('NOx       ')
            n_nox    = m
          CASE ('NO        ')
            n_no     = m
          CASE ('NO2       ')
            n_no2    = m
          CASE ('NO3       ')
            n_no3    = m
          CASE ('N2O5      ')
            n_n2o5   = m
          CASE ('HO2NO2    ')
            n_ho2no2 = m
          CASE ('HONO2     ')
            n_hono2  = m
          CASE ('CH4       ')
            n_ch4    = m
          CASE ('H2O       ')
            n_h2o    = m
          CASE ('H2        ')
            n_h2     = m
          CASE ('RO3       ')
            n_ro3    = m
        END SELECT
      ENDDO

      DO i=1,jpctr
        IF ((advt(i) == 'Ox        ') .OR. (advt(i) == 'O3        ')    &
          .OR. (advt(i) == 'RO3       '))                               &
          n_o3 = i
      END DO
      IF (n_o3 == 0) THEN
        cmessage='Ozone not found among chemical tracers.'
! DEPENDS ON: ereport
        CALL ereport('UKCA_CSPECIES',1,cmessage)
      ENDIF

!     Initialise species numbers

      nn_o3        = 0
      nn_o3s       = 0
      nn_oh        = 0
      nn_ho2       = 0
      nn_no        = 0
      nn_no2       = 0
      nn_o1d       = 0
      nn_o3p       = 0
      nn_meoo      = 0
      nn_meco3     = 0
      nn_etoo      = 0
      nn_etco3     = 0
      nn_nproo     = 0
      nn_nprooh    = 0
      nn_iproo     = 0
      nn_iprooh    = 0
      nn_mecoch2oo = 0
      nn_so2       = 0
      nn_n         = 0
      nn_meo       = 0
      nn_cl        = 0
      nn_clo       = 0
      nn_cl2o2     = 0
      nn_hcl       = 0
      nn_br        = 0
      nn_bro       = 0

!     Find species numbers

      DO m=1,jpspec
        SELECT CASE(speci(m))
          CASE ('O3        ')
            nn_o3 = m
          CASE ('O3S       ')
            nn_o3s = m
          CASE ('OH        ')
            nn_oh = m
          CASE ('HO2       ')
            nn_ho2 = m
          CASE ('NO        ')
            nn_no = m
          CASE ('NO2       ')
            nn_no2 = m
          CASE ('O(1D)     ')
            nn_o1d = m
          CASE ('O(3P)     ')
            nn_o3p = m
          CASE ('MeOO      ')
            nn_meoo = m
          CASE ('MeCO3     ')
            nn_meco3 = m
          CASE ('EtOO      ')
            nn_etoo = m
          CASE ('EtCO3     ')
            nn_etco3 = m
          CASE ('n-PrOO    ')
            nn_nProo = m
          CASE ('n-PrOOH   ')
            nn_nprooh = m
          CASE ('i-PrOO    ')
            nn_iProo = m
          CASE ('i-PrOOH   ')
            nn_iprooh = m
          CASE ('MeCOCH2OO ')
            nn_mecoch2oo = m
          CASE ('NO3       ')
            nn_no3 = m
          CASE ('N2O5      ')
            nn_n2o5 = m
          CASE ('HO2NO2    ')
            nn_ho2no2 = m
          CASE ('HONO2     ')
            nn_hono2 = m
          CASE ('CH4       ')
            nn_ch4 = m
          CASE ('SO2       ')
            nn_so2 = m
          CASE ('N         ')
            nn_n = m
          CASE ('MeO       ')
            nn_meo = m
          CASE ('Cl        ')
            nn_cl = m
          CASE ('ClO       ')
            nn_clo = m
          CASE ('Cl2O2     ')
            nn_cl2o2 = m
          CASE ('HCl       ')
            nn_hcl = m
          CASE ('Br        ')
            nn_br = m
          CASE ('BrO       ')
            nn_bro = m
        END SELECT
      ENDDO
      WRITE(6,*) 'Species indices:'
      WRITE(6,*) 'O3 = ',nn_o3,' OH = ',nn_oh
      WRITE(6,*) 'HO2 = ',nn_ho2,' NO = ',nn_no
      WRITE(6,*) 'NO2 = ',nn_no2,' O(1D) = ',nn_o1d
      WRITE(6,*) 'MeOO = ',nn_meoo,' MeCO3 = ',nn_meco3
      WRITE(6,*) 'EtOO = ',nn_etoo,' EtCO3 = ',nn_etco3
      WRITE(6,*) 'n-Proo = ',nn_nproo,'n-Prooh = ',nn_nprooh
      WRITE(6,*) 'i-Proo = ',nn_iproo,'i-Prooh = ',nn_iprooh
      WRITE(6,*) 'MeCOCH2OO = ',nn_mecoch2oo

      END SUBROUTINE UKCA_CALC_CSPECIES

      END MODULE UKCA_CSPECIES
#endif
