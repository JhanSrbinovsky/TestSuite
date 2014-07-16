#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-j routine for calculating online photolysis rates
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Current Code Owner:       Colin Johsnon/Fiona O'Connor
!                            Oliver Wild
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
       MODULE FASTJ_DATA

       IMPLICIT NONE
       PUBLIC
       SAVE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

!      Setup Blocking of colums for complet fastj model 
!
!      Blocking mode          1: rows
!                             2: whole domain
!                             3: compressed   (Not implemented yet)

       TYPE Block_setup
         INTEGER                               :: Mode=4
         INTEGER                               :: Bl_size=101
       END TYPE Block_setup

       TYPE (Block_setup)                      :: Blocking

!      Parameters
!      ----------
!
!      NB    Number of levels in CTM plus one for above model top
!      NC    Number of levels in the fundamental Fast-J grid
!      NS    Maximum number of species which require J-values calculating
!      NW    Maximum number of wavelength bins that can be used
!      NP    Maximum number of aerosol/cloud types that can be used
!      NH    Maximum number of Herzberg X-sections that can be used
!      MX    Number of aerosol/cloud types supplied from CTM

       INTEGER,PARAMETER             :: NS = 51
       INTEGER,PARAMETER             :: NW = 15
       INTEGER,PARAMETER             :: NP = 25
       INTEGER,PARAMETER             :: NH =  7
       INTEGER,PARAMETER             :: MX =  3 !kk look at set_aer!

       INTEGER                       :: ipcx    !Number of longitudes
       INTEGER                       :: jpcx    !Number of latitudes
       INTEGER                       :: lpar    !Number of levels
       INTEGER                       :: jpcl    !Number of chem. levels
       INTEGER                       :: month   !Number of month (1-12)
       INTEGER                       :: iday    !Day of year
       INTEGER                       :: nslat    !  Latitude of current profile point
       INTEGER                       :: nslon    !  Longitude of current profile point
       INTEGER                       :: NB, NC
       INTEGER                       :: nhz,npdep
       INTEGER                       :: NJVAL,NW1,NW2,NAA,NLBATM
       INTEGER                       :: kpcx, jjpnl, lwepar, lwdpar
       
       INTEGER, ALLOCATABLE               :: jind(:)
       INTEGER,DIMENSION(nh)              :: hzind
       INTEGER,DIMENSION(:,:),pointer     :: nsl
       INTEGER,DIMENSION(mx)              :: MIEDX
       INTEGER,DIMENSION(ns)              :: jpdep
              
!      Logical for degraded CTM resn

       LOGICAL,PARAMETER                  :: ldeg45 = .false.

       CHARACTER(LEN=20),DIMENSION(NP)    :: TITLEA
       CHARACTER(LEN=7), ALLOCATABLE      :: jlabel(:)
       
       REAL                               :: tau     !Time of Day (hours,GMT)
       REAL                               :: RAD,TANHT,ZZHT
       REAL                               :: dtaumax,szamax
       REAL                               :: dtausub,dsubdiv
       REAL                               :: hzo2,hzo3
                     
       REAL,DIMENSION(:),POINTER          :: etaa ! eta(a) for level boundaries
       REAL,DIMENSION(:),POINTER          :: etab ! eta(b) for level boundaries

       REAL,DIMENSION(:,:,:),POINTER      :: p     !Surface pressure
       REAL,DIMENSION(:,:),POINTER        :: sa    !Surface albedo
       REAL,DIMENSION(:,:,:),POINTER      :: T     !Temperature
       REAL,DIMENSION(:,:,:),POINTER      :: OD    !Optical Depth
       REAL,DIMENSION(:,:,:),POINTER      :: ODW   !OD for water
       REAL,DIMENSION(:,:,:),POINTER      :: ODI   !OD for ice
       REAL,DIMENSION(:),POINTER          :: RFLECT,SZA,U0,SZAFAC
       REAL,DIMENSION(:,:),POINTER        :: TJ,PJ,DM,DO3,Z
       REAL,DIMENSION(:,:,:),POINTER      :: AER,AMF
       REAL,DIMENSION(NW+1)               :: WBIN,QRAYL
       REAL,DIMENSION(NW)                 :: WL,FL
       REAL,DIMENSION(NW,3)               :: QO2,QO3,Q1D
       REAL,DIMENSION(3,NS)               :: TQQ
       REAL,DIMENSION(:,:),POINTER        :: DBC
       REAL,DIMENSION(:,:),POINTER        :: FFF,VALJ
       REAL,DIMENSION(NW,2,NS-3)          :: QQQ
       REAL,DIMENSION(4,NP)               :: WAA,QAA,RAA,SSA
       REAL,DIMENSION(8,4,NP)             :: PAA
       REAL,DIMENSION(NW)                 :: QBC
       REAL,DIMENSION(NW,3)               :: zpdep
       REAL,DIMENSION(MX)                 :: ZMIEX
       REAL,DIMENSION(51,18,12)           :: TREF,OREF
       REAL,DIMENSION(51)                 :: BREF,AREF
       REAL,DIMENSION(:,:),POINTER        :: zj
       REAL,DIMENSION(NH)                 :: hztoa
       REAL,DIMENSION(:),POINTER          :: fhz,jfacta
       REAL,DIMENSION(:,:),POINTER        :: SZA_2d,SZAFAC_2d
       REAL,DIMENSION(:,:,:),POINTER      :: sulphur
       REAL,DIMENSION(:,:,:),POINTER      :: rz_3d
       REAL,DIMENSION(:,:,:),POINTER      :: fj_ozone

!      Interface section

       INTERFACE FASTJ_SET_LIMITS
         MODULE PROCEDURE FASTJ_SET_LIMITS
       END INTERFACE FASTJ_SET_LIMITS

       INTERFACE FASTJ_ALLOCATE_MEMORY
         MODULE PROCEDURE  FASTJ_ALLOCATE_MEMORY
       END INTERFACE FASTJ_ALLOCATE_MEMORY

      CONTAINS

! ######################################################################
        SUBROUTINE FASTJ_SET_LIMITS (lon, lat, lev, chemlev)
         IMPLICIT  NONE

         INTEGER,INTENT(IN)            :: lon, lat, lev
         INTEGER,INTENT(IN),OPTIONAL   :: chemlev

         ipcx   = lon
         jpcx   = lat
         lpar   = lev
         IF (PRESENT(chemlev))  THEN
           jpcl = chemlev
         ELSE
           jpcl = lev
         END IF

         lwepar = lpar
         lwdpar = lpar

         NB = LPAR+1
         NC = 2*NB
        
         SELECT CASE (Blocking%Mode)
          CASE(1)
           kpcx  = ipcx                 !Blocking row-wise
           jjpnl  = jpcl*kpcx
          CASE(2)
           kpcx  = ipcx*jpcx            !Blocking domain
           jjpnl  = jpcl*kpcx
          CASE(3)
           kpcx  = Blocking%Bl_size
           jjpnl  = jpcl*kpcx
          CASE(4)
           kpcx  = Blocking%Bl_size
           jjpnl  = jpcl*kpcx
         END SELECT
         
         RETURN
         END SUBROUTINE FASTJ_SET_LIMITS

! ######################################################################
        SUBROUTINE FASTJ_ALLOCATE_MEMORY

         IMPLICIT  NONE

         WRITE(6,*) 'Allocate memory ',ipcx,jpcx,kpcx,lpar

         ALLOCATE(etaa(lpar+1))
         ALLOCATE(etab(lpar+1))

         ALLOCATE(P(ipcx,jpcx,NB))
         ALLOCATE(SA(ipcx,jpcx))
         ALLOCATE(T(ipcx,jpcx,lpar))
         ALLOCATE(OD(ipcx,jpcx,lpar))
         ALLOCATE(ODW(ipcx,jpcx,lpar))
         ALLOCATE(ODI(ipcx,jpcx,lpar))

         ALLOCATE(nsl(2,kpcx))
         ALLOCATE(RFLECT(kpcx))
         ALLOCATE(sza(kpcx))
         ALLOCATE(U0(kpcx))
         ALLOCATE(SZAFAC(kpcx))
         ALLOCATE(TJ(kpcx,NB))
         ALLOCATE(PJ(kpcx,NB+1))
         ALLOCATE(DM(kpcx,NB))
         ALLOCATE(DO3(kpcx,NB))
         ALLOCATE(Z(kpcx,NB))
         ALLOCATE(AER(MX,kpcx,NB))
         ALLOCATE(AMF(kpcx,NB,NB))
         ALLOCATE(DBC(kpcx,NB))
         ALLOCATE(FFF(NW*kpcx,jpcl))
         ALLOCATE(VALJ(kpcx,NS))
         ALLOCATE(zj(jjpnl,jppj))
         ALLOCATE(jfacta(jppj))
         ALLOCATE(fhz(jjpnl))

         ALLOCATE(SZA_2d(ipcx,jpcx))
         ALLOCATE(SZAFAC_2d(ipcx,jpcx))
       
         ALLOCATE(sulphur(ipcx,jpcx,lpar))
         ALLOCATE(rz_3d(ipcx,jpcx,NB))
         ALLOCATE(fj_ozone(ipcx,jpcx,NB))

         ALLOCATE(jind(jppj))
         ALLOCATE(jlabel(jppj))


         RETURN
        END SUBROUTINE FASTJ_ALLOCATE_MEMORY
      END MODULE FASTJ_DATA
#endif
