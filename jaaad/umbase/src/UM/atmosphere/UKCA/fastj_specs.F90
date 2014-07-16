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
      MODULE FASTJ_SPECS
      USE    FASTJ_DATA
      IMPLICIT NONE

      CHARACTER (LEN=78)                          :: TITLE0
      CHARACTER (LEN=7),DIMENSION(3,NS)           :: TITLEJ
      CHARACTER (LEN=7),DIMENSION(NH)             :: hzlab

      INTERFACE FASTJ_RD_JS
        MODULE PROCEDURE FASTJ_RD_JS
      END INTERFACE FASTJ_RD_JS

      INTERFACE FASTJ_RD_TJPL
        MODULE PROCEDURE FASTJ_RD_TJPL
      END INTERFACE FASTJ_RD_TJPL

      PUBLIC FASTJ_RD_JS, FASTJ_RD_TJPL

      CONTAINS

      SUBROUTINE FASTJ_RD_JS(nj1,namfil)


!     Reread the ratj.d file to map photolysis rate to reaction
!     Read in quantum yield 'jfacta' and fastj label 'jlabel'


!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use
!     ipr       Photolysis reaction counter - should total 'jppj'


      INTEGER             :: nj1, ipr, i

      CHARACTER (LEN=6)   ::  namfil
      CHARACTER (LEN=120) :: cline

!     Reread the ratj.d file to map photolysis rate to reaction
!     Read in quantum yield jfacta and fastj label jlabel

      ipr = 0
      OPEN(nj1,FILE=namfil,STATUS='old',FORM='formatted')
 10   READ(nj1,'(a)',ERR=20) cline
      IF (cline(2:5).eq.'9999') THEN
         go to 20
      ELSE IF (cline(1:1).eq.'#') THEN
         go to 10
      ELSE IF (cline(5:5).eq.'$') THEN
         go to 10
      ELSE
         ipr = ipr+1
         READ(cline(90:94),'(f5.1)') jfacta(ipr)
         READ(cline(97:103),'(a7)')   jlabel(ipr)
         jfacta(ipr) = jfacta(ipr)/100.d0
         GO TO 10
      ENDIF
 20   CLOSE(NJ1)

      IF (ipr /= jppj) THEN
         WRITE(6,*) ipr,jppj
!DEPENDS ON: ereport
         CALL EREPORT("FASTJ",1,"rd_js_1")
      ENDIF

      RETURN
      END SUBROUTINE FASTJ_RD_JS

      SUBROUTINE FASTJ_RD_TJPL(NJ1,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh parameters, temperature-
!  dependent cross sections and Rayleigh/aerosol scattering phase functions
!  with temperature dependences. Current data originates from JPL'97
!-----------------------------------------------------------------------

!     NAMFIL   Name of spectral data file (jv_spec.dat)
!     NJ1      Channel number for reading data file
!     NJVAL    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from NW1:NW2
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QBC      Black Carbon absorption extinct. (specific cross-sect.) 
!              (m2/g)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     NK       Number of wavelengths at which functions supplied (set as 4)
!     WAA      Wavelengths for the NK supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SSA      Single scattering albedo
!
!     npdep    Number of pressure dependencies
!     zpdep    Pressure dependencies by wavelength bin
!     jpdep    Index of cross sections requiring pressure dependence
!     lpdep    Label for pressure dependence
!
!-----------------------------------------------------------------------

      INTEGER            :: i, j, k, iw, nk, nqqq, nwww, nj1

      CHARACTER (LEN=7)  :: lpdep(3)
      CHARACTER (LEN=11) :: NAMFIL

      REAL               :: hztmp(nh)

      DO j = 1,ns
        DO k = 1,3
          tqq(k,j) = 0.d0
        ENDDO
      ENDDO

!-------------spectral data---------------------------------------------

      OPEN (NJ1, FILE=NAMFIL)
      READ (NJ1,'(A)') TITLE0
      WRITE(6,'(1X,A)') TITLE0
      READ(NJ1,'(10X,14I5)') njval,nwww,nw1,nw2

      IF (njval > ns) THEN
        WRITE(6,*) njval,ns
!DEPENDS ON: ereport
        CALL EREPORT("FASTJ",1,"RD_TJPL")
      ENDIF

!     NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)

      nqqq = njval-3
      READ(NJ1,'(10X,7E10.3)') (wbin(iw),iw=1,nwww)
      READ(NJ1,'(10X,7E10.3)') (wbin(iw+1),iw=1,nwww)
      READ(NJ1,'(10X,7E10.3)') (wl(iw),iw=1,nwww)
      READ(NJ1,'(10X,7E10.3)') (fl(iw),iw=1,nwww)
      READ(NJ1,'(10X,7E10.3)') (qrayl(iw),iw=1,nwww)
      READ(NJ1,'(10X,7E10.3)') (qbc(iw),iw=1,nwww)   !  From Loiusse et al. 
                                          !  [JGR, 1996]

!     Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)

      DO K=1,3
        READ(NJ1,'(A7,F3.0,7E10.3)') titlej(k,1),tqq(k,1),              &
     &                               (qo2(iw,k),iw=1,nwww)
      ENDDO
      DO K=1,3
        READ(NJ1,'(A7,F3.0,7E10.3)') titlej(k,2),tqq(k,2),              &
     &                               (qo3(iw,k),iw=1,nwww)
      ENDDO
      DO K=1,3
        READ(NJ1,'(A7,F3.0,7E10.3)') titlej(k,3),tqq(k,3),              &
     &                               (q1d(iw,k),iw=1,nwww)
      ENDDO
      DO K=1,3
        WRITE(6,*) 'X-sec ',titlej(1,k),(tqq(i,k),i=1,3)
      ENDDO

!     Read remaining species:  X-sections at 2 T's

      DO J=1,NQQQ
        READ(NJ1,'(A7,F3.0,7E10.3)') titlej(1,j+3),tqq(1,j+3),          &
     &                               (qqq(iw,1,j),iw=1,nwww)
        READ(NJ1,'(A7,F3.0,7E10.3)') titlej(2,j+3),tqq(2,j+3),          &
     &                               (qqq(iw,2,j),iw=1,nwww)
        WRITE(6,*)  'X-sec ',titlej(1,j+3),(tqq(i,j+3),i=1,2)
      ENDDO
      READ(NJ1,'(A)') title0

!     Pressure dependencies

      READ(NJ1,'(13x,i2)') npdep
      DO k = 1,npdep
        READ(NJ1,'(A7,3x,7E10.3)') lpdep(k),(zpdep(iw,k),iw=1,nwww)
        WRITE(6,*) 'pr.dep ',lpdep(k),(zpdep(iw,k),iw=1,nwww)
      ENDDO
      READ(NJ1,'(A)') title0

!     Read aerosol phase functions:

      READ(NJ1,'(A10,I5,/)') title0,naa
      if (naa > np) THEN
        WRITE(6,*) NAA
!DEPENDS ON: ereport
        CALL EREPORT("FASTJ",1,"RD_TJPL_2")
      ENDIF

      NK = 4        ! Fix number of wavelengths at 4
      DO j = 1,NAA
        READ(NJ1,'(3x,a20)') TITLEA(j)
        DO k = 1,NK
          READ(NJ1,*) waa(k,j),qaa(k,j),raa(k,j),ssa(k,j),          &
     &                                 (paa(i,k,j),i=1,8)
        ENDDO
      ENDDO

      WRITE(6,*) 'Aerosol phase functions & wavelengths ',NAA
      DO J = 1,NAA
        WRITE(6,'(1x,A8,I2,A,9F8.1)')                               &
     &                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
        WRITE(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
      ENDDO

      NHZ = 0

!     Zero index arrays

      jind (1:jppj)  = 0
      jpdep(1:njval) = 0
      hzind(1:nh)    = 0      

!     Set mapping index

      DO j = 1,NJVAL
        DO k = 1,jppj
          IF (jlabel(k) == titlej(1,j)) jind(k) = j
        ENDDO
        DO k = 1,npdep
          IF (lpdep(k) == titlej(1,j)) jpdep(j) = k
        ENDDO
      ENDDO

      DO k = 1,jppj
        IF (jfacta(k) == 0.d0)                                    &
     &     WRITE(6,*) 'Not using photolysis reaction ',k
        IF (jind(k) == 0) THEN
          IF (jfacta(k) == 0.d0) THEN
            jind(k)=1
          ELSE
            WRITE(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
!DEPENDS ON: ereport
            CALL EREPORT("FASTJ",1,"RD_TJPL_3")
          ENDIF
        ENDIF
      ENDDO

!     Herzberg index

      i = 0
      DO j=1,nhz
        DO k = 1,jppj
          IF (jlabel(k) == hzlab(j)) THEN
            i        = i+1
            hzind(i) = k
            hztoa(i) = hztmp(j)*jfacta(k)
          ENDIF
        ENDDO
      ENDDO
      nhz = i
      IF (nhz == 0) THEN
        WRITE(6,*) 'Not using Herzberg bin'
      ELSE
        WRITE(6,*) 'Using Herzberg bin for ',                          &
     &                (jlabel(hzind(i)),i=1,nhz)
      ENDIF

      CLOSE(NJ1)

      RETURN
      END SUBROUTINE FASTJ_RD_TJPL
      END MODULE FASTJ_SPECS
#endif
