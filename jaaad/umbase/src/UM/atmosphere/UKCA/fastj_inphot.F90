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
!  Routine to initialise photolysis rate data, called directly from the
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction index.
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
      SUBROUTINE FASTJ_INPHOT

      USE     FASTJ_DATA,         ONLY: kpcx,dtaumax,dtausub,dsubdiv,  &
     &                                     szamax,RAD,ZZHT
      USE     FASTJ_MIE,          ONLY: FASTJ_MIE_ALLOC_MEM
      USE     FASTJ_SPECS,        ONLY: FASTJ_RD_JS, FASTJ_RD_TJPL

      IMPLICIT NONE
      INTEGER                  :: iph  ! Channel for reading data files

      CALL FASTJ_MIE_ALLOC_MEM (kpcx)

!     Use channel 8 to read files at the moment

      iph = 8

!     Defaults & constants

      rad  = 6375.d5
      zzht = 5.d5

!     Calculate new levels if tau exceeds DTAUMAX

      dtaumax = 1.d0

!     Add dsubdiv additional levels in first dtausub of cloud

      dtausub = 1.0d0
      dsubdiv = 10.d0

!     Maximum Zenith Angle  (98 degrees at 63 km; 99 degrees at 80 km)

      szamax = 98.0d0

!     Read in labels of photolysis rates required

      CALL FASTJ_RD_JS(iph,'ratj.d')

!     Read in JPL spectral data set

      CALL FASTJ_RD_TJPL(iph,'jv_spec.dat')

!     Select Aerosol/Cloud types to be used

!DEPENDS ON: fastj_set_aer
      CALL FASTJ_SET_AER

      RETURN
      END SUBROUTINE FASTJ_INPHOT
#endif
