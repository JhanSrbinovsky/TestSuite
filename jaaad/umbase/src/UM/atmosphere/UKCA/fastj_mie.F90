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
      MODULE  FASTJ_MIE
      IMPLICIT NONE
      PUBLIC
      SAVE

      INTEGER, PARAMETER                      :: NL=2700
!                        Max no of levs after insertion of extra Mie levs

      INTEGER, PARAMETER                      :: nlevs_mie=2*NL
!                        No of levs in Mie grid: 2*(2*lpar+2+jaddto(1))+3

      INTEGER, PARAMETER                      :: m_gauss=4
!                        Number of Gauss points used (fixed at 4)

      INTEGER, PARAMETER                      :: M=1
!                        Solve eqn of R.T. only for first-order M=1

      INTEGER, PARAMETER                      :: N=m_gauss
!                        Number of quadrature points (fixed at 4)

      INTEGER, PARAMETER                      :: MFIT=2*m_gauss
!                        Expansion of phase function (fixed at 8)

      INTEGER, PARAMETER                      :: NWB1=1
!                        First wlength bin to use - hardwired for vector'n      

      INTEGER, PARAMETER                      :: NWB2=7
!                        Last wlength bin to use - also use for dimensioning

      INTEGER                                 :: NWB3
      INTEGER                                 :: ND

      INTEGER,DIMENSION(:),POINTER            :: jadsub
      
      REAL                                    :: RADIUS

      REAL,DIMENSION(m_gauss)                 :: wt,emu
      REAL,DIMENSION(m_gauss,2*m_gauss)       :: PM
      REAL,DIMENSION(:),POINTER              :: ZREFL,ZFLUX,ZU0
      REAL,DIMENSION(:,:),POINTER            :: v1,PM0
      REAL,DIMENSION(:,:),POINTER            :: ZTAU,FZ,FJ
      REAL,DIMENSION(:,:,:),POINTER          :: a,c1,h,aa
      REAL,DIMENSION(:,:,:),POINTER          :: cc,s
      REAL,DIMENSION(:,:,:),POINTER          :: POMEGA
      REAL,DIMENSION(:,:,:),POINTER          :: rr
      REAL,DIMENSION(:,:,:,:),POINTER        :: b,dd,w,u1

      INTERFACE FASTJ_MIE_ALLOC_MEM
        MODULE PROCEDURE FASTJ_MIE_ALLOC_MEM
      END INTERFACE FASTJ_MIE_ALLOC_MEM

      CONTAINS

       SUBROUTINE FASTJ_MIE_ALLOC_MEM (kpcx_lo)
         USE        FASTJ_DATA,         ONLY: NC
         IMPLICIT   NONE

         INTEGER,INTENT(IN)                    :: kpcx_lo

         NWB3 = NWB2*kpcx_lo

         ALLOCATE(A(NWB3,m_gauss,nlevs_mie))
         ALLOCATE(B(m_gauss+1,m_gauss+1,NWB3,nlevs_mie))
         ALLOCATE(C1(NWB3,m_gauss,nlevs_mie))
         ALLOCATE(H(NWB3,m_gauss,nlevs_mie))
         ALLOCATE(AA(NWB3,m_gauss,m_gauss))
         ALLOCATE(CC(NWB3,m_gauss,m_gauss))
         ALLOCATE(S(NWB3,m_gauss,m_gauss))
         ALLOCATE(W(NWB3,m_gauss,m_gauss,2))
         ALLOCATE(U1(NWB3,m_gauss,m_gauss,2))
         ALLOCATE(V1(NWB3,m_gauss))
         ALLOCATE(PM0(kpcx_lo,2*m_gauss))
         ALLOCATE(POMEGA(NWB3,2*m_gauss,nlevs_mie))
         ALLOCATE(ZTAU(NWB3,nlevs_mie))
         ALLOCATE(FZ(NWB3,nlevs_mie))
         ALLOCATE(FJ(NWB3,nlevs_mie))
         ALLOCATE(DD(NWB3,m_gauss,m_gauss,nlevs_mie))
         ALLOCATE(RR(NWB3,m_gauss,nlevs_mie))
         ALLOCATE(ZREFL(kpcx_lo))
         ALLOCATE(ZFLUX(NWB3))
         ALLOCATE(ZU0(kpcx_lo))
         ALLOCATE(jadsub(NC))

         RETURN
       END SUBROUTINE FASTJ_MIE_ALLOC_MEM
       
      END MODULE   FASTJ_MIE
#endif
