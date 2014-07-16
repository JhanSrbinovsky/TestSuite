#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Checks component and total average masses MD, MDT [per ptcl] and
!    number concentrations ND for each mode for bad values.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_CHECK_MD_ND(NBOX,PROCESS,JLAT,ND,MD,MDT)
!-----------------------------------------------------------
!
!   Purpose
!   -------
!   Checks cpt and total average masses MD, MDT [per ptcl] and
!   number concentrations ND for each mode for bad values.
!
!   Bad is taken to be:
!
!   MD : greater than 10^20 or less than zero.
!   MDT: greater than 10^20 or less than or equal to zero.
!   ND : greater than 10^9 per cc or less than zero.
!
!   Inputs
!   ------
!   NBOX      : Number of grid boxes
!   PROCESS   : Character string with process just completed
!   JLAT      : Index of latitude array
!   ND        : Aerosol ptcl number density for mode (cm^-3)
!   MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
!   MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
!
!   Outputs
!   -------
!   None
!
!   Inputted by module UKCA_MODE_SETUP
!   ----------------------------------
!   NMODES    : Number of possible aerosol modes
!   NCP       : Number of possible aerosol components
!   MODE      : Logical variable defining which modes are set.
!   COMPONENT : Logical variable defining which cpt are in which dsts
!
!-----------------------------------------------------------
      USE UKCA_MODE_SETUP
      IMPLICIT NONE

! Subroutine interface
      INTEGER :: NBOX
      INTEGER :: JLAT
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      CHARACTER*30 :: PROCESS

! Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: TOTERR
      LOGICAL :: LOGIC1(NBOX)
      LOGICAL :: LOGIC2(NBOX)
      LOGICAL :: LOGIC12(NBOX)
      LOGICAL :: LOGIC3
      CHARACTER(LEN=72) :: cmessage

      TOTERR=0
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        LOGIC1(:)=((ND (:,IMODE) <  1.0E9).AND.                         &
                   (ND (:,IMODE) >= 0.0E0))
        LOGIC2(:)=((MDT(:,IMODE) <  1.0E20).AND.                        &
                   (MDT(:,IMODE) >  0.0E0))
        LOGIC12(:)=((.NOT.LOGIC1(:)).OR.(.NOT.LOGIC2(:)))
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          DO JL=1,NBOX
           LOGIC3=((MD(JL,IMODE,ICP) < 1.0E20).AND.                     &
                 (MD(JL,IMODE,ICP) >= 0.0))
           IF(LOGIC12(JL).OR.(.NOT.LOGIC3)) THEN
            write(6,*) PROCESS
            write(6,'(3i4,4e13.5)') JLAT,IMODE,ICP,                     &
         ND(JL,IMODE),MDT(JL,IMODE),MD(JL,IMODE,ICP),                   &
        (ND(JL,IMODE)*MD(JL,IMODE,ICP))
            TOTERR=TOTERR+1
           ENDIF ! if bad value for either ND,MDT or MD
          ENDDO ! loop over boxes
         ENDIF ! COMPONENT(IMODE,ICP)
        ENDDO
       ENDIF ! MODE(IMODE)
      ENDDO ! loop over modes
!      IF(TOTERR > 0) THEN
!       cmessage='ERRORS!!!!! STOPPING.'
!       write(6,*) cmessage
!!DEPENDS ON: ereport
!       call ereport('UKCA_CHECK_MD_ND',1,cmessage)
!      ENDIF

      RETURN
      END SUBROUTINE UKCA_CHECK_MD_ND
#endif
