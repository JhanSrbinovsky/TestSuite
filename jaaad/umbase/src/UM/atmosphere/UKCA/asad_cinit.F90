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
! Purpose: To initialize variables used in the chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTR_CTL before call to 
!          ASAD_CDRIVE
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Method
!     ------
!     Input arguments are checked and copied to common. Species
!     and reaction data are read in. Other variables used in the
!     chemistry are initialised.
!
!     -- Calculation of peps
!
!     At several places in the ASAD code, we needed a min. value
!     to use to guard against zero divides. We have to compute
!     this to allow for all possible computer hardwares and precisions
!     that ASAD might be run at. [The LAPACK routine slamch.f returns
!     the safe minimum value such that 1/sfmin would not overflow.
!     Since peps is used in ftoy in calcs like y1/y2 where y2 is
!     small we have to ensure that the result never overflows.
!     Therefore, we use the max. likely value of y1 to scale sfmin
!     returned by slamch.] Typically, the max. total no. density is
!     about 1e19 in the troposphere.
!
!     Externals
!     ---------
!     inrats      - Reads species and reaction data.
!     inphot      - Initialises photolysis scheme.
!     ukca_inwdep - Initialises wet deposition data
!                   (user supplied routine)
!     ukca_inddep - Initialises dry deposition data
!                   (user supplied routine)
!     inemit      - Initialises emission data
!                   (user supplied routine).
!     slamch      - BLAS routine to return machine constants.
!                   (NO LONGER USED, see peps)
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_CINIT(p_field)

        USE ASAD_MOD
        IMPLICIT NONE

#include "parparm.h"
#include "typsize.h"
#include "cntlatm.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: p_field

!       Local variables

        INTEGER :: j                      ! Loop variable
        INTEGER :: jc                     ! Loop variable
        INTEGER :: jf                     ! Loop variable
        INTEGER :: jg                     ! Loop variable
        INTEGER :: jl                     ! Loop variable
        INTEGER :: jp                     ! Loop variable
        INTEGER :: jr                     ! Loop variable
        INTEGER :: js                     ! Loop variable
        INTEGER :: jtr                    ! Loop variable
        INTEGER :: jx                     ! Loop variable
        INTEGER :: jpnpx3                 ! Loop variable
        INTEGER, PARAMETER :: nrsteps_max = 200  ! max steps
        INTEGER, PARAMETER :: nit0_max    = 50   ! max
        INTEGER, PARAMETER :: nitfg_max   = 50   ! max
        INTEGER, PARAMETER :: nitnr_max   = 50   ! max

        REAL :: sfmin

        CHARACTER (LEN=72) :: cmessage    ! Error message

!       Initialisations for COMMON variables.

        ljacx = .true.
        
        jpnpx3 = (jpnr/(3*3))+3*3


!       1.  Copy input arguments to COMMON and check them.
!           ---- ----- --------- -- ------ --- ----- -----

!       Logical arguments
        DO j = 1, jpspec
          lemit(j) = .false.
        ENDDO

        IF (nrsteps < 0 .OR. nrsteps > 200) THEN
          cmessage = ' NRSTEPS IS OUT OF RANGE, RESETTING'
          write(6,*) 'NRSTEPS = ',nrsteps,' Reset to: ',nrsteps_max
          nrsteps = nrsteps_max
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_CINIT',-1,cmessage)
        ENDIF

        IF (nit0 < 0 .OR. nit0 > 50) THEN
          cmessage = ' NIT0 IS OUT OF RANGE, RESETTING'
          write(6,*) 'NIT0 = ',nit0,' Reset to: ',nit0_max
          nit0 = nit0_max
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_CINIT',-1,cmessage)
        ENDIF

        IF (nitfg < 0 .OR. nitfg > 50) THEN
          cmessage = ' NITFG IS OUT OF RANGE, RESETTING'
          write(6,*) 'NITFG = ',nitfg,' Reset to: ',nitfg_max
          nitfg = nitfg_max
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_CINIT',-1,cmessage)
       ENDIF

       IF (nitnr < 0 .OR. nitnr > 50 )then
         cmessage = ' NITNR IS OUT OF RANGE, RESETTING'
         write(6,*) 'NITNR = ',nitnr,' Reset to: ',nitnr_max
         nitnr = nitnr_max
! DEPENDS ON: ereport
         CALL EREPORT('ASAD_CINIT',-1,cmessage)
       ENDIF

!      2.  Assign common values
!          ------ ------ ------

        ncsteps = int( dtime / float(kcdt) )
        cdt     = float(kcdt)
        IF (ncsteps*cdt /= dtime) THEN
          ncsteps = ncsteps + 1
          cdt     = dtime / float(ncsteps)
        ENDIF
        kcdt = cdt

!       2.1  Set photolysis frequency.

        IF (kfphot < 0 .AND. abs(kfphot) > dtime) THEN
          write (6,*) '**CINIT WARNING: VALUE OF KFPHOT ',kfphot,      &
          'EXCEEDS THE MODEL TIMESTEP. ROUTINE PHOTOL WILL ONLY',      &
          ' BE CALLED ONCE.'
          nfphot = 0
        ELSEIF ( kfphot > 0 .AND. kfphot > ncsteps ) THEN
          write (6,*) '**CINIT WARNING: FREQUENCY KFPHOT ',kfphot,     &
           ' EXCEEDS THE TOTAL NUMBER OF CHEMICAL SUBSTEPS. ROUTINE ', &
           ' PHOTOL WILL BE CALLED ONCE ONLY.'
          nfphot = 0
        ELSEIF (kfphot < 0) THEN
          nfphot = int( abs(kfphot)/cdt )
        ELSE
          nfphot = kfphot
        ENDIF

!       2.2  Compute minimum safe value (see Method above)

        sfmin = tiny(1.0d0)
        sfmin = 10.0**(int(log10(sfmin))+1)
        peps  = 1.0e19 * sfmin

!       3.  Set fixed vmrs

        fco2 = 350.0e-6
        fh2  = 5.0e-7
        fn2  = 0.78084
        fo2  = 0.20945
        fch4 = 1.76e-6

!       4.  Clear the species arrays

        f      = 0.0
        fdot   = 0.0
        ej     = 0.0
        linfam = .false.
        linfam = .false.

        y    = 0.0
        ydot = 0.0
        prod = 0.0
        slos = 0.0
        dpd  = 0.0
        dpw  = 0.0
        emr  = 0.0

!       5.   Clear the rates and index arrays.
!            ----- --- ----- --- ----- -------


        rk   = 0.0
        prk  = 0.0
        nspi = 0

        DO js = 1, jpspec
          ngrp(js,1)          = 0
          ngrp(js,2)          = 0
          ngrp(js,3)          = 0
          nprdx2(1,js)        = 0
          nprdx2(2,js)        = 0
          nprdx1(js)          = 0
          ngrp(js+jpspec,1)   = 0
          ngrp(js+jpspec,2)   = 0
          ngrp(js+jpspec,3)   = 0
          nprdx2(1,js+jpspec) = 0
          nprdx2(2,js+jpspec) = 0
          nprdx1(js+jpspec)   = 0
          nlall(js)           = 0
          nlstst(js)          = 0
          nlf(js)             = 0
          nlmajmin(js)        = 0
          nldepd(js)          = 0
          nldepw(js)          = 0
          nlemit(js)          = 0
          nldepx(js)          = 0
        ENDDO

        DO js = 1, 2*jpspec
          DO jx = 1, jpnpx3
            nprdx3(1,jx,js) = 0
            nprdx3(2,jx,js) = 0
            nprdx3(3,jx,js) = 0
          ENDDO
        ENDDO

        nbrkx = 0
        ntrkx = 0
        nprkx = 0
        nhrkx = 0

        njacx3(1,:,:) = 0
        njacx3(2,:,:) = 0
        njacx3(3,:,:) = 0

        njcgrp(:,1) = 0
        njcgrp(:,2) = 0
        njcgrp(:,3) = 0
        njacx2(1,:) = 0
        njacx2(2,:) = 0
        njacx1(:)   = 0
        nltrf(:)    = 0
        nltr3(:)    = 0

        DO jc = 1, jpctr
          nmpjac(jc) = 0
          DO jp = 1, jppjac
            npjac1(jp,jc) = 0
          ENDDO
        ENDDO

        frpb  = 0.0
        nfrpx = 0

        ntabfp(:,1) = 0
        ntabfp(:,2) = 0
        ntabfp(:,3) = 0
        nmzjac = 0
        nmsjac = 0
        nzjac1 = 0
        nsjac1 = 0
        ntabpd = 0
        ztabpd = 0.0
        npdfr  = 0

!       6.  Read chemistry data
!           ---- --------- ----

! DEPENDS ON: asad_inrats
        CALL asad_inrats

! Check that deposition and emission is not on for constant species
        DO js = 1, jpspec
          IF ( ldepd(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Dry deposition turned on for constant species'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( ldepw(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Wet deposition turned on for constant species'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( lemit(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Emission turned on for constant species'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
        ENDDO


        IF ( method == 3) THEN   ! For Newton-Raphson solver only
! DEPENDS ON: asad_setsteady
          CALL asad_setsteady    ! Initialize steady-state species
        ENDIF

        IF ( method >= 10 ) THEN
          DO j = 1, jpspec
            IF ( ctype(j)  ==  jpfm .or. ctype(j)  ==  jpif ) THEN
              WRITE(6,*) '*** ASAD ERROR: You cannot use families ',   &
            ' with one of the stiff integrators. If method  >=  10 ',  &
            ' you cannot have species specified as ',jpfm,' or ',jpif
              cmessage = 'ASAD ABORTED'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_CINIT',j,cmessage)
            ENDIF
          ENDDO
        ENDIF

!       7.  Set up the index arrays.
!           --- -- --- ----- -------

! DEPENDS ON: asad_inix
        CALL asad_inix
! DEPENDS ON: asad_inijac
        CALL asad_inijac

!       8.  Initialise photolysis and heterogeneous chemistry
!           ---------- ---------- --- ------------- ---------

! These are dummy routines at vn7.0 of the UM
!! DEPENDS ON: asad_inphot
!        CALL asad_inphot
!! DEPENDS ON: asad_inhet
!        CALL asad_inhet

!       9.  Read deposition and emission data
!           ---- ---------- --- -------- ----

! DEPENDS ON: ukca_inwdep
        IF ( ndepw /= 0 ) CALL UKCA_INWDEP
! DEPENDS ON: ukca_inddep
        IF ( ndepd /= 0 ) CALL UKCA_INDDEP
! DEPENDS ON: asad_inemit
        IF ( nemit /= 0 ) CALL asad_inemit

        RETURN
        END SUBROUTINE ASAD_CINIT
#endif
