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
! Purpose: Reads species and reaction data. Combines reactions into one
!          array and reorders them to put single reactant reactions first
!          to improve code in the prls routine.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
!     Interface
!     ---------
!     Reads the files:
!         chch.d   - Chosen chemistry file, contains information
!                    on species involved in the chemistry and
!                    the families/tracers to which they belong.
!         ratb.d   - Data for bimolecular reactions.
!         ratt.d   - Data for trimolecular reactions.
!         ratj.d   - Data for photolysis reactions.
!         rath.d   - Data for heterogeneous reactions.
!     OR gets info from ukca_chem1 module
!
!     Method
!     ------
!     The chch.d file specifies the species types using 2 letter
!     codes for easier reading.
!
!         chch.d file       Meaning
!             'FM'          Family member
!             'FT'          Tracer but will be put into a family
!                           if lifetime becomes short.
!             'TR'          Tracer, advected by calling model.
!             'SS'          Steady state species.
!             'CT'          Constant species.
!
!     Local variables
!     ---------------
!     iadv         Counter of number of model tracers/families.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INRATS

        USE ASAD_MOD
        USE UKCA_CHEM1_DAT
        IMPLICIT NONE

#include "parvars.h"
#include "typsize.h"
#include "cntlatm.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

!       Local variables

        INTEGER :: ispb(jpbk+1,jpspb)
        INTEGER :: ispt(jptk+1,jpspt)
        INTEGER :: ispj(jppj+1,jpspj)
        INTEGER :: isph(jphk+1,jpsph)
        INTEGER :: ifrpbx(jpbk+1)
        INTEGER :: ilmin(jpspec)
        INTEGER :: ilmaj(jpspec)
        INTEGER :: ixdumt(jptk+1)
        INTEGER :: ixdumj(jppj+1)
        INTEGER :: ixdumh(jphk+1)
        INTEGER :: ierror
        INTEGER :: iadv                  ! Counter
        INTEGER :: imajor                ! Counter
        INTEGER :: iminor                ! Counter
        INTEGER :: ix                    ! Counter
        INTEGER :: icount                ! Counter
        INTEGER :: ifam                  ! Index
        INTEGER :: imaj                  ! Index
        INTEGER :: iflag                 ! Used to test family order
        INTEGER :: idummy                ! Dummy variable
        INTEGER :: istat                 ! Tag for communication
        INTEGER :: j                     ! Loop variable
        INTEGER :: jadv                  ! Loop variable
        INTEGER :: jb                    ! Loop variable
        INTEGER :: jctr                  ! Loop variable
        INTEGER :: jh                    ! Loop variable
        INTEGER :: jj                    ! Loop variable
        INTEGER :: jp                    ! Loop variable
        INTEGER :: jr                    ! Loop variable
        INTEGER :: js                    ! Loop variable
        INTEGER :: jspb                  ! Loop variable
        INTEGER :: jsph                  ! Loop variable
        INTEGER :: jspj                  ! Loop variable
        INTEGER :: jspt                  ! Loop variable
        INTEGER :: jt                    ! Loop variable
        INTEGER :: k                     ! Loop variable

        REAL :: zdumt(1)                 ! Dummy array for trimol file
        REAL :: zdumj(1)                 ! Dummy array for photol file
        REAL :: zdumh(1)                 ! Dummy array for heter file

        CHARACTER (LEN=10) :: cmntb(jpbk+1)
        CHARACTER (LEN=10) :: cmntt(jptk+1)
        CHARACTER (LEN=10) :: cmntj(jppj+1)
        CHARACTER (LEN=10) :: cmnth(jphk+1)
        CHARACTER (LEN=72) :: cmessage        ! Error message

        LOGICAL :: gtype3
        LOGICAL :: gdebug
        LOGICAL :: L_exist
        LOGICAL :: L_fa

!       1.  Read chosen chemistry file and determine chemistry
!           ---- ------ --------- ---- --- --------- ---------


        ispb(:,:) = 0
        ifrpbx(:) = 0
        ispt(:,:) = 0
        ispj(:,:) = 0
        isph(:,:) = 0

!       1.1 Read file.

!       Read in data from PE0 and broardcast to other PEs

!        IF (mype == 0) THEN
!          INQUIRE (FILE=CHCH, EXIST=L_exist)
!          IF (.NOT. L_exist) THEN
!              cmessage = 'INRATS: chch.d file does not exist'
!             CALL EREPORT('INRATS',195,cmessage)
!         ENDIF

!          OPEN( jpcio, file = CHCH, form = 'formatted')
!          DO js = 1, jpspec
!            READ (jpcio,*, IOSTAT=ierror) idummy, speci(js), nodd(js),
!                  ctype(js), family(js)
!            IF (ierror /= 0) THEN
!             cmessage = 'Error reading chch.d file '
!             CALL EREPORT('INRATS',js,cmessage)
!           ENDIF

!            IF ( ctype(js) /= jpfm .and. ctype(js) /= jpif .and.
!              ctype(js) /= jpsp .and. ctype(js) /= jpna .and.
!              ctype(js) /= jpco .and. ctype(js) /= jpcf ) then
!              WRITE(6,*) '*** ERROR in reading ',CHCH,' file'
!              WRITE(6,*) 'Species ',speci(js),' has an invalid type ',
!                     ctype(js)
!              cmessage = 'ERROR reading chch.d file'
!             CALL EREPORT('INRATS',js,cmessage)
!            ENDIF

!          ENDDO
!          CLOSE(jpcio)
!        ENDIF

!        CALL GC_IBCAST (1,jpspec,0,nproc,istat,nodd)
!        CALL GC_CBCAST (2,10*jpspec,0,nproc,istat,speci)
!        CALL GC_CBCAST (3,2*jpspec,0,nproc,istat,ctype)
!        CALL GC_CBCAST (4,10*jpspec,0,nproc,istat,family)

!       1.1.1  Set types and find which species are in families.

! Replace read statement by module:

        IF (size(chch_defs) /= jpspec) THEN
          cmessage='jpspec and chch_defs are inconsistent'
          write(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',1,cmessage)
        ENDIF
        DO k=1,jpspec
          speci(k)  = chch_defs(k)%speci
          nodd(k)   = chch_defs(k)%nodd
          ctype(k)  = chch_defs(k)%ctype
          family(k) = chch_defs(k)%family
        ENDDO

        iadv = 0
        ntrf = 0
        ntr3 = 0
        DO js = 1, jpspec
          gtype3 = .false.
          IF ( ctype(js) /= jpfm .AND. ctype(js) /= jpif )             &
            family(js)='          '
          IF ( ctype(js) == jpsp ) THEN
            iadv = iadv + 1
            ntrf = ntrf + 1
            IF ( iadv > jpctr ) THEN
              WRITE (6,*) '** ASAD ERROR in subroutine inrats'
              WRITE (6,*) '** Parameter jpctr is too low; found',iadv
              WRITE (6,*) '** tracers so far with ',jpspec-js
              WRITE (6,*) '** species to check.'
              cmessage = 'ASAD ERROR: jpctr is too low'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_INRATS',iadv,cmessage)
            ENDIF
            advt(iadv)  = speci(js)
            nltrf(ntrf) = iadv
          ELSE IF( ctype(js) == jpfm .OR. ctype(js) == jpif )THEN
            IF ( ctype(js) == jpif ) THEN
              iadv = iadv + 1
              ntr3 = ntr3 + 1
              IF ( iadv  >   jpctr ) THEN
                WRITE (6,*) '** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '** Parameter jpctr is too low; found',iadv
                WRITE (6,*) '** tracers so far with ',jpspec-js
                WRITE (6,*) '** species to check.'
                cmessage = 'ERROR reading chch.d file'
! DEPENDS ON: ereport
                CALL EREPORT('ASAD_INRATS',js,cmessage)
              ENDIF
              advt(iadv)  = speci(js)
              nltr3(ntr3) = iadv
            ENDIF
            l_fa=.true.
            DO jadv = 1, iadv
              IF ( family(js) == advt(jadv) ) THEN
                l_fa=.false.
                EXIT
              ENDIF
            ENDDO
            IF (l_fa) THEN
              iadv = iadv + 1
              ntrf = ntrf + 1
              IF ( iadv  >   jpctr ) then
                WRITE (6,*) '***** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '***** Param jpctr is too low; found',iadv
                WRITE (6,*) '***** tracers so far with ',jpspec-js
                WRITE (6,*) '***** species to check.'
                cmessage = 'INRATS ERROR : jpctr is too low'
! DEPENDS ON: ereport
                CALL EREPORT('ASAD_INRATS',iadv,cmessage)
              ENDIF
              advt(iadv) = family(js)
              nltrf(ntrf) = iadv
            ENDIF      ! l_fa
          ENDIF
        ENDDO

!       1.2 Find major species of families

        DO jadv = 1, iadv
          DO js = 1, jpspec
            IF( family(js) == advt(jadv) .and. ctype(js) /= jpif )     &
              majors(jadv) = js
            IF( speci(js) == advt(jadv) ) majors(jadv) = js
          ENDDO
        ENDDO

!       1.3 Allocate families to species

        DO js = 1, jpspec
          moffam(js) = 0
          madvtr(js) = 0
          DO jadv = 1, iadv
            IF (family(js) == advt(jadv) ) moffam(js) = jadv
            IF (speci(js)  == advt(jadv) ) madvtr(js) = jadv
            IF (family(js) == advt(jadv) .AND. js > majors(jadv)) THEN
              WRITE (6,*) '** ASAD ERROR: '
              WRITE (6,*) 'RE-ORDER SPECIES FILE SO THAT THE MAJOR '
              WRITE (6,*) 'SPECIES OF A FAMILY OCCURS AFTER THE OTHERS'
              cmessage = 'INRATS ERROR : Order of species is incorrect'
! DEPENDS ON: ereport
              CALL EREPORT('ASAD_INRATS',jadv,cmessage)
            ENDIF
          ENDDO
        ENDDO

!       1.4  Build the list of major and minor species

        nlmajmin(1) = 5
        imajor      = 0
        iminor      = 0
        DO js = 1, jpspec
          ifam = moffam(js)
          imaj = 0
          IF ( ifam /= 0 ) THEN
           imaj = majors(ifam)
           IF ( imaj /= js ) THEN
              iminor = iminor + 1
              ilmin(iminor) = js
            ELSE
              imajor = imajor + 1
              ilmaj(imajor) = js
            ENDIF
          ENDIF
        ENDDO
        nlmajmin(2) = nlmajmin(1) + imajor - 1
        nlmajmin(3) = nlmajmin(1) + imajor
        nlmajmin(4) = nlmajmin(3) + iminor - 1
        DO j = nlmajmin(1), nlmajmin(2)
          nlmajmin(j) = ilmaj(j-nlmajmin(1)+1)
        ENDDO
        DO j = nlmajmin(3), nlmajmin(4)
          nlmajmin(j) = ilmin(j-nlmajmin(3)+1)
        ENDDO

        IF ( iadv /= jpctr ) then
          WRITE (6,*) '** ASAD ERROR: Number of advected tracers',     &
           ' specified in the ',CHCH,' file does not match jpctr'
          WRITE (6,*) 'Found ',iadv,' but expected ',jpctr
          cmessage = 'INRATS ERROR : iadv and jpctr do not match'
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',iadv,cmessage)
        ENDIF

!       2.  Write details of chemistry selection to log file
!           ----- ------- -- --------- --------- -- --- ----

        IF (mype == 0) THEN
          WRITE(6,*)
          WRITE(6,*)'  ***  CHEMISTRY INFORMATION  ***'
          WRITE(6,*)
          WRITE(6,*)'ASAD IS TREATING ADVECTED TRACERS IN THE ORDER:'
          WRITE(6,*)
          WRITE(6,'(5(2x,i2,1x,a10))')(jctr,advt(jctr), jctr= 1,jpctr)
          WRITE(6,*)
          WRITE(6,*)'IF THE TRACERS WERE NOT INITIALISED IN THIS '
          WRITE(6,*)'ORDER THEN THE MODEL RESULTS ARE WORTHLESS '
          WRITE(6,*)

          iflag = 0
          DO jctr = 1, jpctr
            IF (advt(jctr) /= speci(majors(jctr)) ) THEN
              IF (iflag == 0 )then
                WRITE(6,*)'THE MAJOR MEMBER OF EACH OF THE FAMILIES'
                WRITE(6,*)'IS GIVEN BELOW. IF THIS IS NOT ACCEPTABLE,'
                WRITE(6,*)'THEN YOU MUST REORDER THE SPECIES IN ',CHCH
                WRITE(6,*)'SO THE MAJOR SPECIES FOLLOWS THE OTHERS.'
                WRITE(6,*)
                iflag = 1
              ENDIF
              WRITE(6,'(a10,1x,a10)') advt(jctr), speci(majors(jctr))
            ENDIF
          ENDDO
        ENDIF     ! End of IF mype statement

!       3.  Bimolecular ratefile
!           ----------- --------

!        gdebug = .false.
!        CALL READRAT( jpcio,jpbk,jpbk+1,jpspb,jpab,jpfrpb,RATB,
!                     .false.,spb,ab,ifrpbx,frpb,cmntb,gdebug )


!       Get bimolecular rates from module

        IF (size(ratb_defs) /= jpbk) THEN 
          cmessage='size of ratb_defs is inconsistent with jpbk'
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',1,cmessage)
        END IF 
        icount=1
        DO k=1,jpbk
          spb(k,1) = ratb_defs(k)%react1
          spb(k,2) = ratb_defs(k)%react2
          spb(k,3) = ratb_defs(k)%prod1
          spb(k,4) = ratb_defs(k)%prod2
          spb(k,5) = ratb_defs(k)%prod3
          ab(k,1)  = ratb_defs(k)%K0
          ab(k,2)  = ratb_defs(k)%alpha
          ab(k,3)  = ratb_defs(k)%beta
          IF (ratb_defs(k)%pyield1 > 1e-18) THEN
            ifrpbx(k)      = 1
            frpb(icount)   = ratb_defs(k)%pyield1
            frpb(icount+1) = ratb_defs(k)%pyield2
            frpb(icount+2) = ratb_defs(k)%pyield3
            icount         = icount+3
          ENDIF
        ENDDO

        DO jb = 1, jpbk
          DO js = 1, jpspec
            DO jspb = 1, jpspb
              IF ( speci(js) == spb(jb,jspb) ) ispb(jb,jspb) = js
            ENDDO
          ENDDO
        ENDDO

!       4.  Trimolecular ratefile
!           ------------ --------

!        CALL READRAT(jpcio,jptk,jptk+1,jpspt,jpat,1,RATT,.true.,
!                     spt,at,ixdumt,zdumt,cmntt,gdebug )


!       Get trimolecular rates from module for UM version

        IF (size(ratt_defs) /= jptk) THEN
          cmessage='size of ratt_defs is inconsistent with jptk'
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',1,cmessage)
        ENDIF
        DO k=1,jptk
          spt(k,1) = ratt_defs(k)%react1
          spt(k,2) = ratt_defs(k)%react2
          spt(k,3) = ratt_defs(k)%prod1
          spt(k,4) = ratt_defs(k)%prod2
          at(k,1)  = ratt_defs(k)%F
          at(k,2)  = ratt_defs(k)%K1
          at(k,3)  = ratt_defs(k)%alpha1
          at(k,4)  = ratt_defs(k)%beta1
          at(k,5)  = ratt_defs(k)%K2
          at(k,6)  = ratt_defs(k)%alpha2
          at(k,7)  = ratt_defs(k)%beta2
        ENDDO

        DO jt = 1, jptk
          DO js = 1, jpspec
            DO jspt = 1, jpspt
              IF (speci(js) == spt(jt,jspt) ) ispt(jt,jspt) = js
            ENDDO
          ENDDO
        ENDDO


!       5.  Photolysis ratefile
!           ---------- --------

!        CALL READRAT(jpcio,jppj,jppj+1,jpspj,jpaj,1,RATJ,.false.,
!                     spj,aj,ixdumj,zdumj,cmntj,gdebug )


!       use module to get spj

        IF (size(ratj_defs) /= jppj) THEN
          cmessage='size of ratj_defs is not equal to jppj'
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',1,cmessage)
        ENDIF
        DO k=1,jppj
          spj(k,1) = ratj_defs(k)%react1
          spj(k,2) = ratj_defs(k)%prod1
          spj(k,3) = ratj_defs(k)%prod2
          spj(k,4) = ratj_defs(k)%prod3
          spj(k,5) = ratj_defs(k)%prod4
        ENDDO

        DO jj = 1, jppj
          DO js = 1, jpspec
            DO jspj = 1, jpspj
              IF (speci(js) == spj(jj,jspj) ) ispj(jj,jspj) = js
            ENDDO
          ENDDO
        ENDDO

!       6.  Heterogeneous ratefile
!           ------------- --------

        IF ( L_ukca_het_psc ) THEN
!          CALL READRAT( jpcio,jphk,jphk+1,jpsph,jpah,1,RATH,.false.,
!                       sph,ah,ixdumh,zdumh,cmnth,gdebug )

!         use module to get sph

          IF (size(ratj_defs) /= jphk) THEN
            cmessage='size of ratj_defs is not equal to jphk'
! DEPENDS ON: ereport
            CALL EREPORT('ASAD_INRATS',1,cmessage)
          ENDIF
          DO k=1,jphk
            sph(k,1) = rath_defs(k)%react1
            sph(k,2) = rath_defs(k)%react2
            sph(k,3) = rath_defs(k)%prod1
            sph(k,4) = rath_defs(k)%prod2
            sph(k,5) = rath_defs(k)%prod3
          ENDDO

          DO jh = 1, jphk
            DO js = 1, jpspec
              DO jsph = 1, jpsph
                IF (speci(js) == sph(jh,jsph) ) isph(jh,jsph) = js
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!       7.  Reorder reactions, putting single reactants first.
!           ------- ---------- ------- ------ --------- ------

        nuni = 0

!       7.1  Single reactants; scan ratefiles in turn.

!
        DO jr = 1, jpbk
          IF (ispb(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nbrkx(jr) = nuni
            DO jp = 1, jpspb
              nspi(nuni,jp) = ispb(jr,jp)
            ENDDO
            IF ( ifrpbx(jr) /= 0 ) nfrpx(nuni) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) == 0 ) THEN
            nuni = nuni + 1
            ntrkx(jr) = nuni
            DO jp = 1, jpspt
              nspi(nuni,jp) = ispt(jr,jp)
            ENDDO
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF ( ispj(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nprkx(jr) = nuni
            DO jp = 1, jpspj
              nspi(nuni,jp) = ispj(jr,jp)
            ENDDO
          ENDIF
        ENDDO

        IF ( L_ukca_het_psc ) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) == 0 ) THEN
              nuni = nuni + 1
              nhrkx(jr) = nuni
              DO jp = 1, jpsph
                nspi(nuni,jp) = isph(jr,jp)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!       7.2  Two reactants; copy remaining reactions

        ix = nuni
        DO jr = 1, jpbk
          IF ( ispb(jr,2) /= 0 ) THEN
            ix = ix + 1
            nbrkx(jr) = ix
            DO jp = 1, jpspb
              nspi(ix,jp) = ispb(jr,jp)
            ENDDO

            IF ( ifrpbx(jr) /= 0 ) nfrpx(ix) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) /= 0 ) THEN
            ix = ix + 1
            ntrkx(jr) = ix
            DO jp = 1, jpspt
              nspi(ix,jp) = ispt(jr,jp)
            ENDDO
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF (ispj(jr,2) /= 0 ) THEN
            ix = ix + 1
            nprkx(jr) = ix
            DO jp = 1, jpspj
              nspi(ix,jp) = ispj(jr,jp)
            ENDDO
          ENDIF
        ENDDO

        IF ( L_ukca_het_psc ) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) /= 0 ) THEN
              ix = ix + 1
              nhrkx(jr) = ix
              DO jp = 1, jpsph
                nspi(ix,jp) = isph(jr,jp)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

        IF ( ix /= jpnr ) THEN
          WRITE (6,*) '*** INTERNAL ASAD ERROR: Number of reactions',  &
                      ' placed in nspi array does not equal jpnr. '
          WRITE (6,*) '                         Check reaction',       &
                      ' files and code.'
          cmessage = 'No of rxns in nspi array is not equal to jpnr'
! DEPENDS ON: ereport
          CALL EREPORT('ASAD_INRATS',ix,cmessage)
        ENDIF

        RETURN
        END SUBROUTINE ASAD_INRATS
#endif
