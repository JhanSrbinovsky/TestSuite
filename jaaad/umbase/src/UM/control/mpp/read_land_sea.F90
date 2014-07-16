#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM : Reads in the local section of Land-Sea Mask.
!
! Subroutine Interface:
      SUBROUTINE READ_LAND_SEA(NFT,IOSTAT,LOOKUP,LOC_LEN1_LOOKUP,       &
     &                         LOC_LEN2_LOOKUP,FIXHD,LOC_LEN_FIXHD)

#if defined(OASIS3)
      USE oasis3_atm_data_mod
#endif
#if defined(OASIS4)
      USE oasis4_atm_data_mod
#endif

      IMPLICIT NONE
!
! Description:
!  This routine reads the land-sea mask (LSM) from the dump and puts
!  it in a COMMON block defined in IOVARS. It is required for
!  unpacking and packing fields which are stored compressed to
!  land points.
!
! Method:
!  The position of the LSM within the dump is found from examining
!  the LOOKUP headers, it is then read in, and the relevant part
!  of the field sent to each processor. The local number of land
!  points is counted, and the LAND_FIELD variable is reset to this
!  new value.
!  Note : Halos can contain land points - but only those halos
!         which are updated by SWAPBNDS.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    4/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + D.Salmond
!    4.1    18/3/96   Simplified communications    P.Burton
!    4.2    18/11/96  Added *CALL AMAXSIZE for IOVARS
!                     Added atmos_ prefix to landmask fields P.Burton
!    4.2    16/8/96   Add IOSTAT argument to SETPOS_SINGLE and
!                     check return code.                    P.Burton
!    4.2    17/10/96 New name for group of processors in gather_field
!                    P.Burton
!    4.3    11/03/97 Corrected calculation of global LAND_FIELD
!                    Store full global LSM on each PE.  P.Burton
!    4.4    25/04/97 Changes to read well-formed records if the
!                    input dumpfile is in that format (almost PP file
!                    format)
!                      Author: Bob Carruthers, Cray Research
!    4.5    13/01/98 Removed reference to SHMEM COMMON block  P.Burton
!    4.5    15/04/98 Modify output. D. Robinson.
!    5.0    26/07/99 Modified PARVARS variables for 5.0       P.Burton
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!    5.3    05/12/01 Replace reference to MaxFieldSize with the variable
!                    Max2DFieldSize in the error-writing statements -
!                    i.e. tidy up the code.                     S.Cusack
!    5.4    03/09/02 SX now uses this deck - add def UTILIO
!                    Calculate atmos_landmask_local for SX
!                                                   E.Leung
!    5.5    28/01/03 Reverse the above changes          E.Leung
!    5.5    02/08/00 Modification for parallelisation of WAM.
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!                    Includes portability changes allowing for
!                    big_endian I/O on little_endian platforms. P.Dando
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!    6.0    11/09/03 Replaced call to ABORT with one to EREPORT. P.Dando
!    6.2    23/11/05 Removed all references to the wavemodel.
!                    T.Edwards
!
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFT                                                             &
                         ! IN : FORTRAN unit number
     & ,LOC_LEN1_LOOKUP                                                 &
                         ! IN : Dimension of the LOOKUP array
     & ,LOC_LEN2_LOOKUP                                                 &
                         ! IN : Dimension of the LOOKUP array
     & ,LOC_LEN_FIXHD    ! IN : Dimension of the FIXHD array

      INTEGER                                                           &
     &  LOOKUP(LOC_LEN1_LOOKUP,LOC_LEN2_LOOKUP),                        &
!                        ! IN : LOOKUP array from dump header
     &  FIXHD(LOC_LEN_FIXHD) ! IN : FIXHD array from dump header

      REAL                                                              &
     &  IOSTAT           ! OUT : Return code

! Parameters and Common blocks

#include "clookadd.h"
#include "parvars.h"
#include "typsize.h"
#include "cntl_io.h"
#include "c_mdi.h"
#include "atm_lsm.h"
#include "gccom.h"
#include "cntlatm.h"

! Local variables

      INTEGER i,j,k,word_address,ipts,iproc,info,len_io,                &
     &        landpts_local,local_off,global_off

      INTEGER :: GATHER_PE
      INTEGER :: JROWSU, JROWSV
      INTEGER :: OASIS_ROWSU, OASIS_ROWST, OASIS_ROWSV

      LOGICAL::GMASKT(glsize(1,fld_type_p),glsize(2,fld_type_p))

      INTEGER                      ::  ErrorStatus = 0
      CHARACTER (Len=80)           ::  CMessage = ' '

      real io_ret(3)

! Local Parameters
      CHARACTER (Len=*), Parameter ::  RoutineName = 'READ_LAND_SEA'

! --------------------------------------------------------------------

      IOSTAT=-1.0
      ipts=0
      len_io=0

! Find location of LSM in the dump

      IF (mype  ==  0) THEN

        DO i=1,LOC_LEN2_LOOKUP
          IF (LOOKUP(ITEM_CODE,i)  ==  30) GOTO 100
        ENDDO
        write(6,'(/''Error in READ_LAND_SEA_MASK:'',                    &
     &   '' Missing Field of Type '',i4)')                              &
     &   item_code
        ErrorStatus = 1
        CMessage = 'Missing Field or Wrong Dumpfile'
! DEPENDS ON: ereport
        CALL EReport(RoutineName,ErrorStatus,CMessage)

100     CONTINUE

        k=i
        word_address=1
! Old Format dumpfiles
        if((lookup(lbnrec,k) == 0) .or.                                 &
! Prog lookups in dump before vn3.2:
     &    ((lookup(lbnrec,k) == imdi) .and. (fixhd(12) <= 301))) then
! Dump and ancillary files
        word_address=1
        IF (i  >   1) THEN
          DO k=2,i
          IF(MOD((LOOKUP(LBPACK,k-1)),10) == 2) THEN
              ipts=(LOOKUP(LBLREC,k-1)+1)/2
          ELSE
              ipts=(LOOKUP(LBLREC,k-1))
          ENDIF
            word_address=word_address+ipts
          ENDDO
        ENDIF
        word_address=FIXHD(160)+word_address-2
          ipts=lookup(lblrec, i)
        else
! PP type files and new format Dumpfiles (vn4.4 onwards)
          word_address=lookup(lbegin,i)
! Use the stored round-up value
          ipts=lookup(lbnrec,i)
        endif

        CALL SETPOS_SINGLE(NFT,word_address,IOSTAT)
        IF (IOSTAT  /=  0) THEN
          WRITE(6,*) 'READ_LAND_SEA: Error Return from SETPOS_SINGLE',  &
     &               ' Status is ',IOSTAT

          ErrorStatus = 3
          CMessage = 'Error from SETPOS_SINGLE - see output for Status'
! DEPENDS ON: ereport
          CALL EReport(RoutineName,ErrorStatus,CMessage)

        ENDIF

! Read the LSM in to PE 0


!--check that there is space to read the data
        if(ipts >  Max2DFieldSize) then
          write(6,9921) ipts, Max2DFieldSize, lookup(lblrec, i)
9921      format(/'READ_LAND_SEA_MASK: The number of Words',            &
     &     ' to be Read ',i10,' is larger than the Buffer Size ',       &
     &     i10//,'Record length is ',i10/)
#if defined(MPP) && defined(T3E)
          if(my_pe() == 0)                                              &
     &     write(0,9921) ipts, Max2DFieldSize, lookup(lblrec, i)
#endif

          ErrorStatus = 4
          CMessage = 'Insufficient Space for Land Sea Mask'
! DEPENDS ON: ereport
          CALL EReport(RoutineName,ErrorStatus,CMessage)

        endif
!
        call buffin_single(nft,atmos_landmask,ipts,                     &
     &                     len_io,IOSTAT)

      ENDIF   ! (mype == 0)

! Broadcast the I/O Status to the other PE's

      io_ret(1)=len_io
      io_ret(2)=iostat
      io_ret(3)=ipts
!
      call gc_rbcast(99, 3, 0, nproc, info, io_ret)
!
      len_io=nint(io_ret(1))
      iostat=io_ret(2)
      ipts=io_ret(3)

! Check the I/O Status on all PE'e

      IF ((IOSTAT /= -1.0).OR.(LEN_IO /= IPTS)) THEN
        WRITE(6,*)'ERROR READING DUMP ON UNIT ',NFT
! DEPENDS ON: ioerror
        CALL IOERROR('BUFFER IN FROM READ_LAND_SEA_MASK',               &
     &   IOSTAT,LEN_IO,IPTS)
        call abort()
      END IF

! Broadcast the global LSM to all processors

      CALL GC_IBCAST(100,glsize(1,fld_type_p)*glsize(2,fld_type_p),     &
     &               0,nproc,info,atmos_landmask)


      DO i=1,lasize(1,fld_type_p,halo_type_no_halo)*                    &
     &       lasize(2,fld_type_p,halo_type_no_halo)
          atmos_landmask_local(i)=.FALSE.
      ENDDO

! Copy my local part of the full LSM into atmos_landmask_local
#if defined(OASIS3) || defined(OASIS4)

       ALLOCATE (UM_TMASK(                                              &
     & lasize(1,fld_type_p,halo_type_no_halo)                           &
     &,lasize(2,fld_type_p,halo_type_no_halo),1))

       ALLOCATE (UM_UMASK(                                              &
     & lasize(1,fld_type_u,halo_type_no_halo)                           &
     &,lasize(2,fld_type_u,halo_type_no_halo),1))

       ALLOCATE (UM_VMASK(                                              &
     & lasize(1,fld_type_v,halo_type_no_halo)                           &
     &,lasize(2,fld_type_v,halo_type_no_halo),1))

#endif

      DO j=1,blsize(2,fld_type_p)

        local_off=(j-1)*lasize(1,fld_type_p,halo_type_no_halo)
        global_off=(j-1+datastart(2)-1)*glsize(1,fld_type_p)+           &
     &               datastart(1)-1

        DO i=1,blsize(1,fld_type_p)

          atmos_landmask_local(local_off+i)=                            &
     &      atmos_landmask(global_off+i)

#if defined(OASIS3) || defined(OASIS4)
          UM_TMASK(I,J,1) = atmos_landmask_local(local_off+i)
#endif
        ENDDO ! i
      ENDDO ! j


#if defined(OASIS3) || defined(OASIS4)
      ! Always use PE 0 as the OASIS master
      Oasis_CNTLPE= 0
      Oasis_NPROC = NPROC
      Oasis_PE = MYPE

      ! Gather all local domains of our LS mask
      ! to create a global array.
      GATHER_PE = Oasis_CNTLPE

      ! 1st: Gather all onto 1 PE
! DEPENDS ON: gather_field
      CALL GATHER_FIELD(UM_TMASK,GMASKT,                                &
     &  lasize(1,fld_type_p,halo_type_no_halo),                         &
     &  lasize(2,fld_type_p,halo_type_no_halo),                         &
     &  glsize(1,fld_type_p),                                           &
     &  glsize(2,fld_type_p),                                           &
     &  fld_type_p,halo_type_no_halo,                                   &
     &  gather_pe,GC_ALL_PROC_GROUP,info,cmessage)

      ! ALLOCATE space for global U and V masks
      ! We allocate space on all PEs even though
      ! these arrays are only needed on the gather PE
      ! just because it saves mucking about
      ! with if tests when we deallocate them.
      ALLOCATE(GMASKU(glsize(1,fld_type_u),                             &
     &                  glsize(2,fld_type_u)))
      ALLOCATE(GMASKV(glsize(1,fld_type_v),                             &
     &                  glsize(2,fld_type_v)))

      ! Only do the global mask calculations
      ! on the gathering PE.
      IF (MYPE == GATHER_PE) THEN

       GMASKU(:,:) = .FALSE.
       GMASKV(:,:) = .FALSE.

       JROWSU = glsize(2,fld_type_u)
       JROWSV = glsize(2,fld_type_v)

! set up logical land/sea mask on atmosphere U & V grids:
       DO J=1,JROWSU-1
         DO I=1,ROW_LENGTH-1
           GMASKU(I,J)=GMASKT(I,J).OR.                                  &
     &                 GMASKT(I+1,J).OR.                                &
     &                 GMASKT(I,J+1).OR.                                &
     &                 GMASKT(I+1,J+1)
         ENDDO
         GMASKU(ROW_LENGTH,J)=GMASKT(ROW_LENGTH,J).OR.                  &
     &                        GMASKT(1,J).OR.                           &
     &                        GMASKT(ROW_LENGTH,J+1).OR.                &
     &                        GMASKT(1,J+1)
       ENDDO
! N.pole :
       DO I=1,ROW_LENGTH-1
          GMASKU(I,jrowsu)=GMASKT(I,jrowsu).OR.                         &
     &                     GMASKT(I+1,jrowsu)
       ENDDO
       GMASKU(row_length,jrowsu)=GMASKT(row_length,jrowsu)              &
     &                      .OR. GMASKT(1,jrowsu)
!
       DO J=1,JROWSV
         DO I=1,ROW_LENGTH-1
            GMASKV(I,J)=GMASKT(I,J).OR.                                 &
     &                  GMASKT(I+1,J).OR.                               &
     &                  GMASKT(I,J+1).OR.                               &
     &                  GMASKT(I+1,J+1)
         ENDDO
         GMASKV(ROW_LENGTH,J)=GMASKT(ROW_LENGTH,J).OR.                  &
     &                        GMASKT(1,J).OR.                           &
     &                        GMASKT(ROW_LENGTH,J+1).OR.                &
     &                        GMASKT(1,J+1)
       ENDDO

      ENDIF ! GATHER_PE

      IF (L_COUPLE_MASTER) THEN
         IF (MYPE == OASIS_CNTLPE) THEN
            OASIS_ROWST = glsize(2,fld_type_p)
            OASIS_ROWSU = glsize(2,fld_type_u)
            OASIS_ROWSV = glsize(2,fld_type_v)
         ELSE
            ! All non master PEs just get a dimension of 1
            OASIS_ROWST = 1
            OASIS_ROWSU = 1
            OASIS_ROWSV = 1
         ENDIF
      ELSE
         ! Local domain sizes.
         OASIS_ROWST = lasize(2,fld_type_p,halo_type_no_halo)
         OASIS_ROWSU = lasize(2,fld_type_u,halo_type_no_halo)
         OASIS_ROWSV = lasize(2,fld_type_v,halo_type_no_halo)
      ENDIF

      ALLOCATE(OASIS_TMASK                                              &
     &   (lasize(1,fld_type_p,halo_type_no_halo),                       &
     &    OASIS_ROWST,1))

      ALLOCATE(OASIS_UMASK                                              &
     &   (lasize(1,fld_type_u,halo_type_no_halo),                       &
     &    OASIS_ROWSU,1))

      ALLOCATE(OASIS_VMASK                                              &
     &   (lasize(1,fld_type_v,halo_type_no_halo),                       &
     &    OASIS_ROWSV,1))

#if defined(OASIS4)
         IF (L_COUPLE_MASTER) THEN
            IF (MYPE == OASIS_CNTLPE) THEN
               OASIS_TMASK(:,:,1) = GMASKT(:,:)
               OASIS_UMASK(:,:,1) = GMASKU(:,:)
               OASIS_VMASK(:,:,1) = GMASKV(:,:)
            ENDIF
         ELSE


! DEPENDS ON: scatter_field
            CALL SCATTER_FIELD(UM_UMASK,GMASKU,                         &
     &  lasize(1,fld_type_u,halo_type_no_halo),                         &
     &  lasize(2,fld_type_u,halo_type_no_halo),                         &
     &  glsize(1,fld_type_u),                                           &
     &  glsize(2,fld_type_u),                                           &
     &  fld_type_u,halo_type_no_halo,                                   &
     &  gather_pe,GC_ALL_PROC_GROUP,info,cmessage)

! DEPENDS ON: scatter_field
            CALL SCATTER_FIELD(UM_VMASK,GMASKV,                         &
     &  lasize(1,fld_type_v,halo_type_no_halo),                         &
     &  lasize(2,fld_type_v,halo_type_no_halo),                         &
     &  glsize(1,fld_type_v),                                           &
     &  glsize(2,fld_type_v),                                           &
     &  fld_type_v,halo_type_no_halo,                                   &
     &  gather_pe,GC_ALL_PROC_GROUP,info,cmessage)


            ! For OASIS4 purposes, we require that
            ! masks = FALSE over LAND.
            OASIS_TMASK(:,:,:) = .NOT.UM_TMASK(:,:,:)
            OASIS_UMASK(:,:,:) = .NOT.UM_UMASK(:,:,:)
            OASIS_VMASK(:,:,:) = .NOT.UM_VMASK(:,:,:)

         ENDIF ! Coupling through each PE in //
#endif

      ! Deallocate the global U and V masks since we'll
      ! only need the local copies from now on.
      DEALLOCATE (GMASKU)
      DEALLOCATE (GMASKV)

#endif      

! Count the number of global land points

      atmos_number_of_landpts=0
      DO i=1,glsize(1,fld_type_p)*glsize(2,fld_type_p)
        IF (atmos_landmask(i))                                          &
     &      atmos_number_of_landpts=atmos_number_of_landpts+1
      ENDDO

! Do a swap to get land points in halo areas

      landpts_local=0
      DO i=1,lasize(1,fld_type_p,halo_type_no_halo)*                    &
     &       lasize(2,fld_type_p,halo_type_no_halo)
        IF (atmos_landmask_local(i))                                    &
     &    landpts_local=landpts_local+1
      ENDDO


      IF (landpts_local  /=  LAND_FIELD) THEN
        WRITE(6,*) 'PE ',mype,' : LAND_FIELD is being reset from ',     &
     &             LAND_FIELD,' to ',landpts_local
        LAND_FIELD=landpts_local
      ENDIF

      RETURN
      END SUBROUTINE READ_LAND_SEA

#endif
