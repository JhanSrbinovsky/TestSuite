#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM: Select a new decomposition
!
! Subroutine Interface:
      SUBROUTINE CHANGE_DECOMPOSITION(new_decomp,icode)

      IMPLICIT NONE

!
! Description:
! Sets up the PARVARS common blocks with the correct information for
! decomposition new_decomp
!
! Method:
! If new_decomp is already the current decomposition, exit and do
! nothing.
! If decomposition new_decomp has not been initialised, print a
! message and exit, with icode=-1.
! Otherwise, copy the information from the decomp_db arrays in the
! DECOMPDB comdeck into the PARVARS comdecks arrays.
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      22/8/96  New deck created for MPP code.  P.Burton
!  4.3      17/02/97 Changed ICODE to a positive error no. P.Burton
!  5.0      12/04/99 - added new halosize PARVARS array and halo_i/j
!                    - glsize now has extra Nfld_max dimension
!                    - lasize now has extra NHalo_max and Nfld_max
!                      dimensions
!                    - blsizep/u replace by blsize with new
!                      Nfld_max dimension
!                    - attop etc replaced with at_extremity
!                                                         P.Burton
!  5.1      02/02/00 Code for new g_pe_index variables   P.Burton
!  5.2      02/08/00 Code for new g_at_extremity variable  P.Burton
!  5.3      14/09/01 Added sb_model_domain variable    P.Burton
!  5.3      22/11/01 Enable MPP as the only option for
!                    small executables         E.Leung
!  5.5      06/08/00 Modification for parallelisation of WAM
!                    Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  5.5      30/01/03 River routing support. P.Selwood
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!  6.2      23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! Subroutine arguments:

      INTEGER                                                           &
     &  new_decomp                                                      &
                     ! IN : new decomposition to use
     &, icode        ! OUT: return code (-1 is failure)

! Parameters and Common blocks

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"

! Local variables
      INTEGER ineb,idim,ifld,iproc,ihalo,ipt,iside


! ------------------------------------------------------------------


! Check that the new_decomp argument is sensible
      IF ((new_decomp  >   max_decomps) .OR.                            &
     &   ((new_decomp  <   1) .AND. (new_decomp  /=  decomp_unset)))    &
     & THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Error: Cannot change to decomposition ',          &
     &               new_decomp
          WRITE(6,*) 'This decomposition does not exist'
          WRITE(6,*) 'Exiting.'
        ENDIF
        icode=1
        GOTO 999
      ENDIF

! Check if this is already the current decomposition

      IF (new_decomp  ==  current_decomp_type) GOTO 999

! Check to see if setting decomposition to unset

      IF (new_decomp  ==  decomp_unset) THEN
        current_decomp_type = decomp_unset
        GOTO 999
      ENDIF

! Check if this decomposition has been initialised

      IF ( .NOT. decomp_db_set(new_decomp) ) THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Error : Attempt to select uninitialised ',        &
     &               'decomposition ',new_decomp
          WRITE(6,*) 'Exiting.'
        ENDIF
        icode=1
        GOTO 999
      ENDIF

! Now we can copy the information into PARVARS

      first_comp_pe=decomp_db_first_comp_pe(new_decomp)
      last_comp_pe=decomp_db_last_comp_pe(new_decomp)

#if !defined(UTILIO) && !defined(FLDIO) \
 && !defined(FLUXPROC)
      nproc=decomp_db_nproc(new_decomp)
#endif
      nproc_x=decomp_db_gridsize(1,new_decomp)
      nproc_y=decomp_db_gridsize(2,new_decomp)

      sb_model_domain=decomp_db_sb_model_domain(new_decomp)

      DO ihalo=1,NHalo_max
        DO idim=1,Ndim_max
          halosize(idim,ihalo)=                                         &
     &      decomp_db_halosize(idim,ihalo,new_decomp)
        ENDDO
      ENDDO

      Offx=decomp_db_halosize(1,halo_type_single,new_decomp)
      Offy=decomp_db_halosize(2,halo_type_single,new_decomp)

      halo_i=decomp_db_halosize(1,halo_type_extended,new_decomp)
      halo_j=decomp_db_halosize(2,halo_type_extended,new_decomp)

      gc_proc_row_group=decomp_db_gc_proc_row_group(new_decomp)
      gc_proc_col_group=decomp_db_gc_proc_col_group(new_decomp)
      gc_all_proc_group=decomp_db_gc_all_proc_group(new_decomp)

      DO ineb=1,4
        neighbour(ineb)=decomp_db_neighbour(ineb,new_decomp)
      ENDDO

      DO idim=1,Ndim_max
        bound(idim)=decomp_db_bound(idim,new_decomp)
        gridsize(idim)=decomp_db_gridsize(idim,new_decomp)

        DO ifld=1,Nfld_max
          glsize(idim,ifld)=decomp_db_glsize(idim,ifld,new_decomp)
        ENDDO
      ENDDO

      DO iproc=first_comp_pe,last_comp_pe
        DO idim=1,Ndim_max
          g_datastart(idim,iproc)=                                      &
     &      decomp_db_g_datastart(idim,iproc,new_decomp)
          g_gridpos(idim,iproc)=                                        &
     &      decomp_db_g_gridpos(idim,iproc,new_decomp)
        ENDDO
        g_at_extremity(PNorth,iproc)=                                   &
     &   (g_gridpos(2,iproc)  ==  (gridsize(2)-1))
        g_at_extremity(PSouth,iproc)=(g_gridpos(2,iproc)  ==  0)
        g_at_extremity(PEast,iproc)=                                    &
     &   (g_gridpos(1,iproc)  ==  (gridsize(1)-1))
        g_at_extremity(PWest,iproc)=(g_gridpos(1,iproc)  ==  0)

        DO ihalo=1,NHalo_max
          DO ifld=1,Nfld_max
            DO idim=1,Ndim_max
              g_lasize(idim,ifld,ihalo,iproc)=                          &
     &          decomp_db_g_lasize(idim,ifld,ihalo,iproc,new_decomp)
            ENDDO
          ENDDO
        ENDDO

        DO ifld=1,Nfld_max
          DO idim=1,Ndim_max
            g_blsize(idim,ifld,iproc)=                                  &
     &          decomp_db_g_blsize(idim,ifld,iproc,new_decomp)

            g_datastart_f(idim,ifld,iproc)=                             &
     &          decomp_db_g_datastart_f(idim,ifld,iproc,new_decomp)
          ENDDO
        ENDDO
      ENDDO

      DO idim=1,Ndim_max
        DO ifld=1,Nfld_max

          DO ihalo=1,NHalo_max
            lasize(idim,ifld,ihalo)=                                    &
     &        g_lasize(idim,ifld,ihalo,mype)
          ENDDO ! ihalo

          blsize(idim,ifld)=g_blsize(idim,ifld,mype)
          datastart_f(idim,ifld)=g_datastart_f(idim,ifld,mype)

        ENDDO ! ifld

        datastart(idim)=g_datastart(idim,mype)
        gridpos(idim)=g_gridpos(idim,mype)

      ENDDO ! idim


      DO iside=1,4
        at_extremity(iside)=g_at_extremity(iside,mype)
      ENDDO

      DO ipt=1-decomp_db_halosize(1,halo_type_extended,new_decomp),     &
     &       decomp_db_glsize(1,fld_type_p,new_decomp)+                 &
     &         decomp_db_halosize(1,halo_type_extended,new_decomp)
        g_pe_index_EW(ipt)=decomp_db_g_pe_index_EW(ipt,new_decomp)
      ENDDO

      DO ipt=1-decomp_db_halosize(2,halo_type_extended,new_decomp),     &
     &       decomp_db_glsize(2,fld_type_p,new_decomp)+                 &
     &         decomp_db_halosize(2,halo_type_extended,new_decomp)
        g_pe_index_NS(ipt)=decomp_db_g_pe_index_NS(ipt,new_decomp)
      ENDDO

      current_decomp_type=new_decomp

 999  CONTINUE

      RETURN
      END SUBROUTINE CHANGE_DECOMPOSITION

#endif
