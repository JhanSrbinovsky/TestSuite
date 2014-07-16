#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine diagnostics_veg ----------------------------------------
!!!
!!! Purpose : Calculates diagnostics for dynamic vegetation and
!!!           outputs them.
!!!
!!! version for CRAY YMP
!!!
!!!  Model            Modification history:
!!! version  Date
!!!  5.2  15/11/00 New deck.  M. Best
!!!  5.3  26/09/01 Correct code for vegetation diagnostics.   M. Best
!  6.2  01/03/06  new diagnostics for RothC carbon pool stores and
!                 fluxes.                                   C.D. Jones
!!!
!!!END -----------------------------------------------------------------
!
! Subroutine diagnostics_veg

      Subroutine diagnostics_veg(                                       &
     &                       row_length, rows, n_rows                   &
     &,                      global_row_length, global_rows             &
     &,                      DIM_CS1, DIM_CS2                           &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      land_pts                                   &
     &,                      land_index                                 &
     &,                      ntype,npft                                 &
     &,                      c_veg,cv,g_leaf_phen                       &
     &,                      lit_c,lit_c_mn,g_leaf_day                  &
     &,                      lai_phen,g_leaf_dr_out,npp_dr_out          &
     &,                      resp_w_dr_out,resp_s_dr_out,frac_disturb   &
     &,                      frac,lai,ht,cs                             &
     &,                                                                 &
#include "argsts.h"
     & STASHwork                                                        &
     &     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!

      Implicit None

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain     ! indicator as to model type, ie global, lam

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)                                       &
     &, land_pts                                                        &
                 ! No.of land points being processed, can be 0.
     &, ntype                                                           &
                    ! Max. No. of land surface tiles
     &, npft                                                            &
                    ! No. of plant functional types
     &, DIM_CS1, DIM_CS2     ! soil carbon dimensions

      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP


! Primary Arrays used in all models
      Integer                                                           &
     &  land_index(land_pts)      ! set from land_sea_mask


      REAL                                                              &
     & C_VEG(land_pts,NPFT)                                             &
                            ! Total carbon content of vegetation
!                              ! (kg C/m2).
     &,CV(land_pts)                                                     &
                            ! Gridbox mean veg carbon (kg C/m2).
     &,LIT_C(land_pts,NPFT)                                             &
                            ! Carbon Litter (kg C/m2/360days).
     &,LIT_C_MN(land_pts)                                               &
                            ! Gridbox mean carbon litter
!                              ! (kg C/m2/360days)
     &,G_LEAF_DAY(land_pts,NPFT)                                        &
                                   ! Mean leaf turnover rate for
!                                     ! input to PHENOL (/360days).
     &,G_LEAF_PHEN(land_pts,NPFT)                                       &
                                   ! Mean leaf turnover rate over
!                                     ! phenology period (/360days).
     &,G_LEAF_DR_OUT(land_pts,NPFT)                                     &
                                   ! Mean leaf turnover rate for
!                                     ! driving TRIFFID (/360days).
     &,LAI_PHEN(land_pts,NPFT)                                          &
                                   ! LAI of PFTs after phenology.
     &,NPP_DR_OUT(land_pts,NPFT)                                        &
                                   ! Mean NPP for driving TRIFFID
!                                     ! (kg C/m2/360days).
     &,RESP_W_DR_OUT(land_pts,NPFT)                                     &
                                   ! Mean wood respiration for
!                                     ! driving TRIFFID
!                                     ! (kg C/m2/360days).
     &,RESP_S_DR_OUT(land_pts,DIM_CS1+1)                                &
                                         ! Mean soil respiration for
!                                     ! driving TRIFFID
!                                     ! (kg C/m2/360days).
     &,FRAC_DISTURB(land_pts)                                           &
                                   ! Fraction of gridbox in which
!                                     !    vegetation is disturbed.
     &,FRAC(land_pts,NTYPE)                                             &
                                   ! Fractions of surface types.
     &,LAI(land_pts,NPFT)                                               &
                                   ! LAI of plant functional
!                                     !       types.
     &,HT(land_pts,NPFT)                                                &
                                   ! Height of plant functional
!                                     !       types (m).
     &,CS(land_pts,DIM_CS1)   ! Soil carbon content
!                                     !       (kg C/m2).

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace


#include "c_mdi.h"

! Local variables

      LOGICAL                                                           &
     & PLLTYPE(NTYPE)                                                   &
                          ! pseudolevel list for surface types
     &,PLLPFT(NPFT)       ! pseudolevel list for PFTs

      INTEGER                                                           &
     & PSLEVEL                                                          &
                     !  loop counter for pseudolevels
     &,PSLEVEL_OUT   !  index for pseudolevels sent to STASH

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_veg')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)

      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

! ----------------------------------------------------------------------
! Section 1.  Initialisation.
! ----------------------------------------------------------------------

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------

!L ITEM 1: VEGETATION CARBON ON PLANT FUNCTIONAL TYPES

      If (sf(1,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,1,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 1 = c_veg)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = c_veg(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(1,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,1,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1901)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(1,19)


!L ITEM 2: GRIDBOX MEAN VEGETATION CARBON

      IF (SF(2,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cv(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(2,19,im_index)),interp_data,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,2,                                             &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1902)"
           goto 9999
        End if

      END IF     !   sf(2,19)

!L ITEM 3: PHENOLOGICAL LEAF TURNOVER RATE PFTS

      If (sf(3,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,3,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 3 = g_leaf_phen)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_phen(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(3,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,3,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1903)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(3,19)


!L ITEM 4: LITTER CARBON ON PLANT FUNCTIONAL TYPES

      If (sf(4,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,4,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 4 = lit_c)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lit_c(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(4,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,4,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1904)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(4,19)


!L ITEM 5: GRIDBOX MEAN LITTER CARBON

      IF (SF(5,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = lit_c_mn(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(5,19,im_index)),interp_data,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,5,                                             &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1905)"
           goto 9999
        End if

      END IF     !   sf(5,19)

!L ITEM 6: MEAN LEAF TURNOVER RATE ON PFTS FOR PHENOLOGY

      If (sf(6,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,6,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 6 = g_leaf_day)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_day(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(6,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,6,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1906)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(6,19)


!L ITEM 7: LEAF AREA INDEX ON PLANT FUNCTIONAL TYPES AFTER PHENOLOGY

      If (sf(7,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,7,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 7 = lai_phen)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lai_phen(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(7,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,7,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1907)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(7,19)


!L ITEM 8: MEAN LEAF TURNOVER RATE ON PFTS FOR TRIFFID

      If (sf(8,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,8,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 8 = g_leaf_dr_out)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = g_leaf_dr_out(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(8,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,8,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1908)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(8,19)


!L ITEM 9: MEAN NPP ON PFTS FOR TRIFFID

      If (sf(9,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,9,19,im_index)),                        &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 9 = npp_dr_out)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = npp_dr_out(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(9,19,im_index)+(pslevel_out-1)  &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,9,                                         &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1909)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(9,19)


!L ITEM 10: MEAN WOOD RESPIRATION ON PFTS FOR TRIFFID

      If (sf(10,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,10,19,im_index)),                       &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 10 = resp_w_dr_out)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = resp_w_dr_out(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(10,19,im_index)+(pslevel_out-1) &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,10,                                        &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1910)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(10,19)


!L ITEM 11: MEAN SOIL RESPIRATION FOR TRIFFID

      IF (SF(11,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,5)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(11,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,11,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1911)"
           goto 9999
        End if

      END IF     !   sf(11,19)

!L ITEM 12: DISTURBED FRACTION OF VEGETATION

      IF (SF(12,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = frac_disturb(l)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(12,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,12,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1912)"
           goto 9999
        End if

      END IF     !   sf(12,19)

!L ITEM 13: SURFACE TYPE FRACTIONS AFTER TRIFFID

      If (sf(13,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTYPE,LEN_STLIST,                          &
     &       STLIST(1,STINDEX(1,13,19,im_index)),                       &
     &       PLLTYPE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 13 = frac)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NTYPE
          IF (PLLTYPE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = frac(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(13,19,im_index)+(pslevel_out-1) &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,13,                                        &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1913)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(13,19)


!L ITEM 14: LEAF AREA INDEX ON PLANT FUNCTIONAL TYPES AFTER TRIFFID

      If (sf(14,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,14,19,im_index)),                       &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 14 = lai)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = lai(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(14,19,im_index)+(pslevel_out-1) &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,14,                                        &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1914)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(14,19)


!L ITEM 15: CANOPY HEIGHT ON PLANT FUNCTIONAL TYPES AFTER TRIFFID

      If (sf(15,19)) Then

! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NPFT,LEN_STLIST,                           &
     &       STLIST(1,STINDEX(1,15,19,im_index)),                       &
     &       PLLPFT,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "dagveg  : error in set_pseudo_list(item 15 = ht)"
            goto 9999
        END IF
        PSLEVEL_OUT=0

        DO PSLEVEL=1,NPFT
          IF (PLLPFT(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            Do j = 1, rows
              Do i = 1, row_length
                interp_data(i,j) = rmdi
              End Do
            End Do

            Do l = 1, land_pts
              j=(land_index(l)-1)/row_length + 1
              i=land_index(l) - (j-1)*row_length
              interp_data(i,j) = ht(l,pslevel_out)
            End Do

! DEPENDS ON: copydiag
            Call copydiag (STASHwork(si(15,19,im_index)+(pslevel_out-1) &
     &           *row_length*rows),interp_data,                         &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,19,15,                                        &
!#include <argppx/argppx.h>
     &           icode,cmessage)

            If (icode  >   0) then
               cmessage="Error in copydiag( item 1915)"
               goto 9999
            End if
          ENDIF
        ENDDO

      End if     !   sf(15,19)

!L ITEM 16: SOIL CARBON CONTENT AFTER TRIFFID

      IF (SF(16,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,1) + cs(l,2) + cs(l,3) + cs(l,4)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(16,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,16,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1916)"
           goto 9999
        End if

      END IF     !   sf(16,19)


!L ITEM 17-20: MEAN SOIL RESPIRATION FOR TRIFFID, INDIVID. POOLS
! 17: DPM
      IF (SF(17,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,1)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(17,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,17,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1917)"
           goto 9999
        End if
      END IF

! 18: RPM
      IF (SF(18,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,2)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(18,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,18,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1918)"
           goto 9999
        End if
      END IF

! 19: BIO
      IF (SF(19,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,3)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(19,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,19,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1919)"
           goto 9999
        End if
      END IF

! 20: HUM
      IF (SF(20,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = resp_s_dr_out(l,4)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(20,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,20,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1920)"
           goto 9999
        End if
      END IF

!L ITEM 21-24: SOIL CARBON CONTENT AFTER TRIFFID, INDIVID. POOLS
! 21: DPM
      IF (SF(21,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,1)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(21,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,21,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1921)"
           goto 9999
        End if
      END IF

! 22: RPM
      IF (SF(22,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,2)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(22,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,22,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1922)"
           goto 9999
        End if
      END IF

! 23: BIO
      IF (SF(23,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,3)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(23,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,23,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1923)"
           goto 9999
        End if
      END IF

! 24: HUM
      IF (SF(24,19)) THEN
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = rmdi
          End Do
        End Do

        Do l = 1, land_pts
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          interp_data(i,j) = cs(l,4)
        End Do

! DEPENDS ON: copydiag
        Call copydiag (STASHwork(si(24,19,im_index)),interp_data,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,19,24,                                            &
!#include <argppx/argppx.h>
     &       icode,cmessage)

        If (icode  >   0) then
           cmessage="Error in copydiag( item 1924)"
           goto 9999
        End if
      END IF

 9999 continue
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_veg

#endif
