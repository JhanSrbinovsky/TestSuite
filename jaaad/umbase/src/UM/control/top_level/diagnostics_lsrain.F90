#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      Subroutine diagnostics_lsrain(                                    &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      lspice_dim1,lspice_dim2,lspice_dim3        &
     &,                      timestep                                   &
     &,                      at_extremity                               &
     &,                      L_DUST                                     &
     &,                      p_layer_centres                            &
     &,                      T, q, qcl, qcf, cf, cfl, cff               &
     &,                      T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n &
     &,                      ls_rain, ls_snow                           &
     &,                      ls_rain3d,ls_snow3d,rainfrac3d             &
     &,                      RNOUT_TRACER,LSCAV_DUST_ALL,LSCAV_TR       &
     &,                      lscav_nh3                                  &
     &,                      rnout_soot, lscav_soot                     &
     &,                      rnout_bmass, lscav_bmass                   &
     &,                      rnout_ocff, lscav_ocff                     &
     &,                      l_mixing_ratio                             &
     &,                      PSDEP,PSAUT,PSACW,PSACR                    &
     &,                      PSACI,PSMLT,PSMLTEVP                       &
     &,                      PRAUT,PRACW,PREVP                          &
     &,                      PGAUT,PGACW,PGACS,PGMLT                    &
     &,                      PIFRW,PIPRM,PIDEP,PIACW                    &
     &,                      PIACR,PIMLT,PIMLTEVP                       &
     &,                      PIFALL,PSFALL,PRFALL,PGFALL                &
     &,                      PLSET,PLEVPSET                             &
     &     ,                                                            &
#include "argsts.h"
     & STASHwork                                                        &
     &  )

Use ac_diagnostics_mod, Only :  &
    LSRR, LSSR, TINC_PPN

! Purpose:
!          Calculates diagnostics generated from large scale
!          precipitation (UM section 4).
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the large scale
! precipitation routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing, except where indicated.
! NOTE: Although moisture field diagnostics are available from this
! section (and correspond to equivalent variables calculated at UM4.5
! and earlier), the most appropriate place to extract moisture
! diagnostics is after section 9, which is the end of the moist
! processes calculations during the timestep.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (all section 4 )
! ------------------------------------------------------------
! 181 temperature increment across ls precip routines (model levels)
! 182 humidity    increment across ls precip routines (model levels)
! 183 qcl         increment across ls precip routines (model levels)
! 184 qcf         increment across ls precip routines (model levels)
! 192 cf          increment across ls precip routines (model levels)
! 193 cfl         increment across ls precip routines (model levels)
! 194 cff         increment across ls precip routines (model levels)
! 201 large scale rain amount (kg/m2 per timestep)    (surface)
! 202 large scale snow amount (kg/m2 per timestep)    (surface)
! 203 large scale rainfall rate (kg/m2/s)             (surface)
! 204 large scale snowfall rate (kg/m2/s)             (surface)
!   4 temperature           after ls precip           (model levels)
! 205 cloud water (qcl)     after ls precip           (model levels)
! 206 cloud ice (qcf)       after ls precip           (model levels)
! 207 relative humidity     (percent)                 (model levels)
!  10 specific humidity (q) after ls precip           (model levels)
! 220 Large scale rainout of soot (kg/m2/s)           (surface)
! 221 Large scale washout of soot (kg/m2/s)           (surface)
! 223 snowfall rate                                   (model levels)
! 224 supercooled liquid water content                (model levels)
! 225 supercooled rainfall rate                       (model levels)
! 227 rain fraction                                   (model levels)
! 228 Large scale rainout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 229 Large scale washout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 237 Large scale rainout of biomass smoke (kg/m2/s)  (surface)
! 238 Large scale washout of biomass smoke (kg/m2/s)  (surface)
! Microphysical process rate diagnostics all on wet_model_levels
! 240 Homogeneous nucleation rate (kg/kg/s)
! 241 Heterogeneous nucleation rate (kg/kg/s)
! 243 Ice deposition rate (kg/kg/s)
! 245 Snow deposition rate (kg/kg/s)
! 247 Ice collection rate of cloud liquid water (riming) (kg/kg/s)
! 248 Snow collection rate of cloud liquid water (riming) (kg/kg/s)
! 249 Ice collection rate of rain (capture) (kg/kg/s)
! 250 Snow collection rate of rain (capture) (kg/kg/s)
! 251 Evaporation rate of melting ice (kg/kg/s)
! 252 Evaporation rate of melting snow (kg/kg/s)
! 253 Melting rate for ice (kg/kg/s)
! 254 Melting rate for snow (kg/kg/s)
! 255 Snow aggregate autoconversion rate (kg/kg/s)
! 256 Snow collection rate of ice (capture) (kg/kg/s)
! 257 Rain autoconversion rate (kg/kg/s)
! 258 Rain collection rate of cloud liquid water (accretion) (kg/kg/s)
! 259 Evaporation rate of rain (kg/kg/s)
! 260 Graupel autoconversion rate  (kg/kg/s)
! 261 Graupel collection rate of cloud water (accretion) (kg/kg/s)
! 262 Graupel collection rate of snow (capture) (kg/kg/s)
! 263 Melting rate for graupel (kg/kg/s)
! 265 Ice crystal sedimentation rate (kg/kg/s)
! 266 Snow sedimentation rate (kg/kg/s)
! 267 Rain sedimentation rate (kg/kg/s)
! 268 Graupel sedimentation rate (kg/kg/s)
! 269 Droplet settling rate of liquid (kg/kg/s)
! 270 Evaporation rate for settling droplets (kg/kg/s)
!
! History:
! Version   Date     Comment
! ----     -------     -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! 05/04/00  5.1    Provide explicit model increments as STASH output
!                  diagnostics. R Rawlins
! 5.2  21/11/00 Add 3D precipitation field diagnostics. Damian Wilson
! 5.3  04/10/01 Add diagnostics for wet scavenged Sulphur Cycle
!               tracers (SO2 and SO4_DIS)         M Woodage, A Jones
! 5.3  14/08/01 Add (4,207) RH on model levels. R Rawlins
! 5.4  22/07/02 Add PC2 diagnostics.  Damian Wilson
! 5.4  10/07/02 Add (4,220) and (4,221) Rainout and washout of soot.
!                                                           P. Davison
! 5.4  05/09/02  Add diagnostic (4,215) for wet scavenged NH3 (for
!                      Sulphur Cycle).                       M Woodage
! 5.5  05/02/03 Add (4,237) and (4,238) Rainout and washout of biomass
!               smoke aerosol.                          P Davison
!
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
! 6.2  15/02/06 Fix to dust deposition diagnostic.  Stephanie Woodward
!
! 6.2  05/01/06 Add RH wrt water diagnostic (4,208). Richard Forbes.
! 6.2  31/01/06 Add microphysics process rate diags. Richard Forbes
! 6.4  14/08/06 Use mixing ratio version of qsat for 4207. D. Wilson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None
#include "c_dust_ndiv.h"


! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_DUST                                                          &
                        !switch for mineral dust
     &, l_mixing_ratio   ! Use mixing ratio formulation

! Parameters
!
      INTEGER, INTENT(IN) ::                                            &
     &  row_length                                                      &
                            ! number of points on a row
     &, rows                                                            &
                            ! number of rows in a theta field
     &, model_levels                                                    &
                            ! number of model levels
     &, wet_model_levels                                                &
                            ! number of model levels where moisture
     &, lspice_dim1                                                     &
                            ! Dimensions for 3D diagnostic arrays.
     &, lspice_dim2                                                     &
                            !  These are set to 1 in order to save
     &, lspice_dim3         !  memory if the diagnostics are not used.

      REAL, INTENT(IN) ::                                               &
     &  timestep


! Primary Arrays used in all models
      REAL, INTENT(IN) ::                                               &
     &  p_layer_centres(row_length, rows, 0:model_levels)               &
     &, T(row_length, rows, model_levels)                               &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
     &, cf(row_length, rows, wet_model_levels)                          &
     &, cfl(row_length, rows, wet_model_levels)                         &
     &, cff(row_length, rows, wet_model_levels)                         &
! Time level n values for increment diagnostics
     &, T_n  (row_length, rows,     model_levels)                       &
     &, q_n  (row_length, rows, wet_model_levels)                       &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, cf_n  (row_length, rows, wet_model_levels)                      &
     &, cfl_n  (row_length, rows, wet_model_levels)                     &
     &, cff_n  (row_length, rows, wet_model_levels)                     &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, ls_snow3d(lspice_dim1, lspice_dim2, lspice_dim3)                &
     &, rainfrac3d(lspice_dim1, lspice_dim2, lspice_dim3)
! Microphysics Process Rate diagnostic arrays
      REAL, Intent(In) ::                                               &
     &  PSDEP(row_length,rows,wet_model_levels)                         &
     &, PSAUT(row_length,rows,wet_model_levels)                         &
     &, PSACW(row_length,rows,wet_model_levels)                         &
     &, PSACR(row_length,rows,wet_model_levels)                         &
     &, PSACI(row_length,rows,wet_model_levels)                         &
     &, PSMLT(row_length,rows,wet_model_levels)                         &
     &, PSMLTEVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(row_length,rows,wet_model_levels)                         &
     &, PRACW(row_length,rows,wet_model_levels)                         &
     &, PREVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(row_length,rows,wet_model_levels)                         &
     &, PGACW(row_length,rows,wet_model_levels)                         &
     &, PGACS(row_length,rows,wet_model_levels)                         &
     &, PGMLT(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(row_length,rows,wet_model_levels)                         &
     &, PIPRM(row_length,rows,wet_model_levels)                         &
     &, PIDEP(row_length,rows,wet_model_levels)                         &
     &, PIACW(row_length,rows,wet_model_levels)                         &
     &, PIACR(row_length,rows,wet_model_levels)                         &
     &, PIMLT(row_length,rows,wet_model_levels)                         &
     &, PIMLTEVP(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PIFALL(row_length,rows,wet_model_levels)                        &
     &, PSFALL(row_length,rows,wet_model_levels)                        &
     &, PRFALL(row_length,rows,wet_model_levels)                        &
     &, PGFALL(row_length,rows,wet_model_levels)
      REAL, Intent(InOut) ::                                            &
     &  PLSET(row_length,rows,wet_model_levels)                         &
     &, PLEVPSET(row_length,rows,wet_model_levels)


! Used as input and workspace
      REAL, INTENT(INOUT) ::                                            &
     &  ls_rain3d(lspice_dim1, lspice_dim2, lspice_dim3)                &
     &, rnout_tracer(lspice_dim1, lspice_dim2)                          &
     &, LSCAV_DUST_ALL(ROW_LENGTH,ROWS,NDIV)                            &
                                             !scavenged mineral dust

     &, lscav_tr(lspice_dim1, lspice_dim2)                              &
     &, lscav_nh3(lspice_dim1, lspice_dim2)                             &
     &, rnout_soot(lspice_dim1, lspice_dim2)                            &
     &, lscav_soot(lspice_dim1, lspice_dim2)                            &
     &, rnout_bmass(lspice_dim1, lspice_dim2)                           &
     &, lscav_bmass(lspice_dim1, lspice_dim2)                           &
     &, rnout_ocff(lspice_dim1, lspice_dim2)                            &
     &, lscav_ocff(lspice_dim1, lspice_dim2)

#include "csubmodl.h"
#include "typsts.h"
#include "c_0_dg_c.h"

! Diagnostic variables
      REAL, INTENT(INOUT) ::                                            &
     & STASHwork(*)     ! STASH workspace for section 4 (LS precip)

! Local variables
      Integer                                                           &
     & i, j, k, ji                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Integer sect,item    ! STASH section, item no.s
      Parameter (sect = 4) !  for microphysics - large scale rain

      Real work_3d(row_length, rows,model_levels) ! work space

      Character*80  cmessage

      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_lsrain')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     & copydiag, copydiag_3d                                            &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

!L Copy diagnostic information to STASHwork for STASH processing

! increment diagnostics= modified - previous

      item = 181  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = T(i,j,k) - T_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! And set dry level increments to zero explicitly
         Do k = wet_model_levels+1,model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = 0.0
             End do ! i
           End do   ! j
         End do     ! k

        If (.not.allocated(TINC_PPN)) Then
          Allocate ( TINC_PPN(row_length*rows,model_levels) )
        End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

         Do k = 1,model_levels
           Do j = 1,rows
             Do i = 1,row_length
               ji = (j-1)*row_length+i
               TINC_PPN(ji,k) = work_3d(i,j,k) 
             End do ! i
           End do   ! j
         End do     ! k

      End if  !  sf(item,sect)

      item = 182  ! humidity increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k) - q_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 183  ! qcl increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = qcl(i,j,k) - qcl_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 183)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 184  ! qcf increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = qcf(i,j,k) - qcf_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 184)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 192  ! cf increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cf(i,j,k) - cf_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 193  ! cfl increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cfl(i,j,k) - cfl_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 194  ! cff increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = cff(i,j,k) - cff_n(i,j,k)
             End do ! i
           End do   ! j
         End do     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 194)"//cmessage
         End if

      End if  !  sf(item,sect)

! Item 201 Large scale rain

      If (icode <= 0 .and. sf(201,4)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,4,im_index)),ls_rain,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_rain)"
         Endif
! Code to convert rate to amount for a given timestep

         Do i=1,row_length*rows
            STASHwork(si(201,4,im_index)+i-1)=                          &
     &           STASHwork(si(201,4,im_index)+i-1)* timestep
         End do

      End if


! Item 202 Large scale snow

      If (icode <= 0 .and. sf(202,4)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(202,4,im_index)),ls_snow,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,202,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_snow)"
         End if

         Do i=1,row_length*rows
            STASHwork(si(202,4,im_index)+i-1)=                          &
     &           STASHwork(si(202,4,im_index)+i-1)* timestep
        End do

      End if


! Item 203 Large scale rain

      If (icode <= 0 .and. sf(203,4)) Then

      If (.not.allocated(LSRR)) Then
        Allocate ( LSRR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(203,4,im_index)),ls_rain,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,203,                                           &
     &       icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_rain)"
         End if

         Do j = 1,rows
           Do i = 1,row_length
             ji = (j-1)*row_length+i
             LSRR(ji) = ls_rain(i,j)
           End Do
         End Do

      End if


! Item 204 Large scale snow

      If (icode <= 0 .and. sf(204,4)) then

      If (.not.allocated(LSSR)) Then
        Allocate ( LSSR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,4,im_index)),ls_snow,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag(ls_snow)"
         End if

         Do j = 1,rows
           Do i = 1,row_length
             ji = (j-1)*row_length+i
             LSSR(ji) = ls_snow(i,j)
           End Do
         End Do

      End if



! Items 231-236 mineral dust scavenged by LS precip

      IF (L_DUST) THEN

        DO K = 1,NDIV

          IF (ICODE <= 0 .AND. SF(230+K,4)) THEN

!       Convert "per timestep" diagnostic to "per second":
            DO J=1,LSPICE_DIM2
              DO I=1,LSPICE_DIM1
                LSCAV_DUST_ALL(I, J, K)=LSCAV_DUST_ALL(I, J, K)/TIMESTEP
              ENDDO
            ENDDO

! DEPENDS ON: copydiag
            CALL COPYDIAG(STASHWORK(SI(230+K,4,IM_INDEX)),              &
     &         LSCAV_DUST_ALL(1:ROW_LENGTH,1:ROWS,K),                   &
     &         ROW_LENGTH,ROWS,0,0,0,0, AT_EXTREMITY,                   &
     &         ATMOS_IM,4,230+K,                                        &
     &         ICODE,CMESSAGE)

            IF (ICODE  >   0) THEN
              CMESSAGE=": ERROR IN COPYDIAG(LSCAV_DUST_ALL)"
            ENDIF

          ENDIF

        ENDDO !NDIV

      ENDIF !L_DUST


! Item 215 NH3 scavenged by LS precip

      If (icode <= 0 .and. sf(215,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_nh3(i, j)=lscav_nh3(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(215,4,im_index)),lscav_nh3,          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,215,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage="microphysics_ctl  : error in copydiag(lscav_nh3)"
        End if

      End if
!
! Item 216 SO2 scavenged by LS precip

      If (icode <= 0 .and. sf(216,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_tr(i, j)=lscav_tr(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(216,4,im_index)),lscav_tr,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,216,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_tr)"
        End if

      End if

! Item 219 Dissolved SO4 aerosol scavenged by LS precip

      If (icode <= 0 .and. sf(219,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_tracer(i, j)=rnout_tracer(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(219,4,im_index)),rnout_tracer,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,219,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_tracer)"
        End if

      End if
!
! Item 220 Soot scavenged by LS rainout

      If (icode <= 0 .and. sf(220,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_soot(i, j)=rnout_soot(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(220,4,im_index)),rnout_soot,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,220,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_soot)"
        End if

      End if

! Item 221 Soot scavenged by LS washout

      If (icode <= 0 .and. sf(221,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_soot(i, j)=lscav_soot(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(221,4,im_index)),lscav_soot,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,221,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_soot)"
        End if

      End if
!
! Item 228 Fossil-fuel organic carbon scavenged by LS rainout

      If (icode <= 0 .and. sf(228,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_ocff(i, j)=rnout_ocff(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(228,4,im_index)),rnout_ocff,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,228,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_ocff)"
        End if

      End if

! Item 229 Fossil-fuel organic carbon scavenged by LS washout

      If (icode <= 0 .and. sf(229,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_ocff(i, j)=lscav_ocff(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(229,4,im_index)),lscav_ocff,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,229,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_ocff)"
        End if

      End if
!
!L Copy T to STASHwork

      If (icode <= 0 .and. sf(004,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,4,im_index)),T,              &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,004,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,004,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(T)"
         End if

      End if



!L Copy Cloud water to STASHwork

      If (icode <= 0 .and. sf(205,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(205,4,im_index)),qcl,            &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,205,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,205,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(cloud water)"
         End if

      End if

      If (icode <= 0 .and. sf(206,4)) Then

!L Copy Cloud ice to STASHwork

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(206,4,im_index)),qcf,            &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,206,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,206,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(cloud ice)"
         End if

      End if

! -----------------------------
! 207 relative humidity wrt ice (T<0degC) and water (T>0degC) (mdl levs)

      item = 207
      If (icode <= 0 .and. sf(item,sect)) Then
! DEPENDS ON: qsat_mix
         Call qsat_mix(work_3d,T,p_layer_centres(1,1,1),                &
     &             row_length*rows*wet_model_levels,l_mixing_ratio)
         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
               If (work_3d(i,j,k)  <   0.0) Then
                 work_3d(i,j,k) = 0.0
               End if

             End do  ! i
           End do    ! j
         End do      ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 207)"//cmessage
         End if

      End if   ! item 207

! -------------------------------
! 208 relative humidity wrt water (on model levels)

      item = 208
      If (icode  <=  0 .and. sf(item,sect)) Then
         ! q saturation is put in work_3d array
! DEPENDS ON: qsat_wat
         Call qsat_wat(work_3d,T,p_layer_centres(1,1,1),                &
     &             row_length*rows*wet_model_levels )
         Do k = 1,wet_model_levels
           Do j = 1,rows
             Do i = 1,row_length
               work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.
               !  Supersaturation wrt water is limited to =< 100%
               If (work_3d(i,j,k) > 100.0) Then
                 work_3d(i,j,k) = 100.0
               End if
               !  Negative humidity also removed from the diagnostic
               If (work_3d(i,j,k) < 0.0) Then
                 work_3d(i,j,k) = 0.0
               End if

             End do  ! i
           End do    ! j
         End do      ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d (item 208, rhw)"//cmessage
         End if

      End if   ! item 208
!L Copy q  to STASHwork

      If (icode <= 0 .and. sf(010,4)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(010,4,im_index)),q,              &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,010,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,010,                                           &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(q)"
         End if

      End if


! Copy ls_rain3d  to STASHwork
!
      If (icode <= 0 .and. sf(222,4)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(222,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,222,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,222,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(ls_rain3d)"
         End if
!
      End if
!
!
! Copy ls_snoww3d  to STASHwork
!
      If (icode <= 0 .and. sf(223,4)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,4,im_index)),ls_snow3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,223,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,223,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(ls_snow3d)"
         End if
!
      End if
!
!
! Need to produce diagnostic 225 before 224 in order to save memory.
!
      If(icode <= 0 .and. sf(225,4)) Then
!
!  Supercooled 3D rain content. It is equal to
!  the 3D rainrate at T < 0 and equal to 0 at T > 0
!  Alter the array LS_RAIN3D directly
!
          Do k=1,wet_model_levels
            Do j=1,rows
              Do i=1,row_length
                If (T(i,j,k)  >=  zerodegc) THEN
! Warm temperatures
                  ls_rain3d(i,j,k)=0.0
                End if
              End do
            End do
          End do
!
! Copy supercooled rain (now in ls_rain3d)  to STASHwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,225,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,225,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(supercooled rain)"
         End if
!
      End if
!
      IF (icode <= 0 .and. SF(224,4)) THEN
!
!  Supercooled liquid water content at TIME LEVEL N. It is equal to
!  the liquid water content at T < 0 and equal to 0 at T > 0
!  Use LS_RAIN3D as the array in order to save memory
!
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              If (T(i,j,k)  <   zerodegc) Then
! Supercooled temperatures
! Use time level n fields in this diagnostic
                ls_rain3d(i,j,k)=qcl_n(i,j,k)
              Else
! Warm temperatures
                ls_rain3d(i,j,k)=0.0
              End if
            End do
          End do
        End do
!
! Copy supercooled liquid (now in ls_rain3d)  to STASHwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,4,im_index)),ls_rain3d,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,224,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,224,                                           &
     &        icode,cmessage)
!
        If (icode >  0) Then
           cmessage=":error in copydiag_3d(supercooled liq)"
        End if
!
      End if
!
!
      If (icode <= 0 .and. sf(227,4)) Then
!
! Copy rain fraction to stashwork
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(227,4,im_index)),rainfrac3d,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,227,4,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,4,227,                                           &
     &        icode,cmessage)
!
         If (icode >  0) Then
            cmessage=":error in copydiag_3d(rain fraction)"
         End if
!
      End if
!
! Item 237 Biomass scavenged by LS rainout

      If (icode <= 0 .and. sf(237,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            rnout_bmass(i, j)=rnout_bmass(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(237,4,im_index)),rnout_bmass,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,237,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(rnout_bmass)"
        End if

      End if

! Item 238 Biomass scavenged by LS washout

      If (icode <= 0 .and. sf(238,4)) then

!       Convert "per timestep" diagnostic to "per second":

        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_bmass(i, j)=lscav_bmass(i, j)/timestep
          End Do
        End Do

! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(238,4,im_index)),lscav_bmass,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,4,238,                                           &
     &        icode,cmessage)

        If (icode  >   0) Then
          cmessage=": error in copydiag(lscav_bmass)"
        End if

      End if
!

        !---------------------------------------------------------------
        ! Homogeneous nucleation
        !---------------------------------------------------------------
        item = 240

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIFRW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Heterogeneous nucleation
        !---------------------------------------------------------------
        item = 241

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIPRM,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Deposition of ice
        !---------------------------------------------------------------
        item = 243

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIDEP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Deposition of snow aggregates
        !---------------------------------------------------------------
        item = 245

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSDEP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If


        !---------------------------------------------------------------
        ! Ice collection of cloud liquid water (riming)
        !---------------------------------------------------------------
        item = 247

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of cloud liquid water (riming)
        !---------------------------------------------------------------
        item = 248

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Ice collection of rain (capture)
        !---------------------------------------------------------------
        item = 249

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIACR,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of rain (capture)
        !---------------------------------------------------------------
        item = 250

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACR,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If
        !---------------------------------------------------------------
        ! Evaporation of melting ice
        !---------------------------------------------------------------
        item = 251

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIMLTEVP,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporation of melting snow
        !---------------------------------------------------------------
        item = 252

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSMLTEVP,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting ice
        !---------------------------------------------------------------
        item = 253

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting snow
        !---------------------------------------------------------------
        item = 254

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow aggregate autoconversion
        !---------------------------------------------------------------
        item = 255

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow collection of ice (capture)
        !---------------------------------------------------------------
        item = 256

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSACI,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain autoconversion
        !---------------------------------------------------------------
        item = 257

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
        item = 258

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporation of rain
        !---------------------------------------------------------------
        item = 259

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PREVP,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel autoconversion
        !---------------------------------------------------------------
        item = 260

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGAUT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel collection of cloud liquid water (accretion)
        !---------------------------------------------------------------
        item = 261

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGACW,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel collection of snow (capture)
        !---------------------------------------------------------------
        item =262

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGACS,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Melting graupel
        !---------------------------------------------------------------
        item = 263

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGMLT,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Ice crystal sedimentation
        !---------------------------------------------------------------
        item = 265

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PIFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Snow sedimentation
        !---------------------------------------------------------------
        item = 266

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PSFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Rain sedimentation
        !---------------------------------------------------------------
        item = 267

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PRFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Graupel sedimentation
        !---------------------------------------------------------------
        item = 268

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PGFALL,                                                   &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Droplet settling of liquid
        !---------------------------------------------------------------
        item = 269

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PLSET,                                                    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

        !---------------------------------------------------------------
        ! Evaporated settled droplets
        !---------------------------------------------------------------
        item = 270

        If (icode <= 0 .and. SF(item,sect)) Then
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (                                          &
     &        STASHwork(si(item,sect,im_index)),                        &
     &        PLEVPSET,                                                 &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),                  &
     &        len_stlist, stash_levels,num_stash_levels+1,              &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

            If (icode  >   0) Then
              Write(Cmessage,*)                                         &
     &          'ERROR in copydiag_3d for diagnostic ',                 &
     &          'section 4, item ',item
            End If
        End If

! Single point exception handling
      If (icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      End if

      Return
      END SUBROUTINE diagnostics_lsrain

#endif
