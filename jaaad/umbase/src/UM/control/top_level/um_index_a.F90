#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: UM_INDEX-----------------------------------------------
!LL
!LL  Purpose: Calculate addresses of component arrays within a
!LL           series of super arrays, made up of combinations of arrays
!LL           that require dynamic allocation. Lengths of the super
!LL           arrays are calculated and passed into U_MODEL for
!LL           dynamic allocation, reducing the no. of arguments needed
!LL           to be passed between top-level routines.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Introduced as new DECK to allow dynamic allocation
!LL                  of main data arrays in U_MODEL.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  re-dimension PP_XREF and add INDEX_PPXREF.
!LL  3.5   29/03/95  MPP code : Land point fields are allocated P_FIELD
!LL                  amount of space in D1       P.Burton
!LL  3.5   Apr. 95   Sub-Models project.
!LL                  STASH super array modified in accordance with
!LL                  internal model separation scheme for diagnostics.
!LL                  *CALL CSUBMODL introduced.
!LL                  S.J.Swarbrick
!LL  4.0   06/09/95  Added atmos/ocean stash superarrays. K Rogers
!LL  4.1  15/03/96  Introduce Wave sub-model.  RTHBarnes.
!LL  4.1   04/12/95  Increased  A_IXPTR to accomodate 2 extra prognostic
!LL                  arrays J.Smith
!LL  4.1   26/04/96  Increased A_IXPTR to allow for 12 Sulphur Cycle
!LL                  prognostics and ancillaries        MJWoodage
!LL  4.2  11/10/96   Enable atmos-ocean coupling for MPP.
!LL                  (1): Coupled fields. Change to 'global' sizes
!LL                  instead of local. R.Rawlins
!LL  4.2   11/10/96  Enable atmos-ocean coupling for MPP.
!LL                  (2): Swap D1 memory. Add copies of D1 for atmos and
!LL                  ocean. R.Rawlins
!LL  4.3   26/03/97  Added HadCM2 sulphate loading pattern.  Will Ingram
!LL  4.4   01/07/97  Added padding to place the D1 array on a
!LL                  cache line boundary.
!LL                    Author: Bob Carruthers, Cray Research.
!LL  4.4   05/08/97  Allowed prognostic variable CCA to be 3D. JMG
!LL  4.4   11/10/97  Rename AO_D1_MEMORY to L_AO_D1_MEMORY. D Robinson
!LL  4.4   13/10/97  Initialise LEN_A/O/W_SPSTS. D. Robinson.
!LL  4.5   29/07/98  Use U_FIELD_INTFA to set up A_IXINF(21).
!LL                  Compute new A_IXINF(22-25). D. Robinson.
!LL  4.5   04/03/98  Increase IXPTR to allow for 1 new NH3 var in
!LL                  S Cycle and 3 new soot vars.     M Woodage
!LL  4.5   08/05/98  Increase A_IXPTR by 16 to increase maximum number
!LL                  of multi-level user ancillaries to 20
!LL                  Author D.M. Goddard
!LL  4.5   13/05/98  Added RHcrit variable to D1 pointers. S. Cusack
!LL  4.5   15/07/98  Added 3D CO2 to D1 pointers. C.D.Jones
!LL  4.5   17/08/98  Remove pointers for JSOIL_FLDS and JVEG_FLDS.
!LL                  D. Robinson.
!    4.5   30/03/98  Reserve space for land mask in ocean
!                    stash array. R. Hill
!LL  5.0   28/04/99  Remove references to FLOORFLDSA    P.Burton
!LL  5.0   20/05/99 New pointers for new dynamics data arrays added
!                   Redundant pointers removed.
!                   D.M. Goddard
!LL  5.0   21/05/99 Change variable names P.Burton
!    5.0   21/05/99  (1) STASH sizes now held separately.
!                    (2) Atmosphere STASH array changes. R.Rawlins
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!  5.1  28/04/00   Extra pressure level above top theta level
!                  p_top_theta_level removed                A.Malcolm
!    5.1   07/02/00 Extend super array ao_ixcpl() in order to include
!                   the extra arrays xva() & yva()
!    5.2   10/11/00 Extend a_ixinf to include lbc_eta_theta (26) and
!                   lbc_eta_rho (27). Rename MAX_INTF_P_LEVELS to
!                   MAX_INTF_MODEL_LEVELS. D.Robinson
!
!    5.3   19/06/01 Add tropopause-based ozone array. Dave Tan
!    5.3   01/10/01 Extend A_IXCON by 2 arrays for checkerboard
!                   radiation. S Cusack
!    5.3   15/10/01 Include more ocean boundary data. M J Bell
!    5.3   18/07/01 1) Remove unneeded allocation of 3 local array
!                   images of D1.
!                   2) Tidy: remove superfluous #if defined (mpp) and
!                   add test for writing out info messages.
!                   3) Remove innocuous array out of bounds message
!                   for a_spptr. R Rawlins
!    5.4  21/03/02 Remove comment on same line as #include
!                                               S. Carroll
!    5.5   17/02/03 Set pointers for Wave interface code. D.Holmes-Bell
!    5.5   17/02/03 Upgrade Wave model from 4.1 to 5.5, & incorporate
!                   parallel mods. D.Holmes-Bell
!    5.5  05/08/00  Modification for parallelisation of WAM.
!                   Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    5.5  05/02/03  Add biomass aerosol arrays           P Davison
!    5.5  28/02/03  Set up Super array of ATMOS lat. and long.
!                   values for river routing if ATMOS only model
!                                           . C.Bunton
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!    6.0  30/07/03  Include cloud fraction lbcs. Damian Wilson
!    6.1  01/09/04  End of TIMER(UM_SETUP) moved.  R.Barnes
!    6.1  20/08/03  Add STOCHEM pointers.  C.Johnson
!    6.1  18/08/04  Set up A_IXINF(28-33). D Robinson.
!
!  6.1   04/08/04   Add dynamic arrays for diffusion coefficients
!                                                      Terry Davies
!    6.2  23/11/05  Removed all references to the wavemodel.
!                   T.Edwards
!    6.2  02/03/06  Added pointer for direct surface PAR flux.
!                                        M.G. Sanderson
!  6.2   14/02/06   Modify dynamic arrays for diffusion coefficients
!  6.2   14/02/06   Add dynamic arrays for variable resolution
!                                                      Terry Davies
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
#endif
#if defined(CONTROL) && defined(ATMOS)
!LL  Subroutine: UM_INDEX_A---------------------------------------------
!LL
!LL  Purpose: Calculate addresses and lengths of atmosphere super arrays
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Introduced as new DECK to allow dynamic allocation
!LL                  of main data arrays in U_MODEL.
!LL  3.2  22/11/93  Add 3 more A_IXPTR values for aerosol pointers.
!LL                                                   R.T.H.Barnes.
!LL  3.4  29/09/94  Add 6 more A_IXPTR values for multi-level murk and
!LL                  user ancillaries.                R.T.H.Barnes.
!LL  3.4    18/5/94 Add extra u field to A_SPCON (argcona). J Thomson
!LL  4.0  06/09/95  Added atmos stash superarray. K Rogers
!LL  5.0  09/06/99  New addressing for C-P C-grid derived constants for
!LL                 atmosphere, A_IXCON().  M L Gallani.
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.1  26/04/00  Added LBC variables to A_IXPTR variables  P.Burton
!    5.2   01/09/00 - Remove levels from LBC pointers
!                   - Add orography LBC
!                   - Add extra space for LBC LOOKUPs
!                                                  P.Burton
!LL  5.2  13/09/00  Removed 3 aerosol pointers. P.Selwood
!LL  5.2  09/03/01  Replace eta_levels by zseak and Ck arrays for
!LL                 passing into STASH super-array. R Rawlins
!    5.4  02/09/02  Added sea mask as 12th element of a_ixsts K.Williams
!    5.5  03/02/03  Add 3 moisture variables and lbcs. R.M.Forbes
!    6.2  11/05/05  Correct jch4_stoch pointer. M.G. Sanderson.
!    6.2  01/10/04  Add murk aerosol lbc.   R.M.Forbes
!    6.2  05/01/06   Add true_latitude to A_IXCON. Yongming Tang
!    6.2  14/11/04  Add JTR_UKCA UKCA tracer pointer.   R Barnes
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE UM_INDEX_A(                                            &
#include "argszspa.h"
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      USE CSENARIO_MOD
      IMPLICIT NONE
!
!  Subroutines called
!
!  Local variables
!
      INTEGER ICODE             ! Work - Internal return code
      CHARACTER*80  CMESSAGE    ! Work - Internal error message
!
!  Configuration-dependent sizes for dynamic arrays
!
#include "parvars.h"
#include "typsize.h"
!
! Parameters for constants arrays
#include "cmaxsize.h"
#include "cconsts.h"
!  Ancillary file parameters for ancillary length calculations
#include "conanc.h"
!
!  Super array sizes for dynamic allocation in U_MODEL
!
#include "typszspa.h"
!
!
!  Addresses of arrays in super arrays.
!
#include "spindex.h"
!
!L----------------------------------------------------------------------
!L
      ICODE=0

!L
!L 1.0 Calculate size of atmosphere stash array
!L

! [Note that size of a_ixsts must be the same as o_ixsts and any other
!  submodels, since lower level STASH routines assume the existence of
!  _ixsts(i) even if not used.]

      a_ixsts(1) = 1                              ! zseak_rho
      a_ixsts(2) = a_ixsts(1) + model_levels      ! Ck_rho
      a_ixsts(3) = a_ixsts(2) + model_levels      ! zseak_theta
      a_ixsts(4) = a_ixsts(3) + model_levels + 1  ! Ck_theta
      a_ixsts(5) = a_ixsts(4) + model_levels + 1  ! dummy
      a_ixsts(6) = a_ixsts(5) + 1                 ! dummy
      a_ixsts(7) = a_ixsts(6) + 1                 ! p     pointer in D1
      a_ixsts(8) = a_ixsts(7) + 1                 ! pstar pointer in D1
      a_ixsts(9) = a_ixsts(8) + 1                 ! cos_v_latitude
                                                  !         (no halos)
      a_ixsts(10)= a_ixsts(9) + row_length*n_rows ! cos_theta_latitude
                                                  !         (no halos)
      a_ixsts(11)= a_ixsts(10)+ row_length*rows   ! landsea mask
      a_ixsts(12)= a_ixsts(11)+ row_length*rows   ! land fraction
      a_spsts_len= a_ixsts(12)+ row_length*rows
!     Store A_SPSTS_LEN in TYPSIZE
      LEN_A_SPSTS = A_SPSTS_LEN

!L
!L 1.1   DUMP     super array
!L
!L          super array addresses
      A_IXDUM(1) =1
      A_IXDUM(2) =A_IXDUM(1) + LEN_FIXHD
      A_IXDUM(3) =A_IXDUM(2) + A_LEN_INTHD
      A_IXDUM(4) =A_IXDUM(3) + A_LEN_CFI1+1
      A_IXDUM(5) =A_IXDUM(4) + A_LEN_CFI2+1
      A_IXDUM(6) =A_IXDUM(5) + A_LEN_CFI3+1
      A_IXDUM(7) =A_IXDUM(6) + A_LEN_REALHD
      A_IXDUM(8) =A_IXDUM(7) + A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1
      A_IXDUM(9) =A_IXDUM(8) + A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1
      A_IXDUM(10)=A_IXDUM(9) + A_LEN1_COLDEPC*A_LEN2_COLDEPC+1
      A_IXDUM(11)=A_IXDUM(10)+ A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1
      A_IXDUM(12)=A_IXDUM(11)+ A_LEN_EXTCNST+1
      A_IXDUM(13)=A_IXDUM(12)+ LEN_DUMPHIST+1
      A_IXDUM(14)=A_IXDUM(13)+ LEN1_LOOKUP*A_LEN2_LOOKUP

!L
!L          super array length
      A_SPDUM_LEN  =A_IXDUM(14)+ MPP_LEN1_LOOKUP*A_LEN2_LOOKUP
      A_SPDUM_LEN  =A_SPDUM_LEN -1
!L
!L
!L 1.2   Pointers super array
!L
!L          super array addresses

! Comment at end of each line corresponds to the matching
! pointer in ARGPTRA. eg A_IXPTR(3) in ARTPTRA = jtheta in ARGPTRA

! For each line : A_IXPTR(n+1) = A_IXPTR(n) + n_levs
! where n_levs in the no of levels for pointer n.

      A_IXPTR(1) =1                                  ! ju
      A_IXPTR(2) =A_IXPTR(1) + MODEL_LEVELS          ! jv
      A_IXPTR(3) =A_IXPTR(2) + MODEL_LEVELS          ! jw
      A_IXPTR(4) =A_IXPTR(3) + MODEL_LEVELS+1        ! jrho
      A_IXPTR(5) =A_IXPTR(4) + MODEL_LEVELS          ! jtheta
      A_IXPTR(6) =A_IXPTR(5) + MODEL_LEVELS          ! jq
      A_IXPTR(7) =A_IXPTR(6) + WET_LEVELS            ! jqcl
      A_IXPTR(8) =A_IXPTR(7) + WET_LEVELS            ! jqcf
      A_IXPTR(9) =A_IXPTR(8) + WET_LEVELS            ! jexner_rho_lvls
      A_IXPTR(10)=A_IXPTR(9) + MODEL_LEVELS+1         ! ju_adv
      A_IXPTR(11)=A_IXPTR(10)+ MODEL_LEVELS          ! jv_adv
      A_IXPTR(12)=A_IXPTR(11)+ MODEL_LEVELS          ! jw_adv
      A_IXPTR(13)=A_IXPTR(12)+ MODEL_LEVELS+1        ! jp
      A_IXPTR(14)=A_IXPTR(13)+ MODEL_LEVELS+1         ! jp_theta_levels
      A_IXPTR(15)=A_IXPTR(14)+ MODEL_LEVELS          ! jexner_theta_lvls
      A_IXPTR(16)=A_IXPTR(15)+ MODEL_LEVELS          ! jcca
      A_IXPTR(17)=A_IXPTR(16)+ N_CCA_LEV             ! jcf_area
      A_IXPTR(18)=A_IXPTR(17)+ WET_LEVELS            ! jcf_bulk
      A_IXPTR(19)=A_IXPTR(18)+ WET_LEVELS            ! jcf_liquid
      A_IXPTR(20)=A_IXPTR(19)+ WET_LEVELS            ! jcf_frozen
      A_IXPTR(21)=A_IXPTR(20)+ WET_LEVELS            ! j_deep_soil_temp
      A_IXPTR(22)=A_IXPTR(21)+ ST_LEVELS             ! jsmcl
      A_IXPTR(23)=A_IXPTR(22)+ SM_LEVELS             ! jsthu
      A_IXPTR(24)=A_IXPTR(23)+ SM_LEVELS             ! jsthf
      A_IXPTR(25)=A_IXPTR(24)+ SM_LEVELS             ! jsw_incs
      A_IXPTR(26)=A_IXPTR(25)+ MODEL_LEVELS+2        ! jlw_incs
      A_IXPTR(27)=A_IXPTR(26)+ MODEL_LEVELS+1        ! jozone
      A_IXPTR(28)=A_IXPTR(27)+ OZONE_LEVELS          ! jtracer
      A_IXPTR(29)=A_IXPTR(28)+ TR_LEVELS*(TR_VARS+1) ! jmurk
      A_IXPTR(30)=A_IXPTR(29)+ MODEL_LEVELS          ! jmurk_source
      A_IXPTR(31)=A_IXPTR(30)+ MODEL_LEVELS          ! jso2
      A_IXPTR(32)=A_IXPTR(31)+ MODEL_LEVELS          ! jdms
      A_IXPTR(33)=A_IXPTR(32)+ MODEL_LEVELS          ! jso4_aitken
      A_IXPTR(34)=A_IXPTR(33)+ MODEL_LEVELS          ! jso4_accu
      A_IXPTR(35)=A_IXPTR(34)+ MODEL_LEVELS          ! jso4_diss
      A_IXPTR(36)=A_IXPTR(35)+ MODEL_LEVELS          ! jh2o2
      A_IXPTR(37)=A_IXPTR(36)+ MODEL_LEVELS          ! jnh3
      A_IXPTR(38)=A_IXPTR(37)+ MODEL_LEVELS          ! jsoot_new
      A_IXPTR(39)=A_IXPTR(38)+ MODEL_LEVELS          ! jsoot_agd
      A_IXPTR(40)=A_IXPTR(39)+ MODEL_LEVELS          ! jsoot_cld
      A_IXPTR(41)=A_IXPTR(40)+ MODEL_LEVELS          ! jso2_natem
      A_IXPTR(42)=A_IXPTR(41)+ MODEL_LEVELS          ! joh
      A_IXPTR(43)=A_IXPTR(42)+ MODEL_LEVELS          ! jho2
      A_IXPTR(44)=A_IXPTR(43)+ MODEL_LEVELS          ! jh2o2_limit
      A_IXPTR(45)=A_IXPTR(44)+ MODEL_LEVELS          ! jo3_chem
      A_IXPTR(46)=A_IXPTR(45)+ MODEL_LEVELS          ! jhadcm2_so4
      A_IXPTR(47)=A_IXPTR(46)+ NSULPAT               ! jco2
      A_IXPTR(48)=A_IXPTR(47)+ MODEL_LEVELS          ! juser_mult1
      A_IXPTR(49)=A_IXPTR(48)+ MODEL_LEVELS          ! juser_mult2
      A_IXPTR(50)=A_IXPTR(49)+ MODEL_LEVELS          ! juser_mult3
      A_IXPTR(51)=A_IXPTR(50)+ MODEL_LEVELS          ! juser_mult4
      A_IXPTR(52)=A_IXPTR(51)+ MODEL_LEVELS          ! juser_mult5
      A_IXPTR(53)=A_IXPTR(52)+ MODEL_LEVELS          ! juser_mult6
      A_IXPTR(54)=A_IXPTR(53)+ MODEL_LEVELS          ! juser_mult7
      A_IXPTR(55)=A_IXPTR(54)+ MODEL_LEVELS          ! juser_mult8
      A_IXPTR(56)=A_IXPTR(55)+ MODEL_LEVELS          ! juser_mult9
      A_IXPTR(57)=A_IXPTR(56)+ MODEL_LEVELS          ! juser_mult10
      A_IXPTR(58)=A_IXPTR(57)+ MODEL_LEVELS          ! juser_mult11
      A_IXPTR(59)=A_IXPTR(58)+ MODEL_LEVELS          ! juser_mult12
      A_IXPTR(60)=A_IXPTR(59)+ MODEL_LEVELS          ! juser_mult13
      A_IXPTR(61)=A_IXPTR(60)+ MODEL_LEVELS          ! juser_mult14
      A_IXPTR(62)=A_IXPTR(61)+ MODEL_LEVELS          ! juser_mult15
      A_IXPTR(63)=A_IXPTR(62)+ MODEL_LEVELS          ! juser_mult16
      A_IXPTR(64)=A_IXPTR(63)+ MODEL_LEVELS          ! juser_mult17
      A_IXPTR(65)=A_IXPTR(64)+ MODEL_LEVELS          ! juser_mult18
      A_IXPTR(66)=A_IXPTR(65)+ MODEL_LEVELS          ! juser_mult19
      A_IXPTR(67)=A_IXPTR(66)+ MODEL_LEVELS          ! juser_mult20
      A_IXPTR(68)=A_IXPTR(67)+ MODEL_LEVELS          ! j_orog_lbc
      A_IXPTR(69)=A_IXPTR(68)+ 1                     ! ju_lbc
      A_IXPTR(70)=A_IXPTR(69)+ 1                     ! ju_lbc_tend
      A_IXPTR(71)=A_IXPTR(70)+ 1                     ! jv_lbc
      A_IXPTR(72)=A_IXPTR(71)+ 1                     ! jv_lbc_tend
      A_IXPTR(73)=A_IXPTR(72)+ 1                     ! jw_lbc
      A_IXPTR(74)=A_IXPTR(73)+ 1                     ! jw_lbc_tend
      A_IXPTR(75)=A_IXPTR(74)+ 1                     ! jrho_lbc
      A_IXPTR(76)=A_IXPTR(75)+ 1                     ! jrho_lbc_tend
      A_IXPTR(77)=A_IXPTR(76)+ 1                     ! jtheta_lbc
      A_IXPTR(78)=A_IXPTR(77)+ 1                     ! jtheta_lbc_tend
      A_IXPTR(79)=A_IXPTR(78)+ 1                     ! jq_lbc
      A_IXPTR(80)=A_IXPTR(79)+ 1                     ! jq_lbc_tend
      A_IXPTR(81)=A_IXPTR(80)+ 1                     ! jqcl_lbc
      A_IXPTR(82)=A_IXPTR(81)+ 1                     ! jqcl_lbc_tend
      A_IXPTR(83)=A_IXPTR(82)+ 1                     ! jqcf_lbc
      A_IXPTR(84)=A_IXPTR(83)+ 1                     ! jqcf_lbc_tend
      A_IXPTR(85)=A_IXPTR(84)+ 1                     ! jexner_lbc
      A_IXPTR(86)=A_IXPTR(85)+ 1                     ! jexner_lbc_tend
      A_IXPTR(87)=A_IXPTR(86)+ 1                     ! ju_adv_lbc
      A_IXPTR(88)=A_IXPTR(87)+ 1                     ! ju_adv_lbc_tend
      A_IXPTR(89)=A_IXPTR(88)+ 1                     ! jv_adv_lbc
      A_IXPTR(90)=A_IXPTR(89)+ 1                     ! jv_adv_lbc_tend
      A_IXPTR(91)=A_IXPTR(90)+ 1                     ! jw_adv_lbc
      A_IXPTR(92)=A_IXPTR(91)+ 1                     ! jw_adv_lbc_tend
      A_IXPTR(93)=A_IXPTR(92)+ 1                     ! jtracer_lbc
      A_IXPTR(94)=A_IXPTR(93)+ MAX(TR_VARS,1)        ! jtracer_lbc_tend
!     tracer_lbc (93) and tracer_lbc_tend (94) have lengths of TR_VARS
!     which may be zero so use MAX(TR_VARS,1) for lengths.

      A_IXPTR(95)=A_IXPTR(94)+ MAX(TR_VARS,1)        ! jtppsozone

!     tropopause ozone (95) has TPPS_OZONE_LEVELS which may be zero
!     so use MAX(TPPS_OZONE_LEVELS,1) for length.

      A_IXPTR(96)=A_IXPTR(95)+ MAX(TPPS_OZONE_LEVELS,1) ! jbmass_new
      A_IXPTR(97)=A_IXPTR(96)+ MODEL_LEVELS             ! jbmass_agd
      A_IXPTR(98)=A_IXPTR(97)+ MODEL_LEVELS             ! jbmass_cld
      A_IXPTR(99)=A_IXPTR(98)+ MODEL_LEVELS          ! jqcf2
      A_IXPTR(100)=A_IXPTR(99)  + WET_LEVELS         ! jqrain
      A_IXPTR(101)=A_IXPTR(100) + WET_LEVELS         ! jqgraup
      A_IXPTR(102)=A_IXPTR(101) + WET_LEVELS         ! jqcf2_lbc
      A_IXPTR(103)=A_IXPTR(102) + 1                  ! jqcf2_lbc_tend
      A_IXPTR(104)=A_IXPTR(103) + 1                  ! jqrain_lbc
      A_IXPTR(105)=A_IXPTR(104) + 1                  ! jqrain_lbc_tend
      A_IXPTR(106)=A_IXPTR(105) + 1                  ! jqgraup_lbc
      A_IXPTR(107)=A_IXPTR(106) + 1                  ! jqgraup_lbc_tend
      A_IXPTR(108)=A_IXPTR(107) + 1                  ! jdust_div1
      A_IXPTR(109)=A_IXPTR(108) + MODEL_LEVELS       ! jdust_div2
      A_IXPTR(110)=A_IXPTR(109) + MODEL_LEVELS       ! jdust_div3
      A_IXPTR(111)=A_IXPTR(110) + MODEL_LEVELS       ! jdust_div4
      A_IXPTR(112)=A_IXPTR(111) + MODEL_LEVELS       ! jdust_div5
      A_IXPTR(113)=A_IXPTR(112) + MODEL_LEVELS       ! jdust_div6
      A_IXPTR(114)=A_IXPTR(113) + MODEL_LEVELS     ! jcf_bulk_lbc
      A_IXPTR(115)=A_IXPTR(114) + 1                ! jcf_bulk_lbc_tend
      A_IXPTR(116)=A_IXPTR(115) + 1                ! jcf_liquid_lbc
      A_IXPTR(117)=A_IXPTR(116) + 1                ! jcf_liquid_lbc_tend
      A_IXPTR(118)=A_IXPTR(117) + 1                ! jcf_frozen_lbc
      A_IXPTR(119)=A_IXPTR(118) + 1                ! jcf_frozen_lbc_tend
      A_IXPTR(120)=A_IXPTR(119) + 1                ! jch4_stoch
      A_IXPTR(121)=A_IXPTR(120) + MODEL_LEVELS     ! jo3_stoch
      A_IXPTR(122)=A_IXPTR(121) + MODEL_LEVELS     ! jdirpar
      A_IXPTR(123)=A_IXPTR(122) + 1                     ! jtr_ukca
      A_IXPTR(124)=A_IXPTR(123) + TR_LEVELS*(TR_UKCA+1) ! jmurk_lbc
      A_IXPTR(125)=A_IXPTR(124) + 1                  ! jmurk_lbc_tend
      A_IXPTR(126)=A_IXPTR(125) + 1                  ! jOH_UKCA
      A_IXPTR(127)=A_IXPTR(126) + MODEL_LEVELS       ! jHO2_UKCA
      A_IXPTR(128)=A_IXPTR(127) + MODEL_LEVELS       ! jH2O2_UKCA
      A_IXPTR(129)=A_IXPTR(128) + MODEL_LEVELS       ! jO3_UKCA
      A_IXPTR(130)=A_IXPTR(129) + MODEL_LEVELS       ! jlcbase
      A_IXPTR(131)=A_IXPTR(130) + 1                  ! jccw_rad 
      A_IXPTR(132)=A_IXPTR(131) + WET_LEVELS         ! jozone_tracer
      A_IXPTR(133)=A_IXPTR(132) + MODEL_LEVELS       ! jO3_prod_loss
      A_IXPTR(134)=A_IXPTR(133) + MODEL_LEVELS       ! jO3_P_L_VMR
      A_IXPTR(135)=A_IXPTR(134) + MODEL_LEVELS       ! jO3_VMR
      A_IXPTR(136)=A_IXPTR(135) + MODEL_LEVELS       ! jO3_P_L_temp
      A_IXPTR(137)=A_IXPTR(136) + MODEL_LEVELS       ! jO3_temp
      A_IXPTR(138)=A_IXPTR(137) + MODEL_LEVELS       ! jO3_P_L_colO3
      A_IXPTR(139)=A_IXPTR(138) + MODEL_LEVELS       ! jO3_colO3
      A_IXPTR(140)=A_IXPTR(139) + MODEL_LEVELS       ! jarclbiog_bg
      A_IXPTR(141)=A_IXPTR(140) + MODEL_LEVELS       ! jarclbiom_fr
      A_IXPTR(142)=A_IXPTR(141) + MODEL_LEVELS       ! jarclbiom_ag
      A_IXPTR(143)=A_IXPTR(142) + MODEL_LEVELS       ! jarclbiom_ic
      A_IXPTR(144)=A_IXPTR(143) + MODEL_LEVELS       ! jarclblck_fr
      A_IXPTR(145)=A_IXPTR(144) + MODEL_LEVELS       ! jarclblck_ag
      A_IXPTR(146)=A_IXPTR(145) + MODEL_LEVELS       ! jarclsslt_fi
      A_IXPTR(147)=A_IXPTR(146) + MODEL_LEVELS       ! jarclsslt_je
      A_IXPTR(148)=A_IXPTR(147) + MODEL_LEVELS       ! jarclsulp_ac
      A_IXPTR(149)=A_IXPTR(148) + MODEL_LEVELS       ! jarclsulp_ak
      A_IXPTR(150)=A_IXPTR(149) + MODEL_LEVELS       ! jarclsulp_di
      A_IXPTR(151)=A_IXPTR(150) + MODEL_LEVELS       ! jarcldust_b1
      A_IXPTR(152)=A_IXPTR(151) + MODEL_LEVELS       ! jarcldust_b2
      A_IXPTR(153)=A_IXPTR(152) + MODEL_LEVELS       ! jarcldust_b3
      A_IXPTR(154)=A_IXPTR(153) + MODEL_LEVELS       ! jarcldust_b4
      A_IXPTR(155)=A_IXPTR(154) + MODEL_LEVELS       ! jarcldust_b5
      A_IXPTR(156)=A_IXPTR(155) + MODEL_LEVELS       ! jarcldust_b6
      A_IXPTR(157)=A_IXPTR(156) + MODEL_LEVELS       ! jarclocff_fr
      A_IXPTR(158)=A_IXPTR(157) + MODEL_LEVELS       ! jarclocff_ag
      A_IXPTR(159)=A_IXPTR(158) + MODEL_LEVELS       ! jarclocff_ic
      A_IXPTR(160)=A_IXPTR(159) + MODEL_LEVELS       ! jarcldlta_dl
      A_IXPTR(161)=A_IXPTR(160) + MODEL_LEVELS       ! jocff_new
      A_IXPTR(162)=A_IXPTR(161) + MODEL_LEVELS       ! jocff_agd
      A_IXPTR(163)=A_IXPTR(162) + MODEL_LEVELS       ! jocff_cld
      A_IXPTR(164)=A_IXPTR(163) + SM_LEVELS          ! jtsoil_tile 
      A_IXPTR(165)=A_IXPTR(164) + SM_LEVELS          ! jsmcl_tile 
!!! Not used MRD
!!!      A_IXPTR(166)=A_IXPTR(165) + SM_LEVELS          ! jsthu_tile 
      A_IXPTR(166)=A_IXPTR(165) + SM_LEVELS          ! jsthf_tile 
      A_IXPTR(167)=A_IXPTR(166) + 3                  ! jsnow_depth3l 
      A_IXPTR(168)=A_IXPTR(167) + 3                  ! jsnow_mass3l 
      A_IXPTR(169)=A_IXPTR(168) + 3                  ! jsnow_tmp3l 
      A_IXPTR(170)=A_IXPTR(169) + 3                  ! jsnow_rho3l 
      A_IXPTR(171)=A_IXPTR(170) + 1                  ! jsnow_rho1l  
      A_IXPTR(172)=A_IXPTR(171) + 1                  ! jsnow_age 
      A_IXPTR(173)=A_IXPTR(172) + 1                  ! jsnow_flg3l 
! Lestevens Sept 2012 - adding new prognostics
      A_IXPTR(174)=A_IXPTR(173) + 10                 ! jcpool_tile
      A_IXPTR(175)=A_IXPTR(174) + 10                 ! jnpool_tile
      A_IXPTR(176)=A_IXPTR(175) + 12                 ! jppool_tile
      A_IXPTR(177)=A_IXPTR(176) + 1                  ! jsoil_order
      A_IXPTR(178)=A_IXPTR(177) + 1                  ! jnidep
      A_IXPTR(179)=A_IXPTR(178) + 1                  ! jnifix
      A_IXPTR(180)=A_IXPTR(179) + 1                  ! jpwea
      A_IXPTR(181)=A_IXPTR(180) + 1                  ! jpdust
      A_IXPTR(182)=A_IXPTR(181) + 1                  ! jglai
      A_IXPTR(183)=A_IXPTR(182) + 1                  ! jphenphase
!CABLE_LAI
      A_IXPTR(184)=A_IXPTR(183) + 1                  ! jCABLE_LAI 
!     Super array length. If length of last variable could be zero, use
!     max(nlevs,1) to avoid out of bounds warning messages

!CABLE_LAI: updated starting index + model_levels
      A_SPPTR_LEN  =A_IXPTR(184) + MODEL_LEVELS
      A_SPPTR_LEN  =A_SPPTR_LEN -1
!L
!L
!L 1.3   Derived constants super array
!L
!L          super array addresses
!   The size of the increment on each line is actually the size of
!   the array referenced by the previous line,
!   e.g. (model_levels+1) is the size of eta_theta_levels.
      A_IXCON(1) = 1                         ! land_index
      A_IXCON(2) =A_IXCON(1) + land_field  ! land_ice_index
      A_IXCON(3) =A_IXCON(2) + land_field  ! soil_index
      ! other length info is now dealt w/ by ALLOCATE staments in SETCONA

      A_SPCON_LEN = A_IXCON(3) + land_field
!     A_SPCON_LEN  =A_SPCON_LEN -1
! sza: Above statement cause array out side of bounds if land_field = 0
! SWarbrick Sean suggested following, also fix typlndm.h in subroutine initdump  
      IF (A_SPCON_LEN > 1 ) A_SPCON_LEN = A_SPCON_LEN -1
!L
!L
!L 1.4   Interface output (boundary conditions) super array
!L
!L          super array addresses
      A_IXINF(1) =1
      A_IXINF(2) =A_IXINF(1) + LEN_FIXHD*N_INTF_A
      A_IXINF(3) =A_IXINF(2) + PP_LEN_INTHD*N_INTF_A
      A_IXINF(4) =A_IXINF(3) + LEN1_LOOKUP*INTF_LOOKUPSA*N_INTF_A
      A_IXINF(5) =A_IXINF(4) + PP_LEN_REALHD*N_INTF_A
      A_IXINF(6) =A_IXINF(5) +                                          &
     &            MAX_INTF_MODEL_LEVELS * INTF_LEN2_LEVDEPC * N_INTF_A
      A_IXINF(7) =A_IXINF(6) + TOT_LEN_INTFA_P
      A_IXINF(8) =A_IXINF(7) + TOT_LEN_INTFA_P
      A_IXINF(9) =A_IXINF(8) + TOT_LEN_INTFA_U
      A_IXINF(10)=A_IXINF(9) + TOT_LEN_INTFA_U
      A_IXINF(11)=A_IXINF(10)+ TOT_LEN_INTFA_P
      A_IXINF(12)=A_IXINF(11)+ TOT_LEN_INTFA_P
      A_IXINF(13)=A_IXINF(12)+ TOT_LEN_INTFA_P
      A_IXINF(14)=A_IXINF(13)+ TOT_LEN_INTFA_P
      A_IXINF(15)=A_IXINF(14)+ TOT_LEN_INTFA_U
      A_IXINF(16)=A_IXINF(15)+ TOT_LEN_INTFA_U
      A_IXINF(17)=A_IXINF(16)+ TOT_LEN_INTFA_U
      A_IXINF(18)=A_IXINF(17)+ TOT_LEN_INTFA_U
      A_IXINF(19)=A_IXINF(18)+ TOT_LEN_INTFA_U
      A_IXINF(20)=A_IXINF(19)+ TOT_LEN_INTFA_U
      A_IXINF(21)=A_IXINF(20)+ U_FIELD_INTFA
      A_IXINF(22)=A_IXINF(21)+ U_FIELD_INTFA
      A_IXINF(23)=A_IXINF(22)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
      A_IXINF(24)=A_IXINF(23)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
      A_IXINF(25)=A_IXINF(24)+  MAX_INTF_MODEL_LEVELS   *N_INTF_A
      A_IXINF(26)=A_IXINF(25)+  MAX_INTF_MODEL_LEVELS   *N_INTF_A
      A_IXINF(27)=A_IXINF(26)+ (MAX_INTF_MODEL_LEVELS+1)*N_INTF_A
!L
      A_IXINF(28)=A_IXINF(27)+ (MAX_INTF_MODEL_LEVELS*N_INTF_A)
      A_IXINF(29)=A_IXINF(28)+ TOT_LEN_INTFA_U
      A_IXINF(30)=A_IXINF(29)+ TOT_LEN_INTFA_U
      A_IXINF(31)=A_IXINF(30)+ TOT_LEN_INTFA_U
      A_IXINF(32)=A_IXINF(31)+ TOT_LEN_INTFA_U
      A_IXINF(33)=A_IXINF(32)+ TOT_LEN_INTFA_U
      A_IXINF(34)=A_IXINF(33)+ TOT_LEN_INTFA_U
      A_IXINF(35)=A_IXINF(34)+                                          &
     &            MAX_LBCROWS * INTF_LEN2_ROWDEPC * N_INTF_A
!L          super array length
!            [MAX( ,1) to prevent harmless out-of-bounds warnings
!                 when no boundary condition interfaces requested.]
      A_SPINF_LEN  =A_IXINF(35)+                                        &
     &         MAX(MAX_LBCROW_LENGTH * INTF_LEN2_COLDEPC * N_INTF_A ,1)

      A_SPINF_LEN  =A_SPINF_LEN -1

!L
!L 1.5   Ancillary file super array
!L
!L          super array addresses
      A_IXANC(1) =1
      A_IXANC(2) =A_IXANC(1) + LEN_FIXHD*NANCIL_DATASETSA
      A_IXANC(3) =A_IXANC(2) + A_LEN_INTHD*NANCIL_DATASETSA
      A_IXANC(4) =A_IXANC(3) + LEN1_LOOKUP*NANCIL_LOOKUPSA
!L
!L          super array length
      A_SPANC_LEN  =A_IXANC(4)+ A_LEN_REALHD*NANCIL_DATASETSA
      A_SPANC_LEN  =A_SPANC_LEN -1
!L
!L
!L 1.6   Input boundary constants super array
!L
!L          super array addresses
      A_IXBND(1) =1
      A_IXBND(2) =A_IXBND(1) + LEN_FIXHD
      A_IXBND(3) =A_IXBND(2) + A_LEN_INTHD
      A_IXBND(4) =A_IXBND(3) + LEN1_LOOKUP*RIM_LOOKUPSA
      A_IXBND(5) =A_IXBND(4) + LEN1_LBC_COMP_LOOKUP*BOUND_LOOKUPSA
!L
!L          super array length
       A_SPBND_LEN  =A_IXBND(5)+ A_LEN_REALHD
      A_SPBND_LEN  =A_SPBND_LEN -1
!L
!L----------------------------------------------------------------------
!L
      RETURN
      END SUBROUTINE UM_INDEX_A
#endif
#if defined(CONTROL) && defined(OCEAN)
!LL  Subroutine: UM_INDEX_O---------------------------------------------
!LL
!LL  Purpose: Calculate addresses and lengths of ocean super arrays
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history:
!LL version  date
!LL  3.2   30/03/93  Introduced as new DECK to allow dynamic allocation
!LL                  of main data arrays in U_MODEL.
!LL  4.0  06/09/95  Added ocean stash superarray. K Rogers
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!    5.4  02/09/02  Added 12th element of o_ixsts           K.Williams
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
#endif
