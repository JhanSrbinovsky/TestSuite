


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Ozone set in lower
!                                               levels.
!                                               (J. M. Edwards)
!       4.4             26-09-97                Conv. cloud amount on
!                                               model levs allowed for.
!                                               J.M.Gregory
!       4.5             18-05-98                Provision for treating
!                                               extra (H)(C)FCs
!                                               included.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Move HCFCs to a more
!                                               natural place in the
!                                               code.
!                                               (J. M. Edwards)
!       5.1             06-04-00                Remove the explicit
!                                               limit on the
!                                               concentration of
!                                               water vapour.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             29-05-02                Add diagnostic call
!                                               for column-integrated
!                                               cloud droplet number.
!                                               (A. Jones)
!       5.4             25-04-02                Replace land/sea mask
!                                               with land fraction in
!                                               call to NUMBER_DROPLET
!                                               (A. Jones)
!       5.5             17-02-03                Change I_CLIM_POINTER
!                                               for hp compilation
!                                               (M.Hughes)
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!  6.2   03-11-05   Enable HadGEM1 climatological aerosols. C. F. Durman
!                   Reworked to use switch instead of #defined. R Barnes
!  6.2   25/05/05  Convert compilation into a more universally portable
!                  form. Tom Edwards
!  6.2   15/12/05  Set negative specific humidities to zero.
!                                               (J. Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Old formulation over
!                                               sea-ice removed.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Ground temperature
!                                               set equal to that
!                                               in the middle of the
!                                               bottom layer.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.3             25-04-01                Alter the specification
!                                               of temperature on rho
!                                               levels (layer
!                                               boundaries).  S.Cusack
!       5.3             25-04-01   Gather land, sea and
!                                  sea-ice temperatures and
!                                  land fraction. Replace TFS
!                                  with general sea temperature.
!                                       (N. Gedney)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                New flag L_AEROSOL_CCN
!                                               introduced to allow
!                                               inclusion of indirect
!                                               aerosol forcing alone.
!                                               Correction of comments
!                                               for LCCWC1 and LCCWC2.
!                                               Correction of level at
!                                               which temperature for
!                                               partitioning
!                                               convective homogeneously
!                                               mixed cloud is taken.
!                                               (J. M. Edwards)
!       4.4             08-04-97                Changes for new precip
!                                               scheme (qCF prognostic)
!                                               (A. C. Bushell)
!       4.4             15-09-97                A parametrization of
!                                               ice crystals with a
!                                               temperature dependedence
!                                               of the size has been
!                                               added.
!                                               Explicit checking of
!                                               the sizes of particles
!                                               for the domain of
!                                               validity of the para-
!                                               metrization has been
!                                               added.
!                                               (J. M. Edwards)
!       5.0             15-04-98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  A.C.Bushell
!       4.5             18-05-98                New option for
!                                               partitioning between
!                                               ice and water in
!                                               convective cloud
!                                               included.
!                                               (J. M. Edwards)
!       4.5             13/05/98   Changes to R2_SET_CLOUD_FIELD to use
!                                  original sect 9 cloud fraction when
!                                  an extended 'area' cloud fraction is
!                                  used everywhere else in Radiation.
!                                  S. Cusack
!       5.1             04-04-00                Remove obsolete tests
!                                               for convective cloud
!                                               and removal of very
!                                               thin cloud (no longer
!                                               required with current
!                                               solvers, but affects
!                                               bit-comparison).
!                                               (J. M. Edwards)
!       5.1             06-04-00                Correct some comments
!                                               and error messages.
!                                               (J. M. Edwards)
!       5.2             10-11-00                With local partitioning
!                                               of convective cloud
!                                               between water and ice
!                                               force homogeneous
!                                               nucleation at
!                                               -40 Celsius.
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Pass sea-salt variables
!                                               down to R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       5.3             11-10-01                Convert returned
!                                               diagnostics to 2-D
!                                               arrays.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check on small cloud
!                                               fractions and liquid for
!                                               contents for PC2 scheme.
!                                               (D. Wilson)
!       5.5             24-02-03                Addition of new ice
!                                               aggregate
!                                               parametrization.
!                                               (J. M. Edwards)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               R2_RE_MRF_UMIST.
!                                               (A. Jones)
!       6.1             07-04-04                Add variables for
!                                               column-droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Pass through PC2 logical
!                                               (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code to check the
!                                               range of validity of
!                                               parametrizations
!                                               added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Error message for
!                                               ice corrected.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             12-06-96                Code rewritten to
!                                               include two types
!                                               of sulphate provided
!                                               by the sulphur cycle.
!                                               (J. M. Edwards)
!       4.2             08-08-96                Climatological aerosol
!                                               model added.
!                                               (J. M. Edwards)
!       4.4             15-09-97                Code for aerosols
!                                               generalized to allow
!                                               arbitrary combinations.
!                                               (J. M. Edwards)
!       4.5   April 1998   Option to use interactive soot in place
!                          of climatological soot.     Luke Robinson.
!                          (Repositioned more logically at 5.1)
!                                                      J. M. Edwards
!       5.1             11-04-00                The boundary layer
!                                               depth is passed to
!                                               the routine setting
!                                               the climatological
!                                               aerosol.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Set up sea-salt aerosol
!                                               (and deactivate climat-
!                                               ological sea-salt) if
!                                               required.
!                                               (A. Jones)
!       5.3             16-10-01                Switch off the
!                                               climatological
!                                               water soluble aerosol
!                                               when the sulphur
!                                               cycle is on.
!                                               (J. M. Edwards)
!       5.3             04-04-01                Include mesoscale
!                                               aerosols if required.
!                                                            S. Cusack
!       5.4             09-05-02                Use L_USE_SOOT_DIRECT
!                                               to govern the extra-
!                                               polation of aerosol
!                                               soot to the extra top
!                                               layer if required.
!                                               (A. Jones)
!       5.5             05-02-03                Include biomass aerosol
!                                               if required.  P Davison
!       5.5             10-01-03                Revision to sea-salt
!                                               density parameter.
!                                               (A. Jones)
!       5.5             21-02-03                Include mineral dust
!                                               aerosol if required.
!                                               S Woodward
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             29-09-97                Original Code
!                                               very closely based on
!                                               previous versions of
!                                               this scheme.
!                                               (J. M. Edwards)
!  4.5  12/05/98  Swap loop order in final nest of loops to
!                 improve vectorization.  RBarnes@ecmwf.int
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Impose a minimum
!                                               thickness on the
!                                               layer filled by the
!                                               boundary layer
!                                               aerosol and correct
!                                               the dimensioning of T.
!                                               (J. M. Edwards)
!       5.3             17-10-01                Restrict the height
!                                               of the BL to ensure
!                                               that at least one layer
!                                               lies in the free
!                                               troposphere in all
!                                               cases. This is needed
!                                               to allow for changes
!                                               in the range of the
!                                               tropopause.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code added for coherent
!                                               convective cloud.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.4             22-07-02                Check that cloud
!                                               fraction is between 0
!                                               and 1 for PC2 scheme.
!                                               (D. Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.4             15-09-97                Accumulation-mode
!                                               and dissolved sulphate
!                                               passed directly to
!                                               this routine to allow
!                                               the indirect effect to
!                                               be used without
!                                               aerosols being needed
!                                               in the spectral file.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete bounds on
!                                               effective radius
!                                               removed.
!                                               (J. M. Edwards)
!       5.2             14-11-00                Add provision for an
!                                               extra radiative layer
!                                               above the top of the
!                                               model.
!                                               (J. M. Edwards)
!       5.2             15-11-00                Subroutine for droplet
!                                               number concentration
!                                               replaced by new function
!                                               with option to use sea-
!                                               salt to supplement
!                                               sulphate aerosol.
!                                               Treatment of convective
!                                               cloud effective radii
!                                               updated.
!                                               (A. Jones)
!       6.1             07-04-04                Add biomass smoke
!                                               aerosol to call to
!                                               NUMBER_DROPLET.
!                                               (A. Jones)
!       6.1             07-04-04                Add new variables for
!                                               column cloud droplet
!                                               calculation.
!                                               (A. Jones)
!       6.2             24-11-05                Pass Ntot_land and
!                                               Ntot_sea from UMUI.
!                                               (Damian Wilson)
!       6.2             02-03-05                Protect calculations
!                                               from failure in PC2
!                                               (Damian Wilson)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.4             29-05-02                Original Code
!                                               (A. Jones)
!       6.1             07-04-04                Modified in accordance
!                                               with AVHRR retrievals:
!                                               only clouds >273K used,
!                                               convective clouds also
!                                               included.
!                                               (A. Jones)
!
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             04-03-96                Original Code
!                                               (J. M. Edwards)
!                                               Parts of this code are
!                                               rather redundant. The
!                                               form of writing is for
!                                               near consistency with
!                                               HADAM3.
!
!       4.5   April 1998   Check for inconsistencies between soot
!                          spectral file and options used. L Robinson.
!       5.3     04/04/01   Include mesoscale aerosol switch when
!                          checking if aerosols are required.  S. Cusack
!       5.4     09/05/02   Include logical flag for sea-salt aerosol.
!                                                              A. Jones
!       5.5     05/02/03   Include logical for biomass smoke
!                          aerosol.               P Davison
! Description of Code:
!   5.5    21/02/03 Add logical for d mineral dust
!                                                S Woodward
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COMPARE_PROC(IERR, L_PRESENT                        &
     &   , L_RAYLEIGH_PERMITTED, L_GAS_PERMITTED, L_CONTINUUM_PERMITTED &
     &   , L_DROP_PERMITTED, L_AEROSOL_PERMITTED                        &
     &   , L_AEROSOL_CCN_PERMITTED, L_ICE_PERMITTED                     &
     &    ,L_USE_DUST, L_USE_BIOGENIC                                   &
     &   , L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT                     &
     &    ,L_USE_SEASALT_DIRECT                                         &
     &    ,L_USE_SOOT_DIRECT                                            &
     &    ,L_USE_BMASS_DIRECT                                           &
     &   , L_USE_OCFF_DIRECT                                            &
     &   , L_CLIMAT_AEROSOL                                             &
     &   , L_USE_ARCL, N_ARCL_SPECIES                                   &
     &   , L_MURK_RAD                                                   &
     &   , L_RAYLEIGH, L_GAS, L_CONTINUUM                               &
     &   , L_DROP, L_AEROSOL, L_AEROSOL_CCN, L_ICE                      &
     &   , NPD_TYPE                                                     &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
!
!
!     DUMMY ARGUMENTS:
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_TYPE
!             NUMBER OF TYPES OF SPECTRAL DATA
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_PRESENT(0: NPD_TYPE)
!             ARRAY INDICATING BLOCKS OF DATA PRESENT
!             IN THE SPECTRAL FILE.
!
!     PROCESSES PERMITTED WITHIN THE UNIFIED MODEL.
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_RAYLEIGH_PERMITTED                                         &
!             RAYLEIGH SCATTERING PERMITTED IN THE MODEL
     &   , L_GAS_PERMITTED                                              &
!             GASEOUS ABSORPTION PERMITTED IN THE MODEL
     &   , L_CONTINUUM_PERMITTED                                        &
!             CONTINUUM ABSORPTION PERMITTED IN THE MODEL
     &   , L_DROP_PERMITTED                                             &
!             CLOUD DROPLET EXTINCTION PERMITTED IN THE MODEL
     &   , L_AEROSOL_PERMITTED                                          &
!             AEROSOL EXTINCTION PERMITTED IN THE MODEL
     &   , L_AEROSOL_CCN_PERMITTED                                      &
!             DETERMINATION OF CCN FROM AEROSOLS PERMITTED IN THE MODEL
     &   , L_ICE_PERMITTED
!             ICE EXTINCTION PERMITTED IN THE MODEL
!
!     OPTIONS PASSED IN
      LOGICAL                                                           &
     &     L_USE_SULPC_DIRECT                                           &
!             LOGICAL TO USE SULPHUR CYCLE FOR THE DIRECT EFFECT
     &   , L_USE_SULPC_INDIRECT                                         &
!             LOGICAL TO USE SULPHUR CYCLE FOR THE INDIRECT EFFECT
     &    ,L_USE_DUST                                                   &
                      !logical to use direct  effect of mineral dust
     &   , L_USE_BIOGENIC                                               &
!             LOGICAL TO USE BIOGENIC AEROSOL FOR THE DIRECT EFFECT
     &   , L_USE_SEASALT_DIRECT                                         &
!             LOGICAL TO USE SEA-SALT FOR THE DIRECT EFFECT
     &    ,L_USE_SOOT_DIRECT                                            &
!             LOGICAL TO USE DIRECT RADIATIVE EFFECT DUE TO SOOT
     &    ,L_USE_BMASS_DIRECT                                           &
!             LOGICAL TO USE DIRECT RADIATIVE EFFECT OF BIOMASS SMOKE
     &   , L_USE_OCFF_DIRECT                                            &
!             LOGICAL TO USE DIRECT RADIATIVE EFFECT OF FOSSIL-FUEL OC
     &   , L_CLIMAT_AEROSOL                                             &
!             LOGICAL TO USE CLIMATOLOGICAL AEROSOL MODEL
     &   , L_MURK_RAD
!             LOGICAL TO USE MESOSCALE MODEL AEROSOL
!
! arcl_dim.h
!
! Maximum dimensions for the aerosol climatology for NWP
!

      integer, parameter :: NPD_ARCL_SPECIES = 7
      integer, parameter :: NPD_ARCL_COMPNTS = 20

! end of arcl_dim.h
      LOGICAL                                                           &
     &     L_USE_ARCL(NPD_ARCL_SPECIES)
!             LOGICAL TO USE AEROSOLS FROM THE NWP CLIMATOLOGY
      INTEGER                                                           &
     &     N_ARCL_SPECIES
!             NUMBER OF SPECIES USED IN THE NWP AEROSOL CLIMATOLOGY
!             (EQUALS ZERO IF THE NWP CLIMATOLOGY IS NOT USED)

!     PROCESSES TO BE ENABLED IN THE RUN.
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_RAYLEIGH                                                   &
!             RAYLEIGH SCATTERING TO BE ENABLED IN THE RUN
     &   , L_GAS                                                        &
!             GASEOUS ABSORPTION TO BE ENABLED IN THE RUN
     &   , L_CONTINUUM                                                  &
!             CONTINUUM ABSORPTION TO BE ENABLED IN THE RUN
     &   , L_DROP                                                       &
!             CLOUD DROPLET EXTINCTION TO BE ENABLED IN THE RUN
     &   , L_AEROSOL                                                    &
!             AEROSOL EXTINCTION TO BE ENABLED IN THE RUN
     &   , L_AEROSOL_CCN                                                &
!             DETERMINATION OF CCN FROM AEROSOL TO BE ENABLED IN THE RUN
     &   , L_ICE
!             ICE EXTINCTION TO BE ENABLED IN THE RUN
!
!     LOCAL VARIABLES
!
!     ARRAY INDICES FOR AEROSOL SPECIES IN THE CLIMATOLOGY FOR NWP
!
! arcl_ids.h
!
! Constants used by the aerosol climatology for NWP
!
      !
      ! Available species:
      !
      !   1. Sulphate
      !   2. Mineral dust
      !   3. Sea-salt
      !   4. Black-carbon (also named soot)
      !   5. Biomass-burning
      !   6. Fossil-fuel organic carbon
      !   7. Delta aerosol
      !
      ! Note: When adding a species, increase NPD_ARCL_SPECIES in
      !       arcl_dim.h. NPD_ARCL_SPECIES must be equal to the
      !       largest IP_ARCL_????

      integer, parameter :: IP_ARCL_SULP = 1
      integer, parameter :: IP_ARCL_DUST = 2
      integer, parameter :: IP_ARCL_SSLT = 3
      integer, parameter :: IP_ARCL_BLCK = 4
      integer, parameter :: IP_ARCL_BIOM = 5
      integer, parameter :: IP_ARCL_OCFF = 6
      integer, parameter :: IP_ARCL_DLTA = 7
      
      !
      ! List of components.
      !
      !   The number of components depends on the species:
      !
      !   1. Sulphate: 3 (accumulation, Aitken, dissolved)
      !   2. Mineral dust: 6 size bins
      !   3. Sea-salt: 2 (film and jet)
      !   4. Black-carbon: 2 (fresh and aged)
      !   5. Biomass-burning: 3 (fresh, aged, in-cloud)
      !   6. Fossil-fuel organic carbon: 3 (fresh, aged, in-cloud)
      !   7. Delta aerosol: 1
      !
      ! Note: when adding a component, increase NPD_ARCL_COMPNTS
      !       in arcl_dim.h. NPD_ARCL_COMPNTS must be equal to the
      !       largest IP_ARCL_????_??.
      !       Array i_arcl_compnts_per_species in set_arcl_dimensions()
      !       will also have to be updated.
      !
      
      integer, parameter :: IP_ARCL_SULP_AC = 1
      integer, parameter :: IP_ARCL_SULP_AK = 2
      integer, parameter :: IP_ARCL_SULP_DI = 3
      integer, parameter :: IP_ARCL_DUST_B1 = 4
      integer, parameter :: IP_ARCL_DUST_B2 = 5
      integer, parameter :: IP_ARCL_DUST_B3 = 6
      integer, parameter :: IP_ARCL_DUST_B4 = 7
      integer, parameter :: IP_ARCL_DUST_B5 = 8
      integer, parameter :: IP_ARCL_DUST_B6 = 9
      integer, parameter :: IP_ARCL_SSLT_FI = 10
      integer, parameter :: IP_ARCL_SSLT_JT = 11
      integer, parameter :: IP_ARCL_BLCK_FR = 12
      integer, parameter :: IP_ARCL_BLCK_AG = 13
      integer, parameter :: IP_ARCL_BIOM_FR = 14
      integer, parameter :: IP_ARCL_BIOM_AG = 15
      integer, parameter :: IP_ARCL_BIOM_IC = 16
      integer, parameter :: IP_ARCL_OCFF_FR = 17
      integer, parameter :: IP_ARCL_OCFF_AG = 18
      integer, parameter :: IP_ARCL_OCFF_IC = 19
      integer, parameter :: IP_ARCL_DLTA_DL = 20

! end of arcl_ids.h
!
!
!     EACH OPTICAL PROCESS INCLUDED IN THE RADIATION CODE MAY BE
!     PERMITTED OR DENIED IN THE UNIFIED MODEL, DEPENDING ON THE
!     PRESENCE OF SUPPORTING CODE. TO BE ENABLED IN A RUN AN OPTICAL
!     PROCESS MUST BE PERMITTED IN THE UNIFIED MODEL AND HAVE
!     SUITABLE SPECTRAL DATA.
      L_RAYLEIGH=L_RAYLEIGH_PERMITTED.AND.L_PRESENT(3)
      L_GAS=L_GAS_PERMITTED.AND.L_PRESENT(5)
      L_CONTINUUM=L_CONTINUUM_PERMITTED.AND.L_PRESENT(9)
      L_DROP=L_DROP_PERMITTED.AND.L_PRESENT(10)
      L_ICE=L_ICE_PERMITTED.AND.L_PRESENT(12)
!
!     SET THE CONTROLLING FLAG FOR THE DIRECT RADIATIVE EFFECTS OF
!     AEROSOLS.
      IF (L_AEROSOL_PERMITTED) THEN
!        SET THE FLAG AND THEN CHECK THE SPECTRAL FILE.
         L_AEROSOL = L_USE_SULPC_DIRECT                                 &
     &          .OR. L_USE_DUST                                         &
     &          .OR. L_USE_BIOGENIC                                     &
     &          .OR. L_CLIMAT_AEROSOL                                   &
     &          .OR. L_MURK_RAD                                         &
     &          .OR. L_USE_SOOT_DIRECT                                  &
     &          .OR. L_USE_BMASS_DIRECT                                 &
     &          .OR. L_USE_OCFF_DIRECT                                  &
     &          .OR. L_USE_SEASALT_DIRECT                               &
     &          .OR. (N_ARCL_SPECIES > 0)
         IF (L_AEROSOL.AND.(.NOT.L_PRESENT(11))) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** ERROR: THE SPECTRAL FILE CONTAINS NO DATA '         &
     &         //'FOR AEROSOLS.', 'SUCH DATA ARE REQUIRED FOR THE '     &
     &         //'DIRECT EFFECT.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ELSE
!        CHECK THAT AEROSOLS HAVE NOT BEEN REQUESTED
!        WHEN NOT PERMITTED.
         IF (L_USE_SULPC_DIRECT                                         &
     &   .OR.L_USE_DUST                                                 &
     &   .OR.L_USE_BIOGENIC                                             &
     &   .OR.L_CLIMAT_AEROSOL                                           &
     &   .OR.L_USE_SEASALT_DIRECT                                       &
     &   .OR.L_USE_SOOT_DIRECT                                          &
     &   .OR.L_USE_BMASS_DIRECT                                         &
     &   .OR.L_USE_OCFF_DIRECT                                          &
     &   .OR.(N_ARCL_SPECIES > 0)) THEN
            WRITE(IU_ERR, '(/A, /A)')                                   &
     &         '*** ERROR: THE DIRECT EFFECTS AEROSOLS ARE NOT '        &
     &         , 'PERMITTED IN THIS CONFIGURATION OF THE '              &
     &         //'RADIATION CODE.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
      ENDIF
!
!
!     IF THE AEROSOL CLIMATOLOGY FOR NWP AND INTERACTIVE AEROSOL SCHEMES
!     ARE REQUESTED AT THE SAME TIME, WARN THAT THE RADIATION CODE WILL
!     ONLY CONSIDER AEROSOLS FROM THE CLIMATOLOGY
!
      IF (N_ARCL_SPECIES > 0) THEN
      
        IF (L_USE_ARCL(IP_ARCL_SULP) .AND. L_USE_SULPC_DIRECT) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL SULPHATE '    &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
        
        IF (L_USE_ARCL(IP_ARCL_DUST) .AND. L_USE_DUST) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL DUST '        &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
        
        IF (L_USE_ARCL(IP_ARCL_SSLT) .AND. L_USE_SEASALT_DIRECT) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL SEASALT '     &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
        
        IF (L_USE_ARCL(IP_ARCL_BLCK) .AND. L_USE_SOOT_DIRECT) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL SOOT '        &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
        
        IF (L_USE_ARCL(IP_ARCL_BIOM) .AND. L_USE_BMASS_DIRECT) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL BIOMASS '     &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
        
        IF (L_USE_ARCL(IP_ARCL_OCFF) .AND. L_USE_OCFF_DIRECT) THEN
          WRITE(IU_ERR, '(/A, /A)')                                     &
     &       '*** WARNING: INTERACTIVE AND CLIMATOLOGICAL OCFF '        &
     &       //'COEXIST.', 'THE FORMER WILL BE NEGLECTED IN THE '       &
     &       //'RADIATION CODE.'
        ENDIF
      
      ENDIF ! N_ARCL_SPECIES
!
!     SET THE CONTROLLING FLAG FOR THE INDIRECT EFFECTS OF AEROSOLS.
!     AT PRESENT THIS DEPENDS SOLELY ON THE SULPHUR CYCLE.
      L_AEROSOL_CCN=L_USE_SULPC_INDIRECT
!
!
!
      RETURN
      END SUBROUTINE R2_COMPARE_PROC
