


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
      SUBROUTINE R2_COLUMN_DROPLET_CONC(NPD_PROFILE, NPD_LAYER          &
     &   , N_PROFILE, N_LAYER, NCLDS                                    &
     &   , STRAT_LIQ_CLOUD_FRACTION                                     &
     &   , TOTAL_STRAT_LIQ_CLOUD_FRACTION                               &
     &   , CONV_LIQ_CLOUD_FRACTION                                      &
     &   , TOTAL_CONV_LIQ_CLOUD_FRACTION                                &
     &   , N_DROP, D_MASS, DENSITY_AIR                                  &
     &   , NC_DIAG, NC_WEIGHT)
!
!
!
      IMPLICIT NONE
!
!
!
!  Input variables:
!
      INTEGER                                                           &
     &     NPD_PROFILE                                                  &
!             Maximum number of profiles
     &   , NPD_LAYER                                                    &
!             Maximum number of layers
     &   , N_PROFILE                                                    &
!             Number of atmospheric profiles
     &   , N_LAYER                                                      &
!             Number of layers seen in radiation
     &   , NCLDS
!             Number of cloudy layers
!
      REAL                                                              &
     &     STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)             &
!             Stratiform liquid (T>273K) cloud cover in layers
     &   , TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)                  &
!             Total liquid (T>273K) stratiform cloud cover
     &   , CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)              &
!             Convective liquid (T>273K) cloud cover in layers
     &   , TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)                   &
!             Total liquid (T>273K) convective cloud cover
     &   , N_DROP(NPD_PROFILE, NPD_LAYER)                               &
!             Number concentration of cloud droplets (m-3)
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             Mass thickness of layer (kg m-2)
     &   , DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             Air density (kg m-3)
!
!
!  Output variables:
!
      REAL                                                              &
     &     NC_DIAG(NPD_PROFILE)                                         &
!             Column-integrated droplet number diagnostic (m-2)
     &   , NC_WEIGHT(NPD_PROFILE)
!             Weighting factor for column droplet number
!
!
!  Local variables:
!
      INTEGER                                                           &
     &     I, L
!             Loop counters
!
      REAL                                                              &
     &     N_SUM_S                                                      &
     &   , N_SUM_C                                                      &
!             Temporary sums
     &   , NCOL_S                                                       &
     &   , WGT_S                                                        &
     &   , NCOL_C                                                       &
     &   , WGT_C
!             N-column values and weights for stratiform
!             and convective clouds separately.
!
!
!
      DO L=1, N_PROFILE

         N_SUM_S=0.0
         DO I=N_LAYER+1-NCLDS, N_LAYER
            N_SUM_S=N_SUM_S+(STRAT_LIQ_CLOUD_FRACTION(L, I)             &
     &                *N_DROP(L, I)*D_MASS(L, I)/DENSITY_AIR(L, I))
         ENDDO

         N_SUM_C=0.0
         DO I=N_LAYER+1-NCLDS, N_LAYER
            N_SUM_C=N_SUM_C+(CONV_LIQ_CLOUD_FRACTION(L, I)              &
     &                *N_DROP(L, I)*D_MASS(L, I)/DENSITY_AIR(L, I))
         ENDDO

         IF (TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)  >   0.0) THEN
            NCOL_S=N_SUM_S/TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)
            WGT_S=TOTAL_STRAT_LIQ_CLOUD_FRACTION(L)
         ELSE
            NCOL_S=0.0
            WGT_S=0.0
         ENDIF

         IF (TOTAL_CONV_LIQ_CLOUD_FRACTION(L)  >   0.0) THEN
            NCOL_C=N_SUM_C/TOTAL_CONV_LIQ_CLOUD_FRACTION(L)
            WGT_C=TOTAL_CONV_LIQ_CLOUD_FRACTION(L)
         ELSE
            NCOL_C=0.0
            WGT_C=0.0
         ENDIF

         IF ((WGT_S+WGT_C)  >   0.0) THEN
            NC_DIAG(L)=((NCOL_S*WGT_S)+(NCOL_C*WGT_C))/(WGT_S+WGT_C)
            NC_WEIGHT(L)=WGT_S+WGT_C
         ELSE
            NC_DIAG(L)=0.0
            NC_WEIGHT(L)=0.0
         ENDIF

      ENDDO
!
!     Note: weighting is done later in R2_SET_CLOUD_FIELD
!
!
!
      RETURN
      END SUBROUTINE R2_COLUMN_DROPLET_CONC
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
