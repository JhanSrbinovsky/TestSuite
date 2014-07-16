#if defined(A09_2A) || (defined(ATMOS) && defined(CONTROL)) \
                    || (defined(ATMOS) && defined(SCMA))
!
! Parameters required for the PC2 cloud scheme.
!
! Version    Date      Modification
!   5.4    16/08/02    Original Code.            Damian Wilson
!   5.5    13/03/03    Minor fix to defs. P.Selwood.
!   6.1    25/05/04    Change pdf_merge_power and pdf_power and
!                      turbulence in convection parameters.
!                                                Damian Wilson
!   6.2    28/2/05     Adjust cloud fraction tolerances
!                                                Damian Wilson
!   6.2    26/07/06    Change defs to include SCMA. R Barnes
!   6.4    02/02/07    Include parameter to control RH dependent
!                      PC2 erosion of cloud
!
      ! Number of iterations in the initiation
      INTEGER,PARAMETER:: INIT_ITERATIONS=10

      ! Tolerance of cloud fraction for resetting of cloud
      ! Cloud_pc2_tol_2 should always be less than cloud_pc2_tol
      REAL,PARAMETER:: CLOUD_PC2_TOL   = 0.005
      REAL,PARAMETER:: CLOUD_PC2_TOL_2 = 0.001

      ! Tolerance of critical relative humidity for initiation of cloud
      REAL,PARAMETER:: RHCRIT_TOL=0.01

      ! Power that is used to weight the two values of G when they are
      ! merged together
      REAL,PARAMETER:: PDF_MERGE_POWER=0.5

      ! Power that describes the way G varies with s near the boundary
      ! of the cloud probability density function. For a "top-hat" it
      ! is equal to 0, for a "triangular" distribution it is equal to 1.
      REAL,PARAMETER:: PDF_POWER=0.0

      ! Parameters that govern the turbulent decrease in width of the
      ! PDF.  (dbs/dt)/bs = (DBSDTBS_TURB_0 + dQc/dt DBSDTBS_TURB_1)
      !                 * exp( - dbsdtbs_exp Q_c / (a_L qsat(T_L)))
      ! dbsdtbs_turb_0 is now set in the UMUI
!     REAL,PARAMETER:: DBSDTBS_TURB_0 = -2.25E-5
      REAL,PARAMETER:: DBSDTBS_TURB_1 = 0.0
      REAL,PARAMETER:: DBSDTBS_CONV   = 0.0
      Real,Parameter:: dbsdtbs_exp    = 10.05

#endif
