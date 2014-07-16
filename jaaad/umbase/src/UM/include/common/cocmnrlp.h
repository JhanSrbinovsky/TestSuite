!LL COMDECK COCMNRLP ---------------------------------------------
!LL
!LL Contains start point of rigid lid pressure integration for
!LL ocean model
!LL
!LL History:
!LL
!LL   Model    Date     Modification history
!LL  version
!LL   5.2    20/08/00   New COMDECK created. A. Hines
!     6.1    27/07/04   Add MEAN_RLP. A. Hines
!     6.2    13/02/06   Changes for SOR version of RLP.  B. Ingleby
!LL---------------------------------------------------------------
      REAL MEAN_RLP ! Spatial mean rigid lid pressure
      INTEGER MAX_RLP_BASIN
      PARAMETER (MAX_RLP_BASIN=10)
      INTEGER O_RLP_ILON0(MAX_RLP_BASIN),O_RLP_JLAT0(MAX_RLP_BASIN)
      INTEGER O_RLP_ITT  ! Call RLP routine every ITT'th timestep
      REAL RLP_RJAC   ! changes SOR convergence - depends on model grid
      REAL RLP_EPS    ! degree of SOR convergence required
!
      COMMON/MNRLP/MEAN_RLP,O_RLP_ILON0,O_RLP_JLAT0,O_RLP_ITT,          &
     &  RLP_RJAC,RLP_EPS
!LL---------------------------------------------------------------
