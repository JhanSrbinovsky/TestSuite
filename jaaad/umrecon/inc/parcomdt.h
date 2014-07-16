! Data Statements for PARCOMM COMMON block

! History:
!   5.1      22/05/00  New include file created.          P.Burton
!   5.5      07/02/03  Initialise g_at_extremity for SX     E.Leung

#if defined(UTILIO) || defined(FLDIO) || defined(UTILHIST) \
 || defined(FLUXPROC)
      ! As SX uses single processor, PE0 always has data at extremity
      INTEGER, PARAMETER :: NDATA=4*(maxproc+1)
      DATA g_at_extremity / NDATA * .true. /
#endif
      DATA current_decomp_type/-1/  ! set the initial decomposition
!                                   ! to an "unset" value
