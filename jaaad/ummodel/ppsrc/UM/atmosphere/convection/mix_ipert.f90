
! -----------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------
SUBROUTINE MIX_IPERT(npnts, nlev, nbl, ntml, p_layer_boundaries, &
                     exner_layer_centres, dthbydt, dqbydt, flx_init, &
                     thpert, qpert)

! -----------------------------------------------------------------------
! Purpose
! To well-mix convective increments from the initial parcel 
! perturbation in shallow/deep convection throughout the 
! boundary layer.
!
! -----------------------------------------------------------------------

IMPLICIT NONE


!----------------------------------------------------------------------
! Variables with intent in
!----------------------------------------------------------------------

  INTEGER, INTENT(in) :: npnts        ! Number of points

  INTEGER, INTENT(in) :: nlev         ! Number of model levels

  INTEGER, INTENT(in) :: nbl          ! in number of model layers
                                      ! potentially in the boundary layer

  INTEGER, INTENT(in) :: ntml(npnts)  ! number of model layers for which
                                      ! increments are to be well-mixed.

  REAL, INTENT(in)    :: flx_init(npnts) 
                                      ! the initial parcel mass flux
 
  REAL, INTENT(in)    :: thpert(npnts)
                                      ! the initial parcel temperature 
                                      ! perturbation

  REAL, INTENT(in)    :: qpert(npnts) ! in the initial parcel humidity
                                      ! perturbation                      

  REAL, INTENT(in)    :: p_layer_boundaries(npnts,0:nlev)
                                      ! pressure at layer boundaries  

  REAL, INTENT(in)    :: exner_layer_centres(npnts,0:nlev)
                                      ! exner pressure at layer centres

!----------------------------------------------------------------------
! variables with intent inout
!----------------------------------------------------------------------

  REAL, INTENT(inout) :: dthbydt(npnts,nlev)
                                      ! increment to potential
                                      ! temperature due to convection

  REAL, INTENT(inout) :: dqbydt(npnts,nlev)
                                      ! increment to mixing ratio
                                      ! due to convection

!----------------------------------------------------------------------
! Variables that are locally defined
!----------------------------------------------------------------------

  INTEGER :: i,k                      ! loop counters

  REAL    :: delpsum(npnts)           ! summation of model layer thicknesses
                                      ! with height. (pa)

  REAL    :: delpk(npnts,nbl)         ! difference in pressure across a
                                      ! layer (pa)

  REAL    :: delpexsum(npnts)         ! summation of model layer thicknesses
                                      ! multiplied by the exner pressure (pa).

  REAL    :: dthbydt_exdp(npnts)      ! increment to potential temperature
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness and 
                                      ! exner pressure

  REAL    :: dqbydt_dp(npnts)         ! increment to mixing ratio
                                      ! due to intial perturbation at ntml
                                      ! multiplied by the layer thickness

!----------------------------------------------------------------------
! Calculate the layer pressure thickness and sum up.
!----------------------------------------------------------------------

  DO i = 1, npnts
    delpk(i,1)   = p_layer_boundaries(i,0) - p_layer_boundaries(i,1)
    delpsum(i)   = delpk(i,1)
    delpexsum(i) = delpk(i,1) * exner_layer_centres(i,1)
  END DO

  DO k = 2, nbl
    DO i = 1, npnts
      IF (k  <=  ntml(i)) THEN
        delpk(i,k)   = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)
        delpsum(i)   = delpsum(i)   + delpk(i,k)
        delpexsum(i) = delpexsum(i) + delpk(i,k) * exner_layer_centres(i,k)
      ENDIF
    ENDDO
  ENDDO

!----------------------------------------------------------------------
! Calculate the potential temperature and mixing ratio increments due 
! to the initial perturbation multiplied be the appropriate thickness
! nb. the delpk(ntml) in the numerator and denominator cancel
!----------------------------------------------------------------------

  DO i = 1, npnts
    dthbydt_exdp(i) = -flx_init(i) * thpert(i) * exner_layer_centres(i,ntml(i))
    dqbydt_dp(i)    = -flx_init(i) * qpert(i)
    IF (flx_init(i) <= 0.0) THEN
      WRITE(6,*) 'flx_init(i) <= 0.0 ', flx_init(i)
    END IF
  END DO

!----------------------------------------------------------------------
! Mix the increments due to initial perturbation throughout the 
! sub-cloud layer.
!----------------------------------------------------------------------

  DO k = 1, nbl
    DO i = 1, npnts
      IF (k  <=  ntml(i)) THEN
        dthbydt(i,k) = dthbydt(i,k) + dthbydt_exdp(i) / delpexsum(i)
        dqbydt(i,k)  = dqbydt(i,k)  + dqbydt_dp(i)    / delpsum(i)
      ENDIF
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE MIX_IPERT
