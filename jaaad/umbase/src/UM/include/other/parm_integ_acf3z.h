!     ------------------------------------------------------------------
!     Module to set parameters for infinite integrals.
!
!!    NPD_REFINEMENT, NPD_PANEL_POINT and NPD_POINT are not
!!    really independent. ideally, NPD_POINT=1+(NPD_PANEL_POINT-1)
!!    *(2**(NPD_REFINEMENT-1)). Other setting will not make full use of
!!    the declared space.
!
      INTEGER                                                           &
     &    NPD_PANEL_POINT                                               &
!           Number of points in panel
     &  , NPD_INTEGRAL                                                  &
!           Number of integrals available
     &  , NPD_POINT                                                     &
!           Max. no. of points in interval
     &  , NPD_PANEL                                                     &
!           Maximum number of panels
     &  , NPD_REFINEMENT
!           Maximum number of refinements
      REAL  (real64) ::                                                 &
     &    P_PANEL_RATIO
!           Size of each panel
!
      PARAMETER(                                                        &
     &    NPD_PANEL_POINT=5                                             &
     &  , NPD_PANEL=20                                                  &
     &  , NPD_REFINEMENT=12                                             &
     &  , NPD_POINT=16387                                               &
     &  , NPD_INTEGRAL=20                                               &
     &  )
      PARAMETER(                                                        &
     &    P_PANEL_RATIO=3.2                                             &
     &  )
!
!     ------------------------------------------------------------------
