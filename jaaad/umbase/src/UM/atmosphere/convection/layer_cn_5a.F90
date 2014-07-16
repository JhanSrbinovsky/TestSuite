#if defined(A05_5A) 
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates layer dependent constants for model layer k

      SUBROUTINE LAYER_CN(K,npnts,nlev                                 &
      ,                   mdet_on, ent_on                              &
      ,                   ntml,ntpar                                   &
      ,                   l_shallow,l_congestus,l_deep                 &
      ,                   bconv,bwk,bwkp1                              &
      ,                   exner_layer_boundaries                       &
      ,                   exner_layer_centres                          &
      ,                   p_layer_boundaries,p_layer_centres           &
      ,                   recip_pstar,rhum, zk, zkp12, zkp1            &
      ,                   thek, qek,qsek, thekp1,qekp1,qsekp1          &
      ,                   thpk,qpk ,ekm14                              &
      ,                   pkp1,delpkp1,exkp1                           &
      ,                   pk,delpk,delpkp12,exk,delexkp1               &
      ,                   delp_uv_k, delp_uv_kp1                       &
      ,                   ekp14,ekp34,amdetk                           &
                         )

      IMPLICIT NONE
 
!----------------------------------------------------------------------
! Description:
!   Calculates layer dependent constants for layer K i.e.
!            pressure
!            layer thickness
!            entrainment coefficients
!            detrainment coefficients
!
! Method:
!   See Unified Model documentation paper 27.
!
! Current Code Owner: Convection code owner 
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.0 programming standards.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

! Vector lengths and loop counters

      Integer, intent(in) ::   &
       k                       & ! Present model layer 
      ,npnts                   & ! Vector length
      ,nlev                      ! Number of model levels

! Switches

      Integer, intent(in) ::   &
       mdet_on                 & ! flag for adaptive mixing detrainment
                                 !  on = 1, off = 0
      ,ent_on                    ! flag for adaptive entrainment
                                 !  on = 1, off = 0

      Integer, intent(in) ::   &
       ntml(npnts)             & ! Number of levels in the surface-based
                                 ! turbulently mixed layer
      ,ntpar(npnts)              ! Top of initial parcel ascent

! Logical Switches

      Logical, intent(in) ::   &
       l_shallow               & ! indicator all points are shallow
      ,l_congestus             & ! indicator all points are congestus
      ,l_deep                    ! indicator all points are deep

      Logical, intent(in) :: &
       bconv(npnts)       & ! mask for points at which convection takes place
      ,bwk(npnts)         & ! mask for points where condensate is liquid on k
      ,bwkp1(npnts)         ! mask for points where condensate is liquid on k+1

! Field on model levels

      Real,  intent(in) ::                   &
       exner_layer_boundaries(npnts,0:nlev)  &
                                       ! Exner function at layer boundary
                                       ! starting at level k-1/2
      ,exner_layer_centres(npnts,0:nlev)     &
                                       ! Exner function at layer  centre
      ,p_layer_centres(npnts,0:nlev)         &
                                       ! Pressure at layer  centre (Pa)
      ,p_layer_boundaries(npnts,0:nlev)      
                                       ! Pressure at layer  boundary (Pa)
! Field on a level
      Real,  intent(in) :: &
       recip_pstar(npnts)  & ! Reciprocal of pstar array (1/Pa)
      ,rhum(npnts)         & ! Relative humidity at level K
      ,zk(npnts)           & ! height on k
      ,zkp12(npnts)        & ! height on k+1/2
      ,zkp1(npnts)         & ! height on k+1
      ,thek(npnts)         & ! theta for environment on k
      ,qek(npnts)          & ! q for environment on k
      ,qsek(npnts)         & ! q sat for environment on k
      ,thekp1(npnts)       & ! theta for environment on k+1
      ,qekp1(npnts)        & ! q for environment on k+1
      ,qsekp1(npnts)       & ! q sat for environment on k+1
      ,thpk(npnts)         & ! parcel theta on k
      ,qpk(npnts)          & ! parcel q on k
      ,ekm14(npnts)          ! ek14 from previous pass i.e. ek for k-1+1/4 

! Information for layer k+1

      Real,  intent(inout) ::  &
       pkp1(npnts)             & ! pressure at layer K+1 (Pa) 
      ,delpkp1(npnts)          & ! thickness of layer K+1 (Pa) 
      ,exkp1(npnts)              ! Exner function at level K+1

! Information for layer K

      Real,  intent(out) ::    &
       pk(npnts)               & ! pressure at layer K (Pa) 
      ,delpk(npnts)            & ! thickness of layer K (Pa) 
      ,delpkp12(npnts)         & ! thickness between layer K & K+1 (Pa) 
      ,exk(npnts)              & ! Exner function at level K
      ,delexkp1(npnts)         & ! Difference in Exner function between
                                 ! K+3/2 and K+1/2
      ,delp_uv_k(npnts)        & ! thickness of uv layer K (Pa)
      ,delp_uv_kp1(npnts)      & ! thickness of uv layer K+1 (Pa)
      ,ekp14(npnts)            & ! entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
      ,ekp34(npnts)            & ! entrainment coefficient at level k+3/4
                                 ! multiplied by appropriate layer thickness
      ,amdetk(npnts)             ! mixing detrainment coefficient at level k
                                 ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
      Integer :: I             ! Counter for do loop
!
      Real ::                &
       Aekp14,Aekp34         & ! Used in calculation of entrainment rate
      ,amdet_fac             & ! factor to scale down amdetk 'cos eps is bigger
      ,RENT                  & ! parameter controlling fractional entrainment
      ,EL                    & ! latent heat of condensation(/freezing)
      ,gzk                   & ! g*z(k)
      ,cpexnerk              & ! cp*exner(k)
      ,gzkp1                 & ! g*z(k+1)
      ,cpexnerkp1_th         & ! cp*exner(k+1)*theta(k)
      ,hsat_ek(npnts)        &
      ,hsat_ekp1(npnts)      & ! saturated moist static energies
      ,hsat_ekp12(npnts)     &
      ,h_pk(npnts)           & ! moist static energy
      ,h_pkp12(npnts)        &
      ,h_ek(npnts)           &
      ,h_ekp1(npnts)         &
      ,h_ekp12(npnts)        &
      ,delta_hsat_env(npnts) & ! hsat_ekp1 - hsat_ek
      ,dhpbydp(npnts)        &
      ,dhebydp(npnts)        &
      ,epsilonp14(npnts)     & ! fractional entrainment coef calc from hsat
      ,epsilonp34(npnts)     & ! fractional entrainment coef calc from hsat
      ,ekp14_ad(npnts)       & ! EPS * layer thickness for 1/4 layer
      ,ekp34_ad(npnts)       & ! EPS * layer thickness for 3/4 layer
      ,delp_cld(npnts)       & ! thickness of cloud layer (Pa)
      ,delpkp14(npnts)       & ! thickness of layer for k+1/4 (Pa)
      ,delpkp34(npnts)         ! thickness of layer for k+3/4 (Pa)
      
!----------------------------------------------------------------------
! Model constants
!----------------------------------------------------------------------
!
#include "entcnst.h"
#include "c_r_cp.h"
#include "c_lheat.h"
#include "c_g.h"

!----------------------------------------------------------------------
! Set constant ae used in calculation of entrainment and detrainment
! rates depending upon level.
!----------------------------------------------------------------------
!
      aekp14 = ae2
      aekp34 = ae2
 
!---------------------------------------------------------------------
! Initialise various arrays holding layer thicknesses etc.
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Calculate pk and delpk - if K = 1 (lowest model layer) then
! values for previous pass through routine at (k-1)+1 are taken.
! Calculate Exner functions at mid-layers K and K+1, and difference
! of exner function across layer K
!---------------------------------------------------------------------

      If (K == 2) Then
        Do I=1,npnts
          pk(i)    = p_layer_centres(i,k)
          delpk(i) = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)
          exk(i)   = exner_layer_centres(i,k)
        End Do 
      Else
        Do I=1,npnts
          pk(i)    = pkp1(i)
          delpk(i) = delpkp1(i)
          exk(i)   = exkp1(i)
        End Do 
      End If

!---------------------------------------------------------------------
! Calculate pkp1, delpkp1 and delpk+1/2
!---------------------------------------------------------------------
      Do I=1,npnts
        pkp1(i)     = p_layer_centres(i,k+1)
        delpkp1(i)  = p_layer_boundaries(i,k) -                      &
                                       p_layer_boundaries(i,k+1)
        delpkp12(i) = pk(i) - pkp1(i)
        exkp1(i)    = exner_layer_centres(i,k+1)
        delexkp1(i) = exner_layer_boundaries(i,k)-                   &
                                     exner_layer_boundaries(i,k+1)

!---------------------------------------------------------------------
! Calculate delp_uv_k and delp_uv_kp1
! NB: Ensure positive by taking - of difference
!---------------------------------------------------------------------
!
        delp_uv_k(i)   = p_layer_centres(i,k-1) - p_layer_centres(i,k)
        delp_uv_kp1(i) = p_layer_centres(i,k) - p_layer_centres(i,k+1)

! For use later in this routine only 
! Cloud thickness
        delp_cld(i)    = p_layer_boundaries(i,ntml(i))               &
                                   -p_layer_boundaries(i,ntpar(i))
! Thickness of layer for k+1/4
        delpkp14(i)    = pk(i) - p_layer_boundaries(i,k)

! Thickness of layer for k+3/4
        delpkp34(i)    = p_layer_boundaries(i,k) - pkp1(i)
      End Do

!----------------------------------------------------------------------
! Adaptive entrainment option 
! UM documentation paper 27 Section 7.3
!----------------------------------------------------------------------

      If (ent_on  ==  1 ) Then

        rent=0.5

! Initialise everything in sight!,  except ekm14, obviously

        Do I=1, npnts
          hsat_ek(i)=0.0
          hsat_ekp1(i)=0.0
          hsat_ekp12(i)=0.0
          h_pk(i)=0.0
          h_pkp12(i)=0.0
          h_ek(i)=0.0
          h_ekp1(i)=0.0
          h_ekp12(i)=0.0
          delta_hsat_env(i)=0.0
          dhpbydp(i)=0.0
          dhebydp(i)=0.0
          epsilonp14(i)=0.0
          epsilonP34(i)=0.0
          ekp14_ad(i)=0.0
          ekp34_ad(i)=0.0
        End Do
        Do  I=1,npnts

!Calculate hsat in case where adaptive entrainment is switched on

          If(bwk(i)) Then
            EL = LC
          Else
            EL = LC+LF
          End If
          gzk       = g*zk(i)
          cpexnerk  = cp * exner_layer_centres(i,k)
          hsat_ek(i)= cpexnerk * thek(i) + EL * qsek(i) + gzk
          h_ek(i)   = cpexnerk * thek(i) + EL * qek(i)  + gzk
          h_pk(i)   = cpexnerk * thpk(i) + EL * qpk(i)  + gzk

          If(bwkp1(i)) Then
            EL = LC
          Else
            EL = LC+LF
          End If
          gzkp1         = g*zkp1(i)
          cpexnerkp1_th = cp * exner_layer_centres(i,k+1)* thekp1(i)
          hsat_ekp1(i)= cpexnerkp1_th  + EL * qsekp1(i) + gzkp1
          h_ekp1(i)   = cpexnerkp1_th  + EL * qekp1(i)  + gzkp1

          delta_hsat_env(i)=(hsat_ekp1(i)-hsat_ek(i))


!version below has stability criterion, but doesn't really help
!          delta_hsat_env(i)=MIN((hsat_ekp1(i)-hsat_ek(i)), 0.0)
!Not sure precisely what to do with this- if delta_hsat_env
!before correction was positive, convection shouldn't be taking
!place at all, never mind the value of the entrainment coefficient.


          hsat_ekp12(i) = hsat_ek(i)+ delpkp14(i)/delpkp12(i) *         &
                                              delta_hsat_env(i)

       !NB already inside loop over npnts

          dhpbydp(i) = (1.0 - RENT) * delta_hsat_env(i)/                &
                                              (-1.0*delpkp12(i))
          h_pkp12(i) = h_pk(i)  -delpkp14(i)* dhpbydp(i)

          dhebydp(i) = (h_ekp1(i) - h_ek(i))/(-1.0*delpkp12(i)) 

          h_ekp12(i) = h_ek(i)  - delpkp14(i)* dhebydp(i)

        End Do    

! ---------------------------------------------------------------------
!  Calculate entrainment coefficientS multiplies by approppriate
!  layer thickness
! ---------------------------------------------------------------------
!
        If (l_shallow) Then
          Do i=1,npnts
          
            If (K  >   ntml(i) ) Then

              ekp14(i) = delpkp14(i) * 0.03                            &
                        * EXP( -1.0*( p_layer_boundaries(I,ntml(i))    &
                                 - pk(i) )/ delp_cld(i) )              &
                                          / delp_cld(i) 
            Else

              ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *      &
                        recip_pstar(i) * recip_pstar(i)
  
            End If    ! level
            If (K  ==  ntml(i)) ekp14(i) = 0.0
            If ( K  >=  ntml(i) ) Then

              ekp34(i) = delpkp34(i) * 0.03                            &
                         * EXP( -1.0*( p_layer_boundaries(I,ntml(i))   &
                           - p_layer_boundaries(i,k) )/ delp_cld(i) )  &
                                            / delp_cld(i)

            Else ! Deep or mid level (or levels not used)        

              ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *&
                        delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)

            End If   ! type of convection and level
          End Do 

        Else If (L_congestus) Then

          Do i=1,npnts
          
            If (K  >   ntml(i) ) Then
! like shallow
              ekp14(i) = delpkp14(i) * 0.03                            &
                        * EXP( -1.0*( p_layer_boundaries(I,ntml(i))    &
                                 - pk(i) ) / delp_cld(i) )             &
                                           / delp_cld(i) 

! 1/Z
!             ekp14(i) = 1. *2. * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i)) 

            Else

              ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *        &
                        recip_pstar(i) * recip_pstar(i)

            End If    ! type of convection and level
            If (K  ==  ntml(i)) ekp14(i) = 0.0
            If ( K  >=  ntml(i) ) Then

! 1/Z possible alternative
!            ekp34(i) = 1. *2.* (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
! like shallow
              ekp34(i) = delpkp34(i) * 0.03                            &
                         * EXP( -1.0*( p_layer_boundaries(I,ntml(i))   &
                           - p_layer_boundaries(i,k) )/ delp_cld(i) )  &
                                            / delp_cld(i)

            Else 
      
              ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *&
                        delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)

            End If   ! type of convection and level
          End Do

        Else ! Deep or mid level 

          Do i=1,npnts 
            ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *        &
                        recip_pstar(i) * recip_pstar(i)
            ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *  &
                        delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)

          End Do

        End If   ! type of convection 

        !-------------------
        !Calculate epsilon14 and epsilon 34
        !-------------------
        Do i=1,npnts 
          If(bconv(i)) Then
            If((h_pk(i) - h_ek(i))  >   0.0) Then
              epsilonp14(i)= dhpbydp(i)/MAX((h_pk(i) - h_ek(i)),500.0)

          !500 chosen with view to not letting buoyancy drop too much
          !approx 0.2 cp + 0.2/1000 L
            Else
               epsilonp14(i)= 0.0
            End If
            If((h_pkp12(i) - h_ekp12(i))  >   0.0) Then
               epsilonp34(i)= dhpbydp(i)/MAX((h_pkp12(i) - h_ekp12(i)),500.0)
            Else
               epsilonp34(i)= 0.0
            End If

            If(((h_pk(i) - h_ek(i))  >   500.0)                        &
                    .and.(dhpbydp(i)  >   1.0e-9)) Then
              ekp14_ad(i) = MIN((epsilonp14(i) *                       &
                                (pk(i) - p_layer_boundaries(i,k))),1.0)
            Else
              ekp14_ad(i) = 0.0
            End If

            If(((h_pkp12(i) - h_ekp12(i))  >  1.0e-10)                 &
                        .and. (dhpbydp(i)  >  1.0e-10)) Then
              ekp34_ad(i) = MIN((epsilonp34(i) *                       &
                                (p_layer_boundaries(i,k)  - pkp1(i))),1.0)

!Check on sign of difference
              If((((ekp14_ad(i) - ekm14(i))    <   0.0) .and.          &
                  ((ekp34_ad(i) - ekp14_ad(i)) >   0.0))               &
               .OR.                                                    &
                 (((ekp14_ad(i) - ekm14(i))    >   0.0) .and.          &
                  ((ekp34_ad(i) - ekp14_ad(i)) <   0.0))) Then

                ekp34_ad(i) = ekp14_ad(i)

              End If

            Else   
              ekp34_ad(i) = 0.0
            End If     ! test  h gradients

!       reset values of entrainment coefficients to adaptive versions
!        ekp14(i)=MAX(ekp14_ad(i),0.0)
!        ekp34(i)=MAX(ekp34_ad(i),0.0)
!Try an average of the adaptive and original versions
!but using GR as minimum value below which  entrainment doesn't fall
            If(ekp14_ad(i)  >   0.0) Then
              ekp14(i)=0.5*(ekp14_ad(i)+ekp14(i))
            End If
            If(ekp34(i)  >   0.0) Then
              ekp34(i)=0.5*(ekp34_ad(i)+ekp34(i))
            End If
          End If        ! bconv test

        End Do  ! npnts

!----------------------------------------------------------------------

      Else     ! end of adaptive entrainment calculation

!---------------------------------------------------------------------
! Original entrainment calculations.
! Calculate entrainment coefficients multiplied by appropriate
! layer thickness.  
!
! (UM Doc 27 section 2C, equation 14)
!---------------------------------------------------------------------

        If (l_shallow) then   ! shallow scheme

          Do i=1,npnts 
            If ( K  >   ntml(i) ) Then
              ekp14(i) = delpkp14(i)* 0.03                              &
                         * EXP( -1.0*( p_layer_boundaries(i,ntml(i))    &
                                        - pk(i) )/ delp_cld(i) )        &
                                           / delp_cld(i)
            Else
              ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *       &
                         recip_pstar(i) * recip_pstar(i)
            End If
            If (K  ==  ntml(i) ) ekp14(i) = 0.0  
            If (K  >=  ntml(i) ) Then
              ekp34(i) = delpkp34(i)* 0.03                              &
                         * EXP( -1.0*( p_layer_boundaries(i,ntml(i))    &
                         - p_layer_boundaries(i,k) )/ delp_cld(i))      &
                                           / delp_cld(i)
            Else
              ekp34(i) = entcoef * Aekp34 * p_layer_boundaries(i,k) *   &
                         delpkp34(i)*                                   &
                         recip_pstar(i) * recip_pstar(i)
            End If
          End Do      ! npnts

        Else If (l_congestus) Then       ! congestus scheme 
        
          Do I=1,npnts
            If ( K  >  ntml(i) ) Then
              ekp14(i) = delpkp14(i)*0.03          &
                   * EXP( -1.0*( p_layer_boundaries(i,ntml(i))- pk(i) )  &
                            /delp_cld(i))/delp_cld(i)

! 1/z  possible alternative
!              ekp14(i) = 1.0*2.0*(zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i)) 

            Else        ! (k < ntml) levels not used by scheme
              ekp14(i) = 0.0
            End If

            If ( K  >=  ntml(i) ) Then
              ekp34(i) = delpkp34(i)*0.03                               &
                          * EXP( -1.0*( p_layer_boundaries(i,ntml(i))   &
                                    - p_layer_boundaries(i,k) )         &
                            /delp_cld(i))/delp_cld(i)
! 1/z 
!            ekp34(i) = 1.0*2.0*(zkp12(i) - zk(i))/(zkp12(i) + zk(i))

            Else       ! (k < ntml) levels not used by scheme
              ekp34(i) =0.0 
            End If
          End Do

        Else if (l_deep) then      ! deep scheme      

          Do i=1,npnts 
            ekp14(i) = entcoef * Aekp14 * pk(i) *delpkp14(i)*           &
                         recip_pstar(i) * recip_pstar(i)
            ekp34(i) = entcoef * Aekp34 * (p_layer_boundaries(i,k)) *   &
                         delpkp34(i) * recip_pstar(i) * recip_pstar(i)
          End Do

        Else                       ! mid-level scheme

          Do i=1,npnts 
            ekp14(i) = entcoef * Aekp14 * pk(i) *delpkp14(i)*           &
                         recip_pstar(i) * recip_pstar(i)
            ekp34(i) = entcoef * Aekp34 * (p_layer_boundaries(i,k)) *   &
                         delpkp34(i) * recip_pstar(i) * recip_pstar(i)
          End Do

        End If     ! test on type of convection
         
       End If ! entrainment off condition

! ---------------------------------------------------------------------
!  Calculate mixing detrainment coefficient multiplied by appropriate
!  layer thickness.
!
!  UM Documentation paper 27, section (2C) equation(15)
! ---------------------------------------------------------------------
 
      If (l_shallow) Then 
         Do i=1,npnts
           If(K == 1)Then
             amdetk(i) = 0.0
           Else If (K  >=  ntml(i) ) Then
              If (rhum(i)  <=  0.85) Then
                amdetk(i) = (1.0 + 0.3)*(ekp14(i) + ekp34(i))
              Else If (rhum(i)  >   1.0) Then
                amdetk(i) = 1.0*(ekp14(i) + ekp34(i))
              Else
                amdetk(i) = (1.0 + (0.3/0.15)*(1.0-rhum(i)))            &
                                         *(ekp14(i) + ekp34(i))
              End If
           Else
             amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
           End If
         End Do !npnts

      Else If (l_congestus) Then
!-----------------------------------------------------------------------
! Adaptive mixing detrainment allowed as an option.
! See UM documentation paper section 7.2 equation 7.3
!-----------------------------------------------------------------------
         
        If (mdet_on == 0 ) then  ! original mixing detrainment  
          Do i=1,npnts
            If(K == 1)Then
              amdetk(i) = 0.0
            Else If (K  >=  ntml(i) ) Then
              IF (rhum(i)  <=  0.85) Then
                amdetk(i) = (1.0 + 0.1)*(ekp14(i) + ekp34(i))
              Else If (rhum(i)  >   1.0) Then
                amdetk(i) = 1.0*(ekp14(i) + ekp34(i))
              Else
                amdetk(i) = (1.0 + (0.1/0.15)*(1.0-rhum(i)))            &
                                         *(ekp14(i) + ekp34(i))
              End If
            Else
             amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
            End If
          End Do !npnts
        Else             ! adaptive mixing detrainment
          Do i=1,npnts
            If (K == 1) Then
              amdetk(i) = 0.0
            Else If (K  >=  ntml(i) .and. rhum(i) <= 1.0 ) Then
              amdetk(i) = (1.0 -rhum(i))*(ekp14(i) + ekp34(i))
            Else
              amdetk(i) = 0.0
            End If
          End Do !npnts
        End If         ! test on mdet_on  

      Else If (l_deep) Then

!-----------------------------------------------------------------------
! Adaptive mixing detrainment allowed as an option.
! See UM documentation paper section 7.2 equation 7.3
!-----------------------------------------------------------------------

        If (K == 1) Then
          Do i=1,npnts
            amdetk(i) = 0.0
          End Do !npnts
        Else
          If (mdet_on  ==  1) Then

            amdet_fac=1.0   ! used as tuning in adaptive mixing detrainment
                            ! Could be removed as currently set to 1.
            Do i=1,npnts
              If (K  >=  ntml(i) ) Then
                IF (rhum(i)  <=  1.0) Then
                  amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i))*(1-rhum(i))
                Else
                  amdetk(i) = 0.0
                End If
              Else          ! original mixing detrainment option
                amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
              End If
            End Do !npnts

          Else     ! not mdet      
            Do i=1,npnts
              amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
            End Do !npnts
          End If     ! test on mdet_on
        End If   ! test on k=1

      Else      ! mid_level scheme  (no adaptive mixing detrainment)

        If (K == 1) Then
          Do i=1,npnts
            amdetk(i) = 0.0
          End Do !npnts
        Else
          Do i=1,npnts
            amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
          End Do !npnts
        End If

      End If     ! type of convection scheme

!
      RETURN
      END SUBROUTINE LAYER_CN
#endif
