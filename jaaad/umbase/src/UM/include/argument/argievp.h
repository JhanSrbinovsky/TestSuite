!  ARGIEVP Arguments passed around within EVP module
! Current Code Owner: R Hill
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.5   21/02/03  Added this header + Changes necessitated
!                   by removal of SWAP_1D.     R. Hill
!   6.2   30/03/06  Added L_ONOPOLO for use in EVP_SETUP. R. Hill
! Ice model / UM wide variables
     &  L_OCYCLIC,L_ICESSTILT,L_ICYNPOL,L_ONOPOLO,                      &
     &  IMT,IMTM1,IMTM2,JMT,JMTM1,                                      &
     &  ndte,dtts,dte,dter,                                             &
     &  t_dx,t_dxr,t_dy,t_dyr,cs,cst,csr,tng,radius_si,HT_E,HT_N,xymin, &
     &  ocean,icy,icy_uv,icy_uvp,ocean_t_mask,ocean_uv_mask,            &
     &  uice,vice,aice,hice,hsnow,rhoice,rhosnow,rho_water_si,          &
     &  CORIOLIS,UCURRENT,VCURRENT,                                     &
     &  wsx_ice,wsy_ice,isx,isy,                                        &
     &  sstiltx,sstilty,                                                &
     &  waterx,watery,                                                  &
     &  quad_ice_drag,                                                  &
! Internal EVP variables
     &  ecc2,ecci,eccm,eccp,cw,sw,                                      &
     &  eyc, cstar, pstar,                                              &
     &  etan, etas, etae, etaw, prss, umass,fm,                         &
     &  sig11ne, sig11se, sig11sw, sig11nw,                             &
     &  sig12ne, sig12se, sig12sw, sig12nw,                             &
     &  sig22ne, sig22se, sig22sw, sig22nw,                             &
! Temporary EVP variables
     &  a2na, a2sa, a2ea, a2wa, b2n, b2s, b2e, b2w,                     &
     &  t_dx8, t_dy8, edy, edx, eHN, eHE, eHNm, eHEm, HT_N4, HT_E4,     &
     &  h2n, h2s, h2e, h2w, prssn, prsss, prsse, prssw,                 &
     &  Tdter,tdamp2,c2n,d2n,acoef,bcoef,ccoef,dcoef,ecoef,             &
     &  dxtr4,dytr4                                                     &
! ARGIEVP end
