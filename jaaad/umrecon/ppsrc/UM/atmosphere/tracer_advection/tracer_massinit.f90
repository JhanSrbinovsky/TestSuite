!!!Routine to calculate the initial mass of CO2/tracers (global total).
!!!kdcorbin, 02/10

SUBROUTINE TRACER_MASSINIT(    &
    &       row_length, rows, model_levels         &
    &      ,tr_levels, tr_vars                     &
    &      ,halo_i,halo_j,offx,offy,mype,timestep  &
    &      ,r_theta_levels,r_rho_levels,exner_theta_levels   &
    &      ,FV_cos_theta_latitude,delta_lambda,delta_phi     &
    &      ,RHO,q,qcl,qcf                          &
!CO2 and tracer mass flags 
    &      ,L_CO2_MASS,L_TRACER_MASS               &
!CO2 and tracer 3-D fields
    &      ,CO2, tracer                            &
    &      ,tmass                          &
    &      )


    IMPLICIT NONE

!-------------------------------
!INPUT variables
!-------------------------------

     INTEGER :: row_length, rows, model_levels
     INTEGER :: tr_levels,tr_vars
     INTEGER :: halo_i,halo_j,offx,offy,mype,timestep
     REAL    :: delta_lambda,delta_phi
     REAL    :: tmass(21)
     LOGICAL :: L_CO2_MASS,L_TRACER_MASS

     REAL                                                 &
   &  R_THETA_LEVELS(1-halo_i:row_length+halo_i,         &
   &           1-halo_j:rows+halo_j,0:model_levels)      &
   &, R_RHO_LEVELS(1-halo_i:row_length+halo_i,           &
   &           1-halo_j:rows+halo_j,0:model_levels)      &
   &, EXNER_THETA_LEVELS(1-offx:row_length+offx,         &
   &           1-offy:rows+offy,model_levels)            &
   &, FV_cos_theta_latitude(1-offx:row_length+offx,      &
   &           1-offy:rows+offy)                         &
   &, RHO(1-offx:row_length+offx,1-offy:rows+offy,       &
   &           model_levels)                             &
   &, q(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j, &
   &           model_levels)                             &
   &, qcl(1-halo_i:row_length+halo_i,                    &
   &           1-halo_j:rows+halo_j,model_levels)        &
   &, qcf(1-halo_i:row_length+halo_i,                    &
   &           1-halo_j:rows+halo_j,model_levels)        &
   &, CO2(1-offx:row_length+offx,                        &
   &           1-offy:rows+offy,model_levels)            &   
   &, tracer(1-offx:row_length+offx,1-offy:rows+offy,    &
   &           tr_levels,tr_vars)

!-----------------------------
!LOCAL variables
!-----------------------------

    REAL :: temptracer(row_length+2*offx,rows+2*offy,model_levels)
    REAL :: field1mass,field2mass,field3mass,drymass
    CHARACTER*30 :: cmessage

    INTEGER :: n_tr

!------------------------------
!Calculate initial co2 mass
!-----------------------------

IF (L_CO2_MASS) THEN

!DEPENDS ON: tracer_mass
   call tracer_mass(      &
     halo_i,halo_j,offx,offy,    &
     offx,offy,1,             &
     row_length,rows,rows,       &
     model_levels,model_levels,model_levels,  &
     r_theta_levels,r_rho_levels,             &
     FV_cos_theta_latitude,delta_lambda,delta_phi,   &
     rho,q,qcl,qcf,              &
     CO2,CO2,CO2,                &
     exner_theta_levels,"Tracer_Mass: Initial CO2 Mass ",   &
     timestep,0,                 &
     field1mass,field2mass,field3mass,drymass)

   tmass(1)=field1mass/drymass

ELSE

   tmass(1)=0.

ENDIF

!--------------------------------
!Calculate initial tracer mass
!--------------------------------

if (L_TRACER_MASS) then
do n_tr=1,tr_vars

   temptracer(:,:,:) = tracer(:,:,1:model_levels,n_tr)
   write(cmessage,990) "Tracer_Mass: Initial T",n_tr," Mass"
   990 format(a22,i2.2,a6) 


   !DEPENDS ON: tracer_mass
   call tracer_mass(        &
       halo_i,halo_j,offx,offy,      &
       offx,offy,1,               &
       row_length,rows,rows,         &
       model_levels,model_levels,model_levels,   &
       r_theta_levels,r_rho_levels,  &
       FV_cos_theta_latitude,delta_lambda,delta_phi,   &
       rho,q,qcl,qcf,                &
       temptracer,temptracer,temptracer,   &
       exner_theta_levels,cmessage,  &
       timestep,0,                   &
       field1mass,field2mass,field3mass,drymass)

    tmass(n_tr+1) = field1mass/drymass

enddo
endif

RETURN

END SUBROUTINE TRACER_MASSINIT

