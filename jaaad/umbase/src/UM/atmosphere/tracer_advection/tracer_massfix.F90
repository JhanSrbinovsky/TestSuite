SUBROUTINE TRACER_MASSFIX(     &
     &      row_length, rows, model_levels                     &
     &     ,tr_levels, tr_vars, min_tracer   &
     &     ,halo_i, halo_j, offx, offy, mype &
     &     ,timestep,timestep_number         &
     &     ,r_theta_levels, r_rho_levels, exner_theta_levels   &
     &     ,FV_cos_theta_latitude, delta_lambda, delta_phi   &
     &     ,RHO, q, qcl, qcf    &
!CO2 and emission flags
     &     ,L_CO2_MASS, L_TRACER_MASS   &
!CO2 and tracer 3-D fields
     &     ,CO2, tracer   &
!CO2 and tracer fluxes
     &     ,co2emitmass,tracer_flux1,tracer_flux2,tracer_flux3,tracer_flux4   &
     &     ,tracer_flux5,tracer_flux6,tracer_flux7,tracer_flux8             &
     &     ,tracer_flux9,tracer_flux10,tracer_flux11,tracer_flux12          &
     &     ,tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16         &
     &     ,tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20         &
!Previous CO2 and tracer mass
     &     ,prevtmass                        &
     &     )

     IMPLICIT NONE

!--------------------------------
!INPUT variables
!--------------------------------

    INTEGER :: row_length, rows, model_levels
    INTEGER :: tr_levels, tr_vars, min_tracer
    INTEGER :: halo_i, halo_j, offx, offy, mype, timestep_number
    REAL    :: delta_lambda,delta_phi,timestep
    REAL    :: prevtmass(21)
    LOGICAL :: L_CO2_MASS, L_TRACER_MASS             

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

   REAL                                                  &
   &  co2emitmass                                        &
   &, tracer_flux1(row_length,rows), tracer_flux2(row_length,rows)   &
   &, tracer_flux3(row_length,rows), tracer_flux4(row_length,rows)   &
   &, tracer_flux5(row_length,rows), tracer_flux6(row_length,rows)   &
   &, tracer_flux7(row_length,rows), tracer_flux8(row_length,rows)   &
   &, tracer_flux9(row_length,rows), tracer_flux10(row_length,rows)  &
   &, tracer_flux11(row_length,rows),tracer_flux12(row_length,rows)  &
   &, tracer_flux13(row_length,rows),tracer_flux14(row_length,rows)  &
   &, tracer_flux15(row_length,rows),tracer_flux16(row_length,rows)  &
   &, tracer_flux17(row_length,rows),tracer_flux18(row_length,rows)  &
   &, tracer_flux19(row_length,rows),tracer_flux20(row_length,rows) 

!----------------------------
!LOCAL variables
!----------------------------

  REAL :: temptracer(row_length+2*offx,rows+2*offy,model_levels)
  REAL :: tempflux(row_length,rows)

  REAL :: traceremitval,co2emitppm
  REAL :: field1mass,field2mass,field3mass,drymass
  REAL :: badco2mass,newco2mass,co2fixfactor
  REAL :: temitppm(tr_vars)
  REAL :: badtmass(tr_vars),newtmass(tr_vars),tfixfactor(tr_vars)

  INTEGER,parameter :: flux_vars=20

  INTEGER :: i,n_tr
  CHARACTER*2 :: tnumber

!---------------------------
!Fix CO2 Mass
!---------------------------

IF (L_CO2_MASS) THEN

!Calculate tracer mass and dry atmospheric mass
!DEPENDS ON: tracer_mass
  call tracer_mass(       &
       halo_i, halo_j, offx, offy,   &
       offx, offy, 1,             &
       row_length, rows, rows,       &
       model_levels,model_levels,model_levels,    &
       r_theta_levels,r_rho_levels,               &
       FV_cos_theta_latitude,delta_lambda,delta_phi,   &
       rho, q, qcl, qcf,             &
       CO2,CO2,CO2,                  &
       exner_theta_levels,"Tracer_Mass: Fixing Mass      ", &
       timestep_number,-99,                 &
       field1mass,field2mass,field3mass,drymass)
         
!Convert to ppm
   badco2mass=field1mass/drymass
   co2emitppm=co2emitmass/drymass
   newco2mass=prevtmass(1)+co2emitppm

!Calculate fix faxtor
   co2fixfactor=newco2mass/badco2mass

!Change CO2 concentrations:
   CO2(:,:,:)=CO2(:,:,:)*co2fixfactor

!Write out results to output file
if (mype == 0) then
   write(6,*) ''
   write(6,*) 'TRACER_MASSFIX: Re-setting CO2 mass'
   write(6,*) 'Prev CO2 Mass: ',prevtmass(1),' CO2 emit: ', co2emitppm
   write(6,*) 'Bad CO2 Mass: ',badco2mass,' New CO2 Mass: ',newco2mass
   write(6,*) 'CO2 Fix Factor: ',co2fixfactor
endif

!Update previous mass 
   prevtmass(1) = prevtmass(1)+co2emitppm

ENDIF

!----------------------------------
!Fix Tracer Mass
!----------------------------------

if (L_TRACER_MASS) then
if (tr_vars .gt. flux_vars) then
   write(6,*)'ERROR FIXING TRACER MASS: Too many tracers'
   write(6,*)'   Only fixing mass for ',flux_vars,' tracers'
   n_tr=flux_vars
else
   n_tr=tr_vars
endif

do i=min_tracer,n_tr

   Select Case (i)

    Case (1)
       !Tracer Flux 1
       tempflux(:,:) = TRACER_FLUX1(:,:)

    Case (2)
       !Tracer Flux 2
       tempflux(:,:) = TRACER_FLUX2(:,:)

    Case (3)
       !Tracer Flux 3
       tempflux(:,:) = TRACER_FLUX3(:,:)

    Case (4)
       !Tracer Flux 4
       tempflux(:,:) = TRACER_FLUX4(:,:)

    Case (5)
       !Tracer Flux 5
       tempflux(:,:) = TRACER_FLUX5(:,:)

    Case (6)
       !Tracer Flux 6
       tempflux(:,:) = TRACER_FLUX6(:,:)

    Case (7)
       !Tracer Flux 7
       tempflux(:,:) = TRACER_FLUX7(:,:)

    Case (8)
       !Tracer Flux 8
       tempflux(:,:) = TRACER_FLUX8(:,:)

    Case (9)
       !Tracer Flux 9
       tempflux(:,:) = TRACER_FLUX9(:,:)

    Case (10)
       !Tracer Flux 10
       tempflux(:,:) = TRACER_FLUX10(:,:)

    Case (11)
       !Tracer Flux 11
       tempflux(:,:) = TRACER_FLUX11(:,:)

    Case (12)
       !Tracer Flux 12
       tempflux(:,:) = TRACER_FLUX12(:,:)

    Case (13)
       !Tracer Flux 13
       tempflux(:,:) = TRACER_FLUX13(:,:)

    Case (14)
       !Tracer Flux 14
       tempflux(:,:) = TRACER_FLUX14(:,:)

    Case (15)
       !Tracer Flux 15
       tempflux(:,:) = TRACER_FLUX15(:,:)

    Case (16)
       !Tracer Flux 16
       tempflux(:,:) = TRACER_FLUX16(:,:)

    Case (17)
       !Tracer Flux 17
       tempflux(:,:) = TRACER_FLUX17(:,:)

    Case (18)
       !Tracer Flux 18
       tempflux(:,:) = TRACER_FLUX18(:,:)

    Case (19)
       !Tracer Flux 19
       tempflux(:,:) = TRACER_FLUX19(:,:)

    Case (20)
       !Tracer Flux 20
       tempflux(:,:) = TRACER_FLUX20(:,:)
 
    End Select

    temptracer(:,:,:) = tracer(:,:,1:model_levels,i)
    
    !Calculate tracer mass and dry atmospheric mass
    !DEPENDS ON: tracer_mass
    call tracer_mass(       &
       halo_i, halo_j, offx, offy,   &
       offx, offy, 1,             &
       row_length, rows, rows,       &
       model_levels,model_levels,model_levels,    &
       r_theta_levels,r_rho_levels,               &
       FV_cos_theta_latitude,delta_lambda,delta_phi,   &
       rho, q, qcl, qcf,             &
       temptracer,temptracer,temptracer,                  &
       exner_theta_levels,"Tracer_Mass: Fixing Mass      ", &
       timestep_number,-99,                 &
       field1mass,field2mass,field3mass,drymass)
         
     !Calculate global emissions
     !DEPENDS ON: tracer_fluxemit
     call tracer_fluxemit(row_length, rows, mype, timestep,  &
       tempflux,offx,offy,   &
       FV_cos_theta_latitude, delta_lambda, delta_phi, traceremitval)

     !Check to make sure tracer is being used and mass is not 0
     If (traceremitval .ne. 0.) Then
        !Convert to ppm
        badtmass(i)=field1mass/drymass
        temitppm(i)=traceremitval/drymass
        newtmass(i)=prevtmass(i+1)+temitppm(i)

        !Calculate fix faxtor
        tfixfactor(i)=newtmass(i)/badtmass(i)

        !Change tracer concentrations:
        tracer(:,:,1:model_levels,i)=temptracer(:,:,:)*tfixfactor(i)

        !Write out results to output file
        if (mype == 0) then
            write(tnumber,989) i
            write(6,*) ''
            write(6,*) 'TRACER_MASSFIX: Re-setting tracer ',tnumber,' mass'
            write(6,*) '  Prev: ',prevtmass(i+1),' Emit: ', temitppm(i)
            write(6,*) '  Bad Mass: ',badtmass(i),' New: ',newtmass(i)
            write(6,*) '  Fix Factor: ',tfixfactor(i)
         endif

989  format(i2.2)

     !Update previous mass 
     prevtmass(i+1) = prevtmass(i+1)+temitppm(i)

     Endif   !emitval ne 0.

enddo  !loop over tracers
endif !l_tracer_mass

RETURN

END SUBROUTINE TRACER_MASSFIX