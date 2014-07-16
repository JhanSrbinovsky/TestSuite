SUBROUTINE TRACER_MASSPRINT(     &
     &      row_length, rows, model_levels                     &
     &     ,tr_levels, tr_vars   &
     &     ,halo_i, halo_j, offx, offy, mype        &
     &     ,timestep, timestep_number         &
     &     ,r_theta_levels, r_rho_levels, exner_theta_levels   &
     &     ,FV_cos_theta_latitude, delta_lambda, delta_phi   &
     &     ,RHO, q, qcl, qcf, tracer    &
     &     ,lossrate,tracer_flux1,tracer_flux2,tracer_flux3          &
     &     ,tracer_flux4,tracer_flux5,tracer_flux6,tracer_flux7             &
     &     ,tracer_flux8,tracer_flux9,tracer_flux10,tracer_flux11           &
     &     ,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15         &
     &     ,tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19         &
     &     ,tracer_flux20         &
     &     )

     IMPLICIT NONE

!--------------------------------
!INPUT variables
!--------------------------------

    INTEGER :: row_length, rows, model_levels
    INTEGER :: tr_levels, tr_vars
    INTEGER :: halo_i, halo_j, offx, offy, mype
    REAL    :: delta_lambda,delta_phi,timestep
    INTEGER :: timestep_number

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
   &, tracer(1-offx:row_length+offx,1-offy:rows+offy,    &
               tr_levels,tr_vars)                        &
   &, lossrate(row_length,rows,tr_levels,tr_vars)

   REAL                                                  &
   &  tracer_flux1(row_length,rows), tracer_flux2(row_length,rows)   &
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

  REAL :: temptracer(row_length,rows,model_levels)
  REAL :: tempflux(row_length,rows)

  REAL :: field1mass,field2mass,field3mass,drymass,emitval
  REAL :: field1loss,field2loss,field3loss

  INTEGER,parameter :: flux_vars=20
  INTEGER :: i,k,n_tr
  CHARACTER*2 :: tnumber,lnumber
  CHARACTER*30 :: cprint

!----------------------------------
!Print out surface fluxes and loss rates
!----------------------------------

  !Write out mass change from sfc emissions and atmospheric loss
  if (tr_vars .gt. flux_vars) then
     n_tr=flux_vars
  else
     n_tr=tr_vars
  endif

  do i=1,n_tr

        Select Case (i)
           Case (1)
              tempflux(:,:) = TRACER_FLUX1(:,:)

           Case (2)
              tempflux(:,:) = TRACER_FLUX2(:,:)

           Case (3)
              tempflux(:,:) = TRACER_FLUX3(:,:)

           Case (4)
              tempflux(:,:) = TRACER_FLUX4(:,:)

           Case (5)
              tempflux(:,:) = TRACER_FLUX5(:,:)

           Case (6)
              tempflux(:,:) = TRACER_FLUX6(:,:)

           Case (7)
              tempflux(:,:) = TRACER_FLUX7(:,:)

           Case (8)
              tempflux(:,:) = TRACER_FLUX8(:,:)

           Case (9)
              tempflux(:,:) = TRACER_FLUX9(:,:)

           Case (10)
              tempflux(:,:) = TRACER_FLUX10(:,:)

           Case (11)
              tempflux(:,:) = TRACER_FLUX11(:,:)

           Case (12)
              tempflux(:,:) = TRACER_FLUX12(:,:)

           Case (13)
              tempflux(:,:) = TRACER_FLUX13(:,:)

           Case (14)
              tempflux(:,:) = TRACER_FLUX14(:,:)

           Case (15)
              tempflux(:,:) = TRACER_FLUX15(:,:)

           Case (16)
              tempflux(:,:) = TRACER_FLUX16(:,:)

           Case (17)
              tempflux(:,:) = TRACER_FLUX17(:,:)

           Case (18)
              tempflux(:,:) = TRACER_FLUX18(:,:)

           Case (19)
              tempflux(:,:) = TRACER_FLUX19(:,:)

           Case (20)
              tempflux(:,:) = TRACER_FLUX20(:,:)

        End Select

        !Calculate ch4 global surface flux emissions
        !DEPENDS ON: tracer_fluxemit
        call tracer_fluxemit(row_length,rows,mype,timestep,  &
           tempflux,offx,offy,FV_cos_theta_latitude,delta_lambda, &
           delta_phi,emitval)  
           
        temptracer=0.
        temptracer(1:row_length,1:rows,1:model_levels) = &
            lossrate(1:row_length,1:rows,1:model_levels,i)

        !Calculate ch4 global mass loss from oxidation
        !DEPENDS ON: tracer_mass
        call tracer_mass(     &
          halo_i,halo_j,offx,offy,  &
          0,0,1,              &
          row_length,rows,rows,     &
          model_levels,model_levels,model_levels,   &
          r_theta_levels,r_rho_levels,              &
          FV_cos_theta_latitude,delta_lambda,delta_phi,   &
          rho, q, qcl, qcf,         &
          temptracer,temptracer,temptracer,         &
          exner_theta_levels,cprint,  &
          99, 99,            &
          field1loss,field2loss,field3loss,drymass) 

        temptracer=0.
        temptracer(1:row_length,1:rows,1:model_levels) = &
            tracer(1:row_length,1:rows,1:model_levels,i)

        !Calculate ch4 global mass loss from oxidation
        !DEPENDS ON: tracer_mass
        call tracer_mass(     &
          halo_i,halo_j,offx,offy,  &
          0,0,1,              &
          row_length,rows,rows,     &
          model_levels,model_levels,model_levels,   &
          r_theta_levels,r_rho_levels,              &
          FV_cos_theta_latitude,delta_lambda,delta_phi,   &
          rho, q, qcl, qcf,         &
          temptracer,temptracer,temptracer,         &
          exner_theta_levels,cprint,  &
          99, 99,            &
          field1mass,field2mass,field3mass,drymass) 

        if (mype .eq. 0) then
           if (i .eq. 1) then
              write(6,*)
              write(6,*) 'Tracer Flux/Mass/Loss in kg per timestep'
           endif
           write(tnumber,989) i 
           write(6,*) 'Flux T',tnumber,': ',emitval
           write(6,*) 'Total Mass/Loss T',tnumber,': ', &
                   field1mass,field1loss

        endif

      !Calculate ch4 global mass and loss from oxidation by level
      !  once per month
      if (mod(timestep_number,1440) == 0) then
         do k=1, tr_levels
  
            temptracer=0.0
            temptracer(1:row_length,1:rows,k) = &
               tracer(1:row_length,1:rows,k,i)

           !DEPENDS ON: tracer_mass
           call tracer_mass(     &
             halo_i,halo_j,offx,offy,  &
             0,0,1,              &
             row_length,rows,rows,     &
             model_levels,model_levels,model_levels,   &
             r_theta_levels,r_rho_levels,              &
             FV_cos_theta_latitude,delta_lambda,delta_phi,   &
             rho, q, qcl, qcf,         &
             temptracer,temptracer,temptracer,         &
             exner_theta_levels,cprint,  &
             99, 99,            &
             field1mass,field2mass,field3mass,drymass) 

           temptracer=0.0
           temptracer(1:row_length,1:rows,k) = &
              lossrate(1:row_length,1:rows,k,i)

           !DEPENDS ON: tracer_mass
           call tracer_mass(     &
             halo_i,halo_j,offx,offy,  &
             0,0,1,              &
             row_length,rows,rows,     &
             model_levels,model_levels,model_levels,   &
             r_theta_levels,r_rho_levels,              &
             FV_cos_theta_latitude,delta_lambda,delta_phi,   &
             rho, q, qcl, qcf,         &
             temptracer,temptracer,temptracer,         &
             exner_theta_levels,cprint,  &
             99, 99,            &
             field1loss,field2loss,field3loss,drymass) 

             if (mype .eq. 0) then
                write(lnumber,989) k
                write(6,*) '  T',tnumber,' mass/loss for level ', &
                       lnumber,': ', field1mass, field1loss
             endif    
           enddo !1,k
        endif   !mod timestep for once per month

    enddo !1,n_tr

989 format(i2.2)



RETURN

END SUBROUTINE TRACER_MASSPRINT