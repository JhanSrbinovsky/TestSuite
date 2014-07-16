!!Routine to calculate MCF atmospheric loss due to photolysis 
!! and MCF loss from ocean sink
!
!kdcorbin, 04/10

SUBROUTINE TRACER_MCFLOSS(     &
     &      row_length, rows, model_levels,mcf_trnum  &
     &     ,tr_levels, tr_vars, offx, offy,deltat     &
     &     ,halo_i,halo_j,r_theta_levels,r_rho_levels &
     &     ,theta,pressure,density                    &
     &     ,tracer,oh,jrate,oceanloss,lossrate )

     IMPLICIT NONE

!--------------------------------
!INPUT variables
!--------------------------------

    INTEGER :: row_length, rows, model_levels
    INTEGER :: mcf_trnum, tr_levels, tr_vars
    INTEGER :: offx, offy, halo_i, halo_j
    REAL :: deltat,deltaz

    REAL                                                 &
   &  r_theta_levels(1-halo_i:row_length+halo_i,         &
   &       1-halo_j:rows+halo_j,0:model_levels)          &
   &, r_rho_levels(1-halo_i:row_length+halo_i,           &
   &       1-halo_j:rows+halo_j,1:model_levels)          &
   &, theta(1-offx:row_length+offx,1-offy:rows+offy,     &
   &          model_levels)                              &
   &, pressure(1-offx:row_length+offx,1-offy:rows+offy,  &
   &          model_levels)                              &
   &, density(1-offx:row_length+offx,1-offy:rows+offy,   &
   &          model_levels)                              &
   &, tracer(1-offx:row_length+offx,1-offy:rows+offy,    &
   &           tr_levels,tr_vars)                        &
   &, oh(1:row_length,1:rows,1:model_levels)             &
   &, jrate(1:row_length,1:rows,1:model_levels)          &
   &, oceanloss(1:row_length,1:rows,1:model_levels)      &
   &, lossrate(row_length,rows,tr_levels,tr_vars)

!----------------------------
!LOCAL variables
!----------------------------

  REAL, PARAMETER :: avogadro=6.022E23,mass_mcf=133.4
  REAL, PARAMETER :: pref=100000.,rair=287.05,cp=1005.
  REAL, PARAMETER :: gravity=9.8
  REAL, PARAMETER :: kappa=rair/cp

  REAL :: temperature(row_length,rows,model_levels)
  REAL :: koh(row_length,rows,model_levels)
  REAL :: temptracer(row_length,rows,tr_levels)

  INTEGER :: i,j,k,iz

!----------------------------
!Check for MCF field
!----------------------------
if (mcf_trnum .gt. tr_vars) then
   write(6,*) 'MCF not included in run.'
   write(6,*) 'Num tracers: ',tr_vars,' MCF tracer number: ',mcf_trnum
else

!-----------------------------
!Calculate temperature from theta and pressure
!Calculate OH rate constant
!-----------------------------
do i=1,row_length
  do j=1,rows
     do k=1,model_levels
        temperature(i,j,k) =  &
            theta(i,j,k)/((pref/pressure(i,j,k))**kappa)
        koh(i,j,k) = 1.64E-12*exp(-1520/temperature(i,j,k))
     enddo
   enddo
enddo

!----------------------------
!Calculate mcf loss 
!----------------------------

!!!Tracers have MCF mass mixing ratio, with units of kg [MCF] / kg [AIR]
do i=1,row_length
   do j=1,rows
      deltaz=(r_rho_levels(i,j,2)-r_rho_levels(i,j,1))

      if (deltaz .le. 0.) then
         print*,'WARNING: Bad deltaz in tracer_mcfloss.'
         print*,'Setting to 0.'
         deltaz=0.
      endif

      !kdcorbin, 07/10 - do not use loss on top level
      do k=1,tr_levels-1

            temptracer(i,j,k) = tracer(i,j,k,mcf_trnum)

            !Convert MCF units to kg [MCF] / m^3
            temptracer(i,j,k) = temptracer(i,j,k)*density(i,j,k)

            !Convert MCF units to molecules/cm^3
            temptracer(i,j,k) = temptracer(i,j,k)*avogadro/mass_mcf/1000.

            !Calculate MCF rate of change due to OH in molecules/cm^3/s
            lossrate(i,j,k,mcf_trnum) = koh(i,j,k)*oh(i,j,k)*temptracer(i,j,k)

            !Calculate MCF rate of change due to photolysis
            lossrate(i,j,k,mcf_trnum) = lossrate(i,j,k,mcf_trnum) + &
                                jrate(i,j,k)*temptracer(i,j,k)

            !Convert MCF rate to kg/m^3/s
            lossrate(i,j,k,mcf_trnum) = lossrate(i,j,k,mcf_trnum)/  &
                     avogadro*mass_mcf*1000.

            !Calculate loss for entire timestep (kg/m^3/timestep)
            lossrate(i,j,k,mcf_trnum) = lossrate(i,j,k,mcf_trnum)*deltat

            !Convert loss rate back to mixing ratio kg [CH4]/kg [AIR]/timestep
            lossrate(i,j,k,mcf_trnum) = lossrate(i,j,k,mcf_trnum)/ &
                        density(i,j,k)

            !Subtract atmospheric losses from tracer (kg/kg)
            tracer(i,j,k,mcf_trnum) = tracer(i,j,k,mcf_trnum) - &
                       lossrate(i,j,k,mcf_trnum) 

            !Calculate MCF rate of change due to ocean deposition
            if (k .eq. 1) then
               tracer(i,j,k,mcf_trnum) = tracer(i,j,k,mcf_trnum) *  &
                  (exp(-1.*oceanloss(i,j,k)/deltaz*deltat))
           endif

      enddo  !k=1,tr_levels

     !kdcorbin, 08/10 - set top level to next lower level
     tracer(i,j,tr_levels,mcf_trnum) = tracer(i,j,tr_levels-1,mcf_trnum)

   enddo !j=1,rows
enddo !i=1,row_length

endif

RETURN

END SUBROUTINE TRACER_MCFLOSS
