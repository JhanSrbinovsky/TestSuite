!!Routine to calculate CH4 atmospheric loss from OH, O1D and CL 
!
!kdcorbin, 04/10

SUBROUTINE TRACER_METHANELOSS(     &
     &      row_length, rows, model_levels,maxno_ch4tracers    &
     &     ,tr_levels, tr_vars, offx, offy, secs_per_step      &
     &     ,theta,pressure,density                             &
     &     ,tracer,oh,strloss,lossrate )

     IMPLICIT NONE

!--------------------------------
!INPUT variables
!--------------------------------

    INTEGER :: row_length, rows, model_levels
    INTEGER :: maxno_ch4tracers,tr_levels, tr_vars
    INTEGER :: offx, offy
    REAL :: secs_per_step

    REAL                                                 &
   &  theta(1-offx:row_length+offx,1-offy:rows+offy,     &
   &          model_levels)                              &
   &, pressure(1-offx:row_length+offx,1-offy:rows+offy,  &
   &          model_levels)                              &
   &, density(1-offx:row_length+offx,1-offy:rows+offy,   &
   &          model_levels)                              &
   &, tracer(1-offx:row_length+offx,1-offy:rows+offy,    &
   &           tr_levels,tr_vars)                        &
   &, oh(1:row_length,1:rows,1:model_levels)             &
   &, strloss(1:row_length,1:rows,1:model_levels)        &
   &, lossrate(row_length,rows,tr_levels,tr_vars)


!----------------------------
!LOCAL variables
!----------------------------

  REAL, PARAMETER :: avogadro=6.022E23
  REAL, PARAMETER :: pref=100000.,rair=287.05,cp=1005.
  REAL, PARAMETER :: kappa=rair/cp

  REAL :: temperature(row_length,rows,model_levels)
  REAL :: koh(row_length,rows,model_levels)
  REAL :: temptracer(row_length,rows,tr_levels)
  REAL :: ohlossrate(row_length,rows,tr_levels,tr_vars)
  REAL :: strlossrate(row_length,rows,tr_levels,tr_vars)

  INTEGER :: i,j,k,ktemp,tr,n_tr

!-----------------------------
!Calculate temperature from theta and pressure
!Calculate OH rate constant
!-----------------------------
do i=1,row_length
  do j=1,rows
     !kdcorbin, 07/10 - set temp in top level to one below
     do k=1,model_levels  !-1
        temperature(i,j,k) =  &
            theta(i,j,k)/((pref/pressure(i,j,k))**kappa)
        koh(i,j,k) = 2.45E-12*exp(-1775/temperature(i,j,k))
     enddo
     !temperature(i,j,model_levels) = temperature(i,j,model_levels-1)
     !koh(i,j,model_levels) = 2.45E-12*exp(-1775/temperature(i,j,model_levels))
   enddo
enddo


!----------------------------
!Calculate ch4 loss from OH
!----------------------------
if (tr_vars .gt. maxno_ch4tracers) then
   n_tr=maxno_ch4tracers
else
   n_tr=tr_vars
endif


!!!Tracers have methane mass mixing ratio, with units of kg [CH4] / kg [AIR]
 do tr=1,n_tr
   do i=1,row_length
      do j=1,rows

         temptracer(i,j,1:tr_levels) = &
                tracer(i,j,1:tr_levels,tr)

         !kdcorbin, 07/10 - do not use loss on top level
         do k=1,tr_levels-1
   
            !Convert CH4 units to kg [CH4] / m^3
            temptracer(i,j,k) = temptracer(i,j,k)*density(i,j,k)

            !Convert CH4 units to molecules/cm^3
            temptracer(i,j,k) = temptracer(i,j,k)*avogadro/16./1000.

            !Calculate CH4 rate of change due to OH in molecules/cm^3/s
            ohlossrate(i,j,k,tr) = koh(i,j,k)*temptracer(i,j,k)*oh(i,j,k)
            strlossrate(i,j,k,tr) = strloss(i,j,k)*temptracer(i,j,k)
            lossrate(i,j,k,tr) = ohlossrate(i,j,k,tr) + &
                                 strlossrate(i,j,k,tr)

            !Convert CH4 rate to kg/m^3/s
            lossrate(i,j,k,tr) = lossrate(i,j,k,tr)/avogadro*16.*1000.

            !Calculate loss for entire timestep (kg/m^3/timestep)
            lossrate(i,j,k,tr) = lossrate(i,j,k,tr)*secs_per_step

            !Convert loss rate back to mixing ratio kg [CH4]/kg [AIR]/timestep
            lossrate(i,j,k,tr) = lossrate(i,j,k,tr)/density(i,j,k)

            !Subtract loss from tracer (kg/kg)
            tracer(i,j,k,tr) = tracer(i,j,k,tr)-lossrate(i,j,k,tr) 

          enddo  !k=1,tr_levels

          !kdcorbin, 08/10 - set top level to next lower level
          tracer(i,j,tr_levels,tr) = tracer(i,j,tr_levels-1,tr)

      enddo  !j=1,rows
   enddo  !i=1,row_length
 enddo !tr=1,n_tr

RETURN

END SUBROUTINE TRACER_METHANELOSS
