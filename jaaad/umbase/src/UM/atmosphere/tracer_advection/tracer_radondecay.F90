!!Routine to calculate radon decay
!
!kdcorbin, 04/10

SUBROUTINE TRACER_RADONDECAY(     &
     &      row_length, rows, tr_levels, tr_vars,       &
     &      radon_trnum,offx,offy,dtime,tracer )

     IMPLICIT NONE

!--------------------------------
!INPUT variables
!--------------------------------

    INTEGER :: row_length, rows, tr_levels, tr_vars
    INTEGER :: radon_trnum, offx, offy
    REAL :: dtime

    REAL                                                 &
   &  tracer(1-offx:row_length+offx,1-offy:rows+offy,    &
   &           tr_levels,tr_vars)                       


!----------------------------
!LOCAL variables
!----------------------------

  INTEGER :: i,j,k

!-----------------------------
!Calculate radon decay
!-----------------------------

if (radon_trnum .gt. tr_vars) then
   write(6,*) 'Radon not included in run.'
   write(6,*) 'Num tracers: ',tr_vars,' radon tracer number: ',radon_trnum
else
   do i=1,row_length
     do j=1,rows
        do k=1,tr_levels
           tracer(i,j,k,radon_trnum) = &
               exp(-dtime*2.11e-6)*tracer(i,j,k,radon_trnum)
        enddo
      enddo
   enddo
endif

RETURN

END SUBROUTINE TRACER_RADONDECAY