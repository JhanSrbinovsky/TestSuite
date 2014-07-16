#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Writes file domain.dat for the SCM diagnostics

      subroutine write_domain_data(row_length,rows,model_levels,        &
     &     z_top_of_model,first_constant_r_rho_level,                   &
     &     orog,eta_theta,eta_rho)
      implicit none

! Description:
      ! Write file domain.dat which will contain information on the
      ! size and dimensions of the model. This information will be
      ! necessary for making sense of some of the information in the
      ! diagnostic output data files.
! Method:
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      integer row_length,rows    ! IN Horizontal domain
      real orog(row_length*rows) ! IN Orography height
      integer model_levels       ! IN The total no. of model levels
      real z_top_of_model        ! IN The height of the "top of the
                                 ! atmosphere"
      integer first_constant_r_rho_level ! IN Lowest rho level that has
                                         ! constant height
      real eta_theta(model_levels+1), & ! IN The etas of the theta and
     &     eta_rho(model_levels)        ! rho levels

      integer i                 ! Counter
      integer, parameter :: unit=10 ! Temporary unit no. to write to

      open (unit=unit,file='domain.dat')

      ! write out the size of the model domain
      write(unit,'(I4,1X,I4,1X,I4,1X,F9.2,1X,I4)')                      &
     &     row_length,rows,model_levels,                                &
     &     z_top_of_model,first_constant_r_rho_level

      ! write out the eta levels
      write(unit,*)(eta_theta(i),i=1,model_levels+1)
      write(unit,*)(eta_rho(i),i=1,model_levels)

      ! write out the orography data so that PV-wave programs can
      ! calculate heights from eta-levels
      write(unit,*)(orog(i),i=1,row_length*rows)

      close(unit)
      return
      END SUBROUTINE write_domain_data
#endif
