#if defined (A14_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE ENG_MASS_DIAG------------------------------------------
!LL
!LL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!LL            - TO GLOBALLY INTERGATE TOTAL ENERGY AMD MASS OF
!LL              THE ATMOSPHERE
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   3.4   26/05/94  Argument LLINTS added and passed to CALC_RS
!LL                   DEF NOWHBR replaced by LOGICAL LWHITBROM
!LL                                                  S.J.Swarbrick
!LL   5.1    4/02/00  Altered to be consistent with new dynamics
!LL                   and use efficient MPP code. R A Stratton.
!LL                   Corrects mass of atmosphere
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL   5.2   18/09/00  Correct energy corection code. R A Stratton.
!LL   5.3   23/04/01  Correct integrals at polar points and case of
!LL                   wet_model_levels < model levels. R A Stratton.
!     5.4 03/04/02  correction to rho_dry                    A. Malcolm
!  6.2   25/12/05  Variable resolution changes            Yongming Tang
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION :
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE ENG_MASS_DIAG(                                         &
     &                      halo_i, halo_j, off_x, off_y                &
     &,                     global_row_length, proc_row_group           &
     &,                     proc_col_group                              &
     &,                     at_extremity, n_proc, n_procx, n_procy      &
     &,                     neighbour                                   &
     &,                     mype                                        &
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_model_levels               &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     delta_lambda,delta_phi                      &
     &,                     FV_cos_theta_latitude,cos_v_latitude        &
     &,                     cos_theta_longitude,sin_theta_longitude     &
     &,                     theta ,u,v,w, rho_r2,q,qcl,qcf              &
     &,                     exner_theta_levels                          &
     &,                     sum_moist_flux                              &
     &,                     tot_mass_init,tot_m_init                    &
     &,                     Lmass_corr,lqt_corr,lemq_print              &
     &,                     a_energysteps,timestep,tot_energy_final     &
     &,                     tot_mass_final,tot_m_final)

      IMPLICIT NONE

#include "c_r_cp.h"
#include "c_a.h"
#include "c_g.h"
#include "c_lheat.h"
#include "fldtype.h"
!
!----------------------------------------------------------------------
! VECTOR LENGTHS
!----------------------------------------------------------------------
! MPP info
!    halo info
      Integer off_x, off_y, halo_i, halo_j
      Integer                                                           &
     &  global_row_length                                               &
                              ! global row length
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same col
     &, neighbour(4)                                                    &
                       ! Array with the Ids of the four neighbours in
                       ! in the horizontal plane
     &, mype           ! processor number

      Logical                                                           &
     &  at_extremity(4)! Indicates if this processor is at north, south,
                       ! east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! input model dimensional info
      Integer                                                           &
     &   row_length                                                     &
                         ! number of points per row
     &,  rows                                                           &
                         ! number of rows on p grid
     &,  n_rows                                                         &
                         ! number of rows on v grid
     &,  model_domain                                                   &
                         ! global or limited area etc.
     &,  model_levels                                                   &
                         ! number of model levels
     &,  wet_model_levels    ! number of wet model levels

! Input radius of model levels
      REAL                                                              &
     &  R_THETA_LEVELS(1-halo_i:row_length+halo_i,                      &
     &            1-halo_j:rows+halo_j,0:model_levels)                  &
     &, R_RHO_LEVELS(1-halo_i:row_length+halo_i,                        &
     &            1-halo_j:rows+halo_j,model_levels)


! Input model grid spacing in radians

      REAL delta_lambda,delta_phi

! Input trigonometric functions
      Real                                                              &
     &  FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                           1-off_y:rows+off_y)                    &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                    1-off_y:n_rows+off_y)                         &
     &, sin_v_latitude (1-off_x:row_length+off_x,                       &
     &                    1-off_y:n_rows+off_y)                         &
     &, cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

!Input model data
      REAL                                                              &
     &  theta(1-off_x:row_length+off_x,                                 &
     &       1-off_y:rows+off_y, model_levels)                          &
                                                !IN theta
     &, U(1-off_x:row_length+off_x,                                     &
     &       1-off_y:rows+off_y, model_levels)                          &
                                                !IN COMPONENT OF WIND
     &, V(1-off_x:row_length+off_x,                                     &
     &       1-off_y:n_rows+off_y, model_levels)                        &
                                                 !IN COMPONENT OF WIND
     &, w(1-off_x:row_length+off_x,                                     &
     &       1-off_y:rows+off_y,0:model_levels)                         &
                                                 !IN w
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_model_levels)                                           &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)
! IN - sum of moisture fluxes into atmosphere over energy period
      Real                                                              &
     &  sum_moist_flux(row_length,rows)

! In/OUT  model data  Altered in Lmass_corr set
      REAL                                                              &
     &  rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &        model_levels)

      Logical                                                           &
     &  Lmass_corr                                                      &
                        ! if true do a mass correction
     &, lqt_corr                                                        &
                        ! if true do a moisture correction
     &, lemq_print      ! if true print addtional info

!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------

      REAL                                                              &
     &     tot_energy_final                                             &
                             ! TOTAL ENERGY OF ATMOSPHERE
     &,    tot_mass_final                                               &
                             ! total dry mass of atmosphere
     &,    tot_m_final       ! total moisture

!IN only
      INTEGER                                                           &
     &     a_energysteps      ! number of timesteps between energy cor
      REAL                                                              &
     &     tot_mass_init                                                &
                              ! total dry mass of atmosphere at T+0
     &,    tot_m_init                                                   &
                              ! total moisture at previous time
     &,    timestep           ! tiemstep of model in seconds.

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

      INTEGER N_SUMS
      PARAMETER (N_SUMS=10)   ! number of global summations
!
      REAL                                                              &
     & sum_array(row_length,rows,n_sums)                                &
                                            ! array to be summed
     &,dummy_array(row_length,rows,3)                                   &
                                            ! dummy array
     &,SUM_RESULTS(n_sums)                                              &
                                            ! sum of SUM_ARRAY
     &,moist_in(1)                     ! Array to interface with do_sums

      Real                                                              &
     &  cv                                                              &
                                       ! specific heat at constant vol
     &, tot_dry_energy                                                  &
     &, tot_wet_mass                                                    &
     &, tot_dry_mass                                                    &
     &, tot_ke                                                          &
     &, err_mass                                                        &
                                       ! error in dry mass
     &, mass_corr                                                       &
                                       ! correction factor dry mass
     &, factor                                                          &
                                       ! grid resolution factor
     &, tot_moist_in                                                    &
     &, err_moist                                                       &
     &, moist_corr                                                      &
                                       ! moist correction
     &, change_moist                                                    &
                                       ! change in qt over period
     &, per_err                                                         &
     &, q_cld                                                           &
                                       ! total cloud water
     &, weight1, weight2 , weight3                                      &
                                       ! weights
     &, tempw, tempd

      Real                                                              &
     &  rho_dry(row_length,rows,model_levels)                           &
                                                    ! rho dry
     &, dry_to_wet(row_length,rows,model_levels)

      INTEGER I, J, K                ! LOOP COUNTER

! pointers for sum_array
      Integer                                                           &
     & ip_dry_mass                                                      &
                        ! dry mass
     &,ip_wet_mass                                                      &
                        ! wet mass
     &,ip_cvT                                                           &
                        ! cvT
     &,ip_gr                                                            &
                        ! gr
     &,ip_keu                                                           &
                        ! keu
     &,ip_kev                                                           &
                        ! kev
     &,ip_kew                                                           &
                        ! kew
     &,ip_q                                                             &
                        ! q
     &,ip_qcl                                                           &
                        ! qcl
     &,ip_qcf           ! qcf

      Parameter (ip_dry_mass=1, ip_wet_mass=2, ip_cvT=3, ip_gr=4,       &
     &           ip_keu=5, ip_kev=6, ip_kew=7, ip_q=8, ip_qcl=9,        &
     &           ip_qcf=10)
      logical lqflux         ! switch for calculations
      INTEGER ISTAT
      REAL SUM_ROWS(rows)
!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS  -
!----------------------------------------------------------------------
      External                                                          &
     & vert_eng_massq, do_sums


!----------------------------------------------------------------------
! zero summation arrays
!----------------------------------------------------------------------
      Do k=1,n_sums
        sum_results(k)=0.0
      Enddo

! set logical so unnecessary calculation not done

      Lqflux = .false.

!----------------------------------------------------------------------
!   call vert_eng_massq to do vertical integrations
!----------------------------------------------------------------------
! DEPENDS ON: vert_eng_massq
      call vert_eng_massq(                                              &
     &                      halo_i, halo_j, off_x, off_y                &
     &,                     global_row_length, proc_row_group           &
     &,                     at_extremity, n_proc, n_procx, n_procy      &
     &,                     neighbour                                   &
     &,                     mype                                        &
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_model_levels               &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     cos_theta_longitude,sin_theta_longitude     &
     &,                     theta ,u,v,w, rho_r2, q,qcl,qcf             &
     &,                     exner_theta_levels                          &
     &,                     Lqflux                                      &
     &,                     rho_dry,dry_to_wet                          &
     &,                     sum_array, dummy_array )

!----------------------------------------------------------------------
!  multiply all vertical integrals by cos(latitude)
!----------------------------------------------------------------------

      Do k = 1,n_sums
        do j = 1,rows
          do i = 1,row_length
             sum_array(i,j,k) = sum_array(i,j,k)*                       &
     &                                FV_cos_theta_latitude(i,j)
          enddo
        enddo
      enddo


!----------------------------------------------------------------------
! reset sum_array so only one contribution from each pole
! NOTE : FV_cos_theta_latitude at poles holds area for whole pole
!----------------------------------------------------------------------
      If (at_extremity(PSouth)) then
        If (at_extremity(PWest)) then
          do k=1,n_sums
            do i=2,row_length
              sum_array(i,1,k)=0.0
            end do
          end do
        else
          do k=1,n_sums
            do i=1,row_length
              sum_array(i,1,k)=0.0
            end do
          end do
        endif
      endif
      If (at_extremity(PNorth)) then
        If (at_extremity(PWest)) then
          do k=1,n_sums
            do i=2,row_length
              sum_array(i,rows,k)=0.0
            end do
          end do
        else
          do k=1,n_sums
            do i=1,row_length
              sum_array(i,rows,k)=0.0
            end do
          end do
        endif
      endif

!----------------------------------------------------------------------
! Section 3 - Calculation of energy
! Work out the global sums and multiply by correct factor for grid res.
!----------------------------------------------------------------------
      factor=delta_lambda*delta_phi   ! already in radians

! Do the global sums
! Note no halo on input array, grid_type 1, halo_type 3

! DEPENDS ON: do_sums
      call do_sums(sum_array,row_length,rows,0,0,1,3,n_sums,            &
     &             sum_results)

! total dry air energy minus latent heat  held by cloud

      tot_dry_energy =( sum_results(ip_cvT)+ sum_results(ip_gr)         &
     &           + sum_results(ip_keu) + sum_results(ip_kev)            &
     &           + sum_results(ip_kew) )*factor

      tot_energy_final  = tot_dry_energy -                              &
     &( lc*sum_results(ip_qcl)+(lc+lf)*sum_results(ip_qcf) )*factor

      tot_dry_mass =sum_results(ip_dry_mass)*factor
      tot_mass_final = tot_dry_mass
      tot_wet_mass =sum_results(ip_wet_mass)*factor

! Print out info if required

      if (mype == 0.and.lemq_print) then
         tot_ke= (sum_results(ip_keu) + sum_results(ip_kev)             &
     &           + sum_results(ip_kew) )*factor


        write(6,*)'Tot dry mass       ',tot_dry_mass
        write(6,*)'Tot mass           ',tot_wet_mass
        write(6,*)'Tot energy         ',tot_energy_final
        write(6,*)'tot dry energy     ',tot_dry_energy

        write(6,*)'gr( rho cal)       ',sum_results(ip_gr)*factor
        write(6,*)'KE( rho cal)       ',tot_ke
        write(6,*)'cvT( rho cal)      ',sum_results(ip_cvt)*factor
        write(6,*)'lq ( rho cal)      ',sum_results(ip_q)*factor*lc
        write(6,*)'lqcf( rho cal)     ',sum_results(ip_qcl)*factor*lc
        write(6,*)'lqcl( rho cal)     ',sum_results(ip_qcf)*factor      &
     &                                          *(lc+lf)
      endif

!----------------------------------------------------------------------
! Section 4. Mass conservation & correction
! Note this is the best place to do this as rho_dry has already been
! evaluated earlier in subroutine
!----------------------------------------------------------------------
      If (Lemq_print.or.Lmass_corr) then
!----------------------------------------------------------------------
! Work out drift in dry mass and printout results
!----------------------------------------------------------------------

        err_mass  = tot_mass_final - tot_mass_init
        mass_corr = tot_mass_init/tot_mass_final

        if (mype == 0) then
          Write(6,*) 'Final dry mass of atmosphere   = ',               &
     & tot_mass_final, ' KG'
          Write(6,*) 'Initial dry mass of atmosphere = ',               &
     & tot_mass_init, ' KG'
          Write(6,*) 'Correction factor for rho_dry  = ',               &
     & mass_corr
        Endif
      Endif

!----------------------------------------------------------------------
! Apply dry mass correction
!----------------------------------------------------------------------
      if (lmass_corr) then

        DO K=1,model_levels
        Do j = 1, rows
          DO I=1,row_length

! adjust rho_dry in case required for q correction
              rho_dry(i,j,k) = rho_dry(i,j,k)*mass_corr
! use new rho dry to recalculate rho
              rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)

          END DO
        END DO
      END DO

      endif        ! lmass_corr

!----------------------------------------------------------------------
! Section 5. Moisture conservation
!----------------------------------------------------------------------
! total moisture at this time

        tot_m_final = (sum_results(ip_q)+sum_results(ip_qcl)            &
     &                +sum_results(ip_qcf)) *factor

      if (lqt_corr.or.lemq_print) then

! total moisture into atmosphere over energy correction period

        tot_moist_in=0.0
        moist_in(1)=0.0
! DEPENDS ON: do_sums
        call do_sums(sum_moist_flux,row_length,rows,0,0,1,3,1,          &
     &             moist_in)
        tot_moist_in = moist_in(1)*factor

!----------------------------------------------------------------------
! Work out drift in moisture
!----------------------------------------------------------------------

       change_moist = tot_m_final - tot_m_init
       err_moist = change_moist - tot_moist_in
       per_err   = 100.*err_moist/change_moist

        if (mype == 0) then
          Write(6,*) 'Final   moisture               = ',               &
     & tot_m_final, ' KG'
          Write(6,*) 'Initial moisture               = ',               &
     & tot_m_init, ' KG'
          Write(6,*) 'change in moisture             = ',               &
     & change_moist, ' KG'
          Write(6,*) 'Moisture added E-P in period   = ',               &
     & tot_moist_in, ' KG'
          Write(6,*) 'Error in moisture              = ',               &
     & err_moist, ' KG'
          Write(6,*) 'Error as % of change           = ',               &
     & per_err

        write(6,*)'q ( rho cal)      ',sum_results(ip_q)*factor
        write(6,*)'qcf( rho cal)     ',sum_results(ip_qcl)*factor
        write(6,*)'qcl( rho cal)     ',sum_results(ip_qcf)*factor

        endif


!----------------------------------------------------------------------
! Apply moisture correction - highly experimental and not simple.
!----------------------------------------------------------------------
! Not tested as basic conservation of moist very poor when running
! using UM5.0
!----------------------------------------------------------------------
! Note problems as rho and q are related, therefore must alter
! both fields
! At present all error in qt applied to q (no correction of qcl or qcf)
!
!   rho = rho_dry / (1-q-qcl-qcf)
!
!  total global qt = sum qt*rho_dry*r^2 dr dlambda dphi cos(phi)
!                    -------------------------------------------
!                       (1-q-qcl-qcf)
!
!Let (qnew+qcl+qcf)/(1-qnew-qcl-qcf)= f (qold+qcl+qcf)/(1-qold-qcl-qcf)
! then
!     qnew =[ A/(1+A) -qcl -qcf ]
!  where
!     A= f ((qold+qcl+qcf)/(1-qold-qcl-qcf))
!
! where f = 1 - (error in qt)/(total qt old)
!
!  rho new = rho_dry / (1 - qnew-qcl-qcf)
!----------------------------------------------------------------------

        if (lqt_corr) then   ! correct water vapour and wet_rho

! correction factor f
          moist_corr = 1. - err_moist/tot_m_init

! Firstly correct q (on q levels or should this be q_rho ?)

          DO K=1,wet_model_levels
            Do j = 1, rows
              DO I=1,row_length
                q_cld = qcl(i,j,k)+qcf(i,j,k)
!              tempw holds A
                tempw  = moist_corr * (q(i,j,k) + q_cld)/               &
     &                                (1. - q(i,j,k) - q_cld)
                q(i,j,k) = ( tempw/(1.+tempw)  ) - q_cld

              END DO
            END DO
          END DO
! Work out for this new q the dry to wet conversion this order to the
! calculations should ensure dry mass is conserved.
! Using linear interpolation weights.

          k = 1
            Do j = 1, rows
              Do i = 1, row_length
       dry_to_wet(i,j,k)= 1./(1. - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
              End Do
            End Do

          Do k = 2, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
                weight2 = r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3 = r_theta_levels(i,j,k)                         &
     &                                      - r_theta_levels(i,j,k-1)
            tempd = ( weight2 *                                         &
     &              (1. - q(i,j,k)- qcl(i,j,k) - qcf(i,j,k) ) +         &
     &                weight1 *                                         &
     &              (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1) ) )   &
     &               / weight3
                dry_to_wet(i,j,k)= 1./tempd
              End Do
            End Do
          End Do

! Special case of wet_model_levels < model_levels

          If (wet_model_levels <  model_levels) then
            k=wet_model_levels+1
            Do j = 1, rows
              Do i = 1, row_length
                weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
                weight2 = r_rho_levels(i,j,k)-r_theta_levels(i,j,k-1)
                weight3 = r_theta_levels(i,j,k)                         &
     &                                      - r_theta_levels(i,j,k-1)
                tempd = ( weight2 +                                     &
     &                   weight1 *                                      &
     &      (1. - q(i,j,k-1)- qcl(i,j,k-1) - qcf(i,j,k-1)) ) / weight3
                dry_to_wet(i,j,k)= 1./tempd
              End Do
            End Do

          Endif   ! wet_model_levels < model_levels

! convert rho dry back to rho wet
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)
              End Do
            End Do
          End Do

! Case of wet_model_levels < model_levels

          If (wet_model_levels <  model_levels) then
            k=wet_model_levels+1
            Do j = 1, rows
              Do i = 1, row_length
                rho_r2(i,j,k) = rho_dry(i,j,k)*dry_to_wet(i,j,k)
              End Do
            End Do
          endif

        endif ! lqt_corr

      endif ! lqt_corr or print
!----------------------------------------------------------------------

      RETURN
      END SUBROUTINE ENG_MASS_DIAG
#endif
