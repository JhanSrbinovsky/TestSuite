#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for calculating online photolysis rates
!   using fast-j
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Current Code Owner:       Colin Johsnon/Fiona O'Connor
!                            Oliver Wild
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
        SUBROUTINE UKCA_FASTJ(                                         &
#include "arglndm.h"
        dj,                                                            &
        p_rho_levels, p0,                                              &
        p_theta_levels,                                                &
        t, tstar,                                                      &
        so4_aitken, so4_accum,                                         &
        q, qcl, qcf,                                                   &
        conv_cloud_lwp, conv_cloud_top, conv_cloud_base,               &
        conv_cloud_amount,                                             &
        surf_alb,                                                      &
        nd_o3, um_ozone3d,                                             &
        land_frac_ctile,                                               &
        z_top_of_model,                                                &
        l_moses_ii)


      USE FASTJ_DATA,        ONLY: fastj_set_limits,                  &
                                fastj_allocate_memory,                 &
                                Blocking,kpcx,jpcl,                    &
                                jjpnl,jppj,nsl,sa,p,rz_3d,             &
                                tau,month,iday,nslon,nslat,            &
                                sza,sza_2d,SZAFAC,SZAFAC_2d,u0,        &
                                od,odw,odi,sulphur,fj_ozone,           &
                                t_fastj=>t       
      USE FASTJ_MIE,         ONLY: NWB3,NWB2
      USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
      USE trignometric_mod,  ONLY: true_longitude
      USE dyn_coriolis_mod,  ONLY: f3_at_u

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typlndm.h"
#include "chsunits.h"
#include "csubmodl.h"
#include "ccontrol.h"
#include "aercmp3a.h"
#include "aerprm3a.h"
#include "cphyscon.h"
#include "ctime.h"
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
#include "mxsize3a.h"
#include "swspdl3a.h"
#include "swspcm3a.h"
#include "c_v_m.h"


      REAL, DIMENSION(row_length,rows,model_levels,jppj), INTENT(INOUT)&
                                           :: dj
      REAL, DIMENSION(row_length,rows,model_levels+1), INTENT(IN)      &
                                           :: p_rho_levels
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN)        &
                                           :: p_theta_levels
      REAL, DIMENSION(row_length,rows,model_levels)     :: t
      REAL, DIMENSION(row_length,rows,tr_levels), INTENT(IN)           &
                                           :: so4_aitken
      REAL, DIMENSION(row_length,rows,tr_levels), INTENT(IN)           &
                                           :: so4_accum
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN)          &
                                           :: q, qcl, qcf
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN)          &
                                           :: conv_cloud_amount
      REAL, DIMENSION(row_length,rows,ozone_levels), INTENT(IN)        &
                                           :: um_ozone3d ! MMR
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: p0
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: tstar
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: conv_cloud_lwp
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: surf_alb
      REAL, DIMENSION(land_points), INTENT(IN) :: land_frac_ctile
      INTEGER, DIMENSION(row_length,rows), INTENT(IN)                  &
                                           :: conv_cloud_top
      INTEGER, DIMENSION(row_length,rows), INTENT(IN)                  &
                                           :: conv_cloud_base


      REAL, INTENT(IN) :: z_top_of_model
      INTEGER, INTENT(IN) :: nd_o3
      LOGICAL, INTENT(IN) :: L_moses_ii

! Local variables

      INTEGER                           :: i_humidity
      INTEGER                           :: i,j,k,l,ia
      INTEGER                           :: ni,nj,ml
      INTEGER                           :: nr_light_points

      REAL, PARAMETER                   :: pi180=pi/180.
      REAL, PARAMETER                   :: ADIFC = 0.06
      REAL                              :: d_eff
      REAL                              :: delta_humidity
      REAL                              :: f
      REAL                              :: timej
      REAL, DIMENSION(row_length,rows)  :: land_frac
      REAL, DIMENSION(row_length, rows) :: sin_true_latitude
      REAL, DIMENSION(row_length,rows,                                 &
                      model_levels)     :: qsat
      REAL, DIMENSION(row_length,rows,                                 &
                      model_levels)     :: rh


      LOGICAL, SAVE                         :: first=.true.
      LOGICAL, DIMENSION(row_length, rows)  :: is_day

!     interface fastj routines

      INTERFACE 
        SUBROUTINE FASTJ_INPHOT
        END SUBROUTINE FASTJ_INPHOT

        SUBROUTINE FASTJ_SOLAR2(timej,sinlat,longitude,lcal360)
          REAL,INTENT(IN)                :: timej
          REAL,DIMENSION(:,:),INTENT(IN) :: sinlat,longitude
          LOGICAL,INTENT(IN)             :: lcal360
        END SUBROUTINE FASTJ_SOLAR2

        SUBROUTINE FASTJ_PHOTOJ(zpj,timej)
          REAL,DIMENSION(:,:,:,:),INTENT(OUT)         :: zpj
          REAL,INTENT(IN)                             :: timej
        END SUBROUTINE FASTJ_PHOTOJ
      END INTERFACE
!

      ni = row_length
      nj = rows
      ml = model_levels

      IF (first) THEN
        Blocking%Mode = 1
        CALL FASTJ_SET_LIMITS(row_length,rows,model_levels)
        CALL FASTJ_ALLOCATE_MEMORY

! DEPENDS ON: fastj_inphot
        CALL FASTJ_INPHOT

        first = .false.
      ENDIF

      CALL FASTJ_COMPUTE_PHOTJ_VALUES

!     set variables concerning model time

      timej              = secs_per_stepim(a_im)/3600.
      iday               = i_day_number
      month              = int(real(iday)*12.0/365.0)+1    !  Approximately
      tau                = i_hour*1.+i_minute/60.+i_second/3600.

      sa = surf_alb
      if(Blocking%Mode <= 3)   then          !not in load balancing mode
        t_fastj(:,:,1:ml)  = t
        p                  = p_rho_levels/100.0  ! Pa to mbar
      end if
      
!     Convert ozone mmr to vmr and write to fj_ozone array

      fj_ozone(:,:,1:ml) = um_ozone3d/c_o3
      fj_ozone(:,:,ml+1) = fj_ozone(:,:,ml)

! DEPENDS ON: fastj_solar2
      CALL FASTJ_SOLAR2 (timej,sin_true_latitude,true_longitude,       &
                         lcal360)

      WRITE(6,*) 'so4_aitken: ',shape(so4_aitken),                     &
              maxval(so4_aitken),minval(so4_aitken),minloc(so4_aitken)
      WRITE(6,*) 'so4_accum:  ',shape(so4_accum),                      &
                  maxval(so4_accum),minval(so4_accum)
      WRITE(6,*) 'sulphur:  ',shape(sulphur),                          &
                  maxval(sulphur),minval(sulphur)
      WRITE(6,*) 'sa:       ',shape(sa),                               &
                  maxval(sa),minval(sa)
      WRITE(6,*) 'p:        ',shape(p),                                &
                  maxval(p),minval(p),minloc(p)
      WRITE(6,*) 'fj_ozone: ',shape(fj_ozone),                         &
               maxval(fj_ozone),minval(fj_ozone),minloc(fj_ozone)
      WRITE(6,*) 'odw:      ',shape(odw),                              &
                  maxval(odi),maxloc(odw),minval(odw)
      WRITE(6,*) 'odi:      ',shape(odi),                              &
                  maxval(odi),maxloc(odi),minval(odi)
      WRITE(6,*) 'od:       ',shape(od),                               &
                  maxval(od),maxloc(od),minval(od)
      
      nslon = 1
      SELECT CASE (Blocking%Mode)

       CASE(1)
!       Blocking one row

        do j=1,rows
          SZA    = SZA_2d(:,j)
          SZAFAC = SZAFAC_2d(:,j)
          U0     = COS(sza*pi180)

          nslat  = j
          do i=1,row_length
            nsl(1,i)=nslon+(i-1)
            nsl(2,i)=nslat
          enddo

! DEPENDS ON: fastj_photoj
          CALL FASTJ_PHOTOJ (dj,timej)

        end do

       CASE(2)
!       Blocking whole domain

        nslat = 1
        do j=1,rows
          k = (j-1)*row_length+1
          SZA(k:k+rows-1)    = SZA_2d(:,j)
          SZAFAC(k:k+rows-1) = SZAFAC_2d(:,j)
          do i=1,row_length
            l=i+(j-1)*row_length
            nsl(1,l)=nslon+(i-1)
            nsl(2,l)=nslat+(j-1)
          end do
        enddo
        U0 = COS(sza*pi180)

! DEPENDS ON: fastj_photoj
        CALL FASTJ_PHOTOJ (dj,timej)

       CASE(3)
        l     = 0
        nslat = 1
        nr_light_points = 0.0
        do j=1,rows
          do i=1,row_length
            if( SZAFAC_2d(i,j)  > 0.001d0) then
              l = l+1
              SZA(l)    = SZA_2d(i,j)
              SZAFAC(l) = SZAFAC_2d(i,j)
              nsl(1,l)  = nslon+(i-1)
              nsl(2,l)  = nslat+(j-1)
            end if

            if(l == Blocking%BL_size  .OR.                             &
              (j == rows .and. i == row_length .and. l > 0))   then

              nr_light_points = nr_light_points+l
              kpcx = l
              jjpnl = jpcl*kpcx
              NWB3 = NWB2*kpcx
              u0(1:l) = COS(sza(1:l)*pi180)
! DEPENDS ON: fastj_photoj              
              CALL FASTJ_PHOTOJ (dj,timej)
              l = 0
            end if
          end do
        end do
        l = count(SZAFAC_2d > 0.001d0)

       CASE(4)
        nr_light_points = count(SZAFAC_2d > 0.001d0)
        CALL FASTJ_LOADBALANCE_PHOTOJ
      END SELECT

      WRITE(6,*) 'tau: ',tau
      WRITE(6,*) 'dj:        ',sum(dj),maxval(dj),maxloc(dj)

      RETURN

!     internal subroutines

      CONTAINS

      SUBROUTINE FASTJ_COMPUTE_PHOTJ_VALUES
        IMPLICIT NONE
        INTEGER, PARAMETER                    :: sw_band_aer = 4
        INTEGER                               :: i,j,k,l

        REAL,DIMENSION(1:row_length,1:rows)                :: total_mass 
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: d_mass
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_accu
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_aitk
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels+1) :: pj


        DO k=1,model_levels
          CALL qsat_wat(qsat(1:row_length,1:rows,k),                   &
                         t(1:row_length,1:rows,k),                     &
            p_theta_levels(1:row_length,1:rows,k),theta_field_size)
        ENDDO
        rh=q(1:row_length,1:rows,:)*(1.0-qsat)/                        &
         (qsat*(1.0-q(1:row_length,1:rows,:)))
        rh=MIN(rh,0.99) ! limit humidity to 99%


! Expand land_fraction
        land_frac=0.
        DO l=1,land_points
          j=(land_index(l)-1)/row_length+1
          i=land_index(l)-(j-1)*row_length
          land_frac(i,j)=land_frac_ctile(l)
        END DO

        sin_true_latitude = f3_at_u(1:ni,1:nj) / two_Omega

! rz in cm
        rz_3d(:,:,1)    = r_Theta_levels(1:ni,1:nj,0)*100.
        rz_3d(:,:,2:ml) = r_rho_levels(1:ni,1:nj,2:ml)*100.
        rz_3d(:,:,ml+1) = r_Theta_levels(1:ni,1:nj,ml)*100.

        pj(:,:,1)=p0
        pj(:,:,2:model_levels)=p_rho_levels(:,:,2:model_levels)
        pj(:,:,model_levels+1)=0.

!      We multiply by 4.125 to convert from mass mixing ratio of 
!      sulphur atoms to mass mixing ratio of ammonium sulphate.
!      Use d_mass to convert to mass of ammonium sulphate.

        d_mass=(pj(:,:,1:ml)-pj(:,:,2:ml+1))/g
        sulph_accu=so4_accum*4.125*d_mass
        sulph_aitk=so4_aitken*4.125*d_mass

        odw(:,:,1:wet_levels)    = qcl(:,:,1:wet_levels)
        odi(:,:,1:wet_levels)    = qcf(:,:,1:wet_levels)
        odw(:,:,wet_levels+1:ml) = 0.0
        odi(:,:,wet_levels+1:ml) = 0.0

        total_mass = 0.0
        DO k = 1,model_levels
            WHERE ( k <= conv_cloud_top .and. k >= conv_cloud_base) 
              total_mass = total_mass + d_mass(:,:,k)
            ENDWHERE
        END DO

        IF (l_3d_cca) THEN
          do k = 1,n_cca_lev
            DO j = 1, rows
              DO i = 1, row_length
               IF (conv_cloud_top(i,j) > 0) THEN
                 IF (t(i,j,k) > tm)  THEN
                  odw(i,j,k) = odw(i,j,k) +                            &
                   (conv_cloud_lwp(i,j)*conv_cloud_amount(i,j,k))/     &
                   total_mass(i,j)
                 ELSE
                  odi(i,j,k) = odi(i,j,k) +                            &
                   (conv_cloud_lwp(i,j)*conv_cloud_amount(i,j,k))/     &
                   total_mass(i,j)
                 END IF
                END IF
              END DO
            END DO
          END DO
        ELSE
          DO k=1,model_levels
            DO j = 1, rows
              DO i = 1, row_length
               IF ( k <= conv_cloud_top(i,j)                           &
                          .and. k >= conv_cloud_base(i,j))   THEN
                 if(t(i,j,k)>tm)  THEN
                   odw(i,j,k) = odw(i,j,k)+(conv_cloud_lwp(i,j)*       &
                          conv_cloud_amount(i,j,1))/total_mass(i,j)
                 ELSE
                  odi(i,j,k)= odi(i,j,k)+(conv_cloud_lwp(i,j)*         & 
                        conv_cloud_amount(i,j,1))/total_mass(i,j)
                 END IF
                END IF
              END DO
            END DO
          END DO
        END IF

!       Convert mass mixing ratios to column densities.

        odw = odw*d_mass
        odi = odi*d_mass

!       set effective radii for water drops, different for land&sea
!       set effective diameter for ice crystals
!       Formulae from John Edwards to convert from column densities to
!       optical depths

        d_eff=100.0 ! in microns
        DO k=1,model_levels
          WHERE (land_frac>0.5) 
           odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/6.)
          ELSE WHERE
           odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/12.)
          END WHERE
        END DO
        odi = odi*(-2.189E-3+3.311E3/d_eff+3.611/d_eff**2)

        DO ia=1,n_aerosol_sw
          IF (type_aerosol_sw(ia) == ip_accum_sulphate .OR.            &
                 type_aerosol_sw(ia)==ip_aitken_sulphate) THEN
            IF (i_aerosol_parametrization_sw(ia) ==                    &
                   ip_aerosol_param_moist) THEN ! moist aerosol
              delta_humidity=humidities_sw(2,ia)-humidities_sw(1,ia)
              DO k=1,model_levels
               DO j = 1, rows
                DO i = 1, row_length
                 i_humidity = INT(rh(i,j,k)/delta_humidity)+1
                 f = (rh(i,j,k)-humidities_sw(i_humidity,ia))/         &
                      delta_humidity
                 od(i,j,k)=                                            &
                 (aerosol_absorption_sw(i_humidity,ia,sw_band_aer)+    &
                  aerosol_scattering_sw(i_humidity,ia,sw_band_aer))*f+ &
                 (aerosol_absorption_sw(i_humidity+1,ia,sw_band_aer)+  &
                  aerosol_scattering_sw(i_humidity+1,ia,sw_band_aer))  &
                  *(1-f)
                END DO
               END DO
              ENDDO
            ELSE ! dry aerosol
              DO k=1,model_levels
                DO j = 1, rows
                  DO i = 1, row_length
                    od(i,j,k) =                                        &
                       aerosol_absorption_sw(1,ia,sw_band_aer)+        &
                       aerosol_scattering_sw(1,ia,sw_band_aer)
                  END DO
                END DO
              ENDDO
            ENDIF
            IF (type_aerosol_sw(ia) == ip_accum_sulphate)              &
              sulph_accu=sulph_accu*od
            IF (type_aerosol_sw(ia) == ip_aitken_sulphate)             &
              sulph_aitk=sulph_aitk*od
          END IF
        ENDDO

        sulphur = sulph_aitk + sulph_accu

        RETURN
      END SUBROUTINE FASTJ_COMPUTE_PHOTJ_VALUES

      SUBROUTINE FASTJ_LOADBALANCE_PHOTOJ

        USE mpl, ONLY :     &
          MPL_ADDRESS_KIND, &
          MPL_REAL,         &
          MPL_INFO_NULL

        IMPLICIT NONE

!--  local variables

        INTEGER                                :: i,j,l,k,m,ii,istat
        INTEGER                                :: is,ie,ian,kr
        INTEGER                                :: np_new,np_old,me
        INTEGER                                :: vals_per_gp
        INTEGER                                :: vals_per_gp_f
        INTEGER                                :: vals_per_gp_b
        INTEGER                                :: np_trans
        INTEGER                                :: npes
        INTEGER                                :: ll
        INTEGER                                :: my_comm
        INTEGER,SAVE                           :: win_b
        
        INTEGER(kind=MPL_ADDRESS_KIND)         :: winsize
        INTEGER(kind=MPL_ADDRESS_KIND)         :: disp

        INTEGER,DIMENSION(0:nproc-1)           :: npoints
        INTEGER,DIMENSION(0:nproc-1,0:nproc-1) :: sr_matrix
        
        LOGICAL                                :: sender_PE,receiver_PE
        
        REAL,DIMENSION(ni*nj)                  :: sa_lb
        REAL,DIMENSION(ni*nj)                  :: SZA_lb,SZAFAC_lb
        REAL,DIMENSION(ni*nj,ml)               :: t_lb,sulphur_lb
        REAL,DIMENSION(ni*nj,ml)               :: odw_lb,odi_lb
        REAL,DIMENSION(ni*nj,ml+1)             :: p_lb,ozone_lb,rz_lb
        REAL,DIMENSION(ni*nj,ml,jppj)          :: dj_lb

!       Global array for MPI-2 one sided communication
        REAL,DIMENSION(:,:),allocatable        :: rbg

!CDIR GM_ARRAY(rbg)
  
!--     Intercaces 

        INTERFACE FASTJ_CALC_NEW_NUMBER_OF_POINTS
          SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS (                 &
                np,npes,npoints,sr_matrix,me,np_new,sender_PE)

            INTEGER,INTENT(IN)                       :: np
            INTEGER,INTENT(IN)                       :: npes
            INTEGER,INTENT(OUT),DIMENSION(0:npes-1)  :: npoints
            INTEGER,INTENT(OUT),                                       &
                DIMENSION(0:npes-1,0:npes-1)         :: sr_matrix
            INTEGER,INTENT(OUT)                      :: me
            INTEGER,INTENT(OUT)                      :: np_new
            LOGICAL,INTENT(OUT)                      :: sender_PE
          END SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS
        END INTERFACE FASTJ_CALC_NEW_NUMBER_OF_POINTS


        npes = nproc

!       calculate balanced number of gridpoints

! DEPENDS ON: fastj_calc_new_number_of_points
        CALL FASTJ_CALC_NEW_NUMBER_OF_POINTS (nr_light_points,         &
            nproc, npoints,sr_matrix,me,np_new,sender_PE)
        receiver_PE = .not. sender_PE


        vals_per_gp_f = 4*ml+3*(ml+1)+3    !values per pridpoint forward
        vals_per_gp_b = jppj*ml            !backward
        vals_per_gp   = max(vals_per_gp_f,vals_per_gp_b)

!       Number of gridpoints to transport
        np_trans = maxval(abs(sr_matrix))
        np_trans = np_trans*vals_per_gp


        if(np_trans == 0)   then   !No data transfer in this case
           sender_PE   = .false.
           receiver_PE = .false.
           np_trans    = 1     !dummy windows length
        end if

!       Create RealGlobalBuffer for data to bet transported between PEs

        allocate(rbg(np_trans,0:npes-1))

!       Get communicator to use from GCOM (i.e. same communicator
!       as the rest of the model)
        call GC_Get_Communicator(MY_COMM, istat)

!       Create MPI2 windows for one-sided communication
!       to allow remote PEs to access data in window
        winsize = size(rbg)*8  ! in bytes
        call MPL_Win_create (rbg,winsize,8,MPL_INFO_NULL,              &
                             MY_COMM,win_b, istat)

!       Initial setup
        np_old  = nr_light_points


!       First, copy and compress local arrays in large one-dimensional
!       spatial chunks and eliminate night pts

!       2-D arrays

        l = 0
        do j=1,nj
!CDIR NODEP
          do i=1,ni
            if( SZAFAC_2d(i,j)  > 0.001d0) then
              l = l+1
              sa_lb(l)     = sa(i,j)
              SZA_lb(l)    = SZA_2d(i,j)
              SZAFAC_lb(l) = SZAFAC_2d(i,j)
            end if
          end do
        end do
!       3-D arrays
        do k=1,ml
          l = 0
          do j=1,nj
!CDIR NODEP
            do i=1,ni
              if( SZAFAC_2d(i,j)  > 0.001d0) then
                l = l+1
                t_lb(l,k)       = t(i,j,k)
                odw_lb(l,k)     = odw(i,j,k)
                odi_lb(l,k)     = odi(i,j,k)
                sulphur_lb(l,k) = sulphur(i,j,k)
              end if
            end do
          end do
        end do
        do k=1,ml+1
          l = 0
          do j=1,nj
!CDIR NODEP
            do i=1,ni
              if( SZAFAC_2d(i,j)  > 0.001d0) then
                l = l+1
                p_lb(l,k)       = p_rho_levels(i,j,k)/100.0 ! Pa to mbar
                ozone_lb(l,k)   = fj_ozone(i,j,k)
                rz_lb(l,k)      = rz_3d(i,j,k)
              end if
            end do
          end do
        end do

!       synchronise data transfer into window (like barrier)
        call MPL_Win_Fence (0,win_b, istat)

!       Load balancing
!       move data from sender PEs to receiver PEs to establish equal 
!       number of columns on every PE.

        if(sender_PE)   then
!         Move data in one-dimensional buffer
!         All data is move between PEs in ONE MPI call/remote PE

          do i=0,npes-1

            if(sr_matrix(me,i) /= 0) then   !something to do
!             is and ie are start and end points of 1D array
              is  = np_old+sr_matrix(me,i)+1!+ because value is negative
              ie  = np_old
              ian = ie-is+1
   
              kr  = 1

!             Fill buffer
!             2-D arrays
              rbg(kr:kr+ian-1,i) = sa_lb(is:ie) ;          kr = kr+ian
              rbg(kr:kr+ian-1,i) = SZA_lb(is:ie) ;         kr = kr+ian
              rbg(kr:kr+ian-1,i) = SZAFAC_lb(is:ie) ;      kr = kr+ian
!             3-D arrays
              do k=1,ml
                rbg(kr:kr+ian-1,i) = t_lb(is:ie,k) ;       kr = kr+ian
                rbg(kr:kr+ian-1,i) = odw_lb(is:ie,k) ;     kr = kr+ian
                rbg(kr:kr+ian-1,i) = odi_lb(is:ie,k) ;     kr = kr+ian
                rbg(kr:kr+ian-1,i) = sulphur_lb(is:ie,k) ; kr = kr+ian
              end do
              do k=1,ml+1
                rbg(kr:kr+ian-1,i) = p_lb(is:ie,k) ;       kr = kr+ian
                rbg(kr:kr+ian-1,i) = ozone_lb(is:ie,k) ;   kr = kr+ian
                rbg(kr:kr+ian-1,i) = rz_lb(is:ie,k) ;      kr = kr+ian
              end do
              np_old = np_old-ian
            end if
          end do
        end if

        call MPL_Win_Fence (0,win_b, istat)

        if(receiver_PE)   then
          do i=0,npes-1  ! looping over remote/sender PEs
            if(sr_matrix(me,i) /= 0) then   !something to do
              is     = np_old+1
              ie     = np_old+sr_matrix(me,i)
              ian    = ie-is+1
              kr     = 1

              ll   = sr_matrix(me,i)*vals_per_gp_f
              ! me is receiver PE
              disp = me*np_trans ! displacement unit (not in bytes)
                                 ! offset in window 
              CALL MPL_Get(rbg(1,i),ll,MPL_REAL, i,disp,ll,            &
                                       MPL_REAL, win_b,istat)

              
!             get data from  buffer
!             2-D arrays
              sa_lb(is:ie)     = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
              SZA_lb(is:ie)    = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
              SZAFAC_lb(is:ie) = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
!             3-D arrays
              do k=1,ml
                t_lb(is:ie,k)       = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                odw_lb(is:ie,k)     = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                odi_lb(is:ie,k)     = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                sulphur_lb(is:ie,k) = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
              end do
              do k=1,ml+1
                p_lb(is:ie,k)       = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                ozone_lb(is:ie,k)   = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                rz_lb(is:ie,k)      = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
              end do
              np_old = np_old+ian
            end if
          end do
        end if

        call MPL_Win_Fence (0,win_b, istat)

!       Call Photolysis in Block mode
!       data for an equal number of grid points is located in the
!       .._lb arrays. The data is coped in chunks of Blocking%BL_size
!       into aray locted in mo_fastj_data. 
!       photoj is called for thes chunks and the output is copied into
!       dj_lb.

        do m=1,np_new,Blocking%BL_size
           do ii=m,min(np_new,m+Blocking%BL_size-1)
             i = mod(ii-1,row_length)+1
             j = (ii-1)/row_length+1
             l = ii-m+1
             sa(i,j)   = sa_lb(ii)
             SZA(l)    = SZA_lb(ii)
             SZAFAC(l) = SZAFAC_lb(ii)
             nsl(1,l)  = i
             nsl(2,l)  = j
           end do
           do k=1,ml
             do ii=m,min(np_new,m+Blocking%BL_size-1)
               i = mod(ii-1,row_length)+1
               j = (ii-1)/row_length+1
               t_fastj(i,j,k)  = t_lb(ii,k)
               odw(i,j,k)      = odw_lb(ii,k)
               odi(i,j,k)      = odi_lb(ii,k)
               sulphur(i,j,k)  = sulphur_lb(ii,k)
             end do
           end do
           do k=1,ml+1
             do ii=m,min(np_new,m+Blocking%BL_size-1)
               i = mod(ii-1,row_length)+1
               j = (ii-1)/row_length+1
               p(i,j,k)        = p_lb(ii,k)
               fj_ozone(i,j,k) = ozone_lb(ii,k)
               rz_3d(i,j,k)    = rz_lb(ii,k)
             end do
           end do
           u0 = COS(sza*pi180)

           kpcx = l
           jjpnl = jpcl*kpcx
           NWB3 = NWB2*kpcx
! DEPENDS ON: fastj_photoj
           CALL FASTJ_PHOTOJ (dj,timej)

           DO l=1,jppj
             DO k=1,ml
               DO ii=m,min(np_new,m+Blocking%BL_size-1)
                 i = mod(ii-1,row_length)+1
                 j = (ii-1)/row_length+1
                 dj_lb(ii,k,l) = dj(i,j,k,l)
               END DO
             END DO
           END DO
        END DO

        CALL MPL_Win_Fence (0,win_b, istat)

        np_old  = nr_light_points

!       Move results back to origional PEs.

        IF (receiver_PE)   THEN

!         Move data in one-dimensional buffer

          DO i=0,npes-1

            IF(sr_matrix(me,i) /= 0) THEN   !something to do
              is     = np_old+1
              ie     = np_old+sr_matrix(me,i)
              ian    = ie-is+1

              kr  = 1

              DO m=1,jppj
                DO k=1,ml
                  rbg(kr:kr+ian-1,i) = dj_lb(is:ie,k,m) ;   kr = kr+ian
                END DO
              END DO

              np_old = np_old+ian
            END IF
          END DO
        END IF

        CALL MPL_Win_Fence (0,win_b, istat)

        IF (sender_PE)   THEN
          DO i=0,npes-1
            IF(sr_matrix(me,i) /= 0) THEN   !something to do
              is  = np_old+sr_matrix(me,i)+1!+ because value is negative
              ie  = np_old
              ian = ie-is+1

              kr     = 1

              ll   = abs(sr_matrix(me,i))*vals_per_gp_b
              disp = me*np_trans
              CALL MPL_Get(rbg(1,i),ll,MPL_REAL, i,disp,ll,            &
                                       MPL_REAL,  win_b,istat)

              DO m=1,jppj
                DO k=1,ml
                  dj_lb(is:ie,k,m) = rbg(kr:kr+ian-1,i) ;   kr = kr+ian
                END DO
              END DO
              np_old = np_old-ian
            END IF
          END DO
        END IF

!       Move data back from the 1-dimensional (spatial) buffer array
!       into the original output 2-D array.

        dj = 0.0
        DO m=1,jppj
          DO k=1,ml
            l = 0
            DO j=1,nj
              DO i=1,ni
                IF( SZAFAC_2d(i,j)  > 0.001d0) THEN
                  l = l+1
                  dj(i,j,k,m) = dj_lb(l,k,m)
                END IF
              END DO
            END DO
          END DO
        END DO
              
        CALL MPL_Win_Fence (0,win_b, istat)

!       Free MPI-2 recources

        CALL MPL_Win_free(win_b,istat)
        DEALLOCATE (rbg)

        RETURN
      END SUBROUTINE FASTJ_LOADBALANCE_PHOTOJ

      END SUBROUTINE UKCA_FASTJ
#endif
