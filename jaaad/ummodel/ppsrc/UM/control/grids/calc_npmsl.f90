
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_PMSL

      Subroutine Calc_NPMSL(pmsl_local,pstar_local,                     &
     &                 phi_star_local,thetas_local,tm_local,            &
     &                 cos_p_latitude_local,delta_lambda,delta_phi,     &
     &                 row_length_local,rows_local,                     &
     &                 row_length,rows,                                 &
     &                 g, R, Lapse_Rate, earth_radius,                  &
     &                 me, n_proc, n_procx, n_procy,                    &
     &                 off_x, off_y,                                    &
     &                 l_datastart,                                     &
     &                 neighbour, at_extremity,                         &
     &                 all_proc_group,model_domain,                     &
     &                 Cp, npmsl_height)
!
!-------------------------------------------------------------
!

! Purpose:
!          Modifies Pressure at mean sea level as calculated
!          in Calc_PMSL to account for errors due to high
!          orography
!
! Method:
!          As described in section 4.5 of UMDP 80. R. Rawlins
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
!  Version     Date    Comment
!   ----     -------     -------
!    6.3     24/10/06   Calc_NPMSL subroutine split from Calc_PMSL
!                       deck. UM Systems
!    6.4     16/11/06   Modify routine to expect temperature
!                       at MSL instead of Tstar. Paul Earnshaw
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      implicit none
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end


      Integer                                                           &
     &  off_x                                                           &
     &, off_y

      Integer                                                           &
     &  me                                                              &
                            !IN. Processor number
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, l_datastart(3)                                                  &
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, all_proc_group                                                  &
                       ! Group id for all processors
     &, model_domain     ! indicator as to model type, ie global, lam

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Real                                                              &
     & Cp

      Real                                                              &
     &  R                                                               &
     &, g                                                               &
     &, Lapse_Rate                                                      &
     &, earth_radius                                                    &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, npmsl_height ! Orographic height above which relaxation occurs

! --------------------------------------------------
! In the following variable definitions,
!  thetas = theta at Surface
!  thetam = theta at Mean Sea Level
! The same naming convention applies to t and exner
! --------------------------------------------------

!
!  Define variables local to PE
!
      Integer                                                           &
     & row_length_local                                                 &
     &,rows_local

      Real                                                              &
     & pmsl_local (row_length_local, rows_local)                        &
     &,pstar_local(row_length_local, rows_local)                        &
     &,phi_star_local(row_length_local, rows_local)                     &
     &,thetas_local(1-off_x:row_length_local+off_x,                     &
     &                           1-off_y:rows_local+off_y)              &
     &,tm_local(row_length_local, rows_local)                           &
     &,cos_p_latitude_local(row_length_local, rows_local)               &
     &,exnert_local(row_length_local, rows_local)
!
! Define PE0's variables
!
      Integer                                                           &
     & row_length                                                       &
     &,rows

      Real                                                              &
     & pmsl(row_length, rows)                                           &
     &,pstar(row_length, rows)                                          &
     &,phi_star(row_length, rows)                                       &
     &,thetas(row_length, rows)                                         &
     &,tm(row_length, rows)                                             &
     &,cos_p_latitude(row_length, rows)

!
! Define local variables
!
      Real                                                              &
     & exners(row_length, rows)                                         &
     &,exnerm(row_length, rows)                                         &
     &,exnert(row_length, rows,2)                                       &
                                  ! New exner pressure
     &,avexner(row_length, rows)                                        &
     &,orog(row_length, rows)                                           &
     &,thetam(row_length, rows)                                         &
     &,fug(row_length, rows)                                            &
     &,fvg(row_length, rows)                                            &
     &,f(row_length, rows)                                              &
     &,sphlo1(row_length, rows)                                         &
     &,sphlo2(row_length, rows)                                         &
     &,sphla1                                                           &
                ! Latitundinal spherical term
     &,sphla2                                                           &
                ! Latitundinal spherical term squared
     &,niter                                                            &
                        ! Number of iterations
     &,orr                                                              &
                        ! Over relaxation factor
     &,oorr                                                             &
                        ! Over relaxation factor at previous iteration
     &,rho                                                              &
                        ! Spherical radius of convergence
     &,n                ! Typical number of points above npmsl_height

      Real Kappa

      Integer i,j,kk,mm
      Integer ip1,im1,jp1,jm1

!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
!
! Workspace usage
!
      Real aa(row_length, rows),bb(row_length, rows),                   &
     &     pre(row_length, rows),fac

      INTEGER ICODE           ! return code from GATHER_FIELD
      CHARACTER*80 CMESSAGE   ! Error message from GATHER_FIELD
!
!   Constants
!
      n=15
      Kappa=R/Cp
      rho=cos(pi/n)
      niter=20

!


! DEPENDS ON: gather_field
      CALL GATHER_FIELD(pstar_local,pstar,                              &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(pmsl_local,pmsl,                                &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(thetas_local,thetas,                            &
     &                  row_length_local+2*off_x,rows_local+2*off_y,    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_single,                    &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(tm_local,tm,                                    &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(phi_star_local,phi_star,                        &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(cos_p_latitude_local,cos_p_latitude,            &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)



!
!  Calc of new pmsl on PE0
!

      if(me == 0)then

        sphla1=1./(delta_phi*earth_radius)
        sphla2=(delta_phi*earth_radius)**2


      Do j = 1, rows
        Do i = 1, row_length
        exners(i,j)=(pstar(i,j)/100000.)**Kappa
        exnerm(i,j)=(pmsl(i,j)/100000.)**Kappa
        orog(i,j)=phi_star(i,j)/g
        thetam(i,j)=tm(i,j)/exnerm(i,j)
! the check below only applies very close to the poles and the values 
! of the two variables set to arbritrary numbers don't get used in the 
! correction step
        if(abs(cos_p_latitude(i,j)) < 1.0e-8)then
          sphlo1(i,j)=1.0e8
          sphlo2(i,j)=0.0
        else
          sphlo1(i,j)=1./(delta_lambda*earth_radius*cos_p_latitude(i,j))
          sphlo2(i,j)=(delta_lambda*earth_radius*cos_p_latitude(i,j))**2
        end if
        pre(i,j)=(sphlo2(i,j)*sphla2)/((2*sphlo2(i,j))+(2*sphla2))
        exnert(i,j,1)=exnerm(i,j)
        exnert(i,j,2)=exnerm(i,j)
        avexner(i,j)=exnerm(i,j)
        Enddo
      Enddo
!
! Calculate geostrophic wind
!
      Do j = 1, rows
        Do i = 1, row_length
        if (i == 1) then
          ip1=i+1
          if (model_domain  ==  mt_lam) then
            im1 = i
            fac = 1.0
          else
            im1=row_length
            fac = 0.5
          endif

        elseif(i == row_length) then
          if (model_domain  ==  mt_lam) then
            ip1 = i
            fac = 1.0
          else
            ip1=1
            fac = 0.5
          endif
          im1=i-1

        else
          ip1=i+1
          im1=i-1
          fac=0.5
        endif

          fvg(i,j)=fac*sphlo1(i,j)*((phi_star(ip1,j)-phi_star(im1,j))   &
     &          +Cp*thetas(i,j)*(exners(ip1,j)-exners(im1,j)))

        if (j == 1)then
          jp1=j+1
          jm1=j
          fac=1.0
        elseif (j == rows)then
          jp1=j
          jm1=j-1
          fac=1.0
        else
          jp1=j+1
          jm1=j-1
          fac=0.5
        endif

        fug(i,j)=-fac*sphla1*((phi_star(i,jp1)-phi_star(i,jm1))         &
     &          +Cp*thetas(i,j)*(exners(i,jp1)-exners(i,jm1)))

        Enddo
      Enddo
!
! Calculate RHS of equation F
!
      Do j = 1, rows
        Do i = 1, row_length
        aa(i,j)=fvg(i,j)/thetam(i,j)
        bb(i,j)=fug(i,j)/thetam(i,j)
        Enddo
      Enddo

      If (model_domain  ==  mt_lam) Then
          Do j = 2, rows-1
            Do i = 2, row_length-1
              f(i,j)=(0.5/Cp)*sphlo1(i,j)*(aa(i+1,j)-aa(i-1,j))         &
     &           -(0.5/Cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
            Enddo
          Enddo
      Else
        Do j = 2, rows-1
          Do i = 1, row_length
          if (i == 1)then
            ip1=i+1
            im1=row_length
          elseif (i == row_length)then
            ip1=1
            im1=i-1
          else
            ip1=i+1
            im1=i-1
          endif
          f(i,j)=(0.5/Cp)*sphlo1(i,j)*(aa(ip1,j)-aa(im1,j))             &
     &          -(0.5/Cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
          Enddo
        Enddo

      Endif

!
! Relaxation to find exner twiddle
!
      Do kk=1,niter
        if(kk == 1)then
          orr=1.
        elseif(kk == 2)then
          orr=1./(1.-((rho**2)/2.))
        else
         orr=1./(1.-((oorr*rho**2)/4.))
        endif
      If (model_domain  ==  mt_lam) Then
        Do j = 3, rows-2
          Do i = 3, row_length-2
            if(orog(i,j) >  npmsl_height)then
              avexner(i,j)=pre(i,j)*((1./sphla2*(exnert(i+1,j,1)        &
     &            +exnert(i-1,j,2)))                                    &
     &            +(1./sphlo2(i,j)*(exnert(i,j+1,1)                     &
     &            +exnert(i,j-1,2)))-f(i,j))
              exnert(i,j,2)=exnert(i,j,1)+orr*(avexner(i,j)             &
     &            -exnert(i,j,1))
            endif
          Enddo
        Enddo
      Else
        Do j = 3, rows-2
          Do i = 1, row_length
          if(orog(i,j) >  npmsl_height)then
          if(i == 1)then
            ip1=i+1
            im1=row_length
            mm=1
          elseif(i == row_length)then
            ip1=1
            im1=i-1
            mm=2
          else
            ip1=i+1
            im1=i-1
            mm=2
          endif
          avexner(i,j)=pre(i,j)*((1./sphla2*(exnert(ip1,j,1)            &
     &              +exnert(im1,j,mm)))                                 &
     &              +(1./sphlo2(i,j)*(exnert(i,j+1,1)                   &
     &              +exnert(i,j-1,2)))-f(i,j))
          exnert(i,j,2)=exnert(i,j,1)+orr*(avexner(i,j)-exnert(i,j,1))
          endif
          Enddo
        Enddo
      Endif

      Do j = 1, rows
        Do i = 1, row_length
        exnert(i,j,1)=exnert(i,j,2)
        Enddo
      Enddo
        oorr=orr

      ENDDO

      endif
!
!  distribute from exnert on PE0 to exnert_local
!
! DEPENDS ON: scatter_field
      CALL SCATTER_FIELD(exnert_local,exnert,                           &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)

!
! Calculate new pmsl
!

      Do j = 1, rows_local
        Do i = 1, row_length_local
        pmsl_local(i,j)=100000.*(exnert_local(i,j)**(1./Kappa))
        Enddo
      Enddo

! end of routine

      Return
      END SUBROUTINE Calc_NPMSL

