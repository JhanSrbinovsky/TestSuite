#if defined(C92_2A)
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
#include "domtyp.h"


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

#include "c_pi.h"
#include "parparm.h"
!
! Workspace usage
!
      Real aa(row_length, rows),bb(row_length, rows),                   &
     &     pre(row_length, rows),fac

#if defined(T3E)
      Integer add,rem_add
      common/addresses/ add,rem_add

      Real pe0_array(row_length,rows)
      pointer(ptr,pe0_array)
#else
      INTEGER ICODE           ! return code from GATHER_FIELD
      CHARACTER*80 CMESSAGE   ! Error message from GATHER_FIELD
#endif
!
!   Constants
!
      n=15
      Kappa=R/Cp
      rho=cos(pi/n)
      niter=20

!
#if defined(T3E)
!  collect from local arrays to global arrays on PE0
!
      add=loc(pstar)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()


      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  pstar_local(1,j),row_length_local,0)
      enddo

      call barrier()

      add=loc(pmsl)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  pmsl_local(1,j),row_length_local,0)
      enddo

      call barrier()

      add=loc(thetas)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  thetas_local(1,j),row_length_local,0)
      enddo

      call barrier()

      add=loc(tm)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  tm_local(1,j),row_length_local,0)
      enddo

      call barrier()

      add=loc(phi_star)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  phi_star_local(1,j),row_length_local,0)
      enddo

      call barrier()

      add=loc(cos_p_latitude)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_put(                                                   &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  cos_p_latitude_local(1,j),row_length_local,0)
      enddo
      call barrier()
#else

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

#endif

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
#if defined(T3E)
      call barrier()
      add=loc(exnert)
      call barrier()
      call shmem_get(rem_add,add,1,0)
      ptr=rem_add
      call barrier()

      do j=1,rows_local
      call shmem_get(                                                   &
     &  exnert_local(1,j),                                              &
     &  pe0_array(l_datastart(1),l_datastart(2)+(j-1)),                 &
     &  row_length_local,0)
      enddo
      call barrier()
#else
! DEPENDS ON: scatter_field
      CALL SCATTER_FIELD(exnert_local,exnert,                           &
     &                  row_length_local,rows_local,                    &
     &                  row_length,rows,                                &
     &                  fld_type_p,halo_type_no_halo,                   &
     &                  0,all_proc_group,                               &
     &                  ICODE,CMESSAGE)
#endif

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

#endif
