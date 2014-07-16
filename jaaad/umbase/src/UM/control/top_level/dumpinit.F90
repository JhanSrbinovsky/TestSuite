#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!     Subroutine DUMPINIT
!     PURPOSE:-           To initialise primary variables from RESDUMP
!     read in previously from tape
!     PROGRAMMER:-        J. LEAN
!
!
      Subroutine DUMPINIT(                                              &
        row_length, rows, nprimvars, land_points, nlevs, nwet,          &
        nbl_levs, nsoilt_levs, nsoilm_levs, ntrop, n_cca_levels,        &
        land_sea_mask, resdump, u, v, w, t, theta, q, qcl, qcf,         &
        layer_cloud, p, rho, t_deep_soil, smc, canopy_gb, snodep,       &
        tstar, zh, z0msea, cca, rccb, rcct, smcl)

      Implicit none

!     Arguments :
      Integer                                                           &
     &  row_length                                                      &
                                ! IN x dimension
     &  , rows                                                          &
                                ! IN y dimension
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q
                                !  is set.
     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temperature
                                !  levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels
     &  ,ntrop                                                          &
                                ! IN Max number of levels in the
                                !  troposphere
     &  ,nprimvars                                                      &
                                ! IN minimum no. of
                                !  variables required to restart
     &  ,land_points                                                    &
                                ! IN no of land_points            

     &  ,n_cca_levels            ! IN  no of levels of cca


      Integer                                                           &
     &  i,j,k                                                           &
                                ! IN Loop counter
     &  ,icount                 ! IN Counter
      Real                                                              &
     &  resdump(row_length, rows,nprimvars)                             &
                                            ! IN Contains restart dump
     &  ,u(row_length, rows,nlevs)                                      &
                                ! IN Zonal wind at each level
                                !  (m s^-2)
     &  ,v(row_length, rows,nlevs)                                      &
                                ! IN Meridional wind at each level
                                !  (m s^-2)
     &  ,w(row_length, rows,0:nlevs)                                    &
                                     !IN vertical velocity at each
                                ! level
     &  ,p( row_length, rows, nlevs+1)                                  &
     &  ,rho(row_length, rows, nlevs)                                   &
     &  ,t(row_length, rows,nlevs)                                      &
                                ! IN Temperature at each level
     &  ,theta(row_length, rows,nlevs)                                  &
                                ! IN Potential temp. (K)
     &  ,q(row_length, rows,nwet)                                       &
                                ! IN Specific humidity (kg kg^-1)
     &  ,layer_cloud(row_length, rows,nwet)                             &
                                ! IN Layer cloud amount (decima
     &  ,qcl(row_length, rows,nwet)                                     &
                                ! IN Cloud water content (kg kg^-1)
     &  ,qcf(row_length, rows,nwet)                                     &
                                ! IN Cloud ice content (kg kg^-1)
     &  ,pstar(row_length, rows)                                        &
                                ! IN Pressure at earth's surface
                                !  (Pa not HPa)
     &  ,t_deep_soil(land_points, nsoilt_levs)                          &
                                ! IN Deep soil temperatures
     &  ,smc(land_points)                                               &
                                ! IN Soil moisture content
                                !  (kg m^-2)
     &  ,smcl(land_points, nsoilm_levs)                                 &
                                ! IN soil moisture in layers
                                !  (kg m^-2)
     &  ,rccb(row_length, rows)                                         &
                                ! IN Convective cloud base
     &  ,rcct(row_length, rows)                                         &
                                ! IN Convective cloud top
     &  ,canopy_gb(land_points)                                         &
                                ! IN Canopy water content (kg/m2)
     &  ,snodep(row_length, rows)                                       &
                                ! IN Snow depth (kg m^-2)
     &  ,tstar(row_length, rows)                                        &
                                ! IN Surface temperature (K)
     &  ,zh(row_length, rows)                                           &
                                ! IN Height above surface of to
                                ! IN boundary layer (m)
     &  ,z0msea(row_length, rows)                                       &
                                ! IN Sea surface roughness leng
     &  ,cca(row_length, rows, n_cca_levels)
                                ! IN Convective cloud amount

      Logical :: Land_sea_mask(row_length,rows)
                                ! IN True if land point
      ! Local
      Integer :: land_cnt

      
      land_cnt = 0
         
      Do k = 1, rows
        Do j = 1, row_length

        If (land_sea_mask(j,k)) Then
          land_cnt = land_cnt + 1
        End If

        Do i = 1, nlevs
            u(j,k, i) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nlevs-1
            v(j,k, i-icount + 1) = resdump(j,k, i)
          enddo
          icount = i
          Do i = icount, icount + nlevs
            w(j,k, i-icount) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nlevs-1
            t(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nlevs-1
            theta(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nwet-1
            q(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nwet-1
            qcl(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nwet-1
            qcf(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i
        Do i = icount, icount + nwet-1
            layer_cloud(j,k, i-icount + 1) = resdump(j,k, i)
        enddo
        icount = i

        If (land_sea_mask(j,k)) Then
          Do i = icount, icount + nsoilt_levs-1
            t_deep_soil(land_cnt, i-icount + 1) = resdump(j,k, i)
          End Do
          icount = i
        End If

          Do i = icount, icount + nlevs+1
            p(j,k, i-icount + 1) = resdump(j,k, i)
          enddo
          icount = i
          Do i = icount, icount + nlevs
            rho(j,k, i-icount + 1) = resdump(j,k, i)
          enddo
          icount = i

          If (land_sea_mask(j,k)) Then
            smc(land_cnt) = resdump(j,k, icount)
            icount = icount + 1
            canopy_gb(land_cnt) = resdump(j,k, icount)
            icount = icount + 1
          End If

          snodep(j,k) = resdump(j,k, icount)
        icount = icount + 1
          tstar(j,k) = resdump(j,k, icount)
        icount = icount + 1
          zh(j,k) = resdump(j,k, icount)
        icount = icount + 1
          z0msea(j,k) = resdump(j,k, icount)
        icount = icount + 1
          Do i = icount, icount + n_cca_levels
            cca(j,k, i-icount + 1) = resdump(j,k, i)
          enddo
        icount = icount + 1
          rccb(j,k) = resdump(j,k, icount)
        icount = icount + 1
          rcct(j,k) = resdump(j,k, icount)

          If (land_sea_mask(j,k)) Then
            Do i = 1, nsoilm_levs
              smcl(land_cnt, i) = resdump(j,k, icount + i)
            End Do
          End If

        End Do                     ! j
      End Do                     ! k

      Return
      END SUBROUTINE DUMPINIT
#endif
