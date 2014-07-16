#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!     SUBROUTINE RESTART_DUMP
!     PURPOSE:-           To create restart dump for subsequent runs
!
!     PROGRAMMER:-        J. LEAN
!
!     = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
      Subroutine RESTART_DUMP(                                                &
        row_length, rows, nlevs, nwet, nprimvars, land_points,                &
        nbl_levs, nsoilt_levs, nsoilm_levs, n_cca_levels,                     &
        land_sea_mask, resdump, u, v, w, t, theta, q, qcl, qcf,               &
        layer_cloud, p, rho, t_deep_soil, smc, canopy_gb, snodep, tstar,      &
        zh, z0msea, cca, iccb, icct, smcl)

      Implicit none
!
      Integer                                                           &
     &  row_length                                                      &
                                ! IN x dimension of arrays
     &  ,rows                                                           &
                                ! IN y dimension of arrays
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q
                                !   is set.
     &  ,nprimvars                                                      &
                                ! IN minimum no. of
                                !   variables required to restart
                                !   from a dump and is equal to
     &  ,land_points                                                    &
                                ! IN no of land_points            

     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temperature
                                !   levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels
     &,  n_cca_levels      ! IN no of cca levels
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!     Primary model variables + T
!
      Integer                                                           &
     &  iccb(row_length, rows)                                          &
                                ! IN Convective cloud base
     &  ,icct(row_length, rows) ! IN Convective cloud top
      Real                                                              &
     &   canopy_gb(land_points)                                         &
                                ! IN Canopy water content (kg/m2)
     &  ,cca(row_length, rows, n_cca_levels)                            &
                               ! IN Convective cloud amount
     &  ,layer_cloud(row_length, rows,nwet)                             &
                                ! IN Layer cloud amount (decima
     &  ,q(row_length, rows,nwet)                                       &
                                  ! IN Specific humidity (kg kg^-1)
     &  ,qcf(row_length, rows,nwet)                                     &
                                    ! IN Cloud ice content (kg kg^-1)
     &  ,qcl(row_length, rows,nwet)                                     &
                                    ! IN Cloud water content (kg kg^-)
     &  ,smc(land_points)                                               &
                                ! IN Soil moisture content
                                !  (kg m^-2)
     &  ,smcl(land_points, nsoilm_levs)                                 &
                                ! IN Soil moisture in levels
                                !  (kg m^-2)
     &  ,snodep(row_length, rows)                                       &
                                   ! IN Snow depth (kg m^-2)
     &  ,t(row_length, rows,nlevs)                                      &
                                   ! IN Temperature at each level (K)
     &  ,t_deep_soil(land_points, nsoilt_levs)                          &
                                ! IN Deep soil temperatures K
     &  ,theta(row_length, rows,nlevs)                                  &
                                ! IN Potential temperature (K)
     &  ,tstar(row_length, rows)                                        &
                                ! IN Surface temperature (K)
     &  ,rho(row_length, rows, nlevs)                                   &
     &  ,p(row_length, rows,nlevs+1)                                    &
                                     ! IN Pressure (mb)
     &  ,u(row_length, rows,nlevs)                                      &
                                   ! IN Zonal wind at each level
                                !  (m s^-2)
     &  ,v(row_length, rows,nlevs)                                      &
                                   ! IN Meridional wind at each level
                                !  (m s^-2)
     &  ,w(row_length, rows,0:nlevs)                                    &
                                     ! IN vertical vel. at each level
     &  ,zh(row_length, rows)                                           &
                                   ! IN Height above surface of top of
                                !  boundary layer (m)
     &  ,z0msea(row_length, rows)                                       &
                                   ! IN Sea surface roughness length
     &  ,resdump(row_length, rows,nprimvars)
                                   ! IN Contains restart dump

     Logical :: Land_sea_mask(row_length,rows)
                                ! IN True if land point

      Integer                                                           &
        i, j, k                                                         &
                               ! Loop counter
      , icount                 ! Counter

      Integer :: land_cnt      ! land point counter
     

      land_cnt = 0

      Do k = 1, rows
        Do j = 1, row_length

        If (land_sea_mask(j,k)) Then
          land_cnt = land_cnt + 1
        End If

        Do i = 1, nlevs
            resdump(j,k, i) = resdump(j,k,i) + u(j,k,i)
        enddo
        icount = i
        Do i = icount,icount + nlevs-1
            resdump(j,k,i) = resdump(j,k,i) + v(j,k,i-icount + 1)
          enddo
          icount = i
          Do i = icount,icount + nlevs
            resdump(j,k,i) = resdump(j,k,i) + w(j,k,i-icount)
        enddo
        icount = i
        Do i = icount,icount + nlevs-1
            resdump(j,k,i) = resdump(j,k,i) + t(j,k,i-icount + 1)
        enddo
        icount = i
        Do i = icount,icount + nlevs-1
            resdump(j,k,i) = resdump(j,k,i) + theta(j,k,i-icount + 1)
        enddo
        icount = i
        Do i = icount,icount + nwet-1
            resdump(j,k,i) = resdump(j,k,i) + q(j,k,i-icount + 1)
        enddo
        icount = i
        Do i = icount,icount + nwet-1
            resdump(j,k,i) = resdump(j,k,i) + qcl(j,k,i-icount + 1)
        enddo
        icount = i
        Do i = icount,icount + nwet-1
            resdump(j,k,i) = resdump(j,k,i) + qcf(j,k,i-icount + 1)
        enddo
        icount = i
        Do i = icount,icount + nwet-1
            resdump(k, j, i) = resdump(j,k,i) +                         &
     &                layer_cloud(j,k,i-icount + 1)
        enddo
        icount = i

        If (land_sea_mask(j,k)) Then
          Do i = icount,icount + nsoilt_levs-1
              resdump(j,k,i) = resdump(j,k,i)                           &
                             + t_deep_soil(land_cnt,i-icount + 1)
          End Do
          icount = i
        End If

          Do i = icount,icount + nlevs+1
            resdump(j,k,i) = resdump(j,k,i) + p(j,k,i-icount + 1)
          enddo
          icount = i
          Do i = icount, icount + nlevs
            resdump(j,k,i) = resdump(j,k,i) + rho(j,k,i-icount + 1)
          enddo
          icount = i

          If (land_sea_mask(j,k)) Then
            resdump(j,k,icount) = resdump(j,k,icount)                   &
                                + smc(land_cnt)
            icount = icount + 1
            resdump(j,k,icount) = resdump(j,k,icount)                   &
                                + canopy_gb(land_cnt)
            icount = icount + 1
          End If

          resdump(j,k,icount) = resdump(j,k,icount) + snodep(j,k)
        icount = icount + 1
          resdump(j,k,icount) = resdump(j,k,icount) + tstar(j,k)
        icount = icount + 1
          resdump(j,k,icount) = resdump(j,k,icount) + zh(j,k)
        icount = icount + 1
          resdump(j,k,icount) = resdump(j,k,icount) + z0msea(j,k)
        icount = icount + 1
          Do i = icount,icount + n_cca_levels
            resdump(j,k,i) = resdump(j,k,i) + cca(j,k,i-icount + 1)
          enddo
        icount = icount + 1
          resdump(j,k,icount) = resdump(j,k,icount) + real(iccb(j,k))
        icount = icount + 1
          resdump(j,k,icount) = resdump(j,k,icount) + real(icct(j,k))

          If (land_sea_mask(j,k)) Then
            Do i = 1, nsoilm_levs
              resdump(j,k,icount + i) = resdump(j,k,icount + i)         &
                                      + smcl(land_cnt,i)
            End Do
          End If

        End Do                   ! j
      End Do                   ! k

      Return
      END SUBROUTINE RESTART_DUMP
!
#endif
