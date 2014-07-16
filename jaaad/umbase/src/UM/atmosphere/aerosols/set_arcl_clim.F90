! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine for the aerosol climatology for NWP
!
! Purpose:
!  Copy individual fields of the aerosol climatology for NWP into 
!  a single array. Aerosol species that are not requested are not 
!  copied. The array index corresponding to each requested component
!  are stored for latter access to the corresponding field.
!
! Method:
!   Straightforward.
!
! Current owner of code: N. Bellouin
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
subroutine set_arcl_clim(                                               &
                         ! Array dimensions
                         row_length,                                    &
                         rows,                                          &
                         model_levels,                                  &
                         n_arcl_compnts,                                &
                         ! Internal model switches
                         L_USE_ARCL,                                    &
                         ! Climatologies from ancillary files:
                         !    biomass-burning
                         arclbiom_fr, arclbiom_ag, arclbiom_ic,         &
                         !    black-carbon
                         arclblck_fr, arclblck_ag,                      &
                         !    sea-salt
                         arclsslt_fi, arclsslt_jt,                      &
                         !    sulphate
                         arclsulp_ac, arclsulp_ak, arclsulp_di,         &
                         !    mineral dust
                         arcldust_b1, arcldust_b2, arcldust_b3,         &
                         arcldust_b4, arcldust_b5, arcldust_b6,         &
                         !    fossil-fuel organic carbon
                         arclocff_fr, arclocff_ag, arclocff_ic,         &
                         !    delta aerosol
                         arcldlta_dl,                                   &
                         ! Internal climatology array
                         arcl,                                          &
                         ! Component array indices
                         i_arcl_compnts                                 &
                        )

  implicit none

#include "arcl_dim.h"

  !
  ! Arguments with intent(in)
  !

  !
  ! Array dimensions
  !
  integer, intent(in) ::                                                &
                         row_length,                                    &
                         rows,                                          &
                         model_levels,                                  &
                         n_arcl_compnts
  
  !
  ! Internal model switches 
  !
  logical, dimension(NPD_ARCL_SPECIES), intent(in) :: L_USE_ARCL
  
  !
  ! Climatologies from ancillary files
  !
  !    Three components of biomass-burning
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arclbiom_fr,                                                &
            arclbiom_ag,                                                &
            arclbiom_ic
  !    
  !    Two components of black-carbon
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arclblck_fr,                                                &
            arclblck_ag
  !
  !    Two components of sea-salt
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arclsslt_fi,                                                &
            arclsslt_jt
  !
  !    Three components of sulphate
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arclsulp_ac,                                                &
            arclsulp_ak,                                                &
            arclsulp_di
  !
  !    Six components of mineral dust
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arcldust_b1,                                                &
            arcldust_b2,                                                &
            arcldust_b3,                                                &
            arcldust_b4,                                                &
            arcldust_b5,                                                &
            arcldust_b6

  !
  !    Three components of fossil-fuel organic carbon
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arclocff_fr,                                                &
            arclocff_ag,                                                &
            arclocff_ic
  
  !
  !    One component of delta aerosol
  !
  real, dimension(row_length, rows, model_levels), intent(in) ::        &
            arcldlta_dl

  !
  ! Arguments with intent(out)
  !
  
  !
  ! Internal climatology array
  !
  real, dimension(row_length, rows, model_levels, n_arcl_compnts),      &
    intent(out) :: arcl
    
  !
  ! Component array indices
  !
  integer, dimension(NPD_ARCL_COMPNTS), intent(out) :: i_arcl_compnts
    
  !
  ! Local variables
  !
#include "arcl_ids.h"
  
  integer i_cmp
  
  integer i, j, k
  
  i_cmp = 1 ! since this routine has been called, there is at least
            ! one component to process.
  
  !
  ! For each requested species, copy the corresponding components
  ! into the gathering array.
  ! We also keep track of the index where component has been put.
  ! An index of -1 denotes components that are not used.
  !
  if (L_USE_ARCL(IP_ARCL_SULP)) then
    
    i_arcl_compnts(IP_ARCL_SULP_AC) = i_cmp
    i_arcl_compnts(IP_ARCL_SULP_AK) = i_cmp + 1
    i_arcl_compnts(IP_ARCL_SULP_DI) = i_cmp + 2
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arclsulp_ac(i, j, k)
          arcl(i, j, k, i_cmp+1) = arclsulp_ak(i, j, k)
          arcl(i, j, k, i_cmp+2) = arclsulp_di(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 3
    
  else
    
    i_arcl_compnts(IP_ARCL_SULP_AC) = -1
    i_arcl_compnts(IP_ARCL_SULP_AK) = -1
    i_arcl_compnts(IP_ARCL_SULP_DI) = -1
    
  end if 

  if (L_USE_ARCL(IP_ARCL_DUST)) then
    
    i_arcl_compnts(IP_ARCL_DUST_B1) = i_cmp
    i_arcl_compnts(IP_ARCL_DUST_B2) = i_cmp + 1
    i_arcl_compnts(IP_ARCL_DUST_B3) = i_cmp + 2
    i_arcl_compnts(IP_ARCL_DUST_B4) = i_cmp + 3
    i_arcl_compnts(IP_ARCL_DUST_B5) = i_cmp + 4
    i_arcl_compnts(IP_ARCL_DUST_B6) = i_cmp + 5
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arcldust_b1(i, j, k)
          arcl(i, j, k, i_cmp+1) = arcldust_b2(i, j, k)
          arcl(i, j, k, i_cmp+2) = arcldust_b3(i, j, k)
          arcl(i, j, k, i_cmp+3) = arcldust_b4(i, j, k)
          arcl(i, j, k, i_cmp+4) = arcldust_b5(i, j, k)
          arcl(i, j, k, i_cmp+5) = arcldust_b6(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 6
  
  else
    
    i_arcl_compnts(IP_ARCL_DUST_B1) = -1
    i_arcl_compnts(IP_ARCL_DUST_B2) = -1
    i_arcl_compnts(IP_ARCL_DUST_B3) = -1
    i_arcl_compnts(IP_ARCL_DUST_B4) = -1
    i_arcl_compnts(IP_ARCL_DUST_B5) = -1
    i_arcl_compnts(IP_ARCL_DUST_B6) = -1
    
  end if

  if (L_USE_ARCL(IP_ARCL_SSLT)) then
    
    i_arcl_compnts(IP_ARCL_SSLT_FI) = i_cmp
    i_arcl_compnts(IP_ARCL_SSLT_JT) = i_cmp + 1
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arclsslt_fi(i, j, k)
          arcl(i, j, k, i_cmp+1) = arclsslt_jt(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 2
  
  else
  
    i_arcl_compnts(IP_ARCL_SSLT_FI) = -1
    i_arcl_compnts(IP_ARCL_SSLT_JT) = -1
    
  end if 

  if (L_USE_ARCL(IP_ARCL_BLCK)) then
    
    i_arcl_compnts(IP_ARCL_BLCK_FR) = i_cmp
    i_arcl_compnts(IP_ARCL_BLCK_AG) = i_cmp + 1
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arclblck_fr(i, j, k)
          arcl(i, j, k, i_cmp+1) = arclblck_ag(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 2
    
  else
  
    i_arcl_compnts(IP_ARCL_BLCK_FR) = -1
    i_arcl_compnts(IP_ARCL_BLCK_AG) = -1
  
  end if

  if (L_USE_ARCL(IP_ARCL_BIOM)) then
    
    i_arcl_compnts(IP_ARCL_BIOM_FR) = i_cmp
    i_arcl_compnts(IP_ARCL_BIOM_AG) = i_cmp + 1
    i_arcl_compnts(IP_ARCL_BIOM_IC) = i_cmp + 2
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arclbiom_fr(i, j, k)
          arcl(i, j, k, i_cmp+1) = arclbiom_ag(i, j, k)
          arcl(i, j, k, i_cmp+2) = arclbiom_ic(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 3
    
  else
  
    i_arcl_compnts(IP_ARCL_BIOM_FR) = -1
    i_arcl_compnts(IP_ARCL_BIOM_AG) = -1
    i_arcl_compnts(IP_ARCL_BIOM_IC) = -1
  
  end if

  if (L_USE_ARCL(IP_ARCL_OCFF)) then
    
    i_arcl_compnts(IP_ARCL_OCFF_FR) = i_cmp
    i_arcl_compnts(IP_ARCL_OCFF_AG) = i_cmp + 1
    i_arcl_compnts(IP_ARCL_OCFF_IC) = i_cmp + 2
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arclocff_fr(i, j, k)
          arcl(i, j, k, i_cmp+1) = arclocff_ag(i, j, k)
          arcl(i, j, k, i_cmp+2) = arclocff_ic(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 3
  
  else
  
    i_arcl_compnts(IP_ARCL_OCFF_FR) = -1
    i_arcl_compnts(IP_ARCL_OCFF_AG) = -1
    i_arcl_compnts(IP_ARCL_OCFF_IC) = -1  
  
  end if

  if (L_USE_ARCL(IP_ARCL_DLTA)) then
  
    i_arcl_compnts(IP_ARCL_DLTA_DL) = i_cmp
    
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          
          arcl(i, j, k, i_cmp  ) = arcldlta_dl(i, j, k)
          
        end do ! i
      end do ! j
    end do ! k
    
    i_cmp = i_cmp + 1
  
  else
  
    i_arcl_compnts(IP_ARCL_DLTA_DL) = -1
  
  end if

end subroutine set_arcl_clim
