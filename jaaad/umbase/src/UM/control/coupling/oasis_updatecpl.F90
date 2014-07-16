#if defined(OASIS3) || defined(OASIS4)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_updatecpl(                                        &                       
#include "argd1.h" 
#include "arg_atm_fields.h"
  cmessage)


#if defined(OASIS3)
  USE oasis3_atm_data_mod, ONLY : um_co2, oasis_jmt, oasis_jmt_u             &
     ,oasis_jmt_v, oasis_imt
#endif

#if defined(OASIS4)
  USE oasis4_atm_data_mod 
#endif

  IMPLICIT NONE

  !
  ! Description:
  ! Update coupling field data into prognostic fields
  ! to ensure restartability. We do this simply by copying
  ! data from diagnostic (stash) areas at the appropriate timestep.
  !
  ! Author: R. Hill
  ! Current Code Owner : Richard Hill
  !
  !-------------------------------------------------------------------
#include "parvars.h"   
#include "decomptp.h" 
#include "decompdb.h" 
#include "atm_lsm.h"  
#include "cntlatm.h" 

#include "cmaxsize.h" 
#include "typsize.h"  
#include "typd1.h"    
#include "typptra.h"  
#include "ctracera.h" 
#include "typ_atm_fields.h"

#include "nstypes.h"
#include "cruntimc.h"
#include "caoptr.h"   

  INTEGER :: i              ! Loop counter
  INTEGER :: j              ! Loop counter
  INTEGER :: k              ! Ice category counter
  INTEGER :: oasis_error, ft
  INTEGER :: icode          ! Error return code (=0 is OK)
  INTEGER :: info           ! Return code from MPP
  INTEGER :: gather_pe      ! Processor for gathering
  CHARACTER*(*) :: cmessage ! OUT - Error return message

  INTEGER :: pointer_top_d  ! Pointers to D1 for stash data
  INTEGER :: pointer_bot_d  ! multi-cat top/botmelt respectively

  ! Our row lengths and column lengths must be the local domain sizes
  ! for each PE.  


  ! Set up fields from the U grid
  DO j=1,oasis_jmt_u
    DO i=1,oasis_imt
      c_taux(i,j) =                                               &
         D1(ja_taux+i-1+((j-1)*oasis_imt))
      c_u10(i,j)=D1(jc_u10+i-1+((j-1)*oasis_imt))
    END DO
  END DO

!calculate u-comp wind at v point on C-grid
! DEPENDS ON: ub_to_uvc
          CALL  uB_to_uVC(c_u10,                                      &
     &                   row_length,rows,n_rows,1,offx,offy,       &
     &                   c_u10)



  ! Set up fields from the V grid
  DO j=1,oasis_jmt_v
    DO i=1,oasis_imt
      c_tauy(i,j) =                                               &
         D1(ja_tauy+i-1+((j-1)*oasis_imt))
      c_v10(i,j)=D1(jc_v10+i-1+((j-1)*oasis_imt))
    END DO
  END DO

  ! Set up fields from the T grid
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      ! co2
      if (l_co2_interactive) then
        um_co2(i,j) = D1(jc_co2+i-1+((j-1)*oasis_imt))  
        !um_co2(i,j) = co2(i,j,1)
      else
        um_co2(i,j) = co2_mmr
      endif


      ! Copy various heat flux components 
      c_solar(i,j) =                                              &
         D1(ja_solar+i-1+((j-1)*oasis_imt))

      c_blue(i,j) =                                               &
         D1(ja_blue+i-1+((j-1)*oasis_imt))

      c_sublim(i,j) =                                             &
         D1(ja_sublim+i-1+((j-1)*oasis_imt))

      c_longwave(i,j) =                                           &
         D1(ja_longwave+i-1+((j-1)*oasis_imt))

      c_sensible(i,j) =                                           &
         D1(ja_sensible+i-1+((j-1)*oasis_imt))

      ! PME components
      c_evap(i,j) =                                               &
         D1(ja_evap+i-1+((j-1)*oasis_imt))

      c_lsrain(i,j) =                                             &
         D1(ja_lsrain+i-1+((j-1)*oasis_imt))

      c_lssnow(i,j) =                                             &
         D1(ja_lssnow+i-1+((j-1)*oasis_imt))

      c_cvrain(i,j) =                                             &    
         D1(ja_cvrain+i-1+((j-1)*oasis_imt))

      c_cvsnow(i,j) =                                             &
         D1(ja_cvsnow+i-1+((j-1)*oasis_imt))

      ! River runoff
      c_riverout(i,j) =                                           & 
         D1(ja_riverout+i-1+((j-1)*oasis_imt))


      ! Wind mixing energy WME
      c_windmix(i,j) =                                            &
         D1(ja_windmix+i-1+((j-1)*oasis_imt))

#if defined(ACCESS)
! gol124: auscom coupling
      ! Surface pressure
      c_press(i,j) =                                              &
         D1(ja_press+i-1+((j-1)*oasis_imt))
#endif

    END DO
  END DO

  IF (.NOT.L_CTILE) THEN                           
    !          In the non-coastal tiling case we convert sublimation
    !          from a total to a rate. Note that if coastal tiling is 
    !          not used (l_ctile false), field 3231 (sublimation total) is 
    !          set by oasis_inita2o to sublim and needs to be converted to 
    !          a rate here. If coastal tiling is used (l_ctile true), field 
    !          3353 (sublimation rate) is set by oasis_inita2o so no extra 
    !          processing is required.          
    !          (See UM7.0 or earlier swap_a2o for historical perspective
    !          on the way these fields are handled).                        
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_sublim(i,j)=c_sublim(i,j)/                            &
           REAL(OASIS_COUPLE_FREQ*3600)         
      END DO
    END DO
  END IF

  ! Topmelt and Botmelt (multi-category)
  pointer_top_d=ja_topmeltn
  pointer_bot_d=ja_botmeltn

  DO k=1,nice
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_topmeltn(i,j,k)=D1(pointer_top_d)
        c_botmeltn(i,j,k)=D1(pointer_bot_d)
        pointer_top_d=pointer_top_d+1
        pointer_bot_d=pointer_bot_d+1
      END DO
    END DO
  END DO


END SUBROUTINE oasis_updatecpl
#endif
