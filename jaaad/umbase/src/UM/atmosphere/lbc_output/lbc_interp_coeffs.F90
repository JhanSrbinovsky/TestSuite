#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Subroutine Interface:

      Subroutine LBC_Interp_Coeffs (                                    &
     &     lbc_size                                                     &
     &,    src_row_len                                                  &
     &,    src_rows                                                     &
     &,    src_delta_lat                                                &
     &,    src_delta_long                                               &
     &,    src_first_lat                                                &
     &,    src_first_long                                               &
     &,    src_pole_lat                                                 &
     &,    src_pole_long                                                &
     &,    src_cyclic                                                   &
     &,    src_rotated                                                  &
     &,    lbc_row_len                                                  &
     &,    lbc_rows                                                     &
     &,    l_var_lbc                                                    &
     &,    lambda_in                                                    &
     &,    phi_in                                                       &
     &,    lbc_delta_lat                                                &
     &,    lbc_delta_long                                               &
     &,    lbc_first_lat                                                &
     &,    lbc_first_long                                               &
     &,    lbc_pole_lat                                                 &
     &,    lbc_pole_long                                                &
     &,    rimwidth                                                     &
     &,    lbc_halo_x                                                   &
     &,    lbc_halo_y                                                   &
     &,    lbc_index_bl                                                 &
     &,    lbc_index_br                                                 &
     &,    lbc_weights_tr                                               &
     &,    lbc_weights_br                                               &
     &,    lbc_weights_bl                                               &
     &,    lbc_weights_tl                                               &
     &,    lbc_coeff1                                                   &
     &,    lbc_coeff2                                                   &
     &,    i_uv                                                         &
     & )

      Implicit NONE
!
! Description:
!   Calculates the Horizontal Interpolation Coefficients for LBCs.
!
! Method:
!   1. Calculate lat/longs of lbc points on rotated grid.
!   2. Call EQTOLL to get corresponding true lat/longs.
!   3. Call W_COEFF to calculate coefficients to rotate the winds.
!   4. Set up lat/longs for the source grid.
!   5. Call H_INT_CO to calculate the horizontal interpolation
!      coefficients.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    09/10/01   New argument, src_cyclic. Correct calculation of
!                     lat/long for Northern side. Correct w_coeff
!                     argument list. D Robinson
!   5.5    06/08/02   New arguments - src_rotated, src_pole_lat and
!                     src_pole_long. Add call to LLTOEQ. D Robinson
!   6.1    28/07/04   Remove surplus checks on longitude ranges.
!                     D Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"

! Subroutine arguments

      Integer :: lbc_size

!     Source Grid
      Integer :: src_row_len
      Integer :: src_rows
      Real    :: src_delta_lat
      Real    :: src_delta_long
      Real    :: src_first_lat
      Real    :: src_first_long
      Real    :: src_pole_lat
      Real    :: src_pole_long
      Logical :: src_cyclic
      Logical :: src_rotated

!     LBC Grid
      Integer :: lbc_row_len
      Integer :: lbc_rows
      Integer :: lbc_halo_x
      Integer :: lbc_halo_y
      Real    :: lbc_delta_lat
      Real    :: lbc_delta_long
      Real    :: lbc_first_lat
      Real    :: lbc_first_long
      Real    :: lbc_pole_lat
      Real    :: lbc_pole_long
      
      Logical l_var_lbc
! Input VarRes grid info in degrees 
      Real    :: Lambda_in ( 1-lbc_halo_x: lbc_row_len+lbc_halo_x )  
      Real    :: Phi_in ( 1-lbc_halo_y: lbc_rows+lbc_halo_y )
      
      Integer :: RimWidth

      Integer, dimension (lbc_size) :: lbc_index_bl
      Integer, dimension (lbc_size) :: lbc_index_br
      Real,    dimension (lbc_size) :: lbc_weights_tr
      Real,    dimension (lbc_size) :: lbc_weights_br
      Real,    dimension (lbc_size) :: lbc_weights_tl
      Real,    dimension (lbc_size) :: lbc_weights_bl
      Real,    dimension (lbc_size) :: lbc_coeff1
      Real,    dimension (lbc_size) :: lbc_coeff2

      Integer :: i_uv   !  If 1 => Wind field (u or v). Otherwise 0.

! Local parameters:

      Character (Len=*), Parameter :: RoutineName = 'LBC_Interp_Coeffs'

! Local scalars:

      Integer :: ipt         ! LBC point number
      Integer :: row,pt      ! Loop indices for row and point
      Integer :: iside       ! Loop index for LBC sides
      Integer :: lbc_len     ! Computed no of lbc points
      Integer :: ErrorStatus ! Error Code

      Character (Len=80) :: CMessage

! Local dynamic arrays:

      Real, dimension (:), allocatable :: lambda_lbc
      Real, dimension (:), allocatable :: phi_lbc
      Real, dimension (:), allocatable :: lambda_targ
      Real, dimension (:), allocatable :: phi_targ
      Real, dimension (:), allocatable :: lambda_source
      Real, dimension (:), allocatable :: phi_source
!
! Function & Subroutine calls:
      External EqToLL, H_Int_Co, LLToEq, W_Coeff
!
!- End of header

      ErrorStatus = 0
      CMessage = ' '

! --------------------
! Allocate work arrays
! --------------------

      allocate ( lambda_targ (lbc_size) )
      allocate ( phi_targ    (lbc_size) )
      allocate ( lambda_lbc  (lbc_size) )
      allocate ( phi_lbc     (lbc_size) )

! ---------------------------
! Get lat/Longs of LBC points
! ---------------------------
      If(l_var_lbc) Then 
        
        ipt=0  
   
        Do ISide = 1, 4 
 
          If (ISide == PSouth) Then    !  Southern Boundary (SP)  
 
            Do Row = 1 - lbc_halo_y, RimWidth + i_uv  
              Do Pt = 1 - lbc_halo_x, lbc_row_len + lbc_halo_x   
   
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row)  

              End Do  
            End Do 
  
          End If  !  South  
  
          If (ISide == PEast) Then    !  Eastern Boundary   
  
            Do Row = RimWidth + 1, lbc_rows - RimWidth  
              Do Pt = lbc_row_len - i_uv - RimWidth + 1,                &
     &              lbc_row_len + lbc_halo_x   
  
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row)      

              End Do           
            End Do             
   
          End If  ! East   
 
          If (ISide == PNorth) Then   !  Northern Boundary (NP)  
            Do Row = lbc_rows - RimWidth + 1 - i_uv,                    &
     &             lbc_rows + lbc_halo_y  
              Do Pt = 1-lbc_halo_x, lbc_row_len+lbc_halo_x  
  
              ipt = ipt + 1 
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = phi_in(row)             

              End Do              
            End Do          
  
          End If  ! North 
  
          If (ISide == PWest) Then   !  Western Boundary  
  
            Do Row = RimWidth + 1, lbc_rows - RimWidth 
              Do Pt = 1 - lbc_halo_x, RimWidth + i_uv 
    
                ipt = ipt + 1  
                Lambda_lbc(ipt) = Lambda_in(pt) 
                Phi_lbc(ipt)    = Phi_in(row) 

              End Do 
            End Do 
  
          End If  ! West 
  
        End Do  !  ISide 

      Else     ! regular grid
      
      ipt=0

      Do ISide = 1, 4

        If (ISide == PSouth) Then    !  Southern Boundary (SP)

          Do Row = 1 - lbc_halo_y, RimWidth + i_uv
            Do Pt = 1 - lbc_halo_x, lbc_row_len + lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  !  South

        If (ISide == PEast) Then    !  Eastern Boundary

          Do Row = RimWidth + 1, lbc_rows - RimWidth
            Do Pt = lbc_row_len - i_uv - RimWidth + 1,                  &
     &              lbc_row_len + lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! East

        If (ISide == PNorth) Then   !  Northern Boundary (NP)
          Do Row = lbc_rows - RimWidth + 1 - i_uv,                      &
     &             lbc_rows + lbc_halo_y
            Do Pt = 1-lbc_halo_x, lbc_row_len+lbc_halo_x

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc (ipt)   = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! North

        If (ISide == PWest) Then   !  Western Boundary

          Do Row = RimWidth + 1, lbc_rows - RimWidth
            Do Pt = 1 - lbc_halo_x, RimWidth + i_uv

            ipt = ipt + 1
            Lambda_lbc(ipt) = lbc_first_long + (pt-1)  * lbc_delta_long
            Phi_lbc(ipt)    = lbc_first_lat  + (row-1) * lbc_delta_lat

            End Do
          End Do

        End If  ! West

      End Do  !  ISide
      
      End If  ! variable resolution
      
      lbc_len = ipt

! -----------------------------------------------------
! Check no of lbc points computed ; must match lbc_size
! -----------------------------------------------------

      If (lbc_len /= lbc_size) Then
        write (6,*) ' Mismatch in number of LBC points.'
        write (6,*) ' Expected no of lbc points ',lbc_size
        write (6,*) ' Computed no of lbc points ',lbc_len
        write (CMessage,*) ' Mismatch in number of LBC points.'
        ErrorStatus = 10
! DEPENDS ON: ereport
        Call Ereport ( RoutineName, ErrorStatus, CMessage)
      End If

! -------------------------------------
! Get true lat and longs for lbc points
! -------------------------------------

! DEPENDS ON: eqtoll
      Call EqToLL (                                                     &
     &   phi_lbc                                                        &
     &,  lambda_lbc                                                     &
     &,  phi_targ                                                       &
     &,  lambda_targ                                                    &
     &,  lbc_pole_lat                                                   &
     &,  lbc_pole_long                                                  &
     &,  lbc_len                                                        &
     & )

! ----------------------------------------------
! Calculate the coefficients to rotate the winds
! ----------------------------------------------

      If (i_uv == 1) Then

! DEPENDS ON: w_coeff
        Call W_Coeff (                                                  &
     &       lbc_coeff1                                                 &
     &,      lbc_coeff2                                                 &
     &,      lambda_targ                                                &
     &,      lambda_lbc                                                 &
     &,      lbc_pole_lat                                               &
     &,      lbc_pole_long                                              &
     &,      lbc_len                                                    &
     & )

      End If

! -----------------------------------------------------------
! For rotated model grids, get lat/longs w.r.t the model grid
! -----------------------------------------------------------

      If (src_rotated) Then

! DEPENDS ON: lltoeq
        Call LLToEq (                                                   &
     &       phi_targ                                                   &
     &,      lambda_targ                                                &
     &,      phi_targ                                                   &
     &,      lambda_targ                                                &
     &,      src_pole_lat                                               &
     &,      src_pole_long                                              &
     &,      lbc_len                                                    &
     & )

      End If

! -----------------------------------
! Calculate lat/longs for source grid
! -----------------------------------

      allocate ( lambda_source(src_row_len) )
      allocate ( phi_source   (src_rows)    )

      Do pt=1,src_row_len
        Lambda_Source(pt) = src_first_long + src_delta_long * (pt-1)
      Enddo
      Do row=1,src_rows
        Phi_Source(row)   = src_first_lat  + src_delta_lat  * (row-1)
      Enddo

! --------------------------------------------------
! Calculate the horizontal interpolation cofficients
! --------------------------------------------------

! DEPENDS ON: h_int_co
      Call H_Int_Co (                                                   &
     &   lbc_index_bl                                                   &
     &,  lbc_index_br                                                   &
     &,  lbc_weights_tr                                                 &
     &,  lbc_weights_br                                                 &
     &,  lbc_weights_tl                                                 &
     &,  lbc_weights_bl                                                 &
     &,  Lambda_Source                                                  &
     &,  Phi_Source                                                     &
     &,  Lambda_Targ                                                    &
     &,  Phi_Targ                                                       &
     &,  src_row_len                                                    &
     &,  src_rows                                                       &
     &,  lbc_len                                                        &
     &,  src_cyclic                                                     &
     & )

      deallocate (phi_lbc)
      deallocate (lambda_lbc)
      deallocate (Lambda_Source)
      deallocate (Phi_Source)
      deallocate (Lambda_Targ)
      deallocate (Phi_Targ)

      Return
      END SUBROUTINE LBC_Interp_Coeffs
#endif
