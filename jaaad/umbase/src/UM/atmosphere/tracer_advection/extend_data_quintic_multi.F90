#if defined(A11_2A) || defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Extend_data_quintic_multi.

      Subroutine Extend_data_quintic_multi(                             &
     &                             Data_in,                             &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             number_of_inputs, halo_i, halo_j,    &
     &                             Ext_data)

! Purpose:
!          Extends data arrays into bigger arrays to
!          allow use of more efficient cubic interpolation.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date       Comment
! ----     -------     -------
!  6.1     24/06/04    Original deck for super_arrays      A.Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, number_of_inputs !number of fields to interpolate.

      Real                                                              &
     &  Data_in (1-halo_i:dim_i_in+halo_i,                              &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in                    &
     &            ,number_of_inputs)                  ! data to be
                                                      ! interpolated

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,  &
     &            -1:dim_k_in+2, number_of_inputs)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, n     ! Loop indices

! External Routines: None

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Extend input data to bigger area to
!              allow interpolation to be done without having to re-do
!              any end points.
! ----------------------------------------------------------------------

      do n=1,number_of_inputs
        Do k = 1, dim_k_in
          Do j = 1-halo_j, dim_j_in+halo_j
            Do i = 1-halo_i, dim_i_in+halo_i
              Ext_Data (i,j,k,n) = Data_in (i,j,k,n)
            End Do
          End Do
        End Do
      End Do

! Extend in k direction.
      Do n = 1, number_of_inputs
        Do j = 1-halo_j,dim_j_in+halo_j
          Do i = 1-halo_i, dim_i_in+halo_i
            Ext_Data(i,j,-1,n) = Ext_Data(i,j,1,n)
            Ext_Data(i,j,0,n) = Ext_Data(i,j,1,n)
            Ext_Data(i,j,dim_k_in+1,n) = Ext_Data(i,j,dim_k_in,n)
            Ext_Data(i,j,dim_k_in+2,n) = Ext_Data(i,j,dim_k_in,n)
          End Do
        End Do
      End Do

! End of routine.
      return
      END SUBROUTINE Extend_data_quintic_multi

#endif
