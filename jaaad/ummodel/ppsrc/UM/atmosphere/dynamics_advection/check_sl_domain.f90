
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Check_sl_domain.

      Subroutine Check_sl_domain(                                       &
     &                        model_domain, depart_phi, depart_lambda,  &
     &                        row_length, rows_depart, model_levels,    &
     &                        domain_size_x, domain_size_y,             &
     &                        max_lambda, min_lambda,                   &
     &                        max_phi, min_phi, Pi)

! Purpose:
!          Checks data lies inside model domain.
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
! Date     Version     Comment
! ----     -------     -------
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
!LL   5.2   27/09/00  Fix bug for cyclic lam code             A.Malcolm
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                     ! Dimension of u_adv array in i direction.
     &, rows_depart                                                     &
                     ! Dimension of depart arrays in j direction.
     &, model_levels                                                    &
                     ! Dimension of u_adv array in k direction.
     &, model_domain ! holds integer code for model domain

      Real                                                              &
     &  Pi                                                              &
     &, max_lambda                                                      &
     &, max_phi                                                         &
     &, min_lambda                                                      &
     &, min_phi

      Real :: domain_size_x
      Real :: domain_size_y

! Arguments with Intent IN/OUT.

      Real                                                              &
                      ! Departure point co-ordinates.
     &  depart_lambda (row_length, rows_depart, model_levels)           &
     &, depart_phi (row_length, rows_depart, model_levels)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k      ! Loop indices

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
! External Routines: None

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Check points lie inside model domain.
! ----------------------------------------------------------------------

      If (model_domain  ==  mt_Global) Then
! check horizontal
        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = -pi -  depart_phi(i,j,k)
                depart_lambda(i,j,k) = Pi + depart_lambda(i,j,k)
              Else If (depart_phi(i,j,k)  >   max_phi ) Then
                depart_phi(i,j,k) = Pi -  depart_phi(i,j,k)
                depart_lambda(i,j,k) = Pi + depart_lambda(i,j,k)
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = Pi*2. + depart_lambda(i,j,k)
              Else If (depart_lambda(i,j,k)  >=  max_lambda ) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k) - Pi*2.
              End If
            End Do
          End Do
        End Do

      Else If (model_domain  ==  mt_LAM) Then

! check horizontal in the LAM at u points.

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = min_phi
              Else If (depart_phi(i,j,k) >  max_phi)Then
                depart_phi(i,j,k) = max_phi
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = min_lambda
              Else If(depart_lambda(i,j,k) >  max_lambda) Then
                depart_lambda(i,j,k) = max_lambda
              End If
            End Do
          End Do
        End Do

      Else If (model_domain  ==  mt_cyclic_LAM) Then

! check horizontal in the periodic in x LAM at u points.

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi ) Then
                depart_phi(i,j,k) = min_phi
              Else If (depart_phi(i,j,k) >  max_phi)Then
                depart_phi(i,j,k) = max_phi
              End If
              If (depart_lambda(i,j,k)  <   min_lambda ) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               + domain_size_x             
              Else If(depart_lambda(i,j,k) >= max_lambda) Then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               - domain_size_x            
              End If
            End Do
          End Do
        End Do


      Else If (model_domain  ==  mt_bi_cyclic_lam ) Then
! check horizontal in the periodic in both x and y LAM at u points

        Do k = 1, model_levels
          Do j = 1, rows_depart
            Do i = 1, row_length
              If (depart_phi(i,j,k)  <   min_phi) then
                depart_phi(i,j,k) = depart_phi(i,j,k)                   &
     &                            + domain_size_y     
              Elseif (depart_phi(i,j,k)  >   max_phi) then
                depart_phi(i,j,k) = depart_phi(i,j,k)                   &
     &                            - domain_size_y     
              Endif
              If (depart_lambda(i,j,k)  <   min_lambda) then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               + domain_size_x          
              Elseif (depart_lambda(i,j,k)  >   max_lambda) then
                depart_lambda(i,j,k) = depart_lambda(i,j,k)             &
     &                               - domain_size_x          
              Endif
            Enddo
          Enddo
        Enddo

      End if

! End of routine.
      return
      END SUBROUTINE Check_sl_domain

