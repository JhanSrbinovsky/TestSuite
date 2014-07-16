#if defined(ATMOS) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine DERV_INTF_A : Calculates Interface array dimensions.
!
! Subroutine Interface :
!
      SUBROUTINE DERV_INTF_A (TOT_LEN_INTFA_P,TOT_LEN_INTFA_U,          &
     &           MAX_INTF_MODEL_LEVELS,MAX_LBCROW_LENGTH,MAX_LBCROWS,   &
     &           N_INTF_A,U_FIELD,U_FIELD_INTFA)

      IMPLICIT NONE
!
! Description : Calculate array dimensions for boundary data output.
!
! Method : Reads in INTFCNSTA namelist to get grid dimensions of
!          interface area. Calculates array dimensions for boundary
!          data. Also sets dimensions to 1 if no interface areas
!          required to prevent zero dynamic allocation.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.5    03/08/98  Original Code
!   5.2    10/11/00  Cater for 5.x LBCs. D.Robinson
!   6.0    05/09/03  Add new def for makebc. R.Sempers
!   6.0    04/11/03  Correct CPP directives.
!                       A. A. Dickinson
!
! Code Description :
! Language : FORTRAN 77 + common extensions
!
! Declarations :

!     Arguments
      Integer TOT_LEN_INTFA_P   ! OUT  Total length of interface data
                                !      on P grid
      Integer TOT_LEN_INTFA_U   ! OUT  Total length of interface data
                                !      on U grid
      Integer MAX_INTF_MODEL_LEVELS ! OUT  Max no of lbc levels
      Integer MAX_LBCROW_LENGTH     ! OUT  Max no of lbc row_length
      Integer MAX_LBCROWS           ! OUT  Max no of lbc rows
      Integer N_INTF_A          ! IN   No of interface areas
      Integer U_FIELD           ! IN   Dimension of U_FIELD
      Integer U_FIELD_INTFA     ! OUT  Dimension of U_FIELD for dynamic
                                !      allocation.

#include "cmaxsize.h"
#include "cintfa.h"

!  Namelist for atmos interface constants
#include "cnaminfa.h"

!     Local variables
      INTEGER JINTF  !  loop index

!     Read in INTFCNSTA namelist to get output grids for
!     generating boundary data.

      REWIND 5
        LBC_ND (:) = 1
      READ (5,INTFCNSTA)
      REWIND 5

      IF (N_INTF_A >  0) THEN

!       Boundary data to be generated in this run.

        MAX_INTF_MODEL_LEVELS = 0
        MAX_LBCROW_LENGTH     = 0
        MAX_LBCROWS           = 0
        Do JINTF=1,N_INTF_A

          MAX_INTF_MODEL_LEVELS =                                       &
     &        MAX ( MAX_INTF_MODEL_LEVELS , INTF_P_LEVELS(JINTF)+1 )
          MAX_LBCROW_LENGTH =                                           &
     &        MAX ( MAX_LBCROW_LENGTH , INTF_ROW_LENGTH(JINTF)+1 )
          MAX_LBCROWS = MAX ( MAX_LBCROWS , INTF_P_ROWS(JINTF)+1 )  
        ENDDO

!       Check >= 1 to avoid zero dynamic allocation
        MAX_INTF_MODEL_LEVELS = MAX ( MAX_INTF_MODEL_LEVELS , 1)
        MAX_LBCROW_LENGTH     = MAX ( MAX_LBCROW_LENGTH , 1)
        MAX_LBCROWS           = MAX ( MAX_LBCROWS , 1)
        
        TOT_LEN_INTFA_P = 0
        TOT_LEN_INTFA_U = 0
        Do JINTF=1,N_INTF_A

          If (lbc_nd(jintf) == 0) Then  !  Old LBCs.

!         Calculate lengths for interface area JINTF

          LEN_INTFA_P(JINTF) = ( INTF_ROW_LENGTH(JINTF) +               &
     &    INTF_P_ROWS(JINTF) - 2*INTFWIDTHA(JINTF) )                    &
     &    * 2 * INTFWIDTHA(JINTF)
          LEN_INTFA_U(JINTF) = LEN_INTFA_P(JINTF) - 4*INTFWIDTHA(JINTF)

!         Add on to total length

          TOT_LEN_INTFA_P = TOT_LEN_INTFA_P + LEN_INTFA_P(JINTF)
          TOT_LEN_INTFA_U = TOT_LEN_INTFA_U + LEN_INTFA_U(JINTF)

          endif

        ENDDO


!       U_FIELD_INTFA Dimensions COEFF3 & COEFF4 in TYPINFA

        U_FIELD_INTFA = U_FIELD

        if (tot_len_intfa_p == 0 .and. tot_len_intfa_u == 0) Then
          tot_len_intfa_p = 1
          tot_len_intfa_u = 1
          u_field_intfa   = 1
        endif

        write (6,*) ' '
        write (6,*) ' Data lengths calculated in DERV_INTF_A.'
        do jintf=1,n_intf_a
        write (6,*) ' Area no ',jintf,                                  &
     &              ' lbc_nd ',lbc_nd(jintf),                           &
     &              ' len_intfa_p ',len_intfa_p(jintf),                 &
     &              ' len_intfa_u ',len_intfa_u(jintf)
        enddo
        write (6,*) ' n_intf_a ',n_intf_a
        write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
        write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
        write (6,*) ' max_intf_model_levels ',max_intf_model_levels
        write (6,*) ' max_lbcrow_length ',max_lbcrow_length
        write (6,*) ' max_lbcrows ',max_lbcrows
        write (6,*) ' u_field_intfa ',u_field_intfa

      ELSE

!       No boundary conditions to be generated.
!       Initialise to prevent zero length dynamic allocation.

        write (6,*) ' n_intf_a ',n_intf_a

        N_INTF_A = 1
        TOT_LEN_INTFA_P = 1
        TOT_LEN_INTFA_U = 1
        MAX_INTF_MODEL_LEVELS = 1
        MAX_LBCROW_LENGTH = 1
        MAX_LBCROWS = 1
        U_FIELD_INTFA = 1

      write (6,*) ' n_intf_a ',n_intf_a
      write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
      write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
      write (6,*) ' max_intf_model_levels ',max_intf_model_levels
      write (6,*) ' max_lbcrow_length ',max_lbcrow_length
      write (6,*) ' max_lbcrows ',max_lbcrows
      write (6,*) ' u_field_intfa ',u_field_intfa

      ENDIF

      RETURN
      END SUBROUTINE DERV_INTF_A
#endif
