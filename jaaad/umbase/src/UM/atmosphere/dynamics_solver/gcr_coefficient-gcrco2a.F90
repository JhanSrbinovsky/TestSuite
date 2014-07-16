#if defined(A10_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_Coefficient

      Subroutine GCR_Coefficient(                                       &
     &                           first_term, second_term,               &
     &                           row_length, rows, model_levels,        &
     &                           model_domain, coefficient,             &
     &                           at_extremity, n_proc,                  &
     &                           gc_proc_col_group, gc_proc_row_group,  &
     &                           l_datastart                            &
     &                           )

! Purpose:
!          Calculates a coefficient given by
!
!                              < first_term, second_term >
!          inner_product = -   __________________________
!
!                              < second_term, second_term >
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/03/00 5.1        Changed j_start/stop for cyclic LAM  Andy Malcolm
!LL   5.1   10/02/00  Use DOMTYP parameters                    P.Burton
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!     5.3   19/10/01  Use appropriate gcg routines.   S. Cusack
!     6.2   21/10/05  Remove unused external.              P.Selwood
!     6.2   19/10/05  Unroll loop for better vectorisation.
!                     K. Ketelsen/P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels     ! number of model levels.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  first_term(row_length,rows,model_levels)                        &
                                                  ! first term in inner
                                                  ! product definition
     &, second_term(row_length,rows,model_levels)
                                                  ! second term in inner
                                                  ! product definition

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

      Real                                                              &
     &  coefficient

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      Integer                                                           &
     &  l_datastart(3)

      Integer                                                           &
     &  istat                                                           &
     &, gc_proc_col_group                                               &
     &, gc_proc_row_group                                               &
     &, n_proc

! parallel Local arrays
      Real                                                              &
     &  inner_product(2)                                                &
     &, term(row_length,rows,2)                                         &
     &, sum_temp(rows,2)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

#include "parparm.h"
#include "domtyp.h"

#if defined(NEC)
! Variables required for hand-vectorisation. 
      Integer, Parameter :: vec_len = 256   ! vector length
      Integer            :: iii
      Real               :: v1( vec_len )
      Real               :: v2( vec_len )
!CDIR VREG(v1,v2)
#endif


! ----------------------------------------------------------------------
! Section 1.   Calculate Error Norm.
! ----------------------------------------------------------------------

! top term is in term (i,j,1)
! bottom term is in term (i,j,2)
      Do i= 1, row_length
        Do j= 1, rows
          term(i,j,1)=0.0
          term(i,j,2)=0.0
        End Do
      End Do

      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) Then
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
      If (model_domain  ==  mt_global) Then
        if(at_extremity(PSouth))j_start=2
        if(at_extremity(PNorth))j_stop=rows-1
      Endif
      Else If(model_domain  ==  mt_lam)then
! Solve over interior points only.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
        if(at_extremity(PEast)) i_stop = row_length-1
        if(at_extremity(PWest)) i_start = 2
      Else
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
      End If

      If (model_domain  ==  mt_global) Then
! Global model only calculate inner product for one of the polar points.

        If(at_extremity(PSouth).and.(l_datastart(1) == 1))then
          Do k = 1, model_levels
            term(1,1,1)= term(1,1,1) +                                  &
     &                   first_term(1,1,k) * second_term(1,1,k)
            term(1,1,2)= term(1,1,2) +                                  &
     &                   second_term(1,1,k) * second_term(1,1,k)
          End Do
        End If
        if(at_extremity(PNorth).and.(l_datastart(1) == 1))then
          Do k = 1, model_levels
            term(1,rows,1)= term(1,rows,1) + first_term(1,rows,k)       &
     &                      * second_term(1,rows,k)
            term(1,rows,2)= term(1,rows,2) + second_term(1,rows,k)      &
     &                      * second_term(1,rows,k)
          End Do
        End If
      End If

#if defined(NEC)
! If on an NEC we can hand strip-mine the loops to get best
! vector usage. Note that v1 and v2 are vector registers 
! hence limited to 256 words

      Do iii=i_start, i_stop, vec_len
        Do j = j_start , j_stop

! First load up the vector registers with term
!CDIR SHORTLOOP
          Do i = iii,min(iii+vec_len-1,i_stop) 
            v1(i+1-iii)=term(i,j,1)
            v2(i+1-iii)=term(i,j,2)
          End Do

          Do k = 1, model_levels
! Now do the inner products
!CDIR SHORTLOOP
            Do i = iii,min(iii+vec_len-1,i_stop) 
              v1(i+1-iii)=v1(i+1-iii) + first_term(i,j,k)           &
     &                                * second_term(i,j,k)
              v2(i+1-iii)=v2(i+1-iii) + second_term(i,j,k)          &
     &                                * second_term(i,j,k)
            End Do
          End Do

! Finally store the result back in term
!CDIR SHORTLOOP
          Do i = iii,min(iii+vec_len-1,i_stop) 
            term(i,j,1)=v1(i+1-iii)
            term(i,j,2)=v2(i+1-iii)
          End Do
        End Do
      End Do
#else
      Do k = 1, model_levels
        Do j = j_start , j_stop
          Do i = i_start, i_stop
            term(i,j,1) = term(i,j,1) + first_term(i,j,k)               &
     &                                * second_term(i,j,k)
            term(i,j,2) = term(i,j,2) + second_term(i,j,k)              &
     &                                * second_term(i,j,k)
          End Do
        End Do
      End Do
#endif

#if defined(REPROD)
      Call gcg_rvecsumr(row_length, row_length, 1, 2*rows, term,        &
     &                  gc_proc_row_group, istat, sum_temp)
      Call gcg_rvecsumr(rows, rows, 1, 2, sum_temp,                     &
     &                  gc_proc_col_group, istat, inner_product)

#else
! If non-reproducible we can sum each processors contribution
! to the global sum first before a final summation of these
! partial sums to get the true global sum.

!CDIR COLLAPSE
      inner_product(:) = 0.0
      Do j = 1, rows
        Do i = 1, row_length
          inner_product(1) = inner_product(1) + term(i,j,1)
          inner_product(2) = inner_product(2) + term(i,j,2)
        End Do
      End Do

      Call gc_rsum(2,n_proc,istat,inner_product)
#endif

! Calculate coefficient

      coefficient = - inner_product(1) / inner_product(2)

! End of routine
      return
      END SUBROUTINE GCR_Coefficient

#endif
