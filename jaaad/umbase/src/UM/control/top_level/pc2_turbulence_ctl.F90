#if defined(CONTROL) || defined(SCMA)
#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Calculates the effect of changing the width of the PDF used in PC2.

      Subroutine pc2_turbulence_ctl (                                   &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &
!
! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_model_levels                                  &
!
! Model switches
     &, Ltimer, l_mixing_ratio                                          &
!
! Parameters
     &, dbsdtbs_turb_0                                                  &
!
! in time stepping information.
     &, timestep                                                        &
!
! Primary fields passed in
     &, T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                      &
     &, p_layer_centres                                                 &
!
! diagnostic info
     &,                                                                 &
#include "argsts.h"
     &  STASHwork4                                                      &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &
!
! Increment fields passed in/out
     &, T_inc, q_inc, qcl_inc, qcf_inc, cf_inc, cfl_inc, cff_inc        &
!
! error information
     &, Error_code  )

      Implicit None
!
! Description: Condensation and cloud fraction changes in the PC2
!              framework as a result of changing the width of the
!              moisture PDF without changing its shape.
!
! Method:      Uses the equations outlined in the PC2 cloud scheme
!              documentation Annex D.
!
! History:
! Version   Date     Comment.
! -------   ----     --------
! 5.4       16/08/02 Original version. Damian Wilson.
! 6.1       25/05/04 Changes for SCM.  Damian Wilson.
! 6.2       23/01/06 Improvements to Single Column Model Diagnostics
!                    System                          A. Kerr-Munslow
! 6.4       18/08/06 Use mixing ratio formulation.  Damian Wilson
!
! Code Description:
!   Language: Fortran 77 + common extensions
!   This code is written to UMDP3 v6  programming standards.
!
! Declarations:
!
! Arguments with intent IN. ie: input variables.
!
! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number
!
      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
!                          south, east or west of the processor grid
!                          (array index PNorth etc. from parparm.h)
!
! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, wet_model_levels
!
      Logical                                                           &
     &  Ltimer                                                          &
                 ! true then output some timing information
     &, l_mixing_ratio  ! Use mixing ratio formulation
!
! time information for current timestep
      Real                                                              &
     &  timestep
!
! Model parameters
      Real                                                              &
     &  dbsdtbs_turb_0
                        ! IN PC2 erosion rate / s-1
!
! Primary fields passed in
      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, cff_n(row_length, rows, wet_model_levels)
!
      Real                                                              &
     &  p_layer_centres(row_length, rows, 0:model_levels)
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
!
#include "csubmodl.h"
#include "typsts.h"
#include "c_pc2pdf.h"
#if defined(SCMA)
! INOUT SCMop is declared in here
#include "s_scmop.h"
#endif
!
! Diagnostics info
      REAL                                                              &
     & STASHwork4(*)     ! STASH workspace
!
! Additional variables for SCM diagnostics
      INTEGER                                                           &
     &  nSCMDpkgs              ! No of SCM diagnostics packages

      LOGICAL                                                           &
     &  L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages
!
! arguments with intent in/out. ie: input variables changed on output.
      Real                                                              &
     &  T_inc(row_length, rows, model_levels)                           &
     &, q_inc(row_length, rows, wet_model_levels)                       &
     &, qcl_inc(row_length, rows, wet_model_levels)                     &
     &, qcf_inc(row_length, rows, wet_model_levels)                     &
     &, cf_inc(row_length, rows, wet_model_levels)                      &
     &, cfl_inc(row_length, rows, wet_model_levels)                     &
     &, cff_inc(row_length, rows, wet_model_levels)
!
      Integer                                                           &
     &  Error_code
!
! local variables.
#include "fldtype.h"
!
! Local variables
!
      Character*(*), Parameter ::  RoutineName = 'pc2_turbulence_ctl'

      Integer                                                           &
     & i,j,k
! Loop counters

      Real                                                              &
     &  T_work(row_length, rows, model_levels)                          &
     &, q_work(row_length, rows, wet_model_levels)                      &
     &, qcl_work(row_length, rows, wet_model_levels)                    &
     &, cf_work(row_length, rows, wet_model_levels)                     &
     &, cfl_work(row_length, rows, wet_model_levels)                    &
     &, zeros(row_length, rows, wet_model_levels)
!
! External Routines:

!- End of header

!
! Call timer for turbulence code
! DEPENDS ON: timer
      If (Ltimer) Call timer ('PC2 Turbulence',3)
!
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
!
! Work fields are set to time level n fields
!
              q_work(i,j,k)   = q_n(i,j,k)
              qcl_work(i,j,k) = qcl_n(i,j,k)
              cf_work(i,j,k)  = cf_n(i,j,k)
              cfl_work(i,j,k) = cfl_n(i,j,k)
!
! Define zeros array
!
              zeros(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_work(i,j,k)   = T_n(i,j,k)
            End Do
          End Do
        End Do
!
! Call homogenous forcing routine
!
! DEPENDS ON: pc2_homog_plus_turb
        CALL PC2_HOMOG_PLUS_TURB(p_layer_centres(1,1,1),                &
     &      wet_model_levels,                                           &
     &      row_length, rows, timestep, T_work, cf_work,                &
     &      cfl_work,                                                   &
     &      cff_n, q_work, qcl_work, zeros,                             &
     &      zeros, zeros, zeros, dbsdtbs_turb_0, dbsdtbs_turb_1,        &
     &      l_mixing_ratio)
!
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
!
! Redefine work arrays as changes from time level n.
!
              q_work(i,j,k)   = q_work(i,j,k)   - q_n(i,j,k)
              qcl_work(i,j,k) = qcl_work(i,j,k) - qcl_n(i,j,k)
              cf_work(i,j,k)  = cf_work(i,j,k)  - cf_n(i,j,k)
              cfl_work(i,j,k) = cfl_work(i,j,k) - cfl_n(i,j,k)
!
! Update the increment fields
!
              q_inc(i,j,k)   = q_work(i,j,k)   + q_inc(i,j,k)
              qcl_inc(i,j,k) = qcl_work(i,j,k) + qcl_inc(i,j,k)
              cf_inc(i,j,k)  = cf_work(i,j,k)  + cf_inc(i,j,k)
              cfl_inc(i,j,k) = cfl_work(i,j,k) + cfl_inc(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_work(i,j,k)  = T_work(i,j,k) - T_n(i,j,k)
              T_inc(i,j,k)   = T_work(i,j,k) + T_inc(i,j,k)
            End Do
          End Do
        End Do
!
! Call timer for turbulence code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('PC2 Turbulence',4)
!
! ----------------------------------------------------------------------
! Section 2. Output Diagnostics
! ----------------------------------------------------------------------
!
! In the version 5.4 implementation of the PC2 code the diagnostic
! output is not implemented which is why the following code is commented
! out.
!
! Check that turbulence diagnostics requested this timestep
!#if !defined(SCMA)
!      If (sf(0,4)) Then
!        If (Ltimer) Call timer ('Diags   ',3)
!        Call diagnostics_pc2turb(
!     &                       row_length, rows, model_levels
!     &,                      wet_model_levels
!     &,                      me, timestep, at_extremity
!     &,                      T_work,q_work,qcl_work,cf_work,cfl_work
!     &,
!#include <argsts/argsts.h>
!     & STASHwork4
!     &       )
!
!        If (Ltimer) Call timer ('Diags   ',4)
!      End If   ! on error_code and sf(0,4)
!#endif
#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_PC2)) Then

!       Stash 4,161
! DEPENDS ON: scmoutput
        Call SCMoutput(T_work,                                          &
             'dt_pc2turb','Temperature increment PC2 turbulence','K',   &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 4,162
! DEPENDS ON: scmoutput
        Call SCMoutput(q_work,                                          &
             'dq_pc2turb',                                              &
             'Specific humidity increment PC2 turbulence','kg/kg',      &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 4,163
! DEPENDS ON: scmoutput
        Call SCMoutput(qcl_work,                                        &
             'dqcl_pc2turb','QCL increment PC2 turbulence','kg/kg',     &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 4,172
! DEPENDS ON: scmoutput
        Call SCMoutput(cf_work,                                         &
             'dbcf_pc2turb',                                            &
             'Bulk cloud fraction increment PC2 turbulence','Fraction', &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 4,173
! DEPENDS ON: scmoutput
        Call SCMoutput(cfl_work,                                        &
             'dcfl_pc2turb',                                            &
             'Liquid cloud fraction increment PC2 turbulence',          &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_PC2)
#endif
! end of routine turbulence_ctl
      Return
      END SUBROUTINE pc2_turbulence_ctl
#endif
#endif
