
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ A module of variables specifically for mpp usage

Module Rcf_Parvars_Mod

! Description:
!   Variables and parameters for both the decomposition (current)
!   and multi-pe addressing
!
! Derived from UM4.5 code
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   06/01/03   River routing support. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

INTEGER   Ndim_max        ! maximum number of spatial dimensions
PARAMETER (Ndim_max = 3 ) ! 3d data


INTEGER &
&   fld_type_p            &! indicates a grid on P points
&,  fld_type_u            &! indicates a grid on U points
&,  fld_type_v            &! indicates a grid on V points
&,  fld_type_r            &! indicates a river-routing grid
&,  fld_type_unknown     ! indicates a non-standard grid.

PARAMETER ( &
&   fld_type_p=1 &
&,  fld_type_u=2 &
&,  fld_type_v=3 &
&,  fld_type_r=4 &
&,  fld_type_unknown=-1)

INTEGER &
&   local_data &
&,  global_dump_data

PARAMETER ( &
&   local_data=1         &! Used in addressing to indicate if
&,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! For halo types
Integer, Parameter :: halo_type_single   = 1
Integer, Parameter :: halo_type_extended = 2
Integer, Parameter :: halo_type_no_halo  = 3

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================

INTEGER &
&   Maxproc              ! Max number of processors
PARAMETER ( &
&   MAXPROC = 256)

INTEGER &
&   PNorth        &! North processor address in the neighbour array
&,  PEast         &! East  processor address in the neighbour array
&,  PSouth        &! South processor address in the neighbour array
&,  PWest         &! West  processor address in the neighbour array
&,  NoDomain       ! Value in neighbour array if the domain has
                   !  no neighbor in this direction. Otherwise
                   !  the value will be the tid of the neighbor

PARAMETER ( &
&   PNorth   = 1 &
&,  PEast    = 2 &
&,  PSouth   = 3 &
&,  PWest    = 4 &
&,  NoDomain = -1)

INTEGER &
&   BC_STATIC             &! Static boundary conditions
&,  BC_CYCLIC            ! Cyclic boundary conditions

PARAMETER ( &
&   BC_STATIC = 1 &
&,  BC_CYCLIC = 2)


! =======================================================

INTEGER, SAVE :: &
&   first_comp_pe        &! top left pe in LPG
&,  last_comp_pe         &! bottom right pe in LPG
&,  current_decomp_type  &! current decomposition type
&,  Offx                 &! halo size in EW direction
&,  Offy                 &! halo size in NS direction
&,  glsize(Ndim_max)     &! global data size
&,  glsizeu(Ndim_max)    &! global u data size
&,  glsizev(Ndim_max)    &! global v data size
&,  glsizer(Ndim_max)    &! global river-routing data size
&,  lasize(Ndim_max)     &! local data size
&,  blsizep(Ndim_max)    &! personal p data area
&,  blsizeu(Ndim_max)    &! personal u data area
&,  blsizev(Ndim_max)    &! personal v data area
&,  blsizer(Ndim_max)    &! personal river-routing data area
&,  datastart(Ndim_max)  &! position of personal data in global data
                          !   (in terms of standard Fortran array
                          !    notation)
&,  datastartr(Ndim_max) &! position of personal data in global data
                          !  for river routing data
&,  gridsize(Ndim_max)   &! size of the LPG in each dimension
&,  gridpos(Ndim_max)     ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

DATA current_decomp_type/-1/  ! set the initial decomposition
!                                   ! to an "unset" value

LOGICAL, SAVE :: &
&    atSouth             &! process at the bottom of the LPG
&,   atNorth             &! process at the top of the LPG
&,   atWest              &! process at the left of the LPG
&,   atEast               ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

! =======================================================
! Data for Message Passing Software
! =======================================================

INTEGER, SAVE :: &
&  bound(Ndim_max)            &! type of boundary (cyclic or static)
                               !  in each direction
&, g_lasize(Ndim_max,0:maxproc) &
!                                 ! global copy of local data size
&, g_blsizep(Ndim_max,0:maxproc) &
!                                 ! global copy of personal p data area
&, g_blsizeu(Ndim_max,0:maxproc) &
!                                 ! global copy of personal u data area
&, g_blsizev(Ndim_max,0:maxproc) &
!                                 ! global copy of personal v data area
&, g_blsizer(Ndim_max,0:maxproc) &
!                                 ! global copy of personal rr data area
&, g_datastart(Ndim_max,0:maxproc) &
!                                 ! global copy of datastart
&, g_datastartr(Ndim_max,0:maxproc) &
!                                 ! global copy of datastartr
&, g_gridpos(Ndim_max,0:maxproc) &
!                                 ! global copy of gridpos
&, nproc                      &! number of processors in current
!                                 ! decomposition
&, nproc_max                  &! maximum number of processors
&, nproc_x                    &! number of processors in x-direction
&, nproc_y                    &! number of processors in y-direction
&, mype                       &! number of this processor
                               !  (starting from 0)
&, neighbour(4)               &! array with the tids of the four
                               ! neighbours in the horizontal plane
&, gc_proc_row_group          &! GID for procs along a proc row
&, gc_proc_col_group          &! GID for procs along a proc col
&, gc_all_proc_group         ! GID for all procs


End Module Rcf_Parvars_Mod
