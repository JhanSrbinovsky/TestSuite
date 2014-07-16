
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta
      subroutine perturb_theta(                                         &
     &                         theta, row_length, rows, model_levels,   &
     &                         global_row_length, global_rows,          &
     &                         model_domain, at_extremity,              &
     &                         offx, offy, IntRand_Seed, l_datastart    &
     &                        )

! Purpose:
!          Perturb theta at the bit level using a random number  
!
! Method:
!          Is described in ;
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.
      Integer, Intent(In) :: ROW_LENGTH     ! No of points per local row
      Integer, Intent(In) :: ROWS           ! No of local (theta) rows
      Integer, Intent(In) :: MODEL_LEVELS   ! No of model levels
      Integer, Intent(In) :: Offx    ! standard halo size in East-West
      Integer, Intent(In) :: Offy    ! standard halo size in North-South
      Logical, Intent(In) :: At_extremity(4)
      Integer, Intent(In) :: model_domain
      Integer :: IntRand_Seed
      Integer, Intent(In) :: l_datastart(2)
      Integer, Intent(In) :: global_ROW_LENGTH    
      Integer, Intent(In) :: global_ROWs

      Real, Intent (InOut) ::                                           &
     &  theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels) 

! Local Variables
      Integer :: i,j,k    ! loop variables
      Integer :: j_start,j_end
      REAL :: RandomNumber(global_row_length,global_rows,model_levels)
      REAL :: eps_machine
      Integer :: gi,gj

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
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end

!----------------------------------------------------------------------

      if(model_domain==mt_global .and. at_extremity(Psouth))then
        j_start=2
      else 
        j_start=1
      endif
      if(model_domain==mt_global .and. at_extremity(PNorth))then
        j_end=rows-1
      else 
        j_end=rows
      endif

      write(6,*)'PERTURBING THETA FIELD'
      eps_machine=epsilon(1.0)
      write(6,*)' MACHINE epsilon ', eps_machine

! get random number in range 0-1:
! DEPENDS ON: var_randomNumber
      call var_randomNumber(IntRand_seed,RandomNumber,global_row_length,&
     &                      global_rows,model_levels)

      Do k=1,model_levels
        Do j=j_start,j_end
          gj=j+l_datastart(2)-1
          Do i=1, row_length
            gi=i+l_datastart(1)-1
            theta(i,j,k) = theta(i,j,k) +                               &
     &           theta(i,j,k) * (2.0*Randomnumber(gi,gj,k) - 1.0) *     &
     &             eps_machine
          End Do
        End Do
      End Do

      RETURN
      END SUBROUTINE perturb_theta

