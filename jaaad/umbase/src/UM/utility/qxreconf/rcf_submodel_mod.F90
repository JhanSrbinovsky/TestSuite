#if defined(RECON) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! A module to contain information about submodels

Module Rcf_Submodel_Mod

! Description:
!   Data module to contain information about submodels. Largely
!   unnecessary but used for consistency with UM.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

!  1. Maximum internal model/submodel array sizes for this version.

! Max no. of internal models (1 more for expansion according to csubmmax.h)
Integer, Parameter  ::  N_Internal_Model_Max = 2

! Max no. of subm. dump parts (1 more for expansion according to csubmmax.h)
Integer, Parameter  ::  N_Submodel_Partition_Max = 2

! Max value of int. model id
Integer, Parameter  ::  Internal_Id_Max = N_Internal_Model_Max

! Max value of subm. dump id
Integer, Parameter  ::  Submodel_Id_Max = N_Submodel_Partition_Max


!  2. Internal Model identifiers in Long and Short form
Integer, Parameter  ::   Atmos_IM = 1,   & ! Atmosphere IM
                         A_IM     = 1

Integer, Parameter ::    Atmos_SM = 1,   & ! Atmos SM partition
                         A_SM     = 1

!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.

Integer, Save       ::         &
      N_Internal_Model,        & ! No. of internal models
      N_Submodel_Partition,    & ! No. of submodel partitions
      Internal_Model_List(N_Internal_Model_Max), & ! Internal models
      Submodel_For_IM    (N_Internal_Model_Max), & ! Submodel identifier
                                ! for each internal model in list
      Submodel_For_SM(N_Internal_Model_Max) ! Submodel number for
                                            ! each submodel id

! Namelist for information in 3.

Namelist/NSUBMODL/ N_Internal_Model, N_Submodel_Partition,     &
     Internal_Model_List, Submodel_For_IM

!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.

Integer, Save    ::                     &

    ! No of internal models in  each submodel partition
    ! indexed by sm identifier
    N_Internal_For_SM(Submodel_ID_Max), &

    ! List of submodel partition identifiers
    Submodel_Partition_List(N_Submodel_Partition_Max) , &

    ! Submodel partition identifier indexed by internal model identifier
    Submodel_Partition_Index(Internal_ID_Max) , &

    ! Sequence number of internal model indexed by internal model
    ! identifier: required to map from id to STASH internal model
    !sequence
    Internal_Model_Index(Internal_ID_Max)

Logical, Save   :: &
     Last_IM_IN_SM(Internal_ID_Max)   ! Last internal model within
                                      ! a submodel partition if .TRUE.,
                                      ! indexed by internal model id.

Integer, Save  ::  Submodel_Ident

End Module Rcf_Submodel_Mod
#endif
