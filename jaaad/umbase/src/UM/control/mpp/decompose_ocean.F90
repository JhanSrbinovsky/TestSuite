#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM: Perform data decomposition for ocean model
!
! Subroutine Interface:
      SUBROUTINE DECOMPOSE_OCEAN(global_row_len,global_n_rows,          &
     &                           tot_levels,model_type,                 &
     &                           nproc_EW, nproc_NS,                    &
     &                           extended_halo_EW,                      &
     &                           extended_halo_NS,                      &
     &                           local_row_len,local_n_rows,            &
     &                           l_ocyclic)
      IMPLICIT NONE
!
! Desciption:
! This routine currently performs a 1D North-South decomposition on the
! ocean model. nproc_EW is currently ignored.
! The decomposition strategy is much the same as the atmosphere's -
! First try and divide the rows equally between processors, and then
! distribute any left over rows to the processors, starting from the
! top.
!
! Method:
! The local data sizes are calculated and sotred in the COMMON block
! DECOMPDB. The boundary conditions are set (cyclic in East/West
! direction if *DEF,GLOBAL
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/8/96  New deck created for MPP ocean model.  P.Burton
!  4.3      17/02/97 Added initialisation of new ocean decompositon
!                    decomp_nowrap_ocean - which does not include
!                    the wrap around points on the ends of rows.
!                    This requires passing in the l_ocyclic variable
!                    to indicate if these points are present.
!                                                         P.Burton
!  5.0      09/08/99 Reflect changes made to Atmosphere version of
!                    this routine in order to cater for MPP control
!                    code changes introduced in connection with the
!                    new atmosphere dynamics developments.
!                                            R. Hill
!  5.1      27/01/99 Added g_pe_index_EW/NS arrays    P.Burton
!  5.1      19/04/00  Correct lasize for velocity grid.
!                     Include set up of lasize values for
!                     halo_type_no_halo conditions. R. Hill
!  5.3      14/09/01 Added model_domain variable    P.Burton
!  5.3               Remove restriction on odd no of
!                    points E-W. That's not an ocean
!                    limitation.   R. Hill
!  5.4    Apr 2002      Fix bug in "extended halo type". Although
!                       not specifically used by the Ocean, these
!                       values are required in certain generic STASH
!                       routines when processing sub-areas.
!                       Also add g_pe_index_EW/NS arrays for non-wrap
!                       grids to re-enable sub-domain STASH extraction
!                       overlooked by 5.1 MPP developments.    R. Hill
!  5.5      25/03/03 Setup of datastart_f (field based datastart).
!                    P.Selwood.
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine arguments:

      INTEGER                                                           &

     &  global_row_len                                                  &
                        ! IN :  number of E-W points of entire model
     &, global_n_rows                                                   &
                        ! IN :  number of N-S points of entire model
     &, tot_levels                                                      &
                        ! IN :  total number of levels
     &, model_type                                                      &
                        ! IN  : type (Global,LAM etc) of model
     &, nproc_EW                                                        &
                        ! IN :  number of processors to decompose E-W
     &, nproc_NS                                                        &
                        ! IN :  number of processors to decompose N-S
     &, extended_halo_EW                                                &
                         ! IN  : size of extended EW halo
     &, extended_halo_NS                                                &
                         ! IN  : size of extended NS halo
     &, local_row_len                                                   &
                        ! OUT : local number of E-W points
     &, local_n_rows    ! OUT : local number of N-S points
!                       ! local_row_len and local_n_rows include
!                       ! any halos

      LOGICAL                                                           &

     &  l_ocyclic       ! IN : true if extra wrap points are present
!                       !      at the start/ends of rows

! Parameters and Common blocks

#include "parvars.h"
#include "gccom.h"
#include "decomptp.h"
#include "decompdb.h"

! Local variables
      INTEGER                                                           &
     &  iproc                                                           &
     &, iproc_x                                                         &
     &, iproc_y                                                         &
     &, ifld                                                            &
     &, ihalo                                                           &
     &, idim                                                            &
     &, ipt                                                             &
     &, irest                                                           &
     &, jrest                                                           &
     &, start                                                           &
     &, size                                                            &
     &, prow_N                                                          &
     &, prow_S                                                          &
     &, info                                                            &
     &, in_ocn_decomp

      LOGICAL                                                           &
     &  at_north

      INTEGER                                                           &
     &  size_x(0:nproc_EW-1)                                            &
     &, size_y(0:nproc_NS-1)                                            &
     &, start_x(0:nproc_EW-1)                                           &
     &, start_y(0:nproc_NS-1)

! Error reporting
      INTEGER       ErrorStatus ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='DECOMPOSE_OCEAN')


! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

      IF (nproc_EW*nproc_NS  >   Maxproc) THEN
        ErrorStatus=1
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with decomposition ",I3," x ",I3,               &
     &      " (",I3,") processors. ",                                   &
     &      "Maxproc is ",I3," processors.")')                          &
     &       nproc_EW,nproc_NS,nproc_EW*nproc_NS,Maxproc
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF ((nproc_EW  /=  1) .AND. (MOD(nproc_EW,2)  /=  0)) THEN
        ErrorStatus=2
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with an odd (",I3,") number of processors ",    &
     &      "in the East-West direction.")') nproc_EW
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

!      No restriction on odd numbers of E-W points

      IF (extended_halo_EW  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("East-West extended halo size (",I2,") is too large.",      &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_EW,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (extended_halo_NS  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("North-South extended halo size (",I2,") is too large.",    &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_NS,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! ------------------------------------------------------------------
      decomp_db_sb_model_domain(decomp_standard_ocean)=model_type

! Halo Sizes
      ! Reminder: The Ocean only decomposes in the S-N direction so
      ! no halos are required E-W.
      decomp_db_halosize(1,halo_type_single,decomp_standard_ocean) = 0
      decomp_db_halosize(2,halo_type_single,decomp_standard_ocean) = 1
      decomp_db_halosize(3,halo_type_single,decomp_standard_ocean) = 0

      ! As far as the ocean is concerned, we dont worry about
      ! extended halo sizes. Just use the standard settings.
      decomp_db_halosize(1,halo_type_extended,decomp_standard_ocean) =  &
     &  0
      decomp_db_halosize(2,halo_type_extended,decomp_standard_ocean) =  &
     &  1
      decomp_db_halosize(3,halo_type_extended,decomp_standard_ocean) =  &
     &  0

      ! As far as the ocean is concerned, no_halo means exactly that
      ! ... set accordingly.
      decomp_db_halosize(1,halo_type_no_halo,decomp_standard_ocean) = 0
      decomp_db_halosize(2,halo_type_no_halo,decomp_standard_ocean) = 0
      decomp_db_halosize(3,halo_type_no_halo,decomp_standard_ocean) = 0

! Size of global data

      decomp_db_glsize(1,fld_type_p,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_p,decomp_standard_ocean) =            &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_p,decomp_standard_ocean) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_u,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_u,decomp_standard_ocean) =            &
     &  global_n_rows-1
      decomp_db_glsize(3,fld_type_u,decomp_standard_ocean) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_v,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_v,decomp_standard_ocean) =            &
     &  global_n_rows-1
      decomp_db_glsize(3,fld_type_v,decomp_standard_ocean) =            &
     &  tot_levels

! Make sure there's actually enough work for all the processors to do

      IF (nproc_NS  >   global_n_rows) THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Warning : Ocean model has more processors than ', &
     &               'rows. Reducing nproc_y to ',global_n_rows
        ENDIF
        nproc_NS=global_n_rows
      ENDIF

      decomp_db_nproc(decomp_standard_ocean)=nproc_NS
      decomp_db_first_comp_pe(decomp_standard_ocean) = 0
      decomp_db_last_comp_pe(decomp_standard_ocean) =                   &
     &  decomp_db_nproc(decomp_standard_ocean)-1

! Set the size of the Logical Processor Grid (LPG)

      decomp_db_gridsize(1,decomp_standard_ocean) = nproc_EW  ! =1
      decomp_db_gridsize(2,decomp_standard_ocean) = nproc_NS
      decomp_db_gridsize(3,decomp_standard_ocean) = 1

! Calculate processor specific information.

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
!       ! Loop over all processors in this decomposition

! NB : Although the decomposition is currently only N-S, all
! the code is included to allow an E-W decomposition too.
! All that is required is to supply nproc_NS > 1.

! Calculate the position in the LPG:
        decomp_db_g_gridpos(3,iproc,decomp_standard_ocean) = 0
        decomp_db_g_gridpos(2,iproc,decomp_standard_ocean) =            &
     &    iproc / decomp_db_gridsize(1,decomp_standard_ocean)
        decomp_db_g_gridpos(1,iproc,decomp_standard_ocean) =            &
     &    iproc - decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*   &
     &            decomp_db_gridsize(1,decomp_standard_ocean)

! Calculate the local data sizes for processor iproc

! East-West decomposition

        decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_p,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_u,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_u,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_v,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_v,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)


        irest = decomp_db_glsize(1,fld_type_p,decomp_standard_ocean)-   &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)*   &
     &          decomp_db_gridsize(1,decomp_standard_ocean)


        decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =          &
     &    decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)*           &
     & decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)+1

        IF (decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)  <       &
     &      irest) THEN

          decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)= &
     &   decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)+1

          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +      &
     &      decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)

        ELSE
          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +      &
     &      irest
        ENDIF

        decomp_db_g_datastart_f(1,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(1,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(1,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)

        decomp_db_g_lasize(1,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean) +  &
     &  2*decomp_db_halosize(1,halo_type_single,decomp_standard_ocean)
        ! U grid, imt values
        decomp_db_g_lasize(1,fld_type_u,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &   decomp_db_g_lasize(1,fld_type_p,halo_type_single,              &
     &                     iproc,decomp_standard_ocean)

        ! V grid, imt values
        decomp_db_g_lasize(1,fld_type_v,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_lasize(1,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)


!-----------------------------------------
! Set up column-wise lasize for no_halo conditions
!-----------------------------------------
        decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean) +  &
     & 2*decomp_db_halosize(1,halo_type_no_halo,decomp_standard_ocean)

        ! U grid, imt values
        decomp_db_g_lasize(1,fld_type_u,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &   decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,             &
     &                     iproc,decomp_standard_ocean)

        ! V grid, imt values
        decomp_db_g_lasize(1,fld_type_v,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)



! North-South decomposition

        decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) =  &
     &    decomp_db_glsize(2,fld_type_p,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(2,decomp_standard_ocean)

        jrest = decomp_db_glsize(2,fld_type_p,decomp_standard_ocean)-   &
     &   decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)*  &
     &          decomp_db_gridsize(2,decomp_standard_ocean)

        decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =          &
     &    decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*           &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)+1

        IF (decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)  <       &
     &      jrest) THEN

          decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)= &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)+1

          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +      &
     &      decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)

        ELSE
          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +      &
     &      jrest
        ENDIF

        decomp_db_g_datastart_f(2,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(2,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(2,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)

        decomp_db_g_lasize(2,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) +  &
     &    2*decomp_db_halosize(2,halo_type_single,decomp_standard_ocean)


        ! U velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_u,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_single,            &
     &                     iproc,decomp_standard_ocean)

        ! V velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_v,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_single,            &
     &                     iproc,decomp_standard_ocean)
!-----------------------------------------
! Set up row-wise lasize for no_halo conditions
!-----------------------------------------

        decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) +  &
     & 2*decomp_db_halosize(2,halo_type_no_halo,decomp_standard_ocean)


        ! U velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_u,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,           &
     &                     iproc,decomp_standard_ocean)

        ! V velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_v,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,           &
     &                     iproc,decomp_standard_ocean)
! No decomposition in the vertical

        decomp_db_g_datastart(3,iproc,decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)=   &
     &    tot_levels
        DO ihalo=1,NHalo_max
           DO ifld=1,Nfld_max
              decomp_db_g_lasize(3,ifld,ihalo,                          &
     &                iproc,decomp_standard_ocean) = tot_levels
           ENDDO
        ENDDO



! One less U/V row at North. lasize does not reflect this - velocity
! grid sizes are the same as tracer grid sizes, row-wise for MPP
! purposes (mainly in STASH). Effectively, lasize on the northern-most
! PE is over-dimensioned by 1 row.

        decomp_db_g_blsize(1,fld_type_u,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_v,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)

        IF (  decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)        &
     &   ==  (decomp_db_gridsize(2,decomp_standard_ocean)-1)) THEN

          decomp_db_g_blsize(2,fld_type_u,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)-1

          decomp_db_g_blsize(2,fld_type_v,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)-1


        ELSE

          decomp_db_g_blsize(2,fld_type_u,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)

          decomp_db_g_blsize(2,fld_type_v,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)


        ENDIF

        decomp_db_g_blsize(3,fld_type_u,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)

        decomp_db_g_blsize(3,fld_type_v,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)

      ENDDO  ! loop over processors

! Set up the pe_index_EW_array - for each point along a global
! row it indicates the PE index (along the processor
! row) which contains that point

      DO iproc_x=0,nproc_EW-1
        DO ipt=decomp_db_g_datastart(1,iproc_x,decomp_standard_ocean),  &
     &         decomp_db_g_datastart(1,iproc_x,decomp_standard_ocean)+  &
     &         decomp_db_g_blsize(1,fld_type_p,iproc_x,                 &
     &                            decomp_standard_ocean)
          decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=iproc_x
          decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=iproc_x
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_EW,1
        decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=0
        decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=0
      ENDDO

      DO ipt=global_row_len+1,global_row_len+1+extended_halo_EW
        decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=nproc_x-1
        decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=nproc_x-1
      ENDDO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

      DO iproc_y=0,nproc_NS-1
        DO ipt=decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_ocean),            &
     &         decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_ocean)+            &
     &         decomp_db_g_blsize(2,fld_type_p,iproc_y*nproc_EW,        &
     &                            decomp_standard_ocean)
          decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=iproc_y
          decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=iproc_y
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_NS,1
        decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=0
        decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=0
      ENDDO

      DO ipt=global_n_rows+1,global_n_rows+1+extended_halo_NS
        decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=nproc_y-1
        decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=nproc_y-1
      ENDDO
! Set up the boundary types

#if defined(GLOBAL)
      decomp_db_bound(1,decomp_standard_ocean) = BC_CYCLIC
!       ! Cyclic East-West boundaries
#else
      decomp_db_bound(1,decomp_standard_ocean) = BC_STATIC
!       ! No East-West wrap around
#endif
      decomp_db_bound(2,decomp_standard_ocean) = BC_STATIC
!       ! No North-South wrap around
      decomp_db_bound(3,decomp_standard_ocean) = BC_STATIC
!       ! No vertical wrap around

! And set up the neighbour array

! DEPENDS ON: set_neighbour
      CALL SET_NEIGHBOUR(                                               &
     &  decomp_standard_ocean)

! Set up the GCOM groups

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_standard_ocean)  ==  1)          &
     & THEN
       decomp_db_gc_proc_row_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(2,mype,decomp_standard_ocean),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_row_group(decomp_standard_ocean))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_standard_ocean)  ==  1)          &
     & THEN
        decomp_db_gc_proc_col_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(1,mype,decomp_standard_ocean),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_col_group(decomp_standard_ocean))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_standard_ocean)  ==  nproc_max)        &
     & THEN
        decomp_db_gc_all_proc_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        IF ((mype  >=  decomp_db_first_comp_pe(decomp_standard_ocean))  &
     &    .AND.                                                         &
     &     (mype  <=  decomp_db_last_comp_pe(decomp_standard_ocean) ))  &
     &  THEN
          in_ocn_decomp=1
        ELSE
          in_ocn_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_ocn_decomp,info,               &
     &    decomp_db_gc_all_proc_group(decomp_standard_ocean))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_standard_ocean)=.TRUE.

! Initialise decomp_nowrap_ocean which doesn't contain extra wrap
! points at start and end of each row
! Mostly it is a straight copy of the original ocean decomposition

      decomp_db_sb_model_domain(decomp_nowrap_ocean)=                   &
     &  decomp_db_sb_model_domain(decomp_standard_ocean)
      DO idim=1,Ndim_max
        decomp_db_bound(idim,decomp_nowrap_ocean)=                      &
     &    decomp_db_bound(idim,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_p,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_p,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_u,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_u,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_v,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_v,decomp_standard_ocean)


        decomp_db_gridsize(idim,decomp_nowrap_ocean)=                   &
     &    decomp_db_gridsize(idim,decomp_standard_ocean)
        decomp_db_halosize(idim,halo_type_single,decomp_nowrap_ocean)=  &
     &  decomp_db_halosize(idim,halo_type_single,decomp_standard_ocean)
      ENDDO

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
         DO ihalo=1,NHalo_max
            DO ifld=1,Nfld_max
               DO idim=1,Ndim_max
                  decomp_db_g_lasize(idim,ifld,ihalo                    &
     &                       ,iproc,decomp_nowrap_ocean)=               &
     &                    decomp_db_g_lasize(idim,ifld,ihalo            &
     &                         ,iproc,decomp_standard_ocean)
               ENDDO
            ENDDO
         ENDDO

         DO ifld=1,Nfld_max
            DO idim=1,Ndim_max
             decomp_db_g_blsize(idim,ifld,iproc,decomp_nowrap_ocean)=   &
     &       decomp_db_g_blsize(idim,ifld,iproc,decomp_standard_ocean)
            ENDDO
         ENDDO


         DO idim=1,Ndim_max
            decomp_db_g_datastart(idim,iproc,decomp_nowrap_ocean)=      &
     &      decomp_db_g_datastart(idim,iproc,decomp_standard_ocean)

           decomp_db_g_datastart_f(idim,fld_type_p,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_p,iproc,              &
     &                             decomp_standard_ocean)
           decomp_db_g_datastart_f(idim,fld_type_u,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_u,iproc,              &
     &                             decomp_standard_ocean)
           decomp_db_g_datastart_f(idim,fld_type_v,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_v,iproc,              &
     &                             decomp_standard_ocean)

            decomp_db_g_gridpos(idim,iproc,decomp_nowrap_ocean)=        &
     &      decomp_db_g_gridpos(idim,iproc,decomp_standard_ocean)
         ENDDO
      ENDDO

      DO idim=1,4
        decomp_db_neighbour(idim,decomp_nowrap_ocean)=                  &
     &    decomp_db_neighbour(idim,decomp_standard_ocean)
      ENDDO

      decomp_db_first_comp_pe(decomp_nowrap_ocean)=                     &
     &  decomp_db_first_comp_pe(decomp_standard_ocean)
      decomp_db_last_comp_pe(decomp_nowrap_ocean)=                      &
     &  decomp_db_last_comp_pe(decomp_standard_ocean)
      decomp_db_nproc(decomp_nowrap_ocean)=                             &
     &  decomp_db_nproc(decomp_standard_ocean)
      decomp_db_gc_proc_row_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_proc_row_group(decomp_standard_ocean)
      decomp_db_gc_proc_col_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_proc_col_group(decomp_standard_ocean)
      decomp_db_gc_all_proc_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_all_proc_group(decomp_standard_ocean)

      IF (l_ocyclic) THEN
! Make modifications to the decompositions to remove the point at
! the beginning and end of each row
        decomp_db_glsize(1,fld_type_p,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_p,decomp_nowrap_ocean)-2

        decomp_db_glsize(1,fld_type_u,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_u,decomp_nowrap_ocean)-2

        decomp_db_glsize(1,fld_type_v,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_v,decomp_nowrap_ocean)-2

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)          &
     &         ==  0) THEN  ! this processor at left of LPG

            DO ihalo=1,Nhalo_max
               decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
               decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
               decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
            ENDDO


            decomp_db_g_blsize(1,fld_type_p,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_p,iproc,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_u,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_u,iproc,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_v,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_v,iproc,decomp_nowrap_ocean)-1


          ELSE  ! processor not at left of LPG

            decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)=         &
     &        decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_datastart_f(1,fld_type_p,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_p,iproc,                &
     &                               decomp_nowrap_ocean) - 1
            decomp_db_g_datastart_f(1,fld_type_u,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_u,iproc,                &
     &                               decomp_nowrap_ocean) - 1
            decomp_db_g_datastart_f(1,fld_type_v,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_v,iproc,                &
     &                               decomp_nowrap_ocean) - 1

          ENDIF

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)          &
     &        ==  (decomp_db_gridsize(1,decomp_nowrap_ocean)-1) )       &
     &    THEN  ! this processor at right of LPG

            DO ihalo=1,Nhalo_max
               decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1

               decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1

               decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1
            ENDDO

            decomp_db_g_blsize(1,fld_type_p,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_p,iproc             &
     &                                       ,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_u,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_u,iproc             &
     &                              ,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_v,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_v,iproc             &
     &                              ,decomp_nowrap_ocean)-1

          ENDIF

        ENDDO

      ENDIF

! Finally, indicate this decomposition has been initialised

      decomp_db_set(decomp_nowrap_ocean)=.TRUE.
! And return the new horizontal dimensions


       local_row_len=decomp_db_g_lasize(1,fld_type_p,halo_type_single,  &
     &                     mype,decomp_standard_ocean)

       local_n_rows=decomp_db_g_lasize(2,fld_type_p,halo_type_single,   &
     &                     mype,decomp_standard_ocean)

      RETURN
      END SUBROUTINE DECOMPOSE_OCEAN

#endif
#endif
