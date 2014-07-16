
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine COPYDIAG_3D -------------------------------------------
!LL
!LL Purpose : To copy a diagnostic field from secondary space to the
!LL           main data array for stash processing, and to extend the
!LL           data to a full horizontal field. Input data of multilevel
!LL           fields is assumed to be on all model levels. Output data
!LL           is on the levels required.
!LL Service routine
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  4.3     10/02/97 mpp code : Added PPX arguments and modified
!LL                   updating of polar rows.            P.Burton
!LL  4.5     23/10/98  Introduce Single Column Model. JC Thil
!LL  5.0     30/07/99 Re-write to comply to 5.0 mpp decomposition
!LL                   (In and out arrays can have different halos
!LL                   sizes, but they are not copied accross)
!LL                                                      J-C Thil
!LL  5.2     31/01/01 Correction for LAM borders and bypass populating
!                     halos for normal case. R Rawlins
!    5.5    17/04/03   Remove reference to obsolete section
!                      C90_1A. T.White
!    6.2     23/06/06 Allow use for SCM to resolve external.
!                     P.Selwood.
!    6.2     03/02/06 Move to c84_1a. P.Selwood
!LL
!LL Programming Standard :
!LL
!LL Documentation:
!LL
!LLEND --------------------------------------------------------------

!*L Arguments

      Subroutine copydiag_3d(                                           &
     &     diagout,diagin,                                              &
     &     row_length,rows,levels,                                      &
     &     offx_out, offy_out,                                          &
     &     offx_in, offy_in,                                            &
     &     at_extremity,                                                &
     &     stlist,len_stlist,stash_levels,                              &
     &     len_stashlevels,                                             &
     &     im,is,ie,                                                    &
!#include <argppx/argppx.h>
     &     icode,cmessage)








      Implicit None

      Integer                                                           &
     &       row_length,                                                &
                          ! Number of points in a row
     &       rows,                                                      &
                          ! Number of rows
     &       offx_out, offy_out,                                        &
                                 ! offset dimensions of diagout
     &       offx_in, offy_in,                                          &
                                 ! offset dimensions of diagin
     &       levels,                                                    &
                          ! Number of levels in input data
     &       len_stlist,                                                &
                          !
     &       stlist(len_stlist),                                        &
                                 ! Stash list
     &       len_stashlevels,                                           &
                              !
     &       stash_levels(len_stashlevels,*),                           &
                                              ! Stash levels list.
     &       im,is,ie,                                                  &
                          ! Model, section, item
     &       icode        ! Return code =0 Normal exit
!                                       >1 Error message
      LOGICAL    at_extremity(4)
!! logicals indicating if a processor
                             ! is at the edge of the LPG
! argppx arguments:
!#include <csubmodl/csubmodl.h>
!#include <cppxref/cppxref.h>
!#include <ppxlook/ppxlook.h>

      Character*80 cmessage

      Real                                                              &
     &   diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in    &
     &          ,levels)                                                &
                                ! Output data
     &,  diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out&
     &          ,*)             ! Input data

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

!     Local variables

      Integer                                                           &
     &   i,j,                                                           &
                               !
     &   k,                                                             &
                               !
     &   kout                 !
!     &,  gr                   ! grid type of grid
!     &,  fld_type             ! is it a p field or a u field?
!     &,  info                 ! GCOM return code

! Functions called:
!      Integer
!     &   exppxi,get_fld_type


!      Real
!     &   copy_value_start(levels)  ! value to copy into start of array
!     &,  copy_value_end(levels)    ! value to copy into end of array


      Logical                                                           &
     &   list(levels) !

      icode = 0
      cmessage = ""

! Find out the gridtype of the field
!      gr = exppxi(im,is,ie,ppx_grid_type,
!#include <argppx/argppx.h>
!     &            icode,cmessage)

!      If (icode  >   0) goto 9999

! and use this to find the field type (p field or u field)

!      fld_type=get_fld_type(gr)


! Find the values to copy into the start and end of the arrays

!      Do k=1,levels
!        copy_value_start(k)=diagin(start_point,k)
!        copy_value_end(k)=diagin(end_point,k)
!      Enddo

! If this is the Northern processor row - we must make sure we
! get a consistent value over the polar row - so we take
! the value from PE 0 and use that on all processors
!      If (attop) then
!        Call gcg_rbcast(700,levels,first_comp_pe,gc_proc_row_group,
!     &                  info,copy_value_start)
!      Endif

! If this is the Southern processor row - we must make sure we
! get a consistent value over the polar row - so we take
! the value from PE 0 and use that on all processors
!      If (atbase) then
!        Call gcg_rbcast(701,levels,last_comp_pe,gc_proc_row_group,
!     &                  info,copy_value_end)
!      Endif

! DEPENDS ON: set_levels_list
      Call set_levels_list(levels,len_stlist,stlist,list,stash_levels,  &
     &      len_stashlevels,icode,cmessage)
      If(icode >  0) goto 9999

!L Move data from DIAGIN to DIAGOUT at levels requested

      kout = 0
      Do k = 1, levels
         If (list(k)) then
            kout = kout + 1

!           Copy fields
            Do j = 1, rows
               Do i = 1, row_length
                  diagout(i,j,kout) = diagin(i,j,k)
               End do
            End do



         end if
      end do


 9999 continue
      return
      END SUBROUTINE copydiag_3d
