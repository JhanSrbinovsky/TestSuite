
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine COPYDIAG ----------------------------------------------
!LL
!LL Purpose : To copy a single diagnostic field from secondary space to
!LL           the main data array for stash processing, and to extend
!LL           the data to a full horizontal field.
!LL
!LL Service routine
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  4.3     10/02/97 mpp code : Added PPX arguments and modified
!LL                   updating of polar rows.            P.Burton
!LL  5.0     30/07/99 Re-write to comply to 5.0 mpp decomposition
!LL                   (In and out arrays can have different halos
!LL                   sizes, but they are not copied accross)
!LL                                                      J-C Thil
!LL  5.2     31/01/01 Correction for LAM borders and bypass populating
!                     halos for normal case. R Rawlins
!    5.5    17/04/03   Remove reference to obsolete section
!                      C90_1A. T.White
!    6.2     03/02/06 Move to c84_1a. P.Selwood.
!LL
!LL Programming Standard :
!LL
!LL Documentation:
!LL
!LLEND -------------------------------------------------------------
!
!*L Arguments
!
      Subroutine copydiag(                                              &
     &     diagout,diagin,                                              &
     &     row_length,rows,                                             &
     &     offx_out, offy_out,                                          &
     &     offx_in, offy_in,                                            &
     &     at_extremity,                                                &
     &     im,is,ie,                                                    &
!#include <argppx/argppx.h>
     &     icode,cmessage)








      Implicit none

      Integer                                                           &
     &     row_length                                                   &
                                ! Number of points in a row
     &,    rows                                                         &
                                ! Number of rows
     &,    offx_out, offy_out                                           &
                                ! offset dimensions of diagout
     &,    offx_in, offy_in                                             &
                                ! offset dimensions of diagin
     &,    im,is,ie                                                     &
                                ! Model, section, item
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Logical                                                           &
     &     at_extremity(4)      ! Indicates if this processor is at
                                !  north, south east or west of the
                                !  processor grid
      Character*80  cmessage

! ARGPPX arguments:
!#include <csubmodl/csubmodl.h>
!#include <cppxref/cppxref.h>
!#include <ppxlook/ppxlook.h>


      Real                                                              &
     & diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in)     &
                                ! Output field
     &,diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out)
                                ! Input field


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
     &   i,j  ! loop bound
!     &,  gr                     ! grid type of grid
!     &,  fld_type               ! is it a p field or a u field?
!     &,  info                   ! GCOM return code

! Functions called:
!      Integer
!     &   exppxi ,get_fld_type

      icode = 0
      cmessage = ""

! Find out the gridtype of the field
!      gr = exppxi(im,is,ie,ppx_grid_type,
!#include <argppx/argppx.h>
!     &            ICODE,CMESSAGE)

!      If (icode  >   0) goto 9999

! and use this to find the field type (p field, u field or v field)

!      fld_type = get_fld_type(gr)


!     Copy fields
      Do j = 1, rows
         Do i = 1, row_length
            diagout(i,j) = diagin(i,j)
         End do
      End do



 9999 continue
      return
      END SUBROUTINE copydiag
