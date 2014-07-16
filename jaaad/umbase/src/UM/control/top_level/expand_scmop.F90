#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Increase the size of the SCM diagnostic arrays in SCMop

      subroutine expand_SCMop(SCMop)
      implicit none

! Description:
      ! Increases the maximum number of diagnostics entries
      ! allowed in SCMop by reallocating the relevant arrays.
! Method:
      ! On the first timestep, as calls to SCMoutput are being made,
      ! memory has to be allocated to arrays to hold the resulting
      ! information. Since it is not known at the start of the run how
      ! many calls to SCMoutput there will be and what their input
      ! parameters are, no memory is allocated at the outset and the
      ! variable SCMop%maxnentries (the maximum no. of diagnostic
      ! "entries" that the arrays in SCMop can handle before they run
      ! out of space) is zero. In this case most of the
      ! statements in this routine are ignored (since they start with
      ! "if (maxnentries >  0)"), and the arrays are simply allocated
      ! with a size equal to "chunk". On subsequent calls the contents
      ! of the allocatable arrays are copied into temporary arrays,
      ! de-allocated, re-allocated with their original size plus
      ! "chunk", and then the data is copied back in. i.e. the arrays
      ! are re-allocated with large sizes without losing their
      ! information.
      ! At the time of writing documentation is available from
      ! http://www-nwp/~frlj
! Owner: Luke Jones
! History:
! Version  Date      Comment
! =======  ====      =======
! 5.5      06/02/03  Original code (Luke Jones)
! 6.0      05/09/03  Adapt to changes in contents of include files.
!                    Luke Jones.
! 6.0      02/09/03  Add chunks of code to expand new SCMop data
!                    members: ncols, nrows and nlevs. Luke Jones.
! 6.0      02/09/03  Removed all reference to diag_mem. Array is
!                    now (re)allocated in SCMoutput. Luke Jones.
! 6.2      19/01/06  Enabled outputting of SCM diagnostics on multiple
!                    sub-steps. Luke Jones.
!
! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

! SCMop_type is defined in here...
#include "s_scmoptype_defn.h"

      type(SCMop_type) :: SCMop ! INOUT The derived-type structure
                                ! containing all the diagnostic
                                ! information

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
#include "s_scmop_internal.h"

      ! Temporary arrays to hold data while SCMop arrays are being
      ! reallocated
      integer, allocatable :: Ixxx(:)
      logical, allocatable :: Lxxx(:)
      character (len=llname), allocatable :: Cxxx(:)
      type(allocatable_array), allocatable :: Dxxx(:)

      ! The number by which to increment SCMop%maxnentries
      integer, parameter :: chunk=100
      ! Holds SCMop%maxnentries
      integer maxnentries

      if (SCMop%nentries /= SCMop%maxnentries) then
         print*,'expand_SCMop WARNING: SCMop is being expanded '//      &
     &        'before it is full, this could be dangerous'
      endif

      maxnentries=SCMop%maxnentries

      ! Allocate the space required for the temporary arrays
      if (maxnentries >  0) allocate(Cxxx(maxnentries))
      if (maxnentries >  0) allocate(Ixxx(maxnentries))
      if (maxnentries >  0) allocate(Lxxx(maxnentries))
      if (maxnentries >  0) allocate(Dxxx(maxnentries))

      ! Increase the size of all the arrays in SCMop associated
      ! to specific diagnostics...

      if (maxnentries >  0) Cxxx=SCMop%sname
      if (maxnentries >  0) deallocate(SCMop%sname)
      allocate(SCMop%sname(maxnentries+chunk))
      if (maxnentries >  0) SCMop%sname(1:maxnentries)=Cxxx

      if (maxnentries >  0) Cxxx=SCMop%lname
      if (maxnentries >  0) deallocate(SCMop%lname)
      allocate(SCMop%lname(maxnentries+chunk))
      if (maxnentries >  0) SCMop%lname(1:maxnentries)=Cxxx

      if (maxnentries >  0) Cxxx=SCMop%units
      if (maxnentries >  0) deallocate(SCMop%units)
      allocate(SCMop%units(maxnentries+chunk))
      if (maxnentries >  0) SCMop%units(1:maxnentries)=Cxxx

      if (maxnentries >  0) Ixxx=SCMop%domprof
      if (maxnentries >  0) deallocate(SCMop%domprof)
      allocate(SCMop%domprof(maxnentries+chunk))
      if (maxnentries >  0) SCMop%domprof(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%timprof
      if (maxnentries >  0) deallocate(SCMop%timprof)
      allocate(SCMop%timprof(maxnentries+chunk))
      if (maxnentries >  0) SCMop%timprof(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%streams
      if (maxnentries >  0) deallocate(SCMop%streams)
      allocate(SCMop%streams(maxnentries+chunk))
      if (maxnentries >  0) SCMop%streams(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%dump_step
      if (maxnentries >  0) deallocate(SCMop%dump_step)
      allocate(SCMop%dump_step(maxnentries+chunk))
      if (maxnentries >  0) SCMop%dump_step(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%nadd2dump
      if (maxnentries >  0) deallocate(SCMop%nadd2dump)
      allocate(SCMop%nadd2dump(maxnentries+chunk))
      if (maxnentries >  0) SCMop%nadd2dump(1:maxnentries)=Ixxx

      if (maxnentries >  0) Lxxx=SCMop%only_radsteps
      if (maxnentries >  0) deallocate(SCMop%only_radsteps)
      allocate(SCMop%only_radsteps(maxnentries+chunk))
      if (maxnentries >  0) SCMop%only_radsteps(1:maxnentries)=Lxxx

      if (maxnentries >  0) Ixxx=SCMop%ncols
      if (maxnentries >  0) deallocate(SCMop%ncols)
      allocate(SCMop%ncols(maxnentries+chunk))
      if (maxnentries >  0) SCMop%ncols(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%nrows
      if (maxnentries >  0) deallocate(SCMop%nrows)
      allocate(SCMop%nrows(maxnentries+chunk))
      if (maxnentries >  0) SCMop%nrows(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%nlevs
      if (maxnentries >  0) deallocate(SCMop%nlevs)
      allocate(SCMop%nlevs(maxnentries+chunk))
      if (maxnentries >  0) SCMop%nlevs(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%nelements
      if (maxnentries >  0) deallocate(SCMop%nelements)
      allocate(SCMop%nelements(maxnentries+chunk))
      if (maxnentries >  0) SCMop%nelements(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%sname_id
      if (maxnentries >  0) deallocate(SCMop%sname_id)
      allocate(SCMop%sname_id(maxnentries+chunk))
      if (maxnentries >  0) SCMop%sname_id(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%wd
      if (maxnentries >  0) deallocate(SCMop%wd)
      allocate(SCMop%wd(maxnentries+chunk))
      if (maxnentries >  0) SCMop%wd(1:maxnentries)=Ixxx

      if (maxnentries >  0) Ixxx=SCMop%lastencounter
      if (maxnentries >  0) deallocate(SCMop%lastencounter)
      allocate(SCMop%lastencounter(maxnentries+chunk))
      if (maxnentries >  0) SCMop%lastencounter(1:maxnentries)=Ixxx

      if (maxnentries > 0) Ixxx=SCMop%substep
      if (maxnentries > 0) deallocate(SCMop%substep)
      allocate(SCMop%substep(maxnentries+chunk))
      if (maxnentries > 0) SCMop%substep(1:maxnentries)=Ixxx

      if (maxnentries >  0) Dxxx=SCMop%diag
      if (maxnentries >  0) deallocate(SCMop%diag)
      allocate(SCMop%diag(maxnentries+chunk))
      if (maxnentries >  0) SCMop%diag(1:maxnentries)=Dxxx

      if (maxnentries >  0) Ixxx=SCMop%netcdf_id
      if (maxnentries >  0) deallocate(SCMop%netcdf_id)
      allocate(SCMop%netcdf_id(maxnentries+chunk))
      if (maxnentries >  0) SCMop%netcdf_id(1:maxnentries)=Ixxx

      if (maxnentries >  0) deallocate(Cxxx)
      if (maxnentries >  0) deallocate(Ixxx)
      if (maxnentries >  0) deallocate(Lxxx)
      if (maxnentries >  0) deallocate(Dxxx)

      SCMop%maxnentries=SCMop%maxnentries+chunk

      return
      END SUBROUTINE expand_SCMop
#endif
