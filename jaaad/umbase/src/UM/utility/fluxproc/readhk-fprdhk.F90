#if defined(FLUXPROC) || defined(FLXPLPR)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
!
! Author:     J. O. S. Alves
!----------------------------------------------------------------------
! contains routines: readhk
!
! Purpose: Flux processing routine.
!          Reads housekeeping file
!----------------------------------------------------------------------
!*********************************************************************
!**  SUBROUTINE: READHK
!**
!**  AUTHOR:  Oscar Alves                 DATE: 05/11/93
!**
!**  PURPOSE: Reads the date and time from the house keeping file
!**
!*********************************************************************
        subroutine readhk(iunit,hour,day,month,year,icode)

        implicit none

! argument list: intent IN
        integer iunit   ! unit to read

!              intent OUT
        integer hour,day,month,year  ! house keeping time
        integer icode                ! error code

        read (iunit,901,iostat=icode) year,month,day,hour

        if (icode  /=  0) then
          icode = 7
        endif

901     format (i4,3i2)

        return
        END SUBROUTINE readhk
!----------------------------------------------------------------------
#endif
