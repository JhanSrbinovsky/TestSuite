#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!      Subroutine PRINTSUB
!
!     Single Column Unified Model routine to write out T and Q, and
!     increments DT and DQ due to statistical forcing.
!     Will print out any number of columns
!     Jennifer Lean 20/10/90
!
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!
!=====================================================================
!
      Subroutine PRINTSUB(                                              &
!     ! IN
     &  points, nlevs, nwet,                                            &
!     !
     &  title, istep, iday, t, q, dt, dq)
!
      Implicit none
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  points                                                          &
                                ! IN no of model columns.
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                   ! IN no of model levels in which Q is
                                !    set

      Character*60                                                      &
     &  title                   ! Heading
      Integer                                                           &
     &  iday                                                            &
                                ! Day number
     &  ,istep                  ! Timestep
      Real                                                              &
     &  q(points,nlevs)                                                 &
                                ! Specific humidity (Kg Kg^-1)
     &  ,t(points,nlevs)                                                &
                                ! Temperature (K)
     &  ,dt(points,nlevs)                                               &
                                ! Temperature increment(K)
     &  ,dq(points,nlevs)       ! Specific humidity increment
                                !  (Kg Kg^-1)
!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
      Character*28                                                      &
     &  cfmt                    ! Format statement for each row
                                !  of variables Q T DT DQ
      Character*33                                                      &
     &  ctfmt                   ! Format statement for title
                                !  of each row
      Integer                                                           &
     &  i                                                               &
                                ! Write statement loop counter
     &  ,l                                                              &
                                ! Loop counter
     &  ,element                                                        &
                                ! Array element no.
     &  ,lastrow                                                        &
                                ! No. of elements in last row
     &  ,nlevsrows, nlevscount                                          &
                                ! No. of rows and Do Loop counter
     &  ,nwetrows, nwetcount    !
!
!     Set format statements
!
      cfmt = '(''         '',  (1pe10.3,1x))'
      ctfmt = '(''0       '',  (3x,''Level'',i2,1x))'
!
!     Loop over sites :
!
      Do l = 1, points
!
!       Heading
!
        Write (11,200) iday, istep, title
!
!       Write out variables T and DT for NLEVS, maximum of 10
!       variables per row
!
        Write (11,201) nlevs
!
!       Calculate no. of rows and no. of elements in last row
!
        If (mod(nlevs,10)  ==  0) then
          nlevsrows = int(nlevs/10)
          lastrow = 10
        else
          nlevsrows = int(nlevs/10) + 1
          lastrow = mod(nlevs,10)
        endif
        Do nlevscount = 1, nlevsrows
          element = 10*(nlevscount-1)
          If (nlevscount  <   nlevsrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write(11,202) (element+i, i = 1, 10)
            Write(11,203) (t(l,element+i), i = 1,10),                   &
     &        (dt(l,element+i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(13:14),'(i2)') lastrow
            Write (11,ctfmt) (element+i,i=1,lastrow)
            Write (cfmt(14:15),'(i2)') lastrow
            Write (cfmt(4:7),'(''T K '')')
            Write (11,cfmt) (t(l,i+element),i=1,lastrow)
            Write (cfmt(4:7),'(''dT K'')')
            Write (11,cfmt) (dt(l,i+element),i=1,lastrow)
          endif
        enddo
!
!       Write out variables Q and DQ for NWET, maximum of 10
!       variables per row
!
        Write (11,204) nwet
!
!       Calculate no. of rows and no. of elements in last row
!
        If ( mod(nwet,10)  ==  0) then
          nwetrows = int(nwet/10)
          lastrow = 10
        else
          nwetrows = int(nwet/10) + 1
          lastrow = mod(nwet,10)
        endif
        Do nwetcount = 1, nwetrows
          element = 10*(nwetcount-1)
          If (nwetcount  <   nwetrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,202) (element+i, i = 1, 10)
            Write (11,205) (q(l,element+i),i = 1, 10),                  &
     &        (dq(l,element+i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(13:14),'(i2)') lastrow
            Write (11,ctfmt) (element+i,i=1,lastrow)
            Write (cfmt(4:10),'(''Q Kg/Kg'')')
            Write (11,cfmt) (q(l,i+element),i=1,lastrow)
            Write (cfmt(4:11),'(''dQ Kg/Kg'')')
            Write (11,cfmt) (dq(l,i+element),i=1,lastrow)
          endif
        enddo
      enddo                     ! l
!
 200  Format('0Day relative to winter solstice',i5,'; timestep ',i5/    &
     &  A60)
 201  Format('0 Number of atmospheric levels = ',i3)
 202  Format('0          level',i2,9(4x,'level',i2))
 203  Format(' t k     ',10(1pe10.3,1x)/,' dT K    ',10(1pe10.3,1x)/)
 204  Format('0 Number of moist atmospheric levels = ',i3)
 205  Format(' Q Kg/Kg ',10(1pe10.3,1x)/,' dQ Kg/Kg',10(1pe10.3,1x)/)
      Return
      END SUBROUTINE PRINTSUB
#endif
