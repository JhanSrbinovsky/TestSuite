#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine INITSTAT
! Purpose:-           To calculate the initial variables required by
!                     statistical forcing routines used later
!                     and also prints out initial climate datasets
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!  5.3     02/05/01   Change to Charney-Philips grid and hence change
!                     to variables required by the forcing ie p and
!                     rho instead of pstar         Z. Gardner
!  5.4     05/09/02   Corrections to format statements - M.Hughes
!  5.5     06/02/03   Runtime determination of some format statements.
!                     Luke Jones.
!=====================================================================
!
      Subroutine INITSTAT(                                              &
     &  row_length, rows, nlevs, nwet, ntrop,                           &
     &  andayy, dayno, q, t, lat, long,                                 &
     &  p_in, pa, pb, alfada, alfadb, tbara, tbarb,                     &
     &  tsda, tsdb, tgrada, tgradb, dbara, dbarb, dgrada, dgradb,       &
     &  vnbara, vnbarb, vnsda, vnsdb, vpbara, vpbarb, wbara, wbarb,     &
     &  wsda, wsdb, atime, btime, p_theta_levels)

      Implicit none

      Integer                                                           &
     &  row_length                                                      &
                                ! IN x direction dimension
     &  ,rows                                                           &
                                ! IN y direction dimension
     &  ,nlevs                                                          &
                                ! IN no of levels
     &  ,nwet                                                           &
                                ! IN Number of model levels in which
                                !    Q is set.
     &  ,ntrop                  ! IN Max number of levels in the
                                !    troposphere

#include "c_pi.h"
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Real                                                              &
     &  andayy                  ! IN No. of days in 1 year
      Integer                                                           &
     &  dayno                   ! IN  Day number relative to winter
                                !    solstice
      Real                                                              &
     &  alfada(row_length, rows)                                        &
                                ! OUT Amplitude and mean of
     &  ,alfadb(row_length, rows)                                       &
                                 !    seasonal variation of tuning
                                !    hybrid vertical coordinate
     &  ,atime,btime                                                    &
                                ! OUT Constants for calculating
                                !    annual cycle used in eqn 2.33
                                !    in SCM doc.
     &  ,dbara(row_length, rows,nwet)                                   &
                                ! OUT Amplitude and mean of seasonal
     &  ,dbarb(row_length, rows,nwet)                                   &
                                !    variation of mean dew pt.
                                !    depression (K)
     &  ,dgrada(row_length, rows,nwet)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,dgradb(row_length, rows,nwet)                                  &
                                !    variation of gradient of
                                !    dew pt. depression (K/km)
     &  ,lat0                                                           &
                                ! Dummy for I/Os
     &  ,lat(row_length, rows)                                          &
                                ! OUT Latitude and longitude of
     &  ,long(row_length, rows)                                         &
                                !    gridpoint
     &  ,long0                                                          &
                                ! Dummy for I/Os
     &  ,pa(row_length, rows, nlevs+1)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,pb(row_length, rows, nlevs+1)                                  &
                                !    variation of pressure
     &  ,q(row_length, rows,nwet)                                       &
                                 ! INOUT Specific humidity (Kg Kg**-1)
     &  ,p_theta_levels(row_length, rows,nlevs)                         &
                                                  ! IN pressure
                                !  ((HPa or mb))
     &  ,t(row_length, rows,nlevs)                                      &
                                  ! INOUT Temps(K)
     &  ,tbara(row_length, rows,nlevs)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,tbarb(row_length, rows,nlevs)                                  &
                                !    variation of mean temo. (K)
     &  ,tgrada(row_length, rows,nlevs)                                 &
                                ! OUT Amplitude and mean of seasonal
     &  ,tgradb(row_length, rows,nlevs)                                 &
                                !    variation of temp. gradient
                                !    (K km**-1)
     &  ,tsda(row_length, rows,nlevs)                                   &
                                ! OUT Amplitude and mean of seasonal
     &  ,tsdb(row_length, rows,nlevs)                                   &
                                !    variation of SD of temp. (K)
     &  ,tstara(row_length, rows)                                       &
                                ! OUT Amplitude and mean of seasonal
     &  ,tstarb(row_length, rows)                                       &
                                !    variation of surface temp. (K)
     &  ,vnbara(row_length, rows,nlevs)                                 &
                                ! OUT Amplitude and mean of seasonal
     &  ,vnbarb(row_length, rows,nlevs)                                 &
                                !    variation of velocity VN
                                !    (m s**-1)
     &  ,vnsda(row_length, rows,nlevs)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,vnsdb(row_length, rows,nlevs)                                  &
                                !    variation of SD of velocity VN
                                !    (m s**-1)
     &  ,vpbara(row_length, rows,nlevs)                                 &
                                ! OUT Amplitude and mean of seasonal
     &  ,vpbarb(row_length, rows,nlevs)                                 &
                                !    variation of velocity VP
                                !    (m s**-1)
     &  ,wbara(row_length, rows,ntrop)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,wbarb(row_length, rows,ntrop)                                  &
                                !    variation of vert. vel.
                                !    ( mb s**-1)
     &  ,wsda(row_length, rows,ntrop)                                   &
                                ! OUT Amplitude and mean of seasonal
     &  ,wsdb(row_length, rows,ntrop)
                                !    variation of SD of vert. vel.
                                !    (mb s**-1)
!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
!     Variables for printout of climate dataset
!
      Character*29                                                      &
     &  cfmt                    ! Format statement for each row
                                !  of variables
      Character*34                                                      &
     &  ctfmt                   ! Format statement for title
                                !  of each row
      Integer                                                           &
     &  element                                                         &
                                ! Array element no.
     &  ,lastrow                                                        &
                                ! No. of elements in last row
     &  ,nlevsrows, nlevscount                                          &
                                ! No. of rows and Do Loop counter
     &  ,ntroprows, ntropcount                                          &
                                ! elements
     &  ,nwetrows, nwetcount    ! elements
!
      Integer                                                           &
     &  i, j, k, l, m              ! Loop counter
      Real                                                              &
     &  alfad0                                                          &
                                ! Dummy for I/Os
     &  ,alfad1(row_length, rows)                                       &
     &  ,alfad2(row_length, rows)                                       &
                                 ! Tuning factor for Jan and July
     &  ,daynew                                                         &
                                ! Function to calculate SIN arg
                                ! (eqn. 2.33 in SCM doc.)
     &  ,daysol1,daysol2                                                &
                                ! No. of days after winter
     &  ,dbar0(nwet)                                                    &
                                ! dummy var for I/Os
     &  ,dbar1(row_length, rows,nwet)                                   &
                                       ! Mean dew pt. depressions
     &  ,dbar2(row_length, rows,nwet)                                   &
                                       !  for Jan. and July (K)
     &  ,dewpt(row_length, rows,nwet,2)                                 &
                                          ! Dew point (K)
     &  ,dgrad0(nwet)                                                   &
                                ! Dummy var for I/Os
     &  ,dgrad1(row_length, rows,nwet)                                  &
                                ! Gradient dew pt. depressions
     &  ,dgrad2(row_length, rows,nwet)                                  &
                                !  for Jan. and July (K km**-1)
     &  ,p0(nlevs+1)                                                    &
                            ! Dummy for I/Os
     &  ,p1(row_length, rows, nlevs+1)                                  &
                                      ! pressure for
     &  ,p2(row_length, rows, nlevs+1)                                  &
                                      !  Jan. and July
     &  ,p_in(row_length, rows, nlevs+1)                                &
                                        ! Initial pressure
     &  ,qi(row_length, rows,nwet)                                      &
                                   ! Initial specific humidity
                                !  (Kg Kg**-1)
     &  ,rp_theta_levels(row_length, rows,nlevs)                        &
                                                  ! Reciprocal pressure
                                !  ((HPa or mb)**-1)
     &  ,tbar0(nlevs)                                                   &
                                ! dummy var fo I/Os
     &  ,tbar1(row_length, rows,nlevs)                                  &
                                ! Mean temp. for Jan. and July
     &  ,tbar2(row_length, rows,nlevs)                                  &
                                          !  (K)
     &  ,tgrad0(nlevs)                                                  &
                                ! Dummy var for I/Os
     &  ,tgrad1(row_length, rows,nlevs)                                 &
                                        ! Gradient temp.
     &  ,tgrad2(row_length, rows,nlevs)                                 &
                                        !  for Jan. and July (K km**-1)
     &  ,ti(row_length, rows,nlevs)                                     &
                                        ! Initial temps. (K)
     &  ,tsd0(nlevs)                                                    &
                                ! dummy var fo I/Os
     &  ,tsd1(row_length, rows,nlevs)                                   &
                                       ! SD of temp. for Jan. and July
     &  ,tsd2(row_length, rows,nlevs)                                   &
                                          !  (K)
     &  ,tstar0                                                         &
                                ! Dummy for I/Os
     &  ,tstar1(row_length, rows)                                       &
                                          ! Surface temperature for
     &  ,tstar2(row_length, rows)                                       &
                                          !  Jan and July (K)
     &  ,vnbar0(nlevs)                                                  &
                                ! Dummy var for I/Os
     &  ,vnbar1(row_length, rows,nlevs)                                 &
                                        ! Mean horizontal velocity VN
     &  ,vnbar2(row_length, rows,nlevs)                                 &
                                        !  for Jan. and July (m s**-1)
     &  ,vnsd0(nlevs)                                                   &
                                ! Dummy for I/Os
     &  ,vnsd1(row_length, rows,nlevs)                                  &
                                        ! SD horizontal velocity VN
     &  ,vnsd2(row_length, rows,nlevs)                                  &
                                        !  for Jan. and July (m s**-1)
     &  ,vpbar0(nlevs)                                                  &
                                ! Dummy var for I/O.
     &  ,vpbar1(row_length, rows,nlevs)                                 &
                                        ! Mean horizontal velocity VP
     &  ,vpbar2(row_length, rows,nlevs)                                 &
                                        !  for Jan. and July (m s**-1)
     &  ,wbar0(ntrop)                                                   &
                                ! Dummy for I/Os
     &  ,wbar1(row_length, rows,ntrop)                                  &
                                        ! Mean vertical velocity
     &  ,wbar2(row_length, rows,ntrop)                                  &
                                        !  for Jan and July (mb s**-1)
     &  ,wsd0(ntrop)                                                    &
                                ! Dummy for I/Os
     &  ,wsd1(row_length, rows,ntrop)                                   &
                                        ! SD vertical velocity
     &  ,wsd2(row_length, rows,ntrop)                                   &
                                        !  for Jan and July (mb s**-1)
     &  ,xt                     ! Argument of SIN distribution
                                ! eqn. 2.33

      Character (Len=10) :: c1fmt, c2fmt, c3fmt, c4fmt, c5fmt
      Character (Len=200) :: c6fmt, c7fmt, c8fmt, c9fmt, c10fmt
      Character (Len=120) :: c11fmt, c12fmt

      character (len=25) :: txt ! A temporary variable used to test
                                ! the format of the file on unit 25
!  Set up Format statements

      ! The numbers in the files on units 25 and 26 may have been
      ! written with or without decimal points in them. We need to
      ! determine which so we can get the format statements right.
      ! Read the first line of unit 25 into a character variable
      Read (25,'(A)') txt
      ! Does it contain any decimal points?
      if (index(txt,'.') == 0) then
         ! No - use the following format statements
         c1fmt = '(6e8.2)'
         c2fmt = '(10f8.2)'
         c3fmt = '(5e8.2)'
         c4fmt = '(4f7.2)'
         c5fmt = '(f7.6)'
      else
         ! Yes - assume the numbers are one character longer
         c1fmt = '(6e9.2)'
         c2fmt = '(10f9.2)'
         c3fmt = '(5e9.2)'
         c4fmt = '(4f8.2)'
         c5fmt = '(f8.6)'
      endif
      ! Rewind unit 25 in readiness for proper reading below
      rewind(25)

      c6fmt = '("Climate forcing data for july"'                        &
     &  //',/," ________________________________"'                      &
     &  //',/," lat=",f7.2,"  long=",f7.2)'
      c7fmt = '("  tstar (K)  ",f7.2/," tuning factor " ,f7.2,/,'       &
     &  //'" dayno. relative to winter solstice ",f7.2)'
      c8fmt = '("0          level",i2,9(4x,"level",i2))'
      c9fmt = '('                                                       &
     & //'" tmn K    ",10(1e10.3,1x)/," tsd K    ",'                    &
     & //'10(1e10.3,1X)/," tgrd K/km",10(1e10.3,1x)/,'                  &
     & //'" p mb     ",10(1e10.3,1X)/,'                                 &
     & //'" vnmn m/s ",10(1e10.3,1x)/," vpmn m/s ",'                    &
     & //'10(1e10.3,1x)/," vnsd m/s ",10(1e10.3,1x))'
      c10fmt = '('                                                      &
     & //'" dmn K    ",10(1e10.3,1x)/," dgrd K/km",10(1e10.3,1x)/)'
      c11fmt = '('                                                      &
     & //'" wbar mb/s",10(1e10.3,1x)/," wsd mb/s ",10(1e10.3,1x)/)'
      c12fmt = '("1Climate forcing data for january"'                   &
     & //',/," ________________________________"'                       &
     & //',/," lat=",f7.2,"  long=",f7.2)'
!
!     Read climate stats for January and July
!
!     Each column is read in a set of dummy variables,
!     then, the values are copied accross to the real
!     arrays. This is to ensure conistency between
!     one column and multicolumns runs, even though it creates
!     redundancy
      Do m = 1, rows
        Do l = 1, row_length
          Read (25,c4fmt) lat0, long0, tstar0, alfad0
          Read (25,c2fmt) tbar0
          Read (25,c2fmt) tsd0
          Read (25,c2fmt) p0
          Read (25,c2fmt) dbar0
          Read (25,c3fmt) tgrad0
          Read (25,c3fmt) dgrad0
          Read (25,c2fmt) vnbar0
          Read (25,c2fmt) vpbar0
          Read (25,c2fmt) vnsd0
          Read (25,c1fmt) wbar0
          Read (25,c1fmt) wsd0
          Read (25,c5fmt) daysol1
          lat(l,m) = lat0
          long(l,m) = long0
          tstar1(l,m) = tstar0
          alfad1(l,m) = alfad0
        Do k = 1, nlevs
            tbar1(l,m,k) = tbar0(k)
            tsd1(l,m,k) = tsd0(k)
            p1(l,m,k) = p0(k)
            tgrad1(l,m,k) = tgrad0(k)
            vnbar1(l,m,k) = vnbar0(k)
            vpbar1(l,m,k) = vpbar0(k)
            vnsd1(l,m,k) = vnsd0(k)
        enddo
          p1(l,m,nlevs+1) = p0(nlevs+1)
        Do k = 1, nwet
            dbar1(l,m,k) = dbar0(k)
            dgrad1(l,m,k) = dgrad0(k)
        enddo
        Do k = 1, ntrop
            wbar1(l,m,k) = wbar0(k)
            wsd1(l,m,k) = wsd0(k)
        enddo

          Read (26,c4fmt) lat0, long0, tstar0, alfad0
          Read (26,c2fmt) tbar0
          Read (26,c2fmt) tsd0
          Read (26,c2fmt) p0
          Read (26,c2fmt) dbar0
          Read (26,c3fmt) tgrad0
          Read (26,c3fmt) dgrad0
          Read (26,c2fmt) vnbar0
          Read (26,c2fmt) vpbar0
          Read (26,c2fmt) vnsd0
          Read (26,c1fmt) wbar0
          Read (26,c1fmt) wsd0
          Read (26,c5fmt) daysol2
          lat(l,m) = lat0
          long(l,m) = long0
          tstar2(l,m) = tstar0
          alfad2(l,m) = alfad0
        Do k = 1, nlevs
            tbar2(l,m,k) = tbar0(k)
            tsd2(l,m,k) = tsd0(k)
            p2(l,m,k) = p0(k)
            tgrad2(l,m,k) = tgrad0(k)
            vnbar2(l,m,k) = vnbar0(k)
            vpbar2(l,m,k) = vpbar0(k)
            vnsd2(l,m,k) = vnsd0(k)
        enddo
          p2(l,m,nlevs+1) = p0(nlevs+1)
        Do k = 1, nwet
            dbar2(l,m,k) = dbar0(k)
            dgrad2(l,m,k) = dgrad0(k)
        enddo
        Do k = 1, ntrop
            wbar2(l,m,k) = wbar0(k)
            wsd2(l,m,k) = wsd0(k)
        enddo

      enddo                     ! l
      End Do                      ! m
!
!     Calculate amplitude and mean of annual sinusoidal distribution
!     Eqs 10 and 11
!
! DEPENDS ON: abnew
      Call ABNEW( tstar1, tstar2, tstara, tstarb, row_length,           &
     &            rows, 1)
! DEPENDS ON: abnew
      Call ABNEW( p1, p2, pa, pb, row_length,                           &
     &            rows, nlevs+1)
! DEPENDS ON: abnew
      Call ABNEW( alfad1, alfad2, alfada, alfadb, row_length,           &
     &            rows, 1)
! DEPENDS ON: abnew
      Call ABNEW(  tbar1,  tbar2,  tbara,  tbarb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW(   tsd1,   tsd2,   tsda,   tsdb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW( tgrad1, tgrad2, tgrada, tgradb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW(  dbar1,  dbar2,  dbara,  dbarb, row_length,           &
     &            rows, nwet)
! DEPENDS ON: abnew
      Call ABNEW( dgrad1, dgrad2, dgrada, dgradb, row_length,           &
     &            rows, nwet)
! DEPENDS ON: abnew
      Call ABNEW( vnbar1, vnbar2, vnbara, vnbarb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW(  vnsd1,  vnsd2,  vnsda,  vnsdb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW( vpbar1, vpbar2, vpbara, vpbarb, row_length,           &
     &            rows, nlevs)
! DEPENDS ON: abnew
      Call ABNEW(  wbar1,  wbar2,  wbara,  wbarb, row_length,           &
     &            rows, ntrop)
! DEPENDS ON: abnew
      Call ABNEW(   wsd1,   wsd2,   wsda,   wsdb, row_length,           &
     &            rows, ntrop)
!
!     Calculate constants for annual cycle used in eqn. 12
!
      atime = 2. * pi / andayy
      btime = pi * (.5-2.*daysol1)
!
!     Calculate argument of sinusoidal distribution (eqn. 12)
!
! DEPENDS ON: daynew
      xt = DAYNEW (atime, btime, dayno)
!
!     Calculate sinusoidal distribution (eqn. 12)
!
! DEPENDS ON: xnew
      Call XNEW(p_in, pa, pb, row_length, rows,                         &
     &                                          nlevs+1, xt)
! DEPENDS ON: xnew
      Call XNEW(ti, tbara, tbarb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call XNEW(q, dbara, dbarb, row_length, rows, nwet, xt)
!     Calculate default initial profile for Q
!
      Do k = 1, nwet
        Do j = 1, rows
          Do i = 1, row_length
            dewpt(i,j,k,1) = ti(i,j,k) - q(i,j,k)
        enddo
      enddo
      enddo
! DEPENDS ON: qsat
      Call QSAT(qi, dewpt(1,1,1,1), p_theta_levels(1,1,1),              &
     &                 (row_length*rows*nwet))
      Do k = 1, nlevs
        Do j = 1, rows
          Do i = 1, row_length
            t(i,j,k) = ti(i,j,k)
          enddo
        enddo
      enddo
      Do k = 1, nwet
        do j = 1, rows
          Do i = 1, row_length
            q(i,j,k) = qi(i,j,k)
          enddo
        enddo
      enddo                     ! i
!
!*********************************************************************
!     Print out climate datasets for January and July as read in
!     This section of code is very long but is necessary for
!     flexiblity ie to cope with any number of levels
!*********************************************************************
!
      daysol1 = daysol1 * andayy
      daysol2 = daysol2 * andayy

      Do l = 1 , row_length
        Do m = 1, rows
!       Transfer the arrays back to their 1D versions:
          lat0 = lat(l,m)
          long0 = long(l,m)
          tstar0 = tstar1(l,m)
          alfad0 = alfad1(l,m)
        Do k = 1, nlevs
            tbar0(k) = tbar1(l,m,k)
            tsd0(k) = tsd1(l,m,k)
            p0(k) = p1(l,m,k)
            tgrad0(k) = tgrad1(l,m,k)
            vnbar0(k) = vnbar1(l,m,k)
            vpbar0(k) = vpbar1(l,m,k)
            vnsd0(k) = vnsd1(l,m,k)
        enddo
        Do k = 1, nwet
            dbar0(k)= dbar1(l,m,k)
            dgrad0(k)= dgrad1(l,m,k)
        enddo
        Do k = 1, ntrop
            wbar0(k)= wbar1(l,m,k)
            wsd0(k)= wsd1(l,m,k)
        enddo

          if (row_length*rows  >   1) Write (11,*) 'Column no ', l,m
          Write (11,c6fmt) lat0, long0
          Write (11,c7fmt) tstar0, alfad0, daysol1
!
!       Set format statements
!
        cfmt = '(''          '',  (1e10.3,1x))'
        ctfmt = '(''0        '',  (3x,''level'',i2,1x))'
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with nlevs
!
        If ( mod(nlevs,10)  ==  0) then
          nlevsrows = int(nlevs/10)
          lastrow = 10
        else
          nlevsrows = int(nlevs/10) + 1
          lastrow = mod(nlevs,10)
        endif
        Do nlevscount = 1, nlevsrows
          element = 10 * (nlevscount-1)
          If (nlevscount  <   nlevsrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,c8fmt) (element+i, i = 1, 10)
            Write (11,c9fmt) (tbar0(element+i), i = 1, 10),             &
     &        (tsd0(element+i), i = 1, 10),                             &
     &        (tgrad0(element+i), i = 1, 10),                           &
     &          (p0(element+i), i = 1, 10),                             &
     &        (vnbar0(element+i), i = 1, 10),                           &
     &        (vpbar0(element+i), i = 1, 10),                           &
     &        (vnsd0(element+i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)')lastrow
            Write(11,ctfmt)(element+i,i= 1, lastrow)
            Write(cfmt(15:16),'(i2)')lastrow
            Write(cfmt(4:12),'(''tmn k    '')')
            Write(11,cfmt)(tbar0(i+element),i= 1, lastrow)
            Write(cfmt(4:12),'(''tsd k    '')')
            Write(11,cfmt)(tsd0(i+element),i= 1, lastrow)
            Write(cfmt(4:12),'(''tgrd k/km'')')
            Write(11,cfmt)(tgrad0(i+element),i= 1, lastrow)
            Write(cfmt(4:12),'(''P mb '')')
              Write(11,cfmt)(p0(i+element),i= 1, lastrow)
              Write(cfmt(4:12),'(''vnmn m/s '')')
            Write(11,cfmt)(vnbar0(i+element),i= 1, lastrow)
            Write(cfmt(4:12),'(''vpmn m/s '')')
            Write(11,cfmt)(vpbar0(i+element),i= 1, lastrow)
            Write(cfmt(4:12),'(''vnsd m/s '')')
            Write(11,cfmt)(vnsd0(i+element),i= 1, lastrow)
          endif
        enddo
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with NWET
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
            Write (11,c8fmt) (element+i,i = 1, 10)
            Write (11,c10fmt) (dbar0(element+i),i = 1, 10),             &
     &        (dgrad0(element+i),i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!         format statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)') lastrow
            Write (11,ctfmt) (element+i, i = 1, lastrow)
            Write (cfmt(15:16),'(i2)') lastrow
            Write (cfmt(4:12),'(''dmn k    '')')
            Write (11,cfmt) (dbar0(i+element), i = 1, lastrow)
            Write (cfmt(4:12),'(''dgrd k/km'')')
            Write (11,cfmt) (dgrad0(i+element), i = 1, lastrow)
          endif
        enddo
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with NTROP
!
        If ( mod(ntrop,10)  ==  0) then
          ntroprows = int(ntrop/10)
          lastrow = 10
        else
          ntroprows = int(ntrop/10) + 1
          lastrow = mod(ntrop,10)
        endif
        do ntropcount = 1,ntroprows
          element = 10 * (ntropcount-1)
          if ( ntropcount  <   ntroprows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,c8fmt) (element+i,i=1,10)
            Write (11,c11fmt) (wbar0(element+i), i = 1, 10),            &
     &        (wsd0(element+i),i=1,10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)') lastrow
            Write (11,ctfmt) (element+i, i = 1, lastrow)
            Write (cfmt(15:16),'(i2)') lastrow
            Write (cfmt(4:12),'(''wmn mb/s '')')
            Write (11,cfmt) (wbar0(i+element), i = 1, lastrow)
            Write (cfmt(4:12),'(''wsd mb/s '')')
            Write (11,cfmt) (wsd0(i+element), i = 1, lastrow)
          endif
        enddo

!       Transfer the arrays back to their 1D versions:
          lat0 = lat(l,m)
          long0 = long(l,m)
          tstar0 = tstar2(l,m)
          alfad0 = alfad2(l,m)
        Do k = 1, nlevs
            tbar0(k) = tbar2(l,m,k)
            tsd0(k) = tsd2(l,m,k)
            tgrad0(k) = tgrad2(l,m,k)
            p0(k) = p2(l,m,k)
            vnbar0(k) = vnbar2(l,m,k)
            vpbar0(k) = vpbar2(l,m,k)
            vnsd0(k) = vnsd2(l,m,k)
        enddo
        Do k = 1, nwet
            dbar0(k)= dbar2(l,m,k)
            dgrad0(k)= dgrad2(l,m,k)
        enddo
        Do k = 1, ntrop
            wbar0(k)= wbar2(l,m,k)
            wsd0(k)= wsd2(l,m,k)
        enddo


          Write (11,c12fmt) lat0, long0
          Write (11,c7fmt) tstar0, alfad0, daysol2
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with nlevs
!
        If ( mod(nlevs,10) == 0) then
          nlevsrows = int(nlevs/10)
          lastrow = 10
        else
          nlevsrows = int(nlevs/10) + 1
          lastrow = mod(nlevs,10)
        endif
        Do nlevscount = 1, nlevsrows
          element = 10*(nlevscount-1)
          if (nlevscount  <   nlevsrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,c8fmt) (element+i, i = 1, 10)
            Write (11,c9fmt) (tbar0(element+i),i = 1, 10),              &
     &        (tsd0(element+i),i = 1, 10),                              &
     &        (tgrad0(element+i),i = 1, 10),                            &
     &          (p0(element+i),i = 1, 10),                              &
     &        (vnbar0(element+i),i = 1, 10),                            &
     &        (vpbar0(element+i),i = 1, 10),                            &
     &        (vnsd0(element+i),i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           a character string. This will enable a variable format
!           to be created eg NF10.6 where N is the no. of elements in
!           the last row which can be written into the format
!           statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)') lastrow
            Write (11,ctfmt) (element+i,i = 1, lastrow)
            Write (cfmt(15:16),'(i2)') lastrow
            Write (cfmt(4:12),'(''tmn k    '')')
            Write (11,cfmt) (tbar0(i+element), i = 1, lastrow)
            Write (cfmt(4:12),'(''tsd k    '')')
            Write (11,cfmt) (tsd0(i+element),i = 1, lastrow)
            Write (cfmt(4:12),'(''tgrd k/km'')')
            Write (11,cfmt) (tgrad0(i+element),i = 1, lastrow)
            Write (cfmt(4:12),'(''P mb '')')
              Write (11,cfmt) (p0(i+element),i = 1, lastrow)
              Write (cfmt(4:12),'(''vnmn m/s '')')
            Write (11,cfmt) (vnbar0(i+element),i = 1, lastrow)
            Write (cfmt(4:12),'(''vpmn m/s '')')
            Write (11,cfmt) (vpbar0(i+element),i = 1, lastrow)
            Write (cfmt(4:12),'(''vnsd m/s '')')
            Write (11,cfmt) (vnsd0(i+element),i = 1, lastrow)
          endif
        enddo
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with nwet
!
        If ( mod(nwet,10)  ==  0) then
          nwetrows = int(nwet/10)
          lastrow = 10
        else
          nwetrows = int(nwet/10) + 1
          lastrow = mod(nwet,10)
        endif
        Do nwetcount = 1, nwetrows
          element = 10 * (nwetcount-1)
          If (nwetcount  <   nwetrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,c8fmt) (element+i,i = 1, 10)
            Write (11,c10fmt) (dbar0(element+i),i = 1, 10),             &
     &        (dgrad0(element+i),i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)') lastrow
            Write (11,ctfmt) (element+i,i = 1, lastrow)
            Write (cfmt(15:16),'(i2)') lastrow
            Write (cfmt(4:12),'(''dmn k/km '')')
            Write (11,cfmt) (dbar0(i+element),i = 1, lastrow)
            Write (cfmt(4:12),'(''dgrd k/km'')')
            Write (11,cfmt) (dgrad0(i+element),i = 1, lastrow)
          endif
        enddo
!
!       Calculate no. of rows and no. of elements in last row for
!       variables with ntrop
!
        If ( mod(ntrop,10) == 0) then
          ntroprows = int(ntrop/10)
          lastrow = 10
        else
          ntroprows = int(ntrop/10) + 1
          lastrow = mod(ntrop,10)
        endif
        Do ntropcount = 1, ntroprows
          element = 10*(ntropcount-1)
          If (ntropcount  <   ntroprows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (11,c8fmt) (element+i, i = 1, 10)
            Write (11,c11fmt) (wbar0(element+i), i = 1, 10),            &
     &        (wsd0(element+i),i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           formatC to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           format statement via an internal write statement.
!
            Write (ctfmt(14:15),'(i2)') lastrow
            Write (11,ctfmt)(element+i, i = 1, lastrow)
            Write (cfmt(15:16),'(i2)') lastrow
            Write (cfmt(4:12),'(''wmn mb/s '')')
            Write (11,cfmt)(wbar0(i+element), i = 1, lastrow)
            Write (cfmt(4:12),'(''wsd mb/s '')')
            Write (11,cfmt)(wsd0(i+element), i = 1, lastrow)
          endif
        enddo
        enddo                   ! m
      enddo                     ! l
!
!
      Return
      END SUBROUTINE INITSTAT
!
#endif
