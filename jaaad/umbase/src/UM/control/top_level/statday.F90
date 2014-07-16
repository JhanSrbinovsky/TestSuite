#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine statday
! Purpose:-           To calculate statistical forcing required each
!                     CLIM_CHANGE days
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!  5.3     16/05/01   Updated to 5.3          Z. Gardner
!  5.4     05/09/01   Checked for 5.4         M. Hughes
!=====================================================================
!
      Subroutine STATDAY(                                               &
!     ! IN leading dimensions of arrays
     &  row_length, rows, nlevs, nwet, ntrop,                           &
!     !
     &  atime, btime, dayno, deltan, daycount,                          &
     &  tbara, tbarb, tsda, tsdb, dbara, dbarb, vnbara, vnbarb,         &
     &  vnsda, vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb,         &
     &  alfada, alfadb, pa, pb, p, tgrada, tgradb,                      &
     &  dgrada, dgradb, cort, cord, corvn, corw, tdash, ddash,          &
     &  ctbar, ctsd, at, cdbar, cdsd, ad,                               &
     &  cvnbar, cvnsd, avn, cwbar, cwsd, aw, tbar, tsd, dbar, dsd,      &
     &  vnbar, vnsd, vpbar, wbar, wsd)

      Implicit none
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  row_length, rows                                                &
                                ! IN dimensions of arrays
     &  ,nlevs                                                          &
                                ! IN Number of levels of the model.
     &  ,nwet                                                           &
                                ! IN Number of humidity model levels
     &  ,ntrop                                                          &
                                ! IN Max number of levels in the
                                !   troposphere
     &  ,dayno                                                          &
                                ! IN Dayno. relative to winter
                                !    solstice
     &  ,daycount               ! IN Daynumber (ie 1 is 1st january)
      Real                                                              &
     &  alfada(row_length, rows)                                        &
                                ! IN Amplitude and mean of seasonal
     &  ,alfadb(row_length, rows)                                       &
                                 !    variation of tuning factor
     &  ,at(row_length, rows,nlevs-1)                                   &
                                     ! OUT Variable a in eqn 2.22
     &  ,ad(row_length, rows,nwet-1)                                    &
                                     !    used to calculate mean of
                                !    random variable for
                                !    temp. and dew point depression
     &  ,atime,btime                                                    &
                                ! IN Constants for calculating
                                !    annual cycle in eqn. 2.33
     &  ,avn(row_length, rows,nlevs-1)                                  &
                                      ! OUT Variable a in eqn 2.22
     &  ,aw(row_length, rows,ntrop-1)                                   &
                                      !    used to calculate mean of
                                !     random variable for
                                !     horiz. and vert velocity
     &  ,cord(row_length, rows)                                         &
                                ! IN vertical correlation
                                !    coeff. for dew pt. depress(0.9)
     &  ,cort(row_length, rows)                                         &
                                ! IN vertical correlation
                                !    coeff. for temp. (0.9)
     &  ,corvn(row_length, rows)                                        &
                                ! IN vertical correlation
                                !    coeff. for velocity VN (0.5)
     &  ,corw(row_length, rows)                                         &
                                ! IN vertical correlation
                                !    coeff. for velocity W (0.5)
     &  ,cdbar(row_length, rows,nwet)                                   &
                                ! OUT Mean and SD of random variable
     &  ,cdsd(row_length, rows,nwet)                                    &
                                !    for dew pt. depression (eqns
                                !    2.22 and 2.23)
     &  ,ctbar(row_length, rows,nlevs)                                  &
                                       ! OUT Mean and SD of random
     &  ,ctsd(row_length, rows,nlevs)                                   &
                                !    variable for temp. (eqns 2.22
                                !    and 2.23)
     &  ,cvnbar(row_length, rows,nlevs)                                 &
                                ! OUT Mean and SD of random variable
     &  ,cvnsd(row_length, rows,nlevs)                                  &
                                      !    for velocity VN (eqns 2.22
                                !    and 2.23)
     &  ,cwbar(row_length, rows,ntrop)                                  &
                                ! OUT Mean and SD of random variable
     &  ,cwsd(row_length, rows,ntrop)                                   &
                                !    for vertical velocity (eqns 2.22
                                !    and 2.23)
     &  ,dbar(row_length, rows,nwet)                                    &
                                     ! OUT Mean and SD dewpoint
     &  ,dsd(row_length, rows,nwet)                                     &
                                !    depression at daynumber relative
                                !    to winter solstice (K)
     &  ,dbara(row_length, rows,nwet)                                   &
                                ! IN Amplitude and mean of seasonal
     &  ,dbarb(row_length, rows,nwet)                                   &
                                      !    variation of mean dew pt.
                                !    depression (K)
     &  ,ddash(row_length, rows,nwet)                                   &
                                      ! OUT Dew pt. corrections (K)
     &  ,deltan(row_length, rows)                                       &
                                  ! IN Radius of area (m)
     &  ,dgrada(row_length, rows,nwet)                                  &
                                ! IN Amplitude and mean of seasonal
     &  ,dgradb(row_length, rows,nwet)                                  &
                                !    variation of dew pt. depression
                                ! gradient (K km^-1)
     &  ,p(row_length, rows,nlevs+1)                                    &
                                ! OUT Pressure on rho levels (Pa)
     &  ,pa(row_length, rows, nlevs+1)                                  &
                                ! OUT Amplitude and mean of seasonal
     &  ,pb(row_length, rows, nlevs+1)                                  &
                                !    variation of pressure
     &  ,tdash(row_length, rows,nlevs)                                  &
                                      ! Temp. correction (K)
     &  ,tbar(row_length, rows,nlevs)                                   &
                                      ! OUT Mean and SD temperature at
     &  ,tsd(row_length, rows,nlevs)                                    &
                                      !    daycount days from winter
                                !    solstice (K)
     &  ,tbara(row_length, rows,nlevs)                                  &
                                      ! IN Amplitude and mean of
     &  ,tbarb(row_length, rows,nlevs)                                  &
                                !    seasonal variation of temp. (K)
     &  ,tgrada(row_length, rows,nlevs)                                 &
                                ! IN Amplitude and mean of seasonal
     &  ,tgradb(row_length, rows,nlevs)                                 &
                                !    variation of temp. gradient
                                !    (K km^-1)
     &  ,tsda(row_length, rows,nlevs)                                   &
                                ! IN Amplitude and mean of seasonal
     &  ,tsdb(row_length, rows,nlevs)                                   &
                                !    variation of SD of temp. (K)
     &  ,vnbar(row_length, rows,nlevs)                                  &
                                ! OUT Mean and SD velocity VN at
     &  ,vnsd(row_length, rows,nlevs)                                   &
                                      !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vpbar(row_length, rows,nlevs)                                  &
                                       ! OUT Mean  velocity VP at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vnbara(row_length, rows,nlevs)                                 &
                                ! IN Amplitude and mean of seasonal
     &  ,vnbarb(row_length, rows,nlevs)                                 &
                                       !    variation of velocity VN
                                !    (m s^-1)
     &  ,vnsda(row_length, rows,nlevs)                                  &
                                ! IN Amplitude and mean of seasonal
     &  ,vnsdb(row_length, rows,nlevs)                                  &
                                !    variation of SD of velocity VN
                                !    (m s^-1)
     &  ,vpbara(row_length, rows,nlevs)                                 &
                                ! IN Amplitude and mean of seasonal
     &  ,vpbarb(row_length, rows,nlevs)                                 &
                                       !  variation of velocity VP
                                !    (m s^-1)
     &  ,wbar(row_length, rows,ntrop)                                   &
                                      ! OUT Mean and SD vertical
     &  ,wsd(row_length, rows,ntrop)                                    &
                                     !    velocity at daycount days
                                !    from winter solstice (mb s^-1)
     &  ,wbara(row_length, rows,ntrop)                                  &
                                ! IN Amplitude and mean of seasonal
     &  ,wbarb(row_length, rows,ntrop)                                  &
                                !    variation of SD of vert. vel.
                                !    (mb s^-1)
     &  ,wsda(row_length, rows,ntrop)                                   &
                                ! IN Amplitude and mean of seasonal
     &  ,wsdb(row_length, rows,ntrop)
                                !    variation of SD of vert. vel.
                                !    (mb s^-1)
!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
      Real                                                              &
     &  alfad(row_length, rows)                                         &
                                ! Tuning factor at daycount days
                                !  from winter solstice
     &  ,daynew                                                         &
                                ! function to calculate SIN
                                !  of argument (eqn 2.33)
     &  ,dgrad                                                          &
                                ! dew pt. depression gradient
     &  ,tgrad                                                          &
                                ! Temp. gradient
     &  ,xt                     ! Argument of SIN distribution
                                !  (eqn. 2.33)
      Integer                                                           &
     &  i, j, k                 ! Loop counter
!
!     Calculate argument of SIN (in eqn. 12)
!
      If (daycount  ==  1) then
! DEPENDS ON: daynew
        xt = daynew(atime, btime, dayno)
      else
! DEPENDS ON: daynew
        xt = daynew(atime, btime, dayno+1)
      endif
!
!     Calculate sinusoidal distribution (eqn. 12)
!
! DEPENDS ON: xnew
      Call xnew(tbar, tbara, tbarb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call xnew(tsd, tsda, tsdb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call xnew(dbar, dbara, dbarb, row_length, rows, nwet, xt)
! DEPENDS ON: xnew
      Call xnew(vnbar, vnbara, vnbarb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call xnew(vnsd, vnsda, vnsdb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call xnew(vpbar, vpbara, vpbarb, row_length, rows, nlevs, xt)
! DEPENDS ON: xnew
      Call xnew(wbar, wbara, wbarb, row_length, rows, ntrop, xt)
! DEPENDS ON: xnew
      Call xnew(wsd, wsda, wsdb, row_length, rows, ntrop, xt)
! DEPENDS ON: xnew
      Call xnew(alfad, alfada, alfadb, row_length, rows, 1, xt)
! DEPENDS ON: xnew
      Call xnew(p, pa, pb, row_length, rows, 1, xt)

      Do i = 1, row_length
        Do j = 1, rows
        Do k = 1, nlevs
            tgrad = tgrada(i,j,k)*xt + tgradb(i,j,k)
            If (vnbar(i,j,k)  >   0.0) then
            tgrad = -tgrad
          endif
!
!         Calculate corrections for temp. and dew pt. depression
!
            tdash(i,j,k)=deltan(i,j)*tgrad*0.001
        enddo
        Do  k = 1, nwet
            dgrad = dgrada(i,j,k) * xt + dgradb(i,j,k)
            dsd(i,j,k) = alfad(i,j) * tsd(i,j,k)
            If (vnbar(i,j,k)  >   0.0) then
            dgrad = -dgrad
          endif
!
!         Calculate corrections for temp. and dew pt. depression
!
            ddash(i,j,k) = deltan(i,j) * dgrad*0.001
          enddo
        enddo
      enddo

!
!     Calculate mean and SD of random variable (eqns. 6 and 7)
!
! DEPENDS ON: acinit
      Call ACINIT(tbar, tsd, at, ctbar, ctsd, cort, nlevs,              &
     &                                        row_length, rows)
! DEPENDS ON: acinit
      Call ACINIT(dbar, dsd, ad, cdbar, cdsd, cord, nwet,               &
     &                                        row_length, rows)
! DEPENDS ON: acinit
      Call ACINIT(vnbar, vnsd, avn, cvnbar, cvnsd, corvn, nlevs,        &
     &  row_length, rows)
! DEPENDS ON: acinit
      Call ACINIT(wbar, wsd, aw, cwbar, cwsd, corw, ntrop,              &
     &                                        row_length, rows)
      Return
      END SUBROUTINE STATDAY
!
!
#endif
