#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine STATSTEP
! Purpose:-           To calculate statistical forcing required each
!                     timestep
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards + now includes interpolation
!                     of all forcing variables between subsequent days
!                     Last modified on 1/3/93 to include double-
!                     precision versions of NAG random number
!                     generator routines (ie G05) for compatability
!                     with workstation
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!  5.3     10/01      Update to 5.3          Z. Gardner
!  5.4  05/09/02      SCM Stats and Obs forcing changes - M. Hughes
!  5.5  05/09/02      Changed wind increments & initialise t_inc q_inc
!                                                          - M. Hughes

!=====================================================================
!
      Subroutine STATSTEP(                                              &
                                ! IN
     &  row_length, rows, nlevs, nwet, ntrop,                           &
                                !
     &  deltan, px, py, daysteps, stepcount, dayno,                     &
     &  tr, vnr, vpr, qr, wr, tbar, tsd, tdash, dbar, dsd, ddash,       &
     &  vnbar, vpbar, vnsd, wbar, wsd, ctbar, ctsd, at, cdbar,          &
     &  cdsd, ad, cvnbar, cvnsd, avn, cwbar, cwsd, aw,                  &
     &  press, rpress, u, v, w, t, q, prinstat, u_inc, v_inc,           &
     &  w_inc, t_inc, q_inc, daycount, timestep, iv, ntab, iy, idum)
      Implicit none
!


!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  row_length, rows                                                &
                                ! IN no  of model rows/columns
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q is
                                !    set
     &  ,ntrop                                                          &
                                ! IN Max number of levels in the
                                !    troposphere
     &  ,daycount                                                       &
                                ! IN Daynumber (1 represents
                                !    1st january)
     &  ,dayno                                                          &
                                ! IN Daynumber relative to winter
                                !    solstice
     &  ,daysteps                                                       &
                                ! IN No. of timesteps in 1 day
     &  ,stepcount                                                      &
                                ! IN Timestep counter
     &  ,ntab                                                           &
                                ! IN Dimension of array used in
                                !    random generator
     &  ,iv(ntab),iy,idum       ! IN state of number generator saved
                                !    on tape for continuation run

      Logical                                                           &
     &  prinstat                ! T if printout of increments
                                !  due to statistical forcing
                                !  required each timestep
      Real                                                              &
     &  ad(row_length, rows,nwet-1)                                     &
                                ! IN term a of equation 2.22 for dew
     &  ,at(row_length, rows,nlevs-1)                                   &
                                      !    pt depression and temp.
     &  ,avn(row_length, rows,nlevs-1)                                  &
                                ! IN term a of equation 2.22 for
     &  ,aw(row_length, rows,ntrop-1)                                   &
                                !    horiz. and vertical velocities

     &  ,cdbar(row_length, rows,nwet)                                   &
                                ! Mean and SD of random variable
     &  ,cdsd(row_length, rows,nwet)                                    &
                                !    for dew point depression
     &  ,ctbar(row_length, rows,nlevs)                                  &
                                ! IN Mean and SD of random variable
     &  ,ctsd(row_length, rows,nlevs)                                   &
                                      !    for temp.

     &  ,cvnbar(row_length, rows,nlevs)                                 &
                                ! IN Mean and SD of random variable
     &  ,cvnsd(row_length, rows,nlevs)                                  &
                                       !    for velocity VN
     &  ,cwbar(row_length, rows,ntrop)                                  &
     &  ,cwsd(row_length, rows,ntrop)                                   &
                                ! IN Mean and SD of random variable
                                !    for vertical velocity
     &  ,dbar(row_length, rows,nwet)                                    &
                                     ! IN Mean and SD dewpoint
     &  ,dsd(row_length, rows,nwet)                                     &
                                    !    depression at daycount days
                                !    from winter solstice (K)
     &  ,ddash(row_length, rows,nwet)                                   &
                                      ! IN Dew pt. corrections
     &  ,deltan(row_length, rows)                                       &
                                      ! IN Radius of area (m)
     &  ,dq1                                                            &
                                ! IN Spec. humidity differences
     &  ,dq2                                                            &
                                !    (Kg Kg^-1)
     &  ,dt1                                                            &
                                ! IN Temp. differences (K)
     &  ,dt2                                                            &
     &  ,press(row_length, rows,nlevs)                                  &
                                ! IN Pressure coordinates (Pa)
     &  ,px(row_length, rows,ntrop)                                     &
                                ! IN Reciprocal log functions for
     &  ,py(row_length, rows,ntrop-1)                                   &
                                !    calc. of vert. advection
                                !    used in eqns 2.12 and 2.13
     &  ,q(row_length, rows,nwet)                                       &
                                ! OUT Specific humidity (Kg Kg^-1)
     &  ,qr(row_length, rows,nwet,2)                                    &
                                ! INOUT Randomly sampled humidity
                                !    (Kg Kg^-1)
     &  ,rpress(row_length, rows,nlevs)                                 &
                                ! IN Reciprocal pressure for rho
                                !    levels (HPa or mb ^-1)
     &  ,t(row_length, rows,nlevs)                                      &
                                   ! OUT Temp (K)
     &  ,tbar(row_length, rows,nlevs)                                   &
                                ! IN Mean and SD temperature at
                                !    daycount days from
                                !    winter solstice (K)
     &  ,tdash(row_length, rows,nlevs)                                  &
                                      ! IN Temp. corrections (K)
     &  ,timestep                                                       &
                                ! IN model timestep (s)
     &  ,tsd(row_length, rows,nlevs)                                    &
                                ! IN SD of temp. at daycount days
                                !    from winter solstice (K)
     &  ,tr(row_length, rows,nlevs,2)                                   &
                                ! INOUT Randomly sampled temp. (K)
     &  ,u(row_length, rows,nlevs)                                      &
     &  ,v(row_length, rows,nlevs)                                      &
                                ! OUT Zonal and meridional winds
                                !    (m s^-1)
     &  ,vnbar(row_length, rows,nlevs)                                  &
                                ! IN Mean and SD velocity VN at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vnr(row_length, rows,nlevs,2)                                  &
                                ! INOUT Randomly sampled horizontal
                                !    velocity (m s^-1)
     &  ,vnsd(row_length, rows,nlevs)                                   &
                                ! IN Mean and SD velocity VN at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vpbar(row_length, rows,nlevs)                                  &
                                       ! IN Mean  velocity VP at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vpr(row_length, rows,nlevs,2)                                  &
                                ! INOUT Randomly sampled horizontal
                                !    velocity (m s^-1)
     &  ,w(row_length, rows,0:nlevs)                                    &
                                     ! Vertical velocity
     &  ,wbar(row_length, rows,ntrop)                                   &
     &  ,wsd(row_length, rows,ntrop)                                    &
                                      ! IN Mean and SD vertical
                                !    velocity at daycount days
                                !    from winter solstice (mb s^-1)
     &  ,wr(row_length, rows,ntrop,2)                                   &
                                ! INOUT Randomly sampled vertical
                                !    velocity (mb s^-1)
     &  ,u_inc(row_length, rows,nlevs)                                  &
                                       !OUT
                                ! Zonal and meridional wind
     &  ,v_inc(row_length, rows,nlevs)                                  &
                                      !OUT
                                !  increment due to large-scale
                                !  horizontal and vertical
     &  ,w_inc(row_length, rows,nlevs)                                  &
                                      !OUT
     &  ,t_inc(row_length, rows,nlevs)                                  &
                                ! OUT Temp increment due to
                                !    large-scale horizontal and
                                !    vertical advection
                                !    (K s^-1 day^-1)
     &  ,q_inc(row_length, rows,nwet)
                                ! Specific humidity increment
                                !  due to large-scale horizontal
                                !  and vertical advection
                                !  (kg kg^-1 s^-1 day^-1)
!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  i, i1, j,j1, k, l       ! Loop counters
      Real                                                              &
     &  cdd                                                             &
                                ! Randomly sampled variables for
     &  ,ct                                                             &
                                !  temp. and dew pt. depression
     &  ,cvn                                                            &
                                ! Randomly sampled variables for
     &  ,cw                                                             &
                                !  horizontal and vertical velocity
     &  ,d0                                                             &
                                ! Randomly sampled dew pt.
                                !  depression for 1st level
     &  ,dewpt(row_length, rows,nwet,2)                                 &
                                          ! Dew-point
     &  ,dr                                                             &
                                ! Randomly sampled dew pt. depression
                                !   (K)
     &  ,f1                                                             &
                                ! Used in calc of advection term
                                !  in equation 2.12
     &  ,n1, n2                                                         &
                                ! Constants for interpolation
     &  ,qk(row_length, rows,nwet)                                      &
                                   ! Factor to prevent Q becoming
                                !  negative
     &  ,qrint(row_length, rows,nwet)                                   &
                                ! Specific humidity (interpolated
                                !  values) (Kg Kg^-1)
     &  ,rdaysteps, rstepcount                                          &
                                ! Real values of timesteps in day
                                !  and timestep counter
     &  ,t0                                                             &
                                ! Randomly sampled temp. for 1st
                                !  level (K)
     &  ,trint(row_length, rows,nlevs)                                  &
                                ! Temp. (interpolated values) (K)
     &  ,vnrint(row_length, rows,nlevs)                                 &
                                ! Horizontal velocities (linearly
     &  ,vprint(row_length, rows,nlevs)                                 &
                                !  interpolated values) (m s^-1)
     &  ,wrint(row_length, rows,ntrop)                                  &
                                ! Vertical velocity (linearly
                                !  interpolated values)
                                !  (mb s^-1)
     &  ,G05DDE                 ! Function that samples randomly
                                ! from a Gaussian distribution
                                ! with a given mean and SD

!
!
      If (stepcount  ==  1) then
        If (daycount  ==  1) then
          i1 = 1
          j1 = 2
        else
          i1 = 2
          j1 = 2
!
!         Save state of Random Number Generator for continuation
!         STATS run done from tape to allow for the first day of a
!         STATS run, when G05DDE is used twice as many times (to set
!         up 2 profiles) and so the variables after forcing on a
!         continuation run would be different from an unbroken run.
!         - daycount ne 1.
!
! DEPENDS ON: g05cfe
          Call G05CFE(idum,iv,iy)
          Do j = 1, row_length
            Do l = 1, rows
            Do k = 1, nlevs
                tr(j,l,k,j1-1) = tr(j,l,k,j1)
                vnr(j,l,k,j1-1) = vnr(j,l,k,j1)
                vpr(j,l,k,j1-1) = vpr(j,l,k,j1)
            enddo
            Do k = 1, nwet
                qr(j,l,k,j1-1) = qr(j,l,k,j1)
            enddo
            Do k = 1, ntrop
                wr(j,l,k,j1-1) = wr(j,l,k,j1)
            enddo
            End Do
          enddo
        endif
!
!---------------------------------------------------------------------
!       Create new profiles (2 for 1st day and 1 from then on)
!       G05DDE(A,B) returns a pseudo-random real number taken from
!       a normal (Gaussian) distribution with mean A and SD B.
!---------------------------------------------------------------------
!
        Do j = 1, row_length
          Do l = 1, rows
          Do k = i1, j1
! DEPENDS ON: g05dde
              t0 = G05DDE(tbar(j,l,1),tsd(j,l,1))
              tr(j,l,1,k) = t0+tdash(j,l,1)
! DEPENDS ON: g05dde
              d0 = G05DDE(dbar(j,l,1),dsd(j,l,1))
              dr = d0+ddash(j,l,1)
              dewpt(j,l,1,k) = tr(j,l,1,k)-dr
! DEPENDS ON: g05dde
              vnr(j,l,1,k) = G05DDE(vnbar(j,l,1),vnsd(j,l,1))
! DEPENDS ON: g05dde
              wr(j,l,1,k) = G05DDE(wbar(j,l,1),wsd(j,l,1))
            Do i = 1, nlevs-1
! DEPENDS ON: g05dde
                ct = G05DDE(ctbar(j,l,i),ctsd(j,l,i))
                t0 = at(j,l,i)*t0+ct
                tr(j,l,i+1,k) = t0+tdash(j,l,i+1)
! DEPENDS ON: g05dde
                cvn = G05DDE(cvnbar(j,l,i),cvnsd(j,l,i))
                vnr(j,l,i+1,k) = avn(j,l,i)*vnr(j,l,i,k)+cvn
            enddo
            Do i = 1, nwet-1
! DEPENDS ON: g05dde
                cdd = G05DDE(cdbar(j,l,i),cdsd(j,l,i))
                d0 = ad(j,l,i)*d0+cdd
                dr = d0+ddash(j,l,i+1)
                dewpt(j,l,i+1,k) = tr(j,l,i+1,k)-dr
            enddo
            Do i = 1, ntrop-1
! DEPENDS ON: g05dde
                cw = G05DDE(cwbar(j,l,i),cwsd(j,l,i))
                wr(j,l,i+1,k) = aw(j,l,i)*wr(j,l,i,k)+cw
            enddo
            Do i = 1, nlevs
! DEPENDS ON: g05dde
                vpr(j,l,i,k) = G05DDE(vpbar(j,l,i),vnsd(j,l,i))
            enddo
!
!           After the first profile is set up on the first day,
!           save state of Random Number Generator for continuation
!           STATS run done from tape to allow for the first day of a
!           STATS run, when G05DDF is used twice as many times (to set
!           up 2 profiles) and so the variables after forcing on a
!           continuation run would be different from an unbroken run.
!           - daycount EQ 1.
!
            If (k  ==  1 .and. l  ==  1)                                &
! DEPENDS ON: g05cfe
     &        Call G05CFE(idum,iv,iy)
          enddo
        enddo
        enddo
        Do k = i1, j1
! DEPENDS ON: qsat
          Call QSAT(qr(1,1,1,k), dewpt(1,1,1,k), press,                 &
     &            (row_length*rows*nwet))
        enddo
      endif                     !  stepcount=1
!
!     Interpolate between 2 values
!
      rdaysteps = real(daysteps)
      rstepcount = real(stepcount)
      n1 = (rdaysteps-rstepcount+1.) / rdaysteps
      n2 = (rstepcount-1.) / rdaysteps
      Do  k = 1, nlevs
        Do j = 1, rows
          Do l = 1, row_length
            trint(l,j,k) = n1 * tr(l,j,k,1) + n2 * tr(l,j,k,2)
            vnrint(l,j,k) = n1 * vnr(l,j,k,1) + n2 * vnr(l,j,k,2)
            vprint(l,j,k) = n1 * vpr(l,j,k,1) + n2 * vpr(l,j,k,2)
          enddo
        enddo
      enddo
      Do k = 1, nwet
        Do j = 1, rows
          Do l = 1, row_length
            qrint(l,j,k) = n1 * qr(l,j,k,1) + n2 * qr(l,j,k,2)
          enddo
        enddo
      enddo
      Do k = 1, ntrop
        Do j = 1, rows
          Do l = 1, row_length
            wrint(l,j,k) = n1 * wr(l,j,k,1) + n2 * wr(l,j,k,2)
          enddo
        enddo
      enddo


!
!     Set U and V increments
!
      Do k = 1, nlevs
        Do j = 1, rows
          Do l = 1, row_length
            u_inc(l,j,k) = - u(l,j,k) + vnrint(l,j,k)
            v_inc(l,j,k) = - v(l,j,k) + vprint(l,j,k)
          enddo
        enddo
      enddo
      Do k = 1, ntrop
        Do j = 1, rows
          Do l = 1, row_length
            w_inc(l,j,k) = - w(l,j,k) + wrint(l,j,k)

          enddo
        enddo
      enddo




!
!---------------------------------------------------------------------
!     Add vertical advection increments to T and Q (eqns. 2 and 3)
!---------------------------------------------------------------------
!
      Do k = 1, nlevs
        Do j = 1, rows
          Do l = 1, row_length
            t_inc(l,j,k) = 0.0
          enddo
        enddo
      enddo
      Do k = 1, nwet
        Do j = 1, rows
          Do l = 1, row_length
            q_inc(l,j,k) = 0.0
          enddo
        enddo
      enddo



      Do j = 1, rows
        Do l = 1, row_length
          dt1 = t(l,j,2) - t(l,j,1)
          dq1 = q(l,j,2) - q(l,j,1)
        Do i = 2, ntrop
            dt2 = t(l,j,i+1) - t(l,j,i)
            dq2 = q(l,j,i+1) - q(l,j,i)
            f1 = -wrint(l,j,i) * timestep * rpress(l,j,i)
            t_inc(l,j,i) = t_inc(l,j,i) + f1 * (dt2 * px(l,j,i) +       &
     &          dt1 * px(l,j,i-1)                                       &
     &        - (dt1 + dt2) * py(l,j,i-1) - .2856 * t(l,j,i))
            q_inc(l,j,i) = q_inc(l,j,i) + f1 * (dq2 * px(l,j,i) +       &
     &      dq1 * px(l,j,i-1)                                           &
     &        - (dq1+dq2) * py(l,j,i-1))
          dt1 = dt2
          dq1 = dq2
        enddo
      enddo
      enddo                     ! l
!
!     Printout increments due to vert. advection if required
!
! These are commented out for now so that the call can be added in at
! a later date if so required.
!      If (prinstat)
!     &  Call PRINTSUB(
!     ! IN
!     &  row_length, rows, nlevs, nwet,
!     !
!     &  ' Temperature/moisture profiles + incs due to vert advection '
!     &  ,stepcount , dayno, t, q, t_inc, q_inc)
!
!---------------------------------------------------------------------
!     Add horizontal increments to T and Q (eqns. 2 and 3)
!---------------------------------------------------------------------
!
      Do k = 1, nlevs
        Do j = 1, rows
          Do l = 1, row_length
            t_inc(l,j,k) = t_inc(l,j,k) + timestep * abs(vnrint(l,j,k)) &
     &        *    (trint(l,j,k)-t(l,j,k)) / deltan(l,j)
          enddo
        enddo
      enddo
!
!     This section prevents the increment from allowing
!     Q to become negative (lowest value Q can take is 1.E-6)
!
      Do k = 1, nwet
        Do j = 1, rows
          Do l = 1, row_length
            q_inc(l,j,k) = q_inc(l,j,k) + timestep * abs(vnrint(l,j,k)) &
     &        *          (qrint(l,j,k)-q(l,j,k))                        &
     &        /          deltan(l,j)
            qk(l,j,k) = -1 * q(l,j,k) + 1.e-6
            If (q_inc(l,j,k)  <   qk(l,j,k)) then
              q_inc(l,j,k) = qk(l,j,k)
          endif
          enddo
        enddo
      enddo
!
!     Printout increments due to horiz. advection if required
!
      If (prinstat)                                                     &
! DEPENDS ON: printsub
     &  Call PRINTSUB(                                                  &
!     ! IN
     &  row_length*rows, nlevs, nwet,                                   &
!     !
     &  ' Temperature/moisture profiles + incs due to advection'        &
     &  ,stepcount, dayno, t(1,1:1,1:nlevs), q(1,1:1,1:nlevs),          &
     &  t_inc(1,1:1,1:nlevs), q_inc(1,1:1,1:nlevs))
      Return
      END SUBROUTINE STATSTEP
#endif
