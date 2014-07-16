#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SUB_DATA
! Purpose:-           To be used in test runs when detailed
!                     sub-timestep diagnostics are required after
!                     calls to each subroutine or when budget calcs.
!                     are required and so called at start and
!                     end of meaning period
!
! Programmer:-        J. LEAN  15/6/91
!
! Modification History:
!  4.5     07/98    SCM integrated as a standard UM configuration
!                   JC Thil.
!  5.3     03/05/01 Changed to write out p, exner and rho instead of
!                   pstar.  Also extra horizontal dimension added.
!                   Z. Gardner
! 6.2      13/02/06  Removal of MOSES I code   Adrian Lock
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
      Subroutine SUB_DATA(                                              &
!     ! IN leading dimensions of arrays
     &  row_length, rows, model_levels, nwet                            &
     &  ,nfor, nbl_levs, nsoilt_levs, nsoilm_levs, ntrop                &
!     !
     &  ,title1, istep, ayear, aday, atime_string, rday                 &
     &  ,u, v, t, theta, q, qcl, qcf, lca, p,rho, exner, t_deep_soil    &
     &  ,smc,canopy, snodep, tstar, zh, z0msea, cca, iccb, icct, smcl)

      Implicit none
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  row_length                                                      &
                                ! IN x dimension.
     &  ,rows                                                           &
                                ! IN y dimension.
     &  ,model_levels                                                   &
                                       ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q is
                                !   set.
     &  ,nfor                                                           &
                                ! IN Number terms for observational
                                !   forcing
     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temperature
                                !  levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels
     &  ,ntrop                  ! IN Max number of levels in the
                                !  troposphere
      Character*35                                                      &
     &  title1                  ! Table heading
      Character*8                                                       &
     &  atime_string            ! Actual time at end of timestep
      Integer                                                           &
     &  iccb(row_length, rows)                                          &
                                          ! Convective cloud base,top
     &  ,icct(row_length, rows)                                         &
     &  ,aday                                                           &
                                ! Actual day number
     &  ,ayear                                                          &
                                ! Actual year number
     &  ,rday                                                           &
                                ! Runday number
     &  ,istep                  ! Timestep number
      Real                                                              &
     &  canopy(row_length, rows)                                        &
                                ! Surface/canopy water (kg m^-2)
     &  ,cca(row_length, rows)                                          &
                                ! Convective cloud amount
     &  ,lca(row_length, rows,nwet)                                     &
                                   ! Layer cloud amount (decimal
                                !  fraction)
     &  ,P(row_length, rows, model_levels)                              &
                                ! pressure on rho levels (Pa)
     &  ,exner(row_length, rows, model_levels)                          &
                                ! exner pressure on rho levels
     &  ,rho(row_length, rows, model_levels)                            &
                                ! density
     &  ,Q(row_length, rows,nwet)                                       &
                                 ! Specific humidity (kg kg^-1)
     &  ,Qcf(row_length, rows,nwet)                                     &
                                   ! Specific cloud ice (kg kg^-1)
     &  ,Qcl(row_length, rows,nwet)                                     &
                                   ! Specific cloud water (kg kg^-1)
     &  ,smc(row_length, rows)                                          &
                                   ! Soil moisture content (kg m^-2)
     &  ,smcl(row_length, rows,nsoilm_levs)                             &
                                   ! Soil moisture inlayers (kg m^-2)
     &  ,snodep(row_length, rows)                                       &
                                   ! Snow depth (kg m^-2)
     &  ,t_deep_soil(row_length, rows,nsoilt_levs)                      &
                                   ! Soil layer temps (K)
     &  ,t(row_length, rows,model_levels)                               &
                                   ! Temperature (K)
     &  ,theta(row_length, rows,model_levels)                           &
                                   ! Potential temperature (K)
     &  ,Tstar(row_length, rows)                                        &
                                   ! Surface temp.(K)
     &  ,u(row_length, rows,model_levels)                               &
     &  ,v(row_length, rows,model_levels)                               &
                                   ! Zonal,meridional wind (m s^-1)
     &  ,zh(row_length, rows)                                           &
                                   ! Boundary layer depth (m)
     &  ,Z0mSea(row_length, rows)  ! Sea surface roughness length(m)
!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  element                                                         &
                                ! Array element number
     &  ,lastrow                                                        &
                                ! Number of elements in the last row
     &  ,ndeeprows, ndeepcount                                          &
                                ! nbl_levs, nclds, ndeep, model_levels,
     &  ,nlevsrows, nlevscount                                          &
                                ! nwet, nsoil elements
     &  ,nwetrows, nwetcount                                            &
                                !
     &  ,ndeep                                                          &
                                ! No of deep soil layers ie
                                !  surface layer
     &  ,i,k, l     ! Loop counters
      Character*30  c1fmt
      Character*32  c2fmt
      Character (len=150) :: c3fmt, c4fmt, c5fmt, c6fmt, c7fmt, c8fmt
      Character (len=220) :: c9fmt, c10fmt, c11fmt, c12fmt, c13fmt
      Character (len=150) :: c14fmt, c15fmt, c16fmt, c17fmt, c18fmt
      Character (len=150) :: c19fmt
!---------------------------------------------------------------------
      Character*36  ctfmt
      c1fmt = '('' *           '',  (f8.3,3x))'
      c2fmt = '('' *           '',  (1pe10.3,1x))'

      c3fmt = '("Up to actual year ",i6,", day",i6,", time "'           &
     &      //', A8)'
      c4fmt = '(" Runday",i6," timestep",i5, "   primary variables "'   &
     &      //',a24)'
      c5fmt = '('                                                       &
     &//'"*********************************************************"'   &
     &//'"*********************************************************"'   &
     &//'"******")'
      c6fmt = '('                                                       &
     &//'"1********************************************************"'   &
     &//'"*********************************************************"'   &
     &//'"******")'
      c7fmt = '(" * Standard model variables",/,'                       &
     &    //'" + __________________________")'
      c8fmt = '("0* Variable   level",i2,9(4x,"level",i2))'
      c9fmt = '(" * U m s^-1  ",10(F8.3,3X), /,'                        &
     &    //'" * V m s^-1  ",10(F8.3,3X), /,'                           &
     &    //'" * T K       ",10(F8.3,3X), /,'                           &
     &    //'" * theta K   ",10(F8.3,3X), /,'                           &
     &    //'" * rho       ",10(F14.1,3X), /,'                          &
     &    //'" * p Pa      ",10(F8.2,3X), /,'                           &
     &    //'" * exner     ",10(F8.6,3X))'
      c10fmt = '(" * Rain and cloud variables",/,'                      &
     &    //'" + __________________________")'
      c11fmt = '(" * Q kg/kg   ",10(1pe10.3,1x), /,'                    &
     &    //'" * QCL kg/kg ",10(1pe10.3,1x), /,'                        &
     &    //'" * QCF kg/kg ",10(1pe10.3,1x), /,'                        &
     &    //'" * LCA       ",10(1pe10.3,1x))'
      c12fmt = '("0*   CCA",7x,"CCB",6x,"CCT")'
      c13fmt = '(" ",1pe10.3,i6,3x,i6)'
      c14fmt = '(" * Boundary layer and surface variables",/,'          &
     &    //'" + ______________________________________")'
      c15fmt = '(" *Tdeep K    ",10(F10.4))'
      c16fmt = '("0* Tstar",6X,"smc",6X,"canopy",4X,"snodep",'          &
     &    //'6X,"zh",6x,"Z0mSea")'
      c17fmt = '(" *",2X,"(K)",4X,"(kg m^-2)","(kg m^-2)",'             &
     &    //'"(kg m^-2)",3X,"(m)",6X,"(m)")'
      c18fmt = '(" ",f10.3,6f10.4)'
      c19fmt = '(" *SMCL       ",10(F10.4))'
      ctfmt = '(''0* Variable'',  (3x,''level'',i2,1x))'
!
!
!     Write out headings
!
      Write (22,c6fmt)
      Write (22,c4fmt) rday, istep, title1
      Write (22,c3fmt) ayear, aday, atime_string
      Write (22,c5fmt)
      Write (22,c7fmt)
      Write (22,c5fmt)
!
!     Loop on the sites

      Do l = 1, rows
        Do k = 1, row_length
          If (row_length*rows  >   1)                                   &
     &            Write (22,*) " Site No : ", k,l
!
!       Write out variables T theta U and V maximum of 10
!       variables per row
!
          If (mod(model_levels, 10)  ==  0) then
!
!         Calculate no. of rows and no. of elements in last row
!
            nlevsrows = int(model_levels/10)
          lastrow = 10
        else
            nlevsrows = int(model_levels/10) + 1
            lastrow = mod(model_levels,10)
        endif
        Do nlevscount = 1, nlevsrows
          element = 10 * (nlevscount-1)
          If (nlevscount  <   nlevsrows) then
!
!           Write out all complete rows ie of 10 variables per row
!
              Write (22,c8fmt) (element + i,i = 1, 10)
              Write (22,c9fmt) (u(k,l, element + i), i = 1, 10),        &
     &          (v(k,l, element + i), i = 1, 10),                       &
     &          (t(k,l, element + i), i = 1, 10),                       &
     &          (theta(k,l, element + i), i = 1, 10),                   &
     &          (rho(k,l,element + i), i = 1, 10),                      &
     &          (p(k,l,element + i), i = 1, 10),                        &
     &          (exner(k,l,element + i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           Format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           Format statement via an internal write statement.
!
            Write (ctfmt(16:17), '(i2)') lastrow
            Write (22,ctfmt) (element + i, i = 1, lastrow)
            Write (c1fmt(18:19), '(i2)')lastrow
            Write (c1fmt(6:15), '(''U m s^-1  '')')
              Write (22,c1fmt) (u(k,l, i + element), i = 1, lastrow)
            Write (c1fmt(6:15), '(''V m s^-1  '')')
              Write (22,c1fmt) (v(k,l, i + element), i = 1, lastrow)
            Write (c1fmt(6:15), '(''T K       '')')
              Write (22,c1fmt) (t(k,l, i + element), i = 1, lastrow)
            Write (c1fmt(6:15), '(''theta K   '')')
              Write (22,c1fmt) (theta(k,l, i + element), i = 1, lastrow)
              Write (c1fmt(6:15), '(''rho   '')')
              Write (22,c1fmt) (rho(k,l, i + element), i = 1, lastrow)
              Write (c1fmt(6:15), '(''p   '')')
              Write (22,c1fmt) (p(k,l, i + element), i = 1, lastrow)
              Write (c1fmt(6:15), '(''exner'')')
              Write (22,c1fmt) (exner(k,l, i + element), i = 1, lastrow)
          endif
        enddo
        Write (22,c5fmt)
        Write (22,c10fmt)
        Write (22,c5fmt)
!
!       Repeat above section of code for variables on NWET levels
!       Write out variables Q QCL QCF and LCA maximum of 10
!       variables per row
!
        If (mod(nwet,10)  ==  0) then
!
!         Calculate no. of rows and no. of elements in last row
!
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
            Write (22,c8fmt)(element + i,i = 1, 10)
            Write (22,c11fmt)                                           &
     &          (q(k,l, element + i), i = 1, 10),                       &
     &          (qcl(k,l, element + i), i = 1, 10),                     &
     &          (qcf(k,l, element + i), i = 1, 10),                     &
     &          (lca(k,l, element + i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           Format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           Format statement via an internal write statement.
!
            Write (ctfmt(16:17), '(i2)') lastrow
            Write (22,ctfmt) (element + i, i = 1, lastrow)
            Write (c2fmt(18:19), '(i2)') lastrow
            Write (c2fmt(6:15), '(''Q kg/kg   '')')
              Write (22,c2fmt) (q(k,l, i + element), i = 1, lastrow)
            Write (c2fmt(6:15), '(''QCL kg/kg '')')
              Write (22,c2fmt)(qcl(k,l, i + element), i = 1, lastrow)
            Write (c2fmt(6:15), '(''QCF kg/kg '')')
              Write (22,c2fmt)(qcf(k,l, i + element), i = 1, lastrow)
            Write (c2fmt(6:15), '(''LCA       '')')
              Write (22,c2fmt)(lca(k,l, i + element), i = 1, lastrow)
          endif
        enddo
          Write (22,c12fmt)
          Write (22,c13fmt) cca(k,l), iccb(k,l), icct(k,l)
          Write (22,c5fmt)
          Write (22,c14fmt)
          Write (22,c5fmt)
        ndeep = nsoilt_levs
!
!       Repeat above section of code for variables on nsoil levels
!       write out variables Tdeep maximum of 10
!       variables per row
!
        If (mod(ndeep,10)  ==  0) then
!
!         Calculate no. of rows and no. of elements in last row
!
          ndeeprows = int(ndeep/10)
          lastrow = 10
        else
          ndeeprows = int(ndeep/10) + 1
          lastrow = mod(ndeep,10)
        endif
        Do ndeepcount = 1, ndeeprows
          element = 10*(ndeepcount-1)
          if (ndeepcount  <   ndeeprows) then
!
!           Write out all complete rows ie of 10 variables per row
!
            Write (22,c8fmt) (element + i + 1/2, i = 1, 10)
            Write (22,c15fmt)                                           &
     &          (t_deep_soil(k,l, element + i), i = 1, 10)
          else
!
!           Write out last row. Use an internal format statement by
!           creating a character string. This will enable a variable
!           Format to be created eg NF10.6 where N is the no. of
!           elements in the last row which can be written into the
!           Format statement via an internal write statement.
!
            Write (ctfmt(16:17), '(i2)') lastrow
            Write (22,ctfmt) (element + i + 1/2, i = 1, lastrow)
            Write (c1fmt(18:19), '(i2)') lastrow
            Write (c1fmt(6:15), '(''Tdeep K   '')')
            Write (22,c1fmt)                                            &
     &          (t_deep_soil(k,l,  i + element), i = 1, lastrow)
          endif
        enddo
          Write (22,c16fmt)
          Write (22,c17fmt)
          Write (22,c18fmt) tstar(k,l), smc(k,l), canopy(k,l),          &
     &      snodep(k,l), zh(k,l), Z0mSea(k,l)
        enddo                   ! k
      enddo                     ! l

      Return
      END SUBROUTINE SUB_DATA
#endif
