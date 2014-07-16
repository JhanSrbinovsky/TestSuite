!*L --------------------- Comdeck: CENVIRDT ---------------------------
!
!   Purpose: Data statements for character enviroment variables used
!            by portable IO to open/close files (links with CENVIR)
!
!   Author : R A Stratton      Date : 22/10/92
!
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  3.4   1/8/94     Revised Obs file spec + CX files: Stuart Bell
!    3.4   18/5/94    Correct misspelling of CURNTIN and change length
!                     to match. J F Thomson
!    4.0  22/9/95     Units for spectral data added. J. M. Edwards
!LL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1 26/02/96     New env. variables for sulphur ancillary files.
!                     Rename SOURCES to SULPEMIS. D. Robinson.
!      4.3  18/3/97   Add aerosol fcgs of climate change. William Ingram
!    4.4  4/7/97    Add ANLINCR at 108. Chris Jones/Stuart Bell
!LL  4.4 12/9/97      New ancillary file environment variables for
!LL                   initial surface type fracs, initial vegetation
!LL                   state and vegetation disturbance.
!LL                                                  R.A.Betts
!    4.4 28/08/97     Move CACHED from logical unit no 3 to 138 in
!                     order to release a Fortran unit [for OASIS].
!                     R.Rawlins
!    4.5 22/04/98     Add new ancillary file for CO2 emissions:
!                     CO2EMITS - in I/O unit 118. Chris Jones.
!    4.5 22/04/98     Add new ancillary file for soot emissions:
!                     SOOTEMIS - in I/O unit 139. R.Rawlins
!    4.5 29/07/98     Move ALABCOU1/2/3/4 from 101-103 to 140-143.
!                     Add ALABCOU5/6/7/8 to 144-147. D. Robinson.
!    4.5 17/08/98     Remove OLABCOUT from Unit 90. Add OLABCOU1/2/3/4
!                     to 101-103. D. Robinson.
!    5.1 13/04/00     TDF_dump replaces FILE107.
!                     IAU_inc  replaces ANLINCR. Adam Clayton
!    5.2 21/08/00     Add an extra op macro for VAR plus 1 user pp outpu
!                     stream. R Rawlins
!    5.2 18/01/01     Add VERT_LEV. D. Robinson
!    5.3 26/10/01     Add LANDFRAC. Remove OBS06-OBS10, CX01-CX10.
!                     D Robinson
!    5.3 24/10/01     Add IDEALISE.                 Andy Malcolm
!    5.3 19/06/2001   Added file (unit no 119) for TPPSOZON. Dave Tan
!    5.4 29/08/02     Add ALABCIN1/2 to unit no 125/126. Remove ALABCIN
!                     from 95. D Robinson
!    5.4 29/08/02     A
!    5.5 17/02/03     Add Wave model interface & boundary files.
!                                                    D.Holmes-Bell
!    5.5 30/12/02     Add DUSTSOIL/BIOMASS to unit no 75/76.
!                     and RIVSTOR/RIVSEQ/RIVDIR to 77-79. D Robinson
!    6.1 07/04/04     Add DMSCONC to unit no 95.   A. Jones
!    6.1 08/11/04     Alter names of River Routing files. R.Sharp
!    6.2 26/01/06     Add ICFILE to unit no 152. T. Edwards
!    6.2 17/03/06     Add SURFEMIS,AIRCREMS,STRATEMS,EXTRAEMS,RADONEMS
!                     to unit nos 130-134. D.Robinson
!    6.2 24/00/05     Change STRATOUT to PPSMC and MESOUT to PPSCREEN,
!                     files for soil moisture nudging scheme. Clive Jones
!LL
!LL DATA statements for COMDECK CENVIR

      DATA FT_ENVIRON/                                                  &
     &  'PPXREF  ','PPXREFU ','        ','STASHCTL','        ',         &
                                                                !  1- 5
     &  '        ','OUTPUT2 ','        ','        ','XHIST   ',         &
     &  'IHIST   ','THIST   ','        ','ERRFLAG ','CACHE1  ',         &
     &  'CACHE2  ','AOTRANS ','ASWAP   ','OSWAP   ','AINITIAL',         &
     &  'ASTART  ','        ','APSUM1  ','APSTMP1 ','        ',         &
     &  '        ','AOMEAN  ','ATMANL  ','        ','OZONE   ',         &
     &  'SMCSNOWD','DSOILTMP','SOILTYPE','VEGTYPE ','SSTIN   ',         &
     &  'SICEIN  ','PERTURB ','CURNTIN ','        ','OINITIAL',         &
     &  'OSTART  ','        ','OPSUM1  ','OPSTMP1 ','        ',         &
     &  '        ','OCNANL  ','ATRACER ','OTRACER ','WFIN    ',         &
     &  'HFLUXIN ','PMEIN   ','ICEFIN  ','AIRTMP  ','SALINITY',         &
     &  'FLUXCORR','SWSPECTD','BAS_IND ','SLABHCON','PP0     ',         &
     &  'PP1     ','PP2     ','PP3     ','PP4     ','PP5     ',         &
     &  'PP6     ','PP7     ','PP8     ','PP9     ','OBS01   ',         &
                                                               !66-70
     &  'OBS02   ','OBS03   ','OBS04   ','OBS05   ','DUSTSOIL',         &
                                                               !71-75
     &  'BIOMASS ','RIVSTOR ','RIVCHAN ','RIVER2A ','LWSPECTD',         &
                                                               !76-80
     &  'WAVEOUT ','SURGEOUT','PPSCREEN','PPSMC   ','WFOUT   ',         &
     &  'HFLUXOUT','PMEOUT  ','ICFOUT  ','MOSOUT  ','VERT_LEV',         &
     &  'SSTOUT  ','SICEOUT ','CURNOUT ','        ','DMSCONC ',         &
                                                               !91-95
     &  'OROG    ','TRANSP  ','OLABCIN ','OCNDEPTH',                    &
     &  'OLABCOU1','OLABCOU2','OLABCOU3','OLABCOU4','FILE104 ',         &
     &  'FILE105 ','IDEALISE','TDF_dump','IAU_inc ','MURKFILE',         &
     &  'SULPEMIS','USRANCIL','USRMULTI','OUSRANCL','OUSRMULT',         &
     &  'SO2NATEM','CHEMOXID','AEROFCG ','CO2EMITS','TPPSOZON',         &
     &  'LANDFRAC','WLABCOU1','WLABCOU2','WLABCOU3','WLABCOU4',         &
                                                               !120-124
     &  'ALABCIN1','ALABCIN2','        ','OCFFEMIS','HORZGRID',         &
                                                               !125-129
     &  'SURFEMIS','AIRCREMS','STRATEMS','EXTRAEMS','RADONEMS',         &
                                                               !130-134
     &  'FRACINIT','VEGINIT ','DISTURB ','CACHED  ','SOOTEMIS',         &
                                                               !135-139
     &  'ALABCOU1','ALABCOU2','ALABCOU3','ALABCOU4','ALABCOU5',         &
                                                               !140-144
     &  'ALABCOU6','ALABCOU7','ALABCOU8','CARIOLO3','RPSEED  ',         & 
                                                               !145-149
     &  'PPVAR   ','PP10    ','ICFILE  ','VAR_GRID','ARCLBIOG',         &
                                                               !150-154
     &  'ARCLBIOM','ARCLBLCK','ARCLSSLT','ARCLSULP','ARCLDUST',         &
                                                               !155-159
     &  'ARCLOCFF','ARCLDLTA','        ','        ','        ',         &
                                                               !160-164
     &          35*'        '                                           &
     & /
!
      DATA LEN_FT_ENVIR/6,7,0,8,0, 0,7,0,0,5,                           &
                                                    !  1-10
     &                  5,5,0,7,6, 6,7,5,5,8,                           &
                                                    ! 11-20
     &                  6,0,6,7,0, 0,6,6,0,5,                           &
                                                    ! 21-30
     &                  8,8,8,7,5, 6,7,7,0,8,                           &
                                                    ! 31-40
     &                  6,0,6,7,0, 0,6,7,7,4,                           &
                                                    ! 41-50
     &                  7,5,6,6,8, 8,8,7,8,3,                           &
                                                    ! 51-60
     &                  3,3,3,3,3, 3,3,3,3,5,                           &
                                                    ! 61-70
     &                  5,5,5,5,8, 7,7,7,7,8,                           &
                                                    ! 71-80
     &                  7,8,8,5,5, 8,6,6,6,8,                           &
                                                    ! 81-90
     &                  6,7,7,0,7, 4,6,7,8,                             &
                                                    ! 91-99
     &                  8,8,8,8,7, 8,8,8,7,8,                           &
                                                    ! 100-109
     &                  8,8,8,8,8, 8,8,7,8,8,                           &
                                                    ! 110-119
     &                  8,8,8,8,8, 8,8,0,8,8,                           &
                                                    ! 120-129
     &                  8,8,8,8,8, 8,7,7,6,8,                           &
                                                    ! 130-139
     &                  8,8,8,8,8, 8,8,8,8,6,                           &
                                                    ! 140-149
     &                  5,4,6,8,8, 8,8,8,8,8,                           &
                                                    ! 150-159
     &                  8,8,0,0,0, 0,0,0,0,0,                           &
                                                    ! 160-169
     &                  30*0/                       ! 170-199

!End of COMDECK CENVIRDT

!
