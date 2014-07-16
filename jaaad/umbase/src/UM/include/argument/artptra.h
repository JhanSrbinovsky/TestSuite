! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
#if defined(ATMOS)
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
#endif
