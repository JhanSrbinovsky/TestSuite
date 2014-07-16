! ARTATCPL start
! Description: Super Arrays containing gridline coordinates for
! interpolation and area-averaging between atmosphere and
! river-routing grids (Part of ARTAOCPL.hdk)
! Author: C.Bunton 28.02.03
!
! History
! Version  Date    Comment
!  5.5  28/02/03  Original code. C.Bunton
#if defined(ATMOS)
!L --------------- (Atmosphere-TRIP) coupling arrays  -----------
!L ---------------Lat., Long. values of Atmosphere --------------
     &AO_SPCPL(AO_IXCPL( 5)),AO_SPCPL(AO_IXCPL( 6)),                    &
     &AO_SPCPL(AO_IXCPL( 7)),AO_SPCPL(AO_IXCPL( 8)),                    &
     &AO_SPCPL(AO_IXCPL( 9)),AO_SPCPL(AO_IXCPL( 10)),                   &
! END ARTATCPL
#endif
