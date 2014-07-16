#if defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define identifiers for surface tiles.
!
MODULE tileid3z
!
! Description:
!   This module defines identifiers for different surface types
!   as used in the radiation scheme.
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   13/02/06   Original code included into UM Build
!                        (J.-C. Thelen)
!
! Code Description:
!   Language: FORTRAN 90
!
!- End of header
!
!
!
  INTEGER, Parameter :: npd_tile_type  = 3
!   Number of types of tile for which space is allocated
  INTEGER, Parameter :: IP_ocean_tile  = 1
!   Identifier for open sea
  INTEGER, Parameter :: IP_seaice_tile = 2
!   Idenitifer for ice
  INTEGER, Parameter :: IP_land_tile   = 3
!   Identifer for land
!
!
!
END MODULE tileid3z
#endif
