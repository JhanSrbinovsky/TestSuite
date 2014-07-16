#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE INIBOUNDS2
!----------------------------------------------------------------------
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.4    24/09/97  Created.  W.J. Collins
!  5.5    17/07/02  Now uses rectangular grid. W.J. Collins
!  5.5    21/08/02  Data arrays don't cover full longitude. W.J. Collins
!  5.5    13/02/04  Calls HEIGHT_INI to initialise height arrays.
!                   K. Ketelsen.
!  6.1    20/10/04  Minor tidying of code. M.G. Sanderson
!  6.2    20/04/05  Modifications so will work on any PE configuration.
!                   R. Johanni.
!  6.2    06/03/06  Change GSYNC to SSYNC. P.Selwood.
!
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      INTEGER :: info, lndat_next, ltdat_next

      ncol = nproc_x           ! Number of longitude columns of pes
      numrows = nproc_y        ! Number of latitude rows of pes
      nclprc = 2*(ncell/nproc) ! Number of cells per pe

! Number of Met. latitude points per PE. Set to maximum possible
      nlatpe = MAXVAL(g_lasize(2,fld_type_p,halo_type_extended,:))

! Number of Met. longitude points per PE. Same for all fields
      nlonpe = lasize(1,fld_type_p,halo_type_extended)

! Number of Met. latitude points in use on this PE
      rowspe = lasize(2,fld_type_p,halo_type_extended)

      lnbound = datastart(1)-halo_i
      lobound = datastart(2)-halo_j

! Set the decomposition of the STOCHEM grid like the one
! of the Met grid with scaling applied

! First grid point on every PE
      lndat = ((datastart(1)-1)*nlong)/nmetlong+1
      ltdat = ((datastart(2)-1)*mnlat)/(nmetlat-1)+1

! Longitudes/PE
      IF (gridpos(1) == nproc_x-1) THEN
        nlnpe = nlong-lndat+1
      ELSE
        lndat_next = ((g_datastart(1,mype+1)-1)*nlong)/nmetlong+1
        nlnpe = lndat_next-lndat
      END IF

! Latitudes/PE
      IF (gridpos(2) == nproc_y-1) THEN
         nlpe = mnlat-ltdat+1
      ELSE
        ltdat_next =((g_datastart(2,mype+nproc_x)-1)*mnlat) /           &
     &    (nmetlat-1)+1
        nlpe = ltdat_next-ltdat
      END IF

! Get maximum number of gridpoints per PE

      max_stochem_points = nlnpe*nlpe
      CALL GC_SSYNC(nproc,info)
      CALL GC_IMAX(1,nproc,info,max_stochem_points)
      CALL GC_SSYNC(nproc,info)

      WRITE(6,*)
      WRITE(6,*) 'STOCHEM grid (PE',mype,') set to:'
      WRITE(6,*) 'Longitudes/PE   (nlnpe): ',nlnpe
      WRITE(6,*) 'Latitudes/PE    (nlpe) : ',nlpe
      WRITE(6,*) 'First Longitude (lndat): ',lndat
      WRITE(6,*) 'First Latitude  (ltdat): ',ltdat
      WRITE(6,*) 'max_stochem_points:      ',max_stochem_points
      WRITE(6,*) 'UM grid:'
      WRITE(6,*) 'Longitudes/PE (nlonpe) : ',nlonpe
      WRITE(6,*) 'Latitudes/PE  (nlatpe) : ',nlatpe
      WRITE(6,*) 'Latitudes/PE used (rowspe) : ',rowspe
      WRITE(6,*)

      procmap(lndat:lndat+nlnpe-1,ltdat:ltdat+nlpe-1) = mype
! sum procmap over all processors to set all elements.
      CALL GC_SSYNC(nproc,info)
      CALL GC_ISUM(nlong*mnlat,nproc,info,procmap)
      CALL GC_SSYNC(nproc,info)

      END SUBROUTINE INIBOUNDS2
#endif
