#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE IN_STOCHEM_INTF
! Module to contain interfaces
!VVV  V5.0 MODULE IN_STOCHEM_INTF 12/VII/01  ND version
!
! Modification history from model version 6.0:
!  Model
! version  Date     Comment
!
!  6.0   17/09/03   Correct non-standard F90 continuation lines.
!                   Introduce standard UM modification history.
!                                                         P.Dando
!  6.1   20/08/04   Changes to allow vectorisation on SX6.
!                                                   M. Sanderson.
!  6.2   15/11/05  Extra variables added. M.G. Sanderson
!----------------------------------------------------------------------

      INTERFACE
        SUBROUTINE FINDNAME(filetype,filetype_2,submodel_id,prev_h,     &
     &    prev_d,filename,month,year,hour,day)
          INTEGER,       INTENT(IN)    :: prev_h
          INTEGER,       INTENT(IN)    :: prev_d
          INTEGER, INTENT(IN),OPTIONAL :: month       ! month
          INTEGER, INTENT(IN),OPTIONAL :: year        ! year
          INTEGER, INTENT(IN),OPTIONAL :: hour        ! hour
          INTEGER, INTENT(IN),OPTIONAL :: day         ! day
          CHARACTER*1,   INTENT(IN)    :: filetype
          CHARACTER*1,   INTENT(IN)    :: filetype_2
          CHARACTER*1,   INTENT(IN)    :: submodel_id
          CHARACTER*14,  INTENT(OUT)   :: filename
        END SUBROUTINE FINDNAME
      END INTERFACE

      INTERFACE
! DEPENDS ON: st_height
        INTEGER FUNCTION ST_HEIGHT(pos,eta_array)
          CHARACTER(*), INTENT(IN) :: eta_array
          REAL, INTENT(IN) :: pos
        END FUNCTION ST_HEIGHT

! DEPENDS ON: gauss_ran
        FUNCTION GAUSS_RAN(seed,ransize,nfill)
          INTEGER, INTENT(IN) :: ransize
          INTEGER, INTENT(IN) :: nfill
          INTEGER, DIMENSION(ransize), INTENT(IN) :: seed
          REAL, DIMENSION(nfill) :: GAUSS_RAN
        END FUNCTION GAUSS_RAN

! DEPENDS ON: interp_s
        REAL FUNCTION INTERP_S(pos,x,lhalf,lvert,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                     &
     &                        INTENT(IN) :: x
          REAL, DIMENSION(4), INTENT(IN) :: pos      ! Position
          LOGICAL,            INTENT(IN) :: lhalf    ! Horiz=1/2 grids
          LOGICAL,            INTENT(IN) :: lvert    ! .T. for U grid
          LOGICAL, OPTIONAL,  INTENT(IN) :: lu ! .T. for U grid, reqd
                                                     ! only if LHALF=.F.
        END FUNCTION INTERP_S

! DEPENDS ON: hinterp
        REAL FUNCTION HINTERP(x,d1,d2,i,l)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(lnbound:lnbound+nlonpe-1,                     &
     &      lobound:lobound+nlatpe-1), INTENT(IN) :: x
          REAL,    INTENT(IN) :: d1            ! Position
          REAL,    INTENT(IN) :: d2            ! Position
          INTEGER, INTENT(IN) :: i             ! Index
          INTEGER, INTENT(IN) :: l             ! Index
        END FUNCTION HINTERP

! DEPENDS ON: imet
        INTEGER FUNCTION IMET(pos1,lhalf)
          USE IN_STOCHEM_GRD
          REAL,    INTENT(IN) :: pos1
          LOGICAL, INTENT(IN) :: lhalf
        END FUNCTION IMET

! DEPENDS ON: jmet
        INTEGER FUNCTION JMET(pos2,lhalf)
          USE IN_STOCHEM_GRD
          REAL,    INTENT(IN) :: pos2
          LOGICAL, INTENT(IN) :: lhalf
        END FUNCTION JMET

! DEPENDS ON: zen
        REAL FUNCTION ZEN(time,xlong,xlat)
          USE IN_STOCHEM_GRD
          REAL, INTENT(IN) :: time
          REAL, INTENT(IN) :: xlong,xlat
        END FUNCTION ZEN

! DEPENDS ON: secs
        REAL FUNCTION SECS(day,imonth,iyear)
          REAL, INTENT(IN)    :: day
          INTEGER, INTENT(IN) :: imonth,iyear
        END FUNCTION SECS

        REAL FUNCTION LINTERP(pos,x,lhalf,lvert,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev),                     &
     &                    INTENT(IN) :: x        ! Interpolation array
          REAL, DIMENSION(4), INTENT(IN) :: pos      ! Position
          LOGICAL,            INTENT(IN) :: lhalf    ! Horiz=1/2 grids
          LOGICAL,            INTENT(IN) :: lvert    ! .T. for U grid
          LOGICAL, OPTIONAL,  INTENT(IN) :: lu       ! .T. for U grid
        END FUNCTION LINTERP

! DEPENDS ON: etator
        REAL FUNCTION ETATOR(eta,orog,z_top_of_model,                   &
     &    first_constant_r_rho_level)
          USE IN_STOCHEM_GRD
          INTEGER, INTENT(IN) :: first_constant_r_rho_level
          REAL, INTENT(IN)    :: eta
          REAL, INTENT(IN)    :: orog
          REAL, INTENT(IN)    :: z_top_of_model
        END FUNCTION ETATOR

! DEPENDS ON: rtoeta
        REAL FUNCTION RTOETA(r,orog,z_top_of_model,                     &
     &    first_constant_r_rho_level)
          USE IN_STOCHEM_GRD
          INTEGER, INTENT(IN) :: first_constant_r_rho_level
          REAL, INTENT(IN)    :: r
          REAL, INTENT(IN)    :: orog
          REAL, INTENT(IN)    :: z_top_of_model
        END FUNCTION RTOETA

! DEPENDS ON: getmetpoint
        REAL FUNCTION GETMETPOINT(pos,field,Lhalf,lu)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos
          REAL, DIMENSION(nlonpe,nlatpe), INTENT(IN) :: field
          LOGICAL,                        INTENT(IN) :: Lhalf ! T 1/2
          LOGICAL,              OPTIONAL, INTENT(IN) :: lu     ! T for U
        END FUNCTION GETMETPOINT

        SUBROUTINE MET2DATA(yy,x,nlevelsin,nlevelsout,lvert)
          USE IN_STOCHEM_GRD
          INTEGER, INTENT(IN) :: nlevelsin
          INTEGER, INTENT(IN) :: nlevelsout
          LOGICAL, OPTIONAL, INTENT(IN) :: lvert ! True for rho levels
          REAL, DIMENSION(nlonpe,nlatpe,nlevelsin), INTENT(IN) :: x
          REAL, DIMENSION(nlnpe,nlpe,nlevelsout),  INTENT(OUT) :: yy
        END SUBROUTINE MET2DATA

! DEPENDS ON: eta2p
        REAL FUNCTION ETA2P(pos,lnp)
          USE IN_STOCHEM_GRD
          REAL, DIMENSION(4),             INTENT(IN) :: pos ! Pos
          REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: lnp
        END FUNCTION ETA2P

! DEPENDS ON: p2eta
        REAL FUNCTION P2ETA(press,pos,p)
          USE IN_STOCHEM_GRD
          REAL,                                     INTENT(IN) :: press
          REAL, DIMENSION(4),             INTENT(IN) :: pos ! Pos
          REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: P
        END FUNCTION P2ETA
      END INTERFACE

      END MODULE IN_STOCHEM_INTF
#endif
