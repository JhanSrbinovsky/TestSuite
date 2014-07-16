#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Sets up lookup-table for boundary file
!
! Subroutine Interface:

      SUBROUTINE LBC_SetUp_LookUp (                                     &
     &           lbc_lookUp,                                            &
     &           len1_lbc_lookup,                                       &
     &           n_lbc_vars,                                            &
     &           lbc_item_codes,                                        &
     &           jintf,                                                 &
     &           ntime,                                                 &
     &           lbc_fixhd,                                             &
#include "argppx.h"
     &           lbc_rim_size,                                          &
     &           lbc_rim_type,                                          &
     &           lbc_halo_type,                                         &
     &           lbc_fld_type,                                          &
     &           intf_halosize,                                         &
     &           lbc_levels                                             &
     & )

      Implicit None
!
! Description:
!   Set up the lookup headers for lbc variables for one area.
!   Called for each validity time of lbc data.
!
! Method:
!   Set up lookup-table as defined in UMDP F3.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    23/10/01   Correct section number to get packing indicator
!                     from STASHmaster record. Dave Robinson
!   5.3    29/01/02   Set LBSRCE to UM version ID xxxxyyyy,
!                     where xxxx is UM version number e.g. 0503
!                     and yyyy is the model identifier e.g. 1111 for
!                     unified model. D.M. Goddard
!   5.4    09/05/02   Correct initialisation of disk/start addresses in
!                     LookUp for CRUNs. D.Robinson.
!   5.5    24/10/02   Set up grid in Real part of Lookup. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

#include "parvars.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "clookadd.h"
#include "cintfa.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cntl_io.h"
#include "parlbcs.h"
#include "c_model_id.h"

! Subroutine arguments

      Integer  :: len1_lbc_lookup
      Integer  :: n_lbc_vars
      Integer  :: lbc_lookup(len1_lbc_Lookup, n_lbc_vars)
      Integer  :: jintf
      Integer  :: ntime
      Integer  :: lbc_fixhd(256)
      Integer  :: lbc_item_codes(n_lbc_vars)
      Integer  :: lbc_rim_size  (n_lbc_vars)
      Integer  :: lbc_rim_type  (n_lbc_vars)
      Integer  :: lbc_halo_type (n_lbc_vars)
      Integer  :: lbc_fld_type  (n_lbc_vars)
      Integer  :: lbc_levels    (n_lbc_vars)
      Integer  :: intf_halosize (2,Nhalo_max)

! Local parameters:

      Integer,           Parameter :: Sect31   = 31
      Integer,           Parameter :: Sect32   = 32
      Character (Len=*), Parameter :: RoutineName= 'LBC_SetUp_LookUp'

!     Local variables

      Integer :: Start_Address   !  w.r.t start of data
      Integer :: Disk_Address    !  w.r.t start of file
      Integer :: Len_data        !  actual length of data
      Integer :: Disk_length     !  length of data including rounding
                                 !  to sector boundary
      Integer :: N1,N2,N3,N4,N5  !  components of packing indicator
      Integer :: R1R0            !  lbc rimwidth
      Integer :: Y1Y0            !  lbc halo - NS
      Integer :: X1X0            !  lbc halo - EW
      Integer :: Var, Var1       !  loop indexes
      Integer :: lbc_item        !  Item Code

      Real    :: Lat0, Lon0      ! Zeroth Lat & Long

      Integer :: ErrorStatus          !  error code
      Character (Len=80) :: Cmessage  !  error message

! Local dynamic arrays:

!     The lookup table is not stored for previous LBC data times so
!     retain disk and start addresses to avoid unnecessary i/o.


! Function & Subroutine calls:
      Integer Exppxi, Get_um_version_id

!- End of header



!     Initialise disk address and start address

!     - If LBCs already exist in the LBC file then LBC_LOOKUP on entry
!     contains the lookup headers for the last batch of LBCs generated.

!     - At the start of CRUNs LBC_LOOKUP will contain the lookup header
!     corresponding to the last LBC variable before the CRUN restart
!     position. (Lookup Header read in in IN_INTF)

      If (ntime == 1) Then
        disk_address  = lbc_fixhd(160)-1
        start_address = 1
      Else
        disk_address  = LBC_Lookup(lbegin,n_lbc_vars) +                 &
     &                  LBC_Lookup(lbnrec,n_lbc_vars)
        start_address = LBC_Lookup(naddr ,n_lbc_vars) +                 &
     &                  LBC_Lookup(lblrec,n_lbc_vars)
      End If

      If (ntime == 1) Then  !  Include Orography
        Var1 = 1
      Else                  !  Exclude Orography
        Var1 = 2
      Endif

      LBC_Lookup(:,:) = 0

      Do Var = Var1, n_lbc_vars

        lbc_item = lbc_item_codes(var) - 32000

!       Determine length of data for this lbc variable

        len_data =                                                      &
     &  lbc_global_lenrima(lbc_fld_type(var),lbc_halo_type(var) ) *     &
     &  lbc_levels(var)

!       Validity Time
        LBC_Lookup (lbyr ,var) = LBC_VT_Year
        LBC_Lookup (lbmon,var) = LBC_VT_Month
        LBC_Lookup (lbdat,var) = LBC_VT_Day
        LBC_Lookup (lbhr ,var) = LBC_VT_Hour
        LBC_Lookup (lbmin,var) = LBC_VT_Min
        LBC_Lookup (lbday,var) = LBC_VT_DayNo

!       Data Time
        LBC_Lookup (lbyrd ,var) = LBC_DT_Year
        LBC_Lookup (lbmond,var) = LBC_DT_Month
        LBC_Lookup (lbdatd,var) = LBC_DT_Day
        LBC_Lookup (lbhrd ,var) = LBC_DT_Hour
        LBC_Lookup (lbmind,var) = LBC_DT_Min
        LBC_Lookup (lbdayd,var) = LBC_DT_DayNo

!       Length of data
        LBC_Lookup (lblrec,var) = LEN_DATA

!       Grid Type - Rotated regular lat/long
        LBC_Lookup (lbcode,var) = 100+1

!       Hemisphere Indicator : Not used for LBC files
!       Use to store no of lbc levels : 100+(no of levels)
        LBC_Lookup (lbhem,var) = 100+lbc_levels(var)

!       No of rows - not including haloes
        If (lbc_fld_type(var) == fld_type_v) Then  ! v-comp
          LBC_Lookup (lbrow,var) = intf_p_rows(jintf)-1
        Else
          LBC_Lookup (lbrow,var) = intf_p_rows(jintf)
        End If

!       No of points - not including haloes
        If (lbc_fld_type(var) == fld_type_u) Then  ! u-comp
          LBC_Lookup (lbnpt,var) = intf_row_length(jintf)-1
        Else
          LBC_Lookup (lbnpt,var) = intf_row_length(jintf)
        End If


!       Packing Indicator
        Select Case ( intf_pack(jintf) )

          Case (0)   !  No packing
            N1 = 0

          Case (1)   !  32-bit packing
            N1 = 2

          Case (2)   !  Pack according to Stashmaster Record
! DEPENDS ON: exppxi
            N1 = Exppxi(                                                &
     &           atmos_im,Sect32,lbc_item,ppx_dump_packing,             &
#include "argppx.h"
     &           ErrorStatus, cmessage)

        End Select

        N2 = 0   !  Data not compressed
        N3 = 0   !  Compression definition
        N4 = 0   !  Number format
        N5 = 0   !  Not used
        LBC_Lookup (lbpack,var) =                                       &
     &  N5*10000 + N4*1000 +N3*100 + N2*10 + N1

        LBC_Lookup (lbrel,var) = 2

!       Disk Address
        LBC_Lookup (lbegin,var) = disk_address

!       Disk Length

!       Get length of data
        if (mod(lbc_lookup(lbpack, var), 10) == 2) then
          disk_length= (LBC_Lookup(lblrec,var)+1)/2
        else
          disk_length = LBC_Lookup(lblrec,var)
        endif

!       round up to sector boundary
        disk_length=((disk_length+um_sector_size-1)/                    &
     &               um_sector_size)*um_sector_size

!       Disk length : includes rounding-up to sector boundary
        LBC_Lookup (lbnrec,var) = disk_length

!       Level code : 7777 for multi-level field.
        LBC_Lookup (lblev,var) = 7777

! DEPENDS ON: get_um_version_id
        LBC_Lookup (lbsrce,var) = get_um_version_id(model_id)

!       Data Type
        LBC_Lookup (data_type,var) = 1

!       Start Address
        LBC_Lookup (naddr,var) = Start_Address

!       LBUSER(3) stores rimwidth and haloes as R1R0Y1Y0X1X0
!       R1R0 : Rimwidth, Y1Y0 : halo_y, X1X0 : halo_x

        R1R0 = lbc_rim_size(lbc_rim_type(var))
        Y1Y0 = Intf_HaloSize(2,lbc_halo_type(var))
        X1X0 = Intf_HaloSize(1,lbc_halo_type(var))
        LBC_Lookup (lbuser3,var) = R1R0 * 10000 + Y1Y0 * 100 + X1X0

!       Item Code
        LBC_Lookup (item_code,var) = Sect31*1000 + lbc_item

!       Model Code
        LBC_Lookup (model_code,var) = atmos_im

!       REAL section of lookup.

!       Lat & Long of rotated pole
        LBC_Lookup (bplat,var) = Transfer (intf_polelat (jintf),1)
        LBC_Lookup (bplon,var) = Transfer (intf_polelong(jintf),1)

!       Zeroth latitude
        If (lbc_fld_type(var) == fld_type_v) Then
          Lat0 = intf_firstlat(jintf) - 0.5*intf_nsspace(jintf)
        Else
          Lat0 = intf_firstlat(jintf) - intf_nsspace (jintf)
        End If
        LBC_Lookup (bzy,var) = Transfer (Lat0, 1)

!       Latitude Interval
        LBC_Lookup (bdy,var) = Transfer (intf_nsspace(jintf),1)

!       Zeroth longitude
        If (lbc_fld_type(var) == fld_type_u) Then
          Lon0 = intf_firstlong(jintf) - 0.5*intf_ewspace(jintf)
        Else
          Lon0 = intf_firstlong(jintf) - intf_ewspace (jintf)
        End If
        LBC_Lookup (bzx,var) = Transfer (Lon0, 1)

!       Longitude Interval
        LBC_Lookup (bdx,var) = Transfer (intf_ewspace(jintf), 1)

!       Update disk and start address
        disk_address  = disk_address  + LBC_Lookup(lbnrec,var)
        start_address = start_address + LBC_Lookup(lblrec,var)

      End Do
      
      If (intf_ewspace(jintf) < 0) Then
        Do Var = Var1, n_lbc_vars
          LBC_Lookup (bzy,var) = Transfer(intf_ewspace(jintf), 1)
          LBC_Lookup (bzx,var) = Transfer(intf_ewspace(jintf), 1)
        End Do
      End If
      
      return
      END SUBROUTINE LBC_SetUp_LookUp
#endif
