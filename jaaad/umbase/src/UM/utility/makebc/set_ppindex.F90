#if defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!+ Subroutine SET_PPINDEX : Set up array PPINDEX for MAKEBC
!
! Subroutine Interface :
! New argument list new bc parameters that are not used by
! LAMs at present

       subroutine set_ppindex(jorog,ju,jv,jw,                           &
     &           jrho,jtheta,jq,jqcl,jqcf,                              &
     &           jexner,                                                &
     &           ju_adv,jv_adv,jw_adv,                                  &
     &           jqcf2,jqrain,jqgraup,                                  &
     &           jmurk,                                                 &
     &           jcf_bulk,jcf_liquid,jcf_frozen,                        &
     &           jtracer,                                               &
     &           len_ppindex,ppindex,                                   &
     &           len1_lookup,len2_lookup,lookup,                        &
     &           pplengths,lbcdiag_no,                                  &
     &           d1_pointers,                                           &
#include "argppx.h"
! Time_index and next_dump for searching multiple time
! fields files
     &           next_dump,target_time,                                 &
! Add new calculated_length
     &           calc_length,                                           &
     &           halo_size_out,                                         &
     &           l_pc2,l_murk,max_progs,                                &
     &           l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup)

! Type for time
      USE makebc_time_mod, ONLY: &
        time
 
      IMPLICIT NONE

! Description : Initialise PPINDEX array and pointers.
!
! Method : This routine goes through the lookup table and sets up
!          PPINDEX and the pointers to the data in the dump.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.4    10/10/97  Original Code
!
!   6.0    17/7/03   Code updated to handle UM5.5 LBCs and to generate
!                    extra variables/arrays to allow readflds to be used
!                    instead of um_readdump.
!                     R Sempers (frpz)
!   6.1    17/7/03   Declare include file csubmodl.h before PPINDEX so
!                    that N_INTERNAL_MODEL has been declared. P.Dando
!   6.2    14/2/05   Upgrade to allow murk and PC2 lbcs
!                    R. Sempers
!   6.2    23/11/05   Add calculation of length for D1 array
!   6.2    23/11/05   Allow reading from reinitialised fields
!                     files by selecting data on validity time
!   6.2    20/01/06   Change pplengths_pointers to lbcdiag_no in
!                     code not changed by other 6.2 makebc modsets
!                     to make variable name more decriptive of what
!                     it actually is
!                     Remove arguments comented in vn 6.1

!
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
! Global Variables :

! Subroutine arguments
!   Scalar arguments with intent(in) :

#include "parvars.h"
#include "parlbcs.h"

      Integer LEN_PPINDEX   ! Dimension of PP_INDEX
      Integer LEN1_LOOKUP   ! First dimension of Lookup table
      Integer LEN2_LOOKUP   ! Second dimension of Lookup table


!   Array arguments with intent(in) :

      Integer LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)  ! Lookup table

!   Scalar arguments with intent(inout) :

!   Array arguments with intent(inout) :

!   Scalar arguments with intent(out) :

      integer   jorog      !pointer to orography
      integer   ju         !pointer to u wind component
      integer   jv         !pointer to v wind component
      integer   jw         !pointer to w wind component (vertical)
      integer   jrho       !pointer to density
      integer   jtheta     !pointer to theta
      integer   jq         !pointer to specific humidity
      integer   jqcl       !pointer to qcl liquid water
      integer   jqcf       !pointer to qcf frozen water
      integer   jexner     !pointer to exner pressure
      integer   ju_adv     !pointer to u advection
      integer   jv_adv     !pointer to v advection
      integer   jw_adv     !pointer to w advection
! Add pc2 and murk pointers, also max no of prognostics
      integer   jqcf2      !pointer to qcf2 - cloud ice (cry)
      integer   jqrain     !pointer to qrain - rain
      integer   jqgraup    !pointer to qgraup - graupel horiz.
      integer   jmurk
      integer   jcf_bulk   !pointer to cf_bulk
      integer   jcf_liquid !pointer to cf_liquid
      integer   jcf_frozen !pointer to cf_frozen
      integer   max_progs
      integer   jtracer    !pointer to atm. tracer 1

      logical :: l_pc2
      logical :: l_murk
      logical :: l_mcr_qcf2
      logical :: l_mcr_qrain
      logical :: l_mcr_qgraup
      character*80   cmessage   !error message
      integer :: errorstatus

!   Array arguments with intent(out) :

#include "csubmodl.h"

      Integer PPINDEX(LEN_PPINDEX,n_internal_model)
!  Position of fields in LOOKUP


! allows setppindex to generate extra values required
! by makebc to use readflds to read the dump with um_readdump
      Integer pplengths(LEN_PPINDEX,n_internal_model)
!  Length of fields in LOOKUP

! Set length of arrays to max_progs
      integer lbcdiag_no(max_progs)
      integer d1_pointers(max_progs)

#include "decomptp.h"
#include "decompdb.h"
#include "cppxref.h"
#include "ppxlook.h"

      integer :: loop
      integer :: grid_code
      integer :: fld_type
      integer :: halo_type
! Dimension halo_size_out with max_progs
      integer :: halo_size_out(max_progs,3)
      integer :: field_size_x
      integer :: field_size_y
      integer :: field_size_z
      integer :: field_size
      integer :: pointer_next
      integer :: iproc=0
      integer :: num_levels
      integer :: cumil_length_d1

! Target_time is the lbc time required to select lookups
! pointing to the correct data
      type(time) :: target_time
      logical :: next_dump

! Add integer calc_length
      integer :: calc_length

! Declare functions exppxi and get_fld_type
      integer :: exppxi
      integer :: get_fld_type

!   Local parameters :
      character (len=*), parameter :: RoutineName='Set_PPindex'

!   Local scalars :

      Integer I         ! Loop index

!   Local variables for the stashnumbers to avoid 'magic numbers'
!   will make further modification to the subroutine easier.
      integer, parameter :: stash_ju     = 2
      integer, parameter :: stash_jv     = 3
      integer, parameter :: stash_jtheta = 4
      integer, parameter :: stash_jq     = 10
      integer, parameter :: stash_jqcf   = 12
      integer, parameter :: stash_jorog  = 33
      integer, parameter :: stash_jmurk   = 90
      integer, parameter :: stash_jw     = 150
      integer, parameter :: stash_jrho   = 253
      integer, parameter :: stash_jqcl   = 254
      integer, parameter :: stash_jexner = 255
      integer, parameter :: stash_ju_adv = 256
      integer, parameter :: stash_jv_adv = 257
      integer, parameter :: stash_jw_adv = 258
      integer, parameter :: stash_jcf_bulk  =266
      integer, parameter :: stash_jcf_liquid=267
      integer, parameter :: stash_jcf_frozen=268
      integer, parameter :: stash_jqcf2   = 271
      integer, parameter :: stash_jqrain  = 272
      integer, parameter :: stash_jqgraup = 273

! Logical to say whether data is present
!         for target_time
      logical :: data_present
!-  End of Header

! Initialise calc_length to 0 incase it has not been done 
! before
      calc_length=0

! Initialise d1_pointers to -1
      d1_pointers(:)=-1

! Initialise data_present to .false.
      data_present=.false.
      ju          =-1
      jv          =-1
      jtheta      =-1
      jq          =-1
      jqcf        =-1
      jorog       =-1
      jmurk       =-1
      jw          =-1
      jrho        =-1
      jqcl        =-1
      jexner      =-1
      ju_adv      =-1
      jv_adv      =-1
      jw_adv      =-1
      jcf_bulk    =-1
      jcf_liquid  =-1
      jcf_frozen  =-1
      jqcf2       =-1
      jqrain      =-1
      jqgraup     =-1

      jtracer     =-1

! Initialise pplengths such that all elements are set to 0
      pplengths(:,:)=0

! Initialise array as -1
      lbcdiag_no(:)=-1

      lbcdiag_no(1)=2
      lbcdiag_no(2)=3
      lbcdiag_no(3)=4
      lbcdiag_no(4)=10
      lbcdiag_no(5)=12
      lbcdiag_no(6)=33
! Set values for murk in lbcdiag_no only
! if required
      if(l_murk)then
        lbcdiag_no(7)=90
      endif
      lbcdiag_no(8)=150
      lbcdiag_no(9)=253
      lbcdiag_no(10)=254
      lbcdiag_no(11)=255
      lbcdiag_no(12)=256
      lbcdiag_no(13)=257
      lbcdiag_no(14)=258

! Set values for pc2 in lbcdiag_no only
! if required
      if(l_pc2)then
        lbcdiag_no(15)=266
        lbcdiag_no(16)=267
        lbcdiag_no(17)=268
      endif
! Set values for microphysics in lbcdiag_no only
! if required
      if(l_mcr_qcf2)then
        lbcdiag_no(18)=271
      endif
      if(l_mcr_qrain)then
        lbcdiag_no(19)=272
      endif
      if(l_mcr_qgraup)then
        lbcdiag_no(20)=273
      endif

!Get the starting point in each field required, and also set
!pplengths to the length of each field
      DO I=1,LEN2_LOOKUP

! Test for whether the lookup is for the required time
        if(lookup(1,i) == target_time%year  .and.                       &
     &     lookup(2,i) == target_time%month .and.                       &
     &     lookup(3,i) == target_time%day   .and.                       &
     &     lookup(4,i) == target_time%hour  .and.                       &
     &     lookup(5,i) == target_time%min) then

          data_present=.true.

          IF (LOOKUP(42,I) == stash_jorog) THEN
            pplengths(33,1) = pplengths(33,1)+1
            if(jorog == -1)then
              jorog=LOOKUP(40,I)
              PPINDEX(33,1) = I
          endif



          ELSEIF(LOOKUP(42,I) == stash_ju) THEN
            pplengths(2,1) = pplengths(2,1)+1
            if(ju == -1)then
              ju=LOOKUP(40,I)
              PPINDEX(2,1) = I
            endif

          ELSEIF(LOOKUP(42,I) == stash_jv) THEN
            pplengths(3,1) = pplengths(3,1)+1
            if(jv == -1)then
              jv=LOOKUP(40,I)
              PPINDEX(3,1) = I
            endif


          ELSEIF(LOOKUP(42,I) == stash_jw) THEN
            pplengths(150,1) = pplengths(150,1)+1
            if(jw == -1)then
              jw=LOOKUP(40,I)
              PPINDEX(150,1) = I
            endif



          ELSEIF(LOOKUP(42,I) == stash_jrho) THEN
            pplengths(253,1) = pplengths(253,1)+1
            if(jrho == -1)then
              jrho=LOOKUP(40,I)
              PPINDEX(253,1) = I
            endif

          ELSEIF(LOOKUP(42,I) == stash_jtheta) THEN
            pplengths(4,1) = pplengths(4,1)+1
            if(jtheta == -1)then
              jtheta=LOOKUP(40,I)
              PPINDEX(4,1) = I
            endif



          ELSEIF(LOOKUP(42,I) == stash_jq) THEN
            pplengths(10,1) = pplengths(10,1)+1
            if(jq == -1)then
              jq=LOOKUP(40,I)
              PPINDEX(10,1) = I
            endif



          ELSEIF(LOOKUP(42,I) == stash_jqcl) THEN
            pplengths(254,1) = pplengths(254,1)+1
            if (jqcl == -1)then
              jqcl=LOOKUP(40,I)
              PPINDEX(254,1) = I
            endif

          ELSEIF(LOOKUP(42,I) == stash_jqcf) THEN
            pplengths(12,1) = pplengths(12,1)+1
            if(jqcf == -1)then
              jqcf=LOOKUP(40,I)
              PPINDEX(12,1) = I
            endif

          ELSEIF(LOOKUP(42,I) == stash_jexner) THEN
            pplengths(255,1) = pplengths(255,1)+1
            if(jexner == -1)then
              jexner=LOOKUP(40,I)
              PPINDEX(255,1) = I
            endif

          ELSEIF(LOOKUP(42,I) == stash_ju_adv) THEN
            pplengths(256,1) = pplengths(256,1)+1
            if(ju_adv == -1)then
              ju_adv=LOOKUP(40,I)
              PPINDEX(256,1) = I
            endif


          ELSEIF(LOOKUP(42,I) == stash_jv_adv) THEN
            pplengths(257,1) = pplengths(257,1)+1
            if(jv_adv == -1)then
              jv_adv=LOOKUP(40,I)
              PPINDEX(257,1) = I
            endif




        ELSEIF(LOOKUP(42,I) == stash_jcf_bulk) THEN
          pplengths(266,1) = pplengths(266,1)+1
          if(jcf_bulk == -1)then
            jcf_bulk=LOOKUP(40,I)
            PPINDEX(266,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jcf_liquid) THEN
          pplengths(267,1) = pplengths(267,1)+1
          if(jcf_liquid == -1)then
            jcf_liquid=LOOKUP(40,I)
            PPINDEX(267,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jcf_frozen) THEN
          pplengths(268,1) = pplengths(268,1)+1
          if(jcf_frozen == -1)then
            jcf_frozen=LOOKUP(40,I)
            PPINDEX(268,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jqcf2) THEN
          pplengths(271,1) = pplengths(271,1)+1
          if(jqcf2 == -1)then
            jqcf2=LOOKUP(40,I)
            PPINDEX(271,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jqrain) THEN
          pplengths(272,1) = pplengths(272,1)+1
          if(jqrain == -1)then
            jqrain=LOOKUP(40,I)
            PPINDEX(272,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jqgraup) THEN
          pplengths(273,1) = pplengths(273,1)+1
          if(jqgraup == -1)then
            jqgraup=LOOKUP(40,I)
            PPINDEX(273,1) = I
          endif

        ELSEIF(LOOKUP(42,I) == stash_jmurk) THEN
          pplengths(90,1) = pplengths(90,1)+1
          if(jmurk == -1)then
            jmurk=LOOKUP(40,I)
            PPINDEX(90,1) = I
          endif

          ELSEIF(LOOKUP(42,I) == stash_jw_adv) THEN
            pplengths(258,1) = pplengths(258,1)+1
            if(jw_adv == -1)then
              jw_adv=LOOKUP(40,I)
              PPINDEX(258,1) = I
            endif

          ELSEIF((LOOKUP(42,I)  >=  61).AND.(LOOKUP(42,I)  <=  89))THEN
            pplengths(LOOKUP(42,I),1) = pplengths(LOOKUP(42,I),1)+1
            IF (JTRACER == -1) THEN
              JTRACER=LOOKUP(40,I)
              PPINDEX(LOOKUP(42,I),1) = I
              lbcdiag_no(14)=LOOKUP(42,I)
            ENDIF


          ENDIF
        endif !test on target_time
      ENDDO


! If no data found for this time, exit loop
      if(data_present)then
        next_dump=.false.
! Initialise variables prior to entering the loop
      pointer_next=1
      cumil_length_d1=0


! Loop over the maximum number of prognostic fields
      do loop=1,max_progs

! Check that the prognostic has been found before working on it
        if(lbcdiag_no(loop) /= -1)then

! Use exppxi to get the field/grid type from the ppx array
! DEPENDS ON: exppxi
          grid_code=exppxi(1,0,lbcdiag_no(loop),ppx_grid_type,          &
#include "argppx.h"
     &                  errorstatus,cmessage)
          if(errorstatus /= 0)then
! DEPENDS ON: ereport
            call ereport(routinename,errorstatus,cmessage)
          endif
! DEPENDS ON: exppxi
          halo_type=exppxi(1,0,lbcdiag_no(loop),ppx_halo_type,          &
#include "argppx.h"
     &                  errorstatus,cmessage)
          if(errorstatus /= 0)then
! DEPENDS ON: ereport
            call ereport(routinename,errorstatus,cmessage)
          endif
          halo_size_out(loop,1)=decomp_db_halosize                      &
     &                            (1,halo_type,decomp_smexe)
          halo_size_out(loop,2)=decomp_db_halosize                      &
     &                            (2,halo_type,decomp_smexe)
          halo_size_out(loop,3)=decomp_db_halosize                      &
     &                            (3,halo_type,decomp_smexe)

! use get_fld_type to work out which field type it is
! DEPENDS ON: get_fld_type
          fld_type=GET_FLD_TYPE (grid_code)

! try to get the number of levels in the field
          num_levels=lookup(17,lbcdiag_no(loop))-100

          field_size_x=decomp_db_g_blsize(1,fld_type,iproc,decomp_smexe)&
     &               +(2*decomp_db_halosize(1,halo_type,decomp_smexe))
          field_size_y=decomp_db_g_blsize(2,fld_type,iproc,decomp_smexe)&
     &               +(2*decomp_db_halosize(2,halo_type,decomp_smexe))
          field_size_z=pplengths(lbcdiag_no(loop),1)                    &
     &               +(2*decomp_db_halosize(3,halo_type,decomp_smexe))
          field_size=field_size_x*field_size_y*field_size_z

          d1_pointers(loop)=cumil_length_d1+1
          cumil_length_d1=cumil_length_d1 + field_size

! Work out the length required to store each prognostic and add
! up to create a cumilative total
          calc_length=calc_length+field_size




        endif
      enddo


       ju         = d1_pointers(1)
       jv         = d1_pointers(2)
       jtheta     = d1_pointers(3)
       jq         = d1_pointers(4)
       jqcf       = d1_pointers(5)
       jorog      = d1_pointers(6)
       jmurk      = d1_pointers(7)
       jw         = d1_pointers(8)
       jrho       = d1_pointers(9)
       jqcl       = d1_pointers(10)
       jexner     = d1_pointers(11)
       ju_adv     = d1_pointers(12)
       jv_adv     = d1_pointers(13)
       jw_adv     = d1_pointers(14)
       jcf_bulk   = d1_pointers(15)
       jcf_liquid = d1_pointers(16)
       jcf_frozen = d1_pointers(17)
       jqcf2      = d1_pointers(18)
       jqrain     = d1_pointers(19)
       jqgraup    = d1_pointers(20)



      IF (jorog  == -1 .or.                                             &
     &    ju     == -1 .or.                                             &
     &    jv     == -1 .or.                                             &
     &    jw     == -1 .or.                                             &
     &    jrho   == -1 .or.                                             &
     &    jtheta == -1 .or.                                             &
     &    jq     == -1 .or.                                             &
     &    jqcl   == -1 .or.                                             &
     &    jqcf   == -1 .or.                                             &
     &    jexner == -1 .or.                                             &
     &    ju_adv == -1 .or.                                             &
     &    jv_adv == -1 .or.                                             &
     &    jw_adv == -1                                                  &
! Check for pc2, microphysics and murk only if they're required
     &    .or. (jqcf2       == -1 .and. l_mcr_qcf2)                     &
     &    .or. (jqrain      == -1 .and. l_mcr_qrain)                    &
     &    .or. (jqgraup     == -1 .and. l_mcr_qgraup)                   &
     &    .or. (jmurk       == -1 .and. l_murk)                         &
     &    .or. (jcf_bulk    == -1 .and. l_pc2)                          &
     &    .or. (jcf_liquid    == -1 .and. l_pc2)                        &
     &    .or. (jcf_frozen    == -1 .and. l_pc2)                        &
     & ) THEN
        errorstatus = 1
        write (6,*) '  Data missing from input file.'
        CMESSAGE  = '  Data missing from input file.'
        WRITE (6,*) ' jorog = ',jorog
        WRITE (6,*) ' ju=     ',ju
        WRITE (6,*) ' jv=     ',jv
        WRITE (6,*) ' jw=     ',jw
        WRITE (6,*) ' jrho=   ',jrho
        WRITE (6,*) ' jtheta= ',jtheta
        write(6,*) ' jq=      ',jq
        write(6,*) ' jqcl=    ',jqcl
        write(6,*) ' jqcf=    ',jqcf
        write(6,*) ' jexner=  ',jexner
        write(6,*) ' ju_adv=  ',ju_adv
        write(6,*) ' jv_adv=  ',jv_adv
        write(6,*) ' jw_adv=  ',jw_adv
        write(6,*) ' jqcf2=   ',jqcf2
        write(6,*) ' jqrain=  ',jqrain
        write(6,*) ' jqgraup= ',jqgraup
        write(6,*) ' jmurk=   ',jmurk
        write(6,*) ' jcf_bulk=',jcf_bulk
        write(6,*) ' jcf_liquid=',jcf_liquid
        write(6,*) ' jcf_frozen=',jcf_frozen

! DEPENDS ON: ereport
        call ereport(routinename,errorstatus,cmessage)
      ENDIF


!     If no Tracers, reset JTRACER to prevent negative pointer.
      IF (JTRACER  ==  -1) THEN
        WRITE (6,*) ' '
        WRITE (6,*) ' JTRACER = ',JTRACER
        WRITE (6,*) ' No tracer data in this dump.'
        JTRACER = 1
      ENDIF
      else
        next_dump=.true.
      endif  ! if data_present = .true.

! Avoid negative pointers if pc2,microphysics or murk not included
      if(.not. l_murk)then
        write(6,*)'No murk lbcs'
        jmurk =1
        lbcdiag_no(7)=-1
      endif
      if(.not. l_pc2)then
        write(6,*)'No pc2 lbcs'
        jcf_bulk=1
        jcf_liquid=1
        jcf_frozen=1
        lbcdiag_no(15:17)=-1
      endif
      if(.not. l_mcr_qcf2)then
        jqcf2  =1
        lbcdiag_no(18)=-1
      endif
      if(.not. l_mcr_qrain)then
        jqrain =1
        lbcdiag_no(19)=-1
      endif
      if(.not. l_mcr_qgraup)then
        jqgraup=1
        lbcdiag_no(20)=-1
      endif

      RETURN
      END SUBROUTINE set_ppindex
#endif
