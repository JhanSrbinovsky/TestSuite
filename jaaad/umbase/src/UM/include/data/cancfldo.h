!*L--------------------COMDECK  CANCFLDO ---------------------------
!
! Purpose : List of Ancillary Fields - Ocean
!           Stash Codes and Logical file Numbers
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.4   21/03/97  New comdeck. D. Robinson.
!
! -------------------------------------------------------------------
! Type Declarations

      INTEGER                                                           &
     &  ITEM_CODES_ANCIL(NANCIL_FIELDS)                                 &
                                          ! Stash Codes
     &, ANCIL_FILE_NO(NANCIL_FIELDS)      ! Logical file numbers

      DATA ITEM_CODES_ANCIL /                                           &
     &  30, 199, 150, 151, 152, 153, 161, 162, 165, 166,                &
                                                            !  1-10
     & 167, 170, 171, 172, 190, 191, 180, 181, 182, 183,                &
                                                            ! 11-20
     & 185, 186, 187, 331, 332, 333, 334, 335, 336, 337,                &
                                                            ! 21-30
     & 338, 339, 340, 351, 352, 353, 354                                &
                                                            ! 31-40
     & /

      DATA ANCIL_FILE_NO /                                              &
     &   7,   8,   1,   1,   1,   1,   2,   2,   3,   3,                &
                                                            !  1-10
     &   3,   4,   4,   4,   4,   4,   5,   5,   5,   5,                &
                                                            ! 11-20
     &   6,   6,   6,   9,   9,   9,   9,   9,   9,   9,                &
                                                            ! 21-30
     &   9,   9,   9,  10,  10,  10,  10                                &
                                                            ! 31-40
     & /
