! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
#if !defined(RECON) && !defined(UTILIO) && !defined(FLDOP)             \
 && !defined(VAROPSVER)
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
#elif !defined(VAROPSVER)
     &  PPXI,PPXC,ppxRecs,                                              &
#endif
! End of comdeck
