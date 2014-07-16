! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for MPP.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
#if defined(ATMOS) && defined(OCEAN)
     &  SPD1(IXD1( 3)), SPD1(IXD1( 4)),                                 &
#endif
