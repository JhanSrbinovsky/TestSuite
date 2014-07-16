module setnan_m

   implicit none
   private
   public setnan
   interface setnan
      module procedure setnan4, setnan8
   end interface setnan
   integer, parameter :: r4 = selected_real_kind(6)
   integer, parameter :: r8 = selected_real_kind(10)
   integer, parameter :: i4 = selected_int_kind(6)
   integer, parameter :: i8  = selected_int_kind(12)
   integer(kind=i4) inan4
   real(kind=r4) nan4
   equivalence (nan4,inan4)
   ! Signalling NaN value from Kahan's IEEE 754 description.
   data inan4 / Z'7F800001' /

   integer(kind=i4) :: inan8(2)
   real(kind=r8) :: nan8
   equivalence (nan8,inan8(1))
   ! For little endian machine, use this reversed order.
   data inan8 / Z'00000001', Z'7FF00000' /
   
contains
   elemental subroutine setnan4(x)
      real(kind=r4), intent(out) :: x
      x = nan4
      return
   end subroutine setnan4
   
   elemental subroutine setnan8(x)
      real(kind=r8), intent(out) :: x
      x = nan8
   end subroutine setnan8

end module setnan_m

