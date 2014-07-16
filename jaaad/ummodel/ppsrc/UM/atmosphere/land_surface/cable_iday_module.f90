module cable_iday_mod
   implicit none
   public
   integer, pointer, save :: iday_number   
   real, pointer:: Lai_Ma(:,:)

contains

subroutine iday_kick(iday, cable_lai)
   integer,target :: iday
   real, target :: cable_lai(:,:) 
   iday_number => iday  
   Lai_Ma      =>cable_lai
end subroutine iday_kick

end module cable_iday_mod
