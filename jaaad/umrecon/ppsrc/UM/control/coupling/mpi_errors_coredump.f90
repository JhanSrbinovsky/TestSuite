subroutine mpi_errors_coredump (comm, error)

   implicit none

   integer(kind=4) :: comm, error

   external abort

   write (6,*) 'mpi_errors_coredump!'

   call abort

end subroutine mpi_errors_coredump
