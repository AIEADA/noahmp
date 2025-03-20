program test

   use pnetcdf
   use mpi
   implicit none

   integer :: ncid, varid, ierr
   integer :: comm, info
   integer, parameter :: NX = 100, NY = 100
   integer, dimension(:, :), allocatable :: data
   integer :: my_rank, nprocs
   integer :: start(2), count(2)

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

   allocate (data(NX, NY/nprocs))

   data = my_rank+1

   ! Create a file using parallel I/O
   ierr = nf90mpi_create(MPI_COMM_WORLD, "output.nc", NF90_CLOBBER, MPI_INFO_NULL, ncid)

   ! Define dimensions and variable
   ierr = nf90mpi_def_dim(ncid, "x", INT(NX,8), varid)
   ierr = nf90mpi_def_dim(ncid, "y", INT(NY,8), varid)
   ierr = nf90mpi_def_var(ncid, "data", NF90_INT, (/1, 2/), varid)

   ! End definition mode
   ierr = nf90mpi_enddef(ncid)

   ! Select portion to write
   start = (/1, my_rank*(NY/nprocs)+1/)
   count = (/NX, NY/nprocs/)

   ! Write data
   ierr =  nf90mpi_bput_var(ncid, varid, data, NF90_REQ_NULL, start, count, (/1,1/), (/6,1/))

   ! Close file
   ierr =  nf90mpi_close(ncid)

  deallocate (data)
  call MPI_FINALIZE(ierr)

end program test
