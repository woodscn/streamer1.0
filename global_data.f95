module global_data
  implicit none
  save
! Data within this module is shared with many parts of the program. 
  real    :: dxi , deta ! Constant cell width in computational space
  real    :: t = 0. , dt ! Time and physical step size
  integer :: neta , nxi ! Number of computational grid points
end module global_data

