module window_data
  implicit none
  save
! This module contains data on the computational "window" described by Hui. 
!   That is, the min and max parameters defined here describe a box within 
!   which the flow is computed. In reality, a UCS computation will always 
!   extend beyond this box, to some extent.
  real :: xmin , xmax , ymin , ymax  
end module window_data

