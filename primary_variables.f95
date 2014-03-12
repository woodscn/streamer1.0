module primary_variables
  implicit none
  save
  real , allocatable :: E(:,:,:) ! Array of conserved variables. Primary data structure.
  real , allocatable :: x(:,:) , y(:,:) ! Arrays of physical coordinates x , y
  real , allocatable :: h(:,:) ! Grid-angle preserving function h
  real , allocatable :: qold(:,:,:) , hold(:,:)
end module primary_variables

