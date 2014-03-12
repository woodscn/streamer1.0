module init 
! Initialize the flow field. At the moment, this gives an initial guess. If 
!   allowing grid to propogate from upstream, it would be better to 
contains

  subroutine channel_init( nx0 , ny , nx , h0 )
    use global_data
    use primary_variables
    use window_data
    use boundary_conditions_mod
    implicit none
    integer , intent(in) :: nx0 , ny , nx
    real    , intent(in) :: h0
    integer :: i , j
    real :: dx , dy , ymin2
  dx = 0.02 ; dxi = dx ; deta = dx
  nxi = nx ; neta = ny

allocate( &
         x(ny,nx) , y(ny,nx) , qold(ny,nx,2) , hold(ny,nx), &
         E(ny,nx,8) , h(ny,nx) )
    x = 0. ; y = 0. ; qold = 0. ;  hold = 0. ;
    E = 0. ; h = 0. ;! xi = 0. ; eta = 0.
	do i = 1 , nx
       x(:,i) = xmin + dx*(i-1)
       if( x(1,i) .lt. .5 )then
         ymin2 = ymin
       elseif( x(1,i) .lt. 1. )then
         ymin2 = ( x(1,i) - .5 )*tan(15./180.*3.14159265) + ymin
       else
         ymin2 = .5*tan(15./180.*3.14159265) + ymin
       end if
       y(1,i) = ymin2
       dy = (ymax-ymin2)/ny
       do j = 2 , ny
          y(j,i) = y(j-1,i) + dy
       end do
    end do
    h = h0
    E(:,1,:) = upstream_boundary_cons( size(E,1) )
    do i = 2 , nx
       E(:,i,:) = E(:,i-1,:)
    end do


  end subroutine channel_init


  subroutine basic_init( nx0 , ny , nx , h0 ) ! Initialize field with uniform flow
! Inputs:  global_data: E, h, xi, eta, x, y, t, dt, neta, nxi, deta, dxi
!          window_data: xmin, xmax, ymin, ymax
! Outputs: fully initialized E, h, xi, eta, x, y for a given boundary condition.
!          At the moment, this simply applies the upstream boundary condition to 
!          the entire region. It also defines an overall h of 0.25, and a uniform 
!          xi,eta grid that is initially identical to the x,y grid.
!    use vars
    use global_data
    use primary_variables
    use window_data
    use boundary_conditions_mod
    implicit none
    integer , intent(in) :: nx0 , ny , nx
    real    , intent(in) :: h0
    integer :: i , j 
    real , allocatable :: A(:,:) , B(:,:) , L(:,:) , M(:,:)
    real , allocatable :: W(:,:,:)
    real :: dx , dy
    
!    write(*,*) "Allocating..."
!    allocate( &
!         W(ny,nx,4) , x(ny,nx) , y(ny,nx) , qold(ny,nx,2) , hold(ny,nx), &
!         E(ny,nx,8) , h(ny,nx) )
!    W = 0. ; x = 0. ; y = 0. ; qold = 0. ;  hold = 0. ;
!    E = 0. ; h = 0. ;! xi = 0. ; eta = 0.
allocate( &
         x(ny,nx) , y(ny,nx) , qold(ny,nx,2) , hold(ny,nx), &
         E(ny,nx,8) , h(ny,nx) )
    x = 0. ; y = 0. ; qold = 0. ;  hold = 0. ;
    E = 0. ; h = 0. ;! xi = 0. ; eta = 0. 
! Define an initial, uniform mesh in x,y. These arrays give an x,y coordinate
!   for each point in computational space. 
! Determine initial grid spacing given a desired number of computational points
    dx = (xmax-xmin)/(nx0-1) ; dy = (ymax-ymin)/(ny-1) 
!    dx = .01 ; dy = .01
! Use identical spacing, etc. for computational grid.
    dxi = dx ; deta = dy ; neta = ny ; nxi = nx
!! Use a cell-centered mesh
    x(:,1) = xmin! + 0.5*dx
    y(1,:) = ymin! + 0.5*dy
! Initialize x, y.
    do i = 2 , nx
       x(:,i) = x(:,i-1) + dx
    end do
    do i = 2 , ny
       y(i,:) = y(i-1,:) + dy
    end do
! Initialize h.
    h = h0
! Apply upstream boundary condition to entire region.
    E(:,1,:) = upstream_boundary_cons( size(E,1) )
    do i = 2 , nx
       E(:,i,:) = E(:,i-1,:)
    end do
!    deallocate( W )

  end subroutine basic_init

end module init

