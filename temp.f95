module helper_functions
contains
  function area_func( deta , dxi , n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    real                  :: area_func 
    real    , intent(in)  :: deta , dxi
    integer , intent(in)  :: n
! Since we are first advancing in the xi direction, the "area" will be deta. 
!   During the second step of Strang, the opposite is true.
    if (n .ne. 2) then
       area_func = deta
    else
       area_func = dxi
    end if
  end function area_func

  function dlambda_func( dt , n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    real                  :: dlambda_func
    real    , intent(in)  :: dt
    integer , intent(in)  :: n
! In Strang splitting, the solution is  advanced half a time step in one
!   direction, then a full step in the other, and the remaining half in
!   the first direction.
    if (n .eq. 1 .or. n .eq. 3) then
       dlambda_func = dt*0.5
    else
       dlambda_func = dt
    end if
  end function dlambda_func

  function jshift_func( n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    integer :: jshift_func
    integer , intent(in) :: n
    if (n .ne. 2) then
       jshift_func = 0
    else
       jshift_func = 1
    end if
  end function jshift_func

  function ishift_func( n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    integer :: ishift_func
    integer , intent(in) :: n
    if (n .ne. 2) then
       ishift_func = 1
    else
       ishift_func = 0
    end if
  end function ishift_func
  
end module helper_functions

module physical_data 
! This module contains only physical parameters such as the specific gas constant.
  implicit none
  save
  real , parameter :: gamma = 1.4 
end module physical_data

module window_data
  implicit none
  save
! This module contains data on the computational "window" described by Hui. 
!   That is, the min and max parameters defined here describe a box within 
!   which the flow is computed. In reality, a UCS computation will always 
!   extend beyond this box, to some extent.
  real :: xmin , xmax , ymin , ymax  
end module window_data

module global_data
  implicit none
  save
! Data within this module is shared with many parts of the program. 
  real    :: dxi , deta ! Constant cell width in computational space
  real    :: t = 0. , dt ! Time and physical step size
  integer :: neta , nxi ! Number of computational grid points
end module global_data

module boundary_data
  implicit none
  save
  integer :: boundary_switch = 0
end module boundary_data

module primary_variables
  implicit none
  save
  real , allocatable :: E(:,:,:) ! Array of conserved variables. Primary data structure.
  real , allocatable :: x(:,:) , y(:,:) ! Arrays of physical coordinates x , y
  real , allocatable :: h(:,:) ! Grid-angle preserving function h
  real , allocatable :: qold(:,:,:) , hold(:,:)
end module primary_variables

module projection_routines
contains
  subroutine normal_vector(E1 , E2 , n , metric )
    implicit none
    real    , intent(in)  :: E1(:) , E2(:) 
    integer , intent(in)  :: n
    real    , intent(out) :: metric(:,:)
    real                  :: normal(2) , tangential(2) , El(8) = 0. , Er(8) = 0.
    real                  :: deltal , deltar , psi
    logical               :: tick
    
    El(5:8) = E1 ; Er(5:8) = E2
    deltal = ( El(5)*El(8) - El(6)*El(7) )
    deltar = ( Er(5)*Er(8) - Er(6)*Er(7) )
    ! It is better to use the metric components directly, without using
    !   the angle psi; it improves stability.
    tick = .false. ! Use psi-formulation?
    if( tick )then
       ! Compute the angle of rotation, using the average of the two states.
       if( n .ne. 2 )then
          psi = 0.5*( atan2( El(8) , El(7) ) + atan2( Er(8) , Er(7) ) )
       else
          psi = 0.5*( atan2( El(6) , El(5) ) + atan2( Er(6) , Er(5) ) )
       end if
       
       metric(1,:) = (/  sin(psi) , -cos(psi) /)
       metric(2,:) = (/  cos(psi) ,  sin(psi) /)
       
       if( n .eq. 2 ) metric = -1.*metric
    else
       if( n .ne. 2 )then
          normal     = (/  El(8) , -El(7) /)/deltal &
                     + (/  Er(8) , -Er(7) /)/deltar
          
          tangential = (/  El(7) ,  El(8) /)/deltal &
                     + (/  Er(7) ,  Er(8) /)/deltar
       else
          normal     = (/ -El(6) ,  El(5) /)/deltal &
                     + (/ -Er(6) ,  Er(5) /)/deltar
          
          tangential = (/ -El(5) , -El(6) /)/deltal &
                     + (/ -Er(5) , -Er(6) /)/deltar
       end if
       
      
!       if( n .ne. 2 )then
!          normal     = 0.5*( (/  El(8) , -El(7) /)/sqrt( El(8)**2 + El(7)**2 ) &
!                           + (/  Er(8) , -Er(7) /)/sqrt( Er(8)**2 + Er(7)**2 ) )
!          tangential = 0.5*( (/  El(7) ,  El(8) /)/sqrt( El(8)**2 + El(7)**2 ) &
!                           + (/  Er(7) ,  Er(8) /)/sqrt( Er(8)**2 + Er(7)**2 ) )
!       else
!          normal     = 0.5*( (/ -El(6) ,  El(5) /)/sqrt( El(5)**2 + El(6)**2 ) &
!                           + (/ -Er(6) ,  Er(5) /)/sqrt( Er(5)**2 + Er(6)**2 ) )
!          tangential = 0.5*( (/ -El(5) , -El(6) /)/sqrt( El(5)**2 + El(6)**2 ) &
!                           + (/ -Er(5) , -Er(6) /)/sqrt( Er(5)**2 + Er(6)**2 ) )
!       end if
       
       normal     = normal    /sqrt(     normal(1)**2 +     normal(2)**2 )
       tangential = tangential/sqrt( tangential(1)**2 + tangential(2)**2 )
       metric(1,:) = normal
       metric(2,:) = tangential
    end if
       
    
  end subroutine normal_vector
  
  subroutine velocity_projector( W , metric )
    implicit none
    real , intent(inout) :: W(:)
    real , intent(in)    :: metric(:,:)
    real                 :: temp(2)
    temp = (/ W(2) , W(3) /)
    
    temp = matmul(metric,temp)
    W(2) = temp(1) ; W(3) = temp(2)
    
  end subroutine velocity_projector
  
  subroutine velocity_deprojector( W , metric )
    implicit none
    real , intent(in )   :: metric(:,:)
    real , intent(inout) :: W(:)
    real                 :: temp(2) , metric_temp(2,2)
    temp = (/ W(2) , W(3) /)
    ! Invert the metric
    metric_temp(1,:) = (/ metric(2,2) , -metric(1,2) /)
    metric_temp(2,:) = (/-metric(2,1) ,  metric(1,1) /)
    !      metric_temp = metric_temp*1./( metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1) )
    ! Apply the inverse of the metric to the velocity vector
    temp = matmul(metric_temp,temp)
    W(2) = temp(1) ; W(3) = temp(2)
  end subroutine velocity_deprojector
  
end module projection_routines

module vars
! This module contains routines for transforming from conservative variables (e.g. 
!   rho*delta*u) to primitive variables (e.g. u). It uses only the constant gamma,
!   and the arrays of conservative and primitive variables to update either 
!   E(:,:,1:4) or W(:,:,:).
contains
  subroutine primtocons( E , W , geom ) 
! Inputs: E(neta,nxi,5:8) - the array of conservative variables. The last 4 "rows" in
!   the third index are the geometric variables A, B, L, M. The conserved flow 
!   variables will be overwritten based on the  primitive variables. W(neta,nxi,4)
!   - the array of primitive flow variables rho, u, v, p, or density, velocity 
!   components and pressure.

! Outputs: E(neta,nxi,1:4) - the conservative flow variables rho*delta, rho*delta*u,
!   rho*delta*v, and rho*delta*e, where delta = A*M-B*L, and e is the specific 
!   total energy.
    use physical_data
    implicit none
    real , intent(in)    :: W(:,:,:) , geom(:,:,:)
    real , intent(out) :: E(:,:,:)
    real , allocatable :: delta(:,:)
    integer :: neta , nxi ! Recomputed here for the sake of code independence
    neta = size(E,1) ; nxi = size(E,2)
    allocate( delta(neta,nxi) )
    delta = geom(:,:,1)*geom(:,:,4) - geom(:,:,2)*geom(:,:,3) ! AM-BL, see documentation.
! Compute conservative variables based on primitive variables
    E(:,:,1) = W(:,:,1)*delta
    E(:,:,2) = E(:,:,1)*W(:,:,2)
    E(:,:,3) = E(:,:,1)*W(:,:,3)
    E(:,:,4) = delta*(W(:,:,1)*0.5*(W(:,:,2)**2 + W(:,:,3)**2) &
         + W(:,:,4)*(1./(gamma-1.)))
    deallocate( delta )
  end subroutine primtocons

  subroutine constoprim( E , W , geom )
! Inputs: E(neta,nxi,8) - the array of conservative variables. The first 4 "rows"
!   in the third index are conservative flow variables, while the last 4 "rows" 
!   are geometric variables.

! Outputs: W(neta,nxi,4) - the array of primitive flow variables rho, u, v, p.
    use physical_data
    implicit none   
    real , intent(in)  :: E(:,:,:) , geom(:,:,:)
    real , intent(out) :: W(:,:,:)
    real , allocatable :: delta(:,:)
    integer :: neta , nxi
    neta = size(E,1) ; nxi = size(E,2)
    allocate( delta(neta,nxi) )
    delta = geom(:,:,1)*geom(:,:,4) - geom(:,:,2)*geom(:,:,3) ! AM-BL, see documentation.
! Compute primitive variables from conservative variables
    W(:,:,1) = E(:,:,1)/delta
    W(:,:,2) = E(:,:,2)/E(:,:,1)
    W(:,:,3) = E(:,:,3)/E(:,:,1)
    W(:,:,4) = (gamma-1.)*( E(:,:,4)/delta - 0.5 * W(:,:,1)*( W(:,:,2)**2 + W(:,:,3)**2 ) )
    deallocate(delta)
  end subroutine constoprim

  subroutine constoprimvec( E , W , geom )
! Inputs: E(1,1,8) - the array of conservative variables. The first 4 "rows"
!   in the third index are conservative flow variables, while the last 4 "rows" 
!   are geometric variables.

! Outputs: W(1,1,4) - the array of primitive flow variables rho, u, v, p.
    use physical_data
    implicit none   
    real , intent(in)  :: E(:) , geom(:)
    real , intent(out) :: W(:)
    real               :: delta
    delta = geom(1)*geom(4) - geom(2)*geom(3) ! AM-BL
! Compute primitive variables from conservative variables
    W(1) = E(1)/delta
    W(2) = E(2)/E(1)
    W(3) = E(3)/E(1)
    W(4) = (gamma-1.)*(E(4)/delta-0.5*W(1)*(W(2)**2+W(3)**2))
  end subroutine constoprimvec

end module vars

module boundary_conditions_mod
! Apply boundary conditions to the array of conservative variables and to h.
contains
  subroutine boundary_conditions( E , h , flag , h0 )
! Inputs: E - array of conservative variables, h - shape-preserving function,
!   flag - optional flag, indicates whether or not to apply downstream 
!   conditions. It's presence causes the routine to apply upstream conditions
!   only.
    use vars
    implicit none
    real , intent(inout) :: E(:,:,:) , h(:,:)
    integer              :: ny , nx , j , i 
    real , intent(in) , optional :: h0
    logical , intent(in) :: flag 
    real , allocatable   :: Wtemp(:,:,:) , Etemp(:,:,:)
    real                 :: Wl(4) , Wr(4)

    ny = size(E,1) ; nx = size(E,2)
    allocate( Wtemp(ny,1,4) , Etemp(ny,1,8) )
! Define upstream conditions
    call boundary_conditions_prim ( Wtemp)
! Apply zero-derivative condition to h at upstream boundary
    if( size(h,2) .gt. 1 )then
       h(:,1) = h(:,2)
    else
       h(:,1) = h0
    end if

! Define metric components, using a cartesian grid at the inflow boundary.
    do j = 1 , ny
       Etemp(j,1,:) = (/ 0. , 0. , 0. , 0. , 1. , 0. , 0. , 1. /) 
       if( flag ) Etemp(j,1,5) = (h(j,1)*Wtemp(j,1,2))/(h(ny,1)*Wtemp(ny,1,2))
    end do

! Convert Wtemp (primitive variables) to Etemp (conservative variables).
    call primtocons( Etemp(:,:,1:4) , Wtemp , Etemp(:,:,5:8) ) 

! Apply Dirichlet boundary conditions at upsteram, requiring E = Etemp
    E(:,1,:) = Etemp(:,1,:) ! Apply upstream BC 
    deallocate( Wtemp , Etemp )

  end subroutine boundary_conditions
  
  subroutine boundary_conditions_prim( Wtemp )
    use boundary_data
    implicit none
    real , intent(out) :: Wtemp(:,:,:)
    integer              :: ny , nx , i , j , switch
    real                 :: Wl(4) , Wr(4)

    ny = size(Wtemp,1) ; nx = size(Wtemp,2)
    if( boundary_switch .eq. 1 )then
! Basic, uniform flow
       do j = 1 , ny
          Wtemp(j,1,:) = (/ 1. , 1.8*sqrt(1.4*1.00/1.0) , 0. , 1. /) 
       end do
    elseif( boundary_switch .eq. 2 )then
       do j = 1 , ny
! Linear velocity inflow
          Wtemp(j,1,:) = (/ 1. , 1.+real(j)/real(ny) , 0. , 1. /)
       end do
    elseif( boundary_switch .eq. 3 )then
! The Riemann problem of Hui
       Wl = (/ 1.0 , 2.4*sqrt(1.4*1.00/1.0) , 0. , 1.00 /)
       Wr = (/ 0.5 , 7.0*sqrt(1.4*0.25/0.5) , 0. , 0.25 /)
       do j = 1 , ceiling(real(ny)/2)
          Wtemp(j,1,:) = Wl
       end do
       do j = ceiling(real(ny)/2)+1 , ny
          Wtemp(j,1,:) = Wr
       end do
    end if

  end subroutine boundary_conditions_prim

end module boundary_conditions_mod

module init 
! Initialize the flow field. At the moment, this gives an initial guess. If 
!   allowing grid to propogate from upstream, it would be better to 
contains
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
    allocate( &
         W(ny,nx,4) , x(ny,nx) , y(ny,nx) , qold(ny,nx,2) , hold(ny,nx), &
         E(ny,nx,8) , h(ny,nx) )
    W = 0. ; x = 0. ; y = 0. ; qold = 0. ;  hold = 0. ; 
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
    call boundary_conditions( E , h , .false. , h0 )
    do i = 2 , nx
       E(:,i,:) = E(:,i-1,:)
    end do
    deallocate( W )

  end subroutine basic_init

end module init

module write_files_mod
! This module is responsible for outputting data into a form that Matlab can 
!   read for data visualization. 
contains
  subroutine write_files(x,y,E,h,t,dt)
! Inputs: x, y, E (array of conserved variables), t, dt
! Outputs: two_D.dat, an ASCII file that contains the primitive variables together
!          with their coordinate values and the associated time. Data is designed
!          to be read by the matlab file geom_data_reader.m. It should be noted 
!          that the file is opened and headers are written in the main program.
    use vars
    implicit none
    real , intent(in)  :: x(:,:) , y(:,:) , E(:,:,:) , t , dt , h(:,:)
    integer            :: i , j , neta , nxi
    real , allocatable :: W(:,:,:) , u(:,:) , v(:,:) 

    neta = size(E,1) ; nxi = size(E,2)
    allocate( W(neta,nxi,4) , u(neta,nxi) , v(neta,nxi) )
! We output primitive variables
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
    u = E(:,:,2)/E(:,:,1) ; v = E(:,:,3)/E(:,:,1)
    write(1,*) ' ' 
    write(1,*) 't=',t,'nx=',nxi,'ny=',neta,'dt=',dt
    write(1,*) ' '
    write(1,*) 'x= ','y= ','u= ','v= ','rho= ','P= '
    do i = 1 , nxi
       do j = 1 , neta
          write(1,*) ' ', x(j,i) , y(j,i) , u(j,i) , v(j,i) , W(j,i,1) , W(j,i,4) &
               , E(j,i,5) , E(j,i,6) , E(j,i,7) , E(j,i,8) , h(j,i)
       end do
    end do
    deallocate( W , u , v )
  end subroutine write_files
end module write_files_mod
 
module geom_update_mod
! This module exists to advance the x,y coordinates in time, to create/destroy points in order to keep the computation within the computational window, and to solve for updated values of geometric variables.
contains
  subroutine geom_update()
! geom_update acts only as the master routine to direct the action of the three individual subroutines.
    call coords_update ! Move the x,y coordinates with the fluid velocity
    call windower ! Create/destroy points to maintain the computational window
  end subroutine geom_update

  subroutine coords_update()
! This routine advances the values of x,y in time according to the formula dx = h*u*dt and dy = h*v*dt. It also updates t.
    use global_data
    use window_data
    use primary_variables
    implicit none

! Advance spatial coordinates
       x = x + ( hold*qold(:,:,1) + h*E(:,:,2)/E(:,:,1) )*dt*.5
       y = y + ( hold*qold(:,:,2) + h*E(:,:,3)/E(:,:,1) )*dt*.5

  end subroutine coords_update

  subroutine windower()
! Windower is tasked with keeping the solution within a designated computational window. It evaluates whether the grid has overrun the maximum and minimum required values of x, and creates/destroys columns as necessary. 
    use global_data
    use primary_variables
    use window_data
    use boundary_conditions_mod

    implicit none
    logical :: overrun = .false. , shortfall = .false.
    real , allocatable :: x_t(:,:) , y_t(:,:) , E_t(:,:,:) , h_t(:,:) , xi_t(:,:) , eta_t(:,:)
    integer :: j
    logical :: switch1 = .false. , switch2 = .false.
! In contrast to most other subroutines, this one has access to global values that include neta and nxi, so they are not recomputed here.

! Test whether the grid extends too far beyond the right of computational window
    if( sum(x(:,nxi))/real(neta) .gt. xmax )then
       overrun = .true. 
    else
       overrun = .false.
    end if

! Test whether the grid does not extend far enough to the left of the computational window
    if (maxval(x(:,1)) .ge. xmin+dxi) then
       shortfall = .true.
    else
       shortfall = .false.
    end if

! Begin to update coordinates based on the results of the windowing tests above
    if (overrun .and. (.not. shortfall))then
! The grid extends too far in the right (downstream), but not in the left.

! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h 
! Remove a column and reallocate variable arrays
       call deallocate_real
       nxi = nxi - 1
       call allocate_real
! Copy the temporary arrays back into the variable arrays, minus the last, "overrun" column.
       x = x_t(:,1:nxi) ; y  = y_t (:,1:nxi) ; E   = E_t  (:,1:nxi,:) 
       h = h_t(:,1:nxi) 
       call deallocate_temp
    elseif (shortfall .and. (.not. overrun)) then
! The grid does not extend far enough to the left (upstream), but the right is fine.

! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h
! Add a column and reallocate variable arrays
       call deallocate_real
       nxi = nxi + 1
       call allocate_real
! Copy all but the first column from the temporary arrays.
       x(:,2:nxi) = x_t ; y (:,2:nxi) = y_t  ; E  (:,2:nxi,:) = E_t 
       h(:,2:nxi) = h_t 
       call deallocate_temp
! Call boundary_conditions to compute the first column for E,h.
       call upwind_boundary()
    elseif (shortfall .and. overrun) then
! The grid has the appropriate length; it just needs to be shifted upstream.
!       write(*,*) "Both overrun and shortfall."
! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h
! Overwrite all but the first column with data from the temporary arrays
       x(:,2:nxi)   = x_t(:,1:(nxi-1))   ; y(:,2:nxi)   = y_t  (:,1:(nxi-1))
       E(:,2:nxi,:) = E_t(:,1:(nxi-1),:) ; h(:,2:nxi)   = h_t  (:,1:(nxi-1))
       call deallocate_temp
       call upwind_boundary()
    endif
!    write(*,*) "Windower Done"   
    contains
      subroutine upwind_boundary()
! Call boundary_conditions to compute the first column for E,h.
        call boundary_conditions( E , h , .true. )
! Update the first column of the remaining arrays.
        x(:,1) = xmin
        do j = 1 , neta
           y(j,1) = ymin + (real(j)-1.)*deta
        end do
      end subroutine upwind_boundary

      subroutine allocate_temp ! Allocate temporary arrays
        allocate( x_t(neta,nxi) , y_t(neta,nxi) , E_t(neta,nxi,8) &
             , h_t(neta,nxi) )
        x_t = 0. ; y_t = 0. ; E_t = 0. ; h_t = 0.
      end subroutine allocate_temp

      subroutine allocate_real ! Allocate variable arrays
        allocate( x(neta,nxi) , y(neta,nxi) , E(neta,nxi,8)   &
             , h(neta,nxi) , qold(neta,nxi,2) , hold(neta,nxi) )
        x = 0. ; y = 0. ; E = 0. ; h = 0.
        qold = 0. ; hold = 0.
      end subroutine allocate_real

      subroutine deallocate_temp ! Deallocate temporary arrays
        deallocate( x_t , y_t , E_t , h_t )
      end subroutine deallocate_temp

      subroutine deallocate_real ! Deallocate variable arrays
        deallocate( x , y , E , h , qold , hold )
      end subroutine deallocate_real

  end subroutine windower
  
end module geom_update_mod

module h_update_mod
! This module is responsible for solving the PDE for h. It uses an iterative, 
!   first-order, backward difference scheme. Boundary values must be provided 
!   at h(1,:) and h(:,1).
contains

  subroutine h_update( A , B , L , M , h , u , v , xi , eta , h0 )
    use boundary_conditions_mod
    implicit none
    real , intent(in)   :: A(:,:) , B(:,:) ,  L(:,:) ,   M(:,:)
    real , intent(in)   :: u(:,:) , v(:,:) , xi(:,:) , eta(:,:) , h0
    real , intent(inout):: h(:,:)
    real                :: tol  , err , dxi , deta
    real , allocatable  :: dthetadxi(:,:) , dthetadeta(:,:) , temp(:,:)
    real , allocatable  :: alpha(:,:) , beta(:,:) , gamma(:,:) , g(:,:)
    real , allocatable  :: S2(:,:) , T2(:,:) , q(:,:) , theta(:,:) 
    integer             :: ny , nx , i , j , n , k 

    ny = size(h,1) ; nx = size(h,2) ! Define these size parameters
    dxi = xi(2,2)-xi(1,1) ; deta = eta(2,1)-eta(1,1) 
    allocate( g(ny,nx) , q(ny,nx) , theta(ny,nx) , S2(ny,nx) , T2(ny,nx) )
! The PDE we solve to determine h is written in terms of polar velocity components. 
    call vpolar( u , v , q , theta ) ! Transform from cartesian to polar velocity.
    h(ny,:) = 0.25
! The PDE has better continuity properties when written in terms of g rather than h.
    g = log(q*h)
    S2 = L**2 + M**2 ; T2 = A**2 + B**2 ! Define useful parameters.
    tol = 1.0e-14 ! Define convergence criterion.
! We need to compute the gradient of theta (polar velocity angle). We use a 
!   backward-difference method, for consistency with the iterative scheme.
    allocate( dthetadxi(ny,nx) , dthetadeta(ny,nx) )
    do i = 2 , nx
       do j = 2 , ny
          dthetadxi(j,i) = (theta(j,i)-theta(j,i-1))/dxi ;
          dthetadeta(j,i)= (theta(j,i)-theta(j-1,i))/deta;
       end do
    end do
! The PDE can be written in the form alpha*dg/dxi + beta*dg/deta + gamma = 0.
!   We define alpha, beta, and gamma.
    allocate( alpha(ny,nx) , beta(ny,nx) , gamma(ny,nx) )
    alpha = S2*(A*cos(theta)-B*sin(theta))
    beta  = T2*(M*cos(theta)-L*sin(theta))
    gamma =-S2*(B*sin(theta)+M*sin(theta))*dthetadxi &
          + T2*(L*cos(theta)+M*sin(theta))*dthetadeta

! Iteratively solve for g using a backward difference based on BC's at g(1,:)
!   and g(:,1). 
    err = 1 ; 
    allocate( temp(ny,nx) )
    temp = g
    do while (err .gt. tol)
       temp(:,1) = temp(:,2) + deta/beta(:,2) * ( gamma(:,2) + alpha(:,2) * 0 )
       do i = 2 , nx
          do j = 2 , ny
             temp(j,i) = (alpha(j,i)*g(j,i-1)/dxi + beta(j,i)*g(j-1,i) &
                   /deta - gamma(j,i) )/(alpha(j,i)/dxi + beta(j,i)/deta)
          end do
       end do
       temp(:,1) = temp(:,2) ; temp(1,:) = temp(2,:)
       err = maxval(abs(temp-g))
       g   = temp
    end do
! Convert the updated g back to h.
    h = exp(g)/q
    h = h/maxval(h)*h0

    deallocate( temp )
    deallocate( g , q , theta )
    deallocate( dthetadxi , dthetadeta )
    deallocate( alpha , beta , gamma )

  end subroutine h_update

  subroutine vpolar(u,v,r,theta)
! This subroutine computes polar velocity components based on cartesian components.
    real, intent(in) :: u(:,:) , v(:,:)
    real, intent(out):: r(:,:) , theta(:,:)

    r     = sqrt(u**2+v**2)
    theta = atan2(v,u)
  end subroutine vpolar

end module h_update_mod

module muscl_mod
contains
  subroutine MUSCL( Wplus1 , Wplus2 , W , Wminus1 , Wl , Wr )
! This subroutine provides the second-order MUSCL interpolation used by Hui.
!   See his 2-D paper for details.
    real , intent(in)  :: Wplus1(:) , Wplus2(:) , W(:) , Wminus1(:)
    real , intent(out) :: Wl(:) , Wr(:)
    real , allocatable :: rplus(:) , rminus(:) 
    integer :: n , i
    n = size(W,1)
    allocate( rplus(n) , rminus(n) )

    rplus  = (Wplus1-W)/(Wplus2- Wplus1)
    rminus = (Wplus1-W)/(W   -  Wminus1)

    do i = 1 , n
       Wr(i) = Wplus1(i) - 0.5*( Wplus2(i) - Wplus1(i)  )*phi(rplus(i) )
       Wl(i) =    W(i)   + 0.5*(    W(i)   - Wminus1(i) )*phi(rminus(i))
 end do
    deallocate( rplus , rminus )

  contains
    function phi( r )
      real  :: phi
      phi = max(0.,min(1.,r))
    end function phi

  end subroutine MUSCL
end module muscl_mod

module states2
contains
  subroutine left_and_right_states(i,j,n,x,E,h,Er,hr,Wr,El,hl,Wl,h0,ramp)
    use boundary_conditions_mod
    use vars
    use helper_functions
    use projection_routines
    implicit none
    logical , intent(in)  :: ramp
    integer , intent(in)  :: i , j , n
    real    , intent(in)  :: E(:,:,:) , h(:,:) , h0 , x(:,:)
    real    , intent(out) ::  Er(4) , hr , Wr(4) , El(4) , hl , Wl(4) 
    real    :: pi , theta , metric(2,2)
    integer :: is , js , neta , nxi , jshift , ishift
    real , allocatable :: W(:,:,:) , Wtemp(:,:,:) , Etemp(:,:,:) , htemp(:,:) 

    jshift = jshift_func(n)
    ishift = ishift_func(n)
    js = j - jshift
    is = i - ishift

    neta = size(E,1) ; nxi = size(E,2)
    allocate( W(neta,nxi,4) , Wtemp(neta,1,4) , Etemp(neta,1,8) , htemp(neta,1) )
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
    if( n  .ne. 2 )then
       call boundary_conditions_prim( Wtemp )
       call boundary_conditions( Etemp , htemp , .true. , h0 )
    end if

    if( is .eq. 0 )then
       Er = E(j ,i ,5:8)
       hr = h(j ,i     )
       Wr = W(j ,i , : )
       El = Etemp(j,i,5:8)
       hl = hr
       Wl = Wtemp(j,i, : )
    elseif( i .eq. nxi + 1 )then
       El = E(js,is,5:8)
       hl = h(js,is    )
       Wl = W(js,is, : )
       Er = El 
       hr = hl
       Wr = Wl
    elseif( js .eq. 0 )then
       Er = E(j,i,5:8)
       hr = h(j,i    )
       Wr = W(j,i, : )
       El = Er
       hl = hr
       Wl = Wr
       if( x(j,i) .ge. 0.5 .and. x(j,i) .lt. 1. .and. ramp )then
          pi = 3.141592653589793
          theta = pi/12.
          metric(1,:) = (/  cos(theta) , sin(theta) /)
          metric(2,:) = (/ -sin(theta) , cos(theta) /)
          call velocity_projector( Wl , metric )
          Wl(3) = -Wl(3)
          call velocity_deprojector( Wl , metric )
          !      Flat wall segment
       else
          Wl(3) = -Wl(3)
       end if
    elseif( j .eq. neta + 1 )then
       El = E(js,is,5:8)
       hl = h(js,is    )
       Wl = W(js,is, : )
       Er = El
       hr = hl
       Wr = Wl
       Wr(3) = -Wr(3)
    else
       El = E(js,is,5:8)
       hl = h(js,is    )
       Wl = W(js,is, : )
       Er = E(j ,i ,5:8)
       hr = h(j ,i     )
       Wr = W(j ,i , : )
    end if
    deallocate( W , Wtemp , Etemp , htemp )
  end subroutine left_and_right_states
end module states2

module riemann
! This module contains the riemann solver used in my code. It is adapted from 
! the one presented by Toro in his book on Riemann problems in fluid mechanics.
! It is an exact, iterative solver, and it may eventually be beneficial to 
! use a faster, apporoximate solver instead. 
contains
! Inputs: Domlen - the "length" of the domain. Corresponds to either dxi or deta.
!         WinL   - the left state in terms of primitive variables
!         WinR   - the right state in terms of primitive variables
!         Wout   - the output state, used to compute fluxes, in primitive variables
!         Smax   - the speed of the fastest nonlinear wave in the problem. Useful 
!                     for computing time steps as part of a CFL condition.
  subroutine riemann_solve(DOMLEN, WINL, WINR, WOUT, HAVG)

    IMPLICIT NONE

!     Declaration of variables:

    INTEGER :: I, CELLS = 40
    REAL    :: DIAPH = 0., DOMLEN, DS = 0., DX, PM, MPA = 1., PS = 0., S = 0.,&
         &        TIMEOUT = 0.035, UM, US = 0., XPOS, SMAX, WOUT(4), WINL(4), WINR(4)
    REAL    :: GAMMA, G1, G2, G3, G4, G5, G6, G7, G8, &
         &           DL, UL, PL, CL, DR, UR, PR, CR, SL, SR , VL , VR , VS , &
         &           HAVG
    COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
    COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR 

    GAMMA = 1.4
    DIAPH = DOMLEN*0.5
    DL = WINL(1) ; UL = WINL(2) ; VL = WINL(3) ; PL = WINL(4)
    DR = WINR(1) ; UR = WINR(2) ; VR = WINR(3) ; PR = WINR(4)

!     Compute gamma related constants

    G1 = (GAMMA - 1.0)/(2.0*GAMMA)
    G2 = (GAMMA + 1.0)/(2.0*GAMMA)
    G3 = 2.0*GAMMA/(GAMMA - 1.0)
    G4 = 2.0/(GAMMA - 1.0)
    G5 = 2.0/(GAMMA + 1.0)
    G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
    G7 = (GAMMA - 1.0)/2.0
    G8 = GAMMA - 1.0

!   Compute sound speeds
    CL = SQRT(GAMMA*PL/DL)
    CR = SQRT(GAMMA*PR/DR)
    
!   The pressure positivity condition is tested for
    IF(G4*(CL+CR).LE.(UR-UL))THEN
!   The initial data is such that a vacuum is generated
       WRITE(6,*)
       WRITE(6,*)'***Vacuum is generated by data***' 
       WRITE(6,*)'***Program stopped***'
       WRITE(6,*)
       STOP
    ENDIF

!   Exact solution for pressure and velocity in star region is found

    CALL STARPU(PM, UM, MPA)
    
    DX = DOMLEN/REAL(CELLS)

!      Solution at point (X,T) = ( XPOS - DIAPH,TIMEOUT ) is found
         
    S = 0.
    CALL SAMPLE(PM, UM, S, DS, US, PS, VS, HAVG)

    WOUT = (/ DS, US, VS, PS /)
	 
    if(pm.le.pl)then
       SL = ABS( UL - CL )
    else 
       SL = ABS( UL - CL/(2.*GAMMA)*((GAMMA-1.)*PM/PR&
            &               +(GAMMA-1.))**.5 )
    end if

    if(pm.gt.pr)then
       SR = ABS( UR + CR )
    else 
       SR = ABS( UR + CR/(2.*GAMMA)*((GAMMA-1.)*PM/PR&
            &               +(GAMMA-1.))**.5 )
    end if
    smax = max(smax, max(SL,SR))
    !write(*,*) 'smax =  ',smax

  END SUBROUTINE riemann_solve


  SUBROUTINE STARPU(P, U, MPA)
    IMPLICIT NONE

!   Purpose: to compute the solution for pressure and velocity 
!             in the Star Region

!   Declaration of variables

    INTEGER :: I, NRITER
    
    REAL    :: DL, UL, PL, CL, DR, UR, PR, CR,&
         &        CHANGE, FL, FLD, FR, FRD, P, POLD, PSTART,&
         &        TOLPRE, U, UDIFF, MPA, VL, VR

    COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR
    DATA TOLPRE, NRITER/1.0E-14, 10/

!   Guessed value PSTART is computed

    CALL GUESSP(PSTART)

    POLD  = PSTART
    UDIFF = UR - UL

!    WRITE(6,*)'-----------------------------------------'
!    WRITE(6,*)'   Iteration number        Change   '
!    WRITE(6,*)'-----------------------------------------'

    DO 10 I = 1, NRITER
       CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
       CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
       P      = POLD - (FL + FR + UDIFF)/(FLD + FRD)
       CHANGE = 2.0*ABS((P - POLD)/(P + POLD))
       IF(CHANGE.LE.TOLPRE)GOTO 20
       IF(P.LT.0.0)P = TOLPRE
       POLD = P
10     CONTINUE
! I find that this error usually indicates that the solution has gone
! unstable at the flow level. The N-R iteration routine itself works 
! rather well.
       WRITE(6,*)'Divergence in Newton-Raphson iteration'
       STOP
20     CONTINUE

!   Compute velocity in Star Region
       U = 0.5*(UL + UR + FR - FL)

!            WRITE(6,*)'-----------------------------------------'
!            WRITE(6,*)'   Pressure             Velocity'
!            WRITE(6,*)'-----------------------------------------'
!            WRITE(6, 40)P/MPA, U
!            WRITE(6,*)'-----------------------------------------'

30     FORMAT(5X, I5,15X, F12.7)
40     FORMAT(2(F14.6, 5X))

       RETURN
     END SUBROUTINE STARPU


     SUBROUTINE GUESSP(PM)

!   Purpose: to provide a guess value for pressure 
!            PM in the Star Region. The choice is made
!            according to adaptive Riemann solver using
!            the PVRS, TRRS and TSRS approximate
!            Riemann solvers. See Sect. 9.5 of Capht. 9
!            of Ref. 1

       IMPLICIT NONE

!   Declaration of variables

       REAL    :: DL, UL, PL, CL, DR, UR, PR, CR,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,&
            &        CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,&
            &        PTL, PTR, QMAX, QUSER, UM, VL , VR

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR

       QUSER = 2.0

!   Compute guess pressure from PVRS Riemann solver

       CUP  = 0.25*(DL + DR)*(CL + CR)
       PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
       PPV  = AMAX1(0.0, PPV)
       PMIN = AMIN1(PL,  PR)
       PMAX = AMAX1(PL,  PR)
       QMAX = PMAX/PMIN

       IF(QMAX.LE.QUSER.AND.&
            & (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN

!      Selcect PVRS Riemann solver
          PM = PPV
       ELSE
          IF(PPV.LT.PMIN)THEN

!         Select Two-Rarefaction Riemann solver
             PQ  = (PL/PR)**G1
             UM  = (PQ*UL/CL + UR/CR + &
                  &            G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
             PTL = 1.0 + G7*(UL - UM)/CL
             PTR = 1.0 + G7*(UM - UR)/CR
             PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
          ELSE

!      Select Two-Shock Riemann solver with PVRS as estimate

             GEL = SQRT((G5/DL)/(G6*PL + PPV))
             GER = SQRT((G5/DR)/(G6*PR + PPV))
             PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL+GER)
          ENDIF
       ENDIF
       RETURN
     END SUBROUTINE GUESSP


     SUBROUTINE PREFUN(F, FD, P, DK, PK, CK)
!   Purpose: to evaluate the pressure fucntions
!            FL and FR in exact Riemann solver

       IMPLICIT NONE

!   Declaration of variables

       REAL    :: AK, BK, CK, DK, F, FD, P, PK, PRAT, QRT,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       IF(P.LE.PK)THEN

!      Rarefaction wave
          PRAT = P/PK
          F    = G4*CK*(PRAT**G1 - 1.0)
          FD   = (1.0/(DK*CK))*PRAT**(-G2)
       ELSE

!      Shock wave
          AK  = G5/DK
          BK  = G6*PK
          QRT = SQRT(AK/(BK + P))
          F   = (P - PK)*QRT
          FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
       ENDIF
       RETURN
     END SUBROUTINE PREFUN



     SUBROUTINE SAMPLE(PM, UM, S, D, U, P, V, HAVG)

!   Purpose: to sample the solution throughout the wave 
!            pattern. Pressure PM and velocity UM in the 
!            Star Region are known. Sampling is performed
!            in terms of the 'speed' S = X/T. Sampled
!            values are D, U, P

!   Input variables : PM, UM, S, /GAMMAS/, /STATES/
!   Output variables: D, U, P

       IMPLICIT NONE

!   Declaration of variables

       REAL    DL, UL, PL, CL, DR, UR, PR, CR,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,&
            &        C, CML, CMR, D, P, PM, PML, PMR, S,&
            &        SHL, SHR, SL, SR, STL, STR, U, UM, VL, VR, V,&
            &        HAVG

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR
      
       IF(S.LE.UM)THEN

!      Sampling point lies to the left of the contact
!      discontinuity
          V = VL

          IF(PM.LE.PL)THEN

!         Left rarefaction

             SHL = UL - CL

             IF(S.LE.SHL)THEN

!            Sampled point is left data state

                D = DL
                U = UL
                P = PL

             ELSE
                CML = CL*(PM/PL)**G1
                STL = UM - CML

                IF(S.GT.STL)THEN

!               Sampled point is Star Left state

                   D = DL*(PM/PL)**(1.0/GAMMA)
                   U = UM
                   P = PM
                ELSE

!               Sampled point is inside left fan

                   U = G5*(CL + G7*UL + S)
                   C = G5*(CL + G7*(UL - S))
                   D = DL*(C/CL)**G4
                   P = PL*(C/CL)**G3

                   P = PL*(2.*(1.-HAVG)/(GAMMA-2.*HAVG+1.) + (GAMMA-1.)&
                        /((GAMMA-2.*HAVG+1.)*CL)*((1.-HAVG)*UL))**G3
                   D = DL*(P/PL)**(1./GAMMA)
                   U = UL - 2.*CL/(GAMMA-1.)*((P/PL)**((GAMMA-1.)/(2.*GAMMA))-1.)
                   C = SQRT(GAMMA*P/D)

                ENDIF
             ENDIF
          ELSE

!      Left shock

             PML = PM/PL
             SL  = UL - CL*SQRT(G2*PML + G1)

             IF(S.LE.SL)THEN

!            Sampled point is left data state

                D = DL
                U = UL
                P = PL
                
             ELSE

!            Sampled point is Star Left state

                D = DL*(PML + G6)/(PML*G6 + 1.0)
                U = UM
                P = PM
             ENDIF
          ENDIF
       ELSE

!      Sampling point lies to the right of the contact 
!      discontinuity
          V = VR
          IF(PM.GT.PR)THEN
!         Right shock

             PMR = PM/PR
             SR  = UR + CR*SQRT(G2*PMR + G1)

             IF(S.GE.SR)THEN
               
!            Sampled point is right data state
                D = DR
                U = UR
                P = PR
             ELSE

!            Sampled point is Star Right state

                D = DR*(PMR + G6)/(PMR*G6 + 1.0)
                U = UM
                P = PM
             ENDIF
          ELSE

!         Right rarefaction
             SHR = UR + CR

             IF(S.GE.SHR)THEN
!            Sampled point is right data state
                D = DR
                U = UR
                P = PR
             ELSE
                CMR = CR*(PM/PR)**G1
                STR = UM + CMR

                IF(S.LE.STR)THEN

!               Sampled point is Star Right state
                   D = DR*(PM/PR)**(1.0/GAMMA)
                   U = UM
                   P = PM
                ELSE

!               Sampled point is inside right fan

                   U = G5*(-CR + G7*UR + S)
                   C = G5*(CR - G7*(UR - S))
                   D = DR*(C/CR)**G4
                   P = PR*(C/CR)**G3

                   P = PR*(2.*(1.-HAVG)/(GAMMA-2.*HAVG+1.) - (GAMMA-1.)&
                        /((GAMMA-2.*HAVG+1.)*CR)*((1.-HAVG)*UR))**G3
                   D = DR*(P/PR)**(1./GAMMA)
                   U = UR + 2.*CR/(GAMMA-1.)*((P/PR)**((GAMMA-1.)/(2.*GAMMA))-1.)
                   C = SQRT(GAMMA*P/D)

                ENDIF
             ENDIF
          ENDIF
       ENDIF
       
       RETURN
     END SUBROUTINE SAMPLE
     
   end module riemann

  module Wtemp_maker_mod
  contains
    subroutine Wtemp_maker( W , E , h , x , itrue , jtrue , ishift &
         , jshift , n , area , Wtempl , Wtempr , ramp , muscl_update , h0 )
      use riemann
      use muscl_mod
      use boundary_conditions_mod
      use states2
      use projection_routines
      implicit none
      real    , intent(in) :: W(:,:,:) , E(:,:,:) , h(:,:) , area , x(:,:) , h0
      integer , intent(in) :: itrue , jtrue , n , ishift , jshift
      real    , intent(out) :: Wtempl(:) , Wtempr(:)
      logical , intent(in)  :: ramp , muscl_update
      real :: El_met(4) , Er_met(4) , Wl(4) , Wr(4) , hl , hr , metric(2,2)&
           , Wtemp(2,4) 
      real , allocatable    :: Wtemp2(:,:,:) , Etemp(:,:,:) , htemp(:,:)
      integer :: m , i , j , neta , nxi
      
      neta = size(E,1) ; nxi = size(E,2)

      if( n .ne. 2 )then
         allocate( Wtemp2(neta,1,4) , Etemp(neta,1,8) , htemp(neta,1) )
         call boundary_conditions_prim( Wtemp2 )
         call boundary_conditions( Etemp , htemp , .true. , h0) 
      end if

      do m = 1 , 2

         if( m .eq. 1 )then
            i = itrue          ; j = jtrue          
         else
            i = itrue + ishift ; j = jtrue + jshift
         end if

         call left_and_right_states(i,j,n,x,E,h,Er_met,hr,Wr,El_met,hl,Wl,h0,ramp)
         call normal_vector( El_met , Er_met , n , metric )
         call velocity_projector( Wl , metric )
         call velocity_projector( Wr , metric )
         call riemann_solve( area , Wl , Wr , Wtemp(m,:) , 0.5*(hl+hr) )
         call velocity_deprojector( Wtemp(m,:) , metric )
      end do
      Wtempl = Wtemp(1,:) ; Wtempr = Wtemp(2,:)
      if( n .ne. 2 ) deallocate( Wtemp2 )
    
    end subroutine Wtemp_maker
    
  end module Wtemp_maker_mod

module flow_update_mod2
contains
  subroutine flow_update( E , h , dt , dxi , deta , x&
       , qold , ramp , muscl_update , h0 )
    use helper_functions
    use vars
    use window_data
    use riemann
    use Wtemp_maker_mod

    implicit none
    real , intent(in)    :: dt , dxi , deta &
         , h(:,:) , x(:,:) , h0
    real , intent(inout) :: E(:,:,:) 
    real , intent(out)   :: qold(:,:,:)
    logical , intent(in) :: ramp , muscl_update
    integer :: n , i , j , neta , nxi , jshift , ishift , js , is
    real , allocatable   :: W(:,:,:) 
    real :: area , dlambda , Wl(4) , Wr(4) , El_met(4) , Er_met(4)&
         , hl , hr , havg , metric(2,2) , Fl(4) , Fr(4) , vol_inv &
         , theta , pi , Wtempl(4) , Wtempr(4)
    neta = size(E,1) ; nxi = size(E,2) 
    vol_inv = 1./(deta*dxi)
    allocate( W(neta,nxi,4) ) ; W = 0.
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
    qold(:,:,:) = W(:,:,2:3)
    do n = 1 , 3
       do i = 1 , nxi 
          do j = 1 , neta 
             jshift = jshift_func(n) ; ishift = ishift_func(n)
             area = area_func( deta , dxi , n ) 
             dlambda = dlambda_func( dt , n )
             call Wtemp_maker( W , E , h , x , i , j , ishift  &
                  , jshift , n , area , Wtempl , Wtempr , ramp &
                  , muscl_update , h0 )
             call geom_vars_update( E(j,i,5:8) , h(j,i) , Wtempl &
                  , Wtempr , n , dlambda , vol_inv , area )
             call fluxes( E(j,i,5:8) , h(j,i) , Wtempl , n , Fl )
             call fluxes( E(j,i,5:8) , h(j,i) , Wtempr , n , Fr )
             E(j,i,1:4) = E(j,i,1:4) - dlambda*vol_inv*area*(Fr-Fl)
             call constoprimvec( E(j,i,1:4) , W(j,i,:) , E(j,i,5:8) )
          end do
       end do
    end do
    deallocate( W )
  contains
  end subroutine flow_update

  subroutine geom_vars_update( geom , h , Wl , Wr , n , dlambda &
       , vol_inv , area )
    implicit none
    real    , intent(inout) :: geom(:) 
    real    , intent(in)    :: h , Wl(:) , Wr(:) , dlambda , area &
         , vol_inv
    integer , intent(in)    :: n
    
    if( n .ne. 2 )then 
       geom(1) = geom(1) + dlambda*vol_inv*area*h*( Wr(2) - Wl(2) )
       geom(2) = geom(2) + dlambda*vol_inv*area*h*( Wr(3) - Wl(3) )
    else
       geom(3) = geom(3) + dlambda*vol_inv*area*h*( Wr(2) - Wl(2) )
       geom(4) = geom(4) + dlambda*vol_inv*area*h*( Wr(3) - Wl(3) )
    end if
                 
  end subroutine geom_vars_update

  subroutine fluxes( geom , h , W , n , F )
    use physical_data
! This subroutine uses the values from Wtemp to compute the
! intercellular fluxes.
    implicit none
    real    , intent(in)  :: W(:) , geom(:) , h 
    integer , intent(in)  :: n
    real    , intent(out) :: F(:)
    real                  :: I , e

    e = 0.5*( W(2)**2 + W(3)**2 ) + 1./(gamma - 1.)*W(4)/W(1)
    if( n .ne. 2 )then
       I = ( W(2)*geom(4) - W(3)*geom(3) ) 
       F(1) = W(1)*(1.-h)*I
       F(2) = F(1)*W(2) + W(4)*geom(4)
       F(3) = F(1)*W(3) - W(4)*geom(3)
       F(4) = F(1)*e    + W(4)*I
    else
       I = ( W(3)*geom(1) - W(2)*geom(2) ) 
       F(1) = W(1)*(1.-h)*I
       F(2) = F(1)*W(2) - W(4)*geom(2)
       F(3) = F(1)*W(3) + W(4)*geom(1)
       F(4) = F(1)*e    + W(4)*I
    end if

  end subroutine fluxes

end module flow_update_mod2

program two_d

  use window_data
  use global_data
  use primary_variables
  use init
  use write_files_mod
  use geom_update_mod
  use flow_update_mod2
  use h_update_mod 
  use boundary_data

  integer :: nt = 2000 , n = 0 , i , j , skip=100 , a(1) 
  integer :: ny = 100 , nx0 = 60 , nx = 1
  real    :: CFL = .7 , tmax = .5 , h0 = .5
  logical :: ramp = .false. , muscl_update = .false.

! Set window data
  xmin = 0. ; xmax = .6 ; ymin = -0.5 ; ymax =  0.5 ;
! Set boundary data
  boundary_switch = 3 ! {1,2,3} for {uniform,linear,Riemann} inflow, respecitively.


  open(unit = 1, file = 'two_D2.dat')
  write(1,*) ' nt= ' , nt
  write(1,*) ' '
  call basic_init( nx0 , ny , nx , h0 ) ! This initializes the array of conserved variables E.
  call write_files( x , y , E , h , t , dt )

  do while (t .lt. tmax)
     n = n + 1
     if( mod(n,skip) .eq. 0 ) write(*,*) n , t , dt
! Determine dt based on a CFL condition
     dt = timestep()
! Time-advance the flow variables.
     call flow_update( E , h , dt , dxi , deta &
          , x , qold , ramp , muscl_update , h0 )
! Time-advance the coordinates and constrain the solution within a computational window.
     call geom_update() 
! Compute h for next time step.
     hold = h
!     if(n .ne. 1)&
!call h_update( E(:,:,5) , E(:,:,6) , E(:,:,7) , E(:,:,8) , h , E(:,:,2)/E(:,:,1) , E(:,:,3)/E(:,:,1) , xi , eta , h0 )
     t = t + dt
! Write coordinates and flow variables to data file for reading into Matlab.
     if ( mod(n,skip) .eq. 0 )then
        call write_files( x , y , E , h , t , dt )
     end if
  end do
  write(*,*) "Total number of time steps - " , n
  deallocate( x , y , E , h )
contains
  function timestep()
    use physical_data
    use vars
    save
    real :: timestep , specxi , speceta
    integer :: b(1) , m = 0
    logical :: switch
    real , allocatable :: W(:,:,:) , delta(:,:) , a(:,:)
    allocate(W(size(E,1),size(E,2),4) , delta(size(E,1),size(E,2)) &
         , a(size(E,1),size(E,2)) ) 
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
    a = sqrt(gamma*W(:,:,4)/W(:,:,1))
    delta = E(:,:,5)*E(:,:,8) - E(:,:,6)*E(:,:,7)

    specxi = maxval( &
    1./delta*( &
    (1.-h)* abs(   W(:,:,2)*E(:,:,8)      -   W(:,:,3)*E(:,:,7)      ) &
    +   a *sqrt( ( W(:,:,2)*E(:,:,8) )**2 + ( W(:,:,3)*E(:,:,7) )**2 ) &
    ) &
    )

    speceta = maxval( &
    1./delta*( &
    (1.-h)* abs(   W(:,:,2)*E(:,:,6)      -   W(:,:,3)*E(:,:,5)      ) &
    +   a *sqrt( ( W(:,:,2)*E(:,:,6) )**2 + ( W(:,:,3)*E(:,:,5) )**2 ) &
    ) &
    )

    timestep = CFL/( &
         0.5*( specxi/dxi + speceta/deta ) &
         )

    if( maxval(x(:,1)+h(:,1)*E(:,1,2)/E(:,1,1)*dt) .gt. xmin + dxi )then
       b = maxloc(x(:,1)+h(:,1)*E(:,1,2)/E(:,1,1)*dt)
       timestep = (xmin+dxi-x(b(1),1))/(h(b(1),1)*E(b(1),1,2)/E(b(1),1,1))
    elseif( ramp .and. x(1,nxi-m)+h(1,nxi-m)*E(1,nxi-m,2)&
         /E(1,nxi-m,1)*dt .gt. .5 )then
       timestep = (.5 - x(1,nxi-m))*E(1,nxi-m,1)/(h(1,nxi-m)*E(1,nxi-m,2))
       m = m + 1
       switch = .true.
    end if

    deallocate(W,delta,a)
  end function timestep
end program two_d
