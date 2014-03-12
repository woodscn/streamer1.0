module physical_data 
! This module contains only physical parameters such as the specific gas constant.
  real , parameter :: gamma = 1.4 
end module physical_data

module window_data
! This module contains data on the computational "window" described by Hui. 
!   That is, the min and max parameters defined here describe a box within 
!   which the flow is computed. In reality, a UCS computation will always 
!   extend beyond this box, to some extent.
  real , parameter :: xmin = 0. , xmax = 0.6 , ymin = -0.5 , ymax = 0.5 
end module window_data

module global_data
  save
! Data within this module is shared with many parts of the program. 
  real    :: dxi , deta ! Constant cell width in computational space
  real    :: t = 0. , dt ! Time and physical step size
  integer :: neta , nxi ! Number of computational grid points
  real , allocatable :: E(:,:,:) ! Array of conserved variables. Primary data structure.
  real , allocatable :: x(:,:) , y(:,:) ! Arrays of physical coordinates x , y
  real , allocatable :: h(:,:) ! Grid-angle preserving function h
  real , allocatable :: xi(:,:) , eta(:,:) ! Computational coordinates xi , eta
end module global_data

module vars
! This module contains routines for transforming from conservative variables (e.g. 
!   rho*delta*u) to primitive variables (e.g. u). It uses only the constant gamma,
!   and the arrays of conservative and primitive variables to update either 
!   E(:,:,1:4) or W(:,:,:).
contains
  subroutine primtocons( E , W ) 
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
    real , intent(in)    :: W(:,:,:)
    real , intent(inout) :: E(:,:,:)
    real , allocatable :: delta(:,:)
    integer :: neta , nxi ! Recomputed here for the sake of code independence
    neta = size(E,1) ; nxi = size(E,2)
    allocate( delta(neta,nxi) )
    delta = E(:,:,5)*E(:,:,8) - E(:,:,6)*E(:,:,7) ! AM-BL, see documentation.
! Compute conservative variables based on primitive variables
    E(:,:,1) = W(:,:,1)*delta
    E(:,:,2) = E(:,:,1)*W(:,:,2)
    E(:,:,3) = E(:,:,1)*W(:,:,3)
    E(:,:,4) = delta*(W(:,:,1)*0.5*(W(:,:,2)**2 + W(:,:,3)**2) &
         + W(:,:,4)*(1./(gamma-1.)))
    deallocate(delta)
  end subroutine primtocons

  subroutine constoprim( E , W )
! Inputs: E(neta,nxi,8) - the array of conservative variables. The first 4 "rows"
!   in the third index are conservative flow variables, while the last 4 "rows" 
!   are geometric variables.

! Outputs: W(neta,nxi,4) - the array of primitive flow variables rho, u, v, p.
    use physical_data
    implicit none   
    real , intent(in)  :: E(:,:,:) 
    real , intent(out) :: W(:,:,:)
    real , allocatable :: delta(:,:)
    integer :: neta , nxi
    neta = size(E,1) ; nxi = size(E,2)
    allocate( delta(neta,nxi) )
    delta = E(:,:,5)*E(:,:,8) - E(:,:,6)*E(:,:,7) ! AM-BL
! Compute primitive variables from conservative variables
    W(:,:,1) = E(:,:,1)/delta
    W(:,:,2) = E(:,:,2)/E(:,:,1)
    W(:,:,3) = E(:,:,3)/E(:,:,1)
    W(:,:,4) = (gamma-1.)*(E(:,:,4)/delta-0.5*W(:,:,1)*(W(:,:,2)**2+W(:,:,3)**2))
    deallocate(delta)
  end subroutine constoprim
end module vars

module boundary_conditions_mod
! Apply boundary conditions to the array of conservative variables and to h.
contains
  subroutine boundary_conditions( E , h , flag )
! Inputs: E - array of conservative variables, h - shape-preserving function,
!   flag - optional flag, indicates whether or not to apply downstream 
!   conditions. It's presence causes the routine to apply upstream conditions
!   only.
    use vars
    implicit none
    real , intent(inout) :: E(:,:,:) , h(:,:)
    integer              :: ny , nx , j , i 
    logical , intent(in) , optional :: flag 
    real , allocatable   :: Wtemp(:,:,:) , Etemp(:,:,:)
    real                 :: Wl(4) , Wr(4)

    ny = size(E,1) ; nx = size(E,2)
    allocate( Wtemp(ny,1,4) , Etemp(ny,1,8) )
! Define upstream conditions
    do j = 1 , ny
! Define metric components, using a cartesian grid at the inflow boundary.
       Etemp(j,1,:) = (/ 0. , 0. , 0. , 0. , 1. , 0. , 0. , 1. /) 
    end do
    call boundary_conditions_prim ( Wtemp)
! Convert Wtemp (primitive variables) to Etemp (conservative variables).
    call primtocons( Etemp , Wtemp ) 
! Enforce downstream conditions.
    if( .not. present(flag) ) then 
! These are Neumann conditions with zero derivative.
       E(:,nx,:) = E(:,nx-1,:)
    end if
! Apply Dirichlet boundary conditions at upsteram, requiring E = Etemp
    E(:,1,:) = Etemp(:,1,:) ! Apply upstream BC 
! Apply zero-derivative condition to h at upstream boundary
    h(:,1) = h(:,2)
    deallocate( Wtemp , Etemp )
  end subroutine boundary_conditions
  
!          call boundary_conditions_prim( Wtemp(:,1,:) )

  subroutine boundary_conditions_prim( Wtemp )
    real , intent(inout) :: Wtemp(:,:,:)
    integer              :: ny , nx , i , j 
    real                 :: Wl(4) , Wr(4)
    ny = size(Wtemp,1) ; nx = size(Wtemp,2)

!! Basic, uniform flow
!    do j = 1 , ny
!       Wtemp(j,1,:) = (/ 1. , 1. , 0. , 1. /) 
!! Linear velocity inflow
!       Wtemp(j,1,:) = (/ 1. , 1.+real(j)/real(ny) , 0. , 1. /)
!    end do

! The Riemann problem of Hui
    Wl = (/ 1.0 , 2.4*sqrt(1.4*1.00/1.0) , 0. , 1.00 /)
    Wr = (/ 0.5 , 7.0*sqrt(1.4*0.25/0.5) , 0. , 0.25 /)
    do j = 1 , ceiling(real(ny)/2)
       Wtemp(j,1,:) = Wl
    end do
    do j = ceiling(real(ny)/2)+1 , ny
       Wtemp(j,1,:) = Wr
    end do

  end subroutine boundary_conditions_prim

end module boundary_conditions_mod

module init 
! Initialize the flow field. At the moment, this gives an initial guess. If 
!   allowing grid to propogate from upstream, it would be better to 
contains
  subroutine basic_init() ! Initialize field with uniform flow
! Inputs:  global_data: E, h, xi, eta, x, y, t, dt, neta, nxi, deta, dxi
!          window_data: xmin, xmax, ymin, ymax
! Outputs: fully initialized E, h, xi, eta, x, y for a given boundary condition.
!          At the moment, this simply applies the upstream boundary condition to 
!          the entire region. It also defines an overall h of 0.25, and a uniform 
!          xi,eta grid that is initially identical to the x,y grid.
    use vars
    use global_data
    use window_data
    use boundary_conditions_mod
    implicit none
    integer :: i , j , ny = 100 , nx0 = 60 , nx = 60
    real , allocatable :: A(:,:) , B(:,:) , L(:,:) , M(:,:)
    real , allocatable :: W(:,:,:)
    real :: dx , dy
    
!    write(*,*) "Allocating..."
    allocate( &
         W(ny,nx,4) , x(ny,nx) , y(ny,nx) ,&
         E(ny,nx,8) , xi(ny,nx) , eta(ny,nx) , h(ny,nx) )
! Define an initial, uniform mesh in x,y. These arrays give an x,y coordinate
!   for each point in computational space. 
! Determine initial grid spacing given a desired number of computational points
    dx = (xmax-xmin)/(nx0) ; dy = (ymax-ymin)/(ny) 
!    dx = .01 ; dy = .01
! Use identical spacing, etc. for computational grid.
    dxi = dx ; deta = dy ; neta = ny ; nxi = nx
! Use a cell-centered mesh
    x(:,1) = xmin + 0.5*dx
    y(1,:) = ymin + 0.5*dy
! Initialize x, y.
    do i = 2 , nx
       x(:,i) = x(:,i-1) + dx
    end do
    do i = 2 , ny
       y(i,:) = y(i-1,:) + dy
    end do
! Are these necessary?
    xi = x ; eta = y
! Initialize h.
    h = 0.25
! Apply upstream boundary condition to entire region.
    call boundary_conditions( E , h , .true. )
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
    call constoprim( E , W )
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
! This module exists to advance the x,y coordinates in time, to create/destroy
!   points in order to keep the computation within the computational window, 
!   and to solve for updated values of geometric variables.
contains
  subroutine geom_update
! geom_update acts only as the master routine to direct the action of the three
!   individual subroutines.
    call coords_update ! Move the x,y coordinates with the fluid velocity
    call windower ! Create/destroy points to maintain the computational window
  end subroutine geom_update

  subroutine coords_update()
! This routine advances the values of x,y in time according to the formula
!   dx = h*u*dt and dy = h*v*dt. It also updates t.
    use global_data
    implicit none
    integer :: i
    real , allocatable  :: u(:,:) , v(:,:)

! In contrast to most other subroutines, this one has access to global values
!   that include neta and nxi, so they are not recomputed here.
    allocate( u(neta,nxi) , v(neta,nxi) ) 
    u = E(:,:,2)/E(:,:,1);v = E(:,:,3)/E(:,:,1)
    if (maxval(h*u*dt) .ge. dxi .or. maxval(h*v*dt) .ge. deta)then
! If the time step is too large, then x,y components can advance more than one
!   gridpoint forward in a time step. This makes applying BC's difficult.
       write(*,*) "Warning, time step causes coordinates to move too quickly."
       stop
    else
! Advance spatial and temporal coordinates
!       t = t + dt
       x = x + h*u*dt
       y = y + h*v*dt
    end if
    deallocate( u , v )
!! Add a ramp to our simple channel flow at x = 0.
!    do i = 1 , nxi
!       if (x(1,i) .gt. 0.1)then
!          y(1,i) = y(1,1)+(x(1,i)-0.)*0.1
!          E(1,i,5:8) = 1./sqrt(101.)*(/ 10. , 1. , -1. , 10. /)
!       end if
!    end do

  end subroutine coords_update

  subroutine windower()
! Windower is tasked with keeping the solution within a designated computational
!   window. It evaluates whether the grid has overrun the maximum and minimum 
!   required values of x, and creates/destroys columns as necessary. 
    use global_data
    use window_data
    use boundary_conditions_mod
    use vars
    implicit none
    logical :: overrun = .false. , shortfall = .false.
    real , allocatable :: x_t(:,:) , y_t(:,:) , E_t(:,:,:) , h_t(:,:) &
         , xi_t(:,:) , eta_t(:,:) , A(:,:) , B(:,:)
    integer :: j
! In contrast to most other subroutines, this one has access to global values
!   that include neta and nxi, so they are not recomputed here.
    allocate(A(neta,nxi),B(neta,nxi))
    A = E(:,:,5) ; B = E(:,:,6)
! Test whether the grid extends too far beyond the right of computational window
    if(nxi .gt. 1)then
       if (minval(x(:,nxi-1)) .gt. xmax) then
          overrun = .true. 
       else
          overrun = .false.
       end if
    end if
! Test whether the grid does not extend far enough beyond the left of the 
!   computational window
    if(nxi .gt. 1)then
       if (maxval(x(:,2)) .gt. xmin) then
          shortfall = .true.
       else
          shortfall = .false.
       end if
    else
       shortfall = .true.
    end if
! Begin to update coordinates based on the results of the windowing tests above
    if (overrun .and. (.not. shortfall))then
! The grid extends too far in the right (downstream), but not in the left.
!       write(*,*) "Overrun yes, shortfall no."
! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h ; xi_t = xi ; eta_t = eta
! Remove a column and reallocate variable arrays
       call deallocate_real
       nxi = nxi - 1
       call allocate_real
! Copy the temporary arrays back into the variable arrays, minus the last, 
!   "overrun" column.
       x = x_t(:,1:nxi) ; y  = y_t (:,1:nxi) ; E   = E_t  (:,1:nxi,:) 
       h = h_t(:,1:nxi) ; xi = xi_t(:,1:nxi) ; eta = eta_t(:,1:nxi)
       call deallocate_temp
    elseif (shortfall .and. (.not. overrun)) then
! The grid does not extend far enough to the left (upstream), but the 
!   right is fine.
!       write(*,*) "Shortfall yes, overrun no."
! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h ; xi_t = xi ; eta_t = eta
! Add a column and reallocate variable arrays
       call deallocate_real
       nxi = nxi + 1
       call allocate_real
! Copy all but the first column from the temporary arrays.
       x(:,2:nxi) = x_t ; y (:,2:nxi) = y_t  ; E  (:,2:nxi,:) = E_t 
       h(:,2:nxi) = h_t ; xi(:,2:nxi) = xi_t ; eta(:,2:nxi)   = eta_t
       call deallocate_temp
! Call boundary_conditions to compute the first column for E,h.
!       write(*,*) "Writing BC's"
       call boundary_conditions( E , h , .true. )
!       write(*,*) "Done"
! Update the first column of the remaining arrays.
       eta(:,1) = eta(:,2) ; xi(:,1) = xi(:,2)-dxi
       x(:,1) = x(:,2) - 0.5*dxi*( E(:,1,5) + E(:,2,5) )
       y(:,1) = y(:,2) - 0.5*dxi*( E(:,1,6) + E(:,2,6) )
    elseif (shortfall .and. overrun) then
! The grid has the appropriate length; it just needs to be shifted upstream.
!       write(*,*) "Both overrun and shortfall."
! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h ; xi_t = xi ; eta_t = eta
! Overwrite all but the first column with data from the temporary arrays
       x(:,2:nxi)   = x_t(:,1:(nxi-1))   ; y(:,2:nxi)   = y_t  (:,1:(nxi-1))
       E(:,2:nxi,:) = E_t(:,1:(nxi-1),:) ; h(:,2:nxi)   = h_t  (:,1:(nxi-1))
       xi(:,2:nxi)  = xi_t(:,1:(nxi-1))  ; eta(:,2:nxi) = eta_t(:,1:(nxi-1))
       call deallocate_temp
! Call boundary_conditions to compute the first column for E,h.
       call boundary_conditions( E , h , .true. )
! Update the first column of the remaining arrays.
       eta(:,1) = eta(:,2) ; xi(:,1) = xi(:,2)-dxi
       x(:,1) = x(:,2) - 0.5*dxi*( E(:,1,5) + E(:,2,5) )
       y(:,1) = y(:,2) - 0.5*dxi*( E(:,1,6) + E(:,2,6) )
    else
!       write(*,*) "Error in overrun case selection."
!       stop
    endif
    deallocate(A,B)
!    write(*,*) "Windower Done"   
    contains
      subroutine allocate_temp ! Allocate temporary arrays
!        write(*,*) "Allocating temp..."
        allocate( x_t(neta,nxi) , y_t(neta,nxi) , E_t(neta,nxi,8) ,&
             h_t(neta,nxi) , xi_t(neta,nxi) ,eta_t(neta,nxi) )
      end subroutine allocate_temp

      subroutine allocate_real ! Allocate variable arrays
!        write(*,*) "Allocating real..."
        allocate( x(neta,nxi) , y(neta,nxi) , E(neta,nxi,8) ,&
             h(neta,nxi) , xi(neta,nxi) ,eta(neta,nxi) )
      end subroutine allocate_real

      subroutine deallocate_temp ! Deallocate temporary arrays
!       write(*,*) "Deallocating temp..."
        deallocate( x_t , y_t , E_t , h_t , xi_t , eta_t )
      end subroutine deallocate_temp

      subroutine deallocate_real ! Deallocate variable arrays
!        write(*,*) "Deallocating real..."
        deallocate( x , y , E , h , xi , eta )
      end subroutine deallocate_real

  end subroutine windower
  
end module geom_update_mod

module h_update_mod
! This module is responsible for solving the PDE for h. It uses an iterative, 
!   first-order, backward difference scheme. Boundary values must be provided 
!   at h(1,:) and h(:,1).
contains

  subroutine h_update( A , B , L , M , h , u , v , xi , eta )
    use boundary_conditions_mod
    implicit none
    real , intent(in)   :: A(:,:) , B(:,:) ,  L(:,:) ,   M(:,:)
    real , intent(in)   :: u(:,:) , v(:,:) , xi(:,:) , eta(:,:)
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
       err = maxval(abs(temp-g))
       g   = temp
    end do
! Convert the updated g back to h.
    h = exp(g)/q

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

module riemann
contains
  subroutine riemann_solve( WinL , WinR , Wout , h )
    implicit none

    real , intent(in)  :: WinL(:) , WinR(:) , h
    real , intent(out) :: Wout(:)
    real :: eta , xi , gamma_hat , gamma_coeff 
    real :: DL , UL , VL , PL , DR , UR , VR , PR , AL , AR
    real :: Pstar , Ustar , DstarL , DstarR , cL , cR , betaL , betaR
    real :: PsiL , PsiR , UstarL , UstarR , dcdpsi , dbetadpsi 
    real :: dUstarL , dUstarR , temp = 1. , cTL , cHL , cTR , cHR
    real :: Dout , Uout , Vout , Pout , tol = .000001 , gL , gR
    real , parameter :: gamma = 1.4
    integer :: n = 1
    eta = 1.-h
    xi  = sqrt( 1. - h/gamma )
    gamma_hat = gamma - 1.
    gamma_coeff = 2.*eta + gamma_hat

    DL = WinL(1) ; UL = WinL(2) ; VL = WinL(3) ; PL = WinL(4)
    DR = WinR(1) ; UR = WinR(2) ; VR = WinR(3) ; PR = WinR(4)
    AL = sqrt(gamma*PL/DL) ; AR = sqrt(gamma*PR/DR)

! Linearised guess
    Pstar = .5*(PL+PR)-.125*(uR-uL)*(DL+DR)*(aL+aR)
 
!! Two-shock solution
!   gL=sqrt( 2./((gamma_coeff)*DL)/(Pstar+PL*gamma_hat/gamma_coeff))
!   gR=sqrt( 2./((gamma_coeff)*DR)/(Pstar+PR*gamma_hat/gamma_coeff))
!   Pstar = max(tol,&
!         (gL*PL+gR*PR-(UR-UL))/(gL+gR)&
!         )

!! Two-rarefaction solution
!    Pstar = (&
!         (aL+aR-.5*(gamma_hat)*(uR-uL))&
!         /(aL/pL**(gamma_hat/(2.*gamma))+aR/pR**(gamma_hat/(2.*gamma)))&
!         )**(2.*gamma/(gamma_hat))

!    Pstar = .5*( PL + PR )

    do while(abs(temp) .gt. tol)
       write(*,*) n ; n = n + 1
       PsiL = Pstar/PL
       betaL = beta(psiL,gamma_hat,gamma_coeff,eta)
       if( PsiL .gt. 1. )then
          cL = UL - AL*sqrt( betaL*(PsiL-1.) &
         / (gamma*eta*(betaL-1.)) )
          UstarL = cL + (UL - cL)/betaL
          dcdpsi = (cL - UL)/(2.*(PsiL-1.))
          dbetadpsi = ( gamma_coeff**2 - gamma_hat**2 )&
               /( gamma_coeff + gamma_hat*PsiL )**2
          dUstarL = 1./PL*((1.-1./betaL)*dcdpsi-1./betaL**2*(uL-cL)*dbetadpsi)
       else
          UstarL = UL - 2.*AL/(xi*gamma_hat)&
               *(PsiL**(gamma_hat/(2.*(gamma_hat+eta)))-1.)
          dUstarL = -aL*gamma_hat/(PL*xi*gamma*(gamma_hat+eta))&
               *PsiL**(-gamma_coeff/(2.*(gamma_hat+eta)))
       end if

       PsiR = Pstar/PR
       betaR = beta(psiR,gamma_hat,gamma_coeff,eta)
       if( PsiR .gt. 1. )then
          cR = UR + AR*sqrt( betaR*(PsiR-1.) &
         / (gamma*eta*(betaR-1.)) )
          UstarR = cR + (UR - cR)/betaR
          dcdpsi = (cR - UR)/(2.*(PsiR-1.))
          dbetadpsi = (gamma_coeff**2-gamma_hat**2)&
               /(gamma_coeff+gamma_hat*PsiR)**2
          dUstarR = 1./PR*((1.-1./betaR)*dcdpsi-1./betaR**2*(uR-cR)*dbetadpsi)
       else
          UstarR = UR + 2.*AR/(xi*gamma_hat)&
               *(PsiR**(gamma_hat/(2.*(gamma_hat+eta)))-1.)
          dUstarR =  aR*gamma_hat/(PR*xi*gamma*(gamma_hat+eta))&
               * PsiR**(-gamma_coeff/(2.*(gamma_hat+eta)))
       end if
       temp = ( UstarR - UstarL )/( dUstarR - dUstarL )
!       Pstar = max(Pstar - temp,tol)
       Pstar = Pstar - temp
    end do
    Ustar = UstarL
    DstarL = betaL*DL
    DstarR = betaR*DR

!    write(*,*) Pstar , Ustar , DstarL , DstarR

    if( Ustar .gt. 0. )then
       Vout = VL
       if( Pstar .gt. PL )then
          if( cL .gt. 0. )then
             Dout = DL
             Uout = UL
             Pout = PL
          else
             Dout = DstarL
             Uout = Ustar
             Pout = Pstar
          end if
       else
          ctL = sqrt(gamma*PL/DL)
          if( ctL .gt. 0. )then
             Dout = DL
             Uout = UL
             Pout = PL
          else
             chL = sqrt(gamma*Pstar/DstarL)
             if( chL .lt. 0. )then
                Dout = DstarL
                Uout = Ustar
                Pout = Pstar
             else
                Dout = DL*(2./(gamma+1.)-gamma_hat/((gamma+1.)*xi*AL)&
                     *(-eta*UL))**(2./gamma_hat)
                Uout = 2./(gamma+1.)*(gamma_hat/2.*UL+1./eta*(-xi*AL))
                Pout = PL*(2./(gamma+1.)-gamma_hat/((gamma+1)*xi*AL)&
                     *(-eta*UL))**(2.*gamma/gamma_hat)
             end if
          end if
       end if
    else
       Vout = VR
       if( Pstar .gt. PR )then
          if( cR .lt. 0. )then
             Dout = DR
             Uout = UR
             Pout = PR
          else
             Dout = DstarR
             Uout = Ustar
             Pout = Pstar
          end if
       else
          ctR = sqrt(gamma*PR/DR)
          if( ctR .lt. 0. )then
             Dout = DR
             Uout = UR
             Pout = PR
          else
             chR = sqrt(gamma*Pstar/DstarR)
             if( chR .gt. 0. )then
                Dout = DstarR
                Uout = Ustar
                Pout = Pstar
             else
                Dout = DR*(2./(gamma+1.)+gamma_hat/((gamma+1.)*xi*AL)&
                     *(-eta*UL))**(2./gamma_hat)
                Uout = 2./(gamma+1.)*(gamma_hat/2.*UL+1./eta*(+xi*AL))
                Pout = PL*(2./(gamma+1.)+gamma_hat/((gamma+1)*xi*AL)&
                     *(-eta*UL))**(2.*gamma/gamma_hat)
             end if
          end if
       end if
    end if

    Wout = (/ Dout , Uout , Vout , Pout /)

  end subroutine riemann_solve

  function beta(psi,gamma_hat,gamma_coeff,eta)
    implicit none
    real :: psi , gamma_hat , eta , gamma_coeff , beta
    if( psi .gt. 1 )then
       beta = ( gamma_coeff*psi + gamma_hat )&
             /( gamma_coeff + gamma_hat*psi )
    else
       beta = psi**(eta/(gamma_hat+eta))
    end if
  end function beta







module flow_update_mod
! The purpose of this module is to time-advance the flow variables using an 
!   explicit, Godunov scheme. 
contains
  subroutine geom_vars_update2( E , h , Wtemp , n , dlambda , dxi , deta )
    real    , intent(in)    :: h(:,:) , Wtemp(:,:,:) , dlambda , dxi , deta
    integer , intent(in)    :: n
    integer                 :: i , j 
    real    , intent(inout) :: E(:,:,:)
    real    , allocatable   :: flux(:,:,:)
    neta = size(E,1) ; nxi = size(E,2)
    allocate( flux(neta,nxi,2) )
    flux=0.
    if( n .ne. 2 )then
       flux(:,:,1) = h(:,:)*(Wtemp(:,2:nxi+1,2) - Wtemp(:,1:nxi,2))
       flux(:,:,2) = h(:,:)*(Wtemp(:,2:nxi+1,3) - Wtemp(:,1:nxi,3))
       E (:,2:nxi,5:6) = E(:,2:nxi,5:6) + dlambda/dxi *flux(:,2:nxi,:)
       E(:,nxi,5:6) = E(:,nxi-1,5:6)
    elseif( n .eq. 2 )then
       flux(:,:,1) = h(:,:)*(Wtemp(2:neta+1,:,2) - Wtemp(1:neta,:,2))
       flux(:,:,2) = h(:,:)*(Wtemp(2:neta+1,:,3) - Wtemp(1:neta,:,3))
       E (:,2:nxi,7:8) = E(:,2:nxi,7:8) + dlambda/deta*flux(:,2:nxi,:)
       E(:,nxi,7:8) = E(:,nxi-1,7:8)
    end if
!    E(:,nxi,5:8) = E(:,nxi-1,5:8)
    deallocate( flux )
  end subroutine geom_vars_update2
  subroutine geom_vars_update( E , h , Wtemp , n , dlambda , dxi , deta &
       , x , xmin)
    implicit none
    real    , intent(inout) :: E(:,:,:) 
    real    , intent(in)    :: h(:,:) , Wtemp(:,:,:) , x(:,:) &
         , dlambda , dxi , deta , xmin
    integer , intent(in)    :: n
    integer                 :: neta , nxi , i , j
    
    neta = size(E,1) ; nxi = size(E,2)
    do j = 1 , neta
       do i = 1 , nxi
          if( x(j,i) .gt. xmin )then
             if( n .ne. 2 )then 
                E(j,i,5) = E(j,i,5) + dlambda/dxi * h(j,i)* &
                     ( Wtemp(j,i+1,2) - Wtemp(j,i,2) )
                E(j,i,6) = E(j,i,6) + dlambda/dxi * h(j,i)* &
                     ( Wtemp(j,i+1,3) - Wtemp(j,i,3) )
             else
                E(j,i,7) = E(j,i,7) + dlambda/deta * h(j,i)* &
                     ( Wtemp(j+1,i,2) - Wtemp(j,i,2) )
                E(j,i,8) = E(j,i,8) + dlambda/deta * h(j,i)* &
                     ( Wtemp(j+1,i,3) - Wtemp(j,i,3) )
             end if
          end if
       end do
    end do
                 
  end subroutine geom_vars_update

  subroutine flow_update( E , h , xi , eta , dt , dxi , deta , x )
    use window_data
    use vars
    use riemann
    use boundary_conditions_mod
    use muscl_mod
    use helper_functions
    implicit none
    real , intent(in)   :: xi(:,:) , eta(:,:) , dt , dxi , deta , h(:,:) , x(:,:)
    real , intent(inout):: E(:,:,:)
    real , allocatable  :: W(:,:,:) , F(:,:,:) , Wtemp(:,:,:) , delta(:,:) &
         , Etemp(:,:,:) , htemp(:,:)
    integer :: neta , nxi , n , i , j , is , js , ishift , jshift
    real    :: vol_inv , Wl(4) , Wr(4) , havg , temp(2) , Eavg(8) , dlambda , area &
         , metric(2,2) , Fr(4) , Fl(4)
    logical :: tag
    neta = size(E,1) ; nxi = size(E,2) ! Compute grid parameters
    vol_inv = 1./(dxi*deta) ! Compute useful parameter
    allocate( W(neta,nxi,4) )
    do n = 1 , 3 
! The mess of ishifts, jshifts, and other things like area() is used to condense
!   the 3-step Strang splitting technique into a single loop of code.
       jshift = jshift_func(n) ; ishift = ishift_func(n)
       area = area_func( deta , dxi , n ) ; dlambda = dlambda_func( dt , n )
       call constoprim( E , W ) ! The Riemann solver uses primitive variables. 
       allocate( Wtemp(neta+jshift,nxi+ishift,4) )
       Wtemp = 0. 
! March over the interior points.
       do j = 1+jshift , neta
          do i = 1+ishift , nxi
             js = j-jshift ; is = i-ishift
! Use MUSCL interpolation to improve the accuracy away from the edges.
             if( js-jshift .gt. 0 .and. is-ishift .gt. 0 &
                  .and. j+jshift .le. neta .and. i+ishift .le. nxi )then
                call MUSCL( W(j,i,:) , W(j+jshift,i+ishift,:) , W(js,is,:) &
                     , W(js-jshift,is-ishift,:) , Wl , Wr )
             else
                Wl = W(js,is,:) ; Wr = W(j,i,:)
             end if
! Riemann_solve needs input states in the form (/ rho, u_normal, u_tangential, p /)
! Compute the normal vector between the adjacent cells (average of metric quantities).
             call normal_vector( E(js,is,:) , E(j,i,:) , n , metric )
! Project the velocity vector onto normal and tangential components.
             call velocity_projector( Wl , metric )
             call velocity_projector( Wr , metric )
! Call the riemann solver to generate Wtemp.
             havg = 0.5*(h(js,is)+h(j,i))
             call riemann_solve( area , Wl , Wr , Wtemp(j,i,:) , havg )
! De-project the velocity vector
             call velocity_deprojector( Wtemp(j,i,:) , metric )
          end do
       end do

! Apply boundary conditions to the fluxes
       if( n .ne. 2 )then
          call boundary_conditions_prim( Wtemp )
          Wtemp(:,nxi+1,:) = Wtemp(:,nxi,:)
       elseif( n .eq. 2 )then   
          Wtemp(1,:,:)      = Wtemp(2,:,:)    ; Wtemp(1,:,3)      = -Wtemp(1,:,3)
          metric(1,:) = (/ 10. , 1. /) ; metric(2,:) = (/ 1. , 10. /)
          metric = 1./sqrt(101.)*metric
          Wtemp(neta+1,:,:) = Wtemp(neta,:,:) 
          Wtemp(neta+1,:,3) = -Wtemp(neta+1,:,3)
!          do i = 1 , nxi
!             call velocity_projector( W(neta,i,:) , metric )
!             Wtemp(neta+1,i,3) = -Wtemp(neta+1,i,3)
!             call velocity_deprojector( W(neta,i,:) ,metric )
!          end do
       end if
       call geom_vars_update( E , h , Wtemp , n , dlambda , dxi , deta &
            , x , xmin)
! Update the solution       
       do j = 1 , neta
          do i = 2 , nxi
             if( x(j,i) .gt. xmin )then
                call fluxes( Wtemp(j,i,:) , E(j,i,5:8) , h(j,i) , n , Fl )
                call fluxes( Wtemp(j+jshift,i+ishift,:) , E(j,i,5:8) , h(j,i) , n , Fr )
                E(j,i,1:4) = E(j,i,1:4) - dlambda*vol_inv*(Fr-Fl)*&
                     area
             end if
          end do
       end do
       deallocate( Wtemp )
       E(:,nxi,:) = E(:,nxi-1,:)
   end do
   deallocate( W )

  end subroutine flow_update

  subroutine normal_vector(El , Er , n , metric )
    real    , intent(in)  :: El(:) , Er(:) 
    integer , intent(in)  :: n
    real    , intent(out) :: metric(:,:)
    real                  :: normal(2) , tangential(2)
    real                  :: delta , psi
    logical               :: tick

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
          normal = (/ El(8) , -El(7) /)/( El(5)*El(8) - El(6)*El(7) ) &
               + (/ Er(8) , -Er(7) /)/( Er(5)*Er(8) - Er(6)*Er(7) )
          normal = normal/sqrt( normal(1)**2 + normal(2)**2 )
          
          tangential = (/ El(7) , El(8) /)/( El(5)*El(8) - El(6)*El(7) ) &
               + (/ Er(7) , Er(8) /)/( Er(5)*Er(8) - Er(6)*Er(7) )
          tangential = tangential/sqrt( tangential(1)**2 + tangential(2)**2 )
       else
          normal = (/ -El(6) ,  El(5) /)/( El(5)*El(8) - El(6)*El(7) ) &
               + (/ -Er(6) ,  Er(5) /)/( Er(5)*Er(8) - Er(6)*Er(7) )
          normal = normal/sqrt( normal(1)**2 + normal(2)**2 )
          
          tangential = (/ -El(5) , -El(6) /)/( El(5)*El(8) - El(6)*El(7) ) &
               + (/ -Er(5) , -Er(6) /)/( Er(5)*Er(8) - Er(6)*Er(7) )
          tangential = tangential/sqrt( tangential(1)**2 + tangential(2)**2 )
       end if
       
       metric(1,:) = normal
       metric(2,:) = tangential
    end if

  end subroutine normal_vector

  subroutine velocity_projector( W , metric )
    real , intent(inout) :: W(:)
    real , intent(in)    :: metric(:,:)
    real                 :: temp(2)
    temp = (/ W(2) , W(3) /)
    
    temp = matmul(metric,temp)
    W(2) = temp(1) ; W(3) = temp(2)

!    W(2) = matmul(temp , norm)
!    W(3) = matmul(temp , tang)
    
  end subroutine velocity_projector

  subroutine velocity_deprojector( W , metric )
    real , intent(out)   :: metric(:,:)
    real , intent(inout) :: W(:)
    real                 :: temp(2) , metric_temp(2,2)
    temp = (/ W(2) , W(3) /)
! Invert the metric
    metric_temp(1,:) = (/ metric(2,2) , -metric(1,2) /)
    metric_temp(2,:) = (/-metric(2,1) ,  metric(1,1) /)
! Apply the inverse of the metric to the velocity vector
    temp = matmul(metric_temp,temp)
    W(2) = temp(1) ; W(3) = temp(2)
  end subroutine velocity_deprojector

  subroutine fluxes( W , geom , h , n , F )
    use physical_data
! This subroutine uses the values from Wtemp to compute the intercellular fluxes.
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
       F(4) = F(1)*e + W(4)*I
    else
       I = ( W(3)*geom(1) - W(2)*geom(2) ) 
       F(1) = W(1)*(1.-h)*I
       F(2) = F(1)*W(2) - W(4)*geom(2)
       F(3) = F(1)*W(3) + W(4)*geom(1)
       F(4) = F(1)*e + W(4)*I
    end if

  end subroutine fluxes

end module flow_update_mod    

program two_d

  use global_data
  use init
  use write_files_mod
  use geom_update_mod
  use flow_update_mod
  use h_update_mod 

  integer :: nt = 2000 , n = 0 , i , j , skip=10
  real    :: CFL = .7 , tmax = 1.50
  real , allocatable :: u(:,:) , v(:,:)

  open(unit = 1, file = 'two_D.dat')
  write(1,*) ' nt= ' , nt
  write(1,*) ' '
!  write(*,*) "Initializing..."
  call basic_init() ! This initializes the array of conserved variables E.
!  write(*,*) "Initialization complete, entering time-advancement loop"
!     call write_files(x,y,E,t,dt)

  do while (t .lt. tmax)
     n = n + 1
     write(*,*) n !, t
! Determine dt based on a CFL condition
     dt = timestep()
!     write(*,*) "dt = " , dt
!     write(*,*) "Updating flow variables..."
! Time-advance the flow variables.
     call flow_update( E , h , xi , eta , dt , dxi , deta , x )
!     write(*,*) "Updating geometric variables..."
! Time-advance the coordinates and constrain the solution within a computational
! window.
     call geom_update 
!     write(*,*) "Updating h..."
! Compute h for next time step.
     call h_update( E(:,:,5) , E(:,:,6) , E(:,:,7) , E(:,:,8) , h ,&
          E(:,:,2)/E(:,:,1) , E(:,:,3)/E(:,:,1) , xi , eta )
     t = t + dt
!     write(*,*) "Writing files..."
! Write coordinates and flow variables to data file for reading into Matlab.
     if (mod(n,skip) .eq. 0)then
        call write_files(x,y,E,h,t,dt)
     end if
!     write(*,*) "Done"
  end do
  write(*,*) "Total number of time steps - " , n
  deallocate( x , y , xi , eta , E , h )
contains
  function timestep()
    use physical_data
    use vars
    real :: timestep
    real , allocatable :: W(:,:,:)
    allocate(W(size(E,1),size(E,2),4))
    call constoprim(E,W)
    timestep = min(dxi,deta)/(maxval(abs(W(:,:,2)) + sqrt(gamma*W(:,:,4)/W(:,:,1))))*CFL
    deallocate(W)
  end function timestep
end program two_d
