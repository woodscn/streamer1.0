program two_d

  use window_data
  use global_data
  use primary_variables
  use init
  use write_files_mod
  use geom_update_mod
  use flow_update_mod
  use h_update_mod 
  use boundary_conditions_mod

  implicit none

  integer :: nt = 2000 , n = 0 , i , j , skip=100 , a(1)
  integer :: ny = 100 , nx0 = 60 , nx = 1
  real    :: CFL = .7 , tmax = 1.0 , h0 = .25
  logical :: ramp = .false. , muscl_update = .false. , geom_switch = .false. , &
  	h_switch = .true.

! Set window data
  xmin = 0. ; xmax = .9 ; ymin = -0.5 ; ymax =  0.5 ;
! Set boundary data
  boundary_switch = 3 ! {1,2,3} for {uniform,linear,Riemann} inflow, respecitively.

  open(unit = 1, file = 'two_D2.dat')
    open(unit = 2, file = 'walldebug.dat')
  write(1,*) ' nt= ' , nt
  write(1,*) ' '
  call basic_init( nx0 , ny , nx , h0 ) ! This initializes the array of conserved variables E.
!  call channel_init( nx0 , ny , nx , h0 ) ! This initializes the array of conserved variables E.
  call write_files( x , y , E , h , t , dt )
  do while (t .lt. tmax)
     n = n + 1
     if( mod(n,skip) .eq. 0 ) write(*,*) n , t , dt
! Determine dt based on a CFL condition
     dt = .001!max(timestep(),1*10**(-10.))
! Time-advance the flow variables.
     call flow_update( E , h&
          , x , y , qold , ramp , muscl_update , h0 )
! Time-advance the coordinates and constrain the solution within a computational window.
     call geom_update(geom_switch) 
! Compute h for next time step.
     hold = h
!     if(nxi .ne. 1)&

  if(size(E,2) .gt. 1 .and. h_switch )then
     call h_update( E(:,:,5) , E(:,:,6) , E(:,:,7) , E(:,:,8) , h &
        , E(:,:,2)/E(:,:,1) , E(:,:,3)/E(:,:,1) &!, xi , eta
        , h0 )
  end if
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
