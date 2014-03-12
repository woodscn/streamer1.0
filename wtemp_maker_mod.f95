  module Wtemp_maker_mod
    implicit none
  type :: state
     real :: phys(4)
     real :: geom(4)
     real    :: h
  end type state
contains
    subroutine Wtemp_maker( W , E , h , x , itrue , jtrue , ishift &
         , jshift , n , area , Wtempl , Wtempr , ramp , muscl_update , h0 )
      use global_data
      use riemann
      use muscl_mod
      use boundary_conditions_mod
!      use states
!      use projection_routines
      implicit none
      real    , intent(in) :: W(neta,nxi,4) , E(neta,nxi,8) , h(neta,nxi) , area , x(neta,nxi) , h0
      integer , intent(in) :: itrue , jtrue , n , ishift , jshift
      real    , intent(out) :: Wtempl(4) , Wtempr(4)
      logical , intent(in)  :: ramp , muscl_update
      real :: El_met(4) , Er_met(4) , Wl(4) , Wr(4) , hl , hr , metric(2,2)&
           , Wtemp(2,4) 
      real , allocatable    :: Wtemp2(:,:,:) , Etemp(:,:,:) , htemp(:,:)
      integer :: m , i , j 
      type (state)          :: left_state , right_state

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

         call left_and_right_states(i,j,n,x,E,h&
              ,left_state,right_state,h0,ramp,muscl_update)
         Wl = left_state%phys
         El_met = left_state%geom
         hl = left_state%h
         Wr = right_state%phys
         Er_met = right_state%geom
         hr = right_state%h
         call normal_vector( El_met , Er_met , n , metric )
         call velocity_projector( Wl , metric )
         call velocity_projector( Wr , metric )
         call riemann_solve( area , Wl , Wr , Wtemp(m,:) , 0.5*(hl+hr) )
         call velocity_deprojector( Wtemp(m,:) , metric )
      end do
      Wtempl = Wtemp(1,:) ; Wtempr = Wtemp(2,:)
      if( n .ne. 2 ) deallocate( Wtemp2 )
    
    end subroutine Wtemp_maker
    
  subroutine normal_vector(E1 , E2 , n , metric )
    implicit none
    real    , intent(in)  :: E1(4) , E2(4) 
    integer , intent(in)  :: n
    real    , intent(out) :: metric(2,2)
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
    real , intent(inout) :: W(4)
    real , intent(in)    :: metric(2,2)
    real                 :: temp(2)
    temp = (/ W(2) , W(3) /)
    
    temp = matmul(metric,temp)
    W(2) = temp(1) ; W(3) = temp(2)
    
  end subroutine velocity_projector
  
  subroutine velocity_deprojector( W , metric )
    implicit none
    real , intent(in )   :: metric(2,2)
    real , intent(inout) :: W(4)
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
  
  subroutine left_and_right_states(i,j,n,x,E,h,left_state,right_state,h0,ramp,muscl_update)
    use global_data
    use boundary_conditions_mod
    use vars
    use helper_functions
    use projection_routines
    use muscl_mod
    implicit none
    logical , intent(in)  :: ramp , muscl_update
    integer , intent(in)  :: i , j , n
    real    , intent(in)  :: E(neta,nxi,8) , h(neta,nxi) , h0 , x(neta,nxi)
    type (state) , intent(out) :: left_state , right_state
    real    ::  Er(4) , hr , Wr(4) , El(4) , hl , Wl(4) 
    real    :: pi , theta , metric(2,2)
    integer :: is , js , jshift , ishift
    real , allocatable :: W(:,:,:) , Wtemp(:,:,:) , Etemp(:,:,:) , htemp(:,:) 

    jshift = jshift_func(n)
    ishift = ishift_func(n)
    js = j - jshift
    is = i - ishift

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
       if( muscl_update .and. js-1 .gt. 0 .and. is-1 .gt. 0 &
            .and. j+jshift .lt. neta .and. i+ishift .lt. nxi &
            )then
          call MUSCL(  W(j,i,:) , W(j+jshift,i+ishift,:) &
               , W(js,is,:) , W(js-jshift,is-ishift,:) , Wl , Wr )
       else
          Wl = W(js,is, : )
          Wr = W(j ,i , : )
       end if
       El = E(js,is,5:8)
       hl = h(js,is    )
       Er = E(j ,i ,5:8)
       hr = h(j ,i     )
    end if
    deallocate( W , Wtemp , Etemp , htemp )
    left_state%phys = Wl
    left_state%geom = El
    left_state%h = hl
    right_state%phys = Wr
    right_state%geom = Er 
    right_state%h = hr
  end subroutine left_and_right_states

  end module Wtemp_maker_mod
