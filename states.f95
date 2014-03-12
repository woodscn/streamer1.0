module states
  implicit none
  type :: state
     real :: phys(4)
     real :: geom(4)
     real    :: h
  end type state
contains
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
end module states
