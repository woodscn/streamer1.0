module projection_routines
contains
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
  
end module projection_routines
