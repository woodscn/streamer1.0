module flow_update_mod
contains
  subroutine flow_update( E , h &
       , x , y , qold , ramp , muscl_update , h0 )
    use vars
!    use window_data
    use riemann
    use muscl_mod
!	use upstream_boundary_mod
    use global_data
    use physical_fluxes
    use boundary_conditions_mod
    implicit none
    real , intent(in)    :: h(:,:) , x(:,:) , y(:,:) , h0
    real , intent(inout) :: E(:,:,:) 
    real , intent(out)   :: qold(:,:,:)
    logical , intent(in) :: ramp , muscl_update
    integer :: n , i , j , jshift , ishift , jp , ip , js , is , nx , ny , side
    real , allocatable   :: W(:,:,:) , Wtemp(:,:,:) , Wtemp2(:,:,:) &
         , Etemp(:,:,:) , htemp(:,:) , rho_inv(:,:)
    real :: area , dlambda , Wl(4) , Wr(4) , El(4) , Er(4)  , hl , hr , metric(4,4) , wall_rotate_matrix
    real :: Fl(4) , Fr(4) , vol_inv  , theta , pi , Wtempl(4) , Wtempr(4)
    real :: geom1(neta,nxi) , geom2(neta,nxi) , diss
    real , allocatable   :: momx(:,:) , momy(:,:) , kexx(:,:) , kexy(:,:) , keyy(:,:)
    real , allocatable   :: dp(:,:) , rue(:,:) , rve(:,:)
    real , allocatable   :: workx(:,:), worky(:,:) , energy(:,:)
    logical :: tick

    vol_inv = 1./(deta*dxi)
    allocate( W(neta,nxi,4) , rho_inv(neta,nxi) ) ; W = 0.
    qold(:,:,1) = E(:,:,2)/E(:,:,1)
    qold(:,:,2) = E(:,:,3)/E(:,:,1)

    do n = 1 , 3
       call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
       if( n .ne. 2 )then;
          ishift = 1; jshift = 0; area = deta; dlambda = dt*.5
       else
          ishift = 0; jshift = 1; area = dxi ; dlambda = dt
       end if
       allocate( Wtemp(neta+jshift,nxi+ishift,4) )
       if( n  .ne. 2 )then
          allocate( Wtemp2(neta,1,4) , Etemp(neta,1,8) , htemp(neta,1) )
          call upstream_boundary_prim( Wtemp2 )
          Etemp(:,1,:) = upstream_boundary_cons( size(Etemp,1) )
          htemp(:,1) = h(:,1)
       end if
       do i = 1 , nxi +ishift
          do j = 1 , neta +jshift
             js = j - jshift
             is = i - ishift
             if( is .eq. 0 )then
                Er = E(j ,i ,5:8)
                hr = h(j ,i     )
                Wr = W(j ,i , : )
                El = Etemp(j,i,5:8)
                hl = hr
                Wl = Wtemp2(j,i, : )
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
                if( x(j,i) .ge. 0.5 .and. x(j,i) .lt. 1.0 .and. ramp )then
                   pi = 3.141592653589793
                   theta =  pi/180.*15.
                   metric = wall_rotate_matrix( theta )
					Wl = matmul(metric,Wl)
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
                if( muscl_update .and. js-jshift .gt. 0 .and. is-ishift .gt. 0 &
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

			 metric = riemann_rotate_matrix(El,Er,n)

             Wl = matmul(metric,Wl)
             Wr = matmul(metric,Wr)

			 call riemann_solve( Wl , Wr , Wtemp(j,i,:) , 0.5*(hl+hr) )

             Wtemp(j,i,:) = matmul(transpose(metric),Wtemp(j,i,:))

!			 if( j .eq. neta + 1 ) then
!			 	call wall_boundary( Wl , Wtemp(j,i,:) )
!
!		 	write(*,*) Wtemp(j,i,4)
!
!		 	end if
!
!			 if( j .eq. neta + 1 .and. is .gt. 1 .and. is .lt. nxi ) then
!				call wall_boundary2( (/ dxi , deta /) , E(js:js,is-1:is+1,:) )
!
!				write(*,*) Wtemp(j,i,4)
!			end if
          end do
       end do

       if( n .ne. 2 )then
          E(:,:,5) = E(:,:,5) + dlambda/dxi *h(:,:)*( Wtemp(:,2:nxi +1,2)-Wtemp(:,1:nxi ,2) )
          E(:,:,6) = E(:,:,6) + dlambda/dxi *h(:,:)*( Wtemp(:,2:nxi +1,3)-Wtemp(:,1:nxi ,3) )
       else
          E(:,:,7) = E(:,:,7) + dlambda/deta*h(:,:)*( Wtemp(2:neta+1,:,2)-Wtemp(1:neta,:,2) )
          E(:,:,8) = E(:,:,8) + dlambda/deta*h(:,:)*( Wtemp(2:neta+1,:,3)-Wtemp(1:neta,:,3) )
       end if
       
       E(:,:,1:4) = E(:,:,1:4)-dlambda*area*vol_inv*fluxes(Wtemp,E,h,ishift,jshift,nxi,neta,n)

       deallocate( Wtemp )
       if( n .ne. 2 ) deallocate( Wtemp2 , Etemp , htemp )
    end do
    deallocate( W , rho_inv )
    
  end subroutine flow_update

end module flow_update_mod
