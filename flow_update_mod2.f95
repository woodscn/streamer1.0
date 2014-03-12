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
    integer :: n , i , j , jshift , ishift , jp , ip , js , is , nx , ny
    real , allocatable   :: W(:,:,:) , Wtemp(:,:,:) , Wtemp2(:,:,:) &
         , Etemp(:,:,:) , htemp(:,:) , rho_inv(:,:) ,F(:,:,:)
    real :: area , dlambda , Wl(4) , Wr(4) , El(4) , Er(4)  , hl , hr , metric(4,4) !, wall_rotate_matrix
    real :: Fl(4) , Fr(4) , vol_inv  , theta , pi , Wtempl(4) , Wtempr(4)
    real :: geom1(neta,nxi) , geom2(neta,nxi) , diss
    real , allocatable   :: momx(:,:) , momy(:,:) , kexx(:,:) , kexy(:,:) , keyy(:,:)
    real , allocatable   :: dp(:,:) , rue(:,:) , rve(:,:)
    real , allocatable   :: workx(:,:), worky(:,:) , energy(:,:)
    logical :: tick

    vol_inv = 1./(deta*dxi)
    allocate( W(neta,nxi,4) , rho_inv(neta,nxi) ,F(neta,nxi,4) ) ; W = 0.
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
                   theta =  pi/12.
                   !                   metric = wall_rotate_matrix( theta )
                   Wl = matmul(wall_rotate_matrix( theta ),Wl)
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
                Wr = matmul(wall_rotate_matrix( 0. ),Wl)
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
             if( js .ge. 1 .and. is .ge. 1 )then
                if( n .ne. 2 )then
                   E(js,is,5) = E(js,is,5) + dlambda/dxi *h(js,is)*( Wtemp(j,i,2)-Wtemp(js,is,2) )
                   E(js,is,6) = E(js,is,6) + dlambda/dxi *h(js,is)*( Wtemp(j,i,3)-Wtemp(js,is,3) )
                else
                   E(js,is,7) = E(js,is,7) + dlambda/deta*h(js,is)*( Wtemp(j,i,2)-Wtemp(js,is,2) )
                   E(js,is,8) = E(js,is,8) + dlambda/deta*h(js,is)*( Wtemp(j,i,3)-Wtemp(js,is,3) )
                end if
                !			 F = fluxes(Wtemp,E,h,ishift,jshift,nxi,neta,n)
                E(js,is,1:4) = E(js,is,1:4)-dlambda*area*vol_inv&
                     *elemental_fluxes(Wtemp(js,is,:),Wtemp(j,i,:),E(js,is,:),h(js,is),n)

                if( js == 51)then
                   call constoprim(E(:,:,1:4),W,E(:,:,5:8))
                   write(*,*) W(js-1,is,1),W(js,is,1),W(js+1,is,1)
                   write(*,*) "Elemental_fluxes inputs:"
                   write(*,*) "Left interface:"
                   write(*,*) Wtemp(js,is,1:4)
                   write(*,*) "Right interface:"
                   write(*,*) Wtemp(j,i,1:4)
                   write(*,*) "Geometry:"
!                   write(*,*) W(js,is,:)
                   write(*,*) [1.0101010101010101019e-2,&
                        0.,0.,1.0101010101010101019e-2]
                   write(*,*) "Grid motion = "
                   write(*,*) h(js,is)
                   write(*,*) h(js,is)*W(js,is,2:3)
                   write(*,*) "Elemental fluxes debuggery:"
                   write(*,*) "Total change in conservative variables: "
                   write(*,*) "( Mass, Momentum(x), Momentum(y), Momentum(z), Energy )"
                   write(*,*) -1*elemental_fluxes(Wtemp(js,is,:),Wtemp(j,i,:),&
                        [0.,0.,0.,0.,1.0101010101010101019e-2,&
                        0.,0.,1.0101010101010101019e-2],h(js,is),n)
                   read(*,*)
                end if
             end if
             !			 F =
          end do
       end do

       !       if( n .ne. 2 )then
       !          E(:,:,5) = E(:,:,5) + dlambda/dxi *h(:,:)*( Wtemp(:,2:nxi +1,2)-Wtemp(:,1:nxi ,2) )
       !          E(:,:,6) = E(:,:,6) + dlambda/dxi *h(:,:)*( Wtemp(:,2:nxi +1,3)-Wtemp(:,1:nxi ,3) )
       !       else
       !          E(:,:,7) = E(:,:,7) + dlambda/deta*h(:,:)*( Wtemp(2:neta+1,:,2)-Wtemp(1:neta,:,2) )
       !          E(:,:,8) = E(:,:,8) + dlambda/deta*h(:,:)*( Wtemp(2:neta+1,:,3)-Wtemp(1:neta,:,3) )
       !       end if

       !       E(:,:,1:4) = E(:,:,1:4)-dlambda*area*vol_inv*fluxes(Wtemp,E,h,ishift,jshift,nxi,neta,n)

       deallocate( Wtemp )
       if( n .ne. 2 ) deallocate( Wtemp2 , Etemp , htemp )
    end do
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
100 format (f10.8)
    do j = 1, size(W,2)
       write(*,100) W(:,j,4)
    end do
    read(*,*)
    deallocate( W , rho_inv )

  end subroutine flow_update
end module flow_update_mod
