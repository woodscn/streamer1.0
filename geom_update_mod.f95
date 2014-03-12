module geom_update_mod
! This module exists to advance the x,y coordinates in time, to create/destroy points in order to keep the computation within the computational window, and to solve for updated values of geometric variables.
contains
  subroutine geom_update(geom_switch)
    logical :: geom_switch
! geom_update acts only as the master routine to direct the action of the three individual subroutines.
    call coords_update ! Move the x,y coordinates with the fluid velocity
    call windower(geom_switch) ! Create/destroy points to maintain the computational window
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

  subroutine windower(geom_switch)
! Windower is tasked with keeping the solution within a designated computational window. It evaluates whether the grid has overrun the maximum and minimum required values of x, and creates/destroys columns as necessary. 
    use global_data
    use primary_variables
    use window_data
    use boundary_conditions_mod

    implicit none
    logical :: overrun = .false. , shortfall = .false. , geom_switch
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

    if( geom_switch .and. nxi .gt. 10 )then
       overrun = .true.
       geom_switch = .false.
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
! Store the current values in temporary arrays
       call allocate_temp
       x_t = x ; y_t = y ; E_t = E ; h_t = h
! Overwrite all but the first column with data from the temporary arrays
       x(:,2:nxi)   = x_t(:,1:(nxi-1))   ; y(:,2:nxi)   = y_t  (:,1:(nxi-1))
       E(:,2:nxi,:) = E_t(:,1:(nxi-1),:) ; h(:,2:nxi)   = h_t  (:,1:(nxi-1))
       call deallocate_temp
       call upwind_boundary()
    endif

    contains
      subroutine upwind_boundary()
! Call boundary_conditions to compute the first column for E,h.
        E(:,1,:) = upstream_boundary_cons( size(E,1) )
        h(:,1) = h(:,2)
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

