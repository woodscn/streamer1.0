module boundary_conditions_mod

	implicit none
	real , save , allocatable :: theta_bottom(:) , theta_top(:)
	integer , save :: boundary_switch

    contains

		function upstream_boundary_cons( ny , h )
			!Inputs: E - array of conservative variables, h - shape-preserving function,
			use vars
			implicit none
			integer              :: j , i
			integer , intent(in) :: ny
			real , intent(in) , optional :: h(:,:)
			real                 :: upstream_boundary_cons(ny,8)
			real                 :: Wtemp(ny,1,4) , Etemp(ny,1,8)
			real                 :: Wl(4) , Wr(4)

			call upstream_boundary_prim ( Wtemp)

			!Define metric components, using a cartesian grid at the inflow boundary.
			do j = 1 , ny
			Etemp(j,1,:) = (/ 0. , 0. , 0. , 0. , 1. , 0. , 0. , 1. /)
			        !if( flag .and. maxval(h) .gt. 0. )&
			             !Etemp(j,1,5) = (h(j,1)*Wtemp(j,1,2))/(h(ny,1)*Wtemp(ny,1,2))
			end do

			!Convert Wtemp (primitive variables) to Etemp (conservative variables).
			call primtocons( Etemp(:,:,1:4) , Wtemp , Etemp(:,:,5:8) )

			upstream_boundary_cons = Etemp(:,1,:)

		end function upstream_boundary_cons

		subroutine upstream_boundary_prim( Wtemp )

			implicit none
			real , intent(out)   :: Wtemp(:,:,:)
			integer              :: ny , nx , i , j , switch
			real                 :: Wl(4) , Wr(4)

			ny = size(Wtemp,1) ; nx = size(Wtemp,2)

			if( boundary_switch .eq. 1 )then

				!Basic, uniform flow
				do j = 1 , ny
					Wtemp(j,1,:) = (/ 1. , 1.8*sqrt(1.4*1.00/1.0) , 0. , 1. /)
				end do

			elseif( boundary_switch .eq. 2 )then

				!Linear velocity inflow
				do j = 1 , ny
					Wtemp(j,1,:) = (/ 1. , 1.+real(j)/real(ny) , 0. , 1. /)
				end do

			elseif( boundary_switch .eq. 3 )then
				!The Riemann problem of Hui
				Wl = (/ 1.0 , 2.4*sqrt(1.4*1.00/1.0) , 0. , 1.00 /)
				Wr = (/ 0.5 , 7.0*sqrt(1.4*0.25/0.5) , 0. , 0.25 /)

				do j = 1 , ceiling(real(ny)/2)
					Wtemp(j,1,:) = Wl
				end do
				do j = ceiling(real(ny)/2)+1 , ny
					Wtemp(j,1,:) = Wr
				end do
			end if

		end subroutine upstream_boundary_prim

	pure function wall_rotate_matrix( theta )
		implicit none
		real :: wall_rotate_matrix(4,4) , c , s
		real , intent(in) :: theta

		c = cos(2.*theta) ; s = sin(2.*theta)

		wall_rotate_matrix(1,:) = (/ 1. , 0. , 0. , 0. /)
		wall_rotate_matrix(2,:) = (/ 0. , c  , s  , 0. /)
		wall_rotate_matrix(3,:) = (/ 0. , s  ,-c  , 0. /)
		wall_rotate_matrix(4,:) = (/ 0. , 0. , 0. , 1. /)

	end function wall_rotate_matrix

	pure function wall_pressure( rho0 , u0 , p0 , shock  )
	use physical_data
	real , intent(in) :: rho0, u0 , p0
	logical , intent(in) :: shock
	real :: coeff1 , coeff2 , wall_pressure , a0
		a0 = sqrt(gamma*p0/rho0)
		if(shock)then
			coeff1 = 2./((gamma+1.)*rho0)
			coeff2 = (gamma-1.)/(gamma+1.)*p0
			wall_pressure = p0 + u0/(2.*coeff1)*(u0+sqrt(u0**2.+4.*coeff1*(coeff2+p0)))
		else
			wall_pressure = p0*(1+.5*(gamma-1)*u0/a0)**(2*gamma/(gamma-1))
		end if
	end  function wall_pressure

	subroutine wall_boundary2( dxi , Ein )
		use vars
		real , intent(in)  :: Ein(:,:,:) , dxi(:)
		real :: dE(8) , E(8) , dP , W(4) , delta , ddelta , Wtemp(4)
		real , allocatable :: Win(:,:,:)

		allocate( Win(size(Ein,1),size(Ein,2),4) )

		call constoprim( Ein(:,:,1:4) , Win(:,:,:) , Ein(:,:,5:8) )

!		write(*,*) "Ein = "
!		write(*,*) Ein(1,2,1:4)
!		write(*,*) "W = "
!		write(*,*) W  (1,2, : )

		if( size(Ein,2) == 3 )then
			dE = ( Ein(1,3,:) - Ein(1,1,:) )/( 2.*dxi(1) )
			 E =   Ein(1,2,:)
			dP = ( Win(1,3,4) - Win(1,1,4) )/( 2.*dxi(1) )
		elseif( size(Ein,2) == 2 )then
			dE = ( Ein(1,2,:) - Ein(1,1,:) )/( dxi(1) )
			 E =   Ein(1,1,:)
			dP = ( Win(1,2,4) - Win(1,1,4) )/( dxi(1) )
		end if

!		write(*,*) "dE = "
!		write(*,*) dE
!		write(*,*) "E = "
!		write(*,*) E
!		write(*,*) "dP = "
!		write(*,*) dP
!		stop

		delta  =  E(5)* E(8) -  E(6)* E(7)
		ddelta = dE(5)*dE(8) - dE(6)*dE(7)
		Wtemp(1) = W(1)
		Wtemp(2) = W(2)
		Wtemp(3) = W(3)
		Wtemp(4) =										&
			dxi(2)/( E(5)**2 + E(6)**2 )*				& ! dx/(A^2+B^2)
			( W(1)*( W(2)*E(8) - W(3)*E(7) )*			& ! rho(uM-vL)
			( W(3)*( dE(5) - E(5) * ddelta/delta ) 		& ! v(A-A'*del'/del)
			- W(2)*( dE(6) - E(6) * ddelta/delta )	)	& ! u(B-B'*del'/del)
			+ ( E(5)*E(7) + E(6)*E(8) ) * dP )			  ! (AL+BM)P'

		deallocate( Win )

	end subroutine

	subroutine wall_boundary( Wl , Wtemp )
		real , intent(in) :: Wl(:)
		real , intent(inout):: Wtemp(:)
		logical :: shock
		shock = Wtemp(3)>0.
		Wtemp(3) = 0.
		Wtemp(4) = wall_pressure( Wtemp(1) , Wtemp(3) , Wtemp(4) , shock )
!		write(*,*) "Got here"
!		stop
	end  subroutine wall_boundary

end module boundary_conditions_mod
