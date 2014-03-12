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
    integer :: neta , nxi
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

!  subroutine constoprimvec( E , W , geom )
!! Inputs: E(1,1,8) - the array of conservative variables. The first 4 "rows"
!!   in the third index are conservative flow variables, while the last 4 "rows" 
!!   are geometric variables.!

!! Outputs: W(1,1,4) - the array of primitive flow variables rho, u, v, p.
!    use physical_data
!    implicit none   
!    real , intent(in)  :: E(4) , geom(4)
!    real , intent(out) :: W(4)
!    real               :: delta
!    delta = geom(1)*geom(4) - geom(2)*geom(3) ! AM-BL
!! Compute primitive variables from conservative variables
!    W(1) = E(1)/delta
!    W(2) = E(2)/E(1)
!    W(3) = E(3)/E(1)
!    W(4) = (gamma-1.)*(E(4)/delta-0.5*W(1)*(W(2)**2+W(3)**2))
!  end subroutine constoprimvec

end module vars

