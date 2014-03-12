module h_update_mod
! This module is responsible for solving the PDE for h. It uses an iterative, 
!   first-order, backward difference scheme. Boundary values must be provided 
!   at h(1,:) and h(:,1).
contains

  subroutine h_update( A , B , L , M , h , u , v &!, xi , eta 
       , h0 )
    use global_data
    implicit none
    real , intent(in)   :: A(:,:) , B(:,:) ,  L(:,:) ,   M(:,:)
    real , intent(in)   :: u(:,:) , v(:,:) &!, xi(:,:) , eta(:,:) 
         , h0
    real , intent(inout):: h(:,:)
    real                :: tol  , err !, dxi , deta
    real , allocatable  :: dthetadxi(:,:) , dthetadeta(:,:) , temp(:,:)
    real , allocatable  :: alpha(:,:) , beta(:,:) , gamma(:,:) , g(:,:)
    real , allocatable  :: S2(:,:) , T2(:,:) , q(:,:) , theta(:,:) 
    integer             :: ny , nx , i , j , n , k 

    ny = size(h,1) ; nx = size(h,2) ! Define these size parameters
!    dxi = xi(2,2)-xi(1,1) ; deta = eta(2,1)-eta(1,1) 
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
    gamma = S2*(B*sin(theta)+A*cos(theta))*dthetadxi &
          - T2*(L*cos(theta)+M*sin(theta))*dthetadeta

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
       temp(:,1) = temp(:,2) ; temp(1,:) = temp(2,:)
       err = maxval(abs(temp-g))
       g   = temp
    end do
! Convert the updated g back to h.
    h = exp(g)/q
    h = h/maxval(h)*h0

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

