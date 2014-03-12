module riemann
!
!
!program riemann
!implicit none
!real :: WinL(4) , WinR(4) , Wout(4) , h
!
!WinL = (/ 1.0     ,  0.0     , 0.0 ,    1.0    /)
!WinR = (/ 0.125   ,  0.0     , 0.0 ,    0.1    /)

!WinL = (/ 1.0     , -2.0     , 0.0 ,    0.4    /)
!WinR = (/ 1.0     ,  2.0     , 0.0 ,    0.4    /)

!WinL = (/ 1.0     ,  0.0     , 0.0 , 1000.0    /)
!WinR = (/ 1.0     ,  0.0     , 0.0 ,    0.01   /)

!WinL = (/ 1.0     ,  0.0     , 0.0 ,    0.01   /)
!WinR = (/ 1.0     ,  0.0     , 0.0 ,  100.0    /)

!WinL = (/ 5.99924 , 19.5975  , 0.0 ,  460.894  /)
!WinR = (/ 5.99242 , -6.19633 , 0.0 ,   46.0950 /)
!
!h = 1.0
!call riemann_solve( WinL , WinR , Wout , h )


contains

	pure function riemann_rotate_matrix( EINL , EINR , n )
		implicit none
		real :: riemann_rotate_matrix(4,4)
		real , intent(in) :: EINL(4) , EINR(4)
		integer , intent(in) :: n
 		real :: R(4,4)
		real :: geom1, geom2
		real :: S_inv , deltaL_inv , deltaR_inv

		deltaL_inv = 1./(EINL(1)*EINL(4)-EINL(2)*EINL(3))
		deltaR_inv = 1./(EINR(1)*EINR(4)-EINR(2)*EINR(3))

		geom1 = (mod(n,2)*EINL(3) + (mod(n,2)-1)*EINL(1))*deltaL_inv + &
				(mod(n,2)*EINR(3) + (mod(n,2)-1)*EINR(1))*deltaR_inv
		geom2 = (mod(n,2)*EINL(4) + (mod(n,2)-1)*EINL(2))*deltaL_inv + &
				(mod(n,2)*EINR(4) + (mod(n,2)-1)*EINR(2))*deltaR_inv

		S_inv = 1./sqrt(geom1**2+geom2**2)

		R(1,:) = (/ 1. ,           0. ,            0. , 0. /)
		R(2,:) = (/ 0. ,  geom2*S_inv ,  -geom1*S_inv , 0. /)
		R(3,:) = (/ 0. ,  geom1*S_inv ,   geom2*S_inv , 0. /)
		R(4,:) = (/ 0. ,           0. ,            0. , 1. /)

		riemann_rotate_matrix = R

	end function riemann_rotate_matrix

  subroutine riemann_solve( WinL , WinR , Wout , h )
    implicit none

    real , intent(in)  :: WinL(:) , WinR(:) , h
    real , intent(out) :: Wout(:)
    real :: DL , UL , VL , PL , DR , UR , VR , PR , AL , AR
    real :: Pstar , Ustar , DstarL , DstarR , cL , cR , betaL , betaR
    real :: PsiL , PsiR , UstarL , UstarR , dcdpsi , dbetadpsi 
    real :: dUstarL , dUstarR , temp = 1. , cTL , cHL , cTR , cHR
    real :: Dout , Uout , Vout , Pout , tol = .000001 , gL , gR , x
    real :: dfL , dfR , fL , fR
    real , parameter :: gamma = 1.4
    integer :: n = 1 
    integer , parameter :: nx = 1000
    real :: data(nx,5)
    logical :: verbose = .false.

    DL = WinL(1) ; UL = WinL(2) ; VL = WinL(3) ; PL = WinL(4)
    DR = WinR(1) ; UR = WinR(2) ; VR = WinR(3) ; PR = WinR(4)
    AL = sqrt(gamma*PL/DL) ; AR = sqrt(gamma*PR/DR)

! Linearised guess
    Pstar = .5*(PL+PR)-.125*(uR-uL)*(DL+DR)*(aL+aR)
    if(.not.( Pstar .gt. min(PL,PR) .and. Pstar .lt. max(PL,PR) &
         .and. max(PL,PR)/min(PL,PR) .le. 2.0))then
       if(Pstar .lt. min(PL,PR))then
! Two-rarefaction solution
          Pstar = (&
               (aL+aR-.5*(gamma-1.)*(uR-uL))&
               /(aL/pL**((gamma-1.)/(2.*gamma))+aR/pR**((gamma-1.)/(2.*gamma)))&
               )**(2.*gamma/(gamma-1.))
       else
! Two-shock solution
          gL=sqrt( 2./((gamma+1.)*DL)/(Pstar+PL*(gamma-1.)/(gamma+1.)))
          gR=sqrt( 2./((gamma+1.)*DR)/(Pstar+PR*(gamma-1.)/(gamma+1.)))
          Pstar = max(tol,&
               (gL*PL+gR*PR-(UR-UL))/(gL+gR)&
               )
       end if
    end if
!          Pstar = .5*( PL + PR )
	if(verbose)write(*,*) "Initial guess P = " , Pstar
    do while(abs(temp) .gt. tol)
       if(verbose)write(*,*) n , "Pstar = " , Pstar; n = n + 1
!       if(n .gt. 15)then ;  write(*,*) "Failed Convergence" ; stop ; end if

       PsiL = Pstar/PL
       call u_fun(UL,PL,aL,PsiL,gamma,h,fL,dfL)

       PsiR = Pstar/PR
       call u_fun(UR,PR,aR,PsiR,gamma,h,fR,dfR)

       temp = ( UR - UL + fR + fL )/( dfL + dfR )
       Pstar = max( Pstar - temp , tol )
       if(verbose)write(*,*) "fL , fR , Pstar = " , fL , fR , Pstar
    end do

    Ustar = .5*(UR+fR+UL-fL)
    DstarL = beta(PsiL,gamma)*DL
    DstarR = beta(PsiR,gamma)*DR
    if(verbose)write(*,*) Ustar , DstarL , DstarR
    open(unit = 1, file = 'riemann.dat')
    do n = 1 , nx
       x = 10./(real(nx+1))*real(n)-5.
       call sample(x,WinL,WinR,(/Pstar,Ustar,DstarL,DstarR/),h,gamma,Wout)
       data(n,1)   = x
       data(n,2:5) = Wout
       write(1,*) data(n,:)
    end do
!    write(*,*) "Winstar = "
!    write(*,*) Pstar , Ustar , DstarL , DstarR
!    Wout = (/ Dout , Uout , Vout , Pout /)

  end subroutine riemann_solve

  subroutine u_fun(u0,p0,a,psi,gamma,h,f,df)
    real , intent(in)  :: u0 , p0 , a , psi , gamma , h
    real , intent(out) :: f , df
    real :: gamma_hat , gamma_coeff , gamma_fac , xi


    gamma_hat   =  gamma - 1.
    gamma_coeff =  gamma_hat + 2.*eta

    if( psi .gt. 1. )then
       f = 2.*a/(gamma-1.)*(psi**((gamma-1.)/(2.*gamma))-1.)
       df= a/(gamma*p0)*psi**(-(gamma+1.)/(2.*gamma))
    else
       f = a*(psi-1.)/sqrt(.5*gamma*(psi*(gamma+1.)+gamma-1.))
       df= gamma*(-1.+psi+3.*gamma+a*gamma)/&
       	(sqrt(2.)*(gamma*(-1+a+gamma+a*gamma))**1.5)*a/p0
    end if

  end subroutine u_fun


  function beta(psi,gamma)
    implicit none
    real :: psi, beta, gamma
    if( psi .gt. 1 )then
       beta = ( (gamma+1.)*psi + gamma-1. )&
             /( gamma+1. + (gamma-1.)*psi )
    else
       beta = psi**(1./gamma)
    end if
  end function beta

  subroutine sample(x,WinL,WinR,Winstar,h,gamma,Wout)
    implicit none
    real , intent(in)  :: WinL(4) , WinR(4) , Winstar(4) &
         , h , gamma 
    real , intent(out) :: Wout(4)

    real :: DL , UL , VL , PL , DR , UR , VR , PR , Pstar &
         , Ustar , DstarL , DstarR , eta , gamma_hat , PsiL , xi
    real :: PsiR , aL , aR , betaL , betaR , Dout , Uout &
         , Vout , Pout , cL , cR , cLT , cLH , cRT , cRH , term , x
    DL = WinL(1) ; UL = WinL(2) ; VL = WinL(3) ; PL = WinL(4)
    DR = WinR(1) ; UR = WinR(2) ; VR = WinR(3) ; PR = WinR(4)
    Pstar  = Winstar(1) ; Ustar  = Winstar(2) 
    DstarL = Winstar(3) ; DstarR = Winstar(4)

    PsiL = Pstar/PL
    PsiR = Pstar/PR
    aL   = sqrt(gamma*PL/DL)
    aR   = sqrt(gamma*PR/DR)
    betaL= DstarL/DL
    betaR= DstarR/DR

    if( Ustar .gt. x )then
!       write(*,*) "The boundary lies to the left of the contact wave"

       Vout = VL
       if( PsiL .gt. 1. )then
!          write(*,*) " Left shock"
          cL = (1.-h)*UL-aL*sqrt((gamma+1.)/(2.*gamma)*(PsiL-1.)+1.)
!          write(*,*) "Left shock speed = " , cL
          if( cL .gt. x )then
!             write(*,*) " The boundary lies to the left of the shock"
             Pout = PL
             Uout = UL
             Dout = DL
          else
!             write(*,*) " The boundary lies in the left central region"
             Pout = Pstar
             Uout = Ustar
             Dout = DstarL
          end if
       else
!          write(*,*) "Left rarefaction wave"
          cLT = (1.-h)*UL - aL
!          write(*,*) "Left rarefaction tail speed = " , cLT
          cLH = (1.-h)*Ustar - sqrt(gamma*Pstar/DstarL)
!          write(*,*) "Left rarefaction head speed = " , cLH
          if( cLT .gt. x )then
!             write(*,*) " The boundary lies to the left of the wave"
             Pout = PL
             Uout = UL
             Dout = DL

          elseif( cLH .lt. x )then
!             write(*,*) "The boundary lies in the left central region"
             Pout = Pstar
             Uout = Ustar
             Dout = DstarL
          else
!             write(*,*) "The boundary lies within the left expansion wave"
             Pout = PL*(2.*(1.-h)/(gamma-2.*h+1.)+(gamma-1.)/(aL*(gamma-2.*h+1.))&
             	*(1.-h)*UL)**((2.*gamma)/(gamma-1.))
             Uout = UL - 2.*aL/(gamma-1.)*((Pout/PL)**((gamma-1.)/(2.*gamma))-1.)
             Dout = DL*(Pout/PL)**(1./gamma)
          end if
       end if
    else
!       write(*,*) " The boundary lies to the right of the contact wave"
       Vout = VR
       if( PsiR .gt. 1. )then
!          write(*,*) " Right shock"
          cR = (1.-h)*UR+aR*sqrt((gamma+1.)/(2.*gamma)*(PsiR-1.)+1.)
!          write(*,*) " Right shock speed = " , cR
          if( cR .lt. x )then
!             write(*,*) " The boundary lies to the right of the shock"
             Pout = PR
             Uout = UR
             Dout = DR
          else
!             write(*,*) " The boundary lies in the right central region"
             Pout = Pstar
             Uout = Ustar
             Dout = DstarR
          end if
       else
!          write(*,*) " Right rarefaction wave"
          cRT = (1.-h)*UR + aR
!          write(*,*) "Right rarefaction tail speed = " , cRT
          cRH = (1.-h)*Ustar + sqrt(gamma*Pstar/DstarR)
!          write(*,*) "Right rarefaction head speed = " , cRH
          if( cRT .lt. x )then
!             write(*,*) " The boundary lies to the right of the wave"
             Pout = PR
             Uout = UR
             Dout = DR

          elseif( cRH .gt. x )then
!             write(*,*) "The boundary lies in the right central region"
             Pout = Pstar
             Uout = Ustar
             Dout = DstarR
          else
!             write(*,*) " The boundary lies within the right expansion wave"
             Pout = PR*(2.*(1.-h)/(gamma-2.*h+1.)-(gamma-1.)/(aR*(gamma-2.*h+1.))&
             	*(1.-h)*UR)**((2.*gamma)/(gamma-1.))
             Uout = UR + 2.*aR/(gamma-1.)*((Pout/PR)**((gamma-1.)/(2.*gamma))-1.)
             Dout = DR*(Pout/PR)**(1./gamma)
          end if
       end if
    end if

    Wout = (/ Dout , Uout , Vout , Pout /)

  end subroutine sample

!end program riemann
end module riemann
