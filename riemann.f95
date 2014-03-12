module riemann
! This module contains the riemann solver used in my code. It is adapted from 
! the one presented by Toro in his book on Riemann problems in fluid mechanics.
! It is an exact, iterative solver, and it may eventually be beneficial to 
! use a faster, apporoximate solver instead. 

! It also contains a function that creates a matrix to transform a vector of
! primitive variables to a coordinate frame normal to the cell interfaces.
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

  subroutine riemann_solve(WINL, WINR, WOUT, HAVG)

! Inputs: Domlen - the "length" of the domain. Corresponds to either dxi or deta.
!         WinL   - the left state in terms of primitive variables
!         WinR   - the right state in terms of primitive variables
!         Wout   - the output state, used to compute fluxes, in primitive variables
!         Smax   - the speed of the fastest nonlinear wave in the problem. Useful 
!                     for computing time steps as part of a CFL condition.
    IMPLICIT NONE

!     Declaration of variables:

    INTEGER :: I, CELLS = 40
    REAL    :: DIAPH = 0., DOMLEN=1., DS = 0., DX, PM, MPA = 1., PS = 0., S = 0.,&
         &        TIMEOUT = 0.035, UM, US = 0., XPOS, SMAX, WOUT(4), WINL(4), WINR(4)
    REAL    :: GAMMA, G1, G2, G3, G4, G5, G6, G7, G8, &
         &           DL, UL, PL, CL, DR, UR, PR, CR, SL, SR , VL , VR , VS , &
         &           HAVG
    COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
    COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR 

    GAMMA = 1.4
    DIAPH = DOMLEN*0.5
    DL = WINL(1) ; UL = WINL(2) ; VL = WINL(3) ; PL = WINL(4)
    DR = WINR(1) ; UR = WINR(2) ; VR = WINR(3) ; PR = WINR(4)

!     Compute gamma related constants

    G1 = (GAMMA - 1.0)/(2.0*GAMMA)
    G2 = (GAMMA + 1.0)/(2.0*GAMMA)
    G3 = 2.0*GAMMA/(GAMMA - 1.0)
    G4 = 2.0/(GAMMA - 1.0)
    G5 = 2.0/(GAMMA + 1.0)
    G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
    G7 = (GAMMA - 1.0)/2.0
    G8 = GAMMA - 1.0

!   Compute sound speeds
    CL = SQRT(GAMMA*PL/DL)
    CR = SQRT(GAMMA*PR/DR)
    
!   The pressure positivity condition is tested for
    IF(G4*(CL+CR).LE.(UR-UL))THEN
!   The initial data is such that a vacuum is generated
       WRITE(6,*)
       WRITE(6,*)'***Vacuum is generated by data***' 
       WRITE(6,*)'***Program stopped***'
       WRITE(6,*)
       STOP
    ENDIF

!   Exact solution for pressure and velocity in star region is found

    CALL STARPU(PM, UM, MPA)
    
    DX = DOMLEN/REAL(CELLS)

!      Solution at point (X,T) = ( XPOS - DIAPH,TIMEOUT ) is found
         
    S = 0.
    CALL SAMPLE(PM, UM, S, DS, US, PS, VS, HAVG)

    WOUT = (/ DS, US, VS, PS /)
	 
    if(pm.le.pl)then
       SL = ABS( UL - CL )
    else 
       SL = ABS( UL - CL/(2.*GAMMA)*((GAMMA-1.)*PM/PR&
            &               +(GAMMA-1.))**.5 )
    end if

    if(pm.gt.pr)then
       SR = ABS( UR + CR )
    else 
       SR = ABS( UR + CR/(2.*GAMMA)*((GAMMA-1.)*PM/PR&
            &               +(GAMMA-1.))**.5 )
    end if
    smax = max(smax, max(SL,SR))
    !write(*,*) 'smax =  ',smax

  END SUBROUTINE riemann_solve


  SUBROUTINE STARPU(P, U, MPA)
    IMPLICIT NONE

!   Purpose: to compute the solution for pressure and velocity 
!             in the Star Region

!   Declaration of variables

    INTEGER :: I, NRITER
    
    REAL    :: DL, UL, PL, CL, DR, UR, PR, CR,&
         &        CHANGE, FL, FLD, FR, FRD, P, POLD, PSTART,&
         &        TOLPRE, U, UDIFF, MPA, VL, VR

    COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR
    DATA TOLPRE, NRITER/1.0E-14, 10/

!   Guessed value PSTART is computed

    CALL GUESSP(PSTART)

    POLD  = PSTART
    UDIFF = UR - UL

!    WRITE(6,*)'-----------------------------------------'
!    WRITE(6,*)'   Iteration number        Change   '
!    WRITE(6,*)'-----------------------------------------'

    DO 10 I = 1, NRITER
       CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
       CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
       P      = POLD - (FL + FR + UDIFF)/(FLD + FRD)
       CHANGE = 2.0*ABS((P - POLD)/(P + POLD))
       IF(CHANGE.LE.TOLPRE)GOTO 20
       IF(P.LT.0.0)P = TOLPRE
       POLD = P
10     CONTINUE
! I find that this error usually indicates that the solution has gone
! unstable at the flow level. The N-R iteration routine itself works 
! rather well.
       WRITE(6,*)'Divergence in Newton-Raphson iteration'
       STOP
20     CONTINUE

!   Compute velocity in Star Region
       U = 0.5*(UL + UR + FR - FL)

!            WRITE(6,*)'-----------------------------------------'
!            WRITE(6,*)'   Pressure             Velocity'
!            WRITE(6,*)'-----------------------------------------'
!            WRITE(6, 40)P/MPA, U
!            WRITE(6,*)'-----------------------------------------'

30     FORMAT(5X, I5,15X, F12.7)
40     FORMAT(2(F14.6, 5X))

       RETURN
     END SUBROUTINE STARPU


     SUBROUTINE GUESSP(PM)

!   Purpose: to provide a guess value for pressure 
!            PM in the Star Region. The choice is made
!            according to adaptive Riemann solver using
!            the PVRS, TRRS and TSRS approximate
!            Riemann solvers. See Sect. 9.5 of Capht. 9
!            of Ref. 1

       IMPLICIT NONE

!   Declaration of variables

       REAL    :: DL, UL, PL, CL, DR, UR, PR, CR,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,&
            &        CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,&
            &        PTL, PTR, QMAX, QUSER, UM, VL , VR

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR

       QUSER = 2.0

!   Compute guess pressure from PVRS Riemann solver

       CUP  = 0.25*(DL + DR)*(CL + CR)
       PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
       PPV  = AMAX1(0.0, PPV)
       PMIN = AMIN1(PL,  PR)
       PMAX = AMAX1(PL,  PR)
       QMAX = PMAX/PMIN

       IF(QMAX.LE.QUSER.AND.&
            & (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN

!      Selcect PVRS Riemann solver
          PM = PPV
       ELSE
          IF(PPV.LT.PMIN)THEN

!         Select Two-Rarefaction Riemann solver
             PQ  = (PL/PR)**G1
             UM  = (PQ*UL/CL + UR/CR + &
                  &            G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
             PTL = 1.0 + G7*(UL - UM)/CL
             PTR = 1.0 + G7*(UM - UR)/CR
             PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
          ELSE

!      Select Two-Shock Riemann solver with PVRS as estimate

             GEL = SQRT((G5/DL)/(G6*PL + PPV))
             GER = SQRT((G5/DR)/(G6*PR + PPV))
             PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL+GER)
          ENDIF
       ENDIF
       RETURN
     END SUBROUTINE GUESSP


     SUBROUTINE PREFUN(F, FD, P, DK, PK, CK)
!   Purpose: to evaluate the pressure fucntions
!            FL and FR in exact Riemann solver

       IMPLICIT NONE

!   Declaration of variables

       REAL    :: AK, BK, CK, DK, F, FD, P, PK, PRAT, QRT,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       IF(P.LE.PK)THEN

!      Rarefaction wave
          PRAT = P/PK
          F    = G4*CK*(PRAT**G1 - 1.0)
          FD   = (1.0/(DK*CK))*PRAT**(-G2)
       ELSE

!      Shock wave
          AK  = G5/DK
          BK  = G6*PK
          QRT = SQRT(AK/(BK + P))
          F   = (P - PK)*QRT
          FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
       ENDIF
       RETURN
     END SUBROUTINE PREFUN



     SUBROUTINE SAMPLE(PM, UM, S, D, U, P, V, HAVG)

!   Purpose: to sample the solution throughout the wave 
!            pattern. Pressure PM and velocity UM in the 
!            Star Region are known. Sampling is performed
!            in terms of the 'speed' S = X/T. Sampled
!            values are D, U, P

!   Input variables : PM, UM, S, /GAMMAS/, /STATES/
!   Output variables: D, U, P

       IMPLICIT NONE

!   Declaration of variables

       LOGICAL :: tick = .true.
       REAL    DL, UL, PL, CL, DR, UR, PR, CR,&
            &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,&
            &        C, CML, CMR, D, P, PM, PML, PMR, S,&
            &        SHL, SHR, SL, SR, STL, STR, U, UM, VL, VR, V,&
            &        HAVG

       COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
       COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR , VL , VR
      
       IF(S.LE.UM)THEN

!      Sampling point lies to the left of the contact
!      discontinuity
          V = VL

          IF(PM.LE.PL)THEN

!         Left rarefaction

             IF(tick)THEN
                SHL = (1.-HAVG)*UL - CL      
             ELSE
                SHL = UL - CL
             END IF

             IF(S.LE.SHL)THEN

!            Sampled point is left data state

                D = DL
                U = UL
                P = PL

             ELSE
                CML = CL*(PM/PL)**G1
                IF(tick)THEN
                   STL = (1.-HAVG)*UM - CML
                ELSE
                   STL = UM - CML
                END IF

                IF(S.GT.STL)THEN

!               Sampled point is Star Left state

                   D = DL*(PM/PL)**(1.0/GAMMA)
                   U = UM
                   P = PM
                ELSE

!               Sampled point is inside left fan

                   U = G5*(CL + G7*UL + S)
                   C = G5*(CL + G7*(UL - S))
                   D = DL*(C/CL)**G4
                   P = PL*(C/CL)**G3

                   P = PL*(2.*(1.-HAVG)/(GAMMA-2.*HAVG+1.) + (GAMMA-1.)&
                        /((GAMMA-2.*HAVG+1.)*CL)*((1.-HAVG)*UL))**G3
                   D = DL*(P/PL)**(1./GAMMA)
                   U = UL - 2.*CL/(GAMMA-1.)*((P/PL)**((GAMMA-1.)/(2.*GAMMA))-1.)
                   C = SQRT(GAMMA*P/D)

                ENDIF
             ENDIF
          ELSE

!      Left shock

             PML = PM/PL
             IF(tick)THEN
                SL  = (1.-HAVG)*UL - CL*SQRT(G2*PML + G1)
             ELSE
                SL  = UL - CL*SQRT(G2*PML + G1)
             END IF

             IF(S.LE.SL)THEN

!            Sampled point is left data state

                D = DL
                U = UL
                P = PL
                
             ELSE

!            Sampled point is Star Left state

                D = DL*(PML + G6)/(PML*G6 + 1.0)
                U = UM
                P = PM
             ENDIF
          ENDIF
       ELSE

!      Sampling point lies to the right of the contact 
!      discontinuity
          V = VR
          IF(PM.GT.PR)THEN
!         Right shock

             PMR = PM/PR
             IF(tick)then
                SR  = (1.-HAVG)*UR + CR*SQRT(G2*PMR + G1)
             ELSE
                SR  = UR + CR*SQRT(G2*PMR + G1)
             END IF


             IF(S.GE.SR)THEN
               
!            Sampled point is right data state
                D = DR
                U = UR
                P = PR
             ELSE

!            Sampled point is Star Right state

                D = DR*(PMR + G6)/(PMR*G6 + 1.0)
                U = UM
                P = PM
             ENDIF
          ELSE

!         Right rarefaction
             IF(tick)THEN
                SHR = (1.-HAVG)*UR + CR
             ELSE
                SHR = UR + CR
             END IF
             
             IF(S.GE.SHR)THEN
!            Sampled point is right data state
                D = DR
                U = UR
                P = PR
             ELSE
                CMR = CR*(PM/PR)**G1
                IF(tick)THEN
                   STR = (1.-HAVG)*UM + CMR
                ELSE
                   STR = UM + CMR
                END IF
                IF(S.LE.STR)THEN

!               Sampled point is Star Right state
                   D = DR*(PM/PR)**(1.0/GAMMA)
                   U = UM
                   P = PM
                ELSE

!               Sampled point is inside right fan

                   U = G5*(-CR + G7*UR + S)
                   C = G5*(CR - G7*(UR - S))
                   D = DR*(C/CR)**G4
                   P = PR*(C/CR)**G3

                   P = PR*(2.*(1.-HAVG)/(GAMMA-2.*HAVG+1.) - (GAMMA-1.)&
                        /((GAMMA-2.*HAVG+1.)*CR)*((1.-HAVG)*UR))**G3
                   D = DR*(P/PR)**(1./GAMMA)
                   U = UR + 2.*CR/(GAMMA-1.)*((P/PR)**((GAMMA-1.)/(2.*GAMMA))-1.)
                   C = SQRT(GAMMA*P/D)

                ENDIF
             ENDIF
          ENDIF
       ENDIF
       
       RETURN
     END SUBROUTINE SAMPLE
     
   end module riemann
