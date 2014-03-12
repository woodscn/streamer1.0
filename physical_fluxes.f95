module physical_fluxes
contains
  function fluxes(Wtemp,E,h,ishift,jshift,nxi,neta,n)
    implicit none
    integer , intent(in) :: n  , ishift , jshift , nxi , neta
    real , intent(in) :: E(:,:,:) , Wtemp(:,:,:) , h(:,:)
    real :: fluxes(neta,nxi,4)
    integer :: nx , ny
    real , pointer :: geom1(:,:) , geom2(:,:)
    target         :: E
    real , allocatable :: momx(:,:) , momy(:,:) , kexx(:,:) , kexy(:,:) , keyy(:,:) &
         , workx(:,:) , worky(:,:) , dp(:,:) , energy(:,:) , rue(:,:) &
         , rve(:,:)

    if( n .ne. 2 )then
       geom1 => E(:,:,8) ; geom2 => E(:,:,7)
    else
       geom1 => E(:,:,6) ; geom2 => E(:,:,5)
    end if
    !    if( n .ne. 2 )then
    !       geom1 =  E(:,:,8) ; geom2 = -E(:,:,7)
    !    else
    !       geom1 = -E(:,:,6) ; geom2 =  E(:,:,5)
    !    end if

    nx = nxi+ishift ; ny = neta+jshift

    allocate( momx(ny,nx) , momy(ny,nx) , kexx(ny,nx) , kexy(ny,nx) , keyy(ny,nx) &
         , workx(ny,nx) , worky(ny,nx) , dp(neta,nxi) , energy(ny,nx) , rue(ny,nx) &
         , rve(ny,nx) )

    momx = Wtemp(:,:,1)*Wtemp(:,:,2)
    momy = Wtemp(:,:,1)*Wtemp(:,:,3)
    kexx = momx*Wtemp(:,:,2)
    kexy = momx*Wtemp(:,:,3)
    keyy = momy*Wtemp(:,:,3)
    workx= Wtemp(:,:,4)*Wtemp(:,:,2)
    worky= Wtemp(:,:,4)*Wtemp(:,:,3)
    dp   = Wtemp(1+jshift:neta+jshift,1+ishift:nxi+ishift,4)-Wtemp(1:neta,1:nxi,4)
    energy = .5*(Wtemp(:,:,2)**2+Wtemp(:,:,3)**2)+1./.4*Wtemp(:,:,4)/Wtemp(:,:,1)
    rue  = momx*energy
    rve  = momy*energy

    fluxes(:,:,1) = (1.-h)&
         *(geom1*(momx(1+jshift:neta+jshift,1+ishift:nxi+ishift)-momx(1:neta,1:nxi))&
         - geom2*(momy(1+jshift:neta+jshift,1+ishift:nxi+ishift)-momy(1:neta,1:nxi)))
    fluxes(:,:,2) = (1.-h)&
         *(geom1*(kexx(1+jshift:neta+jshift,1+ishift:nxi+ishift)-kexx(1:neta,1:nxi))&
         - geom2*(kexy(1+jshift:neta+jshift,1+ishift:nxi+ishift)-kexy(1:neta,1:nxi)))&
         + geom1*(dp)
    fluxes(:,:,3) = (1.-h)&
         *(geom1*(kexy(1+jshift:neta+jshift,1+ishift:nxi+ishift)-kexy(1:neta,1:nxi))&
         - geom2*(keyy(1+jshift:neta+jshift,1+ishift:nxi+ishift)-keyy(1:neta,1:nxi)))&
         - geom2*(dp)
    fluxes(:,:,4) = (1.-h)&
         *(geom1*(rue(1+jshift:neta+jshift,1+ishift:nxi+ishift)-rue(1:neta,1:nxi))&
         - geom2*(rve(1+jshift:neta+jshift,1+ishift:nxi+ishift)-rve(1:neta,1:nxi)))&
         + geom1*(workx(1+jshift:neta+jshift,1+ishift:nxi+ishift)-workx(1:neta,1:nxi))&
         - geom2*(worky(1+jshift:neta+jshift,1+ishift:nxi+ishift)-worky(1:neta,1:nxi))

    deallocate( momx , momy , kexx , kexy , keyy , workx , worky , dp , energy , rue , rve )
    if( n .eq. 2 ) fluxes = -fluxes

  end function fluxes

  function elemental_fluxes(WL,WR,E,h,n)
    implicit none
    integer , intent(in) :: n
    real , intent(in) :: E(:) , WL(:) , WR(:) , h
    real :: elemental_fluxes(4)
    integer :: nx , ny
    real , pointer :: geom1 , geom2
    target         :: E
    real           :: momx(2) , momy(2) , kexx(2) , kexy(2) , keyy(2) &
         , workx(2) , worky(2) , dp , energy(2) , rue(2) &
         , rve(2) , Wtemp(2,4)

    Wtemp(1,:) = WL ; Wtemp(2,:) = WR;

    if( n .ne. 2 )then
       geom1 => E(8) ; geom2 => E(7)
    else
       geom1 => E(6) ; geom2 => E(5)
    end if
    !    if( n .ne. 2 )then
    !       geom1 =  E(:,:,8) ; geom2 = -E(:,:,7)
    !    else
    !       geom1 = -E(:,:,6) ; geom2 =  E(:,:,5)
    !    end if


    momx = Wtemp(:,1)*Wtemp(:,2)
    momy = Wtemp(:,1)*Wtemp(:,3)
    kexx = momx*Wtemp(:,2)
    kexy = momx*Wtemp(:,3)
    keyy = momy*Wtemp(:,3)
    workx= Wtemp(:,4)*Wtemp(:,2)
    worky= Wtemp(:,4)*Wtemp(:,3)
    dp   = Wtemp(2,4)-Wtemp(1,4)
    energy = .5*(Wtemp(:,2)**2+Wtemp(:,3)**2)+1./.4*Wtemp(:,4)/Wtemp(:,1)
    rue  = momx*energy
    rve  = momy*energy

    elemental_fluxes(1) = (1.-h)&
         *(geom1*(momx(2)-momx(1))&
         - geom2*(momy(2)-momy(1)))
!!$    write(*,*) "Left  density flux = ", (1.-h)*(-geom1*momx(1)+geom2*momy(1))
!!$    write(*,*) "Right density flux = ", (1.-h)*( geom1*momx(2)-geom2*momy(2))

    elemental_fluxes(2) = (1.-h)&
         *(geom1*(kexx(2)-kexx(1))&
         - geom2*(kexy(2)-kexy(1)))&
         + geom1*(dp)
    elemental_fluxes(3) = (1.-h)&
         *(geom1*(kexy(2)-kexy(1))&
         - geom2*(keyy(2)-keyy(1)))&
         - geom2*(dp)
    elemental_fluxes(4) = (1.-h)&
         *(geom1*(rue(2)-rue(1))&
         - geom2*(rve(2)-rve(1)))&
         + geom1*(workx(2)-workx(1))&
         - geom2*(worky(2)-worky(1))

    !		deallocate( momx , momy , kexx , kexy , keyy , workx , worky , dp , energy , rue , rve )
    if( n .eq. 2 ) elemental_fluxes = -elemental_fluxes

  end function elemental_fluxes
end module physical_fluxes
