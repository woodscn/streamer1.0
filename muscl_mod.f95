module muscl_mod
contains
  pure subroutine MUSCL( Wplus1 , Wplus2 , W , Wminus1 , Wl , Wr )
! This subroutine provides the second-order MUSCL interpolation used by Hui.
!   See his 2-D paper for details.
    real , intent(in)  :: Wplus1(:) , Wplus2(:) , W(:) , Wminus1(:)
    real , intent(out) :: Wl(:) , Wr(:)
    real , allocatable :: rplus(:) , rminus(:) 
    integer :: n , i
    n = size(W,1)
    allocate( rplus(n) , rminus(n) )

    rplus  = (Wplus1-W)/(Wplus2- Wplus1)
    rminus = (Wplus1-W)/(W   -  Wminus1)

    do i = 1 , n
       Wr(i) = Wplus1(i) - 0.5*( Wplus2(i) - Wplus1(i)  )*phi(rplus(i) )
       Wl(i) =    W(i)   + 0.5*(    W(i)   - Wminus1(i) )*phi(rminus(i))
 end do
    deallocate( rplus , rminus )

  contains
    pure function phi( r )
      real , intent(in) :: r
      real :: phi
      phi = max(0.,min(1.,r))
    end function phi

  end subroutine MUSCL
end module muscl_mod


