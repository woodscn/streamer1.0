module helper_functions
contains
  function area_func( deta , dxi , n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    real                  :: area_func 
    real    , intent(in)  :: deta , dxi
    integer , intent(in)  :: n
! Since we are first advancing in the xi direction, the "area" will be deta. 
!   During the second step of Strang, the opposite is true.
    if (n .ne. 2) then
       area_func = deta
    else
       area_func = dxi
    end if
  end function area_func

  function dlambda_func( dt , n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    real                  :: dlambda_func
    real    , intent(in)  :: dt
    integer , intent(in)  :: n
! In Strang splitting, the solution is  advanced half a time step in one
!   direction, then a full step in the other, and the remaining half in
!   the first direction.
    if (n .eq. 1 .or. n .eq. 3) then
       dlambda_func = dt*0.5
    else
       dlambda_func = dt
    end if
  end function dlambda_func

  function jshift_func( n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    integer :: jshift_func
    integer , intent(in) :: n
    if (n .ne. 2) then
       jshift_func = 0
    else
       jshift_func = 1
    end if
  end function jshift_func

  function ishift_func( n )
! This function aids in the programming of the time-advancement scheme.
    implicit none
    integer :: ishift_func
    integer , intent(in) :: n
    if (n .ne. 2) then
       ishift_func = 1
    else
       ishift_func = 0
    end if
  end function ishift_func
  
end module helper_functions
