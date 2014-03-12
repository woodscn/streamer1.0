module write_files_mod
! This module is responsible for outputting data into a form that Matlab can 
!   read for data visualization. 
contains
  subroutine write_files(x,y,E,h,t,dt)
! Inputs: x, y, E (array of conserved variables), t, dt
! Outputs: two_D.dat, an ASCII file that contains the primitive variables together
!          with their coordinate values and the associated time. Data is designed
!          to be read by the matlab file geom_data_reader.m. It should be noted 
!          that the file is opened and headers are written in the main program.
    use vars
    implicit none
    real , intent(in)  :: x(:,:) , y(:,:) , E(:,:,:) , t , dt , h(:,:)
    integer            :: i , j , neta , nxi
    real , allocatable :: W(:,:,:) , u(:,:) , v(:,:) 

    neta = size(E,1) ; nxi = size(E,2)
    allocate( W(neta,nxi,4) , u(neta,nxi) , v(neta,nxi) )
! We output primitive variables
    call constoprim( E(:,:,1:4) , W , E(:,:,5:8) )
    u = E(:,:,2)/E(:,:,1) ; v = E(:,:,3)/E(:,:,1)
    write(1,*) ' ' 
    write(1,*) 't=',t,'nx=',nxi,'ny=',neta,'dt=',dt
    write(1,*) ' '
    write(1,*) 'x= ','y= ','u= ','v= ','rho= ','P= '
    do i = 1 , nxi
       do j = 1 , neta
          write(1,*) ' ', x(j,i) , y(j,i) , u(j,i) , v(j,i) , W(j,i,1) , W(j,i,4) &
               , E(j,i,5) , E(j,i,6) , E(j,i,7) , E(j,i,8) , h(j,i)
       end do
    end do
    deallocate( W , u , v )
  end subroutine write_files
end module write_files_mod
 
