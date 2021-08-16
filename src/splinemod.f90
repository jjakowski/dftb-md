module splinemod
 use  utilmod
 implicit none

contains 

 subroutine  gety2tab(N,dx,ytab,y2tab) 
 !
 !  natural splines  with  y"  of the first and last point assumed =0
   implicit none 
   integer , intent (in) :: N
!   real (kind=8), intent(in)   ::   dx, ytab(N)
!   real (kind=8), intent(out)  ::      y2tab(N)
   real (kind=8), intent(in)   ::   dx, ytab(*)
   real (kind=8), intent(out)  ::      y2tab(*)
   real (kind=8)               ::  temp 
   real (kind=8), dimension(:,:), allocatable :: A
   real (kind=8), dimension(:)  , allocatable :: y2
   integer ::  i,j,k , M  

   M = N-2
   allocate (A(M,M), y2(M))
   A =0.0d0
   do  i=1,M-1
      A(i,i+1) = 1.0d0   
      A(i+1,i) = 1.0d0   
      A(i,i)   = 4.0d0
      enddo 
    A(M,M) = 4.0d0
   
   temp  =  6.0d0/(dx*dx) 
   do i=1,M
    y2(i) =  temp * (ytab(i) -2.0d0* ytab(i+1) + ytab(i+2) )
    enddo

! forward loop
   do i=2,M 
    temp = 1.0d0 / A(i-1,i-1) 
     !A(i,:) = A(i,:) - A(i-1,:)*temp
     A(i,i-1:i) = A(i,i-1:i) - A(i-1,i-1:i)*temp
     y2(i) =  y2(i) - y2(i-1)*temp   
     enddo

! backward loop
  do i=M,2,-1
    temp =  1.0d0 / A(i,i)
    y2(i)= y2(i) * temp 
    y2(i-1) = y2(i-1) -  y2(i)   
    enddo 
  y2(1) =   y2(1) /A(1,1)  
  
! finalize second  derivatives    
  y2tab(1) =0.d0
  y2tab(N) = 0.0d0
  y2tab(2:N-1) =  y2

  ! write(*,*) '*** xtab  ';  call  dumpm(xtab)
  ! write(*,*) '*** ytab  ';  call  dumpm(ytab)
  ! write(*,*) '*** y2tab ';  call  dumpm(y2tab)
   deallocate  (A,y2)
 end  subroutine

 subroutine cubics(x,y,dx,y0,ypp0,y1,ypp1)
 !------------------------
 !   A = (x(i+1) - x)  / h 
 !   B = (x   -x(i))   / h
 !   C = (A^3 -A)*h^2  / 6 
 !   D = (B^3 -B)*h^2  / 6
 !!!!
 ! y(x)  =  A*y(i)  + B*y(i+1)  +  C*y2(i) +D*y2(i+1)
 ! y'(x) = (y(i+1)-y(i)) / h   +   (3A^2-1)/6 * h* y2(j)   +   (3B^2 -1)/6 *h *y2(j+1)
 ! y"(x) =  A*y2(j) + B*y2(j+1)
 !-----------------------
    real (kind=8), intent (in)   ::  x,dx,y0,y1, ypp0, ypp1 
    real (kind=8), intent (out)  ::  y 
    real (kind=8)  :: A, B, C, D ,dxi, dx2, sixth 
  
    if (x>dx ) then; write(*,*) 'WARNING: Problem with cubics: x> dx =',dx,x;  endif 
    dxi = 1.0d0/dx 
    sixth = 1.0d0/6.0d0
    dx2 = dx*dx 
    A =  (dx -x)* dxi       !B =  (x  - xtab(i) )* dxi
    B = 1.0d0 - A 
    C = sixth*(A*A*A -A) *dx2 
    D = sixth*(B*B*B -B) *dx2
    y = A*y0 +  B*y1 + C*ypp0  +  D*ypp1 
    end subroutine cubics 

 subroutine dcubics(x,yp,dx,y0,ypp0,y1,ypp1)
 !------------------------------
 ! y' = (y(i+1)-y(i)) / h   +   (3A^2-1)/6 * h* y2(j)   +   (3B^2 -1)/6 *h *y2(j+1)
 ! y" =  A*y2(j) + B*y2(j+1)
 !---------------------------------
    real (kind=8), intent (in)   ::  x,dx,y0,y1, ypp0, ypp1 
    real (kind=8), intent (out)  ::  yp 
    real (kind=8)  :: A, B, dxi, dx6i
  
    if (x>dx ) then; write(*,*) 'Problem with cubics: x> dx =',dx,x; stop ; endif 
    dxi = 1.0d0/dx 
    dx6i = dx/6.0d0
    A =  (dx -x)* dxi       !B =  (x  - xtab(i) )* dxi
    B =  1.0d0 - A 
    !yp  = (y1- y0) *dxi ! +  (3.0d0*A*A -1.0d0)*dx6i * ypp0  +  (3.0d0*B*B-1.0d0)*dx6i * ypp1 
    yp  = 1.0d0*(y1- y0) *dxi    -  (3.0d0*A*A -1.0d0)*dx6i * ypp0  +  (3.0d0*B*B-1.0d0)*dx6i * ypp1 
  end subroutine dcubics 


  subroutine d2cubics(x,ypp,dx,y0,ypp0,y1,ypp1)
  !------------------------------
  ! y" =  A*y2(j) + B*y2(j+1)
  !---------------------------------
    real (kind=8), intent (in)   ::  x,dx,y0,y1, ypp0, ypp1 
    real (kind=8), intent (out)  ::  ypp 
    real (kind=8)  :: A, B, dxi, dx2, sixth 

    if (x>dx ) then; write(*,*) 'Problem with cubics: x> dx =',dx,x; stop ; endif 
    dxi = 1.0d0/dx 
    A =  (dx -x)* dxi       
    B =  1.0d0 - A             ! B =  (x  - xtab(i) )* dxi
    ypp   =  A*ypp0   +B* ypp1 
  end subroutine d2cubics 

end  module splinemod 



