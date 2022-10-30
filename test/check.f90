program check
   use gauss_leg_mod
   implicit none
   real(8)::res
   real(8)::a,b
   a=0.d0
   b=10.d0
   res=gauss_leg(a,b,f1)
   write(*,*)"num(1d)",res
   write(*,*)"exact  ",0.2d0*atan(2.d0)
   res=gauss_leg_nd([4,4], f2, f2bd)
   write(*,*)"num(2d)",res
   write(*,*)"exact  ",7.d0/12
   res=gauss_leg_nd([2,2,2], f3, f3bd)
   write(*,*)"num(3d)",res
   write(*,*)"exact  ",0.25d0*(exp(-1.d0)-exp(1.d0))
contains
   function f1(x) result(res)
      ! 1/(x^2+25)
      real(8),intent(in)::x
      real(8)::res
      res=1.0d0/(x*x+25.d0)
   end function f1

   function f3bd(x,j)result(res)
      ! x<-[-1,1],y<-[-1,0],z<-[0,1]
      real(8),intent(in)::x(:)
      integer,intent(in)::j
      real(8)::res(2)
      if (j==1) then
         res=[-1.d0,1.d0]
      else if (j==2) then
         res=[-1.d0,0.d0]
      else if (j==3) then
         res=[0.d0,1.d0]
      end if
   end function f3bd

   real(8) function f3(x)result(res)
      ! yz*exp(x)
      real(8),intent(in)::x(:)
      res=x(2)*x(3)*exp(x(1))
   end function f3

   function f2bd(x,j)result(res)
      ! x <-[0,1]  y<-[0,x^2+1]
      real(8),intent(in)::x(:)
      integer,intent(in)::j
      real(8)::res(2)
      if (j==1) then
         res=[0.d0,1.d0]
      else if (j==2) then
         res=[0.d0,x(1)**2+1.d0]
      end if
   end function f2bd

   real(8) function f2(x)result(res)
      ! f(x,y)=x*y
      real(8),intent(in)::x(:)
      res=x(1)*x(2)
   end function f2
end program check
