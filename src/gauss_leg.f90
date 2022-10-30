module gauss_leg_mod
  implicit none
  public:: gauss_leg,gauss_leg_nd
  private
  real(8),parameter::xi5(5)=[-0.90617984593866399281d0, -0.53846931010568309105d0,&
     & 0.0d0, 0.53846931010568309105d0, 0.90617984593866399281d0]
  real(8),parameter::wi5(5)=[0.23692688505618908751d0, 0.47862867049936646804d0,&
     & 0.56888888888888888888d0, 0.47862867049936646804d0, 0.23692688505618908751d0]
  real(8),parameter::xi6(6)=[-0.93246951420315202781d0, -0.66120938646626451366d0,&
     & -0.23861918608319690863d0, 0.23861918608319690863d0, 0.66120938646626451366d0, 0.93246951420315202781d0]
  real(8),parameter::wi6(6)=[0.17132449237917034504d0, 0.36076157304813860757d0,&
     & 0.46791393457269104739d0, 0.46791393457269104739d0, 0.36076157304813860757d0, 0.17132449237917034504d0]
  abstract interface
     function func1d(x) result(res)
        real(8),intent(in)::x
        real(8)::res
     end function func1d
     function funcnd(x) result(res)
        real(8),intent(in)::x(:)
        real(8)::res
     end function funcnd
     function funcbd(x,n) result(res)
        real(8),intent(in)::x(:)
        integer,intent(in)::n
        real(8)::res(2)
     end function funcbd
   end interface
   real(8),parameter::eps=1.0d-12
contains
   real(8)function gauss_leg_nd(js,f,fd)result(res)
   procedure(funcnd)::f
   procedure(funcbd)::fd
   integer,intent(in)::js(:)
   integer::is(2,size(js)+1)
   real(8)::d(2,size(js)+1),cc(size(js)+1)
   real(8)::x(size(js))
   real(8)::a(2),p
   integer::m,j,k,n
   n=size(js)
   m=1
   d(1,n+1)=1.0d0
   d(2,n+1)=1.0d0
   do
      do j=m,n
         a=fd(x,j)
         d(1,j)=0.5d0*(a(2)-a(1))/js(j)
         cc(j)=d(1,j)+a(1)
         x(j)=d(1,j)*xi6(1)+cc(j)
         d(2,j)=0.d0
         is(:,j)=1
      end do
      j=n
      do
         k=is(1,j)
         if (j==n) then
            p=f(x)
         else
            p=1.d0
         end if
         d(2,j)=d(2,j+1)*d(1,j+1)*p*wi6(k)+d(2,j)
         is(1,j)=is(1,j)+1
         if (is(1,j)>6) then
            if (is(2,j)>=js(j)) then
               j=j-1
               if(j==0)then
                  res=d(2,1)*d(1,1)
                  return
               end if
               cycle
            end if
            is(2,j)=is(2,j)+1
            cc(j)=cc(j)+d(1,j)*2.d0
            is(1,j)=1
         end if
         k=is(1,j)
         x(j)=d(1,j)*xi6(k)+cc(j)
         if (j/=n) exit
      end do
      m=j+1
   end do
end function gauss_leg_nd


real(8) recursive function gauss_leg(a,b,f)result(res)
   real(8),intent(in)::a,b
   procedure(func1d) :: f
   real(8)::xm
   real(8)::res5,res6
   xm=0.5d0*(a+b)
   res5=gauss_leg_node5(a,xm,f)
   res6=gauss_leg_node6(a,xm,f)
   if(abs(res5-res6)<eps)then
      res=res6
   else
      res=gauss_leg(a,xm,f)
   end if
   res5=gauss_leg_node5(xm,b,f)
   res6=gauss_leg_node6(xm,b,f)
   if(abs(res5-res6)<eps)then
      res=res+res6
   else
      res=res+gauss_leg(xm,b,f)
   end if
end function gauss_leg

function gauss_leg_node5(a,b,f)result(res)
   real(8),intent(in)::a,b
   procedure(func1d) :: f
   real(8)::xi,xm,res,w
   integer::i
   xi=0.5d0*(a+b)
   xm=0.5d0*(b-a)
   res=0.0d0
   do i=1,5
      w=wi5(i)*f(xi+xi5(i)*xm)
      res=res+w
   end do
   res=res*xm
end function gauss_leg_node5

function gauss_leg_node6(a,b,f)result(res)
   real(8),intent(in)::a,b
   procedure(func1d) :: f
   real(8)::xi,xm,res,w
   integer::i
   xi=0.5d0*(a+b)
   xm=0.5d0*(b-a)
   res=0.0d0
   do i=1,6
      w=wi6(i)*f(xi+xi6(i)*xm)
      res=res+w
   end do
   res=res*xm
end function gauss_leg_node6
end module gauss_leg_mod
