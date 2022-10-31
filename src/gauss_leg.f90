module gauss_leg_mod
   use iso_fortran_env,only:rk=>real32
   implicit none
   public:: gauss_leg,gauss_leg_nd
   private
   real(rk),parameter::xi5(5)=[-0.90617984593866399281_rk, -0.53846931010568309105_rk,&
      & 0.0_rk, 0.53846931010568309105_rk, 0.90617984593866399281_rk]
   real(rk),parameter::wi5(5)=[0.23692688505618908751_rk, 0.47862867049936646804_rk,&
      & 0.56888888888888888888_rk, 0.47862867049936646804_rk, 0.23692688505618908751_rk]
   real(rk),parameter::xi6(6)=[-0.93246951420315202781_rk, -0.66120938646626451366_rk,&
      & -0.23861918608319690863_rk, 0.23861918608319690863_rk, 0.66120938646626451366_rk, 0.93246951420315202781_rk]
   real(rk),parameter::wi6(6)=[0.17132449237917034504_rk, 0.36076157304813860757_rk,&
      & 0.46791393457269104739_rk, 0.46791393457269104739_rk, 0.36076157304813860757_rk, 0.17132449237917034504_rk]
   abstract interface
      function func1d(x) result(res)
         import rk
         real(rk),intent(in)::x
         real(rk)::res
      end function func1d
      function funcnd(x) result(res)
         import rk
         real(rk),intent(in)::x(:)
         real(rk)::res
      end function funcnd
      function funcbd(x,n) result(res)
         import rk
         real(rk),intent(in)::x(:)
         integer,intent(in)::n
         real(rk)::res(2)
      end function funcbd
   end interface
   real(rk),parameter::eps=epsilon(1.0_rk)*10
contains
   real(rk)function gauss_leg_nd(js,f,fd)result(res)
      procedure(funcnd)::f
      procedure(funcbd)::fd
      integer,intent(in)::js(:)
      integer::is(2,size(js)+1)
      real(rk)::d(2,size(js)+1),cc(size(js)+1)
      real(rk)::x(size(js))
      real(rk)::a(2),p
      integer::m,j,k,n
      n=size(js)
      m=1
      d(1,n+1)=1.0_rk
      d(2,n+1)=1.0_rk
      do
         do j=m,n
            a=fd(x,j)
            d(1,j)=0.5_rk*(a(2)-a(1))/js(j)
            cc(j)=d(1,j)+a(1)
            x(j)=d(1,j)*xi6(1)+cc(j)
            d(2,j)=0._rk
            is(:,j)=1
         end do
         j=n
         do
            k=is(1,j)
            if (j==n) then
               p=f(x)
            else
               p=1._rk
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
               cc(j)=cc(j)+d(1,j)*2._rk
               is(1,j)=1
            end if
            k=is(1,j)
            x(j)=d(1,j)*xi6(k)+cc(j)
            if (j/=n) exit
         end do
         m=j+1
      end do
   end function gauss_leg_nd


   real(rk) recursive function gauss_leg(a,b,f)result(res)
      real(rk),intent(in)::a,b
      procedure(func1d) :: f
      real(rk)::xm
      real(rk)::res5,res6
      xm=0.5_rk*(a+b)
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
      real(rk),intent(in)::a,b
      procedure(func1d) :: f
      real(rk)::xi,xm,res,w
      integer::i
      xi=0.5_rk*(a+b)
      xm=0.5_rk*(b-a)
      res=0.0_rk
      do i=1,5
         w=wi5(i)*f(xi+xi5(i)*xm)
         res=res+w
      end do
      res=res*xm
   end function gauss_leg_node5

   function gauss_leg_node6(a,b,f)result(res)
      real(rk),intent(in)::a,b
      procedure(func1d) :: f
      real(rk)::xi,xm,res,w
      integer::i
      xi=0.5_rk*(a+b)
      xm=0.5_rk*(b-a)
      res=0.0_rk
      do i=1,6
         w=wi6(i)*f(xi+xi6(i)*xm)
         res=res+w
      end do
      res=res*xm
   end function gauss_leg_node6
end module gauss_leg_mod
