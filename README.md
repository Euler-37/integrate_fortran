# multi dimension integrate

## run

```sh
  fpm test
```

``` sh
 num(1d)  0.22142974355881528
 exact    0.22142974355881809
 num(2d)  0.58333333333333348
 exact    0.58333333333333337
 num(3d) -0.58760059682190058
 exact   -0.58760059682190069
```

## usage
``` toml
[dependencies]
gauss_leg.git="https://github.com/Euler-37/integrate_fortran"
```
### 1d

user-defined function `f(x)`, and $\int_a^b f(x) dx$, `res=gauss_leg(a,b,f)`.
``` fortran
     function func1d(x) result(res)
        real(8),intent(in)::x
        real(8)::res
     end function func1d
```

### nd

user-defined function `f(x)`,and boundary condition `fs(x,k)`,and space slices.

``` fortran
     function funcnd(x) result(res)
        real(8),intent(in)::x(:)
        real(8)::res
     end function funcnd
     function funcbd(x,k) result(res)
        real(8),intent(in)::x(:)
        integer,intent(in)::k
        real(8)::res(2)
     end function funcbd
```

## example

$$\int_0^{10}1/(x^2+25)=0.2atan(2)$$

$$\int_{-1}^{1}\int_{-1}^{0}\int_{0}^{1} yzexp(x)dxdydz=0.25(e^{-1}-e^{1})$$

$$\int_{0}^{1}\int_{0}^{x^2+1} xy dxdy=7/12$$

``` fortran
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
```



[1] https://book.douban.com/subject/20430175/
