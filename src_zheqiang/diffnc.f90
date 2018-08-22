! Difference operator, node to cell
module m_diffnc
implicit none
contains

subroutine diffnc( df, f, i, a, i1, i2, oplevel, bb, x, dx1, dx2, dx3, dx )
real, intent(out) :: df(:,:,:)
real, intent(in) :: f(:,:,:,:), bb(:,:,:,:,:), x(:,:,:,:), &
    dx1(:), dx2(:), dx3(:), dx(3)
integer, intent(in) :: i, a, i1(3), i2(3), oplevel
real :: h
integer :: j, k, l, b, c

if ( any( i1 > i2 ) ) return

select case( oplevel )

! Saved B matrix, flops: 8* 7+
case( 6 )
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
    df(j,k,l) = &
      bb(j,k,l,1,a) * f(j+1,k+1,l+1,i) + f(j,k,l,i) * bb(j,k,l,5,a) &
    + bb(j,k,l,2,a) * f(j+1,k,l,i) + f(j,k+1,l+1,i) * bb(j,k,l,6,a) &
    + bb(j,k,l,3,a) * f(j,k+1,l,i) + f(j+1,k,l+1,i) * bb(j,k,l,7,a) &
    + bb(j,k,l,4,a) * f(j,k,l+1,i) + f(j+1,k+1,l,i) * bb(j,k,l,8,a)
end do
end do
end do

! Constant grid, flops: 1* 7+
case( 1 )
select case( a )
case( 1 )
    h = sign( 0.25 * dx(2) * dx(3), dx(1) )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        + f(j+1,k,l,i) - f(j,k+1,l+1,i) &
        - f(j,k+1,l,i) + f(j+1,k,l+1,i) &
        - f(j,k,l+1,i) + f(j+1,k+1,l,i) )
    end do
    end do
    end do
case( 2 )
    h = sign( 0.25 * dx(3) * dx(1), dx(2) )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        - f(j+1,k,l,i) + f(j,k+1,l+1,i) &
        + f(j,k+1,l,i) - f(j+1,k,l+1,i) &
        - f(j,k,l+1,i) + f(j+1,k+1,l,i) )
    end do
    end do
    end do
case( 3 )
    h = sign( 0.25 * dx(1) * dx(2), dx(3) )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        - f(j+1,k,l,i) + f(j,k+1,l+1,i) &
        - f(j,k+1,l,i) + f(j+1,k,l+1,i) &
        + f(j,k,l+1,i) - f(j+1,k+1,l,i) )
    end do
    end do
    end do
end select
  
! Rectangular grid, flops: 2* 7+
case( 2 )
h = sign( 0.25, product( dx ) )
select case( a )
case( 1 )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * dx2(k) * dx3(l) * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        + f(j+1,k,l,i) - f(j,k+1,l+1,i) &
        - f(j,k+1,l,i) + f(j+1,k,l+1,i) &
        - f(j,k,l+1,i) + f(j+1,k+1,l,i) )
    end do
    end do
    end do
case( 2 )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * dx3(l) * dx1(j) * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        - f(j+1,k,l,i) + f(j,k+1,l+1,i) &
        + f(j,k+1,l,i) - f(j+1,k,l+1,i) &
        - f(j,k,l+1,i) + f(j+1,k+1,l,i) )
    end do
    end do
    end do
case( 3 )
    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)
        df(j,k,l) = h * dx1(j) * dx2(k) * &
        ( f(j+1,k+1,l+1,i) - f(j,k,l,i) &
        - f(j+1,k,l,i) + f(j,k+1,l+1,i) &
        - f(j,k+1,l,i) + f(j+1,k,l+1,i) &
        + f(j,k,l+1,i) - f(j+1,k+1,l,i) )
    end do
    end do
    end do
end select

! Parallelepiped grid, flops: 17* 27+
case( 3 )
h = sign( 0.25, product( dx ) )
b = modulo( a, 3 ) + 1
c = modulo( a + 1, 3 ) + 1
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
df(j,k,l) = h * &
((f(j+1,k+1,l+1,i)-f(j,k,l,i))* &
    (x(j,k+1,l+1,b)*(x(j+1,k,l+1,c)-x(j+1,k+1,l,c)) &
    +x(j+1,k,l+1,b)*(x(j+1,k+1,l,c)-x(j,k+1,l+1,c)) &
    +x(j+1,k+1,l,b)*(x(j,k+1,l+1,c)-x(j+1,k,l+1,c))) &
+(f(j+1,k,l,i)-f(j,k+1,l+1,i))* &
    (x(j,k,l,b)*(x(j+1,k+1,l,c)-x(j+1,k,l+1,c)) &
    +x(j+1,k,l+1,b)*(x(j,k,l,c)-x(j+1,k+1,l,c)) &
    +x(j+1,k+1,l,b)*(x(j+1,k,l+1,c)-x(j,k,l,c))) &
+(f(j,k+1,l,i)-f(j+1,k,l+1,i))* &
    (x(j,k,l,b)*(x(j,k+1,l+1,c)-x(j+1,k+1,l,c)) &
    +x(j+1,k+1,l,b)*(x(j,k,l,c)-x(j,k+1,l+1,c)) &
    +x(j,k+1,l+1,b)*(x(j+1,k+1,l,c)-x(j,k,l,c))) &
+(f(j,k,l+1,i)-f(j+1,k+1,l,i))* &
    (x(j,k,l,b)*(x(j+1,k,l+1,c)-x(j,k+1,l+1,c)) &
    +x(j,k+1,l+1,b)*(x(j,k,l,c)-x(j+1,k,l+1,c)) &
    +x(j+1,k,l+1,b)*(x(j,k+1,l+1,c)-x(j,k,l,c))))
end do
end do
end do

! General grid one-point quadrature, flops: 17* 63+
case( 4 )
h = sign( 0.0625, product( dx ) )
b = modulo( a, 3 ) + 1
c = modulo( a + 1, 3 ) + 1
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
df(j,k,l) = h * &
((f(j+1,k+1,l+1,i)-f(j,k,l,i))* &
    ((x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j,k+1,l,c)-x(j+1,k,l+1,c)-x(j,k,l+1,c)+x(j+1,k+1,l,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j,k,l+1,c)-x(j+1,k+1,l,c)-x(j+1,k,l,c)+x(j,k+1,l+1,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j+1,k,l,c)-x(j,k+1,l+1,c)-x(j,k+1,l,c)+x(j+1,k,l+1,c))) &
+(f(j+1,k,l,i)-f(j,k+1,l+1,i))* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j,k,l+1,c)-x(j+1,k+1,l,c)-x(j,k+1,l,c)+x(j+1,k,l+1,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j+1,k+1,l+1,c)-x(j,k,l,c)-x(j,k,l+1,c)+x(j+1,k+1,l,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j,k+1,l,c)-x(j+1,k,l+1,c)-x(j+1,k+1,l+1,c)+x(j,k,l,c))) &
+(f(j,k+1,l,i)-f(j+1,k,l+1,i))* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j+1,k,l,c)-x(j,k+1,l+1,c)-x(j,k,l+1,c)+x(j+1,k+1,l,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j+1,k+1,l+1,c)-x(j,k,l,c)-x(j+1,k,l,c)+x(j,k+1,l+1,c)) &
    +(x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j,k,l+1,c)-x(j+1,k+1,l,c)-x(j+1,k+1,l+1,c)+x(j,k,l,c))) &
+(f(j,k,l+1,i)-f(j+1,k+1,l,i))* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j,k+1,l,c)-x(j+1,k,l+1,c)-x(j+1,k,l,c)+x(j,k+1,l+1,c)) &
    +(x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j+1,k+1,l+1,c)-x(j,k,l,c)-x(j,k+1,l,c)+x(j+1,k,l+1,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j+1,k,l,c)-x(j,k+1,l+1,c)-x(j+1,k+1,l+1,c)+x(j,k,l,c))))
end do
end do
end do

! General grid exact, flops: 57* 119+
case( 5 )
h = sign( 1.0 / 12.0, product( dx ) )
b = modulo( a, 3 ) + 1
c = modulo( a + 1, 3 ) + 1
do l = i1(3), i2(3)
do k = i1(2), i2(2)
do j = i1(1), i2(1)
df(j,k,l) = h * &
(f(j+1,k+1,l+1,i)* &
    ((x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j+1,k+1,l,c)-x(j+1,k,l+1,c))+x(j,k+1,l+1,b)*(x(j,k,l+1,c)-x(j,k+1,l,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j,k+1,l+1,c)-x(j+1,k+1,l,c))+x(j+1,k,l+1,b)*(x(j+1,k,l,c)-x(j,k,l+1,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j+1,k,l+1,c)-x(j,k+1,l+1,c))+x(j+1,k+1,l,b)*(x(j,k+1,l,c)-x(j+1,k,l,c))) &
+f(j+1,k,l,i)* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j+1,k,l+1,c)-x(j+1,k+1,l,c))+x(j,k,l,b)*(x(j,k+1,l,c)-x(j,k,l+1,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j+1,k+1,l,c)-x(j,k,l,c))+x(j+1,k,l+1,b)*(x(j,k,l+1,c)-x(j+1,k+1,l+1,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j,k,l,c)-x(j+1,k,l+1,c))+x(j+1,k+1,l,b)*(x(j+1,k+1,l+1,c)-x(j,k+1,l,c))) &
+f(j,k+1,l,i)* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j+1,k+1,l,c)-x(j,k+1,l+1,c))+x(j,k,l,b)*(x(j,k,l+1,c)-x(j+1,k,l,c)) &
    +(x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j,k,l,c)-x(j+1,k+1,l,c))+x(j,k+1,l+1,b)*(x(j+1,k+1,l+1,c)-x(j,k,l+1,c)) &
    +(x(j,k,l+1,b)-x(j+1,k+1,l,b))*(x(j,k+1,l+1,c)-x(j,k,l,c))+x(j+1,k+1,l,b)*(x(j+1,k,l,c)-x(j+1,k+1,l+1,c))) &
+f(j,k,l+1,i)* &
    ((x(j+1,k+1,l+1,b)-x(j,k,l,b))*(x(j,k+1,l+1,c)-x(j+1,k,l+1,c))+x(j,k,l,b)*(x(j+1,k,l,c)-x(j,k+1,l,c)) &
    +(x(j+1,k,l,b)-x(j,k+1,l+1,b))*(x(j+1,k,l+1,c)-x(j,k,l,c))+x(j,k+1,l+1,b)*(x(j,k+1,l,c)-x(j+1,k+1,l+1,c)) &
    +(x(j,k+1,l,b)-x(j+1,k,l+1,b))*(x(j,k,l,c)-x(j,k+1,l+1,c))+x(j+1,k,l+1,b)*(x(j+1,k+1,l+1,c)-x(j+1,k,l,c))) &
+f(j,k,l,i)* &
    ((x(j,k+1,l+1,b)-x(j+1,k,l,b))*(x(j,k+1,l,c)-x(j,k,l+1,c))+x(j+1,k,l,b)*(x(j+1,k,l+1,c)-x(j+1,k+1,l,c)) &
    +(x(j+1,k,l+1,b)-x(j,k+1,l,b))*(x(j,k,l+1,c)-x(j+1,k,l,c))+x(j,k+1,l,b)*(x(j+1,k+1,l,c)-x(j,k+1,l+1,c)) &
    +(x(j+1,k+1,l,b)-x(j,k,l+1,b))*(x(j+1,k,l,c)-x(j,k+1,l,c))+x(j,k,l+1,b)*(x(j,k+1,l+1,c)-x(j+1,k,l+1,c))) &
+f(j,k+1,l+1,i)* &
    ((x(j,k,l,b)-x(j+1,k+1,l+1,b))*(x(j,k,l+1,c)-x(j,k+1,l,c))+x(j+1,k+1,l+1,b)*(x(j+1,k+1,l,c)-x(j+1,k,l+1,c)) &
    +(x(j+1,k,l+1,b)-x(j,k+1,l,b))*(x(j+1,k+1,l+1,c)-x(j,k,l+1,c))+x(j,k+1,l,b)*(x(j,k,l,c)-x(j+1,k+1,l,c)) &
    +(x(j+1,k+1,l,b)-x(j,k,l+1,b))*(x(j,k+1,l,c)-x(j+1,k+1,l+1,c))+x(j,k,l+1,b)*(x(j+1,k,l+1,c)-x(j,k,l,c))) &
+f(j+1,k,l+1,i)* &
    ((x(j,k,l,b)-x(j+1,k+1,l+1,b))*(x(j+1,k,l,c)-x(j,k,l+1,c))+x(j+1,k+1,l+1,b)*(x(j,k+1,l+1,c)-x(j+1,k+1,l,c)) &
    +(x(j,k+1,l+1,b)-x(j+1,k,l,b))*(x(j,k,l+1,c)-x(j+1,k+1,l+1,c))+x(j+1,k,l,b)*(x(j+1,k+1,l,c)-x(j,k,l,c)) &
    +(x(j+1,k+1,l,b)-x(j,k,l+1,b))*(x(j+1,k+1,l+1,c)-x(j+1,k,l,c))+x(j,k,l+1,b)*(x(j,k,l,c)-x(j,k+1,l+1,c))) &
+f(j+1,k+1,l,i)* &
    ((x(j,k,l,b)-x(j+1,k+1,l+1,b))*(x(j,k+1,l,c)-x(j+1,k,l,c))+x(j+1,k+1,l+1,b)*(x(j+1,k,l+1,c)-x(j,k+1,l+1,c)) &
    +(x(j,k+1,l+1,b)-x(j+1,k,l,b))*(x(j+1,k+1,l+1,c)-x(j,k+1,l,c))+x(j+1,k,l,b)*(x(j,k,l,c)-x(j+1,k,l+1,c)) &
    +(x(j+1,k,l+1,b)-x(j,k+1,l,b))*(x(j+1,k,l,c)-x(j+1,k+1,l+1,c))+x(j,k+1,l,b)*(x(j,k+1,l+1,c)-x(j,k,l,c))))
end do
end do
end do

case default; stop 'illegal operator'

end select

end subroutine

end module

