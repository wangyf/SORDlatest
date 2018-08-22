! Stress calculation
module m_stress
implicit none
contains

subroutine stress
use m_globals
use m_diffnc
use m_source
use m_util
use m_fieldio
use m_stats
integer :: i1(3), i2(3), i, j, k, l, ic, iid, id, p

if ( verb ) write( 0, * ) 'Stress'

! Modified displacement
do i = 1, 3
    w1(:,:,:,i) = uu(:,:,:,i) + gam * vv(:,:,:,i)
end do
call set_halo( s1, 0.0, i1cell, i2cell )

! Loop over component and derivative direction
doic: do ic  = 1, 3
doid: do iid = 1, 3; id = modulo( ic + iid - 1, 3 ) + 1

! Elastic region: g_ij = (u_i + gamma*v_i),j
i1 = max( i1pml + 1, i1cell )
i2 = min( i2pml - 2, i2cell )
call diffnc( s1, w1, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )

! PML region, non-damped directions: g_ij = u_i,j
do i = 1, 3
if ( id /= i ) then
    i1 = i1cell
    i2 = i2cell
    i2(i) = min( i2(i), i1pml(i) )
    call diffnc( s1, uu, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    i1 = i1cell
    i2 = i2cell
    i1(i) = max( i1(i), i2pml(i) - 1 )
    call diffnc( s1, uu, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
end if
end do

! PML region, damped direction: g'_ij = d_j*g_ij = v_i,j
select case( id )
case( 1 )
    i1 = i1cell
    i2 = i2cell
    i2(1) = min( i2(1), i1pml(1) )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do j = i1(1), i2(1)
        i = j - i1(1) + 1
        p = j + nnoff(1)
        do l = i1(3), i2(3)
        do k = i1(2), i2(2)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g1(i,k,l,ic)
            g1(i,k,l,ic) = s1(j,k,l)
        end do
        end do
    end do
    i1 = i1cell
    i2 = i2cell
    i1(1) = max( i1(1), i2pml(1) - 1 )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do j = i1(1), i2(1)
        i = i2(1) - j + 1
        p = nn(1) - j - nnoff(1)
        do l = i1(3), i2(3)
        do k = i1(2), i2(2)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g4(i,k,l,ic)
            g4(i,k,l,ic) = s1(j,k,l)
        end do
        end do
    end do
case( 2 )
    i1 = i1cell
    i2 = i2cell
    i2(2) = min( i2(2), i1pml(2) )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do k = i1(2), i2(2)
        i = k - i1(2) + 1
        p = k + nnoff(2)
        do l = i1(3), i2(3)
        do j = i1(1), i2(1)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g2(j,i,l,ic)
            g2(j,i,l,ic) = s1(j,k,l)
        end do
        end do
    end do
    i1 = i1cell
    i2 = i2cell
    i1(2) = max( i1(2), i2pml(2) - 1 )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do k = i1(2), i2(2)
        i = i2(2) - k + 1
        p = nn(2) - k - nnoff(2)
        do l = i1(3), i2(3)
        do j = i1(1), i2(1)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g5(j,i,l,ic)
            g5(j,i,l,ic) = s1(j,k,l)
        end do
        end do
    end do
case( 3 )
    i1 = i1cell
    i2 = i2cell
    i2(3) = min( i2(3), i1pml(3) )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do l = i1(3), i2(3)
        i = l - i1(3) + 1
        p = l + nnoff(3)
        do k = i1(2), i2(2)
        do j = i1(1), i2(1)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g3(j,k,i,ic)
            g3(j,k,i,ic) = s1(j,k,l)
        end do
        end do
    end do
    i1 = i1cell
    i2 = i2cell
    i1(3) = max( i1(3), i2pml(3) - 1 )
    call diffnc( s1, vv, ic, id, i1, i2, oplevel, bb, xx, dx1, dx2, dx3, dx )
    do l = i1(3), i2(3)
        i = i2(3) - l + 1
        p = nn(3) - l - nnoff(3)
        do k = i1(2), i2(2)
        do j = i1(1), i2(1)
            s1(j,k,l) = dc2(p) * s1(j,k,l) + dc1(p) * g6(j,k,i,ic)
            g6(j,k,i,ic) = s1(j,k,l)
        end do
        end do
    end do
end select

! Add contribution to potency
i = 6 - ic - id
if ( ic < id ) then
    w2(:,:,:,i) = 0.5 * s1
elseif ( ic > id ) then
    w2(:,:,:,i) = w2(:,:,:,i) + 0.5 * s1
else
    w1(:,:,:,ic) = s1
end if

end do doid
end do doic

! Strain
do i = 1, 3
    w1(:,:,:,i) = w1(:,:,:,i) * vc
    w2(:,:,:,i) = w2(:,:,:,i) * vc
end do

! Add potency source to strain
if ( source == 'potency' ) then
    call finite_source
    call tensor_point_source
end if

! Strain I/O
call fieldio( '<>', 'e11', w1(:,:,:,1) )
call fieldio( '<>', 'e22', w1(:,:,:,2) )
call fieldio( '<>', 'e33', w1(:,:,:,3) )
call fieldio( '<>', 'e23', w2(:,:,:,1) )
call fieldio( '<>', 'e31', w2(:,:,:,2) )
call fieldio( '<>', 'e12', w2(:,:,:,3) )

! Attenuation
!do j = 1, 2
!do k = 1, 2
!do l = 1, 2
!  i = j + 2 * ( k - 1 ) + 4 * ( l - 1 )
!  z1(j::2,k::2,l::2,:) = c1(i) * z1(j::2,k::2,l::2,:) + c2(i) * w1(j::2,k::2,l::2,:)
!  z2(j::2,k::2,l::2,:) = c1(i) * z2(j::2,k::2,l::2,:) + c2(i) * w2(j::2,k::2,l::2,:)
!end do
!end do
!end do

! Hooke's Law: w_ij = lam*g_ij*delta_ij + mu*(g_ij + g_ji)
if ( eplasticity == 'plastic' ) then 
    w1 = w1 - ep1
    w2 = w2 - ep2
end if

s1 = w1(:,:,:,1) + w1(:,:,:,2) + w1(:,:,:,3)

do i = 1, 3
    w1(:,:,:,i) = 2.0 * mu * w1(:,:,:,i) + lam * s1
    w2(:,:,:,i) = 2.0 * mu * w2(:,:,:,i)
end do

! Add moment source to stress
if ( source == 'moment' ) then
    call finite_source
    call tensor_point_source
end if

if ( eplasticity == 'plastic' ) then 

w1 = w1 + si1
w2 = w2 + si2

s1 = w1(:,:,:,1) + w1(:,:,:,2) + w1(:,:,:,3)
s2 = 0.

do i = 1, 3
    z1(:,:,:,i) = w1(:,:,:,i) - s1 / 3.
    z2(:,:,:,i) = w2(:,:,:,i)
    s2 = s2 + z1(:,:,:,i) * z1(:,:,:,i) + 2.0 * z2(:,:,:,i) * z2(:,:,:,i)
end do

s2 = sqrt( 0.5 * s2 )
call invert( s2 )
r3 = - s1 / 3.0 * sin( phi ) + mco * cos( phi )
r3 = r3 * s2

if ( tv > 1.e-6 ) then
   s2 = r3 + ( 1.0 - r3 ) * exp( - dt / tv )
else
   s2 = r3
end if

do i = 1, 3
    r1 = z1(:,:,:,i)
    r2 = z2(:,:,:,i)
    where ( r3 < 1.0 )
        r4 = r1 * s2 
        r5 = r2 * s2
        r1 = r1 - r4
        r2 = r2 - r5
    elsewhere
        r4 = r1
        r5 = r2
        r1 = 0.0
        r2 = 0.0
    end where

    dep1(:,:,:,i) = 0.5 * r1 * mur 
    dep2(:,:,:,i) = 0.5 * r2 * mur 

    w1(:,:,:,i) = s1 / 3.0 + r4
    w2(:,:,:,i) = r5
end do

s2 = 0.0
do i = 1, 3
    s2 = s2 + dep1(:,:,:,i) * dep1(:,:,:,i) + 2.0 * dep2(:,:,:,i) * dep2(:,:,:,i)
end do
epm = epm + sqrt( 2.0 / 3.0 * s2 )

call fieldio( '<>', 'epm', epm )

ep1 = ep1 + dep1
ep2 = ep2 + dep2

w1 = w1 - si1
w2 = w2 - si2

end if

if ( ivols == 'yes' ) then
    z1 = w1 + si1
    z2 = w2 + si2 
 
    call fieldio( '<>', 'w11', z1(:,:,:,1) )
    call fieldio( '<>', 'w22', z1(:,:,:,2) )
    call fieldio( '<>', 'w33', z1(:,:,:,3) )
    call fieldio( '<>', 'w23', z2(:,:,:,1) )
    call fieldio( '<>', 'w31', z2(:,:,:,2) )
    call fieldio( '<>', 'w12', z2(:,:,:,3) )
else
    call fieldio( '<>', 'w11', w1(:,:,:,1) )
    call fieldio( '<>', 'w22', w1(:,:,:,2) )
    call fieldio( '<>', 'w33', w1(:,:,:,3) )
    call fieldio( '<>', 'w23', w2(:,:,:,1) )
    call fieldio( '<>', 'w31', w2(:,:,:,2) )
    call fieldio( '<>', 'w12', w2(:,:,:,3) )
end if

if ( modulo( it, itstats ) == 0 ) then
    call tensor_norm( s1, w1, w2, i1core, i2core, (/ 1, 1, 1 /) )
    call set_halo( s1, -1.0, i1core, i2core )
    wmax = maxval( s1 )
end if
call fieldio( '>', 'wm2', s1  )

end subroutine

end module

