! Rupture boundary condition
module m_rupture
implicit none
contains

! Rupture initialization
subroutine rupture_init
use m_globals
use m_collective
use m_surfnormals
use m_inivolstress
use m_util
use m_fieldio
use m_stats
real :: xhypo(3), xi(3), w  !, rr
integer :: i1(3), i2(3), i, j, k, l

if ( ifn /= 0 ) then
    i1 = i1core
    i2 = i2core
    i1(ifn) = irup
    i2(ifn) = irup
    call nodenormals( nhat, w1, dx, i1, i2, ifn )
    area = sign( 1, faultnormal ) * sqrt( sum( nhat * nhat, 4 ) )
 
    t0 = 0.
end if 

if ( ivols == 'yes' .and. ifn /= 0) then
    call inivolstress
    f1 = area
    call invert( f1 )
    do i = 1, 3
        nhat(:,:,:,i) = nhat(:,:,:,i) * f1
    end do
    t1 = 0.0
    t2 = 0.0
    t3 = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Added by Yongfei Wang May,17 2016
! recalculate strike direction in t1(:,:,:,1:3), along dip t2(:,:,:,1:3)
! slipvector only choose from (1,0,0),(0,1,0),(0,0,1)
! if left-lateral fault, we don't change slipvector, make ts minus
! start:
    t2(:,:,:,1) = nhat(:,:,:,2) * slipvector(3) - nhat(:,:,:,3) * slipvector(2)
    t2(:,:,:,2) = nhat(:,:,:,3) * slipvector(1) - nhat(:,:,:,1) * slipvector(3)
    t2(:,:,:,3) = nhat(:,:,:,1) * slipvector(2) - nhat(:,:,:,2) * slipvector(1)
    f1 = sqrt( sum( t2 * t2, 4 ) )
    call invert( f1 )
    do i = 1, 3
        t2(:,:,:,i) = t2(:,:,:,i) * f1
    end do

    t1(:,:,:,1) = t2(:,:,:,2) * nhat(:,:,:,3) - t2(:,:,:,3) * nhat(:,:,:,2)
    t1(:,:,:,2) = t2(:,:,:,3) * nhat(:,:,:,1) - t2(:,:,:,1) * nhat(:,:,:,3)
    t1(:,:,:,3) = t2(:,:,:,1) * nhat(:,:,:,2) - t2(:,:,:,2) * nhat(:,:,:,1)
    f1 = sqrt( sum( t1 * t1, 4 ) )
    call invert( f1 )
    do i = 1, 3
        t1(:,:,:,i) = t1(:,:,:,i) * f1
    end do
    
    !ts
    t3(:,:,:,1) = sum(t0 * t1, 4)
    !td
    t3(:,:,:,2) = sum(t0 * t2, 4)
    !tn
    t3(:,:,:,3) = sum(t0 * nhat, 4)

    f1 = area
    do i = 1, 3
        nhat(:,:,:,i) = nhat(:,:,:,i) * f1
    end do
    
    call fieldio2d( '>', 'ts',  t3(:,:,:,1) )
    call fieldio2d( '>', 'td',  t3(:,:,:,2) )
    call fieldio2d( '>', 'tn',  t3(:,:,:,3) )

! end;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else if (ivols == 'yes' .and. ifn ==0 ) then
    call inivolstress
end if
if ( ifn == 0 )  return

if ( master ) write( 0, * ) 'Rupture initialization'

! I/O
t1 = 0.0
t2 = 0.0
t3 = 0.0

if ( friction == 'forced' ) then
! we still adopt slip-weakening here
    co = 0.0
    mus = 0.0
    mud = 0.0
    dc = 0.0
    call fieldio2d( '<>', 'co',  co          )
    call fieldio2d( '<>', 'mus', mus         )
    call fieldio2d( '<>', 'mud', mud         )
    call fieldio2d( '<>', 'dc',  dc          )
    
    
    if (nstage < 1) then
    	if (master) write(0,*) 'If forced rupture is used, please set nstage > 0'
    	stop
    end if
    
	tstage(1) = 0.
	rstage(1) = 0.
	do i = 2, nstage + 1
		! compute radius of far-point of each segment
		rstage(i) = rstage(i-1) + sizestage(i-1) 
		! compute analytical time at far-point of each segment
		tstage(i) = tstage(i-1) + sizestage(i-1)/vrupstage(i-1)		
	end do
	rupfrontstage = 1 ! we assume rupture start from 1st stage
end if

if ( friction == 'slipweakening' ) then
    co = 0.0
    mus = 0.0
    mud = 0.0
    dc = 0.0
    call fieldio2d( '<>', 'co',  co          )
    call fieldio2d( '<>', 'mus', mus         )
    call fieldio2d( '<>', 'mud', mud         )
    call fieldio2d( '<>', 'dc',  dc          )   
end if

if ( friction == 'rateandstate' ) then
    af = 0.0 
    bf = 0.0
    v0 = 0.0
    f0 = 0.0
    ll = 0.0
    fw = 0.0
    vw = 0.0
    psi = 0.0

    call fieldio2d( '<>', 'af',  af          )
    call fieldio2d( '<>', 'bf',  bf          )
    call fieldio2d( '<>', 'v0',  v0          )
    call fieldio2d( '<>', 'f0',  f0          )
    call fieldio2d( '<>', 'll',  ll          )
    call fieldio2d( '<>', 'fw',  fw          )
    call fieldio2d( '<>', 'vw',  vw          )
end if

if ( pcdep == 'yes' ) then
    lpc = 0.0
    call fieldio2d( '<>', 'lpc', lpc         )
end if

call fieldio2d( '<>', 's11', t1(:,:,:,1) )
call fieldio2d( '<>', 's22', t1(:,:,:,2) )
call fieldio2d( '<>', 's33', t1(:,:,:,3) )
call fieldio2d( '<>', 's23', t2(:,:,:,1) )
call fieldio2d( '<>', 's31', t2(:,:,:,2) )
call fieldio2d( '<>', 's12', t2(:,:,:,3) )
call fieldio2d( '<', 'ts',  t3(:,:,:,1) )
call fieldio2d( '<', 'td',  t3(:,:,:,2) )
call fieldio2d( '<', 'tn',  t3(:,:,:,3) )

!! Normal traction check
!i1 = maxloc( t3(:,:,:,3) )
!rr = t3(i1(1),i1(2),i1(3),3)
!i1(ifn) = irup
!i1 = i1 + nnoff
!if ( rr > 0.0 ) write( 0, * ) 'warning: positive normal traction: ', rr, i1

! Lock fault in PML region
i1 = i1pml + 1
i2 = i2pml - 1
call set_halo( co, 1e20, i1, i2 )

! Normal vectors
!i1 = i1core
!i2 = i2core
!i1(ifn) = irup
!i2(ifn) = irup
!call nodenormals( nhat, w1, dx, i1, i2, ifn )
!area = sign( 1, faultnormal ) * sqrt( sum( nhat * nhat, 4 ) )

f1 = area

! output area of each node on fault
call fieldio2d( '>', 'area', f1(:,:,:) )

call invert( f1 )
do i = 1, 3
    nhat(:,:,:,i) = nhat(:,:,:,i) * f1
end do

! output nhat with magnitude=1.0
call fieldio2d( '>', 'nhat1', nhat(:,:,:,1) )
call fieldio2d( '>', 'nhat2', nhat(:,:,:,2) )
call fieldio2d( '>', 'nhat3', nhat(:,:,:,3) )

! pore pressure
!pore0 = 0.
!pore = 0.
!call fieldio2d( '<>', 'pore0',  pore0(:,:,:) ) !here input initial pore pressure

do i = 1, 3
    j = modulo( i , 3 ) + 1
    k = modulo( i + 1, 3 ) + 1
    t0(:,:,:,i) = t0(:,:,:,i)   + &
    t1(:,:,:,i) * nhat(:,:,:,i) + &
    t2(:,:,:,j) * nhat(:,:,:,k) + &
    t2(:,:,:,k) * nhat(:,:,:,j)
end do


t2(:,:,:,1) = nhat(:,:,:,2) * slipvector(3) - nhat(:,:,:,3) * slipvector(2)
t2(:,:,:,2) = nhat(:,:,:,3) * slipvector(1) - nhat(:,:,:,1) * slipvector(3)
t2(:,:,:,3) = nhat(:,:,:,1) * slipvector(2) - nhat(:,:,:,2) * slipvector(1)
f1 = sqrt( sum( t2 * t2, 4 ) )
call invert( f1 )
do i = 1, 3
    t2(:,:,:,i) = t2(:,:,:,i) * f1
end do

t1(:,:,:,1) = t2(:,:,:,2) * nhat(:,:,:,3) - t2(:,:,:,3) * nhat(:,:,:,2)
t1(:,:,:,2) = t2(:,:,:,3) * nhat(:,:,:,1) - t2(:,:,:,1) * nhat(:,:,:,3)
t1(:,:,:,3) = t2(:,:,:,1) * nhat(:,:,:,2) - t2(:,:,:,2) * nhat(:,:,:,1)
f1 = sqrt( sum( t1 * t1, 4 ) )
call invert( f1 )
do i = 1, 3
    t1(:,:,:,i) = t1(:,:,:,i) * f1
end do

do i = 1, 3
    t0(:,:,:,i) = t0(:,:,:,i) + &
    t3(:,:,:,1) * t1(:,:,:,i) + &
    t3(:,:,:,2) * t2(:,:,:,i) + &
    t3(:,:,:,3) * nhat(:,:,:,i)
end do

f1 = sum( t0 * nhat, 4 ) 
do i = 1, 3
    ts0(:,:,:,i) = t0(:,:,:,i) - f1 * nhat(:,:,:,i)
end do

if ( ( rcrit > 0.0 .and. (vrup > 0.0 .or. nstage > 0) ) .or. ( rnucl > 0.0 ) ) then
    xhypo = 0.0
    xi = ihypo - nnoff
    i1 = floor( xi )
    if ( all( i1 >= 1 .and. i1 < nm ) ) then
        do l = i1(3), i1(3)+1
        do k = i1(2), i1(2)+1
        do j = i1(1), i1(1)+1
            w = (1.0-abs(xi(1)-j)) * (1.0-abs(xi(2)-k)) * (1.0-abs(xi(3)-l))
            do i = 1, 3
                xhypo(i) = xhypo(i) + w * w1(j,k,l,i)
            end do
        end do
        end do
        end do
    end if
    call rbroadcast1( xhypo, ip2root )
    do i = 1, 3
        select case( ifn )
        case ( 1 ); t2(1,:,:,i) = w1(irup,:,:,i) - xhypo(i)
        case ( 2 ); t2(:,1,:,i) = w1(:,irup,:,i) - xhypo(i)
        case ( 3 ); t2(:,:,1,i) = w1(:,:,irup,i) - xhypo(i)
        end select
    end do
    rhypo = sqrt( sum( t2 * t2, 4 ) )
end if

! Resample mu on to fault plane nodes for moment calculatioin
select case( ifn )
case ( 1 ); lamf(1,:,:) = lam(irup,:,:); muf(1,:,:) = mu(irup,:,:)
case ( 2 ); lamf(:,1,:) = lam(:,irup,:); muf(:,1,:) = mu(:,irup,:)
case ( 3 ); lamf(:,:,1) = lam(:,:,irup); muf(:,:,1) = mu(:,:,irup)
end select
call invert( lamf )
call invert( muf )
j = nm(1) - 1
k = nm(2) - 1
l = nm(3) - 1
if ( ifn /= 1 ) lamf(2:j,:,:) = 0.5 * (lamf(2:j,:,:) + lamf(1:j-1,:,:))
if ( ifn /= 2 ) lamf(:,2:k,:) = 0.5 * (lamf(:,2:k,:) + lamf(:,1:k-1,:))
if ( ifn /= 3 ) lamf(:,:,2:l) = 0.5 * (lamf(:,:,2:l) + lamf(:,:,1:l-1))
if ( ifn /= 1 ) muf(2:j,:,:) = 0.5 * (muf(2:j,:,:) + muf(1:j-1,:,:))
if ( ifn /= 2 ) muf(:,2:k,:) = 0.5 * (muf(:,2:k,:) + muf(:,1:k-1,:))
if ( ifn /= 3 ) muf(:,:,2:l) = 0.5 * (muf(:,:,2:l) + muf(:,:,1:l-1))
call invert( muf )

! Initial state, can be overwritten by read_checkpoint
psv   =  0.0
trup  =  1e9
tarr  =  0.0
efric =  0.0

! Halos
call scalar_swap_halo( co,    nhalo )
call scalar_swap_halo( area,  nhalo )
call scalar_swap_halo( rhypo, nhalo )
call vector_swap_halo( nhat,  nhalo )
call vector_swap_halo( t0,    nhalo )
call vector_swap_halo( ts0,   nhalo )
!call scalar_swap_halo( pore0, nhalo )

if ( friction == 'forced' ) then
    call scalar_swap_halo( mus,   nhalo )
    call scalar_swap_halo( mud,   nhalo )
    call scalar_swap_halo( dc,    nhalo )
end if

if ( friction == 'slipweakening' ) then
    call scalar_swap_halo( mus,   nhalo )
    call scalar_swap_halo( mud,   nhalo )
    call scalar_swap_halo( dc,    nhalo )
end if

if ( friction == 'rateandstate' ) then
    call scalar_swap_halo( af,   nhalo )
    call scalar_swap_halo( bf,   nhalo )
    call scalar_swap_halo( v0,   nhalo )
    call scalar_swap_halo( f0,   nhalo )
    call scalar_swap_halo( ll,   nhalo )    
    call scalar_swap_halo( fw,   nhalo )
    call scalar_swap_halo( vw,   nhalo )
end if

if ( pcdep == 'yes' ) call scalar_swap_halo( lpc,  nhalo )

end subroutine

!------------------------------------------------------------------------------!

! Rupture boundary condition
subroutine rupture
use m_globals
use m_collective
use m_bc
use m_util
use m_fieldio
use m_stats
integer ::  niter, nmax, i1(3), i2(3), i, j, k, l,         &
           j1, k1, l1, j2, k2, l2, j3, k3, l3, j4, k4, l4
real :: tol

if ( ifn == 0 ) return
if ( verb ) write( 0, * ) 'Rupture'

if ( friction == 'rateandstate' ) then
    nmax = 30
    tol = 1.e-3
end if

! Indices
i1 = 1
i2 = nm
i1(ifn) = irup
i2(ifn) = irup
j1 = i1(1); j2 = i2(1)
k1 = i1(2); k2 = i2(2)
l1 = i1(3); l2 = i2(3)
i1(ifn) = irup + 1
i2(ifn) = irup + 1
j3 = i1(1); j4 = i2(1)
k3 = i1(2); k4 = i2(2)
l3 = i1(3); l4 = i2(3)

if ( it == 1) then
    tp = 0.
    tptmp = 0.
    
    if ( tm < tmnucl ) then
        f2 = exp( (tm - tmnucl)**2/(tm*(tm -2*tmnucl)) )
    else
        f2 = 1.
    end if

    f1 = 0.
    where( rhypo < rnucl ) f1 = exp( rhypo**2/(rhypo**2-rnucl**2) ) * f2 * delts

    do i = 1, 3
        tptmp(:,:,:,i) = ts0(:,:,:,i) * f1
    end do
else
    if ( tm < tmnucl ) then
        f2 = exp( (tm - tmnucl)**2/(tm*(tm -2*tmnucl)) )
    else
        f2 = 1.
    end if

    f1 = 0.
    where( rhypo < rnucl ) f1 = exp( rhypo**2/(rhypo**2-rnucl**2) ) * f2 * delts

    do i = 1, 3
        tp(:,:,:,i) = ts0(:,:,:,i) * f1
    end do
end if

if ( it == 1 ) then
    w2 = vv
    if ( maxval(abs(w2(:,:,:,1))) > 1.e-20 ) then
        vv(j3:j4,k3:k4,l3:l4,1) = w2(j3:j4,k3:k4,l3:l4,1)
        vv(j3:j4,k3:k4,l3:l4,2) = &
           - nhat(:,:,:,1)*nhat(:,:,:,2)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1)) * w2(j3:j4,k3:k4,l3:l4,1)
        vv(j3:j4,k3:k4,l3:l4,3) = &
           - nhat(:,:,:,1)*nhat(:,:,:,3)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1)) * w2(j3:j4,k3:k4,l3:l4,1)

        vv(j1:j2,k1:k2,l1:l2,1) = w2(j1:j2,k1:k2,l1:l2,1)
        vv(j1:j2,k1:k2,l1:l2,2) = &
           - nhat(:,:,:,1)*nhat(:,:,:,2)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1)) * w2(j1:j2,k1:k2,l1:l2,1)
        vv(j1:j2,k1:k2,l1:l2,3) = &
           - nhat(:,:,:,1)*nhat(:,:,:,3)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1)) * w2(j1:j2,k1:k2,l1:l2,1)
    end if
    if ( maxval(abs(w2(:,:,:,2))) > 1.e-20 ) then
        vv(j3:j4,k3:k4,l3:l4,1) = vv(j3:j4,k3:k4,l3:l4,1) &
           - nhat(:,:,:,1)*nhat(:,:,:,2)/(1 - nhat(:,:,:,2)*nhat(:,:,:,2)) * w2(j3:j4,k3:k4,l3:l4,2)
        vv(j3:j4,k3:k4,l3:l4,2) = vv(j3:j4,k3:k4,l3:l4,2) + w2(j3:j4,k3:k4,l3:l4,2)
        vv(j3:j4,k3:k4,l3:l4,3) = vv(j3:j4,k3:k4,l3:l4,2) &
           - nhat(:,:,:,2)*nhat(:,:,:,3)/(1 - nhat(:,:,:,2)*nhat(:,:,:,2)) * w2(j3:j4,k3:k4,l3:l4,2)

        vv(j1:j2,k1:k2,l1:l2,1) = vv(j1:j2,k1:k2,l1:l2,1) &
           - nhat(:,:,:,1)*nhat(:,:,:,2)/(1 - nhat(:,:,:,2)*nhat(:,:,:,2)) * w2(j1:j2,k1:k2,l1:l2,2)
        vv(j1:j2,k1:k2,l1:l2,2) = vv(j1:j2,k1:k2,l1:l2,2) + w2(j1:j2,k1:k2,l1:l2,2)
        vv(j1:j2,k1:k2,l1:l2,3) = vv(j1:j2,k1:k2,l1:l2,2) &
           - nhat(:,:,:,2)*nhat(:,:,:,3)/(1 - nhat(:,:,:,2)*nhat(:,:,:,2)) * w2(j1:j2,k1:k2,l1:l2,2)
    end if
    if ( maxval(abs(w2(:,:,:,3))) > 1.e-20 ) then
        vv(j3:j4,k3:k4,l3:l4,1) = vv(j3:j4,k3:k4,l3:l4,1) &
           - nhat(:,:,:,1)*nhat(:,:,:,3)/(1 - nhat(:,:,:,3)*nhat(:,:,:,3)) * w2(j3:j4,k3:k4,l3:l4,3)
        vv(j3:j4,k3:k4,l3:l4,2) = vv(j3:j4,k3:k4,l3:l4,2) &
           - nhat(:,:,:,2)*nhat(:,:,:,3)/(1 - nhat(:,:,:,3)*nhat(:,:,:,3)) * w2(j3:j4,k3:k4,l3:l4,3)
        vv(j3:j4,k3:k4,l3:l4,3) = vv(j3:j4,k3:k4,l3:l4,3) + w2(j3:j4,k3:k4,l3:l4,3)

        vv(j1:j2,k1:k2,l1:l2,1) = vv(j1:j2,k1:k2,l1:l2,1) &
           - nhat(:,:,:,1)*nhat(:,:,:,3)/(1 - nhat(:,:,:,3)*nhat(:,:,:,3)) * w2(j1:j2,k1:k2,l1:l2,3)
        vv(j1:j2,k1:k2,l1:l2,2) = vv(j1:j2,k1:k2,l1:l2,2) &
           - nhat(:,:,:,2)*nhat(:,:,:,3)/(1 - nhat(:,:,:,3)*nhat(:,:,:,3)) * w2(j1:j2,k1:k2,l1:l2,3)
        vv(j1:j2,k1:k2,l1:l2,3) = vv(j1:j2,k1:k2,l1:l2,3) + w2(j1:j2,k1:k2,l1:l2,3)
    end if

!    do i = 2, 3
!    vv(j3:j4,k3:k4,l3:l4,i) = -vv(j3:j4,k3:k4,l3:l4,1) * nhat(:,:,:,i) * &
!                              nhat(:,:,:,1)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1))
!    vv(j1:j2,k1:k2,l1:l2,i) = -vv(j1:j2,k1:k2,l1:l2,1) * nhat(:,:,:,i) * &
!                              nhat(:,:,:,1)/(1 - nhat(:,:,:,1)*nhat(:,:,:,1))
!    end do 
end if
!---------------------------------------------------------------

! Trial traction for zero velocity and zero displacement
! parameter 'c' defined in (5c) in front of normal traction
f3 = area * ( mr(j1:j2,k1:k2,l1:l2) + mr(j3:j4,k3:k4,l3:l4) )
f1 = dt * dt * f3; 
call invert( f1 )

do i = 1, 3
    	t1(:,:,:,i) = t0(:,:,:,i) + tp(:,:,:,i) + f1 * dt * &
        	( vv(j3:j4,k3:k4,l3:l4,i) &
        	- vv(j1:j2,k1:k2,l1:l2,i) &	
        	+ w1(j3:j4,k3:k4,l3:l4,i) * mr(j3:j4,k3:k4,l3:l4) * dt &
        	- w1(j1:j2,k1:k2,l1:l2,i) * mr(j1:j2,k1:k2,l1:l2) * dt )
    	t2(:,:,:,i) = t1(:,:,:,i) + f1 * &
        	( uu(j3:j4,k3:k4,l3:l4,i) - uu(j1:j2,k1:k2,l1:l2,i) )
end do


! Normal and Shear traction  [ZS]
tn = sum( t2 * nhat, 4 )
!if ( faultopening == 1 ) tn = min( 0.0, tn )
tn = min( 0.0, tn )
! Pore pressure comoputation
!pore = pore0 + skepb * (tn - sum(( t0 + tp ) * nhat, 4))
!tn = tn - pore

f1 = sum( t1 * nhat, 4 ) 
do i = 1, 3
    t3(:,:,:,i) = t1(:,:,:,i) - f1 * nhat(:,:,:,i)
end do
ts = sqrt( sum( t3 * t3, 4 ) )



! slip velocity from the previous step
do i = 1, 3     
    t2(:,:,:,i) = vv(j3:j4,k3:k4,l3:l4,i) - vv(j1:j2,k1:k2,l1:l2,i)
end do

if ( friction == 'rateandstate' ) then

svold = sqrt( sum( t2 * t2, 4 ) )
! call scalar_swap_halo( sv0,   nhalo )

if ( it == 1 )  then
   
    psi = af * sngl(dlog( 2.0 * v0 / svold * sinh( dble(ts / (-tn) / af )) ))
!!!!!!!   nucleation to reduce initial psi
    if (psidelts > 0. .and. psidelts < 1.) then    
    if ( tm < tmnucl ) then
        f2 = exp( (tm - tmnucl)**2/(tm*(tm -2*tmnucl)) )
    else
        f2 = 1.
    end if

    f1 = 1.
    where( rhypo < rnucl ) f1 = 1. - exp( rhypo**2/(rhypo**2-rnucl**2) ) * f2 * psidelts
    psi = psi * f1
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!
    f1 = ts
    sv0 = svold
   
    if ( pcdep == 'yes' ) then 
        tnpc = tn 
        tnold = tn
    end if  
else
    where ( svold < sv0 ) svold = sv0
   
    if ( pcdep == 'yes' ) then 
        tnpc = tnpc + dt * svold / lpc *( tnold - tnpc )
        tnold = tn
    end if

    svtrl = ts * f3 * dt  

    f1 = f0 - ( bf - af ) * log( svold / v0 ) 
    
    where ( svold < 1.e-5 )
        f2 = f1
    elsewhere
        f2 = fw + ( f1 - fw ) / ( 1.0 + (svold / vw )**8.0 )**(1.0/8.0)  
    end where
    
    if (any (f2 < 0)) then
        i1 = minloc(f2)
        i1(ifn) = -1
        if( i1(1) /= 1 .and. i1(2) /= 1 .and. i1(3) /= 1 .and. &
            i1(1)/=nm(1) .and. i1(2)/=nm(2) .and. i1(3)/=nm(3)) then
                write(0,*) 'f2 is negative at', i1
                !if (svold(i1(1),i1(2),i1(3)) < 1e-20) write(0,*) 'Reason: sliprate is assigned as 0.0'
                stop 'error in af or vini (af to large) or (vini to small)'
        end if !if svold == 0 at halo, not sure this will genrate NaN data (bug not fully solved)
    end if   
 
    f2 = af * ( log( 2 * v0 ) + sngl(dlog( sinh(dble(f2/af))) ) - log( svold ) ) 

    f4 = exp( - svold * dt / ll )
    
    psi = psi * f4 + ( 1.0 - f4 ) * f2 

    f1 = 0.5 * svtrl / v0 * exp( psi / af ) 
    where ( f1 > 1.e6 )
        f1 = af * ( log(svtrl) - log(v0) ) + psi         
    elsewhere
        f1 = af * log( f1 + sqrt( f1*f1 + 1.0 ) )
    end where

    f2 = 0.5 * svold / v0 * exp( psi / af )  
    where ( f2 > 1.e6 )
        f2 = af * ( log(svold) - log(v0) ) + psi         
    elsewhere
        f2 = af * log( f2 + sqrt( f2*f2 + 1.0 ) )
    end where
    
    f1 = 0.5 * ( f1 + f2 ) 

    if ( pcdep == 'yes' ) then
        f4 = dt * f3 * (-tnpc)
    else
        f4 = dt * f3 * (-tn)
    end if
     
    f2 = psi + 0.82   
    where (f1 > f2)  f1 = f2

    i1 = i1node
    i2 = i2node
    i1(ifn) = 1
    i2(ifn) = 1
    delf = 1.0   

    do l = i1(3), i2(3)
    do k = i1(2), i2(2)
    do j = i1(1), i2(1)

    niter = 0
    
    do while ( abs(delf(j,k,l)) > tol .and. niter <= nmax ) 

        fun(j,k,l) = svtrl(j,k,l) - f4(j,k,l) * f1(j,k,l) - v0(j,k,l) * & 
             ( exp((f1(j,k,l)-psi(j,k,l))/af(j,k,l)) - exp(-(f1(j,k,l)+psi(j,k,l))/af(j,k,l)) )
        dfun(j,k,l) = -f4(j,k,l) - ( exp((f1(j,k,l)-psi(j,k,l))/af(j,k,l)) + &
             exp(-(f1(j,k,l)+psi(j,k,l))/af(j,k,l)) ) *V0(j,k,l) / af(j,k,l)

        delf(j,k,l) = fun(j,k,l) / dfun(j,k,l)
        f1(j,k,l) = f1(j,k,l) - delf(j,k,l)

        niter = niter + 1
!        call set_halo( delf, -1.0, i1core, i2core) 
    end do

    end do
    end do
    end do
    !------------------------------------------  
    
    if ( pcdep == 'yes' ) then 
        f1 = -min( 0.0, tnpc ) * f1
    else
        f1 = -min( 0.0, tn ) * f1
    end if

! Shear traction bounded by friction (need this?)
!    f2 = 1.0
!    where ( ts > f1 ) f2 = f1 / ts
!    do i = 1, 3
!        t3(:,:,:,i) = f2 * t3(:,:,:,i)
!    end do
!    ts = min( ts, f1 )

    f2 = f1 / ts
    do i = 1, 3
        t3(:,:,:,i) = f2 * t3(:,:,:,i)
    end do
    ts = f1
    
end if

end if



if ( friction == 'forced' ) then
	! forced rupture with variable rupture velocity described in each segment
if ( it > 1) then
	f1 = mud
    where ( sl < dc ) f1 = f1 + (1.0 - sl / dc) * (mus - mud)
    f1 = -min( 0.0, tn ) * f1 + co
    
    ! inside the rupture zone
    if ( rcrit > 0.0 .and. vrupstage(1) > 0. .and. sizestage(1) > 0.) then
    
    	! compute analytical forced rupture radius
    	if ( rupfrontstage <=nstage) then
    		rupfront = rstage(rupfrontstage) + ( tm - tstage(rupfrontstage)) * vrupstage(rupfrontstage)
    	else
    		rupfront = rstage(rupfrontstage) + ( tm - tstage(rupfrontstage)) * vrupstage(rupfrontstage-1)
    	end if
    	
    	if (rupfront > rstage(rupfrontstage+1) .and. rupfrontstage <= nstage) rupfrontstage = rupfrontstage + 1 
    	
    	f2 = 1.0
    	if(rrelax > 0.0) f2 = min( (rupfront - rhypo) / rrelax, 1.0 )
    	
    	
    	f2 = (1.0 - f2) * ts + f2 * (-tn * mud + co)
        where ( rhypo < min( rcrit, rupfront ) .and. f2 < f1 ) f1 = f2
        
        !write(0,*) f1(:,1,2)
    end if
    
    ! Shear traction bounded by friction
    f2 = 1.0
    where ( ts > f1 ) f2 = f1 / ts
    do i = 1, 3
        t3(:,:,:,i) = f2 * t3(:,:,:,i)
    end do
    ts = min( ts, f1 )
    
end if
end if

if ( friction == 'slipweakening' ) then
   ! Slip-weakening friction law
if ( it > 1) then    
    f1 = mud
    where ( sl < dc ) f1 = f1 + (1.0 - sl / dc) * (mus - mud)
    f1 = -min( 0.0, tn ) * f1 + co

    ! Nucleation
    if ( rcrit > 0.0 .and. vrup > 0.0 ) then
    
        f2 = 1.0
        if ( trelax > 0.0 ) f2 = min( (tm - rhypo / vrup) / trelax, 1.0 )
        f2 = (1.0 - f2) * ts + f2 * (-tn * mud + co)
        where ( rhypo < min( rcrit, tm * vrup ) .and. f2 < f1 ) f1 = f2
        
        
! user designed
!        f2 = 1.0
!        if ( trelax > 0.0 ) then
!        	f2 = min( (tm - rhypo / vrup) / trelax, 1.0 )
!	  	f2 = (1.0 - f2) * ts + f2 * (-tn * mud + co) 
!	  	where ( rhypo < min( rcrit, tm * vrup ) .and. f2 < f1 ) f1 = f2       	
!        else if (tslope > 0.0 ) then
!            f2 = max(-tn*mud+co,-tn*mud+co+tslope*1.e3*(rhypo-tm*vrup))
!            where ( rhypo < min( rcrit, tm * vrup )) f1 = f2
!        end if
        
    end if

    ! Shear traction bounded by friction
    f2 = 1.0
    where ( ts > f1 ) f2 = f1 / ts
    do i = 1, 3
        t3(:,:,:,i) = f2 * t3(:,:,:,i)
    end do
    ts = min( ts, f1 )
end if

end if

! Total traction
!tn = tn + pore ! tn here is true normal traction
do i = 1, 3
    t1(:,:,:,i) = t3(:,:,:,i) + tn * nhat(:,:,:,i) 
end do


if (friction == 'frictionless' .and. intype == 'sa') then
	f3 = area * ( mr(j1:j2,k1:k2,l1:l2) + mr(j3:j4,k3:k4,l3:l4) )
	call invert( f3 )
	sla = 0.0
	call fieldio( '<', 'sa1', sla(:,:,:,1) )
	call fieldio( '<', 'sa2', sla(:,:,:,2) )
	call fieldio( '<', 'sa3', sla(:,:,:,3) )
	do i = 1, 3
    		t1(:,:,:,i) = t0(:,:,:,i) + tp(:,:,:,i) + f3 * &
        		( -sla(:,:,:,i) &
        		+ w1(j3:j4,k3:k4,l3:l4,i) * mr(j3:j4,k3:k4,l3:l4)  &
        		- w1(j1:j2,k1:k2,l1:l2,i) * mr(j1:j2,k1:k2,l1:l2)  )
	end do
	
	tn = sum( t1 * nhat, 4 )
	tn = min( 0.0, tn )
	f1 = sum( t1 * nhat, 4 ) 
	do i = 1, 3
    		t3(:,:,:,i) = t1(:,:,:,i) - f1 * nhat(:,:,:,i)
	end do
	ts = sqrt( sum( t3 * t3, 4 ) )


elseif(friction == 'frictionless' .and. intype == 'sv') then
	f3 = area * ( mr(j1:j2,k1:k2,l1:l2) + mr(j3:j4,k3:k4,l3:l4) )
	call invert( f3 )
	if(it == 1) then
		slv = 0.0
		call fieldio( '<', 'sv1', slv(:,:,:,1) )
		call fieldio( '<', 'sv2', slv(:,:,:,2) )
		call fieldio( '<', 'sv3', slv(:,:,:,3) )
		slvold = slv
	endif
	
	if (it < nt) then
		it = it + 1
		slv = 0.0
		call fieldio( '<', 'sv1', slv(:,:,:,1) )
		call fieldio( '<', 'sv2', slv(:,:,:,2) )
		call fieldio( '<', 'sv3', slv(:,:,:,3) )
		it = it - 1
	elseif (it == nt) then
		slv = slvold
	endif
	
	do i = 1, 3
		sla(:,:,:,i) = (slv(:,:,:,i)-slvold(:,:,:,i))/dt
    		t1(:,:,:,i) = t0(:,:,:,i) + tp(:,:,:,i) + f3 * &
        		( - sla(:,:,:,i) &
        		+ w1(j3:j4,k3:k4,l3:l4,i) * mr(j3:j4,k3:k4,l3:l4)  &
        		- w1(j1:j2,k1:k2,l1:l2,i) * mr(j1:j2,k1:k2,l1:l2)  )
	end do
	slvold = slv
	
	tn = sum( t1 * nhat, 4 )
	tn = min( 0.0, tn )
	f1 = sum( t1 * nhat, 4 ) 
	do i = 1, 3
    		t3(:,:,:,i) = t1(:,:,:,i) - f1 * nhat(:,:,:,i)
	end do
	ts = sqrt( sum( t3 * t3, 4 ) )
elseif(friction == 'frictionless' .and. intype == 'fsu' .and. tm < period+1.5*dt) then
!elseif(friction == 'frictionless' .and. intype == 'fsu' ) then
	if (it == 1) then
		slu = 0.0
		call fieldio( '<', 'su1', slu(:,:,:,1) )
		call fieldio( '<', 'su2', slu(:,:,:,2) )
		call fieldio( '<', 'su3', slu(:,:,:,3) )
		slvold = 0.
		t1 = t0 + tp
	else
	
	f3 = area * ( mr(j1:j2,k1:k2,l1:l2) + mr(j3:j4,k3:k4,l3:l4) )
	call invert( f3 )
	sla = 0.0

	  do i = 1, 3
		slv(:,:,:,i) = slu(:,:,:,i) * time_function('liu',tm-dt/2.,dt,period)
		sla(:,:,:,i) = (slv(:,:,:,i)-slvold(:,:,:,i))/dt
    		t1(:,:,:,i) = t0(:,:,:,i) + tp(:,:,:,i) + f3 * &
        		( - sla(:,:,:,i) &
        		+ w1(j3:j4,k3:k4,l3:l4,i) * mr(j3:j4,k3:k4,l3:l4)  &
        		- w1(j1:j2,k1:k2,l1:l2,i) * mr(j1:j2,k1:k2,l1:l2)  )
	  end do
	  slvold = slv
	end if
	
!	write(0,*) sla(32,32,1,1)
	tn = sum( t1 * nhat, 4 )
	tn = min( 0.0, tn )
	f1 = sum( t1 * nhat, 4 ) 
	do i = 1, 3
    		t3(:,:,:,i) = t1(:,:,:,i) - f1 * nhat(:,:,:,i)
	end do
	ts = sqrt( sum( t3 * t3, 4 ) )
end if

! Update acceleration (n) step
do i = 1, 3
    f2 = area * ( t1(:,:,:,i) - t0(:,:,:,i) - tp(:,:,:,i) )
    w1(j1:j2,k1:k2,l1:l2,i) = w1(j1:j2,k1:k2,l1:l2,i) + f2
    w1(j3:j4,k3:k4,l3:l4,i) = w1(j3:j4,k3:k4,l3:l4,i) - f2
end do

call vector_bc( w1, bc1, bc2, i1bc, i2bc )

! Output
if (it==1) then
    t1 = t1 + tptmp
    t3 = t3 + tptmp
    ts = sqrt(sum(t3*t3,4))
end if
call fieldio( '>', 't1',  t1(:,:,:,1) )
call fieldio( '>', 't2',  t1(:,:,:,2) )
call fieldio( '>', 't3',  t1(:,:,:,3) )
call fieldio( '>', 'ts1', t3(:,:,:,1) )
call fieldio( '>', 'ts2', t3(:,:,:,2) )
call fieldio( '>', 'ts3', t3(:,:,:,3) )
call fieldio( '>', 'tsm', ts          )
call fieldio( '>', 'tnm', tn          ) 
!call fieldio( '>', 'pore',  pore      )
call fieldio( '>', 'fr',  f1          )
if (it==1) then
    t1 = t1 - tptmp
    t3 = t3 - tptmp
    ts = sqrt(sum(t3*t3,4))
end if

if ( friction == 'rateandstate' ) then
    call fieldio( '>', 'psi', psi         )
end if

call set_halo( ts,      -1.0, i1core, i2core ); tsmax = maxval( ts ) 
call set_halo( tn,  huge(dt), i1core, i2core ); tnmin = minval( tn )
call set_halo( tn, -huge(dt), i1core, i2core ); tnmax = maxval( tn )
call set_halo( tn,       0.0, i1core, i2core )
!call set_halo( pore,       0.0, i1core, i2core )


! Slip acceleration
do i = 1, 3
    t2(:,:,:,i) = &
        w1(j3:j4,k3:k4,l3:l4,i) * mr(j3:j4,k3:k4,l3:l4) - &
        w1(j1:j2,k1:k2,l1:l2,i) * mr(j1:j2,k1:k2,l1:l2)
end do
f2 = sqrt( sum( t2 * t2, 4 ) )

call fieldio( '>', 'sa1', t2(:,:,:,1) )
call fieldio( '>', 'sa2', t2(:,:,:,2) )
call fieldio( '>', 'sa3', t2(:,:,:,3) )
call fieldio( '>', 'sam', f2          )
call set_halo( f2, -1.0, i1core, i2core )
samax = maxval( f2 )

! Friction + fracture energy
!sv(n+1/2) = sv(n-1/2) + dt * sa(n)
! ( sv(n+1/2) + sv(n-1/2) )/2 = sv(n-1/2) + 0.5 * dt * sa(n)
select case( ifn )
case ( 1 )
t2(1,:,:,:) = vv(irup+1,k3:k4,l3:l4,:) - vv(irup,k1:k2,l1:l2,:) + &
                 0.5 * dt * t2(1,k3:k4,l3:l4,:)
case ( 2 )
t2(:,1,:,:) = vv(j3:j4,irup+1,l3:l4,:) - vv(j1:j2,irup,l1:l2,:) + &
                 0.5 * dt * t2(j3:j4,1,l3:l4,:)
case ( 3 )
t2(:,:,1,:) = vv(j3:j4,k3:k4,irup+1,:) - vv(j1:j2,k1:k2,irup,:) + &
                 0.5 * dt * t2(j3:j4,k3:k4,1,:)
end select
f2 = sum( t1 * t2, 4 ) * area
call set_halo( f2, 0.0, i1core, i2core )
efric = efric + dt * sum( f2 )

! Strain energy
t2 = uu(j3:j4,k3:k4,l3:l4,:) - uu(j1:j2,k1:k2,l1:l2,:)
f2 = sum( (t0 + tptmp + t1) * t2, 4 ) * area
call set_halo( f2, 0.0, i1core, i2core )
estrain = 0.5 * sum( f2 )

eradiat = estrain - efric

! Stress drop 
t2 = uu(j3:j4,k3:k4,l3:l4,:) - uu(j1:j2,k1:k2,l1:l2,:)
f2 = sum((t0 + tp - t1) * t2,4) * area
call set_halo( f2, 0.0, i1core, i2core)
strdropint = sum( f2 ) 
f2 = sqrt(sum(t2*t2,4))*area
call set_halo( f2, 0.0, i1core, i2core)
slipint = sum( f2 )

! Moment (negelcts opening lambda contribution)
f2 = muf * area * sqrt( sum( t2 * t2, 4 ) )
call set_halo( f2, 0.0, i1core, i2core )
moment = sum( f2 )

end subroutine

end module

