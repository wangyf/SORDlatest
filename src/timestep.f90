! Time integration
module m_timestep
implicit none
contains

subroutine timestep
use m_globals
use m_util
use m_fieldio
use m_stats

! Status
if ( master ) then
    if ( verb ) then
        write( 0, * ) 'Time step', it
    else
        write( 0, '(a)', advance='no' ) '.'
        if ( modulo( it, 50 ) == 0 .or. it == nt ) write( 0, '(i6)' ) it
    end if
end if

! Save previous slip velocity
if ( ifn /= 0 ) then
    select case( ifn )
    case( 1 ); t2(1,:,:,:) = vv(irup+1,:,:,:) - vv(irup,:,:,:)
    case( 2 ); t2(:,1,:,:) = vv(:,irup+1,:,:) - vv(:,irup,:,:)
    case( 3 ); t2(:,:,1,:) = vv(:,:,irup+1,:) - vv(:,:,irup,:)
    end select
    f2 = sqrt( sum( t2 * t2, 4 ) )
end if

! Velocity time integration
tm = tm0 + dt * ( it - 1 ) - dt * 0.5
vv = vv + dt * w1
call fieldio( '<>', 'v1', vv(:,:,:,1) )
call fieldio( '<>', 'v2', vv(:,:,:,2) )
call fieldio( '<>', 'v3', vv(:,:,:,3) )
if ( modulo( it, itstats ) == 0 ) then
    call vector_norm( s1, vv, i1core, i2core, (/ 1, 1, 1 /) )
    call set_halo( s1, -1.0, i1core, i2core )
    vmax = maxval( s1 )
end if
call fieldio( '>', 'vm2', s1  )

! Displacement time integration
tm = tm0 + dt * ( it - 1 )
uu = uu + dt * vv
call fieldio( '<>', 'u1', uu(:,:,:,1) )
call fieldio( '<>', 'u2', uu(:,:,:,2) )
call fieldio( '<>', 'u3', uu(:,:,:,3) )
if ( modulo( it, itstats ) == 0 ) then
    call vector_norm( s1, uu, i1core, i2core, (/ 1, 1, 1 /) )
    call set_halo( s1, -1.0, i1core, i2core )
    umax = maxval( s1 )
end if
call fieldio( '>', 'um2', s1  )
! note: 
!it = 1 tm = 0
!it = 2 tm = 0
!it = 3 tm = dt ...



! Fault time integration
if ( ifn /= 0 ) then
    select case( ifn )
    case( 1 ); t1(1,:,:,:) = vv(irup+1,:,:,:) - vv(irup,:,:,:)
    case( 2 ); t1(:,1,:,:) = vv(:,irup+1,:,:) - vv(:,irup,:,:)
    case( 3 ); t1(:,:,1,:) = vv(:,:,irup+1,:) - vv(:,:,irup,:)
    end select
    
    ! compute moment rate
    ! \dot M_xx,  \dot M_yy, \dot M_zz
    mort1(:,:,:,1) = lamf(:,:,:)*area(:,:,:)*sum(t1*nhat,4) + &
                     2*muf(:,:,:)*area(:,:,:)*t1(:,:,:,1)*nhat(:,:,:,1)
    mort1(:,:,:,2) = lamf(:,:,:)*area(:,:,:)*sum(t1*nhat,4) + &
                     2*muf(:,:,:)*area(:,:,:)*t1(:,:,:,2)*nhat(:,:,:,2)  
    mort1(:,:,:,3) = lamf(:,:,:)*area(:,:,:)*sum(t1*nhat,4) + &
                     2*muf(:,:,:)*area(:,:,:)*t1(:,:,:,3)*nhat(:,:,:,3) 
    ! \dot M_yz, \dot M_xz, \dot M_xy
    mort2(:,:,:,1) = muf(:,:,:)*area(:,:,:)* &
                     (t1(:,:,:,2)*nhat(:,:,:,3)+t1(:,:,:,3)*nhat(:,:,:,2))
    mort2(:,:,:,2) = muf(:,:,:)*area(:,:,:)* &
                     (t1(:,:,:,1)*nhat(:,:,:,3)+t1(:,:,:,3)*nhat(:,:,:,1))
    mort2(:,:,:,3) = muf(:,:,:)*area(:,:,:)* &
                     (t1(:,:,:,2)*nhat(:,:,:,1)+t1(:,:,:,1)*nhat(:,:,:,2))                
    ! finish moment rate tensor
    
    f1 = sqrt( sum( t1 * t1, 4 ) )
    sl = sl + dt * f1
    psv = max( psv, f1 )
    if ( svtol > 0.0 ) then
        where ( f1 >= svtol .and. trup > 1e8 )
            trup = tm - dt * ( 0.5 + (svtol - f1) / (f2 - f1) )
        end where
        where ( f1 >= svtol )
            tarr = 1e9
        end where
        where ( f1 < svtol .and. f2 >= svtol )
            tarr = tm - dt * ( 0.5 + (svtol - f1) / (f2 - f1) )
        end where
    end if
    select case( ifn )
    case( 1 ); t2(1,:,:,:) = uu(irup+1,:,:,:) - uu(irup,:,:,:)
    case( 2 ); t2(:,1,:,:) = uu(:,irup+1,:,:) - uu(:,irup,:,:)
    case( 3 ); t2(:,:,1,:) = uu(:,:,irup+1,:) - uu(:,:,irup,:)
    end select
    f2 = sqrt( sum( t2 * t2, 4 ) )
    
!    if(master) write(0,*) tm,f2(32,32,1),slu(32,32,1,1)
    
    call fieldio( '>', 'sv1',  t1(:,:,:,1) )
    call fieldio( '>', 'sv2',  t1(:,:,:,2) )
    call fieldio( '>', 'sv3',  t1(:,:,:,3) )
    call fieldio( '>', 'mr11',  mort1(:,:,:,1) )
    call fieldio( '>', 'mr22',  mort1(:,:,:,2) )
    call fieldio( '>', 'mr33',  mort1(:,:,:,3) ) 
    call fieldio( '>', 'mr23',  mort2(:,:,:,1) )
    call fieldio( '>', 'mr31',  mort2(:,:,:,2) )
    call fieldio( '>', 'mr12',  mort2(:,:,:,3) )   
    call fieldio( '>', 'svm',  f1          )
    call fieldio( '>', 'psv',  psv         )
    call fieldio( '>', 'su1',  t2(:,:,:,1) )
    call fieldio( '>', 'su2',  t2(:,:,:,2) )
    call fieldio( '>', 'su3',  t2(:,:,:,3) )
    call fieldio( '>', 'sum',  f2          )
    call fieldio( '>', 'sl',   sl          )
    call fieldio( '>', 'trup', trup        )
    call fieldio( '>', 'tarr', tarr        )
    call set_halo( f1,   -1.0, i1core, i2core )
    call set_halo( f2,   -1.0, i1core, i2core )
    call set_halo( tarr, -1.0, i1core, i2core )
    svmax = maxval( f1 )
    sumax = maxval( f2 )
    slmax = maxval( sl )
    tarrmax = maxval( tarr )
end if

end subroutine
end module

