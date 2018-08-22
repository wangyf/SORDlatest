! Read model parameters
module m_parameters
implicit none
contains

subroutine read_parameters
use m_globals
use m_fieldio
integer :: ios, i
character(12) :: key
character(1) :: op
character(256) :: line

! I/O pointers
allocate( io0 )
io => io0
io%next => io0
io%field = 'head'

open( 1, file='parameters.py', status='old' )

doline: do

! Read line
read( 1, '(a)', iostat=ios ) line
if ( ios /= 0 ) exit doline

! Strip comments and punctuation
str = line
i = scan( str, '#' )
if ( i > 0 ) str(i:) = ' '
do
    i = scan( str, "()[]{}'," )
    if ( i == 0 ) exit
    str(i:i) = ' '
end do

! Read key val pair
if ( str == '' ) cycle doline
read( str, *, iostat=ios ) key

! Select input key
select case( key )
case( 'fieldio', '' )
case( 'nn' );           read( str, *, iostat=ios ) key, op, nn
case( 'nt' );           read( str, *, iostat=ios ) key, op, nt
case( 'dx' );           read( str, *, iostat=ios ) key, op, dx
case( 'dt' );           read( str, *, iostat=ios ) key, op, dt
case( 'tm0' );          read( str, *, iostat=ios ) key, op, tm0
case( 'affine' );       read( str, *, iostat=ios ) key, op, affine
case( 'n1expand' );     read( str, *, iostat=ios ) key, op, n1expand
case( 'n2expand' );     read( str, *, iostat=ios ) key, op, n2expand
case( 'rexpand' );      read( str, *, iostat=ios ) key, op, rexpand
case( 'gridnoise' );    read( str, *, iostat=ios ) key, op, gridnoise
case( 'oplevel' );      read( str, *, iostat=ios ) key, op, oplevel
case( 'rho1' );         read( str, *, iostat=ios ) key, op, rho1
case( 'rho2' );         read( str, *, iostat=ios ) key, op, rho2
case( 'vp1' );          read( str, *, iostat=ios ) key, op, vp1
case( 'vp2' );          read( str, *, iostat=ios ) key, op, vp2
case( 'vs1' );          read( str, *, iostat=ios ) key, op, vs1
case( 'vs2' );          read( str, *, iostat=ios ) key, op, vs2
case( 'gam1' );         read( str, *, iostat=ios ) key, op, gam1
case( 'gam2' );         read( str, *, iostat=ios ) key, op, gam2
case( 'vdamp' );        read( str, *, iostat=ios ) key, op, vdamp
case( 'hourglass' );    read( str, *, iostat=ios ) key, op, hourglass
case( 'bc1' );          read( str, *, iostat=ios ) key, op, bc1
case( 'bc2' );          read( str, *, iostat=ios ) key, op, bc2
case( 'npml' );         read( str, *, iostat=ios ) key, op, npml
case( 'i1pml' );        read( str, *, iostat=ios ) key, op, i1pml
case( 'i2pml' );        read( str, *, iostat=ios ) key, op, i2pml
case( 'ppml' );         read( str, *, iostat=ios ) key, op, ppml
case( 'vpml' );         read( str, *, iostat=ios ) key, op, vpml
case( 'ihypo' );        read( str, *, iostat=ios ) key, op, ihypo
case( 'source' );       read( str, *, iostat=ios ) key, op, source
case( 'timefunction' ); read( str, *, iostat=ios ) key, op, timefunction
case( 'period' );       read( str, *, iostat=ios ) key, op, period
case( 'source1' );      read( str, *, iostat=ios ) key, op, source1
case( 'source2' );      read( str, *, iostat=ios ) key, op, source2
case( 'nsource' );      read( str, *, iostat=ios ) key, op, nsource
case( 'ivols' );        read( str, *, iostat=ios ) key, op, ivols
case( 'eplasticity' );  read( str, *, iostat=ios ) key, op, eplasticity
case( 'tv' );           read( str, *, iostat=ios ) key, op, tv
case( 'friction' );     read( str, *, iostat=ios ) key, op, friction
case( 'pcdep' );        read( str, *, iostat=ios ) key, op, pcdep
case( 'faultnormal' );  read( str, *, iostat=ios ) key, op, faultnormal
case( 'slipvector' );   read( str, *, iostat=ios ) key, op, slipvector
case( 'faultopening' ); read( str, *, iostat=ios ) key, op, faultopening
case( 'vrup' );         read( str, *, iostat=ios ) key, op, vrup
case( 'rcrit' );        read( str, *, iostat=ios ) key, op, rcrit
case( 'trelax' );       read( str, *, iostat=ios ) key, op, trelax
case( 'rnucl' );        read( str, *, iostat=ios ) key, op, rnucl
case( 'tmnucl' );       read( str, *, iostat=ios ) key, op, tmnucl
case( 'delts' );        read( str, *, iostat=ios ) key, op, delts
case( 'svtol' );        read( str, *, iostat=ios ) key, op, svtol
case( 'np3' );          read( str, *, iostat=ios ) key, op, np3
case( 'itstats' );      read( str, *, iostat=ios ) key, op, itstats
case( 'itio' );         read( str, *, iostat=ios ) key, op, itio
case( 'itcheck' );      read( str, *, iostat=ios ) key, op, itcheck
case( 'itstop' );       read( str, *, iostat=ios ) key, op, itstop
case( 'debug' );        read( str, *, iostat=ios ) key, op, debug
case( 'mpin' );         read( str, *, iostat=ios ) key, op, mpin
case( 'mpout' );        read( str, *, iostat=ios ) key, op, mpout
case default
    select case( key(1:1) )
    case( '=', '+' )
        call pappend
        io%ib = -1
        !XXXread( str, *, iostat=ios ) io%mode, io%nc
        read( str, *, iostat=ios ) io%mode, io%nc, io%tfunc, &
            io%period, io%x1, io%x2, io%nb, io%ii, io%filename, &
            io%val, io%field
    case default; ios = 1
    end select
end select

! Error check
if ( ios /= 0 ) then
    if ( master ) write( 0, * ) 'bad input: ', trim( line )
    stop
end if

end do doline

close( 1 )

end subroutine

end module

