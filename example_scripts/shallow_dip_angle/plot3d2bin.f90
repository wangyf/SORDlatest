!	Program plot3d2bin
! 	This program will convert the ascii file of grid file in the format of plot3d
!	into binary file readable to SORD. The grid file is an outcome of Pointwise,inc
!	software, which will generate a single-block structured mesh. During this process,
!	this program can refine the mesh based on higher or similar order interpolation. Plus,
!	this program will have a byproduct of a vtk-format grid file. 
!
!
!	Author: Yongfei Wang
!	  
!   To do:
!		1. refined output using interpolation
program plot3d2bin
integer :: ni,nj,nk,n1f,n2f,n3f,nindex1,nindex2
real,allocatable,dimension(:,:,:) :: x,y,z,xo,yo,zo
integer ::m,nblocks,numfault, nargc
integer,allocatable,dimension(:,:) :: faults
character(len=20) :: filename,arg

! read in parameters
write(0,*) 'Plot3d2bin'
write(0,*) 
write(0,*) 'This program will convert the ascii file of grid file in the format of plot3d'
write(0,*) 'into binary file readable to SORD. The grid file is an outcome of Pointwise,inc'
write(0,*) 'software, which will generate a single-block structured mesh. During this process,'
write(0,*) 'this program can refine the mesh based on higher or similar order interpolation.'
write(0,*) 'The plot3d file should be output with the option of ''unformatted,single and 3d'''
write(0,*) 'in Pointwise.'
write(0,*)
write(0,*) 'Author: Yongfei Wang'
write(0,*)
write(0,*) 'Usage: plot3d2bin test.x faultnumber fault1_normal index fault2_normal index ...'
write(0,*) 'The same fault_normal,index numbers should be ascending e.g., 1 10 1 15 2 10 3 15'
write(0,*)

if (iargc() < 1) then
	stop
end if

nargc = iargc()
! read file name
CALL getarg(1, arg)
filename = trim(arg)
WRITE (0,*) 'Reading ... ',arg

! read fault number
CALL getarg(2, arg)
read(arg,*) numfault
WRITE (0,'(A,I2,A)') 'There are ',numfault ,' faults involved'
write(0,*)

! read fault information
n1f = 0
n2f = 0
n3f = 0
allocate(faults(numfault,2)) !(:,1:2) fault_normal,fault index
do m = 1,numfault
	CALL getarg((m-1)*2+3, arg)
	read(arg,*) faults(m,1)
	CALL getarg((m-1)*2+4, arg)
	read(arg,*) faults(m,2)
	WRITE (0,'(I2,A,I2,A,I4)') m ,' th fault is normal to ',faults(m,1),' direction at ',faults(m,2)
	select case (faults(m,1))
		case (1)
			n1f = n1f + 1
		case (2)
			n2f = n2f + 1
		case (3)
			n3f = n3f + 1
		case default
			write(0,*) 'Fault normal must be 1, 2 or 3'
			stop
	end select
end do

! read in plot3d data (binary unformatted file) 
! Ask for unformatted,single,3D in export option of Pointwise
open ( unit=6, form='unformatted', file=trim(filename), status='old')
!	read(6) nblocks
!	if(nblocks/=1) then
!		write(0,*) 'This file contains ',nblocks,' blocks more than 1'
!		stop
!	end if
    read(6) ni,nj,nk
    write(0,*) 'Mesh Dimension: ',ni,nj,nk
    allocate(x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk))
    
   read(6) &
      ((( x(i,j,k), i=1,ni), j=1,nj), k=1,nk), &
	  ((( y(i,j,k), i=1,ni), j=1,nj), k=1,nk), &
	  ((( z(i,j,k), i=1,ni), j=1,nj), k=1,nk)

close( unit = 6)


! consider split nodes of faults
ni = ni + n1f
nj = nj + n2f
nk = nk + n3f
allocate(xo(ni,nj,nk),yo(ni,nj,nk),zo(ni,nj,nk))

! modify read-in grid coordinate by adding split nodes
! x direction
do m = 1,n1f
	if (m==1) then
		nindex1 = 1
	else
		nindex1 = faults(m-1,2) + m - 1 
	end if
	nindex2 = faults(m,2) + m - 1
	
	xo(nindex1:nindex2,:,:) = x(nindex1 + 1 - m:nindex2 + 1 - m,:,:)
	yo(nindex1:nindex2,:,:) = y(nindex1 + 1 - m:nindex2 + 1 - m,:,:)
	zo(nindex1:nindex2,:,:) = z(nindex1 + 1 - m:nindex2 + 1 - m,:,:)
	
	if (m==n1f) then
		nindex1 = faults(m,2) + m
		nindex2 = ni
		xo(nindex1:nindex2,:,:) = x(faults(m,2):ni-n1f,:,:)
		yo(nindex1:nindex2,:,:) = y(faults(m,2):ni-n1f,:,:)
		zo(nindex1:nindex2,:,:) = z(faults(m,2):ni-n1f,:,:)
	end if
end do

! y direction
do m = 1,n2f
	if (m==1) then
		nindex1 = 1
	else
		nindex1 = faults(m-1 + n1f,2) + m - 1 
	end if
	nindex2 = faults(m+ n1f,2) + m - 1
	
	xo(:,nindex1:nindex2,:) = x(:,nindex1 + 1 - m:nindex2 + 1 - m,:)
	yo(:,nindex1:nindex2,:) = y(:,nindex1 + 1 - m:nindex2 + 1 - m,:)
	zo(:,nindex1:nindex2,:) = z(:,nindex1 + 1 - m:nindex2 + 1 - m,:)
	
	if (m==n2f) then
		nindex1 = faults(m+ n1f,2) + m
		nindex2 = ni
		xo(:,nindex1:nindex2,:) = x(:,faults(m+ n1f,2):ni-n1f,:)
		yo(:,nindex1:nindex2,:) = y(:,faults(m+ n1f,2):ni-n1f,:)
		zo(:,nindex1:nindex2,:) = z(:,faults(m+ n1f,2):ni-n1f,:)
	end if
end do

! z direction
do m = 1,n3f
	if (m==1) then
		nindex1 = 1
	else
		nindex1 = faults(m-1 + n1f +n2f,2) + m - 1 
	end if
	nindex2 = faults(m+ n1f +n2f,2) + m - 1
	
	xo(:,:,nindex1:nindex2) = x(:,:,nindex1 + 1 - m:nindex2 + 1 - m)
	yo(:,:,nindex1:nindex2) = y(:,:,nindex1 + 1 - m:nindex2 + 1 - m)
	zo(:,:,nindex1:nindex2) = z(:,:,nindex1 + 1 - m:nindex2 + 1 - m)
	
	if (m==n2f) then
		nindex1 = faults(m+ n1f +n2f,2) + m
		nindex2 = ni
		xo(:,:,nindex1:nindex2) = x(:,:,faults(m+ n1f +n2f,2):ni-n1f)
		yo(:,:,nindex1:nindex2) = y(:,:,faults(m+ n1f +n2f,2):ni-n1f)
		zo(:,:,nindex1:nindex2) = z(:,:,faults(m+ n1f +n2f,2):ni-n1f)
	end if
end do


! write out as binary file
! By default, the grid file will be x.bin, y.bin and z.bin
open ( unit=1, form='unformatted', file='x.bin', access='direct',recl=4,status='replace' )
open ( unit=2, form='unformatted', file='y.bin', access='direct',recl=4,status='replace' )
open ( unit=3, form='unformatted', file='z.bin', access='direct',recl=4,status='replace' )
m = 0
write(0,*) 'Writing ... x.bin y.bin z.bin'
do k = 1,nk
	do j = 1,nj
		do i = 1,ni
			m = m + 1
			write(1,rec=m) xo(i,j,k)
			write(2,rec=m) yo(i,j,k)
			write(3,rec=m) zo(i,j,k)
		end do
	end do
end do
close (unit = 1)
close (unit = 2)
close (unit = 3)
write(0,*) 'Finish'
end program
