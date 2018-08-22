! Attenuation
module m_attenuation
implicit none
integer, parameter :: ntau=2
real, dimension(ntau,ntau,ntau) :: tautem, tau,tau1,tau2,weights
real,dimension(ntau*ntau*ntau*2) :: coeff
real :: taumin, taumax, w0
complex, parameter:: sqrtm1=(0.0,1.0)
type :: volumeindex
	integer, allocatable, dimension(:,:,:) :: itx
	integer, allocatable, dimension(:,:,:) :: ity
	integer, allocatable, dimension(:,:,:) :: itz
end type
type( volumeindex ), pointer :: vindex
type :: strainstress
	integer :: nm(3)
	integer :: it
	real,allocatable,dimension(:,:,:,:) :: strain1,strain2
	real,allocatable,dimension(:,:,:,:) :: stress1,stress2
end type
type( strainstress ), pointer :: ss_tensor

contains

subroutine attenuaiton_init
use m_globals
real :: tmp, tauu
integer :: j,k,l,i1(3),i2(3)

w0 = 2*pi*fp

tautem = reshape((/1.0,6.0,7.0,4.0,8.0,3.0,2.0,5.0/),shape(tautem))

! assign taumin and taumax
if ((ex .LT. 0.65) .AND. (ex .GT. 0.01)) then   
      taumin=1./(10.*2.*3.1415927)*.2*fac
      taumax=1./(0.01*2.*3.1415927)*1.0*fac
end if

if ((ex .LT. 0.85) .AND. (ex .GT. 0.65)) then
      taumin=1./(12.*2.*3.1415927)*.5*fac
      taumax=1./(0.08*2.*3.1415927)*2.0*fac
end if
      
if ((ex .LT. 0.95) .AND. (ex .GT. 0.85)) then
      taumin=1./(280.*2.*3.1415927)*.1*fac
end if

! constant Q
if ((ex .LT. 0.01)) then
      taumin=1./(10.*2.*3.1415927)*.2*fac
      taumax=1./(0.01*2.*3.1415927)*1.0*fac
end if

! assign tau array (tau_k in paper)
do j = 1, ntau
  do k = 1, ntau
    do l = 1, ntau
    	tmp = tautem(j,k,l)
    	
        tmp = (tmp-0.5)/8
        tmp = 2.0*tmp - 1.0    
        
        tau(j,k,l) = exp(0.5*(log(taumax*taumin) + &
							log(taumax/taumin)*tmp))  	
    end do 
  end do
end do

!assign tau1 and tau2
!c     old way: no exponnetials
!c              tau1(i,j,k)=1/((tauu*dt1)+(1./2.))
!c              tau2(i,j,k)=(tauu*dt1)-(1./2.)
do j = 1, ntau
  do k = 1, ntau
    do l = 1, ntau
    	tauu=tau(j,k,l)
    	tau2(j,k,l)=exp(-dt/tauu)
        tau1(j,k,l)=0.5*(1.-tau2(j,k,l))
    end do
  end do
end do

! assign initial weights and coefficients
! for ex=0.1
if ((ex .LT. 0.15) .AND. (ex .GT. 0.01)) then
! 1,6,7,4,8,3,2,5
         weights(1,1,1)=0.3273
         weights(2,1,1)=1.0434
         weights(1,2,1)=0.044
         weights(2,2,1)=0.9393
         weights(1,1,2)=1.7268
         weights(2,1,2)=0.369
         weights(1,2,2)=0.8478
         weights(2,2,2)=0.4474
         
    coeff = (/7.3781, 4.1655, -83.1627, 13.1326, 69.0839, 0.4981, -37.6966, 5.5263, &
              -51.4056, 8.1934, 13.1865, 3.4775, -36.1049, 7.2107, 12.3809, 3.6117/)
              
!c for ex=0.2
elseif ((ex .LT. 0.25) .AND. (ex .GT. 0.15)) then
	     weights(1,1,1)=0.001
         weights(2,1,1)=1.0349
         weights(1,2,1)=0.0497
         weights(2,2,1)=1.0407
         weights(1,1,2)=1.7245
         weights(2,1,2)=0.2005
         weights(1,2,2)=0.804
         weights(2,2,2)=0.4452
         
    coeff = (/31.8902, 1.6126, -83.2611, 13.0749, 65.485, 0.5118, -42.02, 5.0875, &
              -49.2656, 8.1552, 25.7345, 2.2801, -40.8942, 7.9311, 7.0206, 3.4692/) 
              
!c for ex=0.3
elseif ((ex .LT. 0.35) .AND. (ex .GT. 0.25)) then  
	     weights(1,1,1)=0.001
         weights(2,1,1)=1.0135
         weights(1,2,1)=0.0621
         weights(2,2,1)=1.1003
         weights(1,1,2)=1.7198
         weights(2,1,2)=0.0918
         weights(1,2,2)=0.6143
         weights(2,2,2)=0.4659
         
    coeff = (/43.775, -0.1091, -83.1088, 13.0161, 60.9008, 0.592, -43.4857, 4.5869, &
              -45.3315, 8.0252, 34.3571, 1.199, -41.4422, 8.399, -2.8772, 3.5323/) 

!c for ex=0.4
elseif ((ex .LT. 0.45) .AND. (ex .GT. 0.35)) then
	     weights(1,1,1)=0.001
         weights(2,1,1)=0.9782
         weights(1,2,1)=0.082
         weights(2,2,1)=1.1275
         weights(1,1,2)=1.7122
         weights(2,1,2)=0.001
         weights(1,2,2)=0.4639
         weights(2,2,2)=0.509
         
    coeff = (/41.6858, -0.7344, -164.2252, 14.9961, 103.0301, -0.4199, -41.1157, 3.8266, &
              -73.0432, 8.5857, 38.0868, 0.3937, -43.2133, 8.6747, 5.6362, 3.3287/) 

!c for ex=0.5
elseif ((ex .LT. 0.55) .AND. (ex .GT. 0.45)) then
	     weights(1,1,1)=0.2073
         weights(2,1,1)=0.912
         weights(1,2,1)=0.1186
         weights(2,2,1)=1.081
         weights(1,1,2)=1.6984
         weights(2,1,2)=0.001
         weights(1,2,2)=0.1872
         weights(2,2,2)=0.6016
         
    coeff = (/20.0539, -0.4354, -81.6068, 12.8573, 45.9948, 1.1528, -23.07, 2.6719, &
              -27.8961, 7.1927, 31.4788, -0.0434, -25.1661, 8.245, -45.2178, 4.8476/) 

!c for ex=0.6
elseif ((ex .LT. 0.65) .AND. (ex .GT. 0.55)) then
!c 1,6,7,4,8,3,2,5
         weights(1,1,1)=0.3112
         weights(2,1,1)=0.8339
         weights(1,2,1)=0.1616
         weights(2,2,1)=1.0117
         weights(1,1,2)=1.6821
         weights(2,1,2)=0.001
         weights(1,2,2)=0.001
         weights(2,2,2)=0.7123

 	coeff=(/ 8.0848, -0.1968,-79.9715, 12.7318, 35.7155, 1.68, -13.0365 , 1.8101, &
         -13.2235, 6.3697,25.4548, -0.3947, -10.4478, 7.657, -75.9179 , 6.1791/)

!c for ex=0.7
elseif ((ex .LT. 0.75) .AND. (ex .GT. 0.65)) then
		 weights(1,1,1)=0.1219
         weights(2,1,1)=0.001
         weights(1,2,1)=0.5084
         weights(2,2,1)=0.2999
         weights(1,1,2)=1.2197
         weights(2,1,2)=0.001
         weights(1,2,2)=0.001
         weights(2,2,2)=1.3635

 	coeff=(/ 1.9975, 0.418,-76.6932, 11.3479, 40.7406, 1.9511, -2.7761, 0.5987, &
         0.0, 0.0, 0.0, 0.0, 41.317, 2.1874, -88.8095 , 11.0003/)

!c for ex=0.8
elseif ((ex .LT. 0.85) .AND. (ex .GT. 0.75)) then
!c 1,6,7,4,8,3,2,5
         weights(1,1,1)=0.0462
         weights(2,1,1)=0.001
         weights(1,2,1)=0.4157
         weights(2,2,1)=0.1585
         weights(1,1,2)=1.3005
         weights(2,1,2)=0.001
         weights(1,2,2)=0.001
         weights(2,2,2)=1.4986

	coeff= (/5.1672,0.2129, -46.506,11.7213,-5.8873,1.4279, -8.2448,0.3455, &
         	15.0254,-0.283,0.0,0.0,58.975,0.8131, -108.6828,12.4362/)

!c for ex=0.9
elseif ((ex .LT. 0.95) .AND. (ex .GT. 0.85)) then
!c 1,6,7,4,8,3,2,5
         weights(1,1,1)=0.001
         weights(2,1,1)=0.001
         weights(1,2,1)=0.1342
         weights(2,2,1)=0.1935
         weights(1,1,2)=1.5755
         weights(1,2,2)=0.001
         weights(2,2,2)=1.5297

 	coeff=(/ -0.8151, 0.1621,-61.9333, 12.5014, 0.0358, -0.0006, 0.0 , 0.0, &
         22.0291, -0.4022,0.0, 0.0, 56.0043, 0.7978, -116.9175 , 13.0244/)

! constant Q
elseif ((ex .LT. 0.01)) then
!c 1,6,7,4,8,3,2,5
         weights(1,1,1)=0.8867
         weights(2,1,1)=1.0440
         weights(1,2,1)=0.0423
         weights(2,2,1)=0.8110
         weights(1,1,2)=1.7275
         weights(2,1,2)=0.5615
         weights(1,2,2)=0.8323
         weights(2,2,2)=0.4641

	coeff=(/-27.5089, 7.4177, -82.8803, 13.1952, 72.0312, 0.5298, -34.1779, 6.0293, &
		-52.2607, 8.1754,-1.6270,    4.6858,-27.7770,    6.2852, 14.6295,    3.8839/)
end if ! end of assign weights and coeff array

!global coarse-grain index
allocate(vindex)
allocate(vindex%itx(nm(1),nm(2),nm(3)))
allocate(vindex%ity(nm(1),nm(2),nm(3)))
allocate(vindex%itz(nm(1),nm(2),nm(3)))

vindex%itx = 0
vindex%ity = 0
vindex%itz = 0

i1 = max( i1pml + 1, i1cell )
i2 = min( i2pml - 2, i2cell )

do l = i1(3),i2(3) !here exclude PML
  do k = i1(2),i2(2)
    do j = i1(1),i2(1)
    	vindex%itx(j,k,l) = j + nnoff(1) !global index (from 1 to nn) in material volume
    	vindex%ity(j,k,l) = k + nnoff(2)
    	vindex%itz(j,k,l) = l + nnoff(3)
    	!if (k + nnoff(2) == 1 .and. l + nnoff(3)== 1) write(0,'(5I3)') ipid,j,vindex%itx(j,k,l)
    	!itx ity itz in parallelization are verified
    	vindex%itx(j,k,l) = 2-mod(vindex%itx(j,k,l),2)
    	vindex%ity(j,k,l) = 2-mod(vindex%ity(j,k,l),2)
    	vindex%itz(j,k,l) = 2-mod(vindex%itz(j,k,l),2)
    	!if (k + nnoff(2) == 1 .and. l + nnoff(3)== 1) write(0,'(5I3)') ipid,j,vindex%itx(j,k,l)
    	!itx ity itz in parallelization are verified
    end do
  end do
end do
! global coarse-grain index (1 or 2) is finished (in halo and pml, index = 0)

!Initiate global stress and strain tensor
allocate(ss_tensor)
ss_tensor%nm = nm
!if (master) write(0,*) ss_tensor%it,ss_tensor%nm
allocate(ss_tensor%strain1(nm(1),nm(2),nm(3),3))
allocate(ss_tensor%strain2(nm(1),nm(2),nm(3),3))
allocate(ss_tensor%stress1(nm(1),nm(2),nm(3),3))
allocate(ss_tensor%stress2(nm(1),nm(2),nm(3),3))

end subroutine

!----------------------------------------------------------------------------------------
subroutine dispersion
use m_globals
complex :: value
real :: Mu1, qpinv,qsinv
real,dimension(ntau,ntau,ntau) :: weights_lop,weights_los
integer :: j,k,l,i1(3),i2(3), taul,tauk,tauj,num

i1 = max( i1pml + 1, i1cell )
i2 = min( i2pml - 2, i2cell )

do l = i1(3),i2(3) !here exclude PML
  do k = i1(2),i2(2)
    do j = i1(1),i2(1)

!!!!!!!!!!!!!!!!!!!!!
! P wave attenuation    	
    	if (qp(j,k,l) < 200.0) then ! use Graves and Day 2003
    		num = 0
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	weights_lop(tauj,tauk,taul) = coeff(1+num)/(qp(j,k,l)*qp(j,k,l)) + &
    		    	                              coeff(2+num)/qp(j,k,l)
    		    	num = num + 2
    		    end do
    		  end do
    		end do
    		
    		value = cmplx(0,0)
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	value=value+1./ &
    		    	    (1.-weights_lop(tauj,tauk,taul)/(1+sqrtm1*w0*tau(tauj,tauk,taul))) 
    		    end do
    		  end do
    		end do
			value=1./value
!c calulate the unrelax modulus corresonding to vp and then the velocity at infinite frequncy                                                                                                                                                       
            Mu1=mr(j,k,l)*s1(j,k,l)*s1(j,k,l)/(8.*real(value))  
            s1(j,k,l) = sqrt(Mu1/mr(j,k,l))  	! P wave modified with dispersion		
    	else
    		value = cmplx(0,0)
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	value=value+weights(tauj,tauk,taul) / &
    		    	          ((sqrtm1*w0*tau(tauj,tauk,taul)+1)*qp(j,k,l)*fac**ex)     
    		    end do
    		  end do
    		end do
    		Mu1=mr(j,k,l)*s1(j,k,l)*s1(j,k,l)/(abs(1.-value))
    		s1(j,k,l) = sqrt(Mu1/mr(j,k,l))  	! P wave modified with dispersion
    	end if


!!!!!!!!!!!!!!!!!!!!!
! S wave attenuation
    	if (qs(j,k,l) < 200.0) then ! use Graves and Day 2003
    		num = 0
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	weights_los(tauj,tauk,taul) = coeff(1+num)/(qs(j,k,l)*qs(j,k,l)) + &
    		    	                              coeff(2+num)/qs(j,k,l)
    		    	num = num + 2
    		    end do
    		  end do
    		end do
    		
    		value = cmplx(0,0)
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	value=value+1./ &
    		    	    (1.-weights_los(tauj,tauk,taul)/(1+sqrtm1*w0*tau(tauj,tauk,taul))) 
    		    end do
    		  end do
    		end do
			value=1./value
!c calulate the unrelax modulus corresonding to vp and then the velocity at infinite frequncy                                                                                                                                                       
            Mu1=mr(j,k,l)*s2(j,k,l)*s2(j,k,l)/(8.*real(value))  
            s2(j,k,l) = sqrt(Mu1/mr(j,k,l))  	! S wave modified with dispersion		
    	else
    		value = cmplx(0,0)
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		    	value=value+weights(tauj,tauk,taul) / &
    		    	          ((sqrtm1*w0*tau(tauj,tauk,taul)+1)*qs(j,k,l)*fac**ex)     
    		    end do
    		  end do
    		end do
    		Mu1=mr(j,k,l)*s2(j,k,l)*s2(j,k,l)/(abs(1.-value))
    		s2(j,k,l) = sqrt(Mu1/mr(j,k,l))  	! S wave modified with dispersion
    	end if

!!!!!!!!!!!!!!!!!!
! regenerate a matrix of qp and qs
! KW here I change As and Ap for frequncy dpendent Q
!                tmppq(i,j,k)=tmp1*qpinv/(1.0-tmp2*qpinv)
!                tmpsq(i,j,k)=tmp1*qsinv/(1.0-tmp2*qsinv)

!                tmppq(i,j,k)=qpinv*   (1/Q0)*fac**(-ex)
!                 tmpsq(i,j,k)=qsinv*   (1/Q0)*fac**(-ex)
! Ap and As
	if (qp(j,k,l) <= 0.0) then
    	qpinv=0.0
    	qsinv=0.0
    else
        qpinv=1./qp(j,k,l)
        qsinv=1./qs(j,k,l)
    end if
    qp(j,k,l) = qpinv*fac**(-ex) 
    qs(j,k,l) = qsinv*fac**(-ex) 
    
    end do
  end do
end do

end subroutine

!----------------------------------------------------------------------------------------
subroutine savestrain
use m_globals

!save strain into strainstress_tensor
ss_tensor%it = it
ss_tensor%strain1 = 0.0 !optional
ss_tensor%strain2 = 0.0

ss_tensor%strain1 = w1
ss_tensor%strain2 = w2
if (verb .and. master) then
	write(0,'(A,I10)') 'Save Strain at ', it
	!write(0,*) maxval(ss_tensor%strain1),maxval(ss_tensor%strain2)
end if
end subroutine

!----------------------------------------------------------------------------------------
subroutine savestress
use m_globals

!save strain into strainstress_tensor
!ss_tensor%it = it
ss_tensor%stress1 = 0.0 !optional
ss_tensor%stress2 = 0.0
ss_tensor%stress1 = w1
ss_tensor%stress2 = w2
if (verb .and. master) then
	write(0,'(A,I10)') 'Save Stress at ', it
	!write(0,*) maxval(ss_tensor%stress1),maxval(ss_tensor%stress2)
end if
end subroutine

!----------------------------------------------------------------------------------------
subroutine updatememoryvariable !main part to update stress and memory-variables
use m_globals
use m_fieldio

integer :: i1(3),i2(3),j,k,l,itx,ity,itz,tauj,tauk,taul,num
real,dimension(ntau,ntau,ntau) :: weights_lop,weights_los
real :: weights_p,weights_s,vx1,vx2
real :: vt(6),vh1,vh2,va1

i1 = max( i1pml + 1, i1cell )
i2 = min( i2pml - 2, i2cell )

do l = i1(3),i2(3) !here exclude PML
  do k = i1(2),i2(2)
    do j = i1(1),i2(1)
    	itx = vindex%itx(j,k,l)
    	ity = vindex%ity(j,k,l)
    	itz = vindex%itz(j,k,l)
    	
    	vx1 = tau1(itx,ity,itz)
    	vx2 = tau2(itx,ity,itz)
    	
    	!no harmonic averaging for lam and mu here (not similar in AWP)
        ! P wave
    	if (1./qp(j,k,l) < 200) then !harmonic averaging moduli in Graves and Day, 2003
            num = 0
    		do tauj = 1, ntau
    		  do tauk = 1, ntau
    		    do taul = 1, ntau
    		      weights_lop(tauj,tauk,taul) = coeff(1+num)*(qp(j,k,l)*qp(j,k,l)) + &
                                                  coeff(2+num)*qp(j,k,l)
                  num = num + 2
    		    end do
    		  end do
    		end do
            weights_p=weights_lop(itx,ity,itz)
    	else
    	    weights_p=8*weights(itx,ity,itz)*qp(j,k,l)
    	end if
    	
        ! S wave
        if (1./qs(j,k,l) < 200) then !harmonic averaging moduli in Graves and Day, 2003
            num = 0
            do tauj = 1, ntau
              do tauk = 1, ntau
                do taul = 1, ntau
                  weights_los(tauj,tauk,taul) = coeff(1+num)*(qs(j,k,l)*qs(j,k,l)) + &
                                                  coeff(2+num)*qs(j,k,l)
                  num = num + 2
                end do
              end do
            end do
            weights_s=weights_los(itx,ity,itz)
        else
            weights_s=8*weights(itx,ity,itz)*qs(j,k,l)
        end if       
    	
        !start to modify stress 

        vt(1) = ss_tensor%strain1(j,k,l,1)
        vt(2) = ss_tensor%strain1(j,k,l,2)
        vt(3) = ss_tensor%strain1(j,k,l,3)
        vt(4) = ss_tensor%strain2(j,k,l,1)
        vt(5) = ss_tensor%strain2(j,k,l,2)
        vt(6) = ss_tensor%strain2(j,k,l,3)

        ! first normal stress
        ss_tensor%stress1(j,k,l,1) = ss_tensor%stress1(j,k,l,1) + dt * attr1(j,k,l,1)
        ss_tensor%stress1(j,k,l,2) = ss_tensor%stress1(j,k,l,2) + dt * attr1(j,k,l,2)
        ss_tensor%stress1(j,k,l,3) = ss_tensor%stress1(j,k,l,3) + dt * attr1(j,k,l,3)

        vh1 = -mu(j,k,l)
        va1 = -(lam(j,k,l)+2*mu(j,k,l))*(vt(1)+vt(2)+vt(3))/2

        attr1(j,k,l,1) = vx2 * attr1(j,k,l,1) - (weights_s*vh1*(vt(2)+vt(3))- &
                                                    weights_p*va1) * vx1 *2
        attr1(j,k,l,2) = vx2 * attr1(j,k,l,2) - (weights_s*vh1*(vt(1)+vt(3))- &
                                                    weights_p*va1) * vx1 *2
        attr1(j,k,l,3) = vx2 * attr1(j,k,l,3) - (weights_s*vh1*(vt(1)+vt(1))- &
                                                    weights_p*va1) * vx1 *2

        w1(j,k,l,1) = ss_tensor%stress1(j,k,l,1) + dt * attr1(j,k,l,1)
        w1(j,k,l,2) = ss_tensor%stress1(j,k,l,2) + dt * attr1(j,k,l,2)
        w1(j,k,l,3) = ss_tensor%stress1(j,k,l,3) + dt * attr1(j,k,l,3)

        ! then shear stress
        ss_tensor%stress2(j,k,l,1) = ss_tensor%stress2(j,k,l,1) + dt * attr2(j,k,l,1)
        ss_tensor%stress2(j,k,l,2) = ss_tensor%stress2(j,k,l,2) + dt * attr2(j,k,l,2)
        ss_tensor%stress2(j,k,l,3) = ss_tensor%stress2(j,k,l,3) + dt * attr2(j,k,l,3)

        vh2 = -mu(j,k,l)/2.
        attr2(j,k,l,1) = vx2 * attr2(j,k,l,1) + weights_s*vh2*vt(4)*vx1*2.
        attr2(j,k,l,2) = vx2 * attr2(j,k,l,2) + weights_s*vh2*vt(5)*vx1*2.
        attr2(j,k,l,3) = vx2 * attr2(j,k,l,3) + weights_s*vh2*vt(6)*vx1*2.

        w2(j,k,l,1) = ss_tensor%stress2(j,k,l,1) + dt * attr2(j,k,l,1)
        w2(j,k,l,2) = ss_tensor%stress2(j,k,l,2) + dt * attr2(j,k,l,2)
        w2(j,k,l,3) = ss_tensor%stress2(j,k,l,3) + dt * attr2(j,k,l,3)

        !if(debug > 0 .and. master ) &
        !           write(0,*) j,k,l,vx1,vx2,weights_p,weights_s
    end do
  end do
end do

!if(debug > 0 .and. master) write(0,*) attr1(28,30,30,1)

call fieldio( '>', 'attr11', attr1(:,:,:,1) )
call fieldio( '>', 'attr22', attr1(:,:,:,2) )
call fieldio( '>', 'attr33', attr1(:,:,:,3) )
call fieldio( '>', 'attr23', attr2(:,:,:,1) )
call fieldio( '>', 'attr31', attr2(:,:,:,2) )
call fieldio( '>', 'attr12', attr2(:,:,:,3) )

end subroutine

end module

