!----------------------------------------------------------
!
!----------------------------------------------------------

	subroutine SingleStep

	implicit none
	include 'dpdflow.h'

!   Parameters
!   Locals

    integer k, isample, stepUnEquil
    real    Sampletime

    stepCount = stepCount + 1

    if(stepCount .gt. stepEquil) timeNow = timeNow + deltaT

    if(mod(stepCount,200) .eq. 0) write(*,*) 'stepCount = ', stepCount
      


!    if(stepCount .ge. 800) print*, 'v0'     !test
    if(runId .eq. 3) then
       call Velocity_Verlet0
    else
       call Velocity_Verlet
    endif
    !if(stepCount .ge. 600) print*, 'v0'     !test

    if(stepCount .lt. stepEquil) then
        call EvalProps
        call AccumProps (1)

!       sInitKinEnergy = sInitKinEnergy + kinEnergy
        if(mod(stepCount, stepAvg) .eq. 0) then
            call AccumProps (2)
            call PrintSummary
            call AccumProps (0)
        endif
    endif
     !if(stepCount .ge. 600) print*, 'v1'     !test
   
    if(stepCount .gt. 1000 .and. mod((stepCount - 1000),              &
       stepChainProps) .eq. 0 .and. stepCount .le. stepEquil          &
       .and. nChain .gt.0) call EvalChainProps

    if(nChainconf .gt. 0) then
        if(stepCount .eq. stepEquil .and. nChain .gt.0) call PrintChainConf
        if(stepCount .gt. stepEquil .and. mod((stepCount - stepEquil),    &
           stepChainProps*limitChainProps) .eq. 0 .and. nChain .gt.0)     &
           call PrintChainConf
    endif

    if(nDp .gt. 0) then
       if(stepCount .eq. stepEquil) call PrintDropConf
       if(stepCount .gt. stepEquil .and. mod((stepCount - stepEquil), &
	   (stepChainProps*limitChainProps)) .eq. 0 ) call PrintDropConf
    endif
!if(stepCount .ge. 600) print*, 'v2'     !test
    stepUnEquil = stepCount - stepEquil
    if(stepUnEquil .lt. startSample) then
        isample = 0
    else
        isample = stepUnEquil - stepSample*int(stepUnEquil/stepSample) + 1
    endif
    if(isample .eq. 1) Sampletime = timeNow

    if(isample .gt. 0 .and. isample .le. limitGrid*stepGrid           &
       .and. mod(isample,stepGrid) .eq. 0) then
        countGrid = countGrid + 1
        call GridAverage(1)

        if(mod(countGrid,limitGrid) .eq. 0) then
            call GridAverage(2)
            call PrintLocalProps(Sampletime)
            call GridAverage(0)
        endif
    endif
!if(stepCount .ge. 600) print*, 'v3'     !test
    if(timeNow .ge. timeSteady .and. mod(stepCount,3) .eq. 0          &
       .and. nChain .gt. 0) call GridAvChainProps(1)
! if(stepCount .ge. stepEquil) print*, stepCount           
!if(stepCount .ge. 600) print*, 'v4'     !test
    !return
    end

!-----------------------------------------------------------------

	subroutine EvalProps

	implicit none
	include 'dpdflow.h'

!	parameters
	
!	Locals

	integer	k, n
	real*8	v, vvSum

	vSum = 0.
	vvSum = 0.
	
	do k = 1, NDIM
	   do n = nWallAtom + 1, nAtom
	      v = rv(n,k)
	      vSum = vSum + v
	      vvSum = vvSum + v*v
	   enddo
	enddo

	kinEnergy = 0.5*vvSum / nFreeAtom
	potEnergy = uSum / nFreeAtom
	totEnergy = kinEnergy + potEnergy
	pressure = density*(vvSum + virSum) / (nFreeAtom * NDIM)
	vSum = vSum / nFreeAtom

	return

	end

!---------------------------------------------------------------

	subroutine AccumProps (icord)

	implicit none
	include 'dpdflow.h'

!	Parameters
	
	integer, intent(in):: icord
	
!	Locals
	 
	real*8	d1, d2

	if(icord .eq. 0) then

	   sTotEnergy = 0.
	   ssTotEnergy = 0.
	   sKinEnergy = 0.
	   ssKinEnergy = 0.
	   sPressure = 0.
	   ssPressure = 0.

	elseif(icord .eq. 1) then

	   sTotEnergy = sTotEnergy + totEnergy
	   ssTotEnergy = ssTotEnergy + totEnergy**2
	   sKinEnergy = sKinEnergy + kinEnergy
	   ssKinEnergy = ssKinEnergy + kinEnergy**2
	   sPressure = sPressure + pressure
	   ssPressure = ssPressure + pressure**2

	elseif(icord .eq. 2) then
	   sTotEnergy = sTotEnergy/stepAvg
	   d1 = ssTotEnergy/stepAvg
	   d2 = sTotEnergy**2
	   if(abs(d1 - d2) .lt. 1.e-10) then
	      ssTotEnergy = 0.
	   else
	      ssTotEnergy = sqrt(d1 - d2)
	   endif
	   sKinEnergy = sKinEnergy/stepAvg
	   d1 = ssKinEnergy/stepAvg
	   d2 = sKinEnergy**2
	   if(abs(d1 - d2) .lt. 1.e-10) then
	      ssKinEnergy = 0.
	   else
	      ssKinEnergy = sqrt(d1 - d2)
	   endif
	   sPressure = sPressure/stepAvg
	   d1 = ssPressure/stepAvg
	   d2 = sPressure**2
	   if(abs(d1 - d2) .lt. 1.e-10) then
	      ssPressure = 0.
	   else
	      ssPressure = sqrt(d1 - d2)
	   endif

	endif

	return

	end

!-------------------------------------------------------------------------
	
	subroutine PrintSummary

	implicit none
	include 'dpdflow.h'

	write(lst1,100) stepCount, timeNow, vSum, sTotEnergy, ssTotEnergy,&
                    sKinEnergy, ssKinEnergy, sPressure, ssPressure  

100 format(1x, i8, 2f9.4, f10.4, 3f8.4, f9.4, f8.4)

	return

	end

!------------------------------------------------------------------------

	subroutine GridAverage (opCode)

	implicit none
	include 'dpdflow.h'

!	Parameters

	integer, intent(in):: opCode

!	Locals

	integer	c, j, k, n
	real*8	cbvolm, pSum, ri(NDIM)

	if(opCode .eq. 0) then

	   do j = 1, NHIST
	      do n = 1, hSize
	         histGrid(n,j) = 0.
	      enddo
	   enddo
	   do j = 1, 7
	      do n = 1, hsize
	         strsGrid(n,j) = 0.
	      enddo
	   enddo

	elseif(opCode .eq. 1) then

	   do n = nWallAtom + 1, nAtom
!      do n = nChainEnd + 1, nAtom
!DP	   do n = 1, nAtom

!DP 	  c = int((r(n,3) + regionH(3))*sizeHistGrid(2)/ &
!     	      region(3))*sizeHistGrid(1) + int((r(n,1) + &
!DP     	      regionH(1))*sizeHistGrid(1)/region(1)) + 1
! 	      c = int((r(n,3) + initUcell(3)/2*gap(3))*sizeHistGrid(2)/ &
!     	         initUcell(3)/gap(3))*sizeHistGrid(1) + int((r(n,1) + &
!     	         regionH(1))*sizeHistGrid(1)/region(1)) + 1

          c = (int((r(n,3) + regionH(3))*sizeHistGrid(3)/region(3)) &
              *sizeHistGrid(2) + int((r(n,2) + regionH(2))          &
              *sizeHistGrid(2)/region(2)))*sizeHistGrid(1) +        &
              int((r(n,1) + regionH(1))*sizeHistGrid(1)/region(1)) + 1

!         c = (int((r(n,3) + initUcell(3)/2*gap(3))*sizeHistGrid(3)/  &
!              initUcell(3)/gap(3))*sizeHistGrid(2) +                 &
!              int((r(n,2) + regionH(2))*sizeHistGrid(2)/             &
!              region(2)))*sizeHistGrid(1) +                          &
!              int((r(n,1) + regionH(1))*sizeHistGrid(1)/region(1)) + 1

	      histGrid(c,1) = histGrid(c,1) + 1.0
 	      histGrid(c,2) = histGrid(c,2) + rv(n,1)**2 +                &
     	        		     rv(n,2)**2 + rv(n,3)**2
	      histGrid(c,3) = histGrid(c,3) + rv(n,1)              
	      histGrid(c,4) = histGrid(c,4) + rv(n,2)              
	      histGrid(c,5) = histGrid(c,5) + rv(n,3)              
	      if(n .le. nChainend) histGrid(c,6) = histGrid(c,6) + 1.0

 	      strsGrid(c,1) = strsGrid(c,1) + mass*rv(n,1)**2 + 0.5*rforce(n,1)
 	      strsGrid(c,2) = strsGrid(c,2) + mass*rv(n,2)**2 + 0.5*rforce(n,2)
 	      strsGrid(c,3) = strsGrid(c,3) + mass*rv(n,3)**2 + 0.5*rforce(n,3)
 	      strsGrid(c,4) = strsGrid(c,4) + mass*rv(n,1)*rv(n,2) + 0.5*rforce(n,4)
 	      strsGrid(c,5) = strsGrid(c,5) + mass*rv(n,1)*rv(n,3) + 0.5*rforce(n,5)
 	      strsGrid(c,6) = strsGrid(c,6) + mass*rv(n,2)*rv(n,3) + 0.5*rforce(n,6)

!	contribution from intermolecular forces

!	      strsGrid(c,1) = strsGrid(c,1) + 0.5*rforce(n,1)
!	      strsGrid(c,2) = strsGrid(c,2) + 0.5*rforce(n,2)
!	      strsGrid(c,3) = strsGrid(c,3) + 0.5*rforce(n,3)
!	      strsGrid(c,4) = strsGrid(c,4) + 0.5*rforce(n,4)
!	      strsGrid(c,5) = strsGrid(c,5) + 0.5*rforce(n,5)
!	      strsGrid(c,6) = strsGrid(c,6) + 0.5*rforce(n,6)

!	conritbution from momentium of peculiar velocity

!	      strsGrid(c,1) = strsGrid(c,1) + rv(n,1)**2
!	      strsGrid(c,2) = strsGrid(c,2) + rv(n,2)**2
!	      strsGrid(c,3) = strsGrid(c,3) + rv(n,3)**2
!	      strsGrid(c,4) = strsGrid(c,4) + rv(n,1)*rv(n,2)
!	      strsGrid(c,5) = strsGrid(c,5) + rv(n,1)*rv(n,3)
!	      strsGrid(c,6) = strsGrid(c,6) + rv(n,2)*rv(n,3)

	   enddo
	elseif(opCode .eq. 2) then

!	   pSum = 0.

       cbvolm = limitGrid*binvolm

	   do n = 1, hSize
	      if(histGrid(n,1) .gt. 0.) then
 	         do j = 2, NHIST
	            histGrid(n,j) = histGrid(n,j)/histGrid(n,1)
	         enddo
 	         histGrid(n,2) = (histGrid(n,2) - (histGrid(n,3)**2       &
      	                 + histGrid(n,4)**2 + histGrid(n,5)**2))/3.0
             histGrid(n,2) = histGrid(n,2) * mass

 	         do j = 1, 6
 		        strsGrid(n,j) = strsGrid(n,j)/histGrid(n,1)
 	         enddo

 	         strsGrid(n,1) = strsGrid(n,1) - mass * histGrid(n,3)**2
 	         strsGrid(n,2) = strsGrid(n,2) - mass * histGrid(n,4)**2
 	         strsGrid(n,3) = strsGrid(n,3) - mass * histGrid(n,5)**2
 	         strsGrid(n,4) = strsGrid(n,4) - mass * histGrid(n,3)*histGrid(n,4)
 	         strsGrid(n,5) = strsGrid(n,5) - mass * histGrid(n,3)*histGrid(n,5)
 	         strsGrid(n,6) = strsGrid(n,6) - mass * histGrid(n,4)*histGrid(n,5)

! 	         strsGrid(n,1) = strsGrid(n,1) - histGrid(n,3)**2
! 	         strsGrid(n,2) = strsGrid(n,2) - histGrid(n,4)**2
! 	         strsGrid(n,3) = strsGrid(n,3) - histGrid(n,5)**2
! 	         strsGrid(n,4) = strsGrid(n,4) - histGrid(n,3)*histGrid(n,4)
! 	         strsGrid(n,5) = strsGrid(n,5) - histGrid(n,3)*histGrid(n,5)
! 	         strsGrid(n,6) = strsGrid(n,6) - histGrid(n,4)*histGrid(n,5)
!	         pSum = pSum + histGrid(n,1)
	      endif
	   enddo

!	   pSum = pSum/hSize
 	   do n = 1, hSize
   	      histGrid(n,1) = histGrid(n,1)/cbvolm
	      do j = 1, 6
	         strsGrid(n,j) = -strsGrid(n,j)*histGrid(n,1)
	      enddo
	      strsGrid(n,7) = -(strsGrid(n,1) + strsGrid(n,2) + strsGrid(n,3))/3.0
 	   enddo
	endif

	return

	end

!----------------------------------------------------------------------------------

	subroutine PrintLocalProps(Sampletime)

	implicit none
	include	'dpdflow.h'

!	parameters

	real	Sampletime

!	Locals

	integer	i, j, k, l, n
	real*8	xval, yval, zval

	write(lst1, '(/5x, ''Local Properties of Flow Field'')')
	write(lst2, '(/5x, ''Local Properties on Stress Field'')')
	write(lst1, '(10x, ''time = '', f9.4)') Sampletime
	write(lst2, '(10x, ''time = '', f9.4)') Sampletime
	write(lst1, 100)
	write(lst2, 110)

	do i = 1, sizeHistGrid(1)
	   xval = (i - 0.5)*region(1)/sizeHistGrid(1) - regionH(1)
	   do j = 1, sizeHistGrid(2)
          yval = (j - 0.5)*region(2)/sizeHistGrid(2) - regionH(2)
          do k = 1, sizeHistGrid(3)

!            n = i + (k - 1)*sizeHistGrid(1)
             n = i + ((k - 1)*sizeHistGrid(2) + j - 1) * sizeHistGrid(1)
	         zval = (k - 0.5)*region(NDIM)/sizeHistGrid(3) - &
			         regionH(NDIM)
!	         zval = (k - 0.5)*(initUcell(3)*gap(3))/sizeHistGrid(3) - &
!			        initUcell(3)/2*gap(3)
	         write(lst1,120) n, xval, yval, zval, histGrid(n,2), histGrid(n,3), &
	               histGrid(n,4), histGrid(n,5), histGrid(n,6), histGrid(n,1) 
	         write(lst2,130) n, xval, yval, zval, (strsGrid(n,l), l = 1, 7)
          enddo
	   enddo
	enddo

100	format(4x, 'No', 5x, 'X', 8x, 'Y', 8x, 'Z', 9x, ' T', 9x, 'Vx', 12x, &
	       'Vy', 12x, 'Vz', 6x, 'V_fraction', 2x, 'rho')
110	format(4x, 'No', 6x, 'X', 8x, 'Y', 8x, 'Z', 7x, 'Sxx', 9x, 'Syy', 9x, &
	       'Szz', 9x, 'Sxy', 9x, 'Sxz', 9x, 'Syz',7x,'pressure')
120	format(1x, i7, 3f9.4, f10.4, 3e14.6, 2f10.5)
130	format(1x, i7, 3f9.4, 7e12.4)

	return

	end
