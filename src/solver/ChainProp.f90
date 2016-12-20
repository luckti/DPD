!----------------------------------------------------------
!   Compute chain forces between beads and extract the 
!   chain properties
!----------------------------------------------------------

	subroutine FENEChainForce
	
	implicit none
	include 'dpdflow.h'

!	Parameters

!	Locals

	integer	j, k, m1, m2, n
	real*8 	rr, rij, fenef, f, fr, dr(NDIM), dcc(NDIM)

!	compute the spring forces of the chains

	m1 = nWallAtom
	do n = 1, nChain
	   do j = 1, ChainLen - 1 
	      m1 = m1 + 1
	      m2 = m1 + 1
	      rr = 0.
	      do k = 1, NDIM
!DP          dr(k) = r(m2,k) - r(m1,k)
!DP	         if(k .lt. NDIM) dr(k) = dr(k) - region(k)*anint(dr(k)/region(k))
             if(stepCount .eq. 1) then
                dcc(k) = ncc(m2,k) - ncc(m1,k)
                dr(k) = r(m2,k) + dcc(k)*region(k) - r(m1,k)
                if(dr(k) .gt. regionH(k)) ncc(m2,k) = ncc(m2,k) - 1
                if(dr(k) .lt. - regionH(k)) ncc(m2,k) = ncc(m2,k) + 1
             endif
             dcc(k) = ncc(m2,k) - ncc(m1,k)                       ! differences of crossing times (DP--22/11/11)
             dr(k)  = r(m2,k) + dcc(k)*region(k) - r(m1,k)        ! (DP--22/11/11)
	         rr = rr + dr(k)*dr(k)
	      enddo
	      rij = sqrt(rr)
	      do k = 1, NDIM
	         dr(k) = dr(k)/rij
	      enddo
	      rr = (rij - reqfene)**2

! for FENE chain-------------------------- (DP--21/10/11)
!	      fenef = Hfene*(rij - reqfene)/(1. - rr/rrfene)
!	      if(rij .gt. rmaxfene .or. rij .lt. 2.*reqfene - rmaxfene) then
!	         print*, ' rij = ', rij, m1, m2
!	         print*, 'rv = ', rv(m1,1), rv(m1,2), rv(m1,3)
!	         fenef = 10.*Hfene*(rij - reqfene)
!	      endif
!-----------------------------------------

! for oldroyd-B fluid--------------------- (DP--21/10/11) modified for constant spring (Oldroyd-B model)
          fenef = Hfene*(rij - reqfene)
!-----------------------------------------

	      do k = 1, NDIM
	         f = fenef*dr(k)
	         raSP(m1,k) = raSP(m1,k) + f
	         raSP(m2,k) = raSP(m2,k) - f
!	         fr = 2.*f*rij*dr(k)
		     fr = f*rij*dr(k)
	         rforce(m1,k) = rforce(m1,k) - fr
	         rforce(m2,k) = rforce(m2,k) - fr
	      enddo
	      f = fenef*rij
	      fr = f*dr(1)*dr(2)
	      rforce(m1,4) = rforce(m1,4) - fr 
	      rforce(m2,4) = rforce(m2,4) - fr
	      fr = f*dr(1)*dr(3)
	      rforce(m1,5) = rforce(m1,5) - fr 
	      rforce(m2,5) = rforce(m2,5) - fr
	      fr = f*dr(2)*dr(3)
	      rforce(m1,6) = rforce(m1,6) - fr
	      rforce(m2,6) = rforce(m2,6) - fr
	   enddo
	   m1 = m1 + 1
	enddo

	return

	end

!------------------------------------------------------------------------

	subroutine EvalChainProps

	implicit none
	include	'dpdflow.h'

!	Paramers

!	Locals
	
	integer i, j, k, n, n1
	real*8	a1, a2, a3, dr, ee, sumXY, sumYZ, sumZX, t
	real*8	c(NDIM), g(6), gVal(3), shift(3), sumR(3), sumRR(3)

	countChainProps = countChainProps + 1
	if(countChainProps .eq. 1) then
	   aaDistSq = 0.
	   eeDistSq = 0.
	   radGyrSq = 0.
	   gMomRatio1 = 0.
	   gMomRatio2 = 0.
	endif

	n = nWallAtom
	do i = 1, nChain
	   ee = 0.
	   sumXY = 0.
	   sumYZ = 0.
	   sumZX = 0.
	  do k = 1, NDIM
	     shift(k) = 0.
	     sumR(k) = 0.
	     sumRR(k) = 0.
	  enddo

	  do j = 1, ChainLen
	     n = n + 1
	     if(j .gt. 1) then
	        do k = 1, NDIM
		       dr = r(n,k) - r(n-1,k)
		       if(abs(dr) .gt. regionH(k)) then
		          t = sign(region(k), dr)
		          shift(k) = shift(k) - t
		          dr = dr - t
		       endif
		       aaDistSq = aaDistSq + dr**2          ! distance square of the chain  (DP--23/11/11)
	        enddo
	     else
	       n1 = n
	     endif

	     do k = 1, NDIM
	        c(k) = r(n,k) + shift(k)
	        sumR(k) = sumR(k) + c(k)
	        sumRR(k) = sumRR(k) + c(k)**2
	     enddo

	     sumXY = sumXY + c(1)*c(2)
	     sumYZ = sumYZ + c(2)*c(3)
	     sumZX = sumZX + c(3)*c(1)

	  enddo

	  do k = 1, NDIM
	     ee = ee + (c(k) - r(n1,k))**2
	     sumR(k) = sumR(k) / ChainLen
	     g(k) = sumRR(k) / ChainLen - sumR(k)**2
	  enddo
	  eeDistSq = eeDistSq + ee
	  g(4) = sumXY / ChainLen - sumR(1)*sumR(2)
	  g(5) = sumZX / ChainLen - sumR(3)*sumR(1)
	  g(6) = sumYZ / ChainLen - sumR(2)*sumR(3)
	  a1 = - g(1) - g(2) - g(3)
	  a2 = g(1)*g(2) + g(2)*g(3) + g(3)*g(1) - g(4)**2 - g(5)**2 -g(6)**2
	  a3 = g(1)*g(6)**2 + g(2)*g(5)**2 + g(3)*g(4)**2 - 2.*g(4)*g(5)*g(6) - &
	       g(1)*g(2)*g(3)

	  call SolveCubic(a1, a2, a3, g)

	  gVal(1) = max(g(1), g(2), g(3))                     ! maximum 
	  gVal(3) = min(g(1), g(2), g(3))                     ! minimum
	  gVal(2) = g(1) + g(2) + g(3) - gVal(1) - gVal(3)    ! medium
	  radGyrSq = radGyrSq + gVal(1) + gVal(2) + gVal(3)
	  gMomRatio1 = gMomRatio1 + gVal(2) / gVal(1)
	  gMomRatio2 = gMomRatio2 + gVal(3) / gVal(1)

	enddo

	if(countChainProps .eq. limitChainProps) then
	   aaDistSq = aaDistSq / (nChain*(ChainLen - 1)*limitChainProps)
	   eeDistSq = eeDistSq / (nChain*limitChainProps)
	   radGyrSq = radGyrSq / (nChain*limitChainProps)
	   gMomRatio1 = gMomRatio1 / (nChain*limitChainProps)
	   gMomRatio2 = gMomRatio2 / (nChain*limitChainProps)

	   call PrintChainProps
	   countChainProps = 0
	endif

	return

	end

!-------------------------------------------------------------------

	subroutine PrintChainProps

	implicit none
	include 'dpdflow.h'
	
!	parameters

!	Locals

	integer id
	data 	id / 0 /

	if(id .eq. 0) write(lst3, 100)
	write(lst3, 110) stepCount, sqrt(aaDistSq), eeDistSq, radGyrSq, &
			gMomRatio1, gMomRatio2 
	id = 1

 100	format(1x, 'steps', 8x, '<l>', 7x, '<R^2>', 7x, '<S^2>', 5x, &
		'<g_2/g_1>', 3x, '<g_3/g_1>')
 110	format(i6, 5f12.4)
 
	return

	end

!-------------------------------------------------------------

	subroutine PrintChainConf

	implicit none
	include 'dpdflow.h'

!	parameters

!	Locals

	integer	m, n, i
	do n = 1, nChain, nChainConf
!DP	   m = (n - 1)*ChainLen + nWallAtom
       m = (n - 1)*ChainLen + nDpEnd
	   write(lst3, '(3x, ''Chain No. :'', i4, 2x, ''time = '', f9.2)') n,timeNow
	   do i = 1, ChainLen
	      m = m + 1	
	      write(lst3,'(1x, 3f10.5)') r(m,1),r(m,2),r(m,3)
	   enddo
	enddo

	return

	end

!-------------------------------------------------------------------------------

	subroutine ChainLength

	implicit none
	include 'dpdflow.h'

!	Parameters

!	Locals

	integer	i, j, k, l, m1, m2
	real*8 	rr, ssMax, tt, rc(NDIM), dr(NDIM), shift(NDIM), shifr(NDIM)

	m1 = nWallAtom
	do i = 1, nChain

	   do k = 1, NDIM
	      rc(k) = r(m1+1,k)
	   enddo

	   ssMax = 0.
	   do j = 1, ChainLen - 1
	      do k = 1, NDIM
	         shifr(k) = 0.
	      enddo
	      m1 = m1 + 1
	      do k = 1, NDIM
	         tt = r(m1+j,k) - r(m1+j-1,k)
	         if(k .lt. NDIM) &
		        shifr(k) = shifr(k) - region(k)*anint(tt/region(k))
	         rc(k) = rc(k) + r(m1+j,k) + shifr(k)
	         shift(k) = 0.
	      enddo
	      do l = j + 1, ChainLen
	         m2 = m1 + l - j
	         rr = 0.
	         do k = 1, NDIM
		        tt = r(m2,k) - r(m2-1,k)
		        if(k .lt. NDIM) &
		           shift(k) = shift(k) - region(k)*anint(tt/region(k))
	            dr(k) = r(m2,k) + shift(k) - r(m1,k)
	            rr = rr + dr(k)*dr(k)
	         enddo
	         if(rr .gt. ssMax) ssMax = rr
	      enddo
	   enddo

!	 the length of chain i

	   sChain(i) = sqrt(ssMax) 

!	 the centre position of chain i

	   do k = 1, NDIM
  	      chainCentre(i,k) = rc(k)/ChainLen
	   enddo
	   m1 = m1 + 1
	enddo 

	return

	end
 
!----------------------------------------------------------------------------

	subroutine GridAvChainProps(opcode)
	
	implicit none
	include 'dpdflow.h'

!	Parameters
	
	integer, intent(in) :: opcode

!	Locals

	integer	c, i, n

	if(opCode .eq. 0) then

	   do n = 1, hsize
	      GridChainLen(n,1) = 0.
	      GridChainLen(n,2) = 0.
	   enddo

	elseif(opcode .eq. 1) then

	   call ChainLength 

	   do i = 1, nChain
!DP	      c = int((chainCentre(i,3) + regionH(3))*sizeHistGrid(2)/ &
!		  region(3))*sizeHistGrid(1) + int((chainCentre(i,1) + &
!DP		  regionH(1))*sizeHistGrid(1)/region(1)) + 1

          c = (int((chainCentre(i,3) + initUcell(3)/2*gap(3))*sizeHistGrid(3)/  &
               initUcell(3)/gap(3))*sizeHistGrid(2) +   &
               int((chainCentre(i,2) + regionH(2))*sizeHistGrid(2)/region(2)))* &
               sizeHistGrid(1) + int((chainCentre(i,1) + regionH(1))*sizeHistGrid(1)/region(1)) + 1

	      GridChainLen(c,1) = GridChainLen(c,1) + sChain(i)
	      GridChainLen(c,2) = GridChainLen(c,2) + 1.0
	   enddo

	elseif(opcode .eq. 2) then

	   do n = 1, hsize
	      if(GridChainLen(n,2) .gt. 1.0) &
	         GridChainLen(n,1) = GridChainLen(n,1)/GridChainLen(n,2)
	   enddo

	endif

	return

	end

!---------------------------------------------------------------------------

	subroutine PrintChainLen

	implicit none
	include 'dpdflow.h'

!	Parameters

!	Locals
	
	integer i, j, k, n
	real*8	xval, yval, zval

	lst4 = 14
	write(lst4, '(/10x, ''Chain Length Dstribution'')')
	write(lst4, '(1x, ''No. of Chains = '', i4, 3x, ''No. of Beads per Chain = '', i4)') &
		  nChain, ChainLen
	write(lst4, '(1x, ''No. of Fluid Particles = '', i7)')  nFreeAtom
!DP	write(lst4, '(1x, ''Concentration of Suspesion'', f7.4)')  float(chainLen*nChain/nFreeAtom)
	write(lst4, '(1x, ''Concentration of Suspesion'', f7.4)')  float(chainLen*nChain)/nFreeAtom
	write(lst4,100)

	do i = 1, sizeHistGrid(1)
	   xval = (i - 0.5)*region(1)/sizeHistGrid(1) - regionH(1)
	   do j = 1, sizeHistGrid(2)
          yval = (j - 0.5)*region(2)/sizeHistGrid(2) -regionH(2)
          do k = 1, sizeHistGrid(3)

!	         n = i + (k - 1)*sizeHistGrid(1)
             n = i + ((k - 1)*sizeHistGrid(2) + j - 1) * sizeHistGrid(1)
!DP	         zval = (k - 0.5)*region(NDIM)/sizeHistGrid(2) - regionH(NDIM)
	         zval = (k - 0.5)*(initUcell(3)*gap(3))/sizeHistGrid(3) - &
			         initUcell(3)/2*gap(3)
	         write(lst4,110) n, xval, yval, zval, GridChainLen(n,1), int(GridChainLen(n,2))
          enddo
	   enddo
	enddo

100	format(//4x, 'No.', 5x, 'X', 8x, 'Y', 8x, 'Z', 5x, 'Chain Len', 6x, 'No. of Samples')
110	format(1x, i7, 3f9.4, e13.5, i7)

	return

	end

!-------------------------------------------------------------

	subroutine PrintDropConf

	implicit none
	include 'dpdflow.h'

!	parameters

!	Locals

	integer	m, n, i
   	do n = 1, nDp
	   m = (n - 1)*(nDpEnd - nWallAtom)/nDp + nWallAtom
	   write(lst5, '(3x, ''Drop No. :'', i4, 2x, ''time = '', f9.2)') n,timeNow
	   do i = 1, (nDpEnd - nWallAtom)/nDp
	      m = m + 1	
	      write(lst5,'(1x, 3f10.5, 3i4)') r(m,1:NDIM), ncc(m,1:NDIM)
	   enddo
	enddo

	return

	end

!-----------------------------------------------------------------------------------------------

	subroutine SolveCubic(a1, a2, a3, g)
	
	implicit none

!	Parameters

	real*8	a1, a2, a3, g(3)

!	Locals

	real*8	q1, q2, t, pi, ad3, q12
	data 	pi / 3.14159265358979 /

	ad3 = a1 / 3.0
	q1 = sqrt(a1*a1 - 3.0*a2) / 3.0
	q2 = (2.*a1*a1*a1 - 9.*a1*a2 + 27.*a3) / 54.
	q12 = 2.0*q1
	t = acos(q2 / (q1*q1*q1))

	g(1) = - q12*cos(t/3.0) - ad3
	g(2) = - q12*cos((t + 2.*pi) / 3.0) - ad3
	g(3) = - q12*cos((t + 4.*pi) / 3.0) - ad3
	
	return

	end

