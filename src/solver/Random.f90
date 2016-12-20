!----------------------------------------------------------
!   Generate random number as required
!----------------------------------------------------------

    subroutine ranils(iseed)

!	Choice if iseed: 0<= iseed <= 2E+9;
	
	parameter (in = 214783563, ik = 40014,  & 
     		   iq = 53668,     ir = 12211,  &
     		   ntab = 32)
	integer	iv(ntab)
	common	/ranbls/ idum, idum2, iy, iv

!	Initial seeds for two random number generators

	idum = iseed + 123456789
	idum2 = idum

!	Load the shuffle table (after 8 warm-ups)

	do j = ntab+8, 1, -1
	   k = idum/iq
	   idum = ik*(idum - k*iq) - k*ir
	   if(idum .lt.  0) idum = idum + in
	   if(j .lt. ntab) iv(j) = idum
	enddo

	iy = iv(1)

	return

	end

!-------------------------------------------------------------------

	real function ranuls()

	parameter (in1 = 2147483563, ik1 = 40014, & 
     		   iq1 = 53668,      ir1 = 12211, &
     		   in2 = 2147483399, ik2 = 40692, &
     		   iq2 = 52774,      ir2 = 3791,  &
     	 	   ntab = 32,	     an = 1./in1, &
     		   inm1 = in1 - 1,   ndiv = 1 + inm1/ntab)
	integer	iv(ntab)
	common	/ranbls/ idum, idum2, iy, iv

!	Linear congruential generactor 1

	k = idum/iq1
	idum = ik1*(idum - k*iq1) - k*ir1
	if(idum .lt. 0) idum=idum+in1

!	Linear congruential generator 2

	k = idum2/iq2
	idum2 = ik2*(idum2 - k*iq2) - k*ir2
	if(idum2 .lt. 0) idum2 = idum2 + in2

!	Shuffling and subtracting

	j = 1 + iy/ndiv
	iy = iv(j) - idum2
	iv(j) = idum
	if(iy .lt. 1) iy = iy + inm1
	ranuls = an*iy

	return

	end

!--------------------------------------------------------

	double precision function rangls()

	save	iflag, gauss2
	data	iflag / 0 /

	if(iflag .eq. 0) then

10	   continue

!	Pair of uniform random numbers in [-1,1]x[-1,1]

	   x1 = 2.*ranuls() - 1.
	   x2 = 2.*ranuls() - 1.

!	if not in the unit circle, try again

	   xsq = x1*x1 + x2*x2
	   if(xsq .ge. 1. .or. xsq .eq. 0.) goto 10

!	pair of Gaussian random numbers; return one and 
!	save the other for next time

	   aux = sqrt(-2.*log(xsq)/xsq)
	   rangls = x1*aux
	   gauss2 = x2*aux
	   iflag = 0

	else

	   rangls = gauss2
	   iflag = 0

	endif

	return

	end

!------------------------------------------------------------------

	function RandR (seed)
	
	implicit none

!	Parameters
	
	integer, intent(in out):: seed
	real*4 	RandR

!	Locals

	integer	IMUL, IADD, MASK
	real*8	SCAL
	parameter (IMUL = 314159269, &
     		   IADD = 453806245, &
     		   MASK = 2147483647, &
     		   SCAL = 0.4656612873e-9)

	seed = mod((seed*IMUL + IADD), MASK)
	RandR = seed*SCAL

	return

	end

!----------------------------------------------------------------

	subroutine RandVec3 (p)

	implicit none

!	Parameters

	real*8, intent(out):: p(3)

!	Locals

	real*8	x, y, s, r1, r2
	real*4	ranuls

 	s = 2.0
 	do while (s .gt. 1.)
  	   r1 = ranuls()
  	   r2 = ranuls()
	   x = 2.*r1 - 1.0
	   y = 2.*r2 - 1.0
	   s = x*x + y*y
 	enddo

	p(3) = 1. - 2.*s
	s = 2.*sqrt(1. - s)
	p(1) = s*x
	p(2) = s*y

	return

	end