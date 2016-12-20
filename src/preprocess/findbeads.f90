!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

    subroutine findbeads(m, idom, ndpsb, rc, RdsDp, chainLen, subreg, seeds, fgap, subdom, x)

    implicit none

!   parameters

    integer	MAXCTP, MAXBDS, NDIM
    parameter (MAXCTP = 10, MAXBDS = 90000, NDIM = 3)
    integer	 m, chainLen, nChainend, seeds, ndpsb(MAXCTP)
    real*8	subdom(MAXCTP,8), subreg(MAXCTP,NDIM), x(MAXBDS,NDIM), fgap
    real*8  RdsDp

!	Locals

	integer MAXATOM, MAXWATM, i, k, idom, mc, mp, mk
	real*4 RandR
	real*8	d, s, y(NDIM), shift(NDIM), e(NDIM), rc(MAXCTP, 5000, NDIM)
	parameter (MAXWATM = 50000, MAXATOM = 400000)

	d = 0.5*fgap

!	find a position for the first bead of a Chain

 5 	y(1) = subdom(idom,1) + 0.5*(1.0 + RandR(seeds))*subreg(idom,1)
	y(2) = 0.5*RandR(seeds)*subreg(idom,2) 
 	y(3) = subdom(idom,2) + 0.5*(1.0 + RandR(seeds))*subreg(idom,3)

	do i = 1, m
	   s = 0.
	   do k = 1, 3
	      s = s + (y(k) + shift(k) - x(i,k))**2
	   enddo
	   if(s .lt. 0.25) goto 5
	enddo

    do i = 1, ndpsb(idom)
       s = 0.
       do k = 1, 3
          s = s + (y(k) + shift(k) - rc(idom, i, k))**2
	   enddo
!	   if(s .lt. (RdsDp + 0.25)**2) goto 5
	   if(s .lt. (RdsDp + 0.50)**2) goto 5
    enddo

    if((y(1) - subdom(idom,1)) .lt. 0.2 .or. &
	   (subdom(idom,3) - y(1)) .lt. 0.2) goto 5
	if((y(3) - subdom(idom,2)) .lt. 0.2 .or. &
	   (subdom(idom,8) - y(3)) .lt. 0.2) goto 5

	mp = m
	m = m + 1
	do k = 1, NDIM
	   x(m,k) = y(k)
	enddo

	mc = 1
    mk = 1
	do while(mc .lt. ChainLen)

15	   do k = 1, NDIM
	      shift(k) = 0.
	   enddo

	   call RandVec3(e, seeds)

	   do k = 1, NDIM
	      y(k) = x(m,k) + d*e(k)
	   enddo

	   if(y(2) .lt. -0.5*subreg(idom,2)) then
	      y(2) = y(2) + subreg(idom,2)
	      shift(2) = -subreg(idom,2)
	   endif

	   if(y(2) .gt. 0.5*subreg(idom,2)) then
	      y(2) = y(2) - subreg(idom,2)
	      shift(2) = subreg(idom,2)
	   endif

	   if((y(1) - subdom(idom,1)) .lt. 0.2 .or. &
	      (subdom(idom,3) - y(1)) .lt. 0.2) goto 15
	   if((y(3) - subdom(idom,2)) .lt. 0.2 .or. &
	      (subdom(idom,8) - y(3)) .lt. 0.2) goto 15

	   do i = 1, m
	      s = 0.
	      do k = 1, 3
	         s = s + (y(k) + shift(k) - x(i,k))**2
	      enddo
 	      if(i .le. mp) then
	         if(s .lt. 0.35) goto 15
             mk = mk + 1
             if(mk .gt. 10e6) goto 5
	      else
	         if(s .lt. 0.1) goto 15
             mk = mk + 1
             if(mk .gt. 10e6) goto 5
	      endif
	   enddo

       do i = 1, ndpsb(idom)
          s = 0.
          do k = 1, NDIM
             s = s + (y(k) + shift(k) - rc(idom, i, k))**2
	      enddo
!	      if(s .lt. (RdsDp + 0.25)**2) goto 15
	      if(s .lt. (RdsDp + 0.25)**2) goto 5

! the center of two neighbour beads should locate out of the drop domain (DP--27/12/11)
!          if(mc .gt. 1) then
!             s = 0.
!             do k = 1, NDIM
!                s = s + (0.5*(y(k) + x(m,k)) + shift(k) - rc(idom, i, k))**2
!             enddo
!	         if(s .lt. (RdsDp + 0.25)**2) goto 15
!          endif

       enddo

	   m = m + 1
	   mc = mc + 1
	   do k = 1, NDIM
	      x(m,k) = y(k)
	   enddo
	enddo

	return

	end