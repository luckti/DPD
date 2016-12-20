!----------------------------------------------------------
!   Apply the periodical boundary condition on paritcles
!----------------------------------------------------------

	subroutine ApplyBoundaryCond

	implicit none
	include	'dpdflow.h'

!	parameters
!	Locals

	integer	k, n, vSign
	real*8 	e(NDIM)

!	skip over wall atoms

 	do n = nWallAtom+1, nAtom
       do k = 1, PNDIM
	      if(r(n,k) .gt. region(k) .or. r(n,k) .lt. -region(k)) then
	         print*, 'n = ', n, ' k = ', k, ' stepCount =', stepCount
	         print*, 'r(n,k) = ', r(n,k), '  region(k) = ', region(k)
	         print*, 'rv(n,k) = ', rv(n,k)
	         print*, 'molcule moves too far, reduce time stepCount'
	         pause
	      endif

	      if(r(n,k) .ge. regionH(k)) then
	         r(n,k) = r(n,k) - region(k)

! record the times of the chain particle crossing the boundary (DP--22/11/11)
             ncc(n,k) = ncc(n,k) + 1
                
	      elseif(r(n,k) .lt. - regionH(k)) then         
	         r(n,k) = r(n,k) + region(k)

! same as above, cc: Chain Crossing (DP--22/11/11)
             ncc(n,k) = ncc(n,k) - 1

	      endif
	   enddo
	enddo

!DP	do n = 1, nWallAtom
!	   do k = 1, NDIM - 1
!	      if(r(n,k) .ge. regionH(k)) then
!	         r(n,k) = r(n,k) - region(k)
!	      elseif(r(n,k) .lt. - regionH(k)) then
!	         r(n,k) = r(n,k) + region(k)
!	      endif
!	   enddo
!DP	enddo

	return

	end
