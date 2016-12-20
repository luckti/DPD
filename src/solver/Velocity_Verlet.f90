!----------------------------------------------------------
!   Velocity-Verlet algorithm
!----------------------------------------------------------

	subroutine Velocity_Verlet

	implicit none
	include	'dpdflow.h'

	integer	i, k, n, ni, nk 
	real*8 	temp, vn, e(NDIM), vr(NDIM)
	real*8, allocatable :: rvm(:,:)

	allocate (rvm(nAtom,3))

    do k = 1, NDIM
	   do n = nWallAtom + 1, nAtom
          temp = deltaT*ra(n,k)
	      rvm(n,k) = rv(n,k)
	      r(n,k)   = r(n,k)   + deltaT*(rv(n,k) + 0.5*temp)
	      rv(n,k)  = rv(n,k)  + lambda*temp
	      rvm(n,k) = rvm(n,k) + 0.5*temp
	   enddo
	enddo
    
  	call ApplyBoundaryCond

	do n = nWallAtom + 1, nAtom
	   if(atomID(n) .gt. 0) then
	      call RandVec3(e)
          if(stepCount.le.stepEquil) then
          	 do k = 1, NDIM
				rv(n,k)  = 0.5*vMag*e(k)
				rvm(n,k) = rv(n,k)
			 enddo
		  else
		     rv(n,1)  = 1.2*vMag*e(1) + r(n,3)*shearRate
!            rv(n,1)  = 1.2*vMag*e(1) + initUcell(3)*gap(3)/2.*shearRate
		     rvm(n,1) = rv(n,1)
		     rv(n,2)  = 1.2*vMag*e(2)
		     rvm(n,2) = rv(n,2)
  		     rv(n,3)  = 1.2*vMag*e(3)
		     rvm(n,3) = rv(n,3)
	      endif
       endif
    enddo  
!    if(stepCount .eq. 798) print*, 'cf0'     !test

	call ComputeForces
    call DropIntForce
!    if(stepCount .eq. 798) print*, 'cf0'     !test
    do k = 1, NDIM
       do n = nWallAtom + 1, nAtom
          ra(n,k) = raCV(n,k) + raCR(n,k) + raDP(n,k) + raRD(n,k) +  &
                    raSP(n,k)
       enddo
    enddo

    if(nDp .gt. 0) then
       if(stepCount .eq. stepEquil) call forceRec(0)
       if(stepCount .gt. stepEquil) then
          call forceRec(1)
          if(mod((stepCount - stepEquil),500).eq.0) then
             call forceRec(2)
             call forceRec(0)
          endif
       endif
       do n = nWallAtom + 1, nDpEnd
          i = (n - nWallAtom - 1) / nPDP + 1
          if(DpFzn(i) .eq. 1) ra(n,:) = 0.
       enddo

! ----- froze fluid particles and the droplets particles         -----
! ----- consequently to keep the initial formation of droplet    -----
! ----- (DP--21/07/12)                                           -----
       if (nDp .gt.0)then
           if(stepCount .le. stepEquil/3) then      !froze the fluid particles 
              do n = nDpEnd + 1, nAtom
                 do k = 1, NDIM
                    ra(n,k)  = 0.
                    rv(n,k)  = 0.
                    rvm(n,k) = 0.
                 enddo
              enddo
           elseif(stepCount .le. stepEquil*2/3) then        !froze the droplets particles
              do n = nWallAtom + 1, nDpEnd
                 do k = 1, NDIM
                    ra(n,k)  = 0.
                    rv(n,k)  = 0.
                    rvm(n,k) = 0.
                 enddo
              enddo           
           endif
       endif

!       if(stepCount .ge. stepEquil - 1 .and. stepCount .le. stepFzn) &
!       then
!          do n = nWallAtom + 1, nDpEnd
!             do k = 1, NDIM
!                ra(n,k)  = 0.
!                rv(n,k)  = 0.
!                rvm(n,k) = 0.
!             enddo
!          enddo
!       endif
    endif
    


!DP if(stepCount .ge. stepEquil) call ComputeExternalForce
    call ComputeExternalForce

	do n = nWallAtom + 1, nAtom
       if( atomID(n)>nAtom) atomID(n)=0        
       i = atomID(n)

	   if(i .le. 0) then
	      do k = 1, NDIM
	         rv(n,k) = rvm(n,k) + 0.5*deltaT*ra(n,k)
	      enddo
	   else
	      call RandVec3(e)
	      vn = 0.
          if(stepCount .ge. stepEquil) then                   
!v0
             vr(1) = 1.2*vMag*e(1) + r(n,3)*shearRate          
!v1
!            vr(1) = 1.2*vMag*e(1) + shearRate*initUcell(3)/2*gap(3)
          else                                                 
             vr(1) = 1.3*vMag*e(1)                           
          endif                                             
    
          vn = vn + vr(1)*wn(i,1)
	      vr(2) = 1.2*vMag*e(2)
	      vn = vn + vr(2)*wn(i,2)
  	      vr(3) = 1.2*vMag*e(3)
	      vn = vn + vr(3)*wn(i,3)

	      do k = 1, NDIM
 	         rv(n,k) = vr(k) + (abs(vn) - vn)*wn(i,k)          
	      enddo
	   endif
    enddo
    
	deallocate (rvm)

	return

	end

!---------------------------------------------------------------------------

    subroutine forceRec(icord)
    
    implicit none
    include 'dpdflow.h'

	integer, intent(in):: icord
    real*8 frcTot(10,NDIM), frcCV(10,NDIM), frcDP(10,NDIM), frcRD(10,NDIM)
    integer n, i, k

    if(icord .eq. 0) then
       frcTot(:,:) = 0.
       frcCV(:,:)  = 0.
       frcDP(:,:)  = 0.
       frcRD(:,:)  = 0.
    elseif(icord .eq. 1) then
       do n = nWallAtom + 1, nDpEnd
          i = int((n - nWallAtom - 1) / nPDP) + 1
          do k = 1, NDIM
             frcTot(i,k) = frcTot(i,k) + ra(n,k)
             frcCV(i,k)  = frcCV(i,k)  + raCV(n,k)
             frcDP(i,k)  = frcDP(i,k)  + raDP(n,k)
             frcRD(i,k)  = frcRD(i,k)  + raRD(n,k)
          enddo
       enddo
    elseif(icord .eq. 2) then
       do i = 1, nDp
          do k = 1, NDIM
             frcTot(i,k) = frcTot(i,k) / 500.
             frcCV(i,k)  = frcCV(i,k)  / 500.
             frcDP(i,k)  = frcDP(i,k)  / 500.
             frcRD(i,k)  = frcRD(i,k)  / 500.
             enddo
          write(forc, '(2x, f7.2, 12f12.4)') timeNow, frcTot(i,1:NDIM),&
                frcCV(i,1:NDIM), frcDP(i,1:NDIM), frcRD(i,1:NDIM)
       enddo
    endif
        
    end

!--------------------------------------------------------------------

	subroutine Velocity_Verlet0

	implicit none
	include	'dpdflow.h'

	integer	k, n, i
	real*8 	temp
	real*8, allocatable :: rvm(:,:)

	allocate (rvm(nAtom,3))

	do k = 1, NDIM
	   do n = nWallAtom + 1, nAtom
	      temp = deltaT*ra(n,k)
	      rvm(n,k) = rv(n,k)
	      r(n,k) = r(n,k) + deltaT*(rv(n,k) + 0.5*temp)
	      rv(n,k) = rv(n,k) + lambda*temp
	      rvm(n,k) = rvm(n,k) + 0.5*temp
	   enddo
	enddo

  	call ApplyBoundaryCond

! -----test Poisson ratio ------- (DP--16/03/12)
    if(stepCount .ge. stepEquil) then
       if(nDP .gt. 0) then
! ----- to give the velocity of particles on some sepcified planes --
          if(stepCount .le. stepStop) call VelPlane (0)
       else
          call VelPlane (0)
       endif
    endif

	call ComputeForces

    do k = 1, NDIM
       do n = nWallAtom + 1, nAtom
          ra(n,k) = raCV(n,k) + raCR(n,k) + raDP(n,k) + raRD(n,k) +  &
                    raSP(n,k)
       enddo
    enddo

    if(nDp .gt. 0) then
       if(stepCount .eq. stepEquil) call forceRec(0)
       if(stepCount .gt. stepEquil) then
          call forceRec(1)
          if(mod((stepCount - stepEquil),500).eq.0) then
             call forceRec(2)
             call forceRec(0)
          endif
       endif

       do n = nWallAtom + 1, nDpEnd
          i = (n - nWallAtom - 1) / nPDP + 1
          if(DpFzn(i) .eq. 1) ra(n,:) = 0.
       enddo

! ----- frozen fluid particles or the droplets particles         -----
! ----- consequently to keep the initial formation of droplet    -----
! ----- (DP--21/07/12)                                           -----
       if(stepCount .lt. stepEquil*0.65) then
          do n = nWallAtom + 1, nDpEnd
             do k = 1, NDIM
                ra(n,k) = 0.
                rv(n,k) = 0.
                rvm(n,k)= 0.
             enddo
          enddo
       elseif(stepCount .le. stepEquil) then
          do n = nDpEnd + 1, nAtom
             do k = 1, NDIM
                ra(n,k) = 0.
                rv(n,k) = 0.
                rvm(n,k)= 0.
             enddo
          enddo
       endif

       if(stepCount .ge. stepEquil - 1) then
          do n = nWallAtom + 1, nDpEnd
             i = (n - nWallAtom - 1) / nPDP + 1
             if(DpFzn(i) .eq. 1) then
                do k = 1, NDIM
                   ra(n,k) = 0.
                   rv(n,k) = 0.
                   rvm(n,k)= 0.
                enddo
             else
                if(stepCount .le. stepFzn) then
                   do k = 1, NDIM
                      ra(n,k) = 0.
                      rv(n,k) = 0.
                      rvm(n,k)= 0.
                   enddo
                endif
             endif
          enddo
       endif
    endif
	
! 	if(stepCount .ge. stepEquil) call ComputeExternalForce
    call ComputeExternalForce

	do k = 1, NDIM
	   do n = nWallAtom + 1, nAtom
	      rv(n,k) = rvm(n,k) + 0.5*deltaT*ra(n,k)
	   enddo
	enddo

! -----test Poisson ratio ------- (DP--16/03/12)
    if(stepCount .ge. stepEquil) then
       if(nDP .gt. 0) then
! ----- to give the velocity of particles on some sepcified planes ---
          if(stepCount .le. stepStop) call VelPlane (0)
       else
          call VelPlane (0)
       endif
    endif

	deallocate (rvm)

	return

	end 

!---------------------------------------------------------------------------

    subroutine VelPlane (icord)

    implicit none
    include 'dpdflow.h'

	integer, intent(in):: icord
    integer n, k
    integer npb, npt, npf, npk, npl

    real*8 Px1, Px2, Py1, Py2, Pz1, Pz2
    real*8 deltaX, deltaY, deltaZ

!    real*8 Vx1, Vx2, Vz1, Vz2
    real*8 vpb(NDIM), vpt(NDIM), vpf(NDIM), vpk(NDIM), vpl(NDIM)

    open(31, file = './data/elongCnf.dat', status='unknown')

    read(31, *)
    read(31, *) Px1, Px2, deltaX
    read(31, *) Py1, Py2, deltaY
    read(31, *) Pz1, Pz2, deltaZ

    close(31)

!   Px1 = -30.
!   Px2 =  30.

!   Py1 = -10.
!   Py2 =  10.

!   Pz1 = -20.
!   Pz2 =  20.

!   deltaX = 0.25
!   deltaY = 0.25
!   deltaZ = 0.25

!   Vx1 =  1.0
!   Vx2 = -1.0

!   Vz1 =  1.0
!   Vz2 = -1.0

    do n = nWallAtom + 1, nAtom
!---- the left plane (r=-X+/-delta)----
!       if(r(n,1) .le. Px1 + deltaX .and. r(n,1) .ge. Px1 - deltaX   &
!          .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1    &
!          .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
       if(r(n,1) .le. Px1 + deltaX .and. r(n,1) .ge. Px1 - deltaX   &
          .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
          rv(n,1) =  shearRate*r(n,1)
          rv(n,2) =  0.
          rv(n,3) = -shearRate*r(n,3)
       endif

!---- the right plane (r= X+/-delta)----
!       if(r(n,1) .le. Px2 + deltaX .and. r(n,1) .ge. Px2 - deltaX   &
!          .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1    &
!          .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
       if(r(n,1) .le. Px2 + deltaX .and. r(n,1) .ge. Px2 - deltaX   &
          .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
          rv(n,1) =  shearRate*r(n,1)
          rv(n,2) =  0.
          rv(n,3) = -shearRate*r(n,3)
       endif

!---- the bottom plane (r=-Z+/-delta)----
!       if(r(n,3) .le. Pz1 + deltaZ .and. r(n,3) .ge. Pz1 - deltaZ   &
!          .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1    &
!          .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1 ) then
       if(r(n,3) .le. Pz1 + deltaZ .and. r(n,3) .ge. Pz1 - deltaZ   &
          .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1 ) then
          rv(n,1) =  shearRate*r(n,1)
          rv(n,2) =  0.
          rv(n,3) = -shearRate*r(n,3)
       endif

!---- the top plane (r= Z+/-delta)----
!       if(r(n,3) .le. Pz2 + deltaZ .and. r(n,3) .ge. Pz2 - deltaZ   &
!          .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1    &
!          .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1 ) then
       if(r(n,3) .le. Pz2 + deltaZ .and. r(n,3) .ge. Pz2 - deltaZ   &
          .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1 ) then
          rv(n,1) =  shearRate*r(n,1)
          rv(n,2) =  0.
          rv(n,3) = -shearRate*r(n,3)
       endif
    enddo

! -- to record the particles' velocity on top, bottom, front and back planes
    if(icord .eq. 1) then
       npb = 0
       npt = 0
       npf = 0
       npk = 0
       npl = 0

       vpb(:) = 0.
       vpt(:) = 0.
       vpf(:) = 0.
       vpk(:) = 0.
       vpl(:) = 0.

       do n = nWallAtom + 1, nAtom
!---- the front plane (r=-Y+/-delta)----
          if(r(n,2) .le. Py1 + deltaY .and. r(n,2) .ge. Py1 - deltaY   &
             .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1    &
             .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
             npf = npf + 1
             do k = 1, NDIM
                vpf(k) = vpf(k) + rv(n,k)
             enddo
          endif

!---- the back plane (r= Y+/-delta)----
          if(r(n,2) .le. Py2 + deltaY .and. r(n,2) .ge. Py2 - deltaY   &
             .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1    &
             .and. r(n,3) .le. Pz2 .and. r(n,3) .ge. Pz1 ) then
             npk = npk + 1
             do k = 1, NDIM
                vpk(k) = vpk(k) + rv(n,k)
             enddo
          endif

!---- the bottom plane (r=-Z+/-delta)----
          if(r(n,3) .le. Pz1 + deltaZ .and. r(n,3) .ge. Pz1 - deltaZ   &
             .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1    &
             .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1 ) then
             npb = npb + 1
             do k = 1, NDIM
                vpb(k) = vpb(k) + rv(n,k)
             enddo
          endif

!---- the top plane (r= Z+/-delta)----
          if(r(n,3) .le. Pz2 + deltaZ .and. r(n,3) .ge. Pz2 - deltaZ   &
             .and. r(n,1) .le. Px2 .and. r(n,1) .ge. Px1    &
             .and. r(n,2) .le. Py2 .and. r(n,2) .ge. Py1 ) then
             npt = npt + 1
             do k = 1, NDIM
                vpt(k) = vpt(k) + rv(n,k)
             enddo
          endif
       enddo

       vpb(:) = vpb(:) / npb
       vpt(:) = vpt(:) / npt
       vpf(:) = vpf(:) / npf
       vpk(:) = vpk(:) / npk

!       write(lst6, '(1x, f8.4, 12f10.4)') timeNow, vpf(:), vpk(:), vpb(:), vpt(:)
    endif

    end