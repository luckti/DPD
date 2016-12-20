!----------------------------------------------------------
!   Compute forces acting on particles
!----------------------------------------------------------

    subroutine ComputeForces

    implicit none
    include 'dpdflow.h'

!   Parameters
!   Locals

    integer c, i, j1, j2, k, m1, m2, m1X, m1Y, m1Z, m2X, m2Y, m2Z, n,  &
            offset
    integer	iofX(14), iofY(14), iofZ(14)
    real*8  f, fr, fcVal, fdVal, frVal, ftot, weight, rij, weightc,    &
            weightd, weightr, rr, rdvij, uVal, fstl, rangls, dn,       &
            gamma1, cigama1
    real*8  ri(NDIM), dr(NDIM), invWid(NDIM), shift(NDIM)

    integer, allocatable :: cellList(:),tag(:,:),ntag(:)
    real*8, allocatable :: sqrDist(:),mrho(:),ddr(:,:,:),drij(:,:),weightcd(:,:)
    real*8 fCV, fDP, fRD

!    integer nbSP(MAXATOM), nbDP(MAXATOM)
!    real*8 pct(MAXATOM)

!---MDPD used variables-----------------------------------------------
    integer j, m
    real*8  fcRval, weightp, fCR
    data    pi / 3.14159265358979 /
!---------------------------------------------------------------------

    data    iofX / 0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1 /
    data    iofY / 0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1 /
    data    iofZ / 0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1 /
 
    atomID = 0    
    allocate ( cellList(maxList), tag(nAtom, 150), ntag(nAtom))
    allocate ( sqrDist(nAtom) ,mrho(nAtom), ddr(nAtom, 150, NDIM), drij(nAtom, 150),weightcd(nAtom,150))

    do k = 1, NDIM
        invWid(k) = cells(k) / region(k)
    enddo

    do n = nAtom+1, maxList
        cellList(n) = 0
    enddo

    do n = nStartAtom, nAtom
        do k = 1, NDIM
            ra(n,k)   = 0.
            raCV(n,k) = 0.
            raCR(n,k) = 0.
            raDP(n,k) = 0.
            raRD(n,k) = 0.
            raSP(n,k) = 0.
            rforce(n,k) = 0.
            rforce(n,k+NDIM) = 0.
        enddo
        ntag(n) = 0
        mrho(n) = 0.

!       nbSP(n) = 0
!       nbDP(n) = 0
    enddo

    ! what for ?
    do n = nStartAtom, nAtom
        c = (int((r(n,3) + regionH(3))*invWid(3))*cells(2) +           &
             int((r(n,2) + regionH(2))*invWid(2)))*cells(1) +          &
             int((r(n,1) + regionH(1))*invWid(1)) + nAtom + 1
        
        if(c .lt. 1) then
            print*,regionH(1),regionH(2),regionH(3)
            print*,maxList
            print*, 'in cf, stepCount = ', stepCount
            print*, ' c = ', c, '  n = ', n
            print*, 'r(n,*) = ', (r(n,k), k = 1, 3)
            pause
        endif

        if(c .gt. maxList) then
            print*,regionH(1),regionH(2),regionH(3)
            print*,maxList
            print*, ' c = ', c, '  n = ', n, 'stepCount = ', stepCount
            print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
            pause
        endif

        cellList(n) = cellList(c)
        cellList(c) = n
        
    enddo
  
    
    uSum = 0.
    virSum = 0.

    do m1Z = 1, cells(3)
        do m1Y = 1, cells(2)
            do m1X = 1, cells(1)
                m1 = ((m1Z - 1)*cells(2) + m1Y - 1)*cells(1) + m1X+ nAtom

                do offset = 1, 14
                    m2X = m1X + iofX(offset)
                    shift(1) = 0

                    if(m2X .gt. cells(1)) then
                        m2X = 1
                        shift(1) =  region(1)
                    elseif(m2X .eq. 0) then
                        m2X = cells(1)
                        shift(1) = -region(1)
                    endif

                    m2Y = m1Y + iofY(offset)
                    shift(2) = 0

                    if(m2Y .gt. cells(2)) then
                        m2Y = 1
                        shift(2) =  region(2)
                    elseif( m2Y .eq. 0) then
                        m2Y = cells(2)
                        shift(2) = -region(2)
                    endif

                    m2Z = m1Z + iofZ(offset)
                    shift(3) = 0

                    if(runId .eq. 1) then
                        if(m2Z .eq. 0 .or. m2Z .gt. cells(3)) goto 15
!                    elseif(runId .eq. 4)then
!                        if(m2Z .eq. 0 .or. m2Z .gt. cells(3)) goto 15
                    else
                        if(m2Z .gt. cells(3)) then
                            m2Z = 1
                            shift(3) =  region(3)
                        elseif(m2Z .eq. 0) then
                            m2Z = cells(3)
                            shift(3) = -region(3)
                        endif
                    endif

                    m2 = ((m2Z - 1)*cells(2) + m2Y - 1)*cells(1) + m2X + nAtom

                    j1 = cellList(m1)
                    do while(j1 .gt. 0)
                        j2 = cellList(m2)
                        do while(j2 .gt. 0)
                            if(m1 .ne. m2 .or. j2 .lt. j1) then

                                do k = 1, NDIM
                                    dr(k) = r(j1,k) - r(j2,k) - shift(k)
                                enddo

                                rr = dr(1)*dr(1) + dr(2)*dr(2) +dr(3)*dr(3)

                                if(rr .lt. rrCut) then
                                    rij = sqrt(rr)

                                    rdvij = 0.
                                    do k = 1, NDIM
                                        dr(k) = dr(k)/rij
                                        rdvij = rdvij +dr(k)*(rv(j1,k) - rv(j2,k))
                                    enddo

                                    weight  = 1. - rij/rCut
                                    weightc = weight

! restrict the conservative cut radius equal to 1. (DP--7/12/11)
!                                   weightc = 1. - rij
!                                   weight  = 1. - rij/rcut

! inclusion region for repulsion potential in MDPD
                                    if(rij .lt. rCut2 .and. (j1 .gt. nDpEnd .or.j2 .gt. nDpEnd)) then
! weight function related to local density

                                        weightp  = 15./(2.*pi*r3Cut2)*(1-rij/rCut2)**2
                                        ntag(j1) = ntag(j1) + 1
                                        ntag(j2) = ntag(j2) + 1
                                        tag(j1,ntag(j1)) = j2
                                        tag(j2,ntag(j2)) = j1
                                        weightcd(j1,ntag(j1)) = 1. - rij/rCut2
                                        weightcd(j2,ntag(j2)) = weightcd(j1,ntag(j1))
                                     
                                        do k = 1, NDIM
                                            ddr(j1,ntag(j1),k) =  dr(k)
                                            ddr(j2,ntag(j2),k) = -dr(k)
                                        enddo                                       
                                        drij(j1,ntag(j1)) = rij                                        
                                        drij(j2,ntag(j2)) = rij
! local density parameter in MDPD
                                        mrho(j1) = mrho(j1) + weightp
                                        mrho(j2) = mrho(j2) + weightp
                                      
                                    endif

                                    if(weightc .lt. 0.d0) weightc=0.d0
                                    if(weight  .lt. 0.d0) weight =0.d0

! original version of DPD: W^C(r) = W^R(r) = 1 - r, W^D(r) = (1-r)**2
                                    weightd = weight*weight
                                    weightr = weight

! modified version of DPD: W^C(r) = 1-r, W^R = (1-r)**1/4, W^D(r) = (1-r)**1/2
!                                   weightd = sqrt(weight)
!                                   weightr = sqrt(weightd)

                                    cigama1 = cigamaF
                                    gamma1  = gammaF

! both j1 and j2 are wall particles
                                    fcVal  = 0.
                                    uVal   = 0.

! one of j1 or j2 is a wall particle and another is a chain(or drop) bead or solvent particle
                                    if((j1 .lt. nWallAtom .and.j2 .gt. nWallAtom) .or.(j1 .gt. nWallAtom .and.j2 .lt. nWallAtom)) then

                                        call BellShape(rij, fcVal)
!                                       fcVal   =  alphawf*weightc
!                                       uVal    = -alphawf*(rij - 0.5*rr)
                                        cigama1 =  cigamaw
                                        gamma1  =  gammaw

                                        k  = min(j1,j2)                 ! wall particle (DP--21/11/11)
                                        i  = max(j1,j2)                 ! solvent or bead particle
                                        dn = dr(1)*wn(k,1) +dr(2)*wn(k,2) +dr(3)*wn(k,3)

                                        if(abs(dn) .le. wLayer .and.rr .lt. sqrDist(i)) then
                                            sqrDist(i) = rr
                                            atomID(i)  = k
                                        endif

! both j1 and j2 are droplet particle
                                    elseif(j1 .gt. nWallAtom .and.j2 .gt. nWallAtom .and.j1 .le. nDpEnd .and.j2 .le. nDpEnd) then
                                            goto 12
!                                        if(int((j1-nWallAtom-1)/nPDP) .eq. int((j2-nWallAtom-1)/nPDP))then
!                                            fcVal= 0.                  !switch attractive force to be zero for bubbles
!!                                            fcVal= alphaf*weightc
!                                            uVal =-alphaf*(rij - 0.5*rr)
!!                                           nbSP(j1) = nbSP(j1) + 1
!!                                           nbSP(j2) = nbSP(j2) + 1
!                                        else
!!                                           if(DpSign(j1) .eq. 1 .or. DpSign(j2) .eq. 1) then
!!                                               fcVal =  alphaf*weightc
!!                                               uVal  = -alphaf*(rij - 0.5*rr)
!!                                               nbDP(j1) = nbDP(j1) + 1
!!                                               nbDP(j2) = nbDP(j2) + 1
!!                                           else
!                                            fcVal= alphaDD*weightc
!                                            uVal =-alphaDD*(rij-0.5*rr)
!!                                               nbDP(j1) = nbDP(j1) + 1
!!                                               nbDP(j2) = nbDP(j2) + 1
!!                                           endif
!                                        endif
!
!                                        cigama1 = cigamaD
!                                        gamma1  = gammaD
                                     

! one of j1 and j2 is a drop particle and another is a solvent particle or chain particle
                                    elseif((j1 .gt. nDpEnd .and.j2 .le. nDpEnd) .or.(j2 .gt. nDpEnd .and.j1 .le. nDpEnd)) then

                                        fcVal  = alphaFD*weightc
                                        uVal   =-alphaFD*(rij - 0.5*rr)
                                        cigama1= cigamaFD
                                        gamma1 = gammaFD
                              

! both j1 and j2 are the beads of chains(or drops)
                                    elseif(j1 .gt. nDpEnd.and.j2 .gt. nDpEnd.and.j1 .le. nChainend .and.j2 .le. nChainend) then

                                        if(int((j1-nWallAtom)/ChainLen).eq. int((j2-nWallAtom)/ChainLen))then

!DP--21/12/2011                             if(abs(j1-j2) .le. 4) goto 10
                                            fcVal= alphapp*weightc
                                            uVal =-alphapp*(rij-0.5*rr)
                                        else
                                            fcVal= alphaf *weightc
                                            uVal =-alphaf *(rij-0.5*rr)
                                        endif
               

! one of j1 and j2 is a bead of chains(or drops) and another is a solvent particle
                                    elseif((j1 .gt. nChainend .and.j2 .le. nChainend) .or.(j2 .gt. nChainend.and.j1 .le. nChainend)) then

                                        fcVal =  alphafp*weightc
                                        uVal  = -alphafp*(rij - 0.5*rr)
                            

! both j1 and j2 are solvent particles
                                    elseif(j1 .gt. nChainend .and.j2 .gt. nChainend) then
                                        fcVal =  alphaf*weightc
                                        uVal  = -alphaf*(rij - 0.5*rr)
                                    endif

10                                  fdVal=-gamma1*weightd*rdvij
                                    frVal= cigama1*weightr*rangls()*sdtinv
12                                  ftot = fcVal + fdVal + frVal
                                    fstl = ftot*rij
! exclude the random force when computing stress (DP--7/12/11)
!                                   fstl = (fcVal + fdVal) * rij
                                    do k = 1, NDIM
!                                       f = ftot*dr(k)
                                        fCV = fcVal*dr(k)
                                        fDP = fdVal*dr(k)
                                        fRD = frVal*dr(k)

                                        raCV(j1,k) = raCV(j1,k) +fCV / mass
                                        raCV(j2,k) = raCV(j2,k) -fCV / mass

                                        raDP(j1,k) = raDP(j1,k) +fDP / mass
                                        raDP(j2,k) = raDP(j2,k) -fDP / mass

                                        raRD(j1,k) = raRD(j1,k) +fRD / mass
                                        raRD(j2,k) = raRD(j2,k) -fRD / mass

                                        fr = fstl*dr(k)*dr(k)
                                        rforce(j1,k) = rforce(j1,k) + fr
                                        rforce(j2,k) = rforce(j2,k) + fr
                                    enddo
                                    fr = fstl*dr(1)*dr(2)
                                    rforce(j1,4) = rforce(j1,4) + fr
                                    rforce(j2,4) = rforce(j2,4) + fr

                                    fr = fstl*dr(1)*dr(3)
                                    rforce(j1,5) = rforce(j1,5) + fr
                                    rforce(j2,5) = rforce(j2,5) + fr

                                    fr = fstl*dr(2)*dr(3)
                                    rforce(j1,6) = rforce(j1,6) + fr
                                    rforce(j2,6) = rforce(j2,6) + fr

                                    uSum = uSum + uVal
                                    virSum = virSum + fcVal*rij
                                endif
                            endif
                            j2 = cellList(j2)
                        enddo
                        j1 = cellList(j1)
                    enddo
15              enddo
            enddo
        enddo
    enddo

!   do n = nStartAtom, nAtom
!       m = ntag(n)
!       do i = 1, m
!           j = tag(n,i)
!           fcVal = alphaB*(mrho(n) + mrho(j)) * weightcd(n,i)          ! repulsive potential in MDPD
!           fstl = fcVal * drij(n,i)
!           do k = 1, NDIM
!               fCR = fcVal * ddr(n,i,k)
!               raCR(n,k) = raCR(n,k) + fCR / mass
!
!               fr = fstl*ddr(n,i,k)*ddr(n,i,k)
!               rforce(n,k) = rforce(n,k) + fr
!           enddo
!           fr = fstl*ddr(n,i,1)*ddr(n,i,2)
!           rforce(n,4) = rforce(n,4) + fr
!
!           fr = fstl*ddr(n,i,1)*ddr(n,i,3)
!           rforce(n,5) = rforce(n,5) + fr
!
!           fr = fstl*ddr(n,i,2)*ddr(n,i,3)
!           rforce(n,6) = rforce(n,6) + fr
!       enddo
!   enddo
! --- for Bell Shape potential case -----
    do n = nWallAtom + 1, nAtom
        m = ntag(n)
        do i = 1, m
            j = tag(n,i)
            if(j .le. nWallAtom) goto 20
            if((n.le. nDpEnd).or.(j.gt. nWallAtom .and.j.le. nDpEnd)) goto 20
            fcVal = alphaB*(mrho(n) + mrho(j)) * weightcd(n,i)      ! repulsive potential in MDPD
            fstl  = fcVal * drij(n,i)
            do k = 1, NDIM
                fCR = fcVal * ddr(n,i,k)
                raCR(n,k) = raCR(n,k) + fCR / mass

                fr = fstl*ddr(n,i,k)*ddr(n,i,k)
                rforce(n,k) = rforce(n,k) + fr
            enddo
            fr = fstl*ddr(n,i,1)*ddr(n,i,2)
            rforce(n,4) = rforce(n,4) + fr

            fr = fstl*ddr(n,i,1)*ddr(n,i,3)
            rforce(n,5) = rforce(n,5) + fr

            fr = fstl*ddr(n,i,2)*ddr(n,i,3)
            rforce(n,6) = rforce(n,6) + fr
20          continue
        enddo
    enddo
!    if(runID .eq. 3 .and. stepCount .gt. stepEquil .and.    &
!       mod((stepCount-stepEquil), 20) .eq. 0) then
!       do n = nWallAtom + 1, nDpEnd
!          pct(n) = real(nbDP(n)) / real(nbDP(n) + nbSP(n))
!          if(pct(n) .gt. cpct) DpSign(n) = 1
!       enddo
!    endif
    deallocate ( cellList, tag, ntag)
    deallocate ( sqrDist ,mrho, ddr, drij,weightcd)
!   print*, 'uSum = ', 0.5*uSUm, '  virSum = ', 0.5*virSum


!DPif(stepCount .gt. stepEquil) then        ! added by DP 3/11/11
    call FENEChainForce
!DPendif                                    ! added by DP 3/11/11

    return
	end

!----------------------------------------------------------

    subroutine ComputeExternalForce

    implicit none
    include 'dpdflow.h'

!   Parameters
!   Locals

    integer	n, i, k
    real*8 CDp(10,NDIM), dCDp

    if((runId .eq. 1 .or. runId .eq. 4) .and.                        &
       gravField .ne. 0.                .and.                        &
       StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nAtom
!       do n = nChainEnd + 1, nAtom
	        ra(n,1) = ra(n,1) + gravField
	    enddo
    endif

    if(runId .eq. 2 .and. gravField .ne. 0. .and.  &
       StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nAtom
!       do n = nChainEnd + 1, nAtom
	        if(r(n,3) .gt. 0.) then
                ra(n,1) = ra(n,1) + gravField
            else
                ra(n,1) = ra(n,1) - gravField
            endif
	    enddo
    endif
    
    if((runId .eq. 5) .and.gravField .ne. 0. .and.StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nAtom
!       do n = nChainEnd + 1, nAtom
	        ra(n,3) = ra(n,3) + gravField
	    enddo
    endif
       
    if(((runId .eq. 6).or.(runID .eq. 64)) .and.gravField .ne. 0. .and.StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nDpEnd
	        ra(n,1) = ra(n,1) + gravField
	    enddo
    endif
    
    if(((runId .eq. 7).or.(runID .eq. 74)) .and.gravField .ne. 0. .and.StepCount .ge. StepEquil+10000) then
        do n = nWallAtom + 1, nDpEnd
	        rv(n,1) = rv(n,1) + gravField
	    enddo
    endif
       
       
    do i = 1, nDp
       do k = 1, NDIM
          do n = nWallAtom + (i-1)*nPDP + 1, nWallAtom + i*nPDP
             if(k .ne. NDIM) CDp(i,k) = CDp(i,k) + (r(n,k) + ncc(n,k)*region(k))
             if(k .eq. NDIM) CDp(i,k) = CDp(i,k) + (r(n,k) + ncc(n,k)*initUcell(k)*gap(k))
          enddo
          CDp(i,k) = CDp(i,k) / nPDP
       enddo
    enddo

!    dCDp = (CDp(1,1) - CDp(2,1)) / RdsDp
!
!    if(dCDp .le. 3.0 .and. dCDp .gt. -5.0) PstAdj = 0
!
    !----drop 1-------
!    do n = nWallAtom + 1, nWallAtom + nPDP
!       ra(n,2) = ra(n,2) - HPstAdj*(CDp(1,2)-CDp0(1,2))
!       if(PstAdj .eq. 1) then
!          ra(n,1) = ra(n,1) - HPstAdj*(CDp(1,1)-CDp0(1,1))
!          ra(n,3) = ra(n,3) - HPstAdj*(CDp(1,3)-CDp0(1,3))
!       endif
!    enddo
       
    !----drop 2-------
!    do n = nWallAtom + nPDP + 1, nDpEnd
!      ra(n,2) = ra(n,2) - HPstAdj*(CDp(2,2)-CDp0(2,2))
!       if(PstAdj .eq. 1) then
!!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(2,1)-CDp0(2,1))
!          ra(n,3) = ra(n,3) - HPstAdj*(CDp(2,3)-CDp0(2,3))
!       endif
!    enddo

!   if(PstAdj .eq. 1) then
!   !----drop 1-------
!      do n = nWallAtom + 1, nWallAtom + nPDP
!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(1,1)-CDp0(1,1))
!         ra(n,2) = ra(n,2) - HPstAdj*(CDp(1,2)-CDp0(1,2))
!         ra(n,3) = ra(n,3) - HPstAdj*(CDp(1,3)-CDp0(1,3))
!      enddo
!   !----drop 2-------
!      do n = nWallAtom + nPDP + 1, nDpEnd
!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(2,1)-CDp0(2,1))
!         ra(n,2) = ra(n,2) - HPstAdj*(CDp(2,2)-CDp0(2,2))
!         ra(n,3) = ra(n,3) - HPstAdj*(CDp(2,3)-CDp0(2,3))
!      enddo
!   endif

	return

	end

!---------------------------------------------------------------------

    subroutine BellShape(rwf, fwf)

    implicit none
    include 'dpdflow.h'

    real*8 rwf, fwf, wca, wcr

    if(rwf .lt. rCut/3.) then
        wca = -18*(rCut/3.)**2/r3Cut + 12*(rCut/3.)/rrCut
    elseif(rwf .lt. 0.5*rCut) then
        wca = -18*rwf**2/r3Cut + 12*rwf/rrCut
    else
        wca =  6*rwf**2/r3Cut - 12*rwf/rrCut + 6/rCut
    endif

    if(rwf .lt. rCut2/3.) then
        wcr = -18*(rCut2/3.)**2/r3Cut2 + 12*(rCut2/3.)/rrCut2
    elseif(rwf .lt. 0.5*rCut2) then
        wcr = -18*rwf**2/r3Cut2 + 12*rwf/rrCut2
    elseif(rwf .le. rCut2) then
        wcr =  6*rwf**2/r3Cut2 - 12*rwf/rrCut2 + 6/rCut2
    endif

    fwf = alphaw*wca + alphawB*wcr

    end
    
!----------------------------------------------------------------------

	subroutine DropIntForce
	
	implicit none
	include 'dpdflow.h'

!	Parameters

!	Locals

	integer	c, i, j1, j2, k, m1, m2, m1X, m1Y, m1Z, m2X, m2Y, m2Z, n,&
            offset
	integer	iofX(14), iofY(14), iofZ(14)
	real*8	fr, fcVal, fdVal, frVal, ftot, weight, rij, weightc,  &
            weightd, weightr, rr, rdvij, uVal, fstl, rangls  
	real*8	ri(NDIM), dr(NDIM), invWid(NDIM), shift(NDIM)

	integer, allocatable :: cellListDp(:),tag(:,:),ntag(:)
    real*8,    allocatable :: sqrDist(:),mrho(:),ddr(:,:,:),drij(:,:),weightcd(:,:)
    real*8 fCV, fDP, fRD

!    integer nbSP(MAXATOM), nbDP(MAXATOM)
!    real*8 pct(MAXATOM)

!---MDPD used variables-----------------------------------------------
    integer j, m
    real*8  fcRval, weightp, fCR
    data    pi / 3.14159265358979 /
	!---------------------------------------------------------------------


	data	iofX / 0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1 /
	data	iofY / 0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1 /
	data	iofZ / 0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1 /

    MaxlistDp=1
    rrCutDp = rCutDp*rCutDp
	do k = 1, NDIM
	   cellsDp(k) = region(k) / rCutDp
	   if(cells(k) .lt. 1) cellsDp(k) = 1
	   invWid(k) = cellsDp(k) / region(k)
       MaxlistDp = MaxlistDp*cellsDp(k)
    enddo
    
    MaxlistDp = MaxlistDp + nDpEnd
   
    allocate ( cellListDp(MaxlistDp), tag(nDpEnd, 150), ntag(nDpEnd))
    allocate ( sqrDist(nDpEnd) ,mrho(nDpEnd), ddr(nDpEnd, 150, NDIM), drij(nDpEnd, 150),weightcd(nDpEnd,150))

	do n = nWallAtom + 1,MaxlistDp
	   cellListDP(n) = 0
	enddo

	do n = nWallAtom + 1, nDpEnd
	   c = (int((r(n,3) + regionH(3))*invWid(3))*cellsDp(2) +           &
     	    int((r(n,2) + regionH(2))*invWid(2)))*cellsDp(1) +          &
     	    int((r(n,1) + regionH(1))*invWid(1)) + nDpEnd + 1
       
	   if(c .lt. 1) then
	      print*, 'in cf, stepCount = ', stepCount
	      print*, ' c = ', c, '  n = ', n
	      print*, 'r(n,*) = ', (r(n,k), k = 1, 3)
	      pause
	   endif

	   if(c .gt.MaxlistDp) then
	      print*, ' c = ', c, '  n = ', n, 'stepCount = ', stepCount
	      print*, 'r(n) =  ', (r(n,k),  k = 1, 3) 
	      pause
	   endif

	   cellListDp(n) = cellListDp(c)
	   cellListDp(c) = n
       
        ntag(n) = 0
        mrho(n) = 0.
       
	enddo
    
!	uSum = 0.
!	virSum = 0.

	do m1Z = 1, cellsDp(3)
	   do m1Y = 1, cellsDp(2)
	      do m1X = 1, cellsDp(1)

	         m1 = ((m1Z - 1)*cellsDp(2) + m1Y - 1)*cellsDp(1) + m1X + nDpEnd

             do offset = 1, 14
		        m2X = m1X + iofX(offset)
		        shift(1) = 0

		        if(m2X .gt. cellsDp(1)) then
		           m2X = 1
		           shift(1) = region(1)
	  	        elseif(m2X .eq. 0) then
		           m2X = cellsDp(1)
		           shift(1) = -region(1)
		        endif
		
		        m2Y = m1Y + iofY(offset)
		        shift(2) = 0

		        if(m2Y .gt. cellsDp(2)) then
		           m2Y = 1
		           shift(2) = region(2)
		        elseif( m2Y .eq. 0) then
	 	           m2Y = cellsDp(2)
		           shift(2) = -region(2)
		        endif

		        m2Z = m1Z + iofZ(offset)
		        shift(3) = 0

                if(runId .eq. 1) then
                   if(m2Z .eq. 0 .or. m2Z .gt. cellsDp(3)) goto 15
                elseif(runId .eq. 0 .or. runId .eq. 2 .or. runId .eq. 3   &
                .or. runId .eq. 4) then
                   if(m2Z .gt. cellsDp(3)) then
		              m2Z = 1
                      shift(3) =  region(3)
 		           elseif(m2Z .eq. 0) then
 		              m2Z = cellsDp(3)
 		              shift(3) = -region(3)
 		           endif
                endif
		
		        m2 = ((m2Z - 1)*cellsDp(2) + m2Y - 1)*cellsDp(1) + m2X + nDpEnd

		        j1 = cellListDp(m1)
		        do while(j1 .gt. 0)
                   j2 = cellListDp(m2)               
		           do while(j2 .gt. 0)

                      if(m1 .ne. m2 .or. j2 .lt. j1) then
  		       
                      do k = 1, NDIM
			             dr(k) = r(j1,k) - r(j2,k) - shift(k)
		              enddo
		       
                      rr = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

		              if(rr .lt. rrCutDp) then
			             rij = sqrt(rr)

!                        weight  = 1. - rij/rCut
!                        weightc = weight

! restrict the conservative cut radius equal to 1. (DP--7/12/11)
                         weightc = 1. - rij/rCutDp
                         weight  = weightc
			                
				         if(weightc .lt. 0.d0) weightc = 0.d0
	  	                 if(weight  .lt. 0.d0) weight  = 0.d0

!	original version of DPD: W^C(r) = W^R(r) = 1 - r, W^D(r) = (1-r)**2
!		                 weightd = weight*weight
!			             weightr = weight

!	modified version of DPD: W^C(r) = 1-r, W^R = (1-r)**1/4, W^D(r) = (1-r)**1/2  
		                 weightd = weight*weight
			             weightr = weight

 			             rdvij = 0.
			             do k = 1, NDIM
			                dr(k) = dr(k)/rij
			                rdvij = rdvij + dr(k)*(rv(j1,k) - rv(j2,k))
                         enddo

! inclusion region for repulsion potential in MDPD
                        if(rij .lt. rCutDp2 .and.alphaDB .gt. 0) then
! weight function related to local density
                            weightp  = 15./(2.*pi*(rCutDp2**3))*(1-rij/rCutDp2)**2
                            ntag(j1) = ntag(j1) + 1
                            ntag(j2) = ntag(j2) + 1
                            tag(j1,ntag(j1)) = j2
                            tag(j2,ntag(j2)) = j1
                            weightcd(j1,ntag(j1)) = 1. - rij/rCutDp2
                            weightcd(j2,ntag(j2)) = weightcd(j1,ntag(j1))
                                     
                            do k = 1, NDIM
                                ddr(j1,ntag(j1),k) =  dr(k)
                                ddr(j2,ntag(j2),k) = -dr(k)
                            enddo                                       
                            drij(j1,ntag(j1)) = rij                                        
                            drij(j2,ntag(j2)) = rij
! local density parameter in MDPD
                            mrho(j1) = mrho(j1) + weightp
                            mrho(j2) = mrho(j2) + weightp
                                      
                        endif


!	both j1 and j2 are droplet particle
			             if(j1 .gt. nWallAtom .and. j2 .gt. nWallAtom .and. &
				            j1 .le. nDpEnd    .and. j2 .le. nDpEnd) then

                            if(int((j1-nWallAtom-1)/nPDP) .eq.  &
                               int((j2-nWallAtom-1)/nPDP)) then
                               fcVal   =  alphaDA*weightc
  				               uVal    = -alphaDA*(rij - 0.5*rr)
                            else
                               fcVal   =  alphaDD*weightc
  				               uVal    = -alphaDD*(rij - 0.5*rr)
                              if(rr .lt. 2.0 ) PstAdj = 0
                            endif
!                         else 
 !                           print*, 'wrong'
!                            pause
                     

			             fdVal = -gammaD*weightd*rdvij
 			             frVal =  cigamaD*weightr*rangls()*sdtinv
		                 ftot  =  fcVal + fdVal + frVal
!		                 fstl  =  ftot*rij
! exclude the random force when computing stress (DP--7/12/11)			  
                         fstl  = (fcVal + fdVal) * rij     
                         do k = 1, NDIM
!			                f = ftot*dr(k)
                            fCV = fcVal*dr(k)
                            fDP = fdVal*dr(k)
                            fRD = frVal*dr(k)

			                raCV(j1,k) = raCV(j1,k) + fCV / mass
			                raCV(j2,k) = raCV(j2,k) - fCV / mass

			                raDP(j1,k) = raDP(j1,k) + fDP / mass
			                raDP(j2,k) = raDP(j2,k) - fDP / mass

			                raRD(j1,k) = raRD(j1,k) + fRD / mass
			                raRD(j2,k) = raRD(j2,k) - fRD / mass

			                fr = fstl*dr(k)*dr(k)
			                rforce(j1,k) = rforce(j1,k) + fr
			                rforce(j2,k) = rforce(j2,k) + fr
		  	             enddo

                         fr = fstl*dr(1)*dr(2)
		 	             rforce(j1,4) = rforce(j1,4) + fr 
		 	             rforce(j2,4) = rforce(j2,4) + fr
			             
                         fr = fstl*dr(1)*dr(3)
		 	             rforce(j1,5) = rforce(j1,5) + fr 
		 	             rforce(j2,5) = rforce(j2,5) + fr
			             
                         fr = fstl*dr(2)*dr(3)
		 	             rforce(j1,6) = rforce(j1,6) + fr
		 	             rforce(j2,6) = rforce(j2,6) + fr
 			  
                         uSum = uSum + uVal
 			             virSum = virSum + fcVal*rij
		              endif
		              endif
		              endif
		              j2 = cellListDp(j2)
 		           enddo
 		           j1 = cellListDp(j1)
		        enddo
 
 15	         enddo

	      enddo
	   enddo
    enddo

! --- for Bell Shape potential case -----
if(alphaDB .gt. 0) then
    do n = nWallAtom + 1, nDpEnd
        m = ntag(n)
        do i = 1, m
            j = tag(n,i)
            if(j .le. nWallAtom) goto 20
            fcVal = alphaDB*(mrho(n) + mrho(j)) * weightcd(n,i)      ! repulsive potential in MDPD
            fstl  = fcVal * drij(n,i)
            do k = 1, NDIM
                fCR = fcVal * ddr(n,i,k)
                raCR(n,k) = raCR(n,k) + fCR / mass

                fr = fstl*ddr(n,i,k)*ddr(n,i,k)
                rforce(n,k) = rforce(n,k) + fr
            enddo
            fr = fstl*ddr(n,i,1)*ddr(n,i,2)
            rforce(n,4) = rforce(n,4) + fr

            fr = fstl*ddr(n,i,1)*ddr(n,i,3)
            rforce(n,5) = rforce(n,5) + fr

            fr = fstl*ddr(n,i,2)*ddr(n,i,3)
            rforce(n,6) = rforce(n,6) + fr
20          continue
        enddo
    enddo
endif
    deallocate ( cellListDp, tag, ntag)
    deallocate ( sqrDist ,mrho, ddr, drij,weightcd)


	return

	end