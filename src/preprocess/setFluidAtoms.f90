!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

    subroutine setPolySolution(isbm, nchsb, ndpsb, chainLen, RdsDp,   &
                          nChainend, nDpEnd, seeds, fgap, wLayer, ncc)

    implicit none

!   Parameters

    integer	MAXCTP
    parameter (MAXCTP = 10)
    integer	chainLen, nChainend, nDpEnd, seeds
    integer	isbm(MAXCTP,8), nchsb(MAXCTP), ndpsb(MAXCTP)
    real*8	fgap, wLayer, RdsDp

!   Locals

    integer	MAXATOM, MAXWATM, MAXBDS, NDIM, MAXDP
    parameter (MAXWATM = 50000, MAXATOM = 400000, NDIM = 3,   &
               MAXBDS = 90000,  MAXDP = 30)
    integer	i, j, k, m, m1, m2, nline, n, nX, nY, nZ, nAtom, nWallAtom, ichain, &
            nsubdom, overlap, ni, nClose, nFlag, nDelete, idom, npart, l1x, &
            l1z, l2x, l2z, l3x, l3z, l4x, l4z, iDp, nPDP
	integer	xuc(4,MAXCTP), zuc(4,MAXCTP), ncontrlp(MAXCTP), nsbstart(MAXCTP), &
	 	    nsbend(MAXCTP), initUcell(NDIM), nSelected(100000), isubreg(MAXCTP,NDIM), &
            ncc(MAXBDS,NDIM), ncc0(MAXBDS,NDIM)
	real*8	r(MAXATOM,NDIM), wn(MAXWATM,NDIM), regionH(NDIM), ctsubm(MAXCTP,NDIM), &
		    region(NDIM), gap(NDIM), rb(MAXBDS, NDIM), rb1(MAXBDS, NDIM), rb2(MAXBDS, NDIM), &
            subreg(MAXCTP,NDIM),tempr(MAXBDS,NDIM),lenDpL
	real*8	subdom(MAXCTP,8), c(NDIM), gap1(NDIM), wmingap, widmin
	real*8  ri(NDIM), h1(4), h2(4), h3(4), rc(MAXDP, 5000, NDIM), s
!   real*8  rDP(MAXATOM,NDIM),rF(MAXATOM,NDIM)      !(Irving-Kirkwood interfacial tension test)
!   integer lDP, lF                                 !(Irving-Kirkwood interfacial tension test)

    logical	giveup
	
	data 	h1 / 0., 0., 0.5, 0.5 /
	data 	h2 / 0., 0.5, 0.5, 0. /
	data 	h3 / 0., 0.5, 0., 0.5 /

	common  /intg / nline, nAtom, nWallAtom, nsubdom
	common  /intgar/ ncontrlp, xuc, zuc, initUcell
	common	/realar/ subdom, gap, wmingap
	common	/coords/ r, wn,  region, regionH,tempr,lenDpL
    common /Pdpr/ nPDP

    n = nWallAtom

    npart  = 0
    rb     = 0.
    widmin = 10000.
    m      = 0

    do idom = 1, nsubdom
        subreg(idom,1) = subdom(idom,3) - subdom(idom,1)
        subreg(idom,2) = region(2)
        subreg(idom,3) = max((subdom(idom,8) - subdom(idom,2)), &
                                (subdom(idom,6) - subdom(idom,4)))
        if(subreg(idom,3) .lt. widmin) widmin = subreg(idom,3)
        ctsubm(idom,1) = subdom(idom,1) + 0.5*subreg(idom,1)            ! coordinates of the subdomain center
        ctsubm(idom,2) = 0.
        ctsubm(idom,3) = min(subdom(idom,2), subdom(idom,4)) &
                                + 0.5*subreg(idom,3)

        nsbstart(idom) = m + 1
        m1 = 0
        do iDp = 1, ndpsb(idom)
            call makedrop(m1, idom, RdsDp, iDp, subreg, seeds, fgap, subdom, rc, rb, ncc0)
        enddo

        do i = 1, m1
            m = m + 1
            do k = 1, NDIM
                rb2(m,k) = rb(i,k)
                ncc(m,k) = ncc0(i,k)
            enddo
        enddo

!        nsbstart(idom) = m + 1     ! the start no. of chain bead coordinate in the array of 'rb2'
        m2 = 0
        do ichain = 1, nchsb(idom)
            call findbeads(m2, idom, ndpsb, rc, RdsDp, &
                            chainLen, subreg, seeds, fgap, subdom, rb, ncc0)
        enddo
        do i = 1, m2
            m = m + 1
            do k = 1, NDIM
                rb2(m,k) = rb(i,k)
                ncc(m,k) = ncc0(i,k)
            enddo
        enddo

        nsbend(idom) = m

        isubreg(idom,1) = isbm(idom,3) - isbm(idom,1)
        isubreg(idom,2) = initUcell(2)
        isubreg(idom,3) = max((isbm(idom,8) - isbm(idom,2)), &
                                (isbm(idom,6) - isbm(idom,4)))
        l1x = isbm(idom,3) - isbm(idom,1)
        l1z = isbm(idom,4) - isbm(idom,2)
        l2x = isbm(idom,5) - isbm(idom,3)
        l2z = isbm(idom,6) - isbm(idom,4)
        l3x = isbm(idom,7) - isbm(idom,5)
        l3z = isbm(idom,8) - isbm(idom,6)
        l4x = isbm(idom,1) - isbm(idom,7)
        l4z = isbm(idom,2) - isbm(idom,8)
        npart = npart + 2*(l1x*l2z - l2x*l1z + l3x*l4z - l4x*l3z) &
                            *isubreg(idom,2)
    enddo

    do i = 1, m1
        n = n + 1
        do k = 1, NDIM
            r(n,k) = rb2(i,k)
        enddo
    enddo
    
    do i = 1, nPDP
        n = n + 1
        m = m + 1
        do k = 1, NDIM
            r(n,k) = tempr(i,k)
            ncc(m,k) = 0
        enddo
    enddo
        
    nDpEnd = n

    do i = 1, m2
        n = n + 1
        do k = 1, NDIM
            r(n,k) = rb2(m1+i,k)
        enddo
    enddo

    nChainend = n

    wLayer = min(0.005*widmin,0.5)
    do k = 1, NDIM
        gap(k) = fgap
    enddo

    ni = 0

!    write(*,*)'**', nsubdom,'**',(isubreg(1,i),i=1,3),'**',(ctsubm(1,i),i=1,3),'**',npart,'**'
    
    do idom = 1, nsubdom
        do nZ = 1, isubreg(idom,3)
            c(3) = (nZ - 0.75)*gap(3) - 0.5*subreg(idom,3) + ctsubm(idom,3)
            do nY = 1, isubreg(idom,2)
                c(2) = (nY - 0.75)*gap(2) - 0.5*subreg(idom,2) + ctsubm(idom,2)
                do nX = 1, isubreg(idom,1)
                    c(1) = (nX - 0.75)*gap(1) - 0.5*subreg(idom,1) + ctsubm(idom,1)
                    do j = 1, 4
                        do k = 1, NDIM
                            if(j .eq. k .or. j .eq. 4) then
                                ri(k) = c(k)
                            else
                                ri(k) = c(k) + 0.5*gap(k)
                            endif
                        enddo
                        
                        
                        if(abs(ri(3)) .lt. (regionH(3)-lenDpL) .and. lenDpL .gt. 0) goto 50
                        
                        do i = 1, ndpsb(idom)
                            s = 0
                            do k = 1, NDIM
                                s = s + (ri(k) - rc(idom, i, k))**2
                            enddo
                            if(s .lt. (RdsDp)**2) goto 50
                        enddo

!                        if(ri(NDIM) .gt. 10) goto 50                !¿ÓµùÂð£¿

                        overlap = 0
                        nClose  = 0

                        do i = nsbstart(idom), nsbend(idom)
                            overlap = 1
                            do k = 1, NDIM
                                if(abs(rb2(i,k) - ri(k)) .gt. 0.) then
!                               if(abs(rb2(i,k) - ri(k)) .gt. 0.5) then
                                    overlap = 0
                                    goto 30
                                endif
                            enddo
30                          if(overlap .eq. 1) goto 50
                            nFlag = 1
                            do k = 1, NDIM
                                if(abs(rb2(i,k) - ri(k)) .gt. 0.50) nFlag = 0
!                               if(abs(rb2(i,k) - ri(k)) .gt. 0.75) nFlag = 0
                            enddo
                            if(nFlag .eq. 1) nClose = 1
                        enddo
                        
                        if(overlap .eq. 0) then
                            call checksubdomBD(ri(1), ri(3), idom, subdom, giveup)
                            if(giveup .eqv. .false.) then
                                n = n + 1
                                do k = 1, NDIM
                                    r(n,k) = ri(k)
                                enddo
                                if(nClose .eq. 1) then
                                    ni = ni + 1
                                    nSelected(ni) = n
                                endif
                            endif
                        endif
50                  enddo
                enddo
            enddo
        enddo
    enddo 

    print*, ' n = ', n, ' nWall = ', nWallAtom, '  npart = ', npart
!    nDelete = n + nsbend(nsubdom) - nWallAtom - npart
    nDelete = n - nWallAtom - npart
    print*, ' nDelete = ', nDelete
    if(nDelete .lt. 0) write(*,'(2x,''Density is lower than the given one'')')
    if(nDelete .gt. ni) write(*,'(2x,''Density adjustment cannot be satisfactory'')')

    do i = 1, min(nDelete, ni)
        j = nSelected(i)
        do k = 1, NDIM
            do m = j, n - 1
                r(m,k) = r(m + 1, k)
            enddo
        enddo
        do k = i + 1, ni
            if(nSelected(k) .gt. j) nSelected(k) = nSelected(k) - 1
        enddo
        n = n - 1
    enddo

!--- Irving-Kirkwood interfacial tension test ---
!    lDP = 0
!    lF  = 0
!    do i = nWallAtom + 1, n
!        if((r(i,NDIM) .lt. 0.25*subreg(1,3)) .and.      &
!           (r(i,NDIM) .ge. -0.25*subreg(1,3))) then
!            lDP = lDP + 1
!            do k = 1, NDIM
!                rDP(lDP,k) = r(i,k)
!            enddo
!        else
!            lF = lF + 1
!            do k = 1, NDIM
!                rF(lF,k) = r(i,k)
!            enddo
!        endif
!    enddo
!
!    nDpEnd = nWallAtom + lDP
!    do i = nWallAtom + 1, nDpEnd
!        do k = 1, NDIM
!            r(i,k) = rDP(i-nWallAtom,k)
!        enddo
!    enddo
!
!    nChainEnd = nDpEnd + nsbend(nsubdom)
!    do i = nDpEnd + 1, nChainEnd
!        do k = 1, NDIM
!            r(i,k) = rb2(i-nDpEnd,k)
!        enddo
!    enddo
!
!    n = n + nsbend(nsubdom)
!    do i = nChainEnd + 1, n
!        do k = 1, NDIM
!            r(i,k) = rF(i-nChainEnd, k)
!        enddo
!    enddo
!--------------- END -----------------------

    nAtom = n
    write(*,*) 'natom=',natom

    return
    end