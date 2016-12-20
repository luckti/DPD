!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

    subroutine setWallAtoms(fgap, wgap, n, h3)

    implicit none

!   Parameters

    integer	n
    real*8	wgap, fgap, h3

!   Locals

    integer	MAXCTP, MAXWATM, MAXATOM, NDIM
    parameter (MAXCTP = 10, MAXWATM = 50000, MAXATOM = 400000, NDIM = 3)
    integer	i, j, k, l, n1, n2, ny, nctp, nline, nAtom, nWallAtom, nsubdom
    integer	xuc(4,MAXCTP), zuc(4,MAXCTP), ncontrlp(MAXCTP), initUcell(NDIM)
    real*8	r(MAXATOM,NDIM), wn(MAXWATM,NDIM), regionH(NDIM),  &
            region(NDIM), gap(NDIM)
    real*8	rw(MAXWATM,NDIM), subdom(MAXCTP,8), c(NDIM)
    real*8  length, x1, x2, z1, z2, tx1, tz1, tx2, tz2
    real*8	d, a1, a2, h1, h2, delta, regiony, gap1, gap2, wmingap, det1, det2
    real*8	wx1(4,MAXCTP), wx2(4,MAXCTP), wz1(4,MAXCTP), wz2(4,MAXCTP), &
            wx3(4,MAXCTP), wz3(4,MAXCTP), norml(2,MAXCTP,2)

    common  /intg / nline, nAtom, nWallAtom, nsubdom
    common  /intgar/ ncontrlp, xuc, zuc, initUcell
    common  /realar/ subdom, gap, wmingap
    common  /coords/ r, wn, region, regionH

    d  = -1.0
!   h1 =  0.25*fgap
    h1 =  0.
    h2 =  h1 + 0.5*wgap
    h3 =  h2 + 0.5*wgap
    wmingap = 1000.

    do i = 1, nline
        d = -d
        nctp = ncontrlp(i)
        x1 = xuc(i,1)*fgap
        z1 = zuc(i,1)*fgap
        x2 = xuc(i,2)*fgap
        z2 = zuc(i,2)*fgap
        length = sqrt((x2 - x1)**2 + (z2 - z1)**2)
        tx2 = d*(x2 - x1)/length
        tz2 = d*(z2 - z1)/length
        det2 = d
        wx1(i,1) = x1 + h1*tz2
        wz1(i,1) = z1 - h1*tx2
        wx2(i,1) = x1 + h2*tz2
        wz2(i,1) = z1 - h2*tx2
        wx3(i,1) = x1 + h3*tz2
        wz3(i,1) = z1 - h3*tx2
        do j = 2, nctp - 1
            x1 = x2
            z1 = z2
            tx1 = tx2
            tz1 = tz2
            det1 = det2
            x2 = xuc(i,j+1)*fgap
            z2 = zuc(i,j+1)*fgap
            length = sqrt((x2 - x1)**2 + (z2 - z1)**2)
            tx2 = d*(x2 - x1)/length
            tz2 = d*(z2 - z1)/length
            det2 = tx1*tz2 - tz1*tx2
            if(det2*d .gt. 0.) then
                norml(i,j-1,1) = -0.707106781*(tz1 + tz2)
                norml(i,j-1,2) =  0.707106781*(tx1 + tx2)
            elseif(det1*d .gt. 0) then
                norml(i,j-1,1) = -tz2
                norml(i,j-1,2) =  tx2
            else
                norml(i,j-1,1) = -tz1
                norml(i,j-1,2) =  tx1
            endif

            delta = det2

            a1 = tz1*x1 - tx1*z1
            a2 = tz2*x2 - tx2*z2

            wx1(i,j) = (-(a1 + h1)*tx2 + (a2 + h1)*tx1)/delta
            wz1(i,j) = ( (a2 + h1)*tz1 - (a1 + h1)*tz2)/delta
            wx2(i,j) = (-(a1 + h2)*tx2 + (a2 + h2)*tx1)/delta
            wz2(i,j) = ( (a2 + h2)*tz1 - (a1 + h2)*tz2)/delta
            wx3(i,j) = (-(a1 + h3)*tx2 + (a2 + h3)*tx1)/delta
            wz3(i,j) = ( (a2 + h3)*tz1 - (a1 + h3)*tz2)/delta
        enddo
        wx1(i,nctp) = x2 + h1*tz2
        wz1(i,nctp) = z2 - h1*tx2
        wx2(i,nctp) = x2 + h2*tz2
        wz2(i,nctp) = z2 - h2*tx2
        wx3(i,nctp) = x2 + h3*tz2
        wz3(i,nctp) = z2 - h3*tx2
    enddo

!   do i = 1, nline
!       do j = 1, ncontrlp(i)
!           write(20,'(i3,4f9.4)') j, wx1(i,j), wz1(i,j),wx2(i,j),wz2(i,j)
!       enddo
!   enddo

!   determine the gap and divi

    ny = int(region(2)/wgap)
    if(ny*wgap .lt. region(2) - 0.5*wgap) ny = ny + 1
    gap2 = region(2) / ny
    d = -1.0

    do i = 1, nline
        d = -d
        nctp = ncontrlp(i)
        do j = 1, nctp - 1
            x1 = wx1(i,j)
            z1 = wz1(i,j)
            x2 = wx1(i,j+1)
            z2 = wz1(i,j+1)
            length = sqrt((x2 - x1)**2 + (z2 - z1)**2)
            n1 = int(length/wgap)
!           if(n1*wgap .lt. length) n1 = n1 + 1
            if(j .eq. 1) then
                gap1 = 4.*length/(4.*n1 + 1.0)
                if(nctp .eq. 2) gap1 = length/n1
            elseif(j .eq. ncontrlp(i) - 1) then
                gap1 = 4.*length/(4.*n1 + 3.0)
            else
                gap1 = length/n1
            endif

            if(gap1 .lt. wmingap) wmingap = gap1

            tx1 = (x2 - x1)/length
            tz1 = (z2 - z1)/length
            n2  = n1 + 1

            if(nctp .eq. 2) n2 = n1

            do k = 1, n2
                if(j .eq. 1) then
                    c(1) = x1 + (k - 0.75)*gap1*tx1
                    c(3) = z1 + (k - 0.75)*gap1*tz1
                else
                    c(1) = x1 + (k - 1.0)*gap1*tx1
                    c(3) = z1 + (k - 1.0)*gap1*tz1
                endif

                do l = 1, ny
                    c(2) = (l - 0.75)*gap2 - regionH(2)
                    if(j .eq. 1 .or. k .ne. 1) then
                        n = n + 1
                        rw(n,1) = c(1)
                        rw(n,2) = c(2)
                        rw(n,3) = c(3)

                        if(abs(c(1)-x2) .lt. 1.e-4 .and. abs(c(3)-z2) .lt. 1.e-4) then
                            wn(n,1) = norml(i,j,1)
                            wn(n,2) = 0.
                            wn(n,3) = norml(i,j,2)
                        else
                            wn(n,1) = -d*tz1
                            wn(n,2) = 0.
                            wn(n,3) = d*tx1
                        endif
                    endif

                    if(k .le. n1 .or. j .eq. nctp - 1) then
                        n = n + 1
                        rw(n,1) = c(1) + 0.5*gap1*tx1
                        rw(n,2) = c(2) + 0.5*gap2
                        rw(n,3) = c(3) + 0.5*gap1*tz1
                        wn(n,1) = -d*tz1
                        wn(n,2) = 0.
                        wn(n,3) = d*tx1
                    endif
                enddo
            enddo
        enddo
    enddo

    d = -1.0
    do i = 1, nline
        d = -d
        nctp = ncontrlp(i)
        do j = 1, nctp - 1
            x1 = wx2(i,j)
            z1 = wz2(i,j)
            x2 = wx2(i,j+1)
            z2 = wz2(i,j+1)
            length = sqrt((x2 - x1)**2 + (z2 - z1)**2)
            n1 = int(length/wgap)
!           if(n1*wgap .lt. length) n1 = n1 + 1
            if(j .eq. 1) then
                gap1 = 4.*length/(4.*n1 + 3.0)
                if(nctp .eq. 2) gap1 = length/n1
            elseif(j .eq. nctp - 1) then
                gap1 = 4.*length/(4.*n1 + 1.0)
            else
                gap1 = length/n1
            endif
            tx1 = (x2 - x1)/length
            tz1 = (z2 - z1)/length
            n2  = n1 + 1
            if(nctp .eq. 2) n2 = n1
            do k = 1, n2
                if(j .eq. 1) then
                    c(1) = x1 + (k - 0.25)*gap1*tx1
                    c(3) = z1 + (k - 0.25)*gap1*tz1
                else
                    c(1) = x1 + (k - 1.0)*gap1*tx1
                    c(3) = z1 + (k - 1.0)*gap1*tz1
                endif
                do l = 1, ny
                    c(2) = (l - 0.75)*gap2 - regionH(2)
                    if(j .eq. 1 .or. k .ne. 1) then
                        n = n + 1
                        rw(n,1) = c(1)
                        rw(n,2) = c(2)
                        rw(n,3) = c(3)
                        if(abs(c(1)-x2) .lt. 1.e-4 .and. abs(c(3)-z2) .lt. 1.e-4) then
                            wn(n,1) = norml(i,j,1)
                            wn(n,2) = 0.
                            wn(n,3) = norml(i,j,2)
                        else
                            wn(n,1) = -d*tz1
                            wn(n,2) = 0.
                            wn(n,3) = d*tx1
                        endif
                    endif
                    if(k .le. n1) then
                        if(nctp .eq. 2 .and. k .eq. n1) goto 15
                        n = n + 1
                        rw(n,1) = c(1) + 0.5*gap1*tx1
                        rw(n,2) = c(2) + 0.5*gap2
                        rw(n,3) = c(3) + 0.5*gap1*tz1
                        wn(n,1) = -d*tz1
                        wn(n,2) = 0.
                        wn(n,3) = d*tx1
15                  endif
                    if(j .eq. 1 .and. k .eq. 1) then
                        n = n + 1
                        rw(n,1) = c(1) - 0.5*gap1*tx1
                        rw(n,2) = c(2) + 0.5*gap2
                        rw(n,3) = c(3) - 0.5*gap1*tz1
                        wn(n,1) = -d*tz1
                        wn(n,2) = 0.
                        wn(n,3) = d*tx1
                    endif
                enddo
            enddo
        enddo
    enddo

    d = -1.0
    do i = 1, nline
        d = -d
        nctp = ncontrlp(i)
        do j = 1, nctp - 1
            x1 = wx3(i,j)
            z1 = wz3(i,j)
            x2 = wx3(i,j+1)
            z2 = wz3(i,j+1)
            length = sqrt((x2 - x1)**2 + (z2 - z1)**2)
            n1 = int(length/wgap)
!           if(n1*wgap .lt. length) n1 = n1 + 1
            if(j .eq. 1) then
                gap1 = 4.*length/(4.*n1 + 3.0)
                if(nctp .eq. 2) gap1 = length/n1
            elseif(j .eq. nctp - 1) then
                gap1 = 4.*length/(4.*n1 + 1.0)
            else
                gap1 = length/n1
            endif
            tx1 = (x2 - x1)/length
            tz1 = (z2 - z1)/length
            n2 = n1 + 1
            if(nctp .eq. 2) n2 = n1
            do k = 1, n2
                if(j .eq. 1) then
                    c(1) = x1 + (k - 0.25)*gap1*tx1
                    c(3) = z1 + (k - 0.25)*gap1*tz1
                else
                    c(1) = x1 + (k - 1.0)*gap1*tx1
                    c(3) = z1 + (k - 1.0)*gap1*tz1
                endif
                do l = 1, ny
                    c(2) = (l - 0.75)*gap2 - regionH(2)
                    if(j .eq. 1 .or. k .ne. 1) then
                        n = n + 1
                        rw(n,1) = c(1)
                        rw(n,2) = c(2)
                        rw(n,3) = c(3)
                        if(abs(c(1)-x2) .lt. 1.e-4 .and. abs(c(3)-z2) .lt. 1.e-4) then
                            wn(n,1) = norml(i,j,1)
                            wn(n,2) = 0.
                            wn(n,3) = norml(i,j,2)
                        else
                            wn(n,1) = -d*tz1
                            wn(n,2) = 0.
                            wn(n,3) = d*tx1
                        endif
                    endif
                    if(k .le. n1) then
                        if(nctp .eq. 2 .and. k .eq. n1) goto 25
                        n = n + 1
                        rw(n,1) = c(1) + 0.5*gap1*tx1
                        rw(n,2) = c(2) + 0.5*gap2
                        rw(n,3) = c(3) + 0.5*gap1*tz1
                        wn(n,1) = -d*tz1
                        wn(n,2) = 0.
                        wn(n,3) = d*tx1
25                  endif
                    if(j .eq. 1 .and. k .eq. 1) then
                        n = n + 1
                        rw(n,1) = c(1) - 0.5*gap1*tx1
                        rw(n,2) = c(2) + 0.5*gap2
                        rw(n,3) = c(3) - 0.5*gap1*tz1
                        wn(n,1) = -d*tz1
	                    wn(n,2) = 0.
                        wn(n,3) = d*tx1
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 1, NDIM
        do j = 1, n
            r(j,k) = rw(j,k)
        enddo
    enddo

	wmingap = 0.707106*wmingap

	return

	end

!----------------------------------------------------------------------------