!-------------------------------------------------------------------
!	chainconf.f90
!	generate the configuration of wall, chain and free dpd parti-
!	cels for simulating polymer solutions. Configuration data can
!	be input from the terminal and wtitten into "fort.10" at the 
!	first time and then read from it next time.
!-------------------------------------------------------------------

    implicit none

    integer	MAXCTP, MAXWATM, MAXATOM, NDIM, MAXBDS, MAXDP
    parameter (MAXCTP = 10, MAXWATM = 50000, MAXATOM = 400000, NDIM = 3, &
               MAXBDS = 90000, MAXDP = 30)
    integer	i, j, k, n, nline, ly, nAtom, nWallAtom, nsubdom, &
            nChain, ChainLen, nChainend, nDpEnd, seeds, nBeads1, nBeads2, nDp
    integer	xuc(4,MAXCTP), zuc(4,MAXCTP), ncontrlp(MAXCTP),   &
            initUcell(NDIM), key2, DpFzn(MAXDP),ntemp(2),nPDP,LenDp
    integer	isbm(MAXCTP,8), nchsb(MAXCTP), ndpsb(MAXCTP), ncc(MAXBDS,NDIM)
    real*8	density,dratio, ratio, fgap, wgap, wmingap, h3, wLayer, RdsDp,lenDpL,pdratio
    real*8	r(MAXATOM,NDIM), wn(MAXWATM,NDIM), regionH(NDIM), &
            region(NDIM), gap(NDIM),tempr(MAXBDS,NDIM),tChainLen
    real*8	subdom(MAXCTP,8), Cx(MAXDP), Cy(MAXDP), Cz(MAXDP)
    character key1, lattice*3
    character*30 fdata

    common  /intg / nline, nAtom, nWallAtom, nsubdom
    common  /intgar/ ncontrlp, xuc, zuc, initUcell
    common  /realar/ subdom, gap, wmingap
    common  /coords/ r, wn, region, regionH,tempr,lenDpL
    common  /dp/ Cx, Cy, Cz, density,dratio,pdratio
    common /dpr/ nDp
    common /chain/ nChain
    common /nup/nChainend, nDpEnd
    common /Pdpr/ nPDP

    !write(*,*) 'input data from the terminal (y/n) ?'
    !read(*,'(a)') key1
    key1='n'
!    write(*,*) 'input the name of data file name to read or write on:'
!    read(*,'(a)') fdata
     open(10, file = './data/lint.dat')
!    open(10, file = './data/'//fdata)

    if(key1 .eq. 'y' .or. key1 .eq. 'Y') then
        write(* ,*) 'input the density of fluid molecules'
        write(* ,*) 'and the density ratio of wall/fluid:'
        read (* ,*) density, ratio
        write(10,*) density, ratio
        fgap = (4.0 / density)**(1./3.)
        
        write(* ,*) 'input data of flow geometry:'
        write(* ,*) 'the thickness (in y-direction) in unit cells:'
        read (* ,*) ly
        write(10,*) ly

        write(* ,*) 'number of wall lines  ='
        read (* ,*) nline
        write(10,*) nline

        write(* ,*) 'number of control points for each line:'
        do i = 1, nline
            write(* ,'(''number of control points of line : '', i3)') i
            read (* ,*) ncontrlp(i)
            write(10,*) ncontrlp(i)
        enddo

        write(*,*) 'enter coordinates of controle points (in unit cells)'
        do i = 1, nline
            do j = 1, ncontrlp(i)
                write(* ,'(''line :'', i3, ''  point :'', i3)') i, j
                read (* ,*) xuc(i,j), zuc(i,j)
                write(10,*) xuc(i,j), zuc(i,j)
            enddo
        enddo

        seeds = 27

        write(* ,*) 'enter length of chians and number of chains :'
        read (* ,*) ChainLen, nChain
        write(10,*) ChainLen, nChain
        write(*, *) 'enter x-direction length and density ratio of Droplets:'
        read(*,*) LenDp, dratio
        write(10,*) LenDp, dratio
        write(*, *) 'including droplet or not? (1/0)'
        read (*, *) key2
        write(10,*) key2
        if(key2 .eq. 1) then
            write(* ,*) 'enter radius of droplets ,number and density ratio of droplets :'
            read (* ,*) RdsDp, nDp,dratio
            write(10,*) RdsDp, nDp,dratio
            write(* ,*) 'enter coordinates of droplets centers: (x,y,z)'
            do i = 1, nDp
                read (* ,*) Cx(i), Cy(i), Cz(i)
                write(10,*) Cx(i), Cy(i), Cz(i)
            enddo
        else
            nDp   = 0
            RdsDp = 0.
        endif
    else
        write(* ,*) 'read in data from specified Data file ...'
        read (10,*) density, ratio

        read (10,*) ly
        read (10,*) nline
        do i = 1, nline
            read (10,*) ncontrlp(i)
        enddo
        do i = 1, nline
            do j = 1, ncontrlp(i)
                read (10,*) xuc(i,j), zuc(i,j)
            enddo
        enddo
        read(10,*) ChainLen, nChain
        read(10,*) LenDp, pdratio
        read(10,*) key2
        if(key2 .eq. 1) then
            read(10,*) RdsDp, nDp,dratio
            do i = 1, nDp
!               read(10,*) DpFzn(i)
                read(10,*) Cx(i), Cy(i), Cz(i)
            enddo
        else
            nDp   = 0
            RdsDp = 0.
        endif
    endif
    
    close(10)
    
    fgap = (4.0 / density)**(1./3.)
    
    if (LenDp .gt.0) then

        do i = 1, nline
            do j = 1, ncontrlp(i)
                zuc(i,j)=zuc(i,j)+LenDp*zuc(i,j)/abs(zuc(i,j))
            enddo
        enddo
        tChainLen = ChainLen
        ntemp(1)= nChain
        ntemp(2)= nDp
        
        ChainLen = 0
        nChain = 0
        nDp = 0
        if(pdratio .gt.0)then
            fgap = (4.0*pdratio / density)**(1./3.)
        else 
            fgap = (4.0 / density)**(1./3.)    !build a vacant space
        endif
        
        wgap = (4.0 / density/ratio)**(1./3.)
        gap(1) = fgap
        gap(2) = fgap
        gap(3) = fgap
        initUcell(2) = ly
        region(2) = ly*fgap
        regionH(2) = 0.5*region(2)
        
        i=0
        nPDP=0
        lenDpL=0.
        
        if(pdratio .gt.0)then
            
            lattice = 'FCC'
        
!            call setWallAtoms(fgap, wgap, nWallAtom, h3)
    
	        call setSubDomain(fgap, nChain, nDp, isbm, nchsb, ndpsb)
    
    
	        call setPolySolution(isbm, nchsb, ndpsb, chainLen, RdsDp, &
			                     nChainend, nDpEnd, seeds, fgap, wLayer, ncc)
            
            lenDpL = lenDp*fgap
            do n=nWallAtom+1,nAtom      
                if(abs(r(n,3)) .lt. (regionH(3)-lenDpL))then
                    i=i+1
                    do k=1,NDIM                    
                        tempr(i,k) = r(n,k) 
                    enddo
                endif
            enddo
                
            nPDP = i
    !refresh parameter    
            fgap = (4.0 / density)**(1./3.)
            ly= aint(region(2)/fgap)        
            do i = 1, nline
                do j = 1, ncontrlp(i)
                    xuc(i,j) = aint(regionH(1)/fgap)*xuc(i,j)/abs(xuc(i,j))
                    zuc(i,j) = aint(regionH(3)/fgap)*zuc(i,j)/abs(zuc(i,j))
                enddo
            enddo
        endif        
        
        lenDpL = lenDp*fgap
        ChainLen = tChainLen
        nChain = ntemp(1)
        nDp = ntemp(2)
       
    endif
    
   print*,nPDP,lenDpL
    
    wgap = (4.0 / density/ratio)**(1./3.)
    gap(1) = fgap
    gap(2) = fgap
    gap(3) = fgap
    initUcell(2) = ly
    region(2) = ly*fgap
    regionH(2) = 0.5*region(2)

    write(*,*) 'start to generate the wall molecules ...'
    nWallAtom = 0
    lattice = 'FCC'

    
    call setWallAtoms(fgap, wgap, nWallAtom, h3)

    
	call setSubDomain(fgap, nChain, nDp, isbm, nchsb, ndpsb)
    
    
	call setPolySolution(isbm, nchsb, ndpsb, chainLen, RdsDp, &
			             nChainend, nDpEnd, seeds, fgap, wLayer, ncc)


    
	region(3)  = region(3) + 2.*h3 + 0.25*wgap
	regionH(3) = 0.5*region(3)
    
    
    

   	write(*,*) 'write coordinates of wall molecules into "fort.20" ...'

	nBeads1 = nDpEnd    - nWallAtom 
    nBeads2 = nChainend - nDpEnd
    
    open(20, file = './data/fort.20')
	write(20,'( ''fgap = '', f8.4, 2x, '' wgap = '', f8.4)') fgap, wgap
	write(20,'(/''No. of WallParticles     = '', i7)') nWallAtom
    write(20,'( ''No. of Droplet Particles = '', i7)') nBeads1
	write(20,'( ''No. of Chain Particles   = '', i7)') nBeads2
	write(20,'( ''No. of Total Particles   = '', i7)') nAtom
	write(20,'(/''Volume Concentration(Droplet)  = '', f8.5)') float(nBeads1)/float(nAtom-nWallAtom)
	write(20,'(/''Volume Concentration(Chain)    = '', f8.5)') float(nBeads2)/float(nAtom-nWallAtom)
    if(key2 .eq. 1) then
       write(20,'(/''Droplet info:'')')
	   write(20,'( ''nDrop  = '', i6,  4x, ''Radius   = '', f8.4)') nDp, RdsDp
    endif

    if(nChain .gt. 0) then
       write(20,'(/''Chain info:'')')
	   write(20,'( ''nChain = '', i6,  4x, ''ChainLen = '', i4)') nChain, ChainLen
    endif

 	write(20,'(/''Coordinates of nodes of subdomains'')')
 	do k = 1, nsubdom
 	   write(20,'(''no. of subdomain :'',i3)') k
 	   write(20,'(4f10.4)') (subdom(k,i), i = 1, 7, 2)
 	   write(20,'(4f10.4)') (subdom(k,i), i = 2, 8, 2)
 	enddo
 	write(20,'(''Thickness of Boundary Layer :'', f8.5)') wLayer

	write(20,'(//5x,'' region size: '')')
	write(20,'(2x, 3f9.4)') region(1), region(2), region(3)
	write(20,'(2x, 3f9.4)') regionH(1), regionH(2), regionH(3)

 	write(20,'(//'' Coordinates of Wall Particles'')')
 	write(20,'(2x, ''n'', 6x,''x'', 8x, ''y'', 8x, ''z'', 8x, &
		  ''nx'', 7x, ''ny'', 7x, ''nz'')')
	do n = 1, nWallAtom
	   write(20,'(i6, 6f9.4)') n, (r(n,i), i = 1,3), (wn(n,i),i = 1,3)
	enddo

    if(key2 .eq. 1) then
       write(20,'(//'' Coordinates of Droplet Particles''/)')
       do n = nWallAtom + 1, nDpEnd
	      write(20,'(i6, 3f9.4, 3i4)') n, (r(n,i), i = 1,3), (ncc(n-nWallAtom,i), i = 1,3)
	   enddo
    endif

    if(nChain .gt. 0) then
       write(20,'(//'' Coordinates of Chain Particles''/)')
       do n = nDpEnd + 1, nChainend
	      write(20,'(i6, 3f9.4)') n, (r(n,i), i = 1,3)
	   enddo
    endif

 	write(20,'(//'' Coordinates of Simple Particles''/)')
	do n = nChainend + 1, nAtom
	   write(20,'(i6, 3f9.4)') n, (r(n,i), i = 1,3)
	enddo

	write(*,'(2x,''write binary data to "initdpd.cnf".....'')')
	open(16, file = './data/initdpd.cnf')

	write(16,*) nWallAtom, nDpEnd, nChainend, nAtom, nDp, RdsDp,    &
              DpFzn(1:nDp), nChain, ChainLen, initUcell, region,  &
              regionH, gap, wmingap, wLayer

!	write(16,*) '*********'

    
	do k = 1, NDIM
	   write(16,*) (wn(n,k), n = 1, nWallAtom)
	   write(16,*) (r(n,k), n = 1, nAtom)
    enddo

    if(nDp .gt. 0) then
        do k = 1, NDIM
            write(16,*) (ncc(n,k), n = 1, nBeads1)
        enddo
    endif
    
    if(nChain .gt. 0) then
        do k = 1, NDIM
            write(16,*) (ncc(n,k), n = nBeads1+1, nBeads1+nBeads2)
        enddo
    endif

	close (16)
	close (11)
    
    write(*,*) 'nDp=',nDp,'nChain=',nChain
    !pause
    call ouptParticleLocation
    
    
	stop

	end
	
!--------------------------------------------------------------------------

	subroutine setSubDomain(fgap, nChain, nDp, isbm, nchsb, ndpsb)

	implicit none

!	Parameters

	integer	nChain, MAXCTP, nDp
	parameter (MAXCTP = 10)
	real*8	fgap
	integer	isbm(MAXCTP,8), nchsb(MAXCTP), ndpsb(MAXCTP)

!	Locals

	integer	MAXWATM, MAXATOM, NDIM
	parameter (MAXWATM = 50000, MAXATOM = 400000, NDIM = 3)
	integer	i, j, k, m1, m2, nline, nAtom, nWallAtom, nsubdom, ntot
	integer	ix1, ix2, ix3, ix4, iz1, iz2, iz3, iz4, is
	integer	xuc(4,MAXCTP), zuc(4,MAXCTP), ncontrlp(MAXCTP),&
	 	    initUcell(NDIM)
	real*8	x1, x2, z1, z2, wmingap 
	real*8	r(MAXATOM,NDIM), wn(MAXWATM,NDIM), regionH(NDIM), &
		    region(NDIM), gap(NDIM)
	real*8	subdom(MAXCTP,8)
!DP	character key
	common  /intg / nline, nAtom, nWallAtom, nsubdom
	common  /intgar/ ncontrlp, xuc, zuc, initUcell
	common	/realar/ subdom, gap, wmingap
	common	/coords/ r, wn, region, regionH

	k = 0
	ix2 = xuc(1,1)
	iz2 = zuc(1,1)
 	ix3 = xuc(2,1)
	iz3 = zuc(2,1)
	do i = 1, ncontrlp(1) - 1
	   ix1 = ix2
	   iz1 = iz2
	   ix4 = ix3
	   iz4 = iz3
	   ix2 = xuc(1,i+1)
	   iz2 = zuc(1,i+1)
	   ix3 = xuc(2,i+1)
	   iz3 = zuc(2,i+1)
	   is = (ix3 - ix1)*(iz4 - iz2) + (iz3 - iz1)*(ix2 - ix4)
	   if(iabs(is) .gt. 0) then
	      k = k + 1
	      isbm(k,1) = ix1 
	      isbm(k,2) = iz1 
	      isbm(k,3) = ix2 
	      isbm(k,4) = iz2 
	      isbm(k,5) = ix3 
	      isbm(k,6) = iz3 
	      isbm(k,7) = ix4 
	      isbm(k,8) = iz4
	      nchsb(k) = is
	   endif
	enddo
	
	nsubdom = k
	
    ntot = 0	
	do k = 1, nsubdom
	   ntot = ntot + nchsb(k)
	   do j = 1, 7, 2
	      subdom(k,j) = isbm(k,j)*gap(1)
	      subdom(k,j+1) = isbm(k,j+1)*gap(3)
	   enddo
	enddo
	
!	create size of region

	ix1 =  10000
	iz1 =  10000
	ix2 = -10000
	iz2 = -10000

	m1 = 0
    m2 = 0
	do k = 1, nsubdom
!DP	   if(nChain .eq. 0) then
!DP	      nchsb(k) = 0
!DP	   else
!DP	      nchsb(k) = (nChain*nchsb(k))/ntot + 1
!DP	   endif
	   if(nChain .eq. 0) then                        !---DP
	      nchsb(k) = 0                                 
	   elseif(nsubdom .eq. 1) then
	      nchsb(k) = nChain
	   else
	      nchsb(k) = (nChain*nchsb(k))/ntot + 1
	   endif                                         !---DP

       if(nDp .eq. 0) then                           !---DP
	      ndpsb(k) = 0                                 
	   elseif(nsubdom .eq. 1) then
	      ndpsb(k) = nDp
	   else
	      ndpsb(k) = (nDp*nchsb(k))/ntot + 1
	   endif                                         !---DP

	   do j = 1, 7, 2
	      if(isbm(k,j) .lt. ix1) ix1 = isbm(k,j)
	      if(isbm(k,j) .gt. ix2) ix2 = isbm(k,j)
	      if(isbm(k,j+1) .lt. iz1) iz1 = isbm(k,j+1)
	      if(isbm(k,j+1) .gt. iz2) iz2 = isbm(k,j+1)
	   enddo
	   m1 = m1 + nchsb(k)
       m2 = m2 + ndpsb(k)
	enddo
	nChain = m1
    nDp    = m2
	
	initUcell(1) = (ix2 - ix1)
	initUcell(3) = (iz2 - iz1)
	region(1) = initUcell(1)*gap(1)
	region(3) = initUcell(3)*gap(3)
	regionH(1) = 0.5*region(1)
	regionH(3) = 0.5*region(3)

	return

	end

!-------------------------------------------------------------------------

	subroutine checksubdomBD(x, z, idom, subdom, giveup)

	implicit none
	
!	parameters

	integer	MAXCTP, idom
	parameter (MAXCTP = 10)
	real*8	x, z, subdom(MAXCTP,8)
	logical	giveup

!	Locals

	integer	i, i1, i2
	real*8	dn, t1, t2, sl, delta

	giveup = .false.

	do i = 1, 4
	
!	check if the particles located inside of the subdom

	   i1 = i
	   i2 = i + 1
	   if(i1 .eq. 4) i2 = 1
	   delta = (x - subdom(idom,2*i1-1))*(z - subdom(idom,2*i2)) - &
		       (x - subdom(idom,2*i2-1))*(z - subdom(idom,2*i1)) 	
	   if(delta .lt. 0.) then
	      giveup = .true.
	      return
	   endif

!	check the distance from the particle to the subdom boundary

	   t1 = subdom(idom,2*i2-1) - subdom(idom,2*i1-1)
	   t2 = subdom(idom,2*i2) - subdom(idom,2*i1)
	   sl = sqrt(t1*t1 + t2*t2)
	   dn = (-(x - subdom(idom,2*i1-1))*t2 +  &
		    (z - subdom(idom,2*i1))*t1)/sl
	   if(dn .lt. 0.2) then
	      giveup = .true.
	      return
	   endif

	enddo

	return

	end