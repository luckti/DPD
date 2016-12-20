!! Added by linyuqing

    subroutine ouptParticleLocation
    
    implicit none
        integer	MAXCTP, MAXWATM, MAXATOM, NDIM, MAXBDS, MAXDP
    parameter (MAXCTP = 10, MAXWATM = 50000, MAXATOM = 400000, NDIM = 3, &
               MAXBDS = 90000, MAXDP = 30)
    integer	i, j, k, n, nline, ly, nAtom, nWallAtom, nsubdom, &
            nChain, ChainLen, nChainend, nDpEnd, seeds, nBeads1, nBeads2, nDp
    integer	xuc(4,MAXCTP), zuc(4,MAXCTP), ncontrlp(MAXCTP),   &
            initUcell(NDIM), key2, DpFzn(MAXDP)
    integer	isbm(MAXCTP,8), nchsb(MAXCTP), ndpsb(MAXCTP), ncc(MAXBDS,NDIM)
    real*8	density,dratio, ratio, fgap, wgap, wmingap, h3, wLayer, RdsDp
    real*8	r(MAXATOM,NDIM), wn(MAXWATM,NDIM), regionH(NDIM), &
            region(NDIM), gap(NDIM)
    real*8	subdom(MAXCTP,8), Cx(MAXDP), Cy(MAXDP), Cz(MAXDP)
    
    common  /intg / nline, nAtom, nWallAtom, nsubdom
    common  /coords/ r, wn, region, regionH
    common /dpr/ nDp
    common /chain/ nChain
    common /nup/nChainend, nDpEnd
    !local
    integer step
    character*8 stc
    
    step=0
    
    open(201, file = './data/Particles_0.plt')
    write(201,'(''TITLE = "All Particles Coordinates"'')')
	
    if (nWallAtom.ge. 1) then
        write(201,'(''ZONE'')')
        write(201,'(''T="WallPariticles"'')')
        do n = 1, nWallAtom
	        write(201,'(3f9.4)') (r(n,i), i = 1,3)
        enddo
    endif

    if(nDpEnd>nWallAtom) then
        write(201,'(''ZONE'')')
        write(201,'(''T="BubblePariticles"'')')
       do n = nWallAtom + 1, nDpEnd
	      write(201,'(3f9.4)') (r(n,i), i = 1,3)
       enddo
    endif

    if(nChain .gt. 0) then
        write(201,'(''ZONE'')')
        write(201,'(''T="ChainPariticles"'')')
       do n = nDpEnd + 1, nChainend
	      write(201,'(3f9.4)') (r(n,i), i = 1,3)
       enddo
       close(203)
    endif

    write(201,'(''ZONE'')')
    write(201,'(''T="FluidPariticles"'')')
	do n = nChainend + 1, nAtom
        write(201,'(3f9.4)') (r(n,i), i = 1,3)
    enddo


    write(201,'(''TEXT X=1 Y=90 T="   nAtom    nDp\\n'',2i8,'' "'')') &
                                     nAtom,nDpEnd-nWallAtom
    close(201)
    
 write(*,*) 'We have output the all particles location | stepCount:',step
 !pause
    
    end

!! Added by linyuqing 