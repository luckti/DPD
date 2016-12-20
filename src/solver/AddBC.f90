
!This file is used to add some special Boundary condition
    
    
    
!-----------------------------------------------------
!   Extract a pipe from an original cuboid
!-----------------------------------------------------

    subroutine PipeAdd

    implicit none
    include 'dpdflow.h'

    integer n, k
    integer pnf,nDf,nCf
    integer nwall, nfluid,nDfluid,nCfluid
    real*8 prad, prdw, prad2, prdw2, pryz, pryz2 
     
    real*8,dimension(:,:), allocatable :: trw,twn,trf,trv,tra, tDrf, tDrv, tDra,tCrf, tCrv, tCra
    
    
    pnf=nAtom-nChainend+nWallAtom+1
    nDf=nDpEnd-nWallAtom+1
    nCf=nChainend-nDpEnd+1
    
    allocate(trw(pnf, NDIM), twn(pnf, NDIM), trf(pnf, NDIM), trv(pnf, NDIM), tra(pnf, NDIM))
    allocate(tDrf(nDf, NDIM), tDrv(nDf, NDIM), tDra(nDf, NDIM))
    allocate(tCrf(nCf, NDIM), tCrv(nCf, NDIM), tCra(nCf, NDIM))
    
    prad=(min(region(2),region(3)))/2.-2.*gap(1)

    prdw = prad + 2.0*gap(1)
    
    prad2 = prad**2
    prdw2 = prdw**2
     
    nwall  = 0
    nfluid = 0
    nDfluid = 0
    nCfluid = 0
    
    do n = 1, nAtom
       pryz2 = r(n,2)**2 + r(n,3)**2
       if(pryz2 .ge. prad2 .and. pryz2 .le. prdw2) then
                if(n .gt. nWallAtom .and.n .lt. nDpEnd+1) then
                    print*,'error: the droplet contained in the pipe so that switched into wall  '
                    print*,'  n = ', n
                    print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
                    pause
                elseif(n .gt. nDpEnd .and.n .lt. nChainend+1)then
                    print*,'error: the chain contained in the pipe so that switched into wall  '
                    print*,'  n = ', n
                    print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
                    pause
                endif
              nwall = nwall + 1
              pryz = sqrt(pryz2)
              do k = 1, NDIM
                 trw(nwall,k) = r(n,k)
              enddo
              twn(nWall, 1) = 0.
              twn(nWall, 2) = -trw(nwall, 2)/pryz
              twn(nwall, 3) = -trw(nwall, 3)/pryz
       elseif(pryz2 .lt. prad2.and.(n .gt. nWallAtom .and. n.lt. nDpEnd+1)) then
           nDfluid = nDfluid + 1
          do k = 1, NDIM
             tDrf(nDfluid,k) = r(n,k)
             tDrv(nDfluid,k) = rv(n,k)
             tDra(nDfluid,k) = ra(n,k)
          enddo
       elseif(pryz2 .lt. prad2 .and.(n .gt. nDpEnd .and. n.lt. nChainend+1)) then
           nCfluid = nCfluid + 1
          do k = 1, NDIM
             tCrf(nCfluid,k) = r(n,k)
             tCrv(nCfluid,k) = rv(n,k)
             tCra(nCfluid,k) = ra(n,k)
          enddo    
       elseif(pryz2 .lt. prad2.and. n.gt.nChainend) then
          nfluid = nfluid + 1
          do k = 1, NDIM
             trf(nfluid,k) = r(n,k)
             trv(nfluid,k) = rv(n,k)
             tra(nfluid,k) = ra(n,k)
          enddo
       endif
    enddo

    
    nWallAtom=nwall
    nDpEnd=nWallAtom+nDfluid
    nChainend=nDpEnd+nCfluid
    nAtom=nChainend+nfluid
     

     do n =1, nAtom
        do k = 1, NDIM
             if (n.le.nWallAtom)then
                r(n, k)  = trw(n,k)
                wn(n, k) = twn(n,k)
                rv(n, k) = 0.
                ra(n, k) = 0.
             elseif(n.gt.nWallAtom.and.n.le.nDpEnd)then
                r(n, k)  = tDrf(n-nWallAtom,k)
                rv(n, k) = tDrv(n-nWallAtom,k)
                ra(n, k) = tDra(n-nWallAtom,k)
            elseif(n.gt.nDpEnd.and.n.le.nChainend)then
                r(n, k)  = tCrf(n-nDpEnd,k)
                rv(n, k) = tCrv(n-nDpEnd,k)
                ra(n, k) = tCra(n-nDpEnd,k)
            else
                r(n, k)  = trf(n-nChainend,k)
                rv(n, k) = trv(n-nChainend,k)
                ra(n, k) = tra(n-nChainend,k)
            endif
         enddo
     enddo
     
     nStartAtom = 1
     PNDIM = 1
    
!    call ouptParticleSituation
     
    ! nDpEnd = nWallAtom
    ! nChainend = nWallAtom

    ! write(30,'(//'' Coordinates of Wall Particles'')')
 	! write(30,'(2x, ''n'', 6x,''x'', 8x, ''y'', 8x, ''z'', 8x, &
!		  ''nx'', 7x, ''ny'', 7x, ''nz'')')
!	 do n = 1, nWallAtom
! enddo

 !	 write(30,'(//'' Coordinates of Simple Particles''/)')
	! do n = nChainend + 1, nAtom
	 !   write(30,'(i6, 3f9.4)') n, (r(n,i), i = 1,3)
	 !enddo

    deallocate(trw, twn, trf, trv, tra)
    deallocate(tDrf, tDrv, tDra)
    deallocate(tCrf, tCrv, tCra)
     
    end

!----------------------------------------------------------------------------------------------------------
    subroutine MidPlaneAdd

    implicit none
    include 'dpdflow.h'

    integer n, k
    integer pnf,nDf,nCf
    integer nwall, nfluid,nDfluid,nCfluid
    real*8 prad, prdw, prad2, prdw2, pryz, pryz2, bdp 
     
    real*8,dimension(:,:), allocatable :: trw,twn,trf,trv,tra, tDrf, tDrv, tDra,tCrf, tCrv, tCra
    
    
    pnf=nAtom-nChainend+nWallAtom+1
    nDf=nDpEnd-nWallAtom+1
    nCf=nChainend-nDpEnd+1
    
    allocate(trw(pnf, NDIM), twn(pnf, NDIM), trf(pnf, NDIM), trv(pnf, NDIM), tra(pnf, NDIM))
    allocate(tDrf(nDf, NDIM), tDrv(nDf, NDIM), tDra(nDf, NDIM))
    allocate(tCrf(nCf, NDIM), tCrv(nCf, NDIM), tCra(nCf, NDIM))
    
    bdp=2
    prad=RdsDp+bdp

    prdw = prad + 2.0*gap(3)
    
    prad2 = prad**2
    prdw2 = prdw**2
     
    nwall  = 0
    nfluid = 0
    nDfluid = 0
    nCfluid = 0
    
    do n = 1, nAtom
       pryz2 = r(n,3)**2
       if(pryz2 .ge. prad2 .and. pryz2 .le. prdw2 .and.r(n,3) > 0) then
                if(n .gt. nWallAtom .and.n .lt. nDpEnd+1) then
                    print*,'error: the droplet contained in the pipe so that switched into wall  '
                    print*,'  n = ', n
                    print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
                    pause
                elseif(n .gt. nDpEnd .and.n .lt. nChainend+1)then
                    print*,'error: the chain contained in the pipe so that switched into wall  '
                    print*,'  n = ', n
                    print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
                    pause
                endif
              nwall = nwall + 1
              pryz = sqrt(pryz2)
              do k = 1, NDIM
                 trw(nwall,k) = r(n,k)
              enddo
              twn(nWall, 1) = 0.
              twn(nWall, 2) = -trw(nwall, 2)/pryz
              twn(nwall, 3) = -trw(nwall, 3)/pryz
       elseif(n .gt. nWallAtom .and. n.lt. nDpEnd+1) then
           nDfluid = nDfluid + 1
          do k = 1, NDIM
             tDrf(nDfluid,k) = r(n,k)
             tDrv(nDfluid,k) = rv(n,k)
             tDra(nDfluid,k) = ra(n,k)
          enddo
       elseif(n .gt. nDpEnd .and. n.lt. nChainend+1) then
           nCfluid = nCfluid + 1
          do k = 1, NDIM
             tCrf(nCfluid,k) = r(n,k)
             tCrv(nCfluid,k) = rv(n,k)
             tCra(nCfluid,k) = ra(n,k)
          enddo    
       elseif(n.gt.nChainend) then
          nfluid = nfluid + 1
          do k = 1, NDIM
             trf(nfluid,k) = r(n,k)
             trv(nfluid,k) = rv(n,k)
             tra(nfluid,k) = ra(n,k)
          enddo
       endif
    enddo

    
    nWallAtom=nwall
    nDpEnd=nWallAtom+nDfluid
    nChainend=nDpEnd+nCfluid
    nAtom=nChainend+nfluid
     

     do n =1, nAtom
        do k = 1, NDIM
             if (n.le.nWallAtom)then
                r(n, k)  = trw(n,k)
                wn(n, k) = twn(n,k)
                rv(n, k) = 0.
                ra(n, k) = 0.
             elseif(n.gt.nWallAtom.and.n.le.nDpEnd)then
                r(n, k)  = tDrf(n-nWallAtom,k)
                rv(n, k) = tDrv(n-nWallAtom,k)
                ra(n, k) = tDra(n-nWallAtom,k)
            elseif(n.gt.nDpEnd.and.n.le.nChainend)then
                r(n, k)  = tCrf(n-nDpEnd,k)
                rv(n, k) = tCrv(n-nDpEnd,k)
                ra(n, k) = tCra(n-nDpEnd,k)
            else
                r(n, k)  = trf(n-nChainend,k)
                rv(n, k) = trv(n-nChainend,k)
                ra(n, k) = tra(n-nChainend,k)
            endif
         enddo
     enddo
     
     nStartAtom = 1
     PNDIM = 3
    
!    call ouptParticleSituation
     
    deallocate(trw, twn, trf, trv, tra)
    deallocate(tDrf, tDrv, tDra)
    deallocate(tCrf, tCrv, tCra)
     
    end

!----------------------------------------------------------------------------------------------------------
    subroutine DigCentreBubbleGap

    implicit none
    include 'dpdflow.h'

    integer n, k
    integer pnf,nDf,nCf,dbf
    integer nfluid,nDfluid,nCfluid
    real*8 bubbleR,pxyz2
    logical fexist
    
    real*8,dimension(:,:), allocatable :: trf,trv,tra, tDrf, tDrv, tDra,tCrf, tCrv, tCra
     
    pnf=nAtom-nChainend+1
    nDf=nDpEnd-nWallAtom+1
    nCf=nChainend-nDpEnd+1
    

    allocate(trf(pnf, NDIM), trv(pnf, NDIM), tra(pnf, NDIM))
    allocate(tDrf(nDf, NDIM), tDrv(nDf, NDIM), tDra(nDf, NDIM))
    allocate(tCrf(nCf, NDIM), tCrv(nCf, NDIM), tCra(nCf, NDIM))
    
1   nfluid = 0
    nDfluid = 0
    nCfluid = 0
    dbf  = 41
    
    inquire(file= './data/digBubbleRadius.dat',exist=fexist)
    if(fexist) then
        open(dbf,file= './data/digBubbleRadius.dat')
        read(dbf,*)
        read(dbf,*) bubbleR
        print*,'The bubble radius is:', bubbleR
    else
        write(*,*) 'please input the bubble radius which you want to dig :'
        read(*,*) bubbleR
        open(dbf,file= './data/digBubbleRadius.dat')
        write(dbf,*) 'Bubble radius'
        write(dbf,*) bubbleR
    endif
        
    do n = 1, nAtom
       pxyz2 = r(n,3)**2+r(n,2)**2+r(n,1)**2
      if(n .le. nChainend .and. pxyz2 .le. bubbleR*bubbleR) then
            close(dbf)            
            print*,'The bubble radius so big that contain the no fluid atom.' 
            write(*,*) 'please input the bubble radius which you want to dig :'
            open(dbf,file= './data/digBubbleRadius.dat')
            read(*,*) bubbleR
            write(dbf,*) 'Bubble radius'
            write(dbf,*) bubbleR
            goto 1
      elseif(n .gt. nWallAtom .and. n.lt. nDpEnd+1) then
           nDfluid = nDfluid + 1
          do k = 1, NDIM
             tDrf(nDfluid,k) = r(n,k)
             tDrv(nDfluid,k) = rv(n,k)
             tDra(nDfluid,k) = ra(n,k)
          enddo
       elseif(n .gt. nDpEnd .and. n.lt. nChainend+1) then
           nCfluid = nCfluid + 1
          do k = 1, NDIM
             tCrf(nCfluid,k) = r(n,k)
             tCrv(nCfluid,k) = rv(n,k)
             tCra(nCfluid,k) = ra(n,k)
          enddo    
       elseif(pxyz2 .ge. bubbleR*bubbleR .and. n.gt.nChainend) then
          nfluid = nfluid + 1
          do k = 1, NDIM
             trf(nfluid,k) = r(n,k)
             trv(nfluid,k) = rv(n,k)
             tra(nfluid,k) = ra(n,k)
          enddo
       endif
    enddo

    nDpEnd=nWallAtom+nDfluid
    nChainend=nDpEnd+nCfluid
    nAtom=nChainend+nfluid
     

     do n =1, nAtom
        do k = 1, NDIM
            if(n.gt.nWallAtom.and.n.le.nDpEnd)then
                r(n, k)  = tDrf(n-nWallAtom,k)
                rv(n, k) = tDrv(n-nWallAtom,k)
                ra(n, k) = tDra(n-nWallAtom,k)
            elseif(n.gt.nDpEnd.and.n.le.nChainend)then
                r(n, k)  = tCrf(n-nDpEnd,k)
                rv(n, k) = tCrv(n-nDpEnd,k)
                ra(n, k) = tCra(n-nDpEnd,k)
            elseif(n.gt.nChainend) then
                r(n, k)  = trf(n-nChainend,k)
                rv(n, k) = trv(n-nChainend,k)
                ra(n, k) = tra(n-nChainend,k)
            endif
         enddo
     enddo

     nStartAtom = nWallAtom
     PNDIM = 0
    
!    call ouptParticleSituation
     
    deallocate(trf, trv, tra)
    deallocate(tDrf, tDrv, tDra)
    deallocate(tCrf, tCrv, tCra)
    close(dbf)

    end