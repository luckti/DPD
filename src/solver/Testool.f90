!! Added by linyuqing

    subroutine ouptParticleSituation
    ! To output the coordination, velocity and acceleration of every particle.
    
    implicit none
    include 'dpdflow.h'
    integer ifirst
    character*8 stc
    integer n,i
    
    write(stc,'(i8)')stepCount
    if(stepCount<10) then
        ifirst=index(stc,' ')+7
    elseif(stepCount<100) then
        ifirst=index(stc,' ')+6
    elseif(stepCount<1000) then
        ifirst=index(stc,' ')+5
    elseif(stepCount<10000) then
        ifirst=index(stc,' ')+4
    elseif(stepCount<100000) then
        ifirst=index(stc,' ')+3
    elseif(stepCount<1000000) then
        ifirst=index(stc,' ')+2   
    elseif(stepCount<10000000) then
        ifirst=index(stc,' ')+1 
    endif
    
    
	open(201, file = './data/instantState/Particles_'//stc(ifirst:)//'.plt')    
    write(201,'(''TITLE = "All Particles Coordinates"'')')
    write(201,'(''Variables= "X","Y","Z","Vx","Vy","Vz","ax","ay","az"'')')
	
    if (nWallAtom.ge. 1 .and. runId.ne.0.and. runId.ne.6) then
        write(201,'(''ZONE'')')
        write(201,'(''T="WallPariticles"'')')
    write(201,'(''I= '',i8,'',F = POINT'')') nWallAtom
        do n = 1, nWallAtom
	        write(201,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
        enddo
    endif

    if(nDpEnd>nWallAtom) then
        write(201,'(''ZONE'')')
        write(201,'(''T="BubblePariticles"'')') 
    write(201,'(''I= '',i8,'',F = POINT'')') nDpEnd-nWallAtom
       do n = nWallAtom + 1, nDpEnd
	      write(201,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
       enddo
    endif

    if(nChain .gt. 0) then
        write(201,'(''ZONE'')')
        write(201,'(''T="ChainPariticles"'')')
    write(201,'(''I= '',i8,'',F = POINT'')') nChainend-nDpEnd
       do n = nDpEnd + 1, nChainend
	      write(201,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
       enddo
       close(203)
    endif

    write(201,'(''ZONE'')')
    write(201,'(''T="FluidPariticles"'')')
    write(201,'(''I= '',i8,'',F = POINT'')') nAtom-nChainend
	do n = nChainend + 1, nAtom
        write(201,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
    enddo


    write(201,'(''TEXT X=1 Y=90 T="   nAtom    nDp\\n'',2i8,'' "'')') &
                                     nAtom,nDpEnd-nWallAtom
    close(201)
 write(*,*) 'We have output the all particles location | stepCount:',stepCount
 !pause
    
    end

!! Added by linyuqing

! --------------------------------------------------------------------------------------------------------------------
subroutine ouptParticleGap(icord)
    ! To output the coordination, velocity and acceleration of every particle.
    
    implicit none
    include 'dpdflow.h'
    integer n,i,icord
    
    if (icord .eq. 0) then
	    open(202, file = './data/GapParticles.plt')    
        write(202,'(''TITLE = "All Particles Coordinates"'')')
        write(202,'(''Variables= "X","Y","Z","Vx","Vy","Vz","ax","ay","az"'')')
    endif
    
    write(202,'(''ZONE'')')
    write(202,'(''T="'',i8,''"'')') stepCount
    write(202,'(''I= '',i8,'',F = POINT'')') nAtom
	do n = 1, nAtom
        write(202,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
    enddo

    write(202,'(''TEXT X=1 Y=90 T="   nAtom\\n'',i8,'' "'')') nAtom
    
    if (icord .eq. -1) close(202)

 !pause
    
    end

!! Added by linyuqing

! --------------------------------------------------------------------------------------------------------------------
subroutine ouptSameParticlesAvg(icord)
    ! To output the coordination, velocity and acceleration of every particle.
    
    implicit none
    include 'dpdflow.h'
    integer n,i,icord
    real*8 Cr(3),Crv(3),Cra(3)
    
    if (icord .eq. 0) then 
        open(205, file = './data/Fluid.plt')  
        write(205,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
        write(205,'(''ZONE'')')
        write(205,'(''F = POINT'')') 
        if(nDp.eq. 1.and.nDpEnd>nWallAtom) then
	        open(203, file = './data/Bubble.plt')  
            write(203,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
            write(203,'(''ZONE'')')
            write(203,'(''F = POINT'')') 
        endif
        if(nChain .gt. 0) then
            open(204, file = './data/Chain.plt')  
            write(204,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
            write(204,'(''ZONE'')')
            write(204,'(''F = POINT'')') 
        endif
    endif
    
    do i =1,3
        Cr(i)=0.
        Crv(i)=0.
        Cra(i)=0.
    enddo
	do n = nChainend+1, nAtom
        do i =1,3
            Cr(i)=Cr(i)+r(n,i)
            Crv(i)=Crv(i)+rv(n,i)
            Cra(i)=Cra(i)+ra(n,i)
        enddo
    enddo
    do i =1,3
        Cr(i)=Cr(i)/(nAtom-nChainend)
        Crv(i)=Crv(i)/(nAtom-nChainend)
        Cra(i)=Cra(i)/(nAtom-nChainend)
    enddo
    write(205,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
    if (icord .eq. -1)  close(205)
    
    if(nDp.eq. 1.and.nDpEnd>nWallAtom) then
        do i =1,3
            Cr(i)=0.
            Crv(i)=0.
            Cra(i)=0.
        enddo
	    do n = nWallAtom+1, nDpEnd
            do i =1,3
                Cr(i)=Cr(i)+r(n,i)
                Crv(i)=Crv(i)+rv(n,i)
                Cra(i)=Cra(i)+ra(n,i)
            enddo
        enddo
        do i =1,3
            Cr(i)=Cr(i)/(nDpEnd-nWallAtom)
            Crv(i)=Crv(i)/(nDpEnd-nWallAtom)
            Cra(i)=Cra(i)/(nDpEnd-nWallAtom)
        enddo
        write(203,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
        if (icord .eq. -1)  close(203)
    endif

    if(nChain .gt. 0) then
        do i =1,3
            Cr(i)=0.
            Crv(i)=0.
            Cra(i)=0.
        enddo
	    do n = nDpEnd+1, nChainend
            do i =1,3
                Cr(i)=Cr(i)+r(n,i)
                Crv(i)=Crv(i)+rv(n,i)
                Cra(i)=Cra(i)+ra(n,i)
            enddo
        enddo
        do i =1,3
            Cr(i)=Cr(i)/(nChainend-nDpEnd)
            Crv(i)=Crv(i)/(nChainend-nDpEnd)
            Cra(i)=Cra(i)/(nChainend-nDpEnd)
        enddo
        write(204,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
        if (icord .eq. -1)  close(204)
    endif  
    
 !pause
    
    end


!! Added by linyuqing
! --------------------------------------------------------------------------------------------------------------------
subroutine ouptRDF(icord)
    ! To output the coordination, velocity and acceleration of every particle.
    
    implicit none
    include 'dpdflow.h'
    integer m,n,i,j,icord,Sn,judgeg0(3),judgel0(3),ng0(3),nl0(3),nC
    real*8 Cr(3),rnC(3),vr(3),rr,Crg0(3),Crl0(3),dC,dd,rtemp,Vtemp(3)
    real*8 Ptemp,Sr(6),temp(3,3),SQMatr(3,3),TensorMatr(3,3),&
    		Srv(6),Sra(6),Srforce(6)

    real*8, allocatable :: nG(:,:),nGb(:,:),nGp(:,:)
    integer ic(3),k,Bpout,BnGridH
    
    !deltaQ=0.05
    Sn=int(min(regionH(3),regionH(2),regionH(1))/deltaQ)    
    
    if (icord .eq. 0) then 
        open(301, file = './data/RDF/PDF.plt') 
        open(305, file = './data/RDF/RDF.plt')  

        open(302, file = './data/RDF/BCentre.plt')  
        write(302,'(''Variables= "stepCount","BCrX","BCrY","BCrZ","PCrX","PCrY","PCrZ","nS1","nS2","nS3"'')')
        write(302,'(''ZONE'')')
        write(302,'(''F = POINT'')') 

        open(303, file = './data/RDF/rdflog.dat')  

        open(304 , file = './data/RDF/RDFcell.plt' , status = 'unknown')
	    write(304, *) 'Title=" "'
	    write(304, *) 'Variables= "X","Y","Z","Rho","BRho","FRho"'
	    write(304, *) 'zone'
	    write(304, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', BnGrid, ',J=', BnGrid, ',K=', BnGrid, ',F=POINT'
        
        open(306 , file = './data/RDF/StateCell.plt' , status = 'unknown')
	    write(306, *) 'Title=" "'
	    write(306, *) 'Variables= "X","Y","Z","Vx","Vy","Vz","ax","ay","az"'
	    write(306, *) 'zone'
	    write(306, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', BnGrid, ',J=', BnGrid, ',K=', BnGrid, ',F=POINT'
        
        open(307 , file = './data/RDF/StrCell.plt' , status = 'unknown')
	    write(307, *) 'Title=" "'
	    write(307, *) 'Variables= "X","Y","Z","Sxx","Syy","Szz","Sxy","Sxz","Syz","Pressure"'
	    write(307, *) 'zone'
	    write(307, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', BnGrid, ',J=', BnGrid, ',K=', BnGrid, ',F=POINT'
        
        open(308, file = './data/RDF/PRDF.plt')

        open(309, file = './data/RDF/STRDF.plt') 
!if(stepCount .ge. 600) print*, 'v0'     !test               
        allocate (nS(1:Sn,3),nSB(1:Sn,3))
        allocate (StateSpGrid(1:Sn,1:6),StrSpGrid(1:Sn,1:6))
        do j=1,Sn
        	do k=1,3
            	nS(j,k)=0
            	nSB(j,k)=0
            enddo
            do k=1,6
	            StrSpGrid(j,k) = 0.
	            StateSpGrid(j,k) = 0.
       		enddo 
        enddo

        if (nDp .eq. 1) then
		    BnGrid = int(min(region(3),region(2),region(1))/(cellength))
		    if(mod(BnGrid,2) .eq. 0) BnGrid = BnGrid + 1
		    BnGridH = (BnGrid-1)/2
		    allocate (BmGrid(-BnGridH:BnGridH,-BnGridH:BnGridH,-BnGridH:BnGridH,1:3))
	        allocate (StrGrid(-BnGridH:BnGridH,-BnGridH:BnGridH,-BnGridH:BnGridH,1:6))
	        allocate (StateGrid(-BnGridH:BnGridH,-BnGridH:BnGridH,-BnGridH:BnGridH,1:6))
	    
	        do i = -BnGridH,BnGridH,1
			    do j = -BnGridH,BnGridH,1
				    do k = -BnGridH,BnGridH,1
					    do n = 1,3
						    BmGrid(i,j,k,n) = 0.    			!initial value	                    
	                    enddo
	                    do n=1,6
	                        StrGrid(i,j,k,n) = 0.
	                        StateGrid(i,j,k,n) = 0.
	                    enddo
				    enddo
		        enddo           
	        enddo 
    	endif
!if(stepCount .ge. 600) print*, 'v1'     !test
         
    endif

    BnGridH = (BnGrid-1)/2

        
    if (icord .eq. 1) then
    	judgeg0=0
    	judgel0=0
        do i =1,3
            Cr(i)=0.
            Crg0(i)=0.
            Crl0(i)=0.
            ng0(i)=0
            nl0(i)=0
            ic(i) =0 
        enddo
!if(stepCount .ge. 600) print*, 'v2'     !test
        do n = nWallAtom+1, nDpEnd
            do i =1,3
            	if (r(n,i) .ge.0) then
            		Crg0(i)=Crg0(i)+r(n,i)            		
            		ng0(i)=ng0(i)+1
            		if(r(n,i) .gt. (regionH(i)/10) .and.judgeg0(i) .eq. 0) judgeg0(i)=1
            	else
            		Crl0(i)=Crl0(i)+r(n,i)
            		nl0(i)=nl0(i)+1
            		if(r(n,i) .gt. (regionH(i)/10) .and.judgel0(i) .eq. 0) judgel0(i)=1
            	endif
            enddo
        enddo
!if(stepCount .ge. 600) print*, 'v3'     !test
        do i=1,3
	        if (ng0(i).gt. 0 .and. nl0(i) .gt. 0 .and. (judgeg0(i)+judgel0(i)).eq. 2) then        	
	    		if (ng0(i) .ge. nl0(i)) Crl0(i)=Crl0(i)+nl0(i)*region(i)
	    		if (ng0(i) .lt. nl0(i)) Crg0(i)=Crg0(i)-ng0(i)*region(i) 
	    		write(303,'(''The bubble has acrossed the border'',i5)') stepCount      	
	        endif
	        	Cr(i)=Crg0(i)+Crl0(i)
	            Cr(i)=Cr(i)/(nDpEnd-nWallAtom)
        enddo
!if(stepCount .ge. 600) print*, 'v4'     !test
        ! To find the center particle which is closest to the center point.
        nC= nWallAtom+1
        dC= 0.
        do i=1,3
            dC=dC+sqrt((r(nC,i)-Cr(i))**2)
        enddo
        do n = nWallAtom+2, nDpEnd
            dd=0.
            do i=1,3
                dd=dd+sqrt((r(n,i)-Cr(i))**2)
            enddo
            if(dd .lt. dC) then 
                dC=dd
                nC=n
            endif
        enddo
!if(stepCount .ge. 600) print*, 'v5'     !test 
        do i=1,3
            rnC(i)=r(nC,i)
            if (abs(r(nC,i)) .gt. regionH(i))  rnC(i)=r(nC,i)-(r(nC,i)/abs(r(nC,i)))*region(i) 
        enddo
!if(stepCount .ge. 600) print*, 'v6'     !test       		
        do n = nWallAtom+1,nAtom
!if(stepCount .ge. 600) print*, 'v1',n    !test    
        	!Bubble PDF -------------------------------------------------------------
            rr=0.
            do i =1,3                 
                vr(i) = abs(r(n,i)-rnC(i))
                if (vr(i) .gt. regionH(i)) vr(i)= region(i) - vr(i)
                rr=rr+(vr(i)*vr(i))
            enddo
            
            j=int(sqrt(rr)/deltaQ)+1
            if (j .le. Sn) then 
            	nS(j,1)=nS(j,1)+1   
            	if (n .le. nDpEnd) then
                	nS(j,2)=nS(j,2)+1 
                else 
                	nS(j,3)=nS(j,3)+1 
                endif 
            endif
!if(stepCount .ge. 600) print*, 'v2',n    !test 
            !Bubble RDF --------------------------------------------------------------------
            rr=0.
            do i =1,3                 
                vr(i) = abs(r(n,i)-Cr(i))
                if (vr(i) .gt. regionH(i)) vr(i)= region(i) - vr(i)
                rr=rr+(vr(i)*vr(i))
            enddo
            
            j=int(sqrt(rr)/deltaQ)+1
            if (j .le. Sn) then 
            	nSB(j,1)=nSB(j,1)+1   
            	if (n .le. nDpEnd) then
                	nSB(j,2)=nSB(j,2)+1 
                else 
                	nSB(j,3)=nSB(j,3)+1 
                endif 
				!Spherial coordination stress tensor----------------------------------------------------
				Sr(1) = sqrt(r(n,1)*r(n,1)+r(n,2)*r(n,2)+r(n,3)*r(n,3))
                if (Sr(1).ne.0) then
                    Sr(2) = acos(r(n,3)/Sr(1))
                    if (r(n,1).ne. 0 .and. r(n,2) .ne. 0) then
                        Sr(3) = atan(r(n,2)/r(n,1))
                    else
                        Sr(3) = 0.
                    endif
                else
                    Sr(2) = 0.
                    Sr(3) = 0.
                endif
                
                SQMatr(1,1)=sin(Sr(2))*cos(Sr(3))
                SQMatr(2,1)=sin(Sr(3))*sin(Sr(2))
                SQMatr(3,1)=cos(Sr(2))
                SQMatr(1,2)=cos(Sr(2))*cos(Sr(3))
                SQMatr(2,2)=sin(Sr(3))*cos(Sr(2))
                SQMatr(3,2)=-sin(Sr(2))
                SQMatr(1,2)=-sin(Sr(3))
                SQMatr(2,2)=cos(Sr(3))
                SQMatr(3,2)=0.

                do i=1,3
                	Srv(i) = 0.
                	Sra(i) = 0.
                	do k =1,3 
                		Srv(i) = Srv(i)+SQMatr(k,i)*rv(n,k)
                		Sra(i) = Sra(i)+SQMatr(k,i)*ra(n,k)
					enddo
				enddo

				TensorMatr(1,1) = rforce(n,1)
                TensorMatr(2,2) = rforce(n,2)
                TensorMatr(3,3) = rforce(n,3)
                TensorMatr(1,2) = rforce(n,4)
                TensorMatr(1,3) = rforce(n,5)
                TensorMatr(2,3) = rforce(n,6)
                TensorMatr(2,1) = rforce(n,4)
                TensorMatr(3,1) = rforce(n,5)
                TensorMatr(3,2) = rforce(n,6)
                
                do m=1,3
                    do i=1,3
                        temp(m,i)=SQMatr(i,m)
                    enddo
                enddo
                                    
                TensorMatr =  matmul(SQMatr,TensorMatr)
                TensorMatr =  matmul(TensorMatr,temp)

				Srforce(1) = TensorMatr(1,1)
                Srforce(2) = TensorMatr(2,2)
                Srforce(3) = TensorMatr(3,3)
                Srforce(4) = TensorMatr(1,2)
                Srforce(5) = TensorMatr(1,3)
                Srforce(6) = TensorMatr(2,3)

                StrSpGrid(j,1) = StrSpGrid(j,1) + mass*Srv(1)**2 + 0.5*Srforce(1)
                StrSpGrid(j,2) = StrSpGrid(j,2) + mass*Srv(2)**2 + 0.5*Srforce(2)
                StrSpGrid(j,3) = StrSpGrid(j,3) + mass*Srv(3)**2 + 0.5*Srforce(3)
                StrSpGrid(j,4) = StrSpGrid(j,4) + mass*Srv(1)*Srv(2) + 0.5*Srforce(4)
                StrSpGrid(j,5) = StrSpGrid(j,5) + mass*Srv(1)*Srv(3) + 0.5*Srforce(5)
                StrSpGrid(j,6) = StrSpGrid(j,6) + mass*Srv(2)*Srv(3) + 0.5*Srforce(6)
                !Spherial coordination state velocity and acceleration 
                do i=1,3
                    StateSpGrid(j,i) = StateSpGrid(j,i) + Srv(i)              
                    StateSpGrid(j,i+3) = StateSpGrid(j,i+3) + Sra(i) 
                enddo          

            endif
!if(stepCount .ge. 600) print*, 'v3',n    !test             
            !/*Grid way-------------------------------------------------------------------
           	if (rr .le. 3*((((BnGridH+1)*cellength)**2))) then
           		do i=1,3
	           		rtemp = (vr(i)*abs(r(n,i)-Cr(i)))/((r(n,i)-Cr(i))*cellength)
	           		rtemp = rtemp + (rtemp/(abs(rtemp)))*0.5
	           		ic(i)=int(rtemp)                    
                enddo
           		if(abs(ic(1)) .gt. BnGridH .or. abs(ic(2)) .gt. BnGridH .or. abs(ic(3)) .gt. BnGridH) then
           			if (n .le. nDpEnd) then
	           			write(303,'(''The bubble was two big'')') 
	           			do i =1,3
	           			if(abs(ic(i)) .gt. BnGridH)write(303,'(i5,2f10.4,i4,3i8)') stepCount,Cr(i),r(n ,i), ic(i),i, n,BnGridH
	           			enddo
           			endif
                else
	                BmGrid(ic(1),ic(2),ic(3),1) = BmGrid(ic(1),ic(2),ic(3),1) + 1. 
	                if (n .le. nDpEnd) then
	                	BmGrid(ic(1),ic(2),ic(3),2) = BmGrid(ic(1),ic(2),ic(3),2) + 1. 
	                else 
	                	BmGrid(ic(1),ic(2),ic(3),3) = BmGrid(ic(1),ic(2),ic(3),3) + 1. 
                    endif
                    
                    !stress tensor----------------------------------------------------
                    StrGrid(ic(1),ic(2),ic(3),1) = StrGrid(ic(1),ic(2),ic(3),1) + mass*rv(n,1)**2 + 0.5*rforce(n,1)
                    StrGrid(ic(1),ic(2),ic(3),2) = StrGrid(ic(1),ic(2),ic(3),2) + mass*rv(n,2)**2 + 0.5*rforce(n,2)
                    StrGrid(ic(1),ic(2),ic(3),3) = StrGrid(ic(1),ic(2),ic(3),3) + mass*rv(n,3)**2 + 0.5*rforce(n,3)
                    StrGrid(ic(1),ic(2),ic(3),4) = StrGrid(ic(1),ic(2),ic(3),4) + mass*rv(n,1)*rv(n,2) + 0.5*rforce(n,4)
                    StrGrid(ic(1),ic(2),ic(3),5) = StrGrid(ic(1),ic(2),ic(3),5) + mass*rv(n,1)*rv(n,3) + 0.5*rforce(n,5)
                    StrGrid(ic(1),ic(2),ic(3),6) = StrGrid(ic(1),ic(2),ic(3),6) + mass*rv(n,2)*rv(n,3) + 0.5*rforce(n,6)
                    !Grid state velocity and acceleration 
                    do i=1,3
                        StateGrid(ic(1),ic(2),ic(3),i) = StateGrid(ic(1),ic(2),ic(3),i) + rv(n,i)              
                        StateGrid(ic(1),ic(2),ic(3),i+3) = StateGrid(ic(1),ic(2),ic(3),i+3) + ra(n,i) 
                    enddo                    
                endif
            endif 
            !\*Grid way--------------------------------------------------------------------     
        enddo
if(nS(1,1) .ne. (nS(1,2)+nS(1,3))) write(303,'(''nS'',4i10)') stepCount,(nS(j,1), j = 1,3)    !test 
if(nSB(1,1) .ne. (nSB(1,2)+nSB(1,3))) write(303,'(''nSB'',4i10)') stepCount ,(nSB(j,1), j = 1,3)   !test 
        write(302,'(i10,6f10.4,3i10)') stepCount,(Cr(i), i = 1,3),(rnC(j), j = 1,3),(nS(j,1), j = 1,3)

    endif

    if (icord .eq. 2) then 


        write(309,'(''Variables= "Radius","STRDF","nParticle"'')')
        write(309,'(''ZONE'')')
        write(309,'(''F = POINT'')')
        do j=1,Sn
        	if (nSB(j,1) .gt. 0) then
	        	do i=1,6
	        		!print*,StrSpGrid(j,i),StateSpGrid(j,i)        		
		        	StrSpGrid(j,i) = StrSpGrid(j,i) / nSB(j,1)
	                StateSpGrid(j,i) = StateSpGrid(j,i) / nSB(j,1)         
	            enddo 
        	endif

            StrSpGrid(j,1) = StrSpGrid(j,1) - mass * StateSpGrid(j,1)**2
            StrSpGrid(j,2) = StrSpGrid(j,2) - mass * StateSpGrid(j,2)**2
            StrSpGrid(j,3) = StrSpGrid(j,3) - mass * StateSpGrid(j,3)**2
            StrSpGrid(j,4) = StrSpGrid(j,4) - mass * StateSpGrid(j,1)*StateSpGrid(j,2)
            StrSpGrid(j,5) = StrSpGrid(j,5) - mass * StateSpGrid(j,1)*StateSpGrid(j,3)
            StrSpGrid(j,6) = StrSpGrid(j,6) - mass * StateSpGrid(j,2)*StateSpGrid(j,3)
            
            do i=1,6
                if (j .eq. 1) then
	                StrSpGrid(j,i)=StrSpGrid(j,i)*nSB(j,1)/(RDFstepSample*4*pi*(deltaQ**3)/3)
	            else
	                StrSpGrid(j,i)=StrSpGrid(j,i)*nSB(j,1)/(RDFstepSample*4*pi*((j**3)-((j-1)**3))*(deltaQ**3)/3)
	        	endif         
            enddo 

            write(309,'(2f10.4,i10)') j*deltaQ,(StrSpGrid(j,1)-(StrSpGrid(j,2)-StrSpGrid(j,3))/2),nSB(j,1)
        enddo


        
        allocate (nG(1:Sn,3))
       
        do n =1,3
	        write(301,'(''Variables= "Radius","PDF","nParticle"'')')
	        write(301,'(''ZONE'')')
	        write(301,'(''F = POINT'')')
	        do j=1,Sn
	        	nG(j,n)= 0.
	            if (j .eq. 1) then
	                nG(j,n)=(1.0*nS(j,n)/RDFstepSample)
	                if (n .ne. 3) nG(j,n)=nG(j,n)-1.
	            	nG(j,n)=nG(j,n)/(4*pi*(deltaQ**3)/3)
	            else
	                nG(j,n)=1.0*nS(j,n)/RDFstepSample
	            	nG(j,n)=nG(j,n)/(4*pi*((j**3)-((j-1)**3))*(deltaQ**3)/3)
	        	endif
	            write(301,'(2f10.4,i10)') j*deltaQ,nG(j,n),nS(j,n)
	        enddo
	    enddo

	    allocate (nGB(1:Sn,3))
       
        do n =1,3
	        write(305,'(''Variables= "Radius","RDF","nParticle"'')')
	        write(305,'(''ZONE'')')
	        write(305,'(''F = POINT'')')
	        do j=1,Sn
	        	nGB(j,n)= 0.
	            if (j .eq. 1) then
	                nGB(j,n)=(1.0*nSB(j,n)/RDFstepSample)
	            	nGB(j,n)=nGB(j,n)/(4*pi*(deltaQ**3)/3)
	            else
	                nGB(j,n)=1.0*nSB(j,n)/RDFstepSample
	            	nGB(j,n)=nGB(j,n)/(4*pi*((j**3)-((j-1)**3))*(deltaQ**3)/3)
	        	endif
	            write(305,'(2f10.4,i10)') j*deltaQ,nGB(j,n),nSB(j,n)
	        enddo
        enddo
        
        Sn=int(BnGridH*cellength/deltaQ)
        allocate (nGp(1:Sn,4))
        
        write(308,'(''Variables= "Radius","Pressure","SurTensor","nGrid","xSurTensor"'')')
	    write(308,'(''ZONE'')')
	    write(308,'(''F = POINT'')')    
        do j=1,Sn
	        nGp(j,1)= 0.
            nGp(j,2)= 0.
            nGp(j,3)= 0.
            nGp(j,4)= 0.
        enddo
        do m=1,3
            do n =1,3
                SQMatr(m,n)=0.
                TensorMatr(m,n)=0.
            enddo
        enddo
        
        do i = -BnGridH,BnGridH,1
        	do j = -BnGridH,BnGridH,1
        		do k = -BnGridH,BnGridH,1                    
                    !Average Grid state and stress
                    if (BmGrid(i,j,k,1) .gt. 0) then
                        do n=1,6
                            StateGrid(i,j,k,n) = StateGrid(i,j,k,n)/BmGrid(i,j,k,1)
                            StrGrid(i,j,k,n) = StrGrid(i,j,k,n)/BmGrid(i,j,k,1)
                        enddo                    
                        StrGrid(i,j,k,1) = StrGrid(i,j,k,1) - mass * StateGrid(i,j,k,1)**2
                        StrGrid(i,j,k,2) = StrGrid(i,j,k,2) - mass * StateGrid(i,j,k,2)**2
                        StrGrid(i,j,k,3) = StrGrid(i,j,k,3) - mass * StateGrid(i,j,k,3)**2
                        StrGrid(i,j,k,4) = StrGrid(i,j,k,4) - mass * StateGrid(i,j,k,1)*StateGrid(i,j,k,2)
                        StrGrid(i,j,k,5) = StrGrid(i,j,k,5) - mass * StateGrid(i,j,k,1)*StateGrid(i,j,k,3)
                        StrGrid(i,j,k,6) = StrGrid(i,j,k,6) - mass * StateGrid(i,j,k,2)*StateGrid(i,j,k,3)
                    endif
                    !get the stress
                    do n=1,6
                        StrGrid(i,j,k,n) = -StrGrid(i,j,k,n)*BmGrid(i,j,k,1)/(RDFstepSample*(cellength**3))
                    enddo
                    Ptemp= -(StrGrid(i,j,k,1) + StrGrid(i,j,k,2) + StrGrid(i,j,k,3))/3.0  !get the pressure
                    !get the density
        			do n = 1,3
	        			BmGrid(i,j,k,n) = BmGrid(i,j,k,n)/(RDFstepSample*(cellength**3))
        			enddo
        			write(304,'(3f10.4,3f20.9)') i*cellength,j*cellength,k*cellength,(BmGrid(i,j,k,n), n=1,3)
                    write(306,'(3f10.4,6f20.9)') i*cellength,j*cellength,k*cellength,(StateGrid(i,j,k,n), n=1,6)
                    write(307,'(3f10.4,7f20.9)') i*cellength,j*cellength,k*cellength,(StrGrid(i,j,k,n), n=1,6),Ptemp
                    !get the pressure RDF
!                    print*,stepCount,n
                    Sr(1) = sqrt((i*i+j*j+k*k)*1.0)*cellength
                    if (Sr(1).ne.0) then
                        Sr(2) = acos((k*cellength)/Sr(1))
                        if (j .ne. 0 .and. i .ne. 0) then
                            Sr(3) = atan(j*1.0/i)
                        else
                            Sr(3) = 0.
                        endif
                    else
                        Sr(2) = 0.
                        Sr(3) = 0.
                    endif
                    
                    SQMatr(1,1)=sin(Sr(2))*cos(Sr(3))
                    SQMatr(2,1)=sin(Sr(3))*sin(Sr(2))
                    SQMatr(3,1)=cos(Sr(2))
                    SQMatr(1,2)=cos(Sr(2))*cos(Sr(3))
                    SQMatr(2,2)=sin(Sr(3))*cos(Sr(2))
                    SQMatr(3,2)=-sin(Sr(2))
                    SQMatr(1,2)=-sin(Sr(3))
                    SQMatr(2,2)=cos(Sr(3))
                    SQMatr(3,2)=0.
                    
                    TensorMatr(1,1)=-StrGrid(i,j,k,1)
                    TensorMatr(2,2)=-StrGrid(i,j,k,2)
                    TensorMatr(3,3)=-StrGrid(i,j,k,3)
                    TensorMatr(1,2)=-StrGrid(i,j,k,4)
                    TensorMatr(1,3)=-StrGrid(i,j,k,5)
                    TensorMatr(2,3)=-StrGrid(i,j,k,6)
                    TensorMatr(2,1)=-StrGrid(i,j,k,4)
                    TensorMatr(3,1)=-StrGrid(i,j,k,5)
                    TensorMatr(3,2)=-StrGrid(i,j,k,6)
                    do m=1,3
                        do n=1,3
                            temp(m,n)=SQMatr(n,m)
                        enddo
                    enddo
                                        
                    TensorMatr =  matmul(SQMatr,TensorMatr)
                    TensorMatr =  matmul(TensorMatr,temp)
                    
                    Sr(4) = TensorMatr(1,1)
                    Sr(5) = TensorMatr(2,2)
                    Sr(6) = TensorMatr(3,3)         

                    m=int(Sr(1)/deltaQ)+1
                    if(m .le. Sn) then
                        nGp(m,1)= nGp(m,1)+Ptemp
                        nGp(m,2)= nGp(m,2)+Sr(4)-(Sr(5)+Sr(6))/2
                        nGp(m,3)= nGp(m,3)+1
                        nGp(m,4)= nGp(m,4)+StrGrid(i,j,k,1)-(StrGrid(i,j,k,2)+StrGrid(i,j,k,3))/2
                    endif
!                    print*,nGp(m,2),(Sr(n),n=1,6)
        		enddo
            enddo         
        enddo   
        
        !get the pressure RDF
        do j=1,Sn
            if(nGp(j,3) .gt. 0) then
                nGp(j,1)=nGp(j,1)/nGp(j,3)
                nGp(j,2)=nGp(j,2)/nGp(j,3)
                nGp(j,4)=nGp(j,4)/nGp(j,3)
            endif
	        write(308,'(f10.4,2f16.6,i10,f16.6)') j*deltaQ,nGp(j,1),nGp(j,2),int(nGp(j,3)),nGp(j,4)
	    enddo
            
       
        deallocate (nS,nSB)
        deallocate (nG,nGB,nGp)
        deallocate (BmGrid,StateGrid,StrGrid)
        deallocate (StateSpGrid,StrSpGrid)
        close(301)
        close(302)
        close(303)
        close(304)
        close(305)
        close(306)
        close(307)
        close(308)
        close(309)
    endif

    
 !pause
    
    end

! !! Added by linyuqing -------------------------------------------------------------------------------------------
    
! subroutine DiffusionCompute(icord)
!     ! To output the coordination, velocity and acceleration of every particle.
    
!     implicit none
!     include 'dpdflow.h'

!     integer n,i,icord
!     real * 8 Dmsd,time, DiffStartStep !Dgk,
    
!     real * 8, pointer::ZeroSate(:,:)
!     common / Diff / ZeroSate
!     common / Diffstep / DiffStartStep
    
! !    DiffintStep,DiffStep
    
    
!     if (icord .eq. 0) then
        
!         allocate (ZeroSate(nAtom - nWallAtom,NDIM*2))
!         do n = nWallAtom+1,nAtom
!             do i= 1, NDIM
!                 ZeroSate(n - nWallAtom,i) = r(n,i)+ ncc(n,i) * initUcell(i)*gap(i)
! !                ZeroSate(n - nWallAtom,i+NDIM) = rv(n,i)
!             enddo
!         enddo
!         DiffStartStep = stepCount
!         time = 0.
!         Dmsd = 0.
! !        Dgk = 0.
        
!         open(401, file = './data/CG/Diffusivity.plt')    
!         write(401,'(''TITLE = "Diffusivity('',i6,'')"'')') stepCount
!         write(401,'(''Variables= "t","D(MSD)"'')')  !,"D(GK)"
!         write(401,'(f10.4,2E15.7)') time,Dmsd !,Dgk
!     endif

!     if (icord .eq. 1) then
!         time = (stepCount-DiffStartStep)*deltaT
!         Dmsd = 0.
! !        Dgk = 0.
        
!         do n = nWallAtom+1,nAtom
!             do i= 1, NDIM
!                 Dmsd = Dmsd + (r(n,i)+ ncc(n,i) * initUcell(i)*gap(i) - ZeroSate(n - nWallAtom,i))**2
! !                Dgk = Dgk + ZeroSate(n - nWallAtom,i+NDIM)*rv(n,i)
!             enddo
!         enddo
!         Dmsd = Dmsd /(nAtom - nWallAtom)
! !        Dgk = Dgk /(nAtom - nWallAtom)
!         write(401,'(f10.4,2E15.7)') time,Dmsd !,Dgk
!     endif
    
!     if (icord .eq. 2) then
!         deallocate (ZeroSate)
!         close(401)
!     endif
    
    
!     end

!! Added by linyuqing -------------------------------------------------------------------------------------   

subroutine Barostat
    ! To apply barostat to get constant Pressure
    
    implicit none
    include 'dpdflow.h'

    real*8 mu,P
    
    real*8 alpha,cc,d
    
    integer n,k,i,j,l,rhoSize(3),rhoSizeH(3),ic(3)
    real*8 T,rho,rtemp,rhoV
    
    real*8, allocatable :: rhoGrid(:,:,:,:)
    
    alpha=0.1
    cc=4.16
    d=18.
        
	if(stepCount .eq. stepEquil) then
    open(501, file = './data/Barostat/TrhoP.plt')    
    write(501,'(''TITLE = "EOS"'')')
    write(501,'(''Variables= "time","T","rho","P","mu","tau"'')')
	endif
    !get the liquid temperature    
    T =0.0
    do n = nChainend + 1, nAtom
        do k =1,NDIM
            T= T + rv(n,k)**2
        enddo
    enddo
    T = T/(3.0*(nAtom - nChainend))
    
    !get the liquid density
    do i=1,3
    rhoSize(i) = int(region(i)/(2.*rCut))
	if(mod(rhoSize(i),2) .eq. 0) rhoSize(i) = rhoSize(i) + 1
	rhoSizeH(i) = (rhoSize(i) - 1)/2
    enddo
       
    allocate(rhoGrid(-rhoSizeH(1):rhoSizeH(1),-rhoSizeH(2):rhoSizeH(2),-rhoSizeH(3):rhoSizeH(3),1:3))
    
    do i= -rhoSizeH(1),rhoSizeH(1),1
    	do j= -rhoSizeH(2),rhoSizeH(2),1
    		do k= -rhoSizeH(3),rhoSizeH(3),1
	        do l =1,3
	            rhoGrid(i,j,k,l) =0.0
	        enddo
	        enddo
        enddo       
    enddo
    !pause
    do n = nWallAtom + 1, nAtom
        do i=1,3
       		rtemp = (abs(r(n,i))**2)/(r(n,i)*(2.*rCut))
       		rtemp = rtemp + (rtemp/(abs(rtemp)))*0.5
       		ic(i)=int(rtemp)                    
        enddo

        if (abs(ic(1)) .le. rhoSizeH(1) .and. abs(ic(2)) .le. rhoSizeH(2) .and. abs(ic(3)) .le. rhoSizeH(3)) then
	        rhoGrid(ic(1),ic(2),ic(3),1) =  rhoGrid(ic(1),ic(2),ic(3),1) + 1.
	    
	        if (n .le. nDpEnd) then
	            rhoGrid(ic(1),ic(2),ic(3),2) =  rhoGrid(ic(1),ic(2),ic(3),2) + 1.
	        elseif (n .le. nChainend) then
	            rhoGrid(ic(1),ic(2),ic(3),3) =  rhoGrid(ic(1),ic(2),ic(3),3) + 1.
	        endif 
        endif      
    enddo
    l=0
    rho= 0.
    do i= -rhoSizeH(1)+1,rhoSizeH(1)-1,1
    	do j= -rhoSizeH(2)+1,rhoSizeH(2)-1,1
    		do k= -rhoSizeH(3)+1,rhoSizeH(3)-1,1
	        if (rhoGrid(i,j,k,2) .eq. 0.0 .and. rhoGrid(i,j,k,3) .eq. 0.0 ) then
	            rho = rho + rhoGrid(i,j,k,1)/((2.*rCut)**3)
	            l=l+1
	        endif 
	        enddo
        enddo       
    enddo
    !pause
    rho = rho /l

    !! The Berendsen Barostat
    P = 2*alpha*alphaB*(rCut2**4)*(rho**3)+(alpha*alphaf-2*alpha*alphaB*(rCut2**4)*cc)*(rho**2)+T*rho+2*alpha*alphaB*(rCut2**4)*d
    
    if(stepCount .eq.  (JStep + stepEquil)) P0 = JP0
    mu=1
    if (mod(stepCount,BStep) .eq. 0 .and. (rho0 .ge.0 .or. stepCount .ge.JStep+stepEquil)) then
        !print*,deltaT,tau,rho,P,P0
        mu= 1-((deltaT/tau)*(P0-P))!)**(0.3333333333)

        do while (mu .le. 0.9)
            if ((P0-P) .gt. 0) then
                tau = tau *2
            else
                tau = tau /2
            endif
            print*, "the rise time (tau) is too small, change it to be: ", tau
            mu= 1-((deltaT/tau)*(P0-P))!)**(0.3333333333)
        enddo

        do while (mu .gt. 1.1)
            if ((P0-P) .gt. 0) then
                tau = tau /2
            else
                tau = tau *2
            endif
            print*, "the rise time (tau) is too great, change it to be: ", tau
            mu= 1-((deltaT/tau)*(P0-P))!)**(0.3333333333)
        enddo
        
        !print*,mu,rho,P
        
        !refresh the region
        do n = nChainend + 1, nAtom
            if(abs(r(n,3)) .gt. (regionH(3)*lpercent)) r(n,3) = r(n,3)*mu
        enddo
        

        region(3) = region(3)*mu
        regionH(3) = region(3)/2.

    	cells(3) = int(region(3) / rCut)
    	if(cells(3) .lt. 1) cells(3) = 1

        
        binvolm = region(1)*region(2)*region(3)/hsize
        
        maxList =  nAtom+cells(1)*cells(2)*cells(3)
    
    endif


    write(501,'(f6.2,2f10.4,f16.6,f10.4,f10.4)') timeNow,T,rho,P,mu,tau

    if(stepCount .ge. stepLimit) close(501)

    deallocate(rhoGrid)
    end

!! Added by linyuqing   


subroutine BubbleSize
    ! To get the bubble radius and centre point

    implicit none
    include 'dpdflow.h'
    
    integer i,j,k,l,m,n,Npb
    real*8 BCr(3),theta,Ixx,Iyy,Izz,Ixz,BRadius
    real*8 S1,S2,S3
    
    real*8, allocatable :: BPr(:,:)
    
    Npb = nDpEnd - nWallAtom
    
    allocate(BPr(Npb,3))
    
    do i =1,3
        BCr(i) = 0.
    enddo 
    
    Ixx = 0.
    Iyy = 0.
    Izz = 0.
    Ixz = 0.

    do n = nWallAtom+1, nDpEnd
        do i =1,3
        BCr(i) = BCr(i) + r(n,i) + ncc(n,i)*region(i) 
        BPr(n - nWallAtom,i) = r(n,i) + ncc(n,i)*region(i) 
        enddo
    enddo
    do i =1,3
        BCr(i) = BCr(i) / (nDpEnd - nWallAtom)
        do n =1,Npb
            BPr(n,i) = BPr(n,i) - BCr(i)
        enddo
    enddo   
    
    do n =1,Npb
        Ixx = Ixx + BPr(n,1)**2
        Iyy = Iyy + BPr(n,2)**2
        Izz = Izz + BPr(n,3)**2
        Ixz = Ixz + BPr(n,3)*BPr(n,1)
    enddo
       
    theta = 0.5*atan(2*Ixz/(Ixx-Izz))

    S1 = Ixx*cos(theta)**2 + Izz*sin(theta)**2 + 2*Ixz*sin(theta)*cos(theta)

    S2 = Ixx*sin(theta)**2 + Izz*cos(theta)**2 - 2*Ixz*sin(theta)*cos(theta)

    S3 = Iyy
    
    
    BRadius = (125*S1*S2*S3/(Npb**3))**(0.0666666666666)
    
    if(stepCount .eq. stepEquil) then
    open(601, file = './data/RP.plt')    
    write(601,'(''TITLE = "RP"'')')
    write(601,'(''Variables= "time","Radius","Cx","Cy","Cz"'')')
    endif
    
    !print*, BRadius,(BCr(i),i=1,3)
    
    write(601,'(f6.2,4f10.4)') timeNow,BRadius,(BCr(i),i=1,3)
    
    if(stepCount .ge. stepLimit)  close(601)




end
!! Added by linyuqing   