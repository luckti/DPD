    implicit none
    include 'dpdflow.h' 

    integer	oupt1, oupt2, oupt3, oupt4, oupt5, oupt6, oupt7,oupt,&
            oupt8, oupt9
    character*50 char
    integer i, j, k, n, m, in, nx, ny, nz,Vn,Fn,Bn,Intern,nlp
    real*8 Sampletime
    real*8 FDnst,BDnst,FPress,BPress,InterDnst,InterPress,&
            CBr(3),DBr(3),dspace
    integer ChNum
    !real*8 sliceX(3)
    
    common /drp/ BDnst,InterDnst

    integer,allocatable :: No(:,:,:),SumNBSur(:)
    real*8, allocatable :: Cx(:,:,:), Cy(:,:,:), Cz(:,:,:), Temp(:,:,:), &
                           Vx(:,:,:), Vy(:,:,:), Vz(:,:,:), Dens(:,:,:), &
                           Vof1(:,:,:), Vof2(:,:,:), Vof3(:,:,:),        &
                           Sxx(:,:,:),  Syy(:,:,:),  Szz(:,:,:),         &
                           Sxy(:,:,:),  Sxz(:,:,:),  Syz(:,:,:),         &
                           Press(:,:,:),N1(:,:,:),   N2(:,:,:)
    real*8, allocatable :: VVx(:,:), VVy(:,:), VVz(:,:),     &
                           Strxx(:), Stryy(:), Strzz(:),     &
                           Strxy(:), Strxz(:), Stryz(:),     &
                           Prssr(:), NN1(:),   NN2(:),       &
                           Tmpt(:,:),Dnst(:,:),              &
                           VlmF1(:,:), VlmF2(:,:), VlmF3(:,:)
    real*8, allocatable :: ChX(:,:),ChY(:,:),ChZ(:,:),BSurPress(:)
    
    character*3 str

    lst1 = 11
    lst2 = 12
    lst3 = 13
    lst4 = 14
    lst5 = 16
    
    oupt=2
!   oupt1= 21
!   oupt2= 22
!   oupt3= 23
    oupt4= 24
    oupt5= 25
    oupt6= 26
    oupt7= 27
    oupt8= 28
    oupt9= 29

    open (unit = lst1, file = './data/output/fluid1.dat' , status = 'unknown')
    open (unit = lst2, file = './data/output/fluid2.dat' , status = 'unknown')
    open (unit = lst3, file = './data/output/fluid3.dat' , status = 'unknown')
    open (unit = lst4, file = './data/output/fluid4.dat' , status = 'unknown')

!   call readparam0
    call readparam

!   open(unit = oupt1, file = './data/output/ProfileSumVel.dat', status = 'unknown')
!   open(unit = oupt2, file = './data/output/StressSum.dat',     status = 'unknown')
!   open(unit = oupt3, file = './data/output/ProfileSliVel.dat', status = 'unknown')
    
    if(nChainconf .gt. 0) open(unit = oupt4, file = './data/output/ChainCoordinate.plt', status = 'unknown')

    if(nDp        .gt. 0) then
        open(unit = lst5 , file = './data/output/droplet.dat' , status = 'unknown')
        open(unit = oupt5, file = './data/output/drop_tec.plt', status = 'unknown')
    endif
   
    open(unit = oupt6, file = './data/output/ContourVel.plt',    status = 'unknown')
    open(unit = oupt7, file = './data/output/ContourStress.plt', status = 'unknown')
    open(unit = oupt, file = './data/output/DensityPress.dat',     status = 'unknown')              !output the density and pressure for each phase
    open(unit = oupt+1, file = './data/output/BubbleSPDistr.plt',     status = 'unknown')
    
    !caulate the surface tension
    open(unit = oupt8, file = './data/output/StrCell.plt',    status = 'unknown')
    open(unit = oupt9, file = './data/output/PlaneStr.plt', status = 'unknown')

    in = int(stepEquil/stepAvg)-1
    
    do i=1,in  
        read(lst1,'(1x, i8, 2f9.4, f10.4, 3f8.4, f9.4, f8.4)') stepCount, timeNow, &
            vSum, sTotEnergy, ssTotEnergy, sKinEnergy, ssKinEnergy, &
            sPressure, ssPressure
    enddo

!   write(*, '("Please input the x-coordinates of 3 slices")')
!   read(*, *) sliceX(1:3)

!    stepLimit = 60000
    in = int((stepLimit-stepEquil-startSample)/stepSample)

    do i = 1, in
        oupt1 = 300 + i
        oupt2 = 400 + i
        write(str, '(i3)') i
        open(unit = oupt1, file = './data/output/ProfileSumVel'//str//'.dat', status = 'unknown')
        open(unit = oupt2, file = './data/output/StressSum'//str//'.dat',     status = 'unknown')
    enddo

    nx = sizeHistGrid(1)
    ny = sizeHistGrid(2)
    nz = sizeHistGrid(3)
    
    nlp=int(max(sizeHistGrid(1),sizeHistGrid(2),sizeHistGrid(3))/2.)
    dspace=min(region(1),region(2),region(3))/(nlp*2.)
    print *, nlp,dspace
    allocate(SumNBSur(1:nlp),BSurPress(1:nlp))

    allocate(No(nx,ny,nz),Cx(nx,ny,nz),Cy(nx,ny,nz),Cz(nx,ny,nz),Temp(nx,ny,nz),   &
                Vx(nx,ny,nz),Vy(nx,ny,nz),Vz(nx,ny,nz),Dens(nx,ny,nz),      &
                Vof1(nx,ny,nz),Vof2(nx,ny,nz),Vof3(nx,ny,nz),  &
                Sxx(nx,ny,nz),Syy(nx,ny,nz),Szz(nx,ny,nz),     &
                Sxy(nx,ny,nz),Sxz(nx,ny,nz),Syz(nx,ny,nz),     &
                Press(nx,ny,nz),N1(nx,ny,nz),N2(nx,ny,nz))
    allocate(Strxx(nz),Stryy(nz),Strzz(nz),       &
             Strxy(nz),Strxz(nz),Stryz(nz),       &
             Prssr(nz),NN1(nz),NN2(nz),           &
             Tmpt(nz,4),Dnst(nz,4),               &
             VVx(nz,4),VVy(nz,4),VVz(nz,4),       &
             VlmF1(nz,4),VlmF2(nz,4),VlmF3(nz,4))
    allocate(ChX(nChain,ChainLen),ChY(nChain,ChainLen),ChZ(nChain,ChainLen))   

    do n=1,in
        read(lst1,*)
        read(lst1,*)
        read(lst1,'(10x,A7,f9.4)') char, Sampletime
        read(lst1,*)

        read(lst2,*)
        read(lst2,*)
        read(lst2,'(10x,A7,f9.4)') char, Sampletime
        read(lst2,*)

        Strxx(:) = 0.
        Stryy(:) = 0.
        Strzz(:) = 0.
        Strxy(:) = 0.
        Strxz(:) = 0.
        Stryz(:) = 0.
        Prssr(:) = 0.
        NN1(:)   = 0.
        NN2(:)   = 0.

        Tmpt(:,:) = 0.
        Dnst(:,:) = 0.
        VVx (:,:) = 0.
        VVy (:,:) = 0.
        VVz (:,:) = 0.
        VlmF1(:,:)= 0.
        VlmF2(:,:)= 0.
        
        FDnst=0.
        BDnst=0.
        FPress=0.
        BPress=0.
        InterDnst=0.
        InterPress=0.
        Vn=0
        Fn=0
        Bn=0
        Intern=0
        CBr(1)=0.
        CBr(2)=0.
        CBr(3)=0.
        
        oupt1 = 300 + n
        oupt2 = 400 + n

        do i=1,sizeHistGrid(1)
            do j=1,sizeHistGrid(2)
                do k=1,sizeHistGrid(3) 
!                   read(lst1,'(1x,i7,3f9.4,f8.4,3e14.6,4f10.5)') No(i,j,k), &
!                              cx(i,j,k),cy(i,j,k),cz(i,j,k),Temp(i,j,k), &
!                              Vx(i,j,k),Vy(i,j,k),Vz(i,j,k),             &
!                              Vof1(i,j,k),Vof2(i,j,k),Vof3(i,j,k),dens(i,j,k)
                    read(lst1,'(1x,i7,3f9.4,f10.4,3e14.6,2f10.5)') No(i,j,k), &
                                cx(i,j,k),cy(i,j,k),cz(i,j,k),Temp(i,j,k), &
                                Vx(i,j,k),Vy(i,j,k),Vz(i,j,k),             &
                                Vof1(i,j,k),dens(i,j,k)

                    read(lst2,'(1x,i7,3f9.4,7e12.4)') No(i,j,k),   &
                                cx(i,j,k),cy(i,j,k),cz(i,j,k), &
                             Sxx(i,j,k),Syy(i,j,k),Szz(i,j,k), &
                             Sxy(i,j,k),Sxz(i,j,k),Syz(i,j,k), &
                             Press(i,j,k)
            
                    N1(i,j,k) = Sxx(i,j,k) - Szz(i,j,k)
                    N2(i,j,k) = Szz(i,j,k) - Syy(i,j,k)
            
                    Strxx(k)  = Strxx(k)  + Sxx(i,j,k)
                    Stryy(k)  = Stryy(k)  + Syy(i,j,k)
                    Strzz(k)  = Strzz(k)  + Szz(i,j,k)
                    Strxy(k)  = Strxy(k)  + Sxy(i,j,k)
                    Strxz(k)  = Strxz(k)  + Sxz(i,j,k)
                    Stryz(k)  = Stryz(k)  + Syz(i,j,k)
                    Prssr(k)  = Prssr(k)  + Press(i,j,k)
                    Tmpt(k,1) = Tmpt(k,1) + Temp(i,j,k)
                    Dnst(k,1) = Dnst(k,1) + dens(i,j,k)
                    VVx (k,1) = VVx (k,1) + Vx(i,j,k)
                    VVy (k,1) = VVy (k,1) + Vy(i,j,k)
                    VVz (k,1) = VVz (k,1) + Vz(i,j,k)
                    VlmF1(k,1)= VlmF1(k,1)+ Vof1(i,j,k)
                    VlmF2(k,1)= VlmF2(k,1)+ Vof2(i,j,k)
                    VlmF3(k,1)= VlmF3(k,1)+ Vof3(i,j,k)

! To sum the density and pressure for each phase

                    if(dens(i,j,k)==0)then
                        Vn=Vn+1
                    elseif (Vof1(i,j,k)==1.0) then
                        BDnst = BDnst + dens(i,j,k)
                        BPress = BPress + Press(i,j,k)
                        CBr(1)=CBr(1)+ cx(i,j,k)
                        CBr(2)=CBr(2)+ cy(i,j,k)
                        CBr(3)=CBr(3)+ cz(i,j,k)                        
                        Bn=Bn+1
                    elseif(Vof1(i,j,k)==0.0) then
                        FDnst = FDnst + dens(i,j,k)
                        FPress = FPress + Press(i,j,k)
                        Fn=Fn+1   
                    else
                        InterDnst = InterDnst + dens(i,j,k)
                        InterPress = InterPress + Press(i,j,k)
                        CBr(1)=CBr(1)+ cx(i,j,k)
                        CBr(2)=CBr(2)+ cy(i,j,k)
                        CBr(3)=CBr(3)+ cz(i,j,k)
                        Intern=Intern+1 
                    endif
! To sum the density and pressure for each phase, finished!
                    
 !                  do m = 1, 3
 !                      if(cx(i-1,j,k) .le. sliceX(m) .and. cx(i,j,k) .ge. sliceX(m)) then
 !                          Tmpt(k,m+1) = Tmpt(k,m+1) + Temp(i,j,k)
 !                          Dnst(k,m+1) = Dnst(k,m+1) + dens(i,j,k)
 !                          VVx (k,m+1) = VVx (k,m+1) + Vx(i,j,k)
 !                          VVy (k,m+1) = VVy (k,m+1) + Vy(i,j,k)
 !                          VVz (k,m+1) = VVz (k,m+1) + Vz(i,j,k)
 !                      endif
 !                  enddo
                enddo
            enddo
        enddo

        do k=1,sizeHistGrid(3)
            Strxx(k)  = Strxx(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Stryy(k)  = Stryy(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Strzz(k)  = Strzz(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Strxy(k)  = Strxy(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Strxz(k)  = Strxz(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Stryz(k)  = Stryz(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Prssr(k)  = Prssr(k)  /(sizeHistGrid(1)*sizeHistGrid(2))
            Tmpt(k,1) = Tmpt(k,1) /(sizeHistGrid(1)*sizeHistGrid(2))
            Dnst(k,1) = Dnst(k,1) /(sizeHistGrid(1)*sizeHistGrid(2))
            VVx(k,1)  = VVx(k,1)  /(sizeHistGrid(1)*sizeHistGrid(2))
            VVy(k,1)  = VVy(k,1)  /(sizeHistGrid(1)*sizeHistGrid(2))
            VVz(k,1)  = VVz(k,1)  /(sizeHistGrid(1)*sizeHistGrid(2))
            VlmF1(k,1)= VlmF1(k,1)/(sizeHistGrid(1)*sizeHistGrid(2))
            VlmF2(k,1)= VlmF2(k,1)/(sizeHistGrid(1)*sizeHistGrid(2))
            VlmF3(k,1)= VlmF3(k,1)/(sizeHistGrid(1)*sizeHistGrid(2))
!           do m = 1, 3
!               Tmpt(k,m+1) = Tmpt(k,m+1) /sizeHistGrid(2)
!               Dnst(k,m+1) = Dnst(k,m+1) /sizeHistGrid(2)
!               VVx(k,m+1)  = VVx(k,m+1)  /sizeHistGrid(2)
!               VVy(k,m+1)  = VVy(k,m+1)  /sizeHistGrid(2)
!               VVz(k,m+1)  = VVz(k,m+1)  /sizeHistGrid(2)
!           enddo
          
            NN1(k)    = Strxx(k) - Strzz(k)
            NN2(k)    = Strzz(k) - Stryy(k)
        enddo
        
        write(oupt1,101)
        write(oupt2,102)
!       write(oupt3,103)
        do k=1,sizeHistGrid(3)
            write(oupt1,'(f9.4,3e14.6,f8.4,4f10.5)') Cz(1,1,k), VVx(k,1),VVy(k,1),VVz(k,1),   &
                                                 Tmpt(k,1), Dnst(k,1), VlmF1(k,1), VlmF2(k,1), VlmF3(k,1)
            write(oupt2,'(f9.4,9e12.4)') Cz(1,1,k),Strxx(k),Stryy(k),Strzz(k),   &
                       Strxy(k),Strxz(k),Stryz(k),Prssr(k),NN1(k),NN2(k)
!           write(oupt3,'(f9.4,9e14.6)') Cz(1,1,k), VVx(k,2:4), VVy(k,2:4), VVz(k,2:4)
        enddo
        close(oupt1)
        close(oupt2) 
        
! To get the density and pressure for each phase        
        if(Bn>0)then
            BDnst = BDnst / Bn
            BPress = BPress / Bn            
        endif
        if(Fn>0)then
            FDnst = FDnst / Fn
            FPress = FPress / Fn
        endif
        if(Intern>0) then
            InterDnst = InterDnst / Intern
            InterPress = InterPress / Intern
        endif
        if(n==1) write(oupt,104)
        write(oupt,'(i4,3f10.5,3e12.4,4i8)') n,BDnst,FDnst,InterDnst,BPress,FPress,InterPress,Vn,Bn,Fn,Intern
        
        
! To get the pressure surounding the bubble  
        if((Vn+Bn+Intern) .gt. 0) then
            do m=1,3
                CBr(m)=CBr(m) / (Vn+Bn+Intern)
            enddo
        else
            print*,'we cannot get the centre of bubble!set to origin'
            do m=1,3
                CBr(m)=0.
            enddo
        endif
        write (*,'(i6,3f6.3)') n,(CBr(m),m=1,3)
        SumNBSur(1:nlp)= 0d0
        BSurPress(1:nlp)= 0.d0 
        do i=1,sizeHistGrid(1)
            do j=1,sizeHistGrid(2)
                do k=1,sizeHistGrid(3)
                    DBr(1)=CBr(1)-cx(i,j,k)
                    DBr(2)=CBr(2)-cy(i,j,k)
                    DBr(3)=CBr(3)-cz(i,j,k)
                    do m=1,3
                        if(DBr(m).gt.region(m)/2.) DBr(m)=DBr(m)-region(m)
                        if(DBr(m).lt.-region(m)/2.) DBr(m)=DBr(m)+region(m)
                    enddo 
                    if ((DBr(1)**2+DBr(2)**2+DBr(3)**2) .le. (nlp*dspace)**2) then
                        m=1
                        do while((DBr(1)**2+DBr(2)**2+DBr(3)**2) .gt. (m*dspace)**2)
                                m=m+1
                        enddo
                        if(m>0) then
                            SumNBSur(m)=SumNBSur(m) + 1
                            BSurPress(m)=BSurPress(m) + Press(i,j,k)  
                        endif
                    endif
                enddo
            enddo
        enddo

        if(n==1) write(oupt+1,'(''Variables= "X","Press","N"'')')
        write(oupt+1,'(''ZONE'')')
        write(oupt+1,'(''T="('',f6.3,'','',f6.3,'','',f6.3,'')"'')') (CBr(m),m=1,3)
        write(oupt+1,'(''I= '',i8,'',DATAPACKING = POINT'')') nlp
        do m =1 ,nlp
            if(SumNBSur(m) .gt. 0d0) BSurPress(m) = BSurPress(m)/SumNBSur(m)
            write(oupt+1,'(f10.4,e14.4,i)') (m-0.5)*dspace,BSurPress(m),SumNBSur(m)
        enddo 
                
    enddo
    close(oupt)
    close(oupt+1)

! Finish to ouput the pressure surounding the bubble       

    !enddo
    
    if(nChainconf .gt. 0) then
        read(lst3,*)        
        in = int((stepEquil-1000)/(stepChainProps*limitChainProps))
        if(in .gt. 0) then
            read(lst3,*)
            do i=1,in 
                read(lst3,'(i6,5f12.4)') stepCount,aaDistSq,eeDistSq,radGyrSq,gMomRatio1,gMomRatio2
            enddo
        endif

       in = int((stepLimit-stepEquil)/(stepChainProps*limitChainProps))+1
       do i=1,in
          do n=1,nChain,nChainConf
             read(lst3,'(3x,A11,i4,2x,A7,f9.2)') char,ChNum,char,timeNow
             write(oupt4,'(2x,''zone'' , 2x, ''t="'', f9.2,2x,''ChN='', i4,''"'')') timeNow,ChNum 
             do j=1,ChainLen
                read(lst3,'(1x,3f10.5)') ChX(ChNum,j),ChY(ChNum,j),ChZ(ChNum,j)
                write(oupt4,'(1x,3f10.5)') ChX(ChNum,j),ChY(ChNum,j),ChZ(ChNum,j)
             enddo
          enddo
       enddo
    endif
    
    if(nDp .gt. 0 ) then
       call dropProces(oupt5)
    endif

!    write(oupt1,101)
!    write(oupt2,102)
!!   write(oupt3,103)
!
!    do k=1,sizeHistGrid(3)
!       write(oupt1,'(f9.4,3e14.6,f8.4,4f10.5)') Cz(1,1,k), VVx(k,1),VVy(k,1),VVz(k,1),   &
!                                                 Tmpt(k,1), Dnst(k,1), VlmF1(k,1), VlmF2(k,1), VlmF3(k,1)
!       write(oupt2,'(f9.4,9e12.4)') Cz(1,1,k),Strxx(k),Stryy(k),Strzz(k),   &
!                       Strxy(k),Strxz(k),Stryz(k),Prssr(k),NN1(k),NN2(k)
!!      write(oupt3,'(f9.4,9e14.6)') Cz(1,1,k), VVx(k,2:4), VVy(k,2:4), VVz(k,2:4)
!    enddo

    write(oupt6, *) 'Title=" "'
    write(oupt6, *) 'Variables= "X","Y","Z","U","V","W","T","Rho","VofDP","VofH","VofT"'
    write(oupt6, *) 'zone'
    write(oupt6, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', nx, ',J=', ny, ',K=', nz, ',F=POINT'

    write(oupt7, *) 'Title=" "'
    write(oupt7, *) 'Variables= "X","Y","Z","Sxz","press","N1","N2"'
    write(oupt7, *) 'zone'
    write(oupt7, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', nx, ',J=', ny, ',K=', nz, ',F=POINT'
    
    write(oupt8, *) 'Title=" "'
    write(oupt8, *) 'Variables= "X","Y","Z","Sxx","Syy","Szz","Sxy","Sxz","Syz","Press"'
    write(oupt8, *) 'zone'
    write(oupt8, '(A3, i4, 1x, A4, i4, 1x, A4, i4, A8)') 'I=', nx, ',J=', ny, ',K=', nz, ',F=POINT'

    write(oupt9, *) 'Title=" "'
    write(oupt9, *) 'Variables= "Z","Sxx","Syy","Szz","Sxy","Sxz","Syz","Press","STStr"'
    write(oupt9, *) 'zone'
    write(oupt9, '(A8)') 'F=POINT'

    do k=1,sizeHistGrid(3)
       do j=1,sizeHistGrid(2)
          do i=1,sizeHistGrid(1)
             write(oupt6, '(3f9.4, 3e14.6, f10.4, 4f10.5)')   &
                             cx(i,j,k), cy(i,j,k), cz(i,j,k), &
                             Vx(i,j,k), Vy(i,j,k), Vz(i,j,k), &
                         Temp(i,j,k), Dens(i,j,k), Vof1(i,j,k), Vof2(i,j,k), Vof3(i,j,k)
             write(oupt7, '(3f9.4, 4e14.4)') cx(i,j,k), cy(i,j,k), cz(i,j,k), &
                                             Sxz(i,j,k),Press(i,j,k),N1(i,j,k), N2(i,j,k)
             write(oupt8, '(3f9.4,7e12.4)') cx(i,j,k),cy(i,j,k),cz(i,j,k), &
                                                 Sxx(i,j,k),Syy(i,j,k),Szz(i,j,k), &
                                                 Sxy(i,j,k),Sxz(i,j,k),Syz(i,j,k), &
                                                 Press(i,j,k)
               
            enddo
        enddo
        
         write(oupt9, '(f9.4, 8e12.4)') Cz(1,1,k),Strxx(k),Stryy(k),Strzz(k),   &
            Strxy(k),Strxz(k),Stryz(k),Prssr(k),(Strxx(k)-(Stryy(k)+Strzz(k))/2)   
    enddo
    
    deallocate (SumNBSur, BSurPress)
    deallocate (No,  Cx,  Cy,  Cz,  Temp, Vx,  Vy,    Vz, Vof1, Vof2, Vof3, Dens,  &
                Sxx, Syy, Szz, Sxy, Sxz,  Syz, Press, N1, N2)
    deallocate (Strxx, Stryy, Strzz, Strxy, Strxz, Stryz,   &
                Prssr, NN1,   NN2,   Tmpt,  Dnst,  VlmF1,   &
                VlmF2, VlmF3, VVx,   VVy,   VVz)
    deallocate (ChX, ChY, ChZ) 
    

    close(lst1)
    close(lst2)
    close(lst3)
    close(lst4)



!    close(oupt1)
!    close(oupt2)
!   close(oupt3)

    if(nChainconf .gt. 0) close(oupt4)
    if(nDp .gt. 0) close(oupt5)
    close(oupt6)
    close(oupt7)
    close(oupt8)
    close(oupt9)

101	format(5x, 'Z', 7x, 'Vx', 12x, 'Vy', 12x, 'Vz', 12x, 'T', &
           8x, 'Dens', 6x, 'VoFDP', 5x, 'VoFH', 6x, 'VoFT')
102	format(5x, 'Z', 7x, 'Sxx', 9x, 'Syy', 9x, 'Szz', 9x, 'Sxy', &
           9x, 'Sxz', 9x, 'Syz',9x,'Press',7x,'N1', 10x, 'N2')
103 format(5x, 'Z', 7x, 'Vx1', 11x, 'Vx2', 11x, 'Vx3',  &
                   11x, 'Vy1', 11x, 'Vy2', 11x, 'Vy3',  &
                   11x, 'Vz1', 11x, 'Vz2', 11x, 'Vz3')
104 format(2x, 'in', 1x, 'bubbleRho', 2x, 'fluidRho', 4x, 'MixRho',  &
                   1x, 'bubblePress', 2x, 'fluidPress', 4x, 'MixPress',&
                    2x,'VoidN',2x,'bubbleN',2x,'fluidN',4x,'MixN')

end program

