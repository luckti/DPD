!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

    subroutine makedrop(m , idom, RdsDp, iDp, subreg, seeds, fgap, subdom, rc, x, nx )

	implicit none

!	parameters

	integer	MAXCTP, MAXBDS, NDIM, MAXDP,nDp
	parameter (MAXCTP = 10, MAXBDS = 90000, NDIM = 3, MAXDP = 30)
	integer	 m, nChainend, seeds, iDp
	real*8	subdom(MAXCTP,8), subreg(MAXCTP,NDIM), x(MAXBDS,NDIM),        &
            xx(MAXBDS,NDIM), fgap, density, RdsDp, rc(MAXCTP, 5000, NDIM), RdsDp0, &
            Cx(MAXDP), Cy(MAXDP), Cz(MAXDP), xDp(MAXBDS,NDIM),pi,dratio
    integer nx(MAXBDS, NDIM), nxx(MAXBDS, NDIM)
    parameter (pi=3.1415926536)
!	Locals

	integer MAXATOM, MAXWATM, i, k, idom, mc, mp
	real*4 RandR
	real*8	d,ddensity, y(NDIM), shift(NDIM), e(NDIM)
    real*8  s1, s2, s3, s4, s5, s6, s7, s8, s9, dR
	parameter (MAXWATM = 50000, MAXATOM = 400000)

    logical existfile
    integer sph, npt, nc

    common  /dp/ Cx, Cy, Cz, density,dratio
    common /dpr/ nDp
   

    sph = 999

    inquire(file='./data/sphere/sphere.dat',exist=existfile)
    !if(.not.existfile)call system('mkdir .\data\sphere')
    write(*,*) 'start to generate the droplet molecules ...'
    
    if(dratio.gt.0) then
        ddensity=density/dratio
    elseif(dratio.eq.0.0) then
        npt=0
        goto 120
    else
        ddensity=0.0
    endif
    
    if (existfile)then
        
        open(sph, file='./data/sphere/sphere.dat', status='unknown')
        read(sph, '(1x, i5, f8.5)') npt, RdsDp0 

        if(npt.ne.anint((4./3.*pi*RdsDp**3)*ddensity).or.RdsDp0.ne.RdsDp)then           
            close(sph)
            call sphereDistribution(ddensity,RdsDp)
            open(sph, file='./data/sphere/sphere.dat', status='unknown')
            read(sph, '(1x, i5, f8.5)') npt, RdsDp0
        endif
    else
        call sphereDistribution(ddensity,RdsDp)
        open(sph, file='./data/sphere/sphere.dat', status='unknown')
        read(sph, '(1x, i5, f8.5)') npt, RdsDp0
    endif


    
    if(RdsDp .ne. RdsDp0) then
       write(*,*) 'wrong input droplet data file'
       stop
    endif
    

    do i = 1, npt
       read(sph, *) xDp(i,1), xDp(i,2), xDp(i,3)
    enddo

    close(sph)

120	d = 0.5*fgap
    
! find the center of a droplet
5 continue

  if (nDp .gt. 1)then
!   do k = 1, 200
    y(1) = subdom(idom,1) + 0.5*(1.0 + RandR(seeds))*subreg(idom,1)
    y(2) = 0.5*RandR(seeds)*subreg(idom,2) 
    y(3) = subdom(idom,2) + 0.5*(1.0 + RandR(seeds))*subreg(idom,3)
!   enddo
  else
      y(1)=Cx(1)
      y(2)=Cy(1)
      y(3)=Cz(1)
  endif
  
  
    dR = 0.7

    shift(1) = subreg(idom,1)
    shift(2) = subreg(idom,2)
    shift(3) = 0.

    do i = 1, idp - 1
        s1 = 0.
        do k = 1, NDIM
            s1 = s1 + (y(k) -  rc(idom, i, k))**2
        enddo
        
        s2 = (y(1) + shift(1) - rc(idom,i,1))**2 + &
             (y(2) - rc(idom,i,2))**2 +            &
             (y(3) - rc(idom,i,3))**2

        s3 = (y(1) - shift(1) - rc(idom,i,1))**2 + &
             (y(2) - rc(idom,i,2))**2 +            &
             (y(3) - rc(idom,i,3))**2

        s4 = (y(1) - rc(idom,i,1))**2 +            &
             (y(2) + shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

        s5 = (y(1) - rc(idom,i,1))**2 +            &
             (y(2) - shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

        s6 = (y(1) + shift(1) - rc(idom,i,1))**2 + &
             (y(2) - shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

        s7 = (y(1) + shift(1) - rc(idom,i,1))**2 + &
             (y(2) + shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

        s8 = (y(1) - shift(1) - rc(idom,i,1))**2 + &
             (y(2) - shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

        s9 = (y(1) - shift(1) - rc(idom,i,1))**2 + &
             (y(2) + shift(2) - rc(idom,i,2))**2 + &
             (y(3) - rc(idom,i,3))**2

! the distance between two droplets centers should be larger than 2*radius
        if( s1 .lt. (2*RdsDp + dR)**2 .or.  &
            s2 .lt. (2*RdsDp + dR)**2 .or.  &
            s3 .lt. (2*RdsDp + dR)**2 .or.  &
            s4 .lt. (2*RdsDp + dR)**2 .or.  &
            s5 .lt. (2*RdsDp + dR)**2 .or.  &
            s6 .lt. (2*RdsDp + dR)**2 .or.  &
            s7 .lt. (2*RdsDp + dR)**2 .or.  &
            s8 .lt. (2*RdsDp + dR)**2 .or.  &
            s9 .lt. (2*RdsDp + dR)**2 ) goto 5
    enddo

!    y(1) = subdom(idom,1) + 0.5*subreg(idom,1) + Cx(iDp)*RdsDp
!    y(2) = Cy(iDp)*RdsDp
!    y(3) = subdom(idom,2) + 0.5*subreg(idom,3) + Cz(iDp)*RdsDp

    do i = 1, npt
        do k = 1, NDIM
            xx(i,k) = xDp(i,k) + y(k)
            nxx(i,k) = 0
        enddo

        if(xx(i,1) .lt. -0.5*subreg(idom,1)) then
            xx(i,1) = xx(i,1) + subreg(idom,1)
            nxx(i,1) = -1
        endif
        if(xx(i,1) .ge. 0.5*subreg(idom,1)) then
            xx(i,1) = xx(i,1) - subreg(idom,1)
            nxx(i,1) = 1
        endif

        if(xx(i,2) .lt. -0.5*subreg(idom,2)) then
            xx(i,2) = xx(i,2) + subreg(idom,2)
            nxx(i,2) = -1
        endif
        if(xx(i,2) .ge. 0.5*subreg(idom,2)) then
            xx(i,2) = xx(i,2) - subreg(idom,2)
            nxx(i,2) = 1
        endif

!        if(xx(i,3) .lt. -0.5*subreg(idom,3)) then
!            xx(i,3) = xx(i,3) + subreg(idom,3)
!            nxx(i,3) = -1
!        endif
!        if(xx(i,3) .ge. 0.5*subreg(idom,3)) then
!            xx(i,3) = xx(i,3) - subreg(idom,3)
!            nxx(i,3) = 1
!        endif

!       if((xx(i,1) - subdom(idom,1)) .lt. 0.2 .or. &
!          (subdom(idom,3) - xx(i,1)) .lt. 0.2) then
!          goto 5
!        endif

       if((xx(i,3) - subdom(idom,2)) .lt. 0.2 .or. &
          (subdom(idom,8) - xx(i,3)) .lt. 0.2) then
          goto 5
        endif
    enddo

   write(*,'(2x, A9, i3, 2x, A10)') 'Droplet #', iDp, 'generated!'

    do i = 1, npt
        m = m + 1
        do k = 1, NDIM
            x(m,k) = xx(i,k)
            nx(m,k) = nxx(i,k)
        enddo
    enddo

    do k = 1, NDIM
        rc(idom, idp, k) = y(k)
    enddo
    
    write(* ,*) 'The coordinates of droplets centers: (x,y,z)'
    write(*,'( ''('',f8.4, '','', f8.4, '','', f8.4,'')'')') (y(k), i = 1, NDIM)
    
    
	return

    end
    
    
    
    
    !--------------------------------------------------------------------------
    subroutine sphereDistribution(rho,rdp)
    
    real*8 rho, rdp, volume, length
    integer np, nlayer

    real*8, allocatable :: rad(:), rnp(:), dis(:)
    real*8 :: x(1000,1000), y(1000,1000), z(1000,1000)
    real*8 :: rb(6000,3)

    integer i, ii, jj, n, m1, m2, ni
    parameter (pi=3.1415926536)

!    common / xyz / x, y, z
   ! write(*,*) 'please input "density" and "drop radius"'
   ! read(*,*) rho, rdp

    !rho = density
    !rdp = RdsDp

    volume = 4./3.*pi*rdp**3

    np = anint(volume*rho)

    print*, 'np =', np

    length = (1/rho)**(1./3.)

    nlayer = int(rdp/length)
 
! rad: radius of each layer
! sph: area of each sphere layer
! rnp: number of particle on each sphere
    allocate(dis(0:nlayer+1), rad(nlayer), rnp(0:nlayer+1))

    open(210, file='./data/sphere/sphere_tec.dat')
    open(211, file='./data/sphere/sphere.dat')
    open(200, file='./data/sphere/sphereOrin.dat')

    write(210,*) 'title=""'
	write(210,*) 'VARIABLES="x","y","z" '

    dis(0) = 0.
    rnp(0) = 0.
    n      = 0.

    do i = 1, nlayer
       dis(i) = rdp/nlayer*i
       
       rad(i) = 0.5*(dis(i)+dis(i-1)) 
       
       rnp(i) = anint(4./3.*pi*dis(i)**3/length**3) - n

       if(i .eq. nlayer) rnp(i) = np - n

       m1 = anint(sqrt(rnp(i)))
       m2 = anint(rnp(i)/m1)
!       m1 = 2
!       m2 = int(rnp(i)/m1)

       if ((rnp(i)-m1*m2) .gt. 0.3*max(m1,m2)) m2 = m2 + 1

       write(*,'(2x, f12.6, 2x, 3I8)'), rad(i), int(rnp(i)), m1, m2

       call distribtParticle(m1, m2, x, y, z)

       write(210,*) 'ZONE T="ZONE 01"'
       write(210,*) 'I=', m1,',J=', m2,',K=', 1,',F=POINT'

       do ii = 1, m1
          do jj = 1, m2
             n = n + 1
             rb(n,1) = x(ii,jj)*rad(i)
             rb(n,2) = y(ii,jj)*rad(i)
             rb(n,3) = z(ii,jj)*rad(i)

             write(210, '(2x, 3f12.6)') rb(n,1), rb(n,2), rb(n,3)           
          enddo
       enddo
    enddo

    if(n .lt. np) then
        print*, 'particle is not enough'
!        pause
    elseif(n .gt. np) then
        ni = n - np + 1
    endif

    print*, 'n=', n - ni + 1, ni - 1

    write(211, '(1x, i5, f8.5)') n - ni + 1, rdp
    do i = ni, n
       write(211, '(1x, 3f12.6)') rb(i,1), rb(i,2), rb(i,3)
       write(200, '(3f9.4)') rb(i,1), rb(i,2), rb(i,3)
    enddo

    close(210)
    close(211)
    close(200)
    deallocate (dis, rad, rnp)

    end

!---------------------------------------------------------------------

    subroutine distribtParticle(m, n, xx, yy, zz)

    implicit none

    integer::m,n

    real*8 :: xx(1000,1000), yy(1000,1000), zz(1000,1000)
    real*8, allocatable :: xt1(:,:), yt1(:,:), zt1(:,:)
    real*8, allocatable :: xt2(:,:), yt2(:,:), zt2(:,:), e(:,:,:)

    integer::i,j,k,ii,jj
    real*8::PI, ct

    real*8::alfa(m),beta(n),dalfa,dbeta
    
    real*8::F(m,n,3),Ftemp(3)

    parameter (pi=3.1415926536)
    
    parameter(ct=1e-5)

    real*8::error
    real*8::t1,t2
!    common / xyz / x, y, z

    allocate(xt1(m,n), yt1(m,n), zt1(m,n),  &
             xt2(m,n), yt2(m,n), zt2(m,n),  &
             e(m,n,3))

!    call getCPUtime(t1)

    dalfa=PI/25.0
    dbeta=PI/20.0
    
    do i=-m/2,-1
       alfa(i+m/2+1)=i*dalfa
    end do

    do i=1,m/2
       alfa(i+m/2)=i*dalfa
    end do

    do j=0,n-1
       beta(j+1)=j*dbeta
    end do

    do i=1,m
       do j=1,n
          call angletoCartesian(alfa(i),beta(j),xx(i,j),yy(i,j),zz(i,j))
       end do
    end do

    do while(.true.)
       F=0
       do i=1,m
          do j=1,n
             do ii=1,m
                do jj=1,n
                   if((ii/=i).or.(jj/=j)) then
                      call Ftwopoints(xx(ii,jj),yy(ii,jj),zz(ii,jj),  &
                                      xx(i,j),yy(i,j),zz(i,j),Ftemp)
	                  F(i,j,:)=F(i,j,:)+Ftemp(:)
	               end if
                end do
             end do

             xt1(i,j)=xx(i,j)+F(i,j,1)*ct
             yt1(i,j)=yy(i,j)+F(i,j,2)*ct
             zt1(i,j)=zz(i,j)+F(i,j,3)*ct

             call pointtosphere(xt1(i,j),yt1(i,j),zt1(i,j),   &
                                xt2(i,j),yt2(i,j),zt2(i,j))

             e(i,j,1)=abs(xx(i,j)-xt2(i,j))
             e(i,j,2)=abs(yy(i,j)-yt2(i,j))
             e(i,j,3)=abs(zz(i,j)-zt2(i,j))

	         xx(i,j)=xt2(i,j)
	         yy(i,j)=yt2(i,j)
	         zz(i,j)=zt2(i,j)
          end do
       end do

       error=maxval(e)
       if(error<1e-4) exit

!    write(*,*)"error=",error
    end do
!    call getCPUtime(t2)
!    write(*,*)t2-t1

!    open(unit=210, file='points360.dat')
!    do i=1,m
!       do j=1,n
!          write(210,*)x(i,j),y(i,j),z(i,j)
!       end do
!    end do

    deallocate(xt1, yt1, zt1,  &
               xt2, yt2, zt2,  &
               e)

    end
     
!--------------------------------------------------------------------

    subroutine angletoCartesian(angle1,angle2,x1,y1,z1)

    implicit none

    real*8,intent(in)::angle1,angle2
    real*8,intent(out)::x1,y1,z1
    real*8::rxy

    rxy=cos(angle1)
    z1=sin(angle1)
    x1=rxy*cos(angle2)
    y1=rxy*sin(angle2)
    
    end

!--------------------------------------------------------------------

    subroutine Ftwopoints(x1,y1,z1,x2,y2,z2,F2)

    implicit none

    real*8,intent(in)::x1,y1,z1,x2,y2,z2
    real*8,intent(out)::F2(3)
    real*8::dx,dy,dz,d,dd,f
    real*8::xx
    dx=x2-x1
    dy=y2-y1
    dz=z2-z1

    dd=dx**2+dy**2+dz**2
    d=sqrt(dd)

    f=1.0/dd
!    f = 100000/dd

    F2(1)=f*dx/d
    F2(2)=f*dy/d
    F2(3)=f*dz/d
    
    end

!--------------------------------------------------------------------

    subroutine pointtosphere(x0,y0,z0,x1,y1,z1)

    implicit none

    real*8,intent(in)::x0,y0,z0
    real*8,intent(out)::x1,y1,z1
    real*8::r

    r=sqrt(x0**2+y0**2+z0**2)
    x1=x0/r
    y1=y0/r
    z1=z0/r

    end

!--------------------------------------------------------------------

    subroutine getCPUtime(seconds)

    implicit none

!   Parameters

    integer	mclock
    real*4, intent(out):: seconds

!   Locals

    real*4	tarray(2), etime
    external etime

!   seconds= mclock()*0.01

    seconds= etime(tarray)
    
    return
    
    end

!--------------------------------------------------------------------