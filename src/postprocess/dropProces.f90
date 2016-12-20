    subroutine dropProces(oupt51)

    implicit none
    include 'dpdflow.h'

    character*50 char
    integer oupt51, oupt52
    integer DpNum,npt
    real*8 S1, S2, S3, LABC,BDnst,InterDnst
    real*8, allocatable :: DpX(:,:),DpY(:,:),DpZ(:,:)
    real*8, allocatable :: sumDpX(:),sumDpY(:),sumDpZ(:), &
                           sqDpX(:), sqDpY(:), sqDpZ(:),  &
                           sqDpXY(:),sqDpXZ(:),sqDpYZ(:), &
                           theta(:), LA(:), LB(:), LC(:), &
                           volume(:), deform(:)
    integer i,j,k,n,in,ChNum
    character*3 str

	data 	pi / 3.14159265358979 /
    real  rho
    
    common /drp/ BDnst,InterDnst
    rho = BDnst
    shearRate = 0.05

    read(lst5, *)
    read(lst5, *)
    read(lst5, *)
    read(lst5, '(1x,A25,i4)') char, nDp
    read(lst5, '(1x,A31,i4)') char, npt

    allocate(DpX(nDp,npt),DpY(nDp,npt),DpZ(nDp,npt))
    allocate(sumDpX(nDp), sumDpY(nDp), sumDpZ(nDp), &
              sqDpX(nDp),  sqDpY(nDp),  sqDpZ(nDp), &
              sqDpXY(nDp), sqDpXZ(nDp), sqDpYZ(nDp),&
              theta(nDp), LA(nDp), LB(nDp), LC(nDp),&
              volume(nDp), deform(nDp) )
       
    do i = 1, nDP
       oupt52 = 90 + i
       write(str, '(i3)') i
       open(unit = oupt52, file = './data/output/drop_deform'//str//'.dat', status = 'unknown')
    enddo

    in = int((stepLimit-stepEquil)/(stepChainProps*limitChainProps))+1

    do i = 1, in
       sumDpX(:) = 0
       sumDpY(:) = 0
       sumDpZ(:) = 0

       sqDpX(:)  = 0
       sqDpY(:)  = 0
       sqDpZ(:)  = 0

       sqDpXY(:) = 0
       sqDpXZ(:) = 0
       sqDpYZ(:) = 0

       do n = 1, nDp
          read(lst5, '(3x,A10,i4,2x,A7,f9.2)') char, ChNum, char, timeNow
          write(oupt51,'(2x,''zone'')') !, 2x, ''"t='', f9.2,2x,''DpN='', i4,''"'')') timeNow,ChNum 
          do j = 1, npt
             read(lst5,'(1x,3f10.5,3i4)')   DpX(ChNum,j),DpY(ChNum,j),DpZ(ChNum,j),ncc(j,1:NDIM)
             write(oupt51,'(1x,3f10.5)') DpX(ChNum,j),DpY(ChNum,j),DpZ(ChNum,j)


             DpX(ChNum,j) = DpX(ChNum,j) + ncc(j,1)*region(1) 
             DpY(ChNum,j) = DpY(ChNum,j) + ncc(j,2)*region(2)
             DpZ(ChNum,j) = DpZ(ChNum,j) + ncc(j,3)*initUcell(3)*gap(3)
             sumDpX(ChNum) = sumDpX(ChNum) + DpX(ChNum,j)
             sumDpY(ChNum) = sumDpY(ChNum) + DpY(ChNum,j)
             sumDpZ(ChNum) = sumDpZ(ChNum) + DpZ(ChNum,j)
!100         continue           
          enddo
!         npt=npt-1
          sumDpX(ChNum) = sumDpX(ChNum)/npt
          sumDpY(ChNum) = sumDpY(ChNum)/npt
          sumDpZ(ChNum) = sumDpZ(ChNum)/npt

          do j = 1, npt
             sqDpX(ChNum)  = sqDpX(ChNum)  + (DpX(ChNum,j)-sumDpX(ChNum))**2
             sqDpY(ChNum)  = sqDpY(ChNum)  + (DpY(ChNum,j)-sumDpY(ChNum))**2
             sqDpZ(ChNum)  = sqDpZ(ChNum)  + (DpZ(ChNum,j)-sumDpZ(ChNum))**2
!            sqDpXY(ChNum) = sqDpXY(ChNum) + (DpX(ChNum,j)-sumDpX(ChNum))*  &
!                                            (DpY(ChNum,j)-sumDpY(ChNum))
             sqDpXZ(ChNum) = sqDpXZ(ChNum) + (DpX(ChNum,j)-sumDpX(ChNum))*  &
                                             (DpZ(ChNum,j)-sumDpZ(ChNum))
!            sqDpYZ(ChNum) = sqDpYZ(ChNum) + (DpY(ChNum,j)-sumDpY(ChNum))*  &
!                                            (DpZ(ChNum,j)-sumDpZ(ChNum))
          enddo

          theta(ChNum) = 0.5*atan(2*sqDpXZ(ChNum)/(sqDpX(ChNum)-sqDpZ(ChNum)))
!         theta(ChNum) = 0.

          S1 = sqDpX(ChNum)*cos(theta(ChNum))**2 + &
               sqDpZ(ChNum)*sin(theta(ChNum))**2 + &
               2*sqDpXZ(ChNum)*sin(theta(ChNum))*cos(theta(ChNum))

          S2 = sqDpX(ChNum)*sin(theta(ChNum))**2 + &
               sqDpZ(ChNum)*cos(theta(ChNum))**2 - &
               2*sqDpXZ(ChNum)*sin(theta(ChNum))*cos(theta(ChNum))

          S3 = sqDpY(ChNum)

          LABC = (S1*S2*S3/(4*pi/15)**3/rho**3)**0.25
!          LABC = (S1*S2*S3/(4*pi/15)**3/rho**4)**0.25

          LA(ChNum) = (S1/rho/(4*pi/15)/LABC)**0.5
          LB(ChNum) = (S2/rho/(4*pi/15)/LABC)**0.5
          LC(ChNum) = (S3/rho/(4*pi/15)/LABC)**0.5

          volume(ChNum) = 4./3.*pi*LA(ChNum)*LB(ChNum)*LC(ChNum)

          deform(ChNum) = (LA(ChNum)-LB(ChNum))/(LA(ChNum)+LB(ChNum))   

       enddo

       do n = 1, nDP
          oupt52 = 90 + n
          write(oupt52, '(1x, f8.2, f9.4, 3f12.4, 3f9.4, 2f12.5 )') timeNow, theta(n)*180./pi, &
                sumDpX(n), sumDpY(n), sumDpZ(n), &
                LA(n), LB(n), LC(n), volume(n), deform(n)
       enddo
    enddo

    
    end