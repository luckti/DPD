!----------------------------------------------------------

subroutine readparam

implicit none
include 'dpdflow.h' 

character*50 char

read(lst1,*)
read(lst1,'(1x,A16, i7  )') char,runId
read(lst1,'(1x,A16,3i7  )') char,initUcell
read(lst1,'(1x,A16,4f7.2)') char,alphaf, alphaFD, alphaDD, alphaB
read(lst1,'(1x,A16,2f7.2)') char, alphaDA, alphaDB
read(lst1,'(1x,A16, f7.2)') char,alphaw
read(lst1,'(1x,A16,2f7.2)') char,rCut,rCut2
read(lst1,'(1x,A16,2f7.2)') char,rCutDp,rCutDp2
read(lst1,'(1x,A16,3f7.2)') char,gammaF,gammaD,gammaFD
read(lst1,'(1x,A16, f7.2)') char,gammaw
read(lst1,'(1x,A16, f7.2)') char,lambda
read(lst1,'(1x,A16, f7.2)') char,density
read(lst1,'(1x,A16, f7.2)') char,temperature
read(lst1,'(1x,A16, f7.2)') char,gravField
read(lst1,'(1x,A16, f7.2)') char,WLadjust
read(lst1,'(1x,A16, i7  )') char,nDp
read(lst1,'(1x,A16, f7.2)') char,RdsDp
read(lst1,'(1x,A16, i7  )') char,nChain
read(lst1,'(1x,A16, i7  )') char,ChainLen
read(lst1,'(1x,A16,3f7.2)') char,Hfene,rmaxfene,reqfene
read(lst1,'(1x,A16, f7.4)') char,deltaT
read(lst1,'(1x,A16, i7  )') char,stepAvg
read(lst1,'(1x,A16, i7  )') char,stepEquil
read(lst1,'(1x,A16, i7  )') char,startSample
read(lst1,'(1x,A16, i7  )') char,stepSample
read(lst1,'(1x,A16, i7  )') char,stepLimit
read(lst1,'(1x,A16,3i7  )') char,sizeHistGrid
read(lst1,'(1x,A16, i7  )') char,stepGrid
read(lst1,'(1x,A16, i7  )') char,limitGrid
read(lst1,'(1x,A16, i7  )') char,stepChainProps
read(lst1,'(1x,A16, i7  )') char,limitChainProps
read(lst1,'(1x,A16, i7  )') char,nChainConf
read(lst1,'(1x,A16, f7.2)') char,timeSteady

read(lst1,'(//1x, A24)')    char
read(lst1,'(1x,A19,3f9.4)') char,region
read(lst1,'(1x,A19,3(3x,i5  ))')  char,cells
read(lst1,'(1x,A19,3(3x,f7.4))')  char,gap
read(lst1,'(1x,A19, 5x, i8)') char,maxList
read(lst1,'(1x,A19, 5x, i8)') char,nFreeAtom
read(lst1,'(1x,A19, 5x, i8)') char,nWallAtom
read(lst1,'(1x,A19, 5x, i8)') char,nAtom
read(lst1,*)
read(lst1,*)
read(lst1,*)

read(lst2,*)
read(lst2,'(1x,A16, i7  )') char,runId
read(lst2,'(1x,A16,3i7  )') char,initUcell
read(lst2,'(1x,A16,4f7.2)') char,alphaf, alphaFD, alphaDD, alphaB
read(lst2,'(1x,A16, f7.2)') char,alphaw
read(lst2,'(1x,A16,2f7.2)') char,rCut,rCut2
read(lst2,'(1x,A16,3f7.2)') char,gammaF,gammaD,gammaFD
read(lst2,'(1x,A16, f7.2)') char,gammaw
read(lst2,'(1x,A16, f7.2)') char,lambda
read(lst2,'(1x,A16, f7.2)') char,density
read(lst2,'(1x,A16, f7.2)') char,temperature
read(lst2,'(1x,A16, f7.2)') char,gravField
read(lst2,'(1x,A16, f7.2)') char,WLadjust
read(lst2,'(1x,A16, i7  )') char,nDp
read(lst2,'(1x,A16, f7.2)') char,RdsDp
read(lst2,'(1x,A16, i7  )') char,nChain
read(lst2,'(1x,A16, i7  )') char,ChainLen
read(lst2,'(1x,A16,3f7.2)') char,Hfene, rmaxfene, reqfene
read(lst2,'(1x,A16, f7.4)') char,deltaT
read(lst2,'(1x,A16, i7  )') char,stepAvg
read(lst2,'(1x,A16, i7  )') char,stepEquil
read(lst2,'(1x,A16, i7  )') char,startSample
read(lst2,'(1x,A16, i7  )') char,stepSample
read(lst2,'(1x,A16, i7  )') char,stepLimit
read(lst2,'(1x,A16,3i7  )') char,sizeHistGrid
read(lst2,'(1x,A16, i7  )') char,stepGrid
read(lst2,'(1x,A16, i7  )') char,limitGrid
read(lst2,'(1x,A16, i7  )') char,stepChainProps
read(lst2,'(1x,A16, i7  )') char,limitChainProps
read(lst2,'(1x,A16, i7  )') char,nChainConf
read(lst2,'(1x,A16, f7.2)') char,timeSteady

read(lst3,'(//, A24)') char
read(lst3,'(1x,A19,3f9.4)') char,region
read(lst3,'(1x,A19,3(3x,i5))') char,cells
read(lst3,'(1x,A19,3(3x,f7.4))') char,gap
read(lst3,'(1x,A19,5x,i8)') char,maxlist
read(lst3,'(1x,A19,5x,i8)') char,nFreeAtom
read(lst3,'(1x,A19,5x,i8)') char,nWallAtom
read(lst3,'(1x,A19,5x,i8)') char,nAtom
read(lst3,'(/3x,A27)') char

end

!!----------------------------------------------------------
!
!subroutine readparam0
!
!implicit none
!include 'dpdflow.h' 
!
!character*50 char
!
!read(lst1,*)
!read(lst1,'(1x,A16,i3)')    char,runId
!read(lst1,'(1x,A16,3i5)')   char,initUcell
!read(lst1,'(1x,A16,f7.2)')  char,alphaf
!read(lst1,'(1x,A16,f7.2)')  char,alphaw
!read(lst1,'(1x,A16,f7.4)')  char,rCut
!read(lst1,'(1x,A16,f6.2)')  char,cigama
!read(lst1,'(1x,A16,f6.2)')  char,cigamaw
!read(lst1,'(1x,A16,f6.2)')  char,lambda
!read(lst1,'(1x,A16,f7.4)')  char,density
!read(lst1,'(1x,A16,f5.2)')  char,temperature
!read(lst1,'(1x,A16,f5.2)')  char,gravField
!read(lst1,'(1x,A16,f7.4)')  char,WLadjust
!read(lst1,'(1x,A16,i4)')    char,nChain
!read(lst1,'(1x,A16,i4)')    char,ChainLen
!read(lst1,'(1x,A16,3f8.2)') char,Hfene, rmaxfene, reqfene
!read(lst1,'(1x,A16,f7.4)')  char,deltaT
!read(lst1,'(1x,A16,i5)')    char,stepAvg
!read(lst1,'(1x,A16,i5)')    char,stepEquil
!read(lst1,'(1x,A16,i7)')    char,startSample
!read(lst1,'(1x,A16,i7)')    char,stepSample
!read(lst1,'(1x,A16,i7)')    char,stepLimit
!read(lst1,'(1x,A16,3i7)')   char,sizeHistGrid
!read(lst1,'(1x,A16,i7)')    char,stepGrid
!read(lst1,'(1x,A16,i7)')    char,limitGrid
!read(lst1,'(1x,A16,i7)')    char,stepChainProps
!read(lst1,'(1x,A16,i7)')    char,limitChainProps
!read(lst1,'(1x,A16,i3)')    char,nChainConf
!read(lst1,'(1x,A16,f9.2)')  char,timeSteady
!read(lst1,'(//1x, A24)')    char
!read(lst1,'(1x,A19,3f9.4)') char,region
!read(lst1,'(1x,A19,3(3x,i5))')  char,cells
!read(lst1,'(1x,A19,3(3x,f7.4))')    char,gap
!read(lst1,'(1x,A19,5x,i8)') char,maxList
!read(lst1,'(1x,A19,5x,i8)') char,nFreeAtom
!read(lst1,'(1x,A19,5x,i8)') char,nWallAtom
!read(lst1,'(1x,A19,5x,i8)') char,nAtom
!read(lst1,*)
!read(lst1,*)
!read(lst1,*)
!
!read(lst2,*)
!read(lst2,'(1x,A16,i3)')    char,runId
!read(lst2,'(1x,A16,3i5)')   char,initUcell
!read(lst2,'(1x,A16,f7.2)')  char,alphaf
!read(lst2,'(1x,A16,f7.2)')  char,alphaw
!read(lst2,'(1x,A16,f7.4)')  char,rCut
!read(lst2,'(1x,A16,f6.2)')  char,cigama
!read(lst2,'(1x,A16,f6.2)')  char,cigamaw
!read(lst2,'(1x,A16,f6.2)')  char,lambda
!read(lst2,'(1x,A16,f7.4)')  char,density
!read(lst2,'(1x,A16,f5.2)')  char,temperature
!read(lst2,'(1x,A16,f5.2)')  char,gravField
!read(lst2,'(1x,A16,f7.4)')  char,WLadjust
!read(lst2,'(1x,A16,i4)')    char,nChain
!read(lst2,'(1x,A16,i4)')    char,ChainLen
!read(lst2,'(1x,A16,3f8.2)') char,Hfene, rmaxfene, reqfene
!read(lst2,'(1x,A16,f7.4)')  char,deltaT
!read(lst2,'(1x,A16,i5)')    char,stepAvg
!read(lst2,'(1x,A16,i5)')    char,stepEquil
!read(lst2,'(1x,A16,i7)')    char,startSample
!read(lst2,'(1x,A16,i7)')    char,stepSample
!read(lst2,'(1x,A16,i7)')    char,stepLimit
!read(lst2,'(1x,A16,3i7)')   char,sizeHistGrid
!read(lst2,'(1x,A16,i7)')    char,stepGrid
!read(lst2,'(1x,A16,i7)')    char,limitGrid
!read(lst2,'(1x,A16,i7)')    char,stepChainProps
!read(lst2,'(1x,A16,i7)')    char,limitChainProps
!read(lst2,'(1x,A16,i3)')    char,nChainConf
!read(lst2,'(1x,A16,f9.2)')  char,timeSteady
!
!read(lst3,'(//, A24)') char
!read(lst3,'(1x,A19,3f9.4)') char,region
!read(lst3,'(1x,A19,3(3x,i5))') char,cells
!read(lst3,'(1x,A19,3(3x,f7.4))') char,gap
!read(lst3,'(1x,A19,5x,i8)') char,maxlist
!read(lst3,'(1x,A19,5x,i8)') char,nFreeAtom
!read(lst3,'(1x,A19,5x,i8)') char,nWallAtom
!read(lst3,'(1x,A19,5x,i8)') char,nAtom
!read(lst3,'(/3x,A27)') char
!
!end