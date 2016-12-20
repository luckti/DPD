!-----------------------------------------------------------------------
!	dpdflow.fi
!	include file of dpdflow.f90
!-----------------------------------------------------------------------
!	Parameters to difine the size of the program

	  integer NDIM, NHIST   !, MAXHSIZE, MAXWTM, MAXCHAIN, MAXATOM
	  parameter (NDIM = 3, NHIST = 6)
	  logical filexist
	  integer inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
  	  integer moreCycles, nAtom, runId, stepAvg, stepCount, nWallAtom, &
	          nFreeAtom, nStartAtom, PNDIM, stepEquil, stepLimit, nPDP,&
			  randSeed, maxList, hSize, startSample, stepSample, nDp,  &
			  nChain, ChainLen, limitChainProps, stepChainProps,       &
			  countChainProps, nChainend, nChainConf, nDpEnd, stepFzn, &
			  stepStop,MaxlistDp, cellsDp(NDIM),RDFstepSample,BnGrid,  &
			  DiffintStep, DiffStep, JStep, BStep
!             nflow, latticeMod, forceMod, stepAdjustTemp,             &
!             stepInitlzTemp
	  integer countGrid, limitGrid, stepGrid, nsolventend, PstAdj
	  integer initUcell(NDIM), cells(NDIM), DpFzn(10),sizeHistGrid(NDIM)
	
	  real*8 deltaT, density, kinEnergy, pi, potEnergy, pressure, rCut,&
	         rrCut, sKinEnergy, sPressure, sTotEnergy, ssKinEnergy,    &
     		 ssPressure, ssTotEnergy, temperature, timeNow, gravField, &
     		 totEnergy, uSum, vMag, vSum, virSum, sInitKinEnergy,      &
			 wmingap, binvolm, alphaf, alphafp, alphapp, WLadjust,     &
		     cigamaw, gammaw, lambda, sdtinv, alphaw, alphawf, rcutw,  &
			 wLayer, aaDistSq, eeDistSq, radGyrSq, gMomRatio1, mass,   &
			 gMomRatio2, Hfene, rrfene, rmaxfene, reqfene, timeSteady, &
			 shearRate, RdsDp, alphafd, alphaDD, cigamaF, cigamaD,     &
			 gammaF, gammaD, gammaFD, cigamaFD, HPstAdj, r3Cut, cpct,  &
			 rCutDp,rrCutDp,rCutDp2, alphaDA, alphaDB,deltaQ,cellength
!            epslon, nEval, gravity(NDIM) 	
	real * 8 P0, tau, JP0, lpercent,rho0
     		 
	  real*8 region(NDIM), regionH(NDIM), gap(NDIM), CDp0(10,NDIM)
!!-----------allocatable---------------------------------------------------------
!	integer, allocatable :: atomID(:),DpSign(:),ncc(:,:)
!	real*8 , dimension(:,:), allocatable ::  r,rv, ra,raCV,raDP,raRD,raSP,&
!		raCR,rforce,rw,wn,chainCentre,strsGrid,GridChainLen,histGrid !,profileV,profileT,flowvel
!	real*8 , allocatable :: sChain(:)
!!-------------------------------------------------------------------------------------------------------

! ---- MDPD	------------------------------------------------------------	
      real*8 rCut2, alphaB, rrCut2, r3Cut2, alphawB
! ----------------------------------------------------------------------
		
      common / int1 / inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
	  common / int2 / moreCycles, nAtom, nWallAtom, nFreeAtom, runId,  &
		              stepAvg, stepCount, nStartAtom, PNDIM, stepFzn,  &
     		          stepEquil, stepLimit, randSeed, maxList, hSize,  &
		              startSample, stepSample, nChain, ChainLen, nDp,  &
					  nChainend, nChainConf, nDpEnd, nPDP, stepStop
!                     latticeMod, forceMod, stepAdjustTemp, nflow,     &
!                     stepInitlzTemp ,MaxlistDp, cellsDp(NDIM)
	  common / int3 / countGrid, limitGrid, stepGrid, limitChainProps, &
		              stepChainProps, countChainProps, nsolventend,    &
					  PstAdj, RDFstepSample, BnGrid,                   &
					  DiffintStep, DiffStep, JStep, BStep
	  common / int4 /  initUcell, cells, sizeHistGrid,DpFzn
	                  
	  common / real1 / deltaT, density, kinEnergy, pi, potEnergy, vMag,&
                         pressure, rCut, rrCut, sKinEnergy, sPressure, &
					   sTotEnergy, ssKinEnergy, ssPressure, virSum,    &
					   ssTotEnergy, temperature, timeNow, totEnergy,   &
					   uSum, vSum, gravField, sInitKinEnergy, WLadjust,&
		               wmingap, binvolm, alphaf, alphafp, alphapp,     &
					   cigamaw, gammaw, lambda, sdtinv, alphaw,        &
					   alphawf, rcutw, Hfene, rrfene, wLayer, aaDistSq,&
					   eeDistSq, radGyrSq, gMomRatio1, gMomRatio2,     &
					   rmaxfene, reqfene, timeSteady, RdsDp, alphaFD,  &
					   alphaDD, cigamaF, cigamaD, cigamaFD, gammaF,    &
					   gammaD, gammaFD, mass, HPstAdj, r3Cut, cpct ,   &
					   rCutDp, rrCutDp,rCutDp2, alphaDA, alphaDB,	&
					   deltaQ, cellength, P0, tau, JP0, lpercent,rho0
!                      epslon, nEval, gravity 
	  common / real2 /  region, regionH, histGrid, gap,shearRate, CDp0
!                      profileV, profileT, flowvel
! ----- MDPD -----------------------------------------------------------
      common / realm / rCut2, alphaB,rrCut2, r3Cut2, alphawB

!-----------allocatable---------------------------------------------------------
	common / int5 / atomID, ncc,DpSign
	common / real3 / r, rv, ra, rw,raCV, raDP, raRD, raSP, raCR,rforce,&
                      strsGrid,chainCentre, sChain, GridChainLen, wn
	common / RDF / nS,nSB
	common / RDFGrid / BmGrid, StrGrid, StateGrid
	common / SGrid / StateSpGrid,StrSpGrid
	integer, pointer :: atomID(:),DpSign(:),ncc(:,:),nS(:,:),nSB(:,:)
	real*8 , dimension(:,:), pointer ::  r,rv, ra,raCV,raDP,raRD,raSP,raCR,&
		    rforce,rw,wn,chainCentre,strsGrid,GridChainLen,histGrid !,profileV,profileT,flowvel
	real * 8, pointer::sChain(:), BmGrid(:, : , : , : ), StrGrid(:, : , : , : ), &
						StateGrid(:, : , : , : ),StateSpGrid(:,:),StrSpGrid(:,:)