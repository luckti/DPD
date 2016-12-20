!-----------------------------------------------------------------------
!	dpdpolym.h
!	include file of dpdpolym.f90
!-----------------------------------------------------------------------
!	Parameters to difine the size of the program

	  integer NDIM, MAXATOM, NHIST, MAXHSIZE, MAXWTM, MAXCHAIN
	  parameter (NDIM = 3, NHIST = 6,                                  &
     		     MAXATOM = 600000, MAXCHAIN = 100000,                  &
		         MAXHSIZE = 3000000, MAXWTM = 40000)

	  integer inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
  	  integer moreCycles, nAtom, runId, stepAvg, stepCount, nWallAtom, &
	          nFreeAtom, nStartAtom, PNDIM, stepEquil, stepLimit, nPDP,&
			  randSeed, maxList, hSize, startSample, stepSample, nDp,  &
			  nChain, ChainLen, limitChainProps, stepChainProps,       &
			  countChainProps, nChainend, nChainConf, nDpEnd, stepFzn, &
			  stepStop
!             nflow, latticeMod, forceMod, stepAdjustTemp,             &
!             stepInitlzTemp
	  integer countGrid, limitGrid, stepGrid, nsolventend, PstAdj
	  integer atomID(MAXATOM), initUcell(NDIM), cells(NDIM), DpFzn(150),&
	          sizeHistGrid(NDIM), ncc(MAXATOM,NDIM)
	
	  real*8 deltaT, density, kinEnergy, pi, potEnergy, pressure, rCut,&
	         rrCut, sKinEnergy, sPressure, sTotEnergy, ssKinEnergy,    &
     		 ssPressure, ssTotEnergy, temperature, timeNow, gravField, &
     		 totEnergy, uSum, vMag, vSum, virSum, sInitKinEnergy,      &
			 wmingap, binvolm, alphaf, alphafp, alphapp, WLadjust,     &
		     cigamaw, gammaw, lambda, sdtinv, alphaw, alphawf, rcutw,  &
			 wLayer, aaDistSq, eeDistSq, radGyrSq, gMomRatio1, mass,   &
			 gMomRatio2, Hfene, rrfene, rmaxfene, reqfene, timeSteady, &
			 shearRate, RdsDp, alphafd, alphaDD, cigamaF, cigamaD,     &
			 gammaF, gammaD, gammaFD, cigamaFD, HPstAdj,rCut2,rCutDp,  &
			 rCutDp2, alphaB, alphaDA, alphaDB
!            epslon, nEval, gravity(NDIM)
	  real*8 r(MAXATOM,NDIM), rv(MAXATOM,NDIM), ra(MAXATOM,NDIM),      &
     		 rw(MAXWTM,NDIM), chainCentre(MAXCHAIN,NDIM),              &
			 sChain(MAXCHAIN), raCV(MAXATOM,NDIM), raDP(MAXATOM,NDIM), &
			 raRD(MAXATOM,NDIM), raSP(MAXATOM,NDIM)
	  real*8 region(NDIM), regionH(NDIM), histGrid(MAXHSIZE,NHIST),    &
		     strsGrid(MAXHSIZE,7), GridChainLen(MAXHSIZE,2), gap(NDIM),&
			 rforce(MAXATOM,6), wn(MAXWTM,NDIM), CDp0(150,NDIM)
!            profileV(MAXHSIZE), profileT(MAXHSIZE), flowvel(MAXHSIZE,3)	

	  integer DpsizeHistGrid(NDIM), Dphsize, drp2
	  real*8 CDp(10,NDIM), DpregionH(NDIM), DphistGrid(MAXHSIZE,NHIST),&
	         Dpbinvolm

      common / int1 / inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
	  common / int2 / moreCycles, nAtom, nWallAtom, nFreeAtom, runId,  &
		              stepAvg, stepCount, nStartAtom, PNDIM, stepFzn,  &
     		          stepEquil, stepLimit, randSeed, maxList, hSize,  &
		              startSample, stepSample, nChain, ChainLen, nDp,  &
					  nChainend, nChainConf, nDpEnd, nPDP, stepStop
!                     latticeMod, forceMod, stepAdjustTemp, nflow,     &
!                     stepInitlzTemp 
	  common / int3 / countGrid, limitGrid, stepGrid, limitChainProps, &
		              stepChainProps, countChainProps, nsolventend,    &
					  PstAdj
	  common / int4 / atomID, initUcell, cells, sizeHistGrid, ncc, DpFzn

	  common / real1 / deltaT, density, kinEnergy, pi, potEnergy, vMag,&
     		           pressure, rCut, rrCut, sKinEnergy, sPressure,   &
					   sTotEnergy, ssKinEnergy, ssPressure, virSum,    &
					   ssTotEnergy, temperature, timeNow, totEnergy,   &
					   uSum, vSum, gravField, sInitKinEnergy, WLadjust,&
		               wmingap, binvolm, alphaf, alphafp, alphapp,     &
					   cigamaw, gammaw, lambda, sdtinv, alphaw,        &
					   alphawf, rcutw, Hfene, rrfene, wLayer, aaDistSq,&
					   eeDistSq, radGyrSq, gMomRatio1, gMomRatio2,     &
					   rmaxfene, reqfene, timeSteady, RdsDp, alphaFD,  &
					   alphaDD, cigamaF, cigamaD, cigamaFD, gammaF,    &
					   gammaD, gammaFD, mass, HPstAdj,rCut2,rCutDp,    &
					   rCutDp2, alphaB, alphaDA, alphaDB
!                      epslon, nEval, gravity 
	  common / real2 / r, rv, ra, rw, region, regionH, histGrid, gap,  &
	                   rforce, strsGrid, chainCentre, sChain,          &
					   GridChainLen, wn, shearRate, raCV, raDP, raRD,  &
					   raSP, CDp0
!                      profileV, profileT, flowvel 

      common / dropi / DpsizeHistGrid, Dphsize, drp2
	  common / dropr / CDp, DpregionH, DphistGrid, Dpbinvolm
!-----------------------------------------------------------------------