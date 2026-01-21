// Empirical Convolution model for a general Compton shoulder shape
// emcomp.cxx: based on the htsmooth.cxx
//
// ver. date author
// 0.1 2024 Aug 17 Kenji Hamaguchi original (use calcLinesLM_v02.cxx)
// 0.2 2024 Sept 5 Kenji Hamaguchi add parameter 3 for high-E cut (use calcLinesLM_v03.cxx)
// 0.3 2024 Sept 15 Kenji Hamaguchi Change the parameter sets
// 0.4 2025 Dec 11 Kenji Hamaguchi, change the function name from compsh to emcomp
// 
// Parameters:
//     0    Ratio of the Compton shoulder component to the original model
//     1    y value at x=0 of the linear redistribution function
//     2    y value at x=1 of the linear redistribution function
//     3    High-E cut (E ratio)

#include <xsTypes.h>
#include <functionMap.h>

// debug
#include <iostream>

// Uses the calcManyLines function from calcLines.cxx

void calcManyLinesLM(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);

extern "C" void emcomp(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  size_t nE = energyArray.size();
  float Ee_keV = 510.997876;
  size_t nparams = params.size();
  // std::cout << "Number of Params: " << nparams << std::endl;
  
  // construct the array of dE (maximum energy shift realized with the 180 degree backscattering)
  // It's 1 bin smaller than energyArray because it's between the energyArray
  std::vector<RealArray> dE(nE-1+nparams);
  RealArray centerEnergy(nE-1);
  for (size_t i=0; i<nE-1; i++) {
    dE[i].resize(1);
    centerEnergy[i] = 0.5*(energyArray[i]+energyArray[i+1]);
    dE[i][0] = -1.0*centerEnergy[i]/(1 + Ee_keV/(2*centerEnergy[i]));
  }

  // std::cout << "compsh:ProcA" << std::endl;
  
  // add the parameter values at the end of the dE array
  for (size_t i=0; i<nparams; i++) {
    // std::cout << "params: " << i << " " << params[i] << std::endl;
    dE[nE-1+i].resize(1);
    dE[nE-1+i][0] = params[i];
    // std::cout << "dE[" << nE-1+i << "][0]= " << dE[nE-1+i][0] << std::endl;
  }
  
  // qspeedy is false as the calculation doesn't take computer time.
  
  RealArray inputFlux(flux);
  flux = 0.0;
  calcManyLinesLM(energyArray, centerEnergy, dE, inputFlux, (Real)1.0e-6, 4,
		false, flux);
}
