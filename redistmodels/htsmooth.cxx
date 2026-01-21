// Convolution model with a half triangle kernel
// htsmooth.cxx: based on the gsmooth.cxx
//
// ver.0.1 2024 Aug 14 Kenji Hamaguchi
// 
// Parameters:
//     0    Largest velocity shift in km/s (negative for the blueshift)

#include <xsTypes.h>
#include <functionMap.h>

// debug
#include <iostream>

// Uses the calcManyLines function from calcLines.cxx

void calcManyLinesLM(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);

extern "C" void htsmooth(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  size_t nE = energyArray.size();
  float clight = 2.99792458e5;

  // construct the array of dE (between the half triangle tip and the center)
  std::vector<RealArray> dE(nE-1);
  RealArray centerEnergy(nE-1);
  for (size_t i=0; i<nE-1; i++) {
    dE[i].resize(1);
    centerEnergy[i] = 0.5*(energyArray[i]+energyArray[i+1]);
    dE[i][0] = -1.0*centerEnergy[i]*params[0]/clight;
  }

  // do the smoothing by just treating the input flux as the sum of a lot of
  // half triangles with the left side at the center energy and the normalization based
  // on the flux in that bin
  // qspeedy is false as the calculation doesn't take computer time.
  
  RealArray inputFlux(flux);
  flux = 0.0;
  calcManyLinesLM(energyArray, centerEnergy, dE, inputFlux, (Real)1.0e-6, 3,
		false, flux);
}
