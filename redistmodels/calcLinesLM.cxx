// Generic code for calculating lines with various shapes
// calcLinesLM.cxx: changed from calcLines.cxx for a local model input
//
// ver.0.1 2024 Aug 14 Kenji Hamaguchi
// ver 0.2 2024 Aug 17 Kenji Hamaguchi  Add lines for compsh.cxx
//                                      compsh.cxx adds line profile params after the maximum energy shift arrays.
// ver 0.3 2024 Sept 5 Kenji Hamaguchi  Change the formula of the comptonshoulderLM and add one more parameter for high-E cut
// ver 0.4 2024 Sept 15 Kenji Hamaguchi  Fixing a bug of the comptonshoulderLM function
// ver 0.5 2025 Sept 19 Kenji Hamaguchi  Adding the BScut parameter to handle the forward scattering case as well.
//                                       The function returns zero if BScut is lower than FScut
//                                       (though the parameter range should handle it.)
//                                       Changed the parameter names used in comptonshoulderLM as well, so it's not confused.
// ver 0.6 2025 Dec 11 Kenji Hamaguchi Remove the functions in the original calcLines.cxx

#define HALFTRIANGLE 3
#define COMPTONSHOULDER 4

#define SQRT2 1.4142135623730950488

#include <xsTypes.h>
#include <XSstreams.h>
#include <XSUtil/Numerics/AdaptiveIntegrate.h>
#include <XSUtil/Numerics/BinarySearch.h>
#include <XSUtil/Numerics/Faddeeva.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <sstream>
#include <iostream>

using std::atan;
using namespace Numerics;

void calcManyLinesLM(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);
void calcLineLM(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

Real lineFractionLM(const int lineShape, const Real energy, const Real ecenter, 
		  const RealArray& lineParams, const bool qspeedy);
Real halftriangleFractionLM(const Real delta);
Real comptonshoulderLM(const Real delta, const Real y0, const Real y1, const Real E0, const Real E1);


// to calculate line shapes. If qspeedy is true then uses a tabulation
// otherwise does a proper calculation. Lines are calculated down to
// crtlevel*(flux in center bin)
// Arguments:
//      energyArray                 model energy bins
//      ecenterArray                array of line center energies
//      lineParamsArray             vector of line parameter arrays
//      linefluxArray               array of line fluxes
//      crtLevelIn                  critical level down which to calculate lines
//                                    can be overridden by xset LINECRITLEVEL
//      lineShape                   0==Gauss, 1==Lorentz, 2==Voigt 3==HalfTriangle
//      qspeedy                     if true then do fast but less accurate calculation
//      fluxArray                   input/output fluxes - note that this routine
//                                    does not initialize input flux to zero but uses +=

void calcManyLinesLM(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevelIn, 
		   const int inLineShape, const bool qspeedy, RealArray& fluxArray)
{
  // Find out whether we want to override the crtLevel

  Real crtLevel(crtLevelIn);
  string pname = "LINECRITLEVEL";
  string pvalue = FunctionUtility::getModelString(pname);
  if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(pvalue);
    Real tmpVal(crtLevel);
    if (!(iss >> tmpVal) || !iss.eof()) {
      std::ostringstream oss;
      oss << "Failed to read value from LINECRITLEVEL - assuming "
	  << crtLevel << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
    } else {
      crtLevel = tmpVal;
    }
  }
  
  // Find the bins containing the line centers. This assumes that the ecenter
  // array is in increasing order of energy. If the center energy falls outside
  // the energy array then lineCenterBin will be set to < 0. Otherwise
  // energyArray[lineCenterBin] < ecenterArray <= energyArray[lineCenterBin+1]

  IntegerVector lineCenterBin(ecenterArray.size());
  lineCenterBin = BinarySearch(energyArray, ecenterArray);

  Real efirst = energyArray[0];
  int nE = energyArray.size();
  Real elast = energyArray[nE-1];

  // useful to have a non-const lineShape variable so we can change it
  // for Voigt special cases
  int lineShape = inLineShape;

  // Loop around the lines ignoring those with zero line flux

  for (size_t iline=0; iline<ecenterArray.size(); iline++) {

    Real lineflux = linefluxArray[iline];
    if ( lineflux == 0.0 ) continue;

    int icen = lineCenterBin[iline];
    Real ecenter = ecenterArray[iline];

    // handle the special cases for Voigt where we can use Gauss or Lorentz
    

    RealArray lineParams;
    if ( lineShape == COMPTONSHOULDER ) {
      // 5 is the dE param (1) + the number of the compsh function params
      // std::cout << "calcLineLM:calcManyLinesLM:lineParamsSet" << std::endl;
      lineParams.resize(5);
      lineParams[0] = lineParamsArray[iline][0];
      lineParams[1] = lineParamsArray[nE-1][0];
      lineParams[2] = lineParamsArray[nE][0];
      lineParams[3] = lineParamsArray[nE+1][0];
      lineParams[4] = lineParamsArray[nE+2][0];
      lineParams[5] = lineParamsArray[nE+3][0];
      // Show the parameter for iline == 10 for debugging
      //std::cout << "iline:" << iline << std::endl;
      //std::cout << "lineParams0:" << lineParams[0] << std::endl;
      //std::cout << "lineParams1:" << lineParams[1] << std::endl;
      //std::cout << "lineParams2:" << lineParams[2] << std::endl;
      //std::cout << "lineParams3:" << lineParams[3] << std::endl;
    } else {
      lineParams.resize(lineParamsArray[iline].size());
      lineParams = lineParamsArray[iline];
    }

    // first do case of zero width line. assume for the moment that this occurs
    // if all the lineParams are zero

    bool deltaFunction(true);
    for (size_t i=0; i<lineParams.size(); i++) {
      if (lineParams[i] != 0.0) {
	deltaFunction = false;
	break;
      }
    }

    if ( deltaFunction ) {
      if ( icen >= 0 ) fluxArray[icen] += lineflux;
      continue;
    }

    // If the line center is below first bin then don't calculate the 
    // lower part of the line. If line center is above the first bin
    // then just calculate part of line within the energy range.

    if ( ecenter >= efirst ) {
      int ielow = icen;
      Real alow = 0.0;
      Real fractionInsideRange = lineFractionLM(lineShape, efirst, ecenter, lineParams, qspeedy);
      if ( ecenter > elast ) {
	ielow = nE-2;
	alow = lineFractionLM(lineShape, elast, ecenter, lineParams, qspeedy);
	fractionInsideRange -= alow;
      }

      // Do the low energy part of the line

      Real lineSum = 0.0;
      Real ahi;
      while ( ielow >= 0 ) {
	ahi = lineFractionLM(lineShape, energyArray[ielow], ecenter, lineParams, qspeedy);
	Real fract = ahi-alow;

	fluxArray[ielow] += fract*lineflux;
	lineSum += fract;
	if ( (fractionInsideRange-lineSum) < crtLevel ) {
	  // Too many sigma away so stop now and add the rest of the line into this
	  // bin. Not strictly correct but shouldn't matter and ensures that the total
	  // flux is preserved.
	  fluxArray[ielow] += (fractionInsideRange-lineSum)*lineflux;
	  ielow = 0;
	}
	alow = ahi;
	ielow -= 1;
      }
    }

    // If the line center is above the last bin then don't calculate 
    // the upper part of the line. If line center is below first bin then 
    // just calculate the part of the line within energy range.

    if ( ecenter <= elast ) {
      Real alow = 0.0;
      int ielow = icen;
      Real fractionInsideRange = lineFractionLM(lineShape, elast, ecenter, lineParams, qspeedy);
      if ( ecenter <  efirst ) {
	ielow = 0;
	alow = lineFractionLM(lineShape, energyArray[ielow], ecenter, lineParams, qspeedy);
	fractionInsideRange -= alow;
      }

      // Do the high energy part of the line

      Real lineSum = 0.0;
      Real ahi;
      while ( ielow < nE-1 ) {
	ahi = lineFractionLM(lineShape, energyArray[ielow+1], ecenter, lineParams, qspeedy);
	Real fract = ahi-alow;
	fluxArray[ielow] += fract*lineflux;
	lineSum += fract;
	if ( (fractionInsideRange-lineSum) < crtLevel ) {
	  // Too many sigma away so stop now and add the rest of the Gaussian into this
	  // bin. Not strictly correct but shouldn't matter and ensures that the total
	  // flux is preserved.
	  fluxArray[ielow] += (fractionInsideRange-lineSum)*lineflux;
	  ielow = nE;
	}
	alow = ahi;
	ielow += 1;
      }
    }
	  
  }

  return;
}

void calcLineLM(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray)
{
  RealArray ecenterArray(ecenter,1);
  std::vector<RealArray> lineParamsArray(1,lineParams);
  RealArray linefluxArray(lineflux,1);
  calcManyLinesLM(energyArray, ecenterArray, lineParamsArray, linefluxArray, 
		crtLevel, lineShape, qspeedy, fluxArray);
  return;
}

// Calculates the fraction of the line between energy and ecenter.
// If lineShape==HALFTRIANGLE then lineParams[0] is maximum blueshift in energy

Real lineFractionLM(const int lineShape, const Real energy, const Real ecenter,
		  const RealArray& lineParams, const bool qspeedy)
{
  if ( lineShape == HALFTRIANGLE ) {
    return halftriangleFractionLM((energy-ecenter)/lineParams[0]);
  } else if ( lineShape == COMPTONSHOULDER ) {
    return comptonshoulderLM((energy-ecenter)/lineParams[0], lineParams[2], lineParams[3], lineParams[4], lineParams[5]) * lineParams[1];
  } else {
    return 0.0;
  }
}

Real halftriangleFractionLM(const Real delta)
{
  if ( delta < 0.0 ){
    return 0.0;
  } else if ( delta > 1.0 ){
    return 1.0;
  } else {
    return (2.0 - delta) * delta;
  }
}

Real comptonshoulderLM(const Real delta, const Real E0, const Real y0, const Real E1, const Real y1)
{
  // if the E0 gets lower than the E1, the function always return zero.
  if ( E0 >= E1 ){
    return 0.0;
  } else if ( delta < E0 ){
    return 0.0;
  } else if ( delta > E1 ){
    return 1.0;
  } else {
    Real ydelta = (y0 * (E1 - delta) + y1 * (delta - E0))/(E1-E0);
    return ((y0 + ydelta) * (delta - E0))/((y0 + y1) * (E1 - E0));
  }
}

