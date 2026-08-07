// A multiplicative absorption version of the voigt line shape.
// number of model parameters:3
//       0       EL      line energy (in energy units, e.g. keV)
//       1       sigma   gaussian line width (in energy units)
//       2       gamma   lorentzian line width (in energy units)
//       3       Tau     line optical depth
// intrinsic energy range: none
// algorithm:
//       M(E) = exp(-tau*Re[w(z)]/(sigma*sqrt(2pi)))
//          where z = (E-El + i*gamma)/(sigma*sqrt(2))
//          and w(z) = exp(-z^2) erfc(-iz)
//          is the Faddeeva function calculated in XSUtil/Numerics/Faddeeva.{cxx,h}

#include <XSFunctions/functionMap.h>
#include <xsTypes.h>
#include <XSstreams.h>
#include <cmath>
#include <XSFunctions/Utilities/FunctionUtility.h>

// In calcAbsorptionLines.cxx
void calcAbsorptionLine(const RealArray& energyArray, const RealArray& lineParams,
			const Real crtLevel, const int lineShape, RealArray& fluxArray);

extern "C"
void voigtAbsorptionLine(const RealArray& energyArray, const RealArray& params, 
			 int spectrumNumber, RealArray& fluxArray, 
			 RealArray& fluxErrArray, const string& initString)
{
  const Real crtLevel = 1.0e-8;

  int nE = energyArray.size();
  fluxArray.resize(nE-1);
  fluxErrArray.resize(0);
  fluxArray = 1.0;

  calcAbsorptionLine(energyArray, params, crtLevel, 2, fluxArray);

  return;
}
