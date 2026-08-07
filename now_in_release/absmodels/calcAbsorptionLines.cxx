// Generic code for calculating absorption lines with various shapes

#define GAUSS 0
#define LORENTZ 1
#define VOIGT 2

#define SQRT2 1.4142135623730950488
#define SQRT2PI sqrt(2*M_PI)

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

// from calcLines.cxx
Real voigtFraction(const Real energy, const Real ecenter, const Real sigma, const Real gamma);

// local routines
void calcAbsorptionLine(const RealArray& energyArray, const RealArray& lineParams,
			const Real crtLevel, const int lineShape, RealArray& fluxArray);

Real absorptionLineBin(const int lineShape, const Real elow, const Real ehigh,
		       const RealArray& lineParams);
RealArray gaussAbsorptionLineIntegrand(const RealArray& x, void *p);
RealArray lorentzAbsorptionLineIntegrand(const RealArray& x, void *p);
RealArray voigtAbsorptionLineIntegrand(const RealArray& x, void *p);


// to calculate line shape. Line is calculated down to (1-crtLevel)
// Arguments:
//      energyArray                 model energy bins
//      inLineParams                line parameter array
//      crtLevelIn                  critical level down which to calculate lines
//                                    can be overridden by xset LINECRITLEVEL
//      lineShape                   0==Gauss, 1==Lorentz, 2==Voigt
//      fluxArray                   input/output fluxes - this routine initializes to 1.0

void calcAbsorptionLine(const RealArray& energyArray, const RealArray& inLineParams,
			const Real crtLevelIn, const int inLineShape, RealArray& fluxArray)
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

  // Find the bin containing the line center. This assumes that the ecenter
  // array is in increasing order of energy. If the center energy falls outside
  // the energy array then icen will be set to < 0. Otherwise
  // energyArray[icen] < ecenter <= energyArray[icen+1]

  Real ecenter = inLineParams[0];
  int icen = BinarySearch(energyArray, ecenter);

  int nE = energyArray.size();
  fluxArray.resize(nE-1);
  fluxArray = 1.0;

  // useful to have non-const lineShape and lineParams variables so we can change them
  // for Voigt special cases
  int lineShape = inLineShape;
  RealArray lineParams(inLineParams);

  // handle the special cases for Voigt where we can use Gauss or Lorentz

  //  if ( (lineShape == VOIGT) && (lineParams[1] <= 0.0) ) {
  //    lineShape = LORENTZ;
  //    lineParams.resize(3);
  //    lineParams[0] = inLineParams[0];
  //    lineParams[1] = inLineParams[2];
  //    lineParams[2] = inLineParams[3];
  //  } else if ( (lineShape == VOIGT) && (lineParams[1] <= 0.0) ) {
  //    lineShape = GAUSS;
  //    lineParams.resize(3);
  //    lineParams[0] = inLineParams[0];
  //    lineParams[1] = inLineParams[1];
  //    lineParams[2] = inLineParams[3];
  //  }

  Real tau;
  if ( lineShape == VOIGT ) {
    tau = lineParams[3];
  } else {
    tau = lineParams[2];
  }

  if ( tau == 0.0 ) return;
  
  // first do case of zero width line. assume for the moment that this occurs
  // if all the lineParams are zero

  bool deltaFunction(true);
  if ( lineShape == VOIGT ) {
    if ( lineParams[1] != 0.0 || lineParams[2] != 0.0 ) deltaFunction = false;
  } else {
    if ( lineParams[1] != 0.0 ) deltaFunction = false;
  }

  if ( deltaFunction ) {
    if ( icen >= 0 ) fluxArray[icen] = exp(-tau);
    return;
  }

  // start at the line center bin and work down in energy till the absorption line
  // fraction differs from 1.0 by less than the critical value

  int ie = icen;
  while ( ie >= 0 ) {
    Real fract = absorptionLineBin(lineShape, energyArray[ie], energyArray[ie+1], lineParams);
    fluxArray[ie] = fract;
    if ( (1.0 - fract) < crtLevel ) ie = 0;
    ie--;
  }

  // repeat for the upper part of the line

  ie = icen+1;
  while ( ie <= nE-2 ) {
    Real fract = absorptionLineBin(lineShape, energyArray[ie], energyArray[ie+1], lineParams);
    fluxArray[ie] = fract;
    if ( (1.0 - fract) < crtLevel ) ie = nE-2;
    ie++;
  }

  return;
}

// Calculates the line between elow and ehigh.
// For all lineShape values lineParams[0] is the center energy
// If lineShape==GAUSS then lineParams[1] is sigma
// If lineShape==LORENTZ then lineParams[1] is width
// If lineShape==VOIGT then lineParams[1] is sigma and lineParams[2] is gamma
// The final element in the lineParams array is tau.

 Real absorptionLineBin(const int lineShape, const Real elow, const Real ehigh,
			const RealArray& lineParams)
{
  static bool first(true);
  static Real saveGaussEcenter, saveGaussWidth, gaussNorm;
  static Real saveLorentzEcenter, saveLorentzWidth, lorentzNorm;
  static Real saveVoigtEcenter, saveVoigtSigma, saveVoigtGamma, voigtNorm;
  static int ModelEvaluations(0);
  Real Precision(1.0e-6);
  Real Integral, IntegralError;

  Real ecenter = lineParams[0];
  
  if ( lineShape == GAUSS ) {
    Real sigma = lineParams[1];
    if ( first || ecenter != saveGaussEcenter || sigma != saveGaussWidth ) {
      saveGaussEcenter = ecenter;
      saveGaussWidth = sigma;
      // line normalization is set so that G(0,inf) = 1.
      gaussNorm = 2.0/(sigma*SQRT2PI*(1.0-erf(-ecenter/(sigma*SQRT2))));
      first = false;
    }
    RealArray p(4);
    p[0] = lineParams[2];
    p[1] = ecenter;
    p[2] = sigma;
    p[3] = gaussNorm;

    ModelEvaluations +=
      AdaptiveIntegrate<gaussAbsorptionLineIntegrand>(elow, ehigh, &p, Precision,
						      Integral, IntegralError);
    Real ebin = fabs(ehigh - elow);
    if ( ebin > 0.0 ) Integral /= ebin;

    return Integral;
    
  } else if ( lineShape == LORENTZ ) {
    Real gamma = lineParams[1];
    if ( first || ecenter != saveLorentzEcenter || gamma != saveLorentzWidth ) {
      saveLorentzEcenter = ecenter;
      saveLorentzWidth = gamma;
      // line normalization is set so that L(0,inf) = 1.
      lorentzNorm = 1.0/(M_PI/2.0 - atan(-2.0*ecenter/gamma));
      first = false;
    }
    RealArray p(4);
    p[0] = lineParams[2];
    p[1] = ecenter;
    p[2] = gamma;
    p[3] = lorentzNorm;

    ModelEvaluations +=
      AdaptiveIntegrate<lorentzAbsorptionLineIntegrand>(elow, ehigh, &p, Precision,
						     Integral, IntegralError);
    Real ebin = fabs(ehigh - elow);
    if ( ebin > 0.0 ) Integral /= ebin;

    return Integral;

  } else if ( lineShape == VOIGT ) {
    Real sigma = lineParams[1];
    Real gamma = lineParams[2];
    if ( first || ecenter != saveVoigtEcenter || sigma != saveVoigtSigma ||
	 gamma != saveVoigtGamma ) {
      saveVoigtEcenter = ecenter;
      saveVoigtSigma = sigma;
      saveVoigtGamma = gamma;
      // line normalization is set so that V(0,inf) = 1.
      voigtNorm = 1.0/(SQRT2PI*sigma);
      Real voigtFracBelowZero = 0.5 - voigtNorm*voigtFraction(0.0, ecenter, sigma, gamma);
      if ( voigtFracBelowZero > 0.0 ) voigtNorm /= (1.0 - voigtFracBelowZero);
      first = false;
    }
    RealArray p(5);
    p[0] = lineParams[3];
    p[1] = ecenter;
    p[2] = sigma;
    p[3] = gamma;
    p[4] = voigtNorm;

    ModelEvaluations +=
      AdaptiveIntegrate<voigtAbsorptionLineIntegrand>(elow, ehigh, &p, Precision,
						     Integral, IntegralError);
    Real ebin = fabs(ehigh - elow);
    if ( ebin > 0.0 ) Integral /= ebin;

    return Integral;

  } else {
    return 0.0;
  }
}

RealArray gaussAbsorptionLineIntegrand(const RealArray& x, void *p)
{

  // Function to evaluate the gaussian absorption line integrand with p used
  // to pass in the parameter values.

  RealArray OutArray(x.size());

  RealArray params(4);
  params = *(static_cast<RealArray*>(p));

  Real tau    = params[0];
  Real center = params[1];
  Real width  = params[2];
  Real gnorm  = params[3];

  for (size_t i=0; i<x.size(); i++) {
    Real gauss = gnorm * exp(-0.5*(x[i]-center)*(x[i]-center)/width/width);
    OutArray[i] = exp(-tau * gauss);
  }

  return OutArray;
}


RealArray lorentzAbsorptionLineIntegrand(const RealArray& x, void *p)
{

  // Function to evaluate the lorentzian absorption line integrand with p used
  // to pass in the parameter values.

  RealArray OutArray(x.size());

  RealArray params(4);
  params = *(static_cast<RealArray*>(p));

  Real tau    = params[0];
  Real center = params[1];
  Real width  = params[2];
  Real lnorm  = params[3];

  Real halfwidth = width/2.0;

  for (size_t i=0; i<x.size(); i++) {
    Real delta = x[i] - center;
    Real lorentz = lnorm * halfwidth / (delta*delta + halfwidth*halfwidth);
    OutArray[i] = exp(-tau * lorentz);
  }

  return OutArray;
}

// evaluate the Voigt absorption line integrand with p used to pass in the parameter values

RealArray voigtAbsorptionLineIntegrand(const RealArray& x, void *p)
{
  RealArray OutArray(x.size());

  RealArray params(5);
  params = *(static_cast<RealArray*>(p));

  Real tau    = params[0];
  Real center = params[1];
  Real sigma  = params[2];
  Real gamma  = params[3];
  Real vnorm  = params[4];

  Real znorm = 1.0/(sigma*SQRT2);
  for (size_t i = 0; i<x.size(); i++) {
    // note that this use of the Faddeeva function assumes that the Lorentzian width
    // is the HWHM so divide the FWHM by 2.
    std::complex<double> z(x[i]-center,gamma/2.0);
    z *= znorm;
    std::complex<double> wofz = Faddeeva::w(z);
    OutArray[i] = exp(-tau * vnorm * wofz.real());
  }

  return OutArray;
}
