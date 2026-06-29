#include "TTransition.h"

#include <iostream>
#include <cmath>

#include "TString.h"

EVerbosity TTransition::fVerbosity = EVerbosity::kQuiet;

TTransition::TTransition(TLevel* level, double energy, double energyUncertainty, double intensity, double intensityUncertainty, double mixingRatio, double mixingRatioUncertainty, double conversionCoeff, double conversionCoeffUncertainty, double totalIntensity, double totalIntensityUncertainty)
   : fLevel(level), fEnergy(energy), fEnergyUncertainty(energyUncertainty), fIntensity(intensity), fIntensityUncertainty(intensityUncertainty),
   fMixingRatio(mixingRatio), fMixingRatioUncertainty(mixingRatioUncertainty), fConversionCoeff(conversionCoeff), fConversionCoeffUncertainty(conversionCoeffUncertainty),
   fTotalIntensity(totalIntensity), fTotalIntensityUncertainty(totalIntensityUncertainty)
{
}

TTransition::TTransition(TLevel* level, double energy, double energyUncertainty, double intensity, double intensityUncertainty, double totalIntensity, double totalIntensityUncertainty)
   : fLevel(level), fEnergy(energy), fEnergyUncertainty(energyUncertainty), fIntensity(intensity), fIntensityUncertainty(intensityUncertainty),
   fTotalIntensity(totalIntensity), fTotalIntensityUncertainty(totalIntensityUncertainty)
{
}

TTransition::TTransition(double energy, double energyUncertainty, double intensity, double intensityUncertainty, double mixingRatio, double mixingRatioUncertainty, double conversionCoeff, double conversionCoeffUncertainty, double totalIntensity, double totalIntensityUncertainty)
   : fEnergy(energy), fEnergyUncertainty(energyUncertainty), fIntensity(intensity), fIntensityUncertainty(intensityUncertainty),
   fMixingRatio(mixingRatio), fMixingRatioUncertainty(mixingRatioUncertainty), fConversionCoeff(conversionCoeff), fConversionCoeffUncertainty(conversionCoeffUncertainty),
   fTotalIntensity(totalIntensity), fTotalIntensityUncertainty(totalIntensityUncertainty)
{
}

void TTransition::Clear(Option_t*)
{
   fEnergy                     = 0.;
   fEnergyUncertainty          = std::numeric_limits<double>::quiet_NaN();
   fIntensity                  = 0.;
   fIntensityUncertainty       = std::numeric_limits<double>::quiet_NaN();
   fMixingRatio                = 0.;
   fMixingRatioUncertainty     = std::numeric_limits<double>::quiet_NaN();
   fConversionCoeff            = 0.;
   fConversionCoeffUncertainty = std::numeric_limits<double>::quiet_NaN();
   fTotalIntensity             = 0.;
   fTotalIntensityUncertainty  = std::numeric_limits<double>::quiet_NaN();
}

void TTransition::Print(Option_t*) const
{
   std::cout << "Energy:    " << fEnergy;
   if(!std::isnan(fEnergyUncertainty)) {
      std::cout << " +/- " << fEnergyUncertainty;
   }
   std::cout << std::endl;
   if(!std::isnan(fIntensity)) {
      std::cout << "\tIntensity: " << fIntensity;
      if(!std::isnan(fIntensityUncertainty)) {
         std::cout << " +/- " << fIntensityUncertainty;
      }
      std::cout << std::endl;
   }
   if(!std::isnan(fTotalIntensity)) {
      std::cout << "\tTotal intensity (including conversion electrons): " << fTotalIntensity;
      if(!std::isnan(fTotalIntensityUncertainty)) {
         std::cout << " +/- " << fTotalIntensityUncertainty;
      }
      std::cout << std::endl;
   }
   if(!std::isnan(fMixingRatio)) {
      std::cout << "\tMixing ratio: " << fMixingRatio;
      if(!std::isnan(fMixingRatioUncertainty)) {
         std::cout << " +/- " << fMixingRatioUncertainty;
      }
      std::cout << std::endl;
   }
   if(!std::isnan(fConversionCoeff)) {
      std::cout << "\tConversion coefficent: " << fConversionCoeff;
      if(!std::isnan(fConversionCoeffUncertainty)) {
         std::cout << " +/- " << fConversionCoeffUncertainty;
      }
      std::cout << std::endl;
   }
}

std::string TTransition::PrintToString() const
{
   std::string toString;
   toString.append(Form("%f\t", fEnergy));
   toString.append(Form("%f\t", fEnergyUncertainty));
   toString.append(Form("%f\t", fIntensity));
   toString.append(Form("%f\t", fIntensityUncertainty));

   return toString;
}

int TTransition::Compare(const TObject* obj) const
{
   if(fCompareIntensity) {
      return CompareIntensity(obj);
   }
   return CompareEnergy(obj);
}

int TTransition::CompareIntensity(const TObject* obj) const
{
   /// Compares the intensities of the TTransitions
   if(fIntensity > static_cast<const TTransition*>(obj)->fIntensity) {
      return -1;
   }
   if(fIntensity == static_cast<const TTransition*>(obj)->fIntensity) {
      return 0;
   }
   return 1;
}

int TTransition::CompareEnergy(const TObject* obj) const
{
   /// Compares the energies of the TTransitions
   if(fEnergy < static_cast<const TTransition*>(obj)->fEnergy) {
      return -1;
   }
   if(fEnergy == static_cast<const TTransition*>(obj)->fEnergy) {
      return 0;
   }
   return 1;
}
