#ifndef TTRANSITION_H
#define TTRANSITION_H

/** \addtogroup Fitting Fitting & Analysis
 *  @{
 */

#include <cstdio>
#include <limits>

#include "TObject.h"

#include "Globals.h"

/////////////////////////////////////////////////////////////////
///
/// \class TTransition
///
/// This Class contains the information about a nuclear
/// transition. These transitions are a part of a TNucleus
/// and are typically set within the TNucleus framework
///
/////////////////////////////////////////////////////////////////

class TNucleus;
class TLevel;

class TTransition : public TObject {
   friend class TNucleus;

public:
   TTransition()                       = default;
   TTransition(double energy, double energyUncertainty, double intensity, double intensityUncertainty, double mixingRatio, double mixingRatioUncertainty, double conversionCoeff, double conversionCoeffUncertainty, double totalIntensity, double totalIntensityUncertainty);
   TTransition(TLevel* level, double energy, double energyUncertainty, double intensity, double intensityUncertainty, double mixingRatio, double mixingRatioUncertainty, double conversionCoeff, double conversionCoeffUncertainty, double totalIntensity, double totalIntensityUncertainty);
   TTransition(TLevel* level, double energy, double energyUncertainty, double intensity, double intensityUncertainty, double totalIntensity, double totalIntensityUncertainty);
   TTransition(const TTransition&)     = default;
   TTransition(TTransition&&) noexcept = default;
   ~TTransition()                      = default;

   TTransition& operator=(const TTransition&)     = default;
   TTransition& operator=(TTransition&&) noexcept = default;

   bool IsSortable() const override { return true; }
   int  Compare(const TObject* obj) const override;
   int  CompareIntensity(const TObject* obj) const;
   int  CompareEnergy(const TObject* obj) const;

   void SetEnergy(const double& tmpenergy) { fEnergy = tmpenergy; }
   void SetEnergyUncertainty(const double& tmperror) { fEnergyUncertainty = tmperror; }
   void SetIntensity(const double& tmpintens) { fIntensity = tmpintens; }
   void SetIntensityUncertainty(const double& tmpinterror) { fIntensityUncertainty = tmpinterror; }
   void SetCompareIntensity(const bool& val) { fCompareIntensity = val; }

   double GetEnergy() const { return fEnergy; }
   double GetEnergyUncertainty() const { return fEnergyUncertainty; }
   double GetIntensity() const { return fIntensity; }
   double GetIntensityUncertainty() const { return fIntensityUncertainty; }

   void Energy(double val) { fEnergy = val; }
   void EnergyUncertainty(double val) { fEnergyUncertainty = val; }
   void Intensity(double val) { fIntensity = val; }
   void IntensityUncertainty(double val) { fIntensityUncertainty = val; }
   void MixingRatio(double val) { fMixingRatio = val; }
   void MixingRatioUncertainty(double val) { fMixingRatioUncertainty = val; }
   void ConversionCoefficent(double val) { fConversionCoeff = val; }
   void ConversionCoefficentUncertainty(double val) { fConversionCoeffUncertainty = val; }
   void TotalIntensity(double val) { fTotalIntensity = val; }
   void TotalIntensityUncertainty(double val) { fTotalIntensityUncertainty = val; }
   void FinalLevel(double val) { fFinalLevel = val; }
   void UncertainPlacement(bool val) { fUncertainPlacement = val; }
                                                                                   
   double Energy() const { return fEnergy; }
   double EnergyUncertainty() const { return fEnergyUncertainty; }
   double Intensity() const { return fIntensity; }
   double IntensityUncertainty() const { return fIntensityUncertainty; }
   double MixingRatio() const { return fMixingRatio; }
   double MixingRatioUncertainty() const { return fMixingRatioUncertainty; }
   double ConversionCoefficent() const { return fConversionCoeff; }
   double ConversionCoefficentUncertainty() const { return fConversionCoeffUncertainty; }
   double TotalIntensity() const { return fTotalIntensity; }
   double TotalIntensityUncertainty() const { return fTotalIntensityUncertainty; }
   double FinalLevel() const { return fFinalLevel; }
   bool UncertainPlacement() const { return fUncertainPlacement; }
                                                                                   
   void Clear(Option_t* opt = "") override;
   void Print(Option_t* opt = "") const override;

   std::string PrintToString() const;

   bool operator>(const TTransition& rhs) const { return GetEnergy() > rhs.GetEnergy(); }
   bool operator<(const TTransition& rhs) const { return GetEnergy() < rhs.GetEnergy(); }

   static void Verbosity(EVerbosity val) { fVerbosity = val; }
   static EVerbosity Verbosity() { return fVerbosity; }

private:
   static EVerbosity fVerbosity;                                                              ///< verbosity of the class
   TLevel*           fLevel{nullptr};                                                         ///< The level this transition drains
   double            fEnergy{0.};                                                             ///< Energy of the transition
   double            fEnergyUncertainty{std::numeric_limits<double>::quiet_NaN()};            ///< Uncertainty in the energy of the transition
   double            fIntensity{0.};                                                          ///< Intensity of the transition
   double            fIntensityUncertainty{std::numeric_limits<double>::quiet_NaN()};         ///< Uncertainty in the intensity
   double            fMixingRatio{0.};                                                        ///< Mixing ratio of the transition
   double            fMixingRatioUncertainty{std::numeric_limits<double>::quiet_NaN()};       ///< Uncertainty in the mixing ratio
   double            fConversionCoeff{0.};                                                    ///< Conversion coefficent of the transition
   double            fConversionCoeffUncertainty{std::numeric_limits<double>::quiet_NaN()};   ///< Uncertainty in the conversion coefficent
   double            fTotalIntensity{0.};                                                     ///< Total intensity of the transition (includes conversion electrons for gamma-rays)
   double            fTotalIntensityUncertainty{std::numeric_limits<double>::quiet_NaN()};    ///< Uncertainty in the total intensity
   double            fFinalLevel{std::numeric_limits<double>::quiet_NaN()};                   ///< Final level for this transition
   bool              fUncertainPlacement{false};                                              ///< Flag to indicate placement of this transition in the level scheme is uncertain

   bool fCompareIntensity{true};   ///< Whether to sort by intensity or energy

   /// \cond CLASSIMP
   ClassDefOverride(TTransition, 0)   // NOLINT(readability-else-after-return)
   /// \endcond
};
/*! @} */
#endif
