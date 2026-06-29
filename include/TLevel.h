#ifndef TLEVEL_H
#define TLEVEL_H

/////////////////////////////////////////////////////////////////
///
/// \class TLevel
///
/// This class represents a level in a nucleus. Probably not of
/// much use on it's own, but it is part of the framework to
/// provide information about a nucleus (see also TNucleus).
///
/////////////////////////////////////////////////////////////////

#include <map>

#include "TObject.h"

#include "Globals.h"
#include "TTransition.h"

class TNucleus;

class TLevel : public TObject {
public:
   explicit TLevel(TNucleus* nucleus = nullptr, const double& energy = -1., std::string label = "");
   TLevel(TNucleus* nucleus, const double& energy, const double& energyUncertainty = 0.);
   TLevel(const TLevel& rhs);
   TLevel(TLevel&& rhs) noexcept = default;
   TLevel& operator=(const TLevel& rhs);
   TLevel& operator=(TLevel&& rhs) noexcept = default;
   ~TLevel()                                = default;

   TTransition* AddTransition(double levelEnergy, double br = 100., double ts = 1.);
   TTransition* AddTransition(double levelEnergy, double energyUncertainty, double br = 100., double ts = 1.);
   TLevel*      AddTransition(TTransition* transition, bool quiet = false);

   void Energy(const double val) { fEnergy = val; }
   void EnergyUncertainty(const double val) { fEnergyUncertainty = val; }
   void Label(const char* val) { fLabel = val; }
   void AddFeeding() { ++fNofFeeding; }

   double      Energy() const { return fEnergy; }
   double      EnergyUncertainty() const { return fEnergyUncertainty; }
   std::string Label() const { return fLabel; }

   std::pair<double, double> GetMinMaxTransition() const;
   size_t                    NofDrainingTransitions() const { return fTransitions.size(); }
   size_t                    NofFeedingTransitions() const { return fNofFeeding; }

   std::map<double, TTransition>::iterator       begin() { return fTransitions.begin(); }
   std::map<double, TTransition>::iterator       end() { return fTransitions.end(); }
   std::map<double, TTransition>::const_iterator begin() const { return fTransitions.begin(); }
   std::map<double, TTransition>::const_iterator end() const { return fTransitions.end(); }

   // comparison operators (level-level, level-double, and double-level)
   friend bool operator<(const TLevel& lhs, const TLevel& rhs) { return lhs.fEnergy < rhs.fEnergy; }
   friend bool operator>(const TLevel& lhs, const TLevel& rhs) { return rhs < lhs; }
   friend bool operator<=(const TLevel& lhs, const TLevel& rhs) { return !(rhs < lhs); }
   friend bool operator>=(const TLevel& lhs, const TLevel& rhs) { return !(lhs < rhs); }

   friend bool operator<(const TLevel& lhs, const double& rhs) { return lhs.fEnergy < rhs; }
   friend bool operator>(const TLevel& lhs, const double& rhs) { return rhs < lhs; }
   friend bool operator<=(const TLevel& lhs, const double& rhs) { return !(rhs < lhs); }
   friend bool operator>=(const TLevel& lhs, const double& rhs) { return !(lhs < rhs); }

   friend bool operator<(const double& lhs, const TLevel& rhs) { return lhs < rhs.fEnergy; }
   friend bool operator>(const double& lhs, const TLevel& rhs) { return rhs < lhs; }
   friend bool operator<=(const double& lhs, const TLevel& rhs) { return !(rhs < lhs); }
   friend bool operator>=(const double& lhs, const TLevel& rhs) { return !(lhs < rhs); }

   void Print(Option_t* option = "") const override;

   static void       Verbosity(EVerbosity val) { fVerbosity = val; }
   static EVerbosity Verbosity() { return fVerbosity; }

private:
   static EVerbosity             fVerbosity;               ///< verbosity of the class
   double                        fEnergy{0.};              ///< energy of this level
   double                        fEnergyUncertainty{0.};   ///< energy uncertainty of this level
   std::string                   fLabel;                   ///< label for this level
   std::map<double, TTransition> fTransitions;             ///< transitions draining this level, each pointing to a level
   size_t                        fNofFeeding{0};           ///< counter for gammas feeding this level
   TNucleus*                     fNucleus{nullptr};        ///< pointer to the nucleus this level belongs to

   /// \cond CLASSIMP
   ClassDefOverride(TLevel, 1)   // NOLINT(readability-else-after-return)
   /// \endcond
};

#endif
