#include "TLevel.h"
#include "TNucleus.h"

#include <iostream>

EVerbosity TLevel::fVerbosity = EVerbosity::kDefault;

TLevel::TLevel(TNucleus* nucleus, const double& energy, std::string label)
   : fEnergy(energy), fLabel(std::move(label)), fNucleus(nucleus)
{
   if(nucleus == nullptr) {
      std::cout << "Warning created new level w/o a parent nucleus!" << std::endl;
   }
}
TLevel::TLevel(TNucleus* nucleus, const double& energy, const double& energyUncertainty, const char& identifier)
   : fEnergy(energy), fEnergyUncertainty(energyUncertainty), fNucleus(nucleus), fIdentifier(identifier)
{
}

TLevel::TLevel(const TLevel& rhs)
   : TObject(rhs), fEnergy(rhs.fEnergy), fEnergyUncertainty(rhs.fEnergyUncertainty), fLabel(rhs.fLabel),
     fTransitions(rhs.fTransitions), fNofFeeding(rhs.fNofFeeding), fNucleus(rhs.fNucleus)
{
   if(fNucleus == nullptr) {
      std::cout << "Copying level w/o a parent nucleus!" << std::endl;
   }
}

TLevel& TLevel::operator=(const TLevel& rhs)
{
   if(this != &rhs) {
      TObject::operator=(rhs);
      fEnergy            = rhs.fEnergy;
      fEnergyUncertainty = rhs.fEnergyUncertainty;
      fLabel             = rhs.fLabel;
      fTransitions       = rhs.fTransitions;
      fNofFeeding        = rhs.fNofFeeding;
      fNucleus           = rhs.fNucleus;
   }

   return *this;
}

TTransition* TLevel::AddTransition(const double levelEnergy, double br, double ts)
{
   /// MENU function to add gamma from this level to another level at "levelEnergy"
   /// Adds a new gamma ray with specified branching ratio (default 100.), and transition strength (default 1.).
   /// Returns gamma if it doesn't exist yet and was successfully added, null pointer otherwise.
   return AddTransition(levelEnergy, 0., br, ts);
}

TTransition* TLevel::AddTransition(const double levelEnergy, const double energyUncertainty, double br, double ts)
{
   /// Function to add gamma from this level to another level at "levelEnergy +- energyUncertainty".
   /// Adds a new gamma ray with specified branching ratio (default 100.), and transition strength (default 1.).
   /// Returns gamma if it doesn't exist yet and was successfully added, null pointer otherwise.

   // we expect the final level to have the same identifier as this level
   auto* level = fNucleus->FindLevel(levelEnergy, energyUncertainty, fIdentifier);
   if(level == nullptr) {
      std::cerr << DRED << "Failed to find level at " << levelEnergy << " keV, can't add transition!" << RESET_COLOR << std::endl;
      return nullptr;
   }

   // now that we found a level, we will use its energy instead of the levelEnergy parameter
   if(fTransitions.count(level->Energy()) == 1) {
      if(fVerbosity > EVerbosity::kQuiet) {
         std::cout << "Already found gamma ray from ";
         Print();
         std::cout << " to " << levelEnergy << "/" << level->Energy() << " keV";
         level->Print();
      }
      return nullptr;
   }

   // create the transition and set its properties
   // energy of the transition is calculated from the energy of this level and the level it populates

   fTransitions.emplace(std::piecewise_construct, std::forward_as_tuple(level->Energy()), std::forward_as_tuple(this, Energy() - levelEnergy, energyUncertainty, br, std::numeric_limits<double>::quiet_NaN(), ts, std::numeric_limits<double>::quiet_NaN()));
   level->AddFeeding();
   //// re-calculate the branching ratios in % for all gammas
   //double sum = std::accumulate(fTransitions.begin(), fTransitions.end(), 0., [](double r, std::pair<const double, TTransition>& g) { r += g.second.Intensity(); return r; });
   //for(auto& [en, gamma] : fTransitions) {
   //   gamma.RelativeIntensity(gamma.Intensity() / sum);
   //   gamma.RelativeIntensityUncertainty(gamma.IntensityUncertainty() / sum);
   //}
   level->AddFeeding();

   if(fVerbosity > EVerbosity::kQuiet) { Print(); }

   return &(fTransitions[level->Energy()]);
}

TLevel* TLevel::AddTransition(TTransition* transition, bool quiet)
{
   /// Adds the provided transition to the map of transitions.
   // simple search for the level using this levels energy, the transition energy and the sum of the uncertainties
   // this will return the level that matches the most closely in energy, so it's probably not necessary to change the uncertainty to sqrt of sum of squares?
   TLevel* level = nullptr;
   // we expect the final level to have the same identifier as this level
   if(std::isnan(transition->FinalLevel())) {
      level = fNucleus->FindLevel(fEnergy - transition->Energy(), fEnergyUncertainty + transition->EnergyUncertainty(), fIdentifier);
      if(fVerbosity >= EVerbosity::kBasicFlow) { std::cout << "Finding level in \"" << fNucleus->GetName() << "\" with energy " << fEnergy - transition->Energy() << " +- " << fEnergyUncertainty + transition->EnergyUncertainty() << " returned " << level << std::endl; }
   } else {
      // not sure what the right uncertainty is here, if a level is given can we assume it's exact?
      // for now we just use the transition's energy uncertainty and the function should return the best matching one anyways
      level = fNucleus->FindLevel(transition->FinalLevel(), transition->EnergyUncertainty(), transition->FinalLevelIdentifier());   //(std::isalpha(transition->FinalLevelIdentifier()) != 0) ? transition->FinalLevelIdentifier() : fIdentifier);
      if(fVerbosity >= EVerbosity::kBasicFlow) { std::cout << "Finding final level in \"" << fNucleus->GetName() << "\" with energy " << transition->FinalLevel() << " +- " << transition->EnergyUncertainty() << " returned " << level << std::endl; }
   }
   if(level == nullptr) {
      if(!quiet) {
         std::cerr << DRED << "Failed to find level at " << fEnergy - transition->Energy() << " keV, can't add " << transition->Energy() << " keV transition!" << RESET_COLOR << std::endl;
      }
      return nullptr;
   }
   if(this == level) {
      std::cout << DYELLOW << "Warning, adding transition of " << transition->Energy() << " +- " << transition->EnergyUncertainty() << " keV to level at " << fEnergy << " +- " << fEnergyUncertainty << " keV feeds the same level? Should be at " << fEnergy - transition->Energy() << " +- " << fEnergyUncertainty + transition->EnergyUncertainty() << RESET_COLOR << std::endl;
      if(fNucleus != nullptr) { fNucleus->Print(); }
   }
   // check if we have this transition already
   // if there is no final level in the transition, based on the energy of this level and the transiton energy, otherwise based on the energy of the final level
   if((std::isnan(transition->FinalLevel()) && fTransitions.find(fEnergy - transition->Energy()) != fTransitions.end()) || (!std::isnan(transition->FinalLevel()) && fTransitions.find(transition->FinalLevel()) != fTransitions.end())) {
      std::cout << DYELLOW << "Warning, can't add new transition " << transition << " to level at " << Energy() << " keV in \"" << (fNucleus != nullptr ? fNucleus->GetName() : "N/A") << "\" because we already have a transition with " << transition->Energy() << " keV to level at " << (std::isnan(transition->FinalLevel()) ? fEnergy - transition->Energy() : transition->FinalLevel()) << " keV (final level of transition is " << transition->FinalLevel() << " keV)" << RESET_COLOR << std::endl;
      return nullptr;
   }
   level->AddFeeding();

   fTransitions.emplace(level->Energy(), *transition);

   return level;
}

std::pair<double, double> TLevel::GetMinMaxTransition() const
{
   if(fTransitions.empty()) { return std::make_pair(-1., -1.); }
   auto result = std::make_pair(fTransitions.begin()->second.Intensity(), fTransitions.begin()->second.Intensity());
   for(const auto& [energy, gamma] : fTransitions) {
      if(gamma.Intensity() < result.first) { result.first = gamma.Intensity(); }
      if(gamma.Intensity() > result.second) { result.second = gamma.Intensity(); }
   }

   return result;
}

void TLevel::Print(Option_t* opt) const
{
   /// Print level information to stdout (with an ending newline).
   /// If opt is "g" also print list of all gamma rays draining this level.
   std::cout << "Level \"" << fLabel << "\" (" << this << ") at " << fEnergy << " keV has " << fTransitions.size() << " draining gammas and " << fNofFeeding << " feeding gammas" << std::endl;
   if(fVerbosity > EVerbosity::kQuiet || strcmp(opt, "g") == 0) {
      for(const auto& [level, gamma] : fTransitions) {
         std::cout << gamma.Energy() << " keV gamma to level at " << level << " keV" << std::endl;
      }
   }
}
