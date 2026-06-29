#ifndef TNUCLEUS_H
#define TNUCLEUS_H

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>

#include "Globals.h"
#include "TTransition.h"
#include "TLevel.h"

#include "TNamed.h"
#include "TSortedList.h"

/////////////////////////////////////////////////////////////////
///
/// \class TNucleus
///
/// This class builds a nucleus and sets all the basic information
/// (mass, Z, symbol, radius, etc.)
///
/////////////////////////////////////////////////////////////////

class TNucleus : public TNamed {
public:
   TNucleus() = default;                                                                   ///< Should not be used, here so we can write things to a root file.
   explicit TNucleus(const char* name, bool loadTransitions = true, bool quiet = false);   ///< Creates a nucleus based on symbol and sets all parameters from mass.dat
   TNucleus(int charge, int neutrons, double mass, const char* symbol);                    ///< Creates a nucleus with Z, N, mass, and symbol
   TNucleus(int charge, int neutrons, const char* MassFile = nullptr);                     ///< Creates a nucleus with Z, N using mass table (default MassFile = "mass.dat")
   TNucleus(const TNucleus&)                = delete;
   TNucleus(TNucleus&&) noexcept            = delete;
   TNucleus& operator=(const TNucleus&)     = delete;
   TNucleus& operator=(TNucleus&&) noexcept = delete;

   ~TNucleus();

   enum class EFlag : uint8_t { kDefault,
                                kUnobserved,
                                kInferred,
                                kTentative,
                                kObserved };

   void Observed(bool val)
   {
      if(val) {
         fFlag = EFlag::kObserved;
      } else {
         fFlag = EFlag::kDefault;
      }
   }
   void Unobserved(bool val)
   {
      if(val) {
         fFlag = EFlag::kUnobserved;
      } else {
         fFlag = EFlag::kDefault;
      }
   }
   void Inferred(bool val)
   {
      if(val) {
         fFlag = EFlag::kInferred;
      } else {
         fFlag = EFlag::kDefault;
      }
   }
   void Tentative(bool val)
   {
      if(val) {
         fFlag = EFlag::kTentative;
      } else {
         fFlag = EFlag::kDefault;
      }
   }

   bool Observed() { return fFlag == EFlag::kObserved; }
   bool Unobserved() { return fFlag == EFlag::kUnobserved; }
   bool Inferred() { return fFlag == EFlag::kInferred; }
   bool Tentative() { return fFlag == EFlag::kTentative; }

   static void ParseName(const char* name, std::string& symbol, int& number, std::string& element)
   {
      ParseName(std::string(name), symbol, number, element);
   }
   static void        ParseName(std::string name, std::string& symbol, int& number, std::string& element);
   static std::string SortName(const char* input)
   {
      std::string tmp(input);
      return SortName(tmp);
   }
   static std::string SortName(std::string input);
   void               SetZ(int);                ///< Sets the Z (# of protons) of the nucleus
   void               SetN(int);                ///< Sets the N (# of neutrons) of the nucleus
   void               SetMassExcess(double);    ///< Sets the mass excess of the nucleus (in MeV)
   void               SetMass(double);          ///< Sets the mass manually (in MeV)
   void               SetMass();                ///< Sets the mass based on the A and mass excess of nucleus (in MeV)
   void               SetSymbol(const char*);   ///< Sets the atomic symbol for the nucleus

   TLevel* AddLevel(Double_t energy, Double_t energyUncertainty);
   TLevel* FindLevel(Double_t levelEnergy, Double_t energyUncertainty, int index = -1);

   void AddTransition(Double_t energy, Double_t intensity, Double_t energy_uncertainty = 0.0, Double_t intensity_uncertainty = 0.0);
   void AddTransition(TTransition* tran);

   int         GetZ() const { return fZ; }                     ///< Gets the Z (# of protons) of the nucleus
   int         GetN() const { return fN; }                     ///< Gets the N (# of neutrons) of the nucleus
   int         GetA() const { return fN + fZ; }                ///< Gets the A (Z + N) of the nucleus
   double      GetMassExcess() const { return fMassExcess; }   ///< Gets the mass excess of the nucleus (in MeV)
   double      GetMass() const { return fMass; }               ///< Gets the mass of the nucleus (in MeV)
   const char* GetSymbol() const { return fSymbol.c_str(); }   ///< Gets the atomic symbol of the nucleus

   void QValue(const double& val) { fQValue = val; }
   void QValueUncertainty(const double& val) { fQValueUncertainty = val; }
   void AlphaQValue(const double& val) { fAlphaQValue = val; }
   void AlphaQValueUncertainty(const double& val) { fAlphaQValueUncertainty = val; }
   void NeutronSeparation(const double& val) { fNeutronSeparation = val; }
   void NeutronSeparationUncertainty(const double& val) { fNeutronSeparationUncertainty = val; }
   void ProtonSeparation(const double& val) { fProtonSeparation = val; }
   void ProtonSeparationUncertainty(const double& val) { fProtonSeparationUncertainty = val; }

   double QValue() const { return fQValue; }
   double QValueUncertainty() const { return fQValueUncertainty; }
   double AlphaQValue() const { return fAlphaQValue; }
   double AlphaQValueUncertainty() const { return fAlphaQValueUncertainty; }
   double NeutronSeparation() const { return fNeutronSeparation; }
   double NeutronSeparationUncertainty() const { return fNeutronSeparationUncertainty; }
   double ProtonSeparation() const { return fProtonSeparation; }
   double ProtonSeparationUncertainty() const { return fProtonSeparationUncertainty; }

   // Returns total kinetic energy in MeV
   double GetEnergyFromBeta(double beta) const;
   double GetBetaFromEnergy(double energy_MeV) const;

   TTransition* GetTransition(Int_t idx) { return GetTransitionByIntensity(idx); }
   TTransition* GetTransitionByIntensity(Int_t idx);
   TTransition* GetTransitionByEnergy(Int_t idx);

   Int_t  NTransitions() const { return fTransitionListByIntensity.GetSize(); };
   Int_t  GetNTransitions() const { return fTransitionListByIntensity.GetSize(); };
   double GetRadius() const;
   int    GetZfromSymbol(char* symbol);
   int    GetZFromSymbol(std::string symbol);

   void Print(Option_t* opt = "") const override;
   void WriteSourceFile(const std::string& outfilename = "");

   const TSortedList* GetTransitionList() const { return GetTransitionListByIntensity(); }
   const TSortedList* GetTransitionListByIntensity() const { return &fTransitionListByIntensity; }
   const TSortedList* GetTransitionListByEnergy() const { return &fTransitionListByEnergy; }

   bool operator==(const TNucleus& rhs) const { return ((fN == rhs.fN) && (fZ == rhs.fZ)); }
   bool operator!=(const TNucleus& rhs) const { return !(*this == rhs); }

   static std::string& SourceDirectory();   /// < Returns the directory with the .sou files and the mass file

   static void       Verbosity(EVerbosity val) { fVerbosity = val; }
   static EVerbosity Verbosity() { return fVerbosity; }

private:
   bool LoadTransitionFile();

   static EVerbosity fVerbosity;   ///< verbosity of the class

   static std::string& Massfile();   ///< Returns the massfile to be used, which includes Z, N, atomic symbol, and mass excess
   void                SetName(const char* name = "") override;

   static std::string fSourceDirectory;          //!<! path of directory with .sou files
   static bool        fSourceDirectoryChecked;   //!<! flag to indicate whetehr the source directory path has been checked

   int         fN{0};                    ///< Number of neutrons (N)
   int         fZ{0};                    ///< Number of protons (Z)
   double      fMass{0.};                ///< Mass (in MeV)
   double      fMassExcess{0.};          ///< Mass excess (in MeV)
   std::string fSymbol;                  ///< Atomic symbol (ex. Ba, C, O, N)
   EFlag       fFlag{EFlag::kDefault};   ///< Flag indicating if nucleus is observed, unobserved, inferred, or tentative
   double      fQValue{std::numeric_limits<double>::quiet_NaN()};
   double      fQValueUncertainty{std::numeric_limits<double>::quiet_NaN()};
   double      fAlphaQValue{std::numeric_limits<double>::quiet_NaN()};
   double      fAlphaQValueUncertainty{std::numeric_limits<double>::quiet_NaN()};
   double      fNeutronSeparation{std::numeric_limits<double>::quiet_NaN()};
   double      fNeutronSeparationUncertainty{std::numeric_limits<double>::quiet_NaN()};
   double      fProtonSeparation{std::numeric_limits<double>::quiet_NaN()};
   double      fProtonSeparationUncertainty{std::numeric_limits<double>::quiet_NaN()};

   std::map<double, std::vector<TLevel>> fLevels;

   TSortedList fTransitionListByIntensity;
   TSortedList fTransitionListByEnergy;

   /// \cond CLASSIMP
   ClassDefOverride(TNucleus, 2)   // NOLINT(readability-else-after-return)
   /// \endcond
};

#endif
