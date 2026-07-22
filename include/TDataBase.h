#ifndef TDATABASE_H
#define TDATABASE_H

#include <string>
#include <map>
#include <algorithm>

#include "Globals.h"
#include "TNucleus.h"

class TDataBase {
public:
   enum class EUncertainty : std::uint8_t { kLessThan, kGreaterThan, kLessOrEqual, kGreaterOrEqual, kApproximate, kCalculated, kSystematic, kDefault };

   explicit TDataBase(const std::string& path);

   bool LoadEnsdf(const std::string& path);
   bool LoadAbundance(const std::string& path);
   bool LoadMasses(const std::string& path);

   void ListNuclei(bool print = false) const;

   TNucleus* Nucleus(std::string label) const { std::transform(label.begin(), label.end(), label.begin(), ::toupper); return fNuclei.at(label); }

   static void Verbosity(EVerbosity val) { fVerbosity = val; }
   static EVerbosity Verbosity() { return fVerbosity; }

private:
   bool ParseEnsdfFile(const std::string& fileName);
   bool ReadEnergy(std::istringstream& str, double& energy, char& identifier);
   bool ReadUncertainty(std::istringstream& str, double& uncertainty, EUncertainty& uncertaintyLabel);

   static EVerbosity fVerbosity;
   std::map<std::string, TNucleus*> fNuclei;
};

#endif
