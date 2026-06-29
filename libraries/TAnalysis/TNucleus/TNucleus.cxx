// g++ -c -fPIC Nucleus.cc -I./ `root-config --cflags`

#include "TNucleus.h"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>

static double amu = 931.494043;

EVerbosity TNucleus::fVerbosity = EVerbosity::kQuiet;

std::string& TNucleus::Massfile()
{
   static std::string output = std::string(getenv("GRSISYS")) + "/libraries/TAnalysis/SourceData/mass.dat";
   return output;
}

TNucleus::TNucleus(const char* name, bool loadTransitions, bool quiet)
{
   /// Creates a nucleus based on symbol (ex. 26Na OR Na26) and sets all parameters from mass.dat
   int           number = 0;
   std::string   symbol;
   std::string   element;
   int           z = 0;
   int           n = 0;
   std::string   sym_name;
   double        mass  = 0.;
   bool          found = false;
   std::string   line;
   std::ifstream infile;
   std::string   massFile;
   ParseName(name, symbol, number, element);
   SetSymbol(symbol.c_str());
   z = GetZFromSymbol(symbol);
   n = number - z;
   SetZ(z);
   SetN(n);
   SetName();
   try {
      massFile = Massfile();
      infile.open(massFile.c_str());
      while(!getline(infile, line).fail()) {
         if(line.empty()) {
            continue;
         }
         std::istringstream ss(line);
         ss >> n;
         ss >> z;
         ss >> sym_name;
         ss >> mass;
         if(strcasecmp(element.c_str(), sym_name.c_str()) == 0) {
            found = true;
            break;
         }
      }
   } catch(std::out_of_range&) {
      std::cout << "Could not parse element " << name << std::endl
                << "Nucleus not set!" << std::endl;
      return;
   }
   if(fVerbosity >= EVerbosity::kBasicFlow) {
      std::cout << "After we " << (found ? "found":"didn't find") << " symbol \"" << symbol << "\"/\"" << GetSymbol() << "\", we got Z = " << z << ", N = " << n << ", and set name to \"" << GetName() << "\" (A = " << GetA() << ")" << std::endl;
   }
   if(!found) {
      if(!quiet) {
         std::cout << "Warning: Element " << element << " not found in the mass table " << massFile << ", Z, N, and mass excess not set!" << std::endl;
      }
      return;
   }
   infile.close();
   SetZ(z);
   SetN(n);
   SetMassExcess(mass / 1000.0);
   SetMass();
   if(loadTransitions) {
      LoadTransitionFile();
   }
}

TNucleus::TNucleus(int charge, int neutrons, double mass, const char* symbol)
   : fN(neutrons), fZ(charge), fMass(mass), fSymbol(symbol)
{
   /// Creates a nucleus with Z, N, mass, and symbol
   SetName();
   LoadTransitionFile();
}

TNucleus::TNucleus(int charge, int neutrons, const char* MassFile)
   : fN(neutrons), fZ(charge)
{
   /// Creates a nucleus with Z, N using mass table (default MassFile = "mass.dat")
   if(MassFile == nullptr) {
      MassFile = Massfile().c_str();
   }
   int           i     = 0;
   int           n     = 0;
   int           z     = 0;
   double        emass = 0.;
   std::string   tmp;
   std::ifstream mass_file;
   mass_file.open(MassFile, std::ios::in);
   while(!mass_file.bad() && !mass_file.eof() && i < 3008) {
      mass_file >> n;
      mass_file >> z;
      mass_file >> tmp;
      mass_file >> emass;
      if(n == fN && z == fZ) {
         fMassExcess = emass / 1000.;
         fSymbol     = tmp;
#ifdef debug
         std::cout << "Symbol " << fSymbol << " tmp " << tmp << std::endl;
#endif
         SetMass();
         SetSymbol(fSymbol.c_str());
         break;
      }
      i++;
      mass_file.ignore(256, '\n');
   }

   mass_file.close();
   std::string name   = fSymbol;
   std::string number = name.substr(0, name.find_first_not_of("0123456789 "));

   name = name.substr(name.find_first_not_of("0123456789 "));
   name.append(number);

   SetName();
   LoadTransitionFile();
}

void TNucleus::SetName(const char*)
{
   std::string name = GetSymbol();
   name.append(std::to_string(GetA()));
   TNamed::SetName(name.c_str());
}

TNucleus::~TNucleus()
{
   fTransitionListByIntensity.Delete();
   fTransitionListByEnergy.Delete();
}

void TNucleus::ParseName(std::string name, std::string& symbol, int& number, std::string& element)
{
   /// Strips any non-alphanumeric character from input, parses rest as number and symbol.
   /// E.g. turns "na26" into "Na" and 26 or "  152  EU... " into "Eu" and 152.
   /// Only uses the first number and first letter, so "na26   eu152" would be turned into "Na" and 26,
   /// but "na26  152eu" would be turned into "Na" and 26152, and "26na  152eu" would be turned into "Na152eu" and 26.
   /// Special inputs are "p" for "H" and 1, "d" for "H" and 2, "t" for "H" and 3, and "a" for "He" and 4.
   /// element simply is the combined string of number and element.

   // remove any characters that aren't alphanumeric
   name.erase(std::remove_if(name.begin(), name.end(), [](char c) { return std::isalnum(c) == 0; }), name.end());

   // special single letter cases
   if(name.length() == 1) {
      switch(name[0]) {
      case 'p':
         symbol.assign("H");
         return;
      case 'd':
         symbol.assign("H");
         return;
      case 't':
         symbol.assign("H");
         return;
      case 'a':
         symbol.assign("He");
         return;
      default:
         std::cout << "error, type numbersymbol, or symbolnumber, i.e. 30Mg oder Mg30, not " << name << std::endl;
         return;
      };
   }

   size_t firstDigit  = name.find_first_of("0123456789");
   size_t firstLetter = name.find_first_not_of("0123456789");
   if(firstDigit == std::string::npos || firstLetter == std::string::npos) {
      std::cout << "Failed to find either first digit (" << firstDigit << ") or first letter (" << firstLetter << ") in name \"" << name << "\", using name as symbol!" << std::endl;
      element = name;
      return;
   }
   if(firstDigit > firstLetter) {
      number = atoi(name.substr(firstDigit).c_str());
      symbol.append(name.substr(firstLetter, firstDigit - firstLetter));
   } else {
      number = atoi(name.substr(firstDigit, firstLetter - firstDigit).c_str());
      symbol.append(name.substr(firstLetter));
   }
   // make certain the symbol starts with upper case and rest is lower case by first transforming everything to lower case and then the first character to upper case.
   std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::tolower);
   symbol[0] = toupper(symbol[0]);
   element.append(std::to_string(number));
   element.append(symbol);

   if(fVerbosity >= EVerbosity::kBasicFlow) {
      std::cout << "parsed \"" << name << "\" to symbol \"" << symbol << "\", element \"" << element << "\", number " << number << std::endl;
   }
}

std::string TNucleus::SortName(std::string input)
{
   /// Strips any non-alphanumeric character from input, parses rest as number and symbol.
   /// E.g. turns "na26" into "26Na" or "  152  EU... " into "152Eu".
   /// Special inputs are "p" for "1H", "d" for "2H", "t" for "2H", and "a" for "4He".
   int         number = 0;
   std::string symbol;
   std::string element;
   ParseName(std::move(input), symbol, number, element);

   return element;
}

void TNucleus::SetZ(int charge)
{
   // Sets the Z (# of protons) of the nucleus
   fZ = charge;
}
void TNucleus::SetN(int neutrons)
{
   // Sets the N (# of neutrons) of the nucleus
   fN = neutrons;
}
void TNucleus::SetMassExcess(double mass_ex)
{
   // Sets the mass excess of the nucleus (in MeV)
   fMassExcess = mass_ex;
}

void TNucleus::SetMass(double mass)
{
   // Sets the mass manually (in MeV)
   fMass = mass;
}

void TNucleus::SetMass()
{
   // Sets the mass based on the A and mass excess of nucleus (in MeV)
   fMass = amu * GetA() + GetMassExcess();
}

void TNucleus::SetSymbol(const char* symbol)
{
   // Sets the atomic symbol for the nucleus
   fSymbol = symbol;
}

int TNucleus::GetZfromSymbol(char* symbol)
{
   // Figures out the Z of the nucleus based on the atomic symbol
   std::array<std::array<char, 3>, 118> symbols = {{{"H"}, {"HE"}, {"LI"}, {"BE"}, {"B"}, {"C"}, {"N"}, {"O"}, {"F"}, {"NE"}, {"NA"}, {"MG"}, {"AL"}, {"SI"}, {"P"}, {"S"}, {"CL"}, {"AR"}, {"K"}, {"CA"}, {"SC"}, {"TI"}, {"V"}, {"CR"}, {"MN"}, {"FE"}, {"CO"}, {"NI"}, {"CU"}, {"ZN"}, {"GA"}, {"GE"}, {"AS"}, {"SE"}, {"BR"}, {"KR"}, {"RB"}, {"SR"}, {"Y"}, {"ZR"}, {"NB"}, {"MO"}, {"TC"}, {"RU"}, {"RH"}, {"PD"}, {"AG"}, {"CD"}, {"IN"}, {"SN"}, {"SB"}, {"TE"}, {"F"}, {"XE"}, {"CS"}, {"BA"}, {"LA"}, {"CE"}, {"PR"}, {"ND"}, {"PM"}, {"SM"}, {"EU"}, {"GD"}, {"TB"}, {"DY"}, {"HO"}, {"ER"}, {"TM"}, {"YB"}, {"LU"}, {"HF"}, {"TA"}, {"W"}, {"RE"}, {"OS"}, {"IR"}, {"PT"}, {"AU"}, {"HG"}, {"TI"}, {"PB"}, {"BI"}, {"PO"}, {"AT"}, {"RN"}, {"FR"}, {"RA"}, {"AC"}, {"TH"}, {"PA"}, {"U"}, {"NP"}, {"PU"}, {"AM"}, {"CM"}, {"BK"}, {"CF"}, {"ES"}, {"FM"}, {"MD"}, {"NO"}, {"LR"}, {"RF"}, {"DB"}, {"SG"}, {"BH"}, {"HS"}, {"MT"}, {"DS"}, {"RG"}, {"CN"}, {"NH"}, {"FL"}, {"MC"}, {"LV"}, {"TS"}, {"OG"}}};
   size_t                               length  = strlen(symbol);
   auto*                                search  = new char[length + 1];
   for(size_t i = 0; i < length; i++) {
      search[i] = toupper(symbol[i]);   // make sure symbol is in uppercase
   }
   search[length] = '\0';
   for(int i = 0; i < 105; i++) {
      if(strcmp(search, symbols[i].data()) == 0) {
         delete[] search;
         SetZ(i + 1);
         return i + 1;
      }
   }

   delete[] search;
   SetZ(0);
   return 0;
}

int TNucleus::GetZFromSymbol(std::string symbol)
{
   /// Return the Z of the nucleus based on the atomic symbol
   std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::toupper);
   std::array<std::string, 118> symbols = {{{"H"}, {"HE"}, {"LI"}, {"BE"}, {"B"}, {"C"}, {"N"}, {"O"}, {"F"}, {"NE"}, {"NA"}, {"MG"}, {"AL"}, {"SI"}, {"P"}, {"S"}, {"CL"}, {"AR"}, {"K"}, {"CA"}, {"SC"}, {"TI"}, {"V"}, {"CR"}, {"MN"}, {"FE"}, {"CO"}, {"NI"}, {"CU"}, {"ZN"}, {"GA"}, {"GE"}, {"AS"}, {"SE"}, {"BR"}, {"KR"}, {"RB"}, {"SR"}, {"Y"}, {"ZR"}, {"NB"}, {"MO"}, {"TC"}, {"RU"}, {"RH"}, {"PD"}, {"AG"}, {"CD"}, {"IN"}, {"SN"}, {"SB"}, {"TE"}, {"F"}, {"XE"}, {"CS"}, {"BA"}, {"LA"}, {"CE"}, {"PR"}, {"ND"}, {"PM"}, {"SM"}, {"EU"}, {"GD"}, {"TB"}, {"DY"}, {"HO"}, {"ER"}, {"TM"}, {"YB"}, {"LU"}, {"HF"}, {"TA"}, {"W"}, {"RE"}, {"OS"}, {"IR"}, {"PT"}, {"AU"}, {"HG"}, {"TI"}, {"PB"}, {"BI"}, {"PO"}, {"AT"}, {"RN"}, {"FR"}, {"RA"}, {"AC"}, {"TH"}, {"PA"}, {"U"}, {"NP"}, {"PU"}, {"AM"}, {"CM"}, {"BK"}, {"CF"}, {"ES"}, {"FM"}, {"MD"}, {"NO"}, {"LR"}, {"RF"}, {"DB"}, {"SG"}, {"BH"}, {"HS"}, {"MT"}, {"DS"}, {"RG"}, {"CN"}, {"NH"}, {"FL"}, {"MC"}, {"LV"}, {"TS"}, {"OG"}}};

   auto* it = std::find(symbols.begin(), symbols.end(), symbol);

   if(it != symbols.end()) {
      return std::distance(symbols.begin(), it) + 1;
   }

   return 0;
}

double TNucleus::GetRadius() const
{
   // Gets the radius of the nucleus (in fm).
   // The radius is calculated using 1.12*A^1/3 - 0.94*A^-1/3
   return 1.12 * pow(GetA(), 1. / 3.) - 0.94 * pow(GetA(), -1. / 3.);
}

void TNucleus::AddTransition(Double_t energy, Double_t intensity, Double_t energy_uncertainty,
                             Double_t intensity_uncertainty)
{
   auto* tran = new TTransition();
   tran->SetEnergy(energy);
   tran->SetEnergyUncertainty(energy_uncertainty);
   tran->SetIntensity(intensity);
   tran->SetIntensityUncertainty(intensity_uncertainty);

   AddTransition(tran);
}

void TNucleus::AddTransition(TTransition* tran)
{
   fTransitionListByIntensity.Add(tran);
   auto* tran2 = new TTransition(*tran);
   tran2->SetCompareIntensity(false);
   fTransitionListByEnergy.Add(tran2);
}

TTransition* TNucleus::GetTransitionByIntensity(Int_t idx)
{
   auto* tran = static_cast<TTransition*>(fTransitionListByIntensity.At(idx));
   if(tran == nullptr) {
      std::cout << "Out of Range" << std::endl;
   }

   return tran;
}

TTransition* TNucleus::GetTransitionByEnergy(Int_t idx)
{
   auto* tran = static_cast<TTransition*>(fTransitionListByEnergy.At(idx));
   if(tran == nullptr) {
      std::cout << "Out of Range" << std::endl;
   }

   return tran;
}

void TNucleus::Print(Option_t*) const
{
   // Prints out the Name of the nucleus, as well as the numerated transition list
   std::cout << "Nucleus: " << GetName() << std::endl;
   std::cout << "\tA = " << GetA() << ", N = " << fN << ", Z = " << fZ << std::endl;
   std::cout << "\tMass = " << fMass << " MeV, mass excess = " << fMassExcess << " MeV" << std::endl;
   if(!std::isnan(fQValue)) {
      std::cout << "\tQ value = " << fQValue << " +- " << fQValueUncertainty << std::endl;
   }
   if(!std::isnan(fAlphaQValue)) {
      std::cout << "\tAlpha Q value = " << fAlphaQValue << " +- " << fAlphaQValueUncertainty << std::endl;
   }
   if(!std::isnan(fNeutronSeparation)) {
      std::cout << "\tS_n = " << fNeutronSeparation << " +- " << fNeutronSeparationUncertainty << std::endl;
   }
   if(!std::isnan(fProtonSeparation)) {
      std::cout << "\tS_p = " << fProtonSeparation << " +- " << fProtonSeparationUncertainty << std::endl;
   }
   if(!fLevels.empty()) {
      std::cout << fLevels.size() << " levels:" << std::endl;
      for(const auto& [energy, level] : fLevels) {
         std::cout << energy << " keV: ";
         if(level.size() == 1) {
            level[0].Print("g");
         } else {
            std::cout << level.size() << " levels at this energy" << std::endl;
            for(const auto& l : level) {
               l.Print("g");
            }
         }
      }
   }
   TIter next(&fTransitionListByIntensity);
   int   counter = 0;
   while(auto* tran = static_cast<TTransition*>(next())) {
      std::cout << "\t" << counter++ << "\t";
      tran->Print();
   }
}

void TNucleus::WriteSourceFile(const std::string& outfilename)
{
   if(!outfilename.empty()) {
      std::ofstream sourceout;
      sourceout.open(outfilename.c_str());
      for(int i = 0; i < fTransitionListByIntensity.GetSize(); i++) {
         std::string transtr = (static_cast<TTransition*>(fTransitionListByIntensity.At(i)))->PrintToString();
         sourceout << transtr.c_str();
         sourceout << std::endl;
      }
      sourceout << std::endl;
      sourceout.close();
   }
}

bool TNucleus::LoadTransitionFile()
{
   if(fTransitionListByIntensity.GetSize() != 0) {
      return false;
   }
   std::string filename;
   filename           = std::string(getenv("GRSISYS")) + "/libraries/TAnalysis/SourceData/";
   std::string symbol = GetSymbol();
   std::transform(symbol.begin(), symbol.end(), symbol.begin(), ::tolower);
   filename.append(symbol);
   filename.append(std::to_string(GetA()));
   filename.append(".sou");
   std::ifstream transfile;
   transfile.open(filename.c_str());
   if(!transfile.is_open()) {
      std::cout << "failed to open source file: " << filename.c_str() << std::endl;
      return false;
   }

   std::string line;

   while(!getline(transfile, line).fail()) {
      if(line.compare(0, 2, "//") == 0) {
         continue;
      }
      if(line.compare(0, 1, "#") == 0) {
         continue;
      }
      double             temp = 0.;
      auto*              tran = new TTransition;
      std::istringstream str(line);
      int                counter = 0;
      while(!(str >> temp).fail()) {
         counter++;
         if(counter == 1) {
            tran->SetEnergy(temp);
         } else if(counter == 2) {
            tran->SetEnergyUncertainty(temp);
         } else if(counter == 3) {
            tran->SetIntensity(temp);
         } else if(counter == 4) {
            tran->SetIntensityUncertainty(temp);
         } else {
            break;
         }
      }
      AddTransition(tran);
   }

   fTransitionListByIntensity.Sort();

   return true;
}

double TNucleus::GetEnergyFromBeta(double beta) const
{
   double gamma = 1 / std::sqrt(1 - beta * beta);
   return fMass * (gamma - 1);
}

double TNucleus::GetBetaFromEnergy(double energy_MeV) const
{
   double gamma = energy_MeV / fMass + 1;
   double beta  = std::sqrt(1 - 1 / (gamma * gamma));
   return beta;
}

TLevel* TNucleus::AddLevel(Double_t energy, Double_t energyUncertainty)
{
   /// Add a new level at provided energy. Since we have a vector of levels for each energy, we can have multiple levels with the same energy.
   /// Returns the newly created level.
   // The [] operator either returns the existing vector for that key, or (if that key doesn't exist yet), creates a new key-vector pair and returns the empty vector.
   // Either way, we simply add the given level to this key.
   fLevels[energy].emplace_back(this, energy, energyUncertainty);
   return &(fLevels[energy].back());
}

TLevel* TNucleus::FindLevel(Double_t levelEnergy, Double_t energyUncertainty, int index)
{
   /// Returns level at the provided energy +- the uncertainty. If there isn't any level in that range a null pointer is returned.
   /// If there are multiple levels in the range the one with the smallest energy difference is returned.
   auto low  = fLevels.lower_bound(levelEnergy - energyUncertainty);
   auto high = fLevels.upper_bound(levelEnergy + energyUncertainty);
   if(low == high) {
      if(fVerbosity >= EVerbosity::kBasicFlow) {
         std::cout << this << ": failed to find level with " << levelEnergy << " +- " << energyUncertainty << " keV" << std::endl;
         Print();
      }
      return nullptr;
   }
   // if we found multiple matching levels, use the one with the smallest energy difference
   // if there is only one matching level, the loop won't executre and level is already set correctly
   double en    = low->first;
   auto   level = low->second;
   for(auto& it = ++low; it != high; ++it) {
      if(std::fabs(levelEnergy - en) > std::fabs(levelEnergy - it->first)) {
         en    = it->first;
         level = it->second;
      }
   }
   if(index < 0) {
      if(level.size() > 1) {
         std::cout << "Warning, found " << level.size() << " levels with energy " << en << " keV (best match for requested energy " << levelEnergy << " keV), but no index has been provided, going to return the first one!" << std::endl;
      }
      index = 0;
   }
   return &((low->second)[index]);
}

