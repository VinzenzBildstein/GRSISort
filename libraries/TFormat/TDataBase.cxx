#include "TDataBase.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "Globals.h"
#include "TGRSIUtilities.h"
#include "TTransition.h"

EVerbosity TDataBase::fVerbosity = EVerbosity::kDefault;

TDataBase::TDataBase(const std::string& path)
{
   // we need the ensdf information later, so if we fail to load it here we can't proceed
   if(!LoadEnsdf(path)) {
      std::ostringstream str;
      str << "Failed to load ENSDF data base from path \"" << path << "\", please check earlier error message for more details" << std::endl;
      throw std::runtime_error(str.str().c_str());
   }
   if(fVerbosity > EVerbosity::kQuiet) {
      std::cout << "Successfully loaded ENSDF from path \"" << path << "\"" << std::endl;
   }
   // these could be just error messages instead of exceptions if there are use cases with only the ensdf data being of interest
   if(!LoadAbundance(path)) {
      std::ostringstream str;
      str << "Failed to load ENSDF data base from path \"" << path << "\", please check earlier error message for more details" << std::endl;
      throw std::runtime_error(str.str().c_str());
   }
   if(fVerbosity > EVerbosity::kQuiet) {
      std::cout << "Successfully loaded abundance from path \"" << path << "\"" << std::endl;
   }
   if(!LoadMasses(path)) {
      std::ostringstream str;
      str << "Failed to load ENSDF data base from path \"" << path << "\", please check earlier error message for more details" << std::endl;
      throw std::runtime_error(str.str().c_str());
   }
   if(fVerbosity > EVerbosity::kQuiet) {
      std::cout << "Successfully loaded masses from path \"" << path << "\"" << std::endl;
   }
}

bool TDataBase::LoadEnsdf(const std::string& path)
{
   int                fileIndex = 1;
   std::ostringstream fileName;
   fileName << path << "/ensdf." << std::setw(3) << std::setfill('0') << fileIndex;

   while(FileExists(fileName.str().c_str())) {
      if(fVerbosity >= EVerbosity::kLoops) {
         std::cout << "Trying to read file \"" << fileName.str() << "\"" << std::endl;
      }

      if(!ParseEnsdfFile(fileName.str())) {
         return false;
      }
      ++fileIndex;
      fileName.clear();
      fileName.str("");
      fileName << path << "/ensdf." << std::setw(3) << std::setfill('0') << fileIndex;
   }
   if(fVerbosity > EVerbosity::kQuiet) {
      std::cout << "Successfully read files up to file \"" << fileName.str() << "\"" << std::endl;
   }

   return true;
}

bool TDataBase::ReadEnergy(std::istringstream& str, double& energy, char& identifier)
{
   /// Reads an energy field as defined in V.18 of the ensdf manual.
   /// The field can have on of four forms:
   /// 1. NUM - unsigned number, can be integer or real (including integer exponents like E+3)
   /// 2. NUM+A or A+NUM, where A is a single upper case character in the order X, Y, Z, U, V, W, A, B, ...
   /// 3. SN+NUM or SP+NUM for resonances
   /// 4. A with A as defined in 2.
   /// Parentheses are allowed for numbers that have been deduced or taken from other experiments, not sure how to deal with these yet.

   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Reading energy from str \"" << str.str() << "\", starting with '" << str.peek() << "' = '" << static_cast<char>(str.peek()) << "', isdigit " << std::isdigit(str.peek()) << std::endl;
   }
   // try and read the energy (this skips any leading whitespace)
   str >> energy;
   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Read energy " << energy << ", string status is " << (str.fail() ? "fail" : "not fail") << ", " << (str.good() ? "good" : "not good") << ", and " << (str.eof() ? "eof" : "not eof") << std::endl;
   }
   if(str.eof()) {   // reached the end of the string, which implies we read a number and there's nothing else but blank space (case 1.)
      str.clear();   // clear the eof flag
      if(fVerbosity >= EVerbosity::kSubroutines) {
         std::cout << "Case 1 (NUM)" << std::endl;
      }
   } else if(str.good()) {   // check if we succeeded in reading a number (case 2.) and still have something left to read
      str.clear();
      // there are a few instances where there is a space between the number and the '+' sign, so we skip all whitespace
      while(std::isspace(str.peek()) != 0) {
         str.get();
      }
      if(fVerbosity >= EVerbosity::kSubroutines) {
         std::cout << "after skipping whitespace, string status is " << (str.fail() ? "fail" : "not fail") << ", " << (str.good() ? "good" : "not good") << ", and " << (str.eof() ? "eof" : "not eof") << std::endl;
      }
      // check if we stopped because we reached the end
      if(str.eof()) {
         if(fVerbosity >= EVerbosity::kSubroutines) {
            std::cout << "stringstream reached end-of-file after skipping whitespace, returning" << std::endl;
         }
         return !str.fail();
      }
      // check if the number is followed by a '+' (case 2.)
      if(str.peek() == '+') {
         // discard plus sign and read the identifier
         str.get();
         str.get(identifier);
         if(fVerbosity >= EVerbosity::kSubroutines) {
            std::cout << "Case 2a (NUM+A)" << std::endl;
         }
      } else if(fVerbosity >= EVerbosity::kSubroutines) {
         // this is not a real problem, sometimes we just have a number that we read successfully but still do not reach eof?
         std::cout << "Read energy " << energy << " from string \"" << str.str() << "\", but we don't have blank space or a + sign next: '" << str.peek() << "', string status is " << (str.fail() ? "fail" : "not fail") << ", " << (str.good() ? "good" : "not good") << ", and " << (str.eof() ? "eof" : "not eof") << std::endl;
      }
   } else {
      // failed to read a number so clear fail bits and try the other cases that don't start with a number
      // do we need to skip leading whitespace?
      str.clear();
      if(str.str().substr(0, 3) == "SN+" || str.str().substr(0, 3) == "SP+") {   // check if there is a plus sign in third position (case 3.)
                                                                                 // not sure how to properly handle this, for now we set the identifier to either n or p (lowercase!)
         str.get();                                                              // discard the 'S'
         str.get(identifier);                                                    // read the 'N' or 'P'
         identifier = std::tolower(identifier);                                  // convert to lower case to distinguish it from a normal identifier
         str.get();                                                              // discard the '+'
         str >> energy;
         if(fVerbosity >= EVerbosity::kSubroutines) {
            std::cout << "Case 3 (SN+NUM or SP+NUM)" << std::endl;
         }
      } else {   // case 4.: read the identifier and set the energy to zero, could also be 2b (A+NUM)?
         str.get(identifier);
         if(str.peek() == '+') {
            str.get();
            str >> energy;
            if(fVerbosity >= EVerbosity::kSubroutines) {
               std::cout << "Case 2b (A+NUM), read identifier '" << identifier << "' and energy " << energy << std::endl;
            }
         } else {
            energy = 0.;
            if(fVerbosity >= EVerbosity::kSubroutines) {
               std::cout << "Case 4 (A)" << std::endl;
            }
         }
      }
   }

   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Read energy " << energy << ", identifier '" << identifier << "' (" << static_cast<int>(identifier) << "), returning " << (!str.fail() ? "true" : "false") << std::endl;
   }

   return !str.fail();
}

bool TDataBase::ReadUncertainty(std::istringstream& str, double& uncertainty, EUncertainty& uncertaintyLabel)
{
   // The uncertainty is a two character field that is either blank, an integer, or two letters:
   // LT - less than, GT - greater than, LE - less or equal, GE - greater or equal, AP - approximate, CA - calculated, or SY - systematics
   // so we try and read the number, if that fails, we check if we have the two letters
   str >> uncertainty;
   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Reading uncertainty from \"" << str.str() << "\", got uncertainty " << uncertainty << ", " << (str.fail() ? "fail" : "not fail") << ", and " << (str.eof() ? "eof" : "not eof") << std::endl;
   }
   if(str.fail() && !str.eof()) {
      // reset the uncertainty to NaN (the default value)
      uncertainty = std::numeric_limits<double>::quiet_NaN();
      // not sure what the best way is to check the two letters, so for now we just read them and do if/else?
      str.clear();
      auto first  = str.get();
      auto second = str.get();
      if(fVerbosity >= EVerbosity::kSubroutines) {
         std::cout << "Got first '" << first << "' and second '" << second << "'" << std::endl;
      }
      if(first == 'L') {
         if(second == 'T') {
            uncertaintyLabel = EUncertainty::kLessThan;
            return true;
         }
         if(second == 'E') {
            uncertaintyLabel = EUncertainty::kLessOrEqual;
            return true;
         }
         return false;
      }
      if(first == 'G') {
         if(second == 'T') {
            uncertaintyLabel = EUncertainty::kGreaterThan;
            return true;
         }
         if(second == 'E') {
            uncertaintyLabel = EUncertainty::kGreaterOrEqual;
            return true;
         }
         return false;
      }
      if(first == 'A' && second == 'P') {
         uncertaintyLabel = EUncertainty::kApproximate;
         return true;
      }
      if(first == 'C' && second == 'A') {
         uncertaintyLabel = EUncertainty::kCalculated;
         return true;
      }
      if(first == 'S' && second == 'Y') {
         uncertaintyLabel = EUncertainty::kSystematic;
         return true;
      }
      return false;
   }
   return true;
}

bool TDataBase::ParseEnsdfFile(const std::string& fileName)
{
   /// This function parses the given file assuming it's ENSDF formatted.
   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Reading ENSDF file \"" << fileName << "\"" << std::endl;
   }

   std::ifstream file(fileName);

   if(!file.is_open()) {
      std::cerr << DRED << "Failed to open \"" << fileName << "\"" << RESET_COLOR << std::endl;
      return false;
   }

   std::istringstream str;
   std::string        line;
   int                subSection = 0;
   //int firstQLine = 1;
   std::string tmpString;
   double      tmpDouble = 0.;
   //std::vector<std::string> values;
   TNucleus*    currentNucleus = nullptr;
   TLevel*      currentLevel   = nullptr;
   TTransition* currentGamma   = nullptr;
   int          lineNumber     = 0;
   bool         rejectedGamma  = false;

   while(file.good()) {
      std::getline(file, line);
      ++lineNumber;
      if(fVerbosity == EVerbosity::kAll) {
         std::cout << lineNumber << ": read line \"" << line << "\"" << std::endl;
      }

      if(line.empty()) {
         ++subSection;
         //firstQLine = 1;
         if(fVerbosity == EVerbosity::kAll) {
            std::cout << lineNumber << ": empty line, increased sub section to " << subSection << std::endl;
         }
         continue;
      }

      // general identifier format of line (first 8 characters)
      // 1-3 mass
      // 4-5 element symbol (except for reference record where it's blank)
      // 6 blank/'1' for primary, any ascii character (but not '1) for continuation
      // 7 blank or c, C, d, D, t, or T for comment
      // 8 R - reference, X - cross reference, H - history, Q - q-value, P - parent, N - normalization, L - level, B - beta, E - EC, A - alpha, D - delayed particle, G - gamma

      // first line should be identification record (1-5 nuclide, 10-39 data set ident., 40-65 refs., 66-74 publ. inf., 75-80 date
      auto label = line.substr(0, 5);
      trimWS(label);
      // check special label case
      if(label == "1NN") {
         label = "1n";
      }
      if(line.substr(9, 14) == "ADOPTED LEVELS") {
         if(fVerbosity >= EVerbosity::kLoops) {
            std::cout << lineNumber << ": found adopted levels for \"" << label << "\"" << std::endl;
         }
         // this is the beginning of the entry of this nuclide, so we re-set some values
         subSection = 0;
         // check if we already have this nucleus (we should never get here?)
         if(fNuclei.find(label) != fNuclei.end()) {
            std::cerr << RED << R"(Error, found line with "ADOPTED LEVELS" for nucleus ")" << label << "\" (\" " << line << "\"), but we already have this nucleus:" << std::endl;
            fNuclei[label]->Print();
            return false;
         }
         fNuclei[label] = new TNucleus(label.c_str(), false, true);   // create nucleus with given label, but do not load transition file (.sou file), and be quiet
         currentNucleus = fNuclei[label];
         if(fVerbosity >= EVerbosity::kQuiet) {
            std::cout << lineNumber << ": created new nucleus \"" << label << "\"" << std::endl;
         }
         // check flags for unobserved, inferred, or tentative nuclei
         auto flag = line.substr(23, 7);
         if(flag == ":NOT OB" || flag == ":UNOBSE") {
            currentNucleus->Unobserved(true);
         } else if(flag == ":INFERR") {
            currentNucleus->Inferred(true);
         } else if(flag == ":TENTAT") {
            currentNucleus->Tentative(true);
         } else {
            currentNucleus->Observed(true);
            // maybe check and record if this is the min/max N for this Z and vice versa?
         }

         // if the identification record is continued, 6 will not be blank
         // for now we simly skip these lines
         do {
            std::getline(file, line);
            ++lineNumber;
            if(fVerbosity == EVerbosity::kAll) {
               std::cout << lineNumber << ": skipping line \"" << line << "\"" << std::endl;
            }
         } while(file.good() && (line.length() < 6 || line[5] != ' '));

         // now we have the first line that is not part of the identification record, so we check it's format and continue reading lines as long as we can
         do {
            if(fVerbosity >= EVerbosity::kLoops) {
               std::cout << lineNumber << ": processing line \"" << line << "\"" << std::endl;
            }
            switch(line[7]) {
            case 'H':   // history record (1-5 nuclide, 6 blank or anything but '1' for cont., 7 blank, 8 'H', 9 blank, 10-80 history)
                        // ignored
               break;
            case 'Q':   // q-value record (1-5 nuclide, 6-9 "  Q ", 10-19 beta- q-value, 20-21 uncert., 22-29 S_n, 30-31 uncert., 32-39 S_p, 40-41 uncert., 42-49 alpha q-value, 50-51 uncert., 56-80 refs.)
               if(line.substr(5, 4) == "  Q ") {
                  // we only read the beta- q-value and the neutron separation energy
                  str.str(line.substr(9, 10));
                  str.clear();
                  str >> tmpDouble;
                  currentNucleus->QValue(tmpDouble);
                  str.str(line.substr(19, 2));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  str >> tmpDouble;
                  currentNucleus->QValueUncertainty(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Q-value: " << currentNucleus->QValue() << " +- " << currentNucleus->QValueUncertainty() << std::endl; }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading S_n old stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str.str(line.substr(21, 8));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->NeutronSeparation(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => S_n " << currentNucleus->NeutronSeparation() << std::endl; }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading old stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str.str(line.substr(29, 2));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->NeutronSeparationUncertainty(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => delta S_n " << currentNucleus->NeutronSeparationUncertainty() << std::endl; }
                  str.str(line.substr(31, 8));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->ProtonSeparation(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => S_p " << currentNucleus->ProtonSeparation() << std::endl; }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading old stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str.str(line.substr(39, 2));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->ProtonSeparationUncertainty(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => delta S_p " << currentNucleus->ProtonSeparationUncertainty() << std::endl; }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading old stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str.str(line.substr(41, 8));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->AlphaQValue(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => Q_alpha " << currentNucleus->AlphaQValue() << std::endl; }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading old stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str.str(line.substr(49, 6));
                  str.clear();
                  tmpDouble = std::numeric_limits<double>::quiet_NaN();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new stream \"" << str.str() << "\" @ " << str.tellg() << (str.fail() ? " failed" : " good"); }
                  str >> tmpDouble;
                  currentNucleus->AlphaQValueUncertainty(tmpDouble);
                  if(fVerbosity >= EVerbosity::kLoops) {
                     std::cout << std::endl
                               << "Found q-value line: " << currentNucleus->QValue() << " +- " << currentNucleus->QValueUncertainty() << ", S_n " << currentNucleus->NeutronSeparation() << " +- " << currentNucleus->NeutronSeparationUncertainty() << std::endl
                               << "\"" << line.substr(9, 10) << "\", \"" << line.substr(19, 2) << "\", \"" << line.substr(21, 8) << "\", \"" << line.substr(29, 2) << "\"" << std::endl;
                  }
               } else if(fVerbosity >= EVerbosity::kLoops) {   // if 7 is not blank this is a comment for the q-value and we ignore it
                  std::cout << "Ignoring q-value comment \"" << line << "\"" << std::endl;
               }
               break;
            case 'X':   // cross-reference record (1-5 nuclide, 6-7 blank, 8 'X', 9 identifier, 10-39 DSID used, 40-80 blank)
                        // ignored
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring cross-reference \"" << line << "\"" << std::endl; }
               break;
            case 'P':   // parent record (1-5 nuclide, 6-7 blank, 8 'P', 9 blank or integer, ...), seems to be for decay data sets only?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring parent \"" << line << "\"" << std::endl; }
               break;
            case 'N':   // this can be one of two records, depending on whether 7 is blank or 'P'
                        // normalization record (1-5 nuclide, 6-7 blank, 8 'N', 9 blank or integer, ...), seems to be for decay data sets only?
                        // production normalization record (1-5 nuclide, 6-9 " PN ", ...) ignored
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring normalization \"" << line << "\"" << std::endl; }
               break;
            case 'L':   // level record (1-5 nuclide, 6 blank or not '1' for cont., 7-9 " L ", 10-19 energy keV, 20-21 uncert., 22-39 spin and parity, 40-49 half life with units, 50-55 uncert.
                        // 56-74 angular momentum transfer in reaction, 75-76 uncert., 77 comment flag, 78-79 metastable as "M ", or "M1", "M2", etc., 80 '?' denotes uncertain level and 'S' denotes neutron, proton, alpha sep. en.)
               if(line[5] == ' ' && line[6] == ' ') {
                  if(currentGamma != nullptr) {
                     // we still have the last gamma for the previous level, so first add that to the level (or nucleus)
                     if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Adding gamma at " << currentGamma->Energy() << " keV to " << std::flush; }
                     if(currentLevel == nullptr) {
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "nucleus \"" << currentNucleus->GetName() << "\"" << std::endl;
                           currentNucleus->Print();
                        }
                        currentNucleus->AddTransition(currentGamma);
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "done:" << std::endl;
                           currentNucleus->Print();
                        }
                     } else {
                        // try to (quietly) add the transition to the current level, if it fails, add it to the nucleus
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "level at " << currentLevel->Energy() << " keV to " << std::endl;
                           currentLevel->Print();
                        }
                        if(currentLevel->AddTransition(currentGamma, true) == nullptr) {
                           currentNucleus->AddTransition(currentGamma);
                        }
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "done:" << std::endl;
                           currentLevel->Print();
                        }
                     }
                     currentGamma = nullptr;
                  }
                  // read energy
                  double energy     = 0.;
                  char   identifier = '\0';
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << "reading level energy old \"" << str.str() << "\""; }
                  str.str(line.substr(9, 10));
                  str.clear();
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << ", new \"" << str.str() << "\""; }
                  if(!ReadEnergy(str, energy, identifier)) {
                     if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Failed to read energy (got " << energy << " from \"" << str.str() << "\"/\"" << line.substr(9, 10) << "\"), skipping this line" << std::endl; }
                     break;
                  }
                  if(fVerbosity >= EVerbosity::kLoops) { std::cout << " => energy " << energy << std::endl; }
                  // read energy uncertainty
                  double       energyUncertainty = 0.;
                  EUncertainty uncertaintyLabel  = EUncertainty::kDefault;
                  str.str(line.substr(19, 2));
                  str.clear();
                  if(!ReadUncertainty(str, energyUncertainty, uncertaintyLabel)) {
                     std::cout << "Failed to read uncertainty for level energy, should have been blank, integer, or one of LT, GT, LE, GE, AP, CA, or SY, but we got \"" << line.substr(19, 2) << "\"" << std::endl;
                  }
                  // read string for spin and parity, half-life, and half-life uncertainty seperately
                  // reading a string from stringstream means we discard all leading whitespace and stop the moment we encounter more whitespace
                  // i.e. we read only one word
                  auto spinParity = line.substr(21, 18);
                  trimWS(spinParity);
                  auto halfLife = line.substr(39, 10);   // include the unit!
                  trimWS(halfLife);
                  auto halfLifeUncertainty = line.substr(49, 6);
                  trimWS(halfLifeUncertainty);
                  if(fVerbosity >= EVerbosity::kLoops) {
                     std::cout << "Adding new level " << energy << " +- " << energyUncertainty << ", " << spinParity << ", half-life " << halfLife << " +- " << halfLifeUncertainty << std::endl;
                     std::cout << "\"" << line.substr(9, 10) << "\", \"" << line.substr(19, 2) << "\", \"" << line.substr(21, 18) << "\", \"" << line.substr(39, 10) << "\", \"" << line.substr(49, 6) << "\"" << std::endl;
                  }
                  // create new level and set it's energy uncertainty
                  currentLevel = currentNucleus->AddLevel(energy, energyUncertainty, identifier);
               } else if(fVerbosity >= EVerbosity::kLoops) {
                  std::cout << "Ignoring level \"" << line << "\"" << std::endl;
               }
               break;
            case 'G':   // gamma record (1-5 nuclide, 6 blank or not '1' for cont., 7-9 " G ", 10-19 energy in keV, 20-21 uncert. 22-29 relative photon intens., 30-31 uncert., 32-41 multipol.
                        // 42-49 mixing ratio, 50-55 uncert., 56-62 conversion coeff., 63-64 uncert., 65-74 relative total intens., 75-76 uncert., 77 comment flag, 78 'C' confirmed coincidence, '?' questionable coincidence,
                        // 79 blank, 80 '?' questionable placement, 'S' expected but unobserved)
               if(line[5] == ' ' && line[6] == ' ') {
                  // read energy and uncertainty
                  double energy     = 0.;
                  char   identifier = '\0';   // not really used for gamma records?
                  str.str(line.substr(9, 10));
                  str.clear();
                  // sometime the next field is an X (with or without preceding whitespace), so check if the last read failed and if so, break
                  if(!ReadEnergy(str, energy, identifier)) {
                     if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Failed to read energy (got " << energy << " from \"" << str.str() << "\"/\"" << line.substr(9, 10) << "\"), skipping this line" << std::endl; }
                     rejectedGamma = true;
                     break;
                  }
                  if(energy == 0. && identifier != '\0') {
                     if(fVerbosity >= EVerbosity::kBasicFlow) {
                        std::cout << DYELLOW << "Found gamma ray with unknown energy (identifier = '" << identifier << "' = " << static_cast<int>(identifier) << "), ignoring it since we do not have a mechanism for this (yet?)" << RESET_COLOR << std::endl;
                     }
                     rejectedGamma = true;
                     break;
                  }
                  double       energyUncertainty = 0.;
                  EUncertainty uncertaintyLabel  = EUncertainty::kDefault;
                  str.str(line.substr(19, 2));
                  str.clear();
                  if(!ReadUncertainty(str, energyUncertainty, uncertaintyLabel)) {
                     std::cout << "Failed to read uncertainty for level energy, should have been blank, integer, or one of LT, GT, LE, GE, AP, CA, or SY, but we got \"" << line.substr(19, 2) << "\"" << std::endl;
                  }
                  // read relative photon intensity and uncertainty
                  double photonIntensity            = 0.;
                  double photonIntensityUncertainty = 0.;
                  str.str(line.substr(21, 8));
                  str.clear();
                  str >> photonIntensity;
                  str.str(line.substr(29, 2));
                  str.clear();
                  str >> photonIntensityUncertainty;
                  // read multipolarity
                  std::string multipolarity = line.substr(31, 10);
                  trimWS(multipolarity);
                  // read mixing ratio and uncertainty
                  double mixingRatio            = 0.;
                  double mixingRatioUncertainty = 0.;
                  str.str(line.substr(41, 8));
                  str.clear();
                  str >> mixingRatio;
                  str.str(line.substr(49, 6));
                  str.clear();
                  str >> mixingRatioUncertainty;
                  // read conversion coeff. and uncertainty
                  double conversionCoeff            = 0.;
                  double conversionCoeffUncertainty = 0.;
                  str.str(line.substr(55, 7));
                  str.clear();
                  str >> conversionCoeff;
                  str.str(line.substr(62, 2));
                  str.clear();
                  str >> conversionCoeffUncertainty;
                  // read relative total intensity and uncertainty
                  double totalIntensity            = 0.;
                  double totalIntensityUncertainty = 0.;
                  str.str(line.substr(64, 10));
                  str.clear();
                  str >> totalIntensity;
                  str.str(line.substr(74, 2));
                  str.clear();
                  str >> totalIntensityUncertainty;
                  if(currentGamma != nullptr) {
                     // we still have the previous gamma (which might have gotten additional information from a continuation record), so first add that to the level (or nucleus)
                     if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Adding gamma at " << currentGamma->Energy() << " keV to " << std::flush; }
                     if(currentLevel == nullptr) {
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "nucleus \"" << currentNucleus->GetName() << "\"" << std::endl;
                           currentNucleus->Print();
                        }
                        currentNucleus->AddTransition(currentGamma);
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "done:" << std::endl;
                           currentNucleus->Print();
                        }
                     } else {
                        // try to (quietly) add the transition to the current level, if it fails, add it to the nucleus
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "level at " << currentLevel->Energy() << " keV to " << std::endl;
                           currentLevel->Print();
                        }
                        if(currentLevel->AddTransition(currentGamma, true) == nullptr) {
                           currentNucleus->AddTransition(currentGamma);
                        }
                        if(fVerbosity >= EVerbosity::kLoops) {
                           std::cout << "done:" << std::endl;
                           currentLevel->Print();
                        }
                     }
                  }
                  currentGamma  = new TTransition(energy, energyUncertainty, photonIntensity, photonIntensityUncertainty, mixingRatio, mixingRatioUncertainty, conversionCoeff, conversionCoeffUncertainty, totalIntensity, totalIntensityUncertainty);
                  rejectedGamma = false;
                  if(fVerbosity >= EVerbosity::kLoops) {
                     std::cout << "Created gamma with energy " << energy << " +- " << energyUncertainty << ", " << photonIntensity << " +- " << photonIntensityUncertainty << ", " << multipolarity << ", " << mixingRatio << " +- " << mixingRatioUncertainty << ", " << conversionCoeff << " +- " << conversionCoeffUncertainty << ", " << totalIntensity << " +- " << totalIntensityUncertainty << std::endl;
                     std::cout << "\"" << line.substr(9, 10) << "\", \"" << line.substr(19, 2) << "\", \"" << line.substr(21, 8) << "\", \"" << line.substr(29, 2) << "\", \"" << line.substr(31, 10) << "\", \"" << line.substr(41, 8) << "\", \"" << line.substr(49, 6) << "\", \"" << line.substr(55, 7) << "\", \"" << line.substr(62, 2) << "\", \"" << line.substr(64, 10) << "\", \"" << line.substr(74, 2) << "\"" << std::endl;
                  }
               } else if(line[5] != '1' && line[6] == ' ') {
                  // continuation record
                  // currently we only check if a final level is given
                  if(line.substr(9, 3) == "FL=" && line[12] != '?') {
                     if(currentGamma == nullptr) {
                        // if we rejected the previous gamma ray (mainly because the energy is only "X"), it makes sense that we do not have a gamma ray now, so no warning needed
                        if(!rejectedGamma) {
                           std::cout << "Found final level information (\"" << line.substr(12, 68) << "\"), but have no current gamma ray to add this to? Nucleus is " << (currentNucleus == nullptr ? "also unknown" : currentNucleus->GetName()) << std::endl;
                        }
                        break;
                     }
                     // we have a final level with a known energy (might have a ? or a +x after it though)
                     // could also be FL=X or FL=X+energy
                     double finalLevel = 0.;
                     char   identifier = '\0';
                     str.str(line.substr(12, 68));
                     str.clear();
                     if(!ReadEnergy(str, finalLevel, identifier)) {
                        std::cout << DRED << "Error, failed to read energy of final level, got final level energy " << finalLevel << ", identifier '" << identifier << "' (" << static_cast<int>(identifier) << ") but the stringstream failed: \"" << str.str() << "\", going to ignore line " << lineNumber << ": \"" << line << "\"" << RESET_COLOR << std::endl;
                        break;
                     }
                     char next = '\0';
                     str.get(next);
                     if(next == '?') {
                        currentGamma->UncertainPlacement(true);
                     } else if(next == '+') {
                        str.get(identifier);
                     }
                     currentGamma->FinalLevel(finalLevel, identifier);
                     if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Added final level " << finalLevel << ", identifier '" << identifier << "' (" << static_cast<int>(identifier) << ") to gamma of " << currentGamma->Energy() << " keV" << std::endl; }
                  }
               } else if(fVerbosity >= EVerbosity::kLoops) {
                  std::cout << "Ignoring gamma \"" << line << "\"" << std::endl;
               }
               break;
            case 'B':   // beta record (1-5 nuclide, 6 blank or not '1' for cont., 7-9 " B ", ...), seems to be for decay data sets only?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring beta \"" << line << "\"" << std::endl; }
               break;
            case 'E':   // EC record (1-5 nuclide, 6 blank or not '1' for cont., 7-9 " E ", ...), seems to be for decay data sets only?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring EC \"" << line << "\"" << std::endl; }
               break;
            case 'A':   // alpha record (1-5 nuclide, 6 blank or not '1' for cont., 7-9 " A ", ...), seems to be for decay data sets only?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring alpha \"" << line << "\"" << std::endl; }
               break;
            case 'D':   // delayed particle record (1-5 nuclide, 6 blank or not '1' for cont., 7-8 " D", 9 particle N, P, or A, ...), seems to be for decay data sets only?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring delayed \"" << line << "\"" << std::endl; }
               break;
            case 'R':   // reference record (1-3 mass, 4-7 blank, 8 'R', 9 blank, ...) not present?
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring reference \"" << line << "\"" << std::endl; }
               break;
            case ' ':   // comment (if line[6] is 'c' or maybe 'C', 'd', 'D', 't', or 'T')
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Ignoring comment \"" << line << "\"" << std::endl; }
               break;
            default:
               // probably a comment record (1-5 nuclide, 6 blank or character not '1' for cont., 7 'c', 'D', 'T', or 't', plus other stuff we ignore)
               if(fVerbosity >= EVerbosity::kLoops) { std::cout << "Skipping unknown character " << line[7] << " from line \"" << line << "\"" << std::endl; }
               break;
            };
            ++lineNumber;
         } while(std::getline(file, line) && !line.empty() && file.good());
      } else {
         // not adopted levels
      }
   }

   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Done reading \"" << fileName << "\"" << std::endl;
   }

   file.close();

   return true;
}

bool TDataBase::LoadAbundance(const std::string& path)
{
   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Loading abundances from path \"" << path << "\"" << std::endl;
   }

   return true;
}

bool TDataBase::LoadMasses(const std::string& path)
{
   if(fVerbosity >= EVerbosity::kSubroutines) {
      std::cout << "Loading masses from path \"" << path << "\"" << std::endl;
   }

   return true;
}

void TDataBase::ListNuclei(bool print) const
{
   /// List all nuclei in the data base, if print flag is true, print each of them.

   std::cout << "Found " << fNuclei.size() << " nuclei:" << std::endl;
   for(const auto& [label, nucleus] : fNuclei) {
      std::cout << label << " = " << nucleus << std::endl;
      if(print) {
         nucleus->Print();
      }
   }
}
