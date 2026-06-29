#ifndef ExampleTreeHelper_h
#define ExampleTreeHelper_h

/** \addtogroup Helpers
 * @{
 */

////////////////////////////////////////////////////////////////////////////////
/// \class ExampleTreeHelper
///
/// This selector shows how to create a tree with selected events (beta-tagged
/// with gamma multiplicities of at least three), and selected information
/// (suppressed addback energies, and beta-gamma timing differences).
///
////////////////////////////////////////////////////////////////////////////////

// Header file for the classes stored in the TTree if any.
#include "TGRSIHelper.h"
#include "TGriffin.h"
#include "TGriffinBgo.h"
#include "TSceptar.h"
#include "TZeroDegree.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class ExampleTreeHelper : public TGRSIHelper, public ROOT::Detail::RDF::RActionImpl<ExampleTreeHelper> {
public:
   explicit ExampleTreeHelper(TList* list)
      : TGRSIHelper(list)
   {
      Prefix("ExampleTree");
      Setup();
   }
   // These functions are expected to exist
   ROOT::RDF::RResultPtr<std::map<std::string, TList>> Book(ROOT::RDataFrame* d) override
   {
      auto columnNames = ColumnNames();
      std::cout << "Using column names";
      for(const auto& name : columnNames) {
         std::cout << " \"" << name << "\"";
      }
      std::cout << " we get:" << std::endl;
      if(std::none_of(columnNames.begin(), columnNames.end(), [](const std::string& name) { return name == "TGriffin"; }) ||
         std::none_of(columnNames.begin(), columnNames.end(), [](const std::string& name) { return name == "TGriffinBgo"; })) {
         throw std::runtime_error("Missing either TGriffin or TGriffinBgo in list of column names?");
      }
      fZds     = std::any_of(columnNames.begin(), columnNames.end(), [](const std::string& name) { return name == "TZeroDegree"; });
      fSceptar = std::any_of(columnNames.begin(), columnNames.end(), [](const std::string& name) { return name == "TSceptar"; });
      std::cout << "Zds " << (fZds ? "present" : "missing") << ", and Sceptar " << (fSceptar ? "present" : "missing") << std::endl;
      if(fZds && fSceptar) {
         return d->Book<TGriffin, TGriffinBgo, TZeroDegree, TSceptar>(std::move(*this), {"TGriffin", "TGriffinBgo", "TZeroDegree", "TSceptar"});
      }
      if(fZds) {
         return d->Book<TGriffin, TGriffinBgo, TZeroDegree>(std::move(*this), {"TGriffin", "TGriffinBgo", "TZeroDegree"});
      }
      if(fSceptar) {
         return d->Book<TGriffin, TGriffinBgo, TSceptar>(std::move(*this), {"TGriffin", "TGriffinBgo", "TSceptar"});
      }
      return d->Book<TGriffin, TGriffinBgo>(std::move(*this), {"TGriffin", "TGriffinBgo"});
   }
   void CreateHistograms(unsigned int slot) override;
   void Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TZeroDegree& zds, TSceptar& scep);
   void Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TZeroDegree& zds)
   {
      TSceptar scep;
      Exec(slot, grif, grifBgo, zds, scep);
   }
   void Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo, TSceptar& scep)
   {
      TZeroDegree zds;
      Exec(slot, grif, grifBgo, zds, scep);
   }
   void Exec(unsigned int slot, TGriffin& grif, TGriffinBgo& grifBgo)
   {
      TZeroDegree zds;
      TSceptar    scep;
      Exec(slot, grif, grifBgo, zds, scep);
   }

private:
   // branches of output trees
   // all of these need to be maps for the different slots/workers, we're using maps as they don't need to be resized and once accessed
   // the address of key stays the same (important to be able to create branches)
   std::map<unsigned int, double*>             fSuppressedAddback2;    ///< vector of suppressed addback energies
   std::map<unsigned int, std::vector<double>> fSuppressedAddback;     ///< vector of suppressed addback energies
   std::map<unsigned int, std::vector<double>> fBetaGammaTiming;       ///< vector of beta-gamma timing
   std::map<unsigned int, int>                 fGriffinMultiplicity;   ///< multiplicity of suppressed addback energies
   bool                                        fZds{false};            ///< flag whether ZDS is present in the data
   bool                                        fSceptar{false};        ///< flag whether SCEPTAR is present in the data
};

extern "C" ExampleTreeHelper* CreateHelper(TList* list) { return new ExampleTreeHelper(list); }

extern "C" void DestroyHelper(TGRSIHelper* helper) { delete helper; }

/*! @} */
#endif
