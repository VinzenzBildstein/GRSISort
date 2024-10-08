//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 13:18:27 2016 by ROOT version 5.34/24
// from TTree FragmentTree/FragmentTree
// found on file: fragment07844_000.root
//////////////////////////////////////////////////////////

#ifndef CrossTalk_h
#define CrossTalk_h

#include "TChain.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

// Header file for the classes stored in the TTree if any.
#include "TGriffin.h"
#include "TSceptar.h"
#include "TGRSISelector.h"
#include "TGriffinBgo.h"
// Fixed size dimensions of array or collections stored in the TTree if any.

class CrossTalk : public TGRSISelector {

public:
   TGriffin*    fGrif;
   TSceptar*    fScep;
   TGriffinBgo* fGriffinBgo;

   explicit CrossTalk(TTree* /*tree*/ = 0) : TGRSISelector(), fGrif(0), fScep(0), fGriffinBgo(0) { SetOutputPrefix("Crosstalk"); }
   virtual ~CrossTalk() = default;
   virtual Int_t Version() const { return 2; }
   void          CreateHistograms();
   void          FillHistograms();
   void          InitializeBranches(TTree* tree);

   ClassDef(CrossTalk, 2);
};

#endif

#ifdef CrossTalk_cxx
void CrossTalk::InitializeBranches(TTree* tree)
{
   if(!tree) return;
   if(tree->SetBranchAddress("TGriffin", &fGrif) == TTree::kMissingBranch) {
      fGrif = new TGriffin;
   }
   if(tree->SetBranchAddress("TSceptar", &fScep) == TTree::kMissingBranch) {
      fScep = new TSceptar;
   }
   if(tree->SetBranchAddress("TGriffinBgo", &fGriffinBgo) == TTree::kMissingBranch) {
      fGriffinBgo = new TGriffinBgo;
   }
}

#endif   // #ifdef CrossTalk_cxx
