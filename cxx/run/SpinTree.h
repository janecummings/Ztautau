#ifndef SpinTree_h
#define SpinTree_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>

using namespace std;

class SpinTree: public TObject
{
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  UInt_t          EventNumber;

  Int_t           mc_n;
  vector<int>     *mc_pdgId;
  vector<int>     *mc_status;
  vector<vector<int> > *mc_child_index;
  vector<vector<int> > *mc_parent_index;
  vector<vector<int> > *mc_parents;
  vector<float>   *mc_pt;
  vector<float>   *mc_eta;
  vector<float>   *mc_phi;
  vector<float>   *mc_m;

  Int_t           mcevt_n;
  vector<int>     *mcevt_signal_process_id;
  vector<int>     *mcevt_event_number;
  vector<double>  *mcevt_event_scale;
  vector<double>  *mcevt_alphaQCD;
  vector<double>  *mcevt_alphaQED;
  vector<int>     *mcevt_pdf_id1;
  vector<int>     *mcevt_pdf_id2;
  vector<double>  *mcevt_pdf_x1;
  vector<double>  *mcevt_pdf_x2;
  vector<double>  *mcevt_pdf_scale;
  vector<double>  *mcevt_pdf1;
  vector<double>  *mcevt_pdf2;
  vector<vector<double> > *mcevt_weight;
  vector<int>     *mcevt_nparticle;
  vector<short>   *mcevt_pileUpType;

  double          polari;
  Int_t           tau_mode;
  Int_t           tau2_mode;

  //
  Bool_t          EF_mu24i_tight;
  
  
  TBranch        *b_EventNumber;   //!
  TBranch        *b_EF_mu24i_tight;

  TBranch        *b_mc_n;   //!
  TBranch        *b_mc_pt;   //!
  TBranch        *b_mc_m;   //!
  TBranch        *b_mc_eta;   //!
  TBranch        *b_mc_phi;   //!
  TBranch        *b_mc_status;   //!
  TBranch        *b_mc_pdgId;   //!
  TBranch        *b_mc_child_index;   //!
  TBranch        *b_mc_parent_index;   //!
  TBranch        *b_mc_parents;   //!

  TBranch        *b_mcevt_n;   //!
  TBranch        *b_mcevt_signal_process_id;   //!
  TBranch        *b_mcevt_event_number;   //!
  TBranch        *b_mcevt_event_scale;   //!
  TBranch        *b_mcevt_alphaQCD;   //!
  TBranch        *b_mcevt_alphaQED;   //!
  TBranch        *b_mcevt_pdf_id1;   //!
  TBranch        *b_mcevt_pdf_id2;   //!
  TBranch        *b_mcevt_pdf_x1;   //!
  TBranch        *b_mcevt_pdf_x2;   //!
  TBranch        *b_mcevt_pdf_scale;   //!
  TBranch        *b_mcevt_pdf1;   //!
  TBranch        *b_mcevt_pdf2;   //!
  TBranch        *b_mcevt_weight;   //!
  TBranch        *b_mcevt_nparticle;   //!
  TBranch        *b_mcevt_pileUpType;   //!



  TBranch        *b_polari;
  TBranch        *b_tau_mode;
  TBranch        *b_tau2_mode;

  //Constructor
  SpinTree();
  //Destructor
  ~SpinTree();
  
  
  void Init(TTree *tree);
  int GetEntry(int entry);
  int GetNEn();
  int LoadTree(int entry);
  
  
 private:
  
  ClassDef(SpinTree,1)
    };
#endif
