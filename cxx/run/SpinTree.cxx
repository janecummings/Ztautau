#include "SpinTree.h"
#include "TObject.h"
#include "math.h"
#include "TRandom.h"


ClassImp(SpinTree) //integrate class into Root
SpinTree::SpinTree()
{
  //access to leafs
}
SpinTree::~SpinTree()
{
   if(!fChain) return;
   delete fChain->GetCurrentFile();
}
void SpinTree::Init(TTree *tree)
{
 
  EventNumber = 0;
  EF_mu24i_tight = 0;

  mc_n = 0;
  mc_pdgId = 0;
  mc_status = 0;
  mc_child_index = 0;
  mc_parents = 0;//
  mc_parent_index = 0;
  mc_pt = 0;
  mc_eta = 0;
  mc_phi = 0;
  mc_m = 0;

  mcevt_n = 0;
  mcevt_signal_process_id = 0;
  mcevt_event_number = 0;
  mcevt_event_scale = 0; 
  mcevt_alphaQCD = 0;
  mcevt_alphaQED = 0;
  mcevt_pdf_id1 = 0;
  mcevt_pdf_id2 = 0;
  mcevt_pdf_x1 = 0;
  mcevt_pdf_x2 = 0;
  mcevt_pdf_scale = 0;
  mcevt_pdf1 = 0;
  mcevt_pdf2 = 0;
  mcevt_weight = 0;
  mcevt_nparticle = 0;
  mcevt_pileUpType = 0;


  polari = 0;
  tau_mode = 0;
  tau2_mode = 0;

  

  if (!tree) return;
  fChain = tree;
  fCurrent = -1;

  fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);

  fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
  fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
  fChain->SetBranchAddress("mc_child_index", &mc_child_index, &b_mc_child_index);
  fChain->SetBranchAddress("mc_parent_index", &mc_parent_index, &b_mc_parent_index);
  fChain->SetBranchAddress("mc_parents", &mc_parents, &b_mc_parents);
  fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
  fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
  fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
  fChain->SetBranchAddress("mc_m", &mc_m, &b_mc_m);

  fChain->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
  fChain->SetBranchAddress("mcevt_signal_process_id", &mcevt_signal_process_id, &b_mcevt_signal_process_id);
  fChain->SetBranchAddress("mcevt_event_number", &mcevt_event_number, &b_mcevt_event_number);
  fChain->SetBranchAddress("mcevt_event_scale", &mcevt_event_scale, &b_mcevt_event_scale);
  fChain->SetBranchAddress("mcevt_alphaQCD", &mcevt_alphaQCD, &b_mcevt_alphaQCD);
  fChain->SetBranchAddress("mcevt_alphaQED", &mcevt_alphaQED, &b_mcevt_alphaQED);
  fChain->SetBranchAddress("mcevt_pdf_id1", &mcevt_pdf_id1, &b_mcevt_pdf_id1);
  fChain->SetBranchAddress("mcevt_pdf_id2", &mcevt_pdf_id2, &b_mcevt_pdf_id2);
  fChain->SetBranchAddress("mcevt_pdf_x1", &mcevt_pdf_x1, &b_mcevt_pdf_x1);
  fChain->SetBranchAddress("mcevt_pdf_x2", &mcevt_pdf_x2, &b_mcevt_pdf_x2);
  fChain->SetBranchAddress("mcevt_pdf_scale", &mcevt_pdf_scale, &b_mcevt_pdf_scale);
  fChain->SetBranchAddress("mcevt_pdf1", &mcevt_pdf1, &b_mcevt_pdf1);
  fChain->SetBranchAddress("mcevt_pdf2", &mcevt_pdf2, &b_mcevt_pdf2);
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight, &b_mcevt_weight);
  fChain->SetBranchAddress("mcevt_nparticle", &mcevt_nparticle, &b_mcevt_nparticle);
  fChain->SetBranchAddress("mcevt_pileUpType", &mcevt_pileUpType, &b_mcevt_pileUpType);

  //
  fChain->SetBranchAddress("polari",&polari,&b_polari);
  fChain->SetBranchAddress("tau_mode",&tau_mode,&b_tau_mode);
  fChain->SetBranchAddress("tau2_mode",&tau2_mode,&b_tau2_mode);

  //selection
  fChain->SetBranchAddress("EF_mu24i_tight",&EF_mu24i_tight,&b_EF_mu24i_tight);

}

int SpinTree::LoadTree(int entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   int centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) 
     {
       fCurrent = chain->GetTreeNumber();
       Notify();
     }
   return centry;
}
int SpinTree::GetEntry(int entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
int SpinTree::GetNEn()
{
   if (fChain == 0) return 0;
   fChain->GetEntries(); 
   return int(fChain->GetEntriesFast());
}
