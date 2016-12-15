#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TFile.h"
#include "math.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "SpinTree.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinnerTool.h"

using namespace std;
using namespace TauSpinner;


float GeV = 1000.;
float kPI = TMath::Pi();

bool DEBUG = false; 
bool addwts = true;
//TH1F* weights;


const int nkin = 8;
const int nhand = 2;
const int nmodes = 4;

TString kinematics[nkin] = {"Zpt","Zmass","upsilon","taupt","ptlep","vistaupt","taueta","vistaueta"};
TString modes[nmodes] = {"0","1","3","4"};

TH1F* hKinematics[nkin][nmodes][nhand];
TH1F* hWeighted[nkin][nmodes][nhand];

int kbins[nkin] = {200,200, 120,  100,  100, 100, 120, 120 };
double lbins[nkin] = {0, 0, -1.2, 0,     0,     0, -6,   -6};
double hbins[nkin] = {200,200, 1.2, 100, 100, 100, 6,   6};

const int nsamples = 6;

int main(int argc, char *argv[])
{
  if (argc<2) {
    //cout << " usage : " << endl;
    return 0;
  }


  //if (addwts){
  TFile* rwf = new TFile("ZpTWeights.root");
  TH1F* weights = (TH1F*)rwf->Get("ZPt_weights");


  //}

  //TString outFileName = argv[1];
  TString sample = argv[1];
  TString outFileName = "ZptReweight2016.root";
  vector<std::string> fileList;

  bool DEBUG = false;//true;
  
  std::string argStr;
  ifstream files(argv[2]);
  std::getline(files,argStr);
  for (size_t i=0,n; i <= argStr.length(); i=n+1) {
    n = argStr.find_first_of(',',i);
    if (n == std::string::npos) n = argStr.length();
    std::string tmp = argStr.substr(i,n-i);
    fileList.push_back(tmp);
  }


  TChain *chain = new TChain("tau");
  Int_t addedFiles = 0;
  for (unsigned int iFile=0; iFile<fileList.size(); iFile++) {
    cout<<"open "<< fileList[iFile].c_str() << endl;
    addedFiles += chain->Add(fileList[iFile].c_str());
  }
  cout<<"added "<<addedFiles<<" files to the TChain"<<endl;

  int isam = -1;
  if (sample == "Np0") isam = 0;
  else if (sample == "Np1") isam = 1;
  else if (sample == "Np2") isam = 2;
  else if (sample == "Np3") isam = 3;
  else if (sample == "Np4") isam = 4;
  else if (sample == "Np5") isam = 5;

  //Set Branches
  //SetBranches(chain);
  


  TFile* out = new TFile("../output/" + outFileName,"UPDATE");

  TString hist_name;
  for(int ihk=0; ihk < nkin; ihk++){
    for( int im = 0; im < nmodes; im++){
      hist_name = sample + "_" + kinematics[ihk] + "_" + modes[im];
      hKinematics[ihk][im][0] = new TH1F( hist_name + "_LH", hist_name+"_LH", kbins[ihk], lbins[ihk], hbins[ihk]);
      hKinematics[ihk][im][1] = new TH1F( hist_name + "_RH", hist_name+"_RH", kbins[ihk], lbins[ihk], hbins[ihk]);
      hKinematics[ihk][im][0]->Sumw2();
      hKinematics[ihk][im][1]->Sumw2();
      hWeighted[ihk][im][0] = new TH1F( hist_name + "_LH_weighted", hist_name+"_LH_weighted", kbins[ihk], lbins[ihk], hbins[ihk]);
      hWeighted[ihk][im][1] = new TH1F( hist_name + "_RH_weighted", hist_name+"_RH_weighted", kbins[ihk], lbins[ihk], hbins[ihk]);
      hWeighted[ihk][im][0]->Sumw2();
      hWeighted[ihk][im][1]->Sumw2();
    }
  }

  gROOT->ProcessLine("#include <vector>");
  SpinTree* spin = new SpinTree(); 
  spin->Init(chain);
  
  
  TauSpinnerHelpers::D3PD_MC* d3pd_mc = new TauSpinnerHelpers::D3PD_MC();
  //following HSG4 implementation but need Ipol = 1, not 2 (pol = <0> for higgs)
  Tauolapp::Tauola::initialize();
  double CMSENE = 8000.0; // center of mass system energy
  bool Ipp = true; //for pp collsions
  int Ipol = 1; // events generated with all spin effects         
  int nonSM2 = 0; //turn nonSM calculations off
  int nonSMN = 0; //turn on if calculation nonSM weight for shapes only
  TauSpinner::initialize_spinner(Ipp,Ipol,nonSM2,nonSMN,CMSENE);
  //Initialize PDFs
  std::string pdf="MRSTMCal.LHgrid";
  LHAPDF::initPDFSetByName(pdf);
  TauSpinnerHelpers::TauSpinnerTool tool;
  
  int numEvents = spin->GetNEn();
  int skippedEvents = 0;
  int nbytes = 0, nb = 0;

  cout<<"Sample: "<<sample<<" Number of Events: "<<numEvents<<endl;

  
  out->cd();

  double hand,weight;
  int mode, mode_2;

  int itau = 0;

  //loop over events
  for(int jentry = 0; jentry < numEvents; jentry++){
    int ientry = spin->LoadTree(jentry);
    if(ientry<0) break;
    nb = spin->GetEntry(jentry); nbytes += nb;

    //pass branches to d3pd_mc interface
    d3pd_mc->set_addresses(
         spin->mc_pt, 
         spin->mc_eta, 
         spin->mc_phi, 
         spin->mc_m, 
         spin->mc_pdgId, 
         spin->mc_status, 
         spin->mc_child_index
         ); 
    //read in event with TauSpinnerTool
    tool.read_event(*d3pd_mc);
    //get weight -- polari gets filled so can used getTauSpin from tau_reweight_lib
    weight = tool.get_spin_weight();
    hand = tool.get_helicity();
    TauSpinnerHelpers::DecayParticles dp = tool.get_decay_particles();
    mode = TauSpinnerHelpers::get_decay_mode(dp.tau,dp.tau_daughters);
    mode_2 = TauSpinnerHelpers::get_decay_mode(dp.tau2,dp.tau2_daughters);
    int m;
    //cout<<dp.tau.pdgid()<<" "<<dp.tau2.pdgid()<<endl;
    TLorentzVector tau_had;
    TLorentzVector nu_tau;
    vector<TauSpinner::SimpleParticle> tau_daughters;
    vector<TauSpinner::SimpleParticle> nu_tau_daughters;
    vector<TauSpinner::Particle> lf_tau_daughters;
    TauSpinner::Particle lf_tau,lf_nu_tau;

    if (mode == 1 || mode == 2) {
      m = mode_2;
      tau_had = TauSpinnerHelpers::simp2tlv(dp.tau2);
      nu_tau = TauSpinnerHelpers::simp2tlv(dp.tau);
      tau_daughters = dp.tau2_daughters;
      nu_tau_daughters = dp.tau_daughters;
      lf_tau = TauSpinnerHelpers::simp2pcle(dp.tau2);
      lf_nu_tau = TauSpinnerHelpers::simp2pcle(dp.tau);
    }
    else {
      m = mode;
      tau_had = TauSpinnerHelpers::simp2tlv(dp.tau);
      nu_tau = TauSpinnerHelpers::simp2tlv(dp.tau2);
      tau_daughters = dp.tau_daughters;
      nu_tau_daughters = dp.tau2_daughters;
      lf_tau = TauSpinnerHelpers::simp2pcle(dp.tau);
      lf_nu_tau = TauSpinnerHelpers::simp2pcle(dp.tau2);
    }
    //cout<<"**************************"<<endl;
    //cout<<m<<endl;
    // mode filter
    int ncharged = 0;
    if (m==5){
      for( int p=0; p < tau_daughters.size(); p++){
        if (fabs(tau_daughters[p].pdgid()) == 211) ncharged++;
        else if (fabs(tau_daughters[p].pdgid()) == 321) ncharged++;
      }
    }
    int im = 0;
    if (ncharged == 1) im = 1; // single prong
    if (m==3 || m==6) im = 2; // single charged pion
    if (m==4 || m == 7) im = 3; // 1p1n


    // get visible tau
    TLorentzVector tau_vis;
    vector<TauSpinner::SimpleParticle>::iterator id;
    TLorentzVector pi0,pi;
    //loop over daughters and leave out neutrinos
    for (id=tau_daughters.begin(); id!=tau_daughters.end();id++){
      if ( fabs(id->pdgid()) == 16 || fabs(id->pdgid()) == 22 ) continue;
      tau_vis+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
      if( fabs(id->pdgid()) == 211 || fabs(id->pdgid())==321){
        pi += TauSpinnerHelpers::simp2tlv(*id );
      }
      if( fabs(id->pdgid()) == 111 || fabs(id->pdgid()) == 130 || fabs(id->pdgid()) == 310){
        pi0 += TauSpinnerHelpers::simp2tlv(*id);
      }
    }
    //TString kinematics[nkin] = {"Zpt","Zmass","upsilon","taupt","ptlep","vistaupt","taueta","vistaueta"};

    vector<double> fillvalues(nkin);

    TLorentzVector Zbos = TauSpinnerHelpers::simp2tlv(dp.X);
    fillvalues[0] = Zbos.Pt();
    fillvalues[1] = Zbos.M();
    fillvalues[2] = (pi.Pt() - pi0.Pt())/(pi.Pt() + pi0.Pt()); //upsilon

    fillvalues[3] = tau_had.Pt();
    fillvalues[4] = nu_tau.Pt();
    fillvalues[6] = tau_had.Eta();
    fillvalues[5] = tau_vis.Pt();
    fillvalues[7] = tau_vis.Eta();

    rwf->cd();
    // // Get event weight
    float evtwt = 1.0;
    if (addwts){
      int zbin = weights->FindBin(Zbos.Pt());
      if (zbin > 31) {
        evtwt = weights->GetBinContent(31);
      }
      else {
        evtwt = weights->GetBinContent(zbin);
      }

    }
    //evtwt = 1.;
    //if (jentry>100) break;
  


    //fillvalues[4] = ( tau_vis.Pt() / tau_had.Pt() );
    int ih = hand;
    if (ih == -1) ih = 0;
    itau++;

    out->cd();

    for (int ikk=0; ikk < nkin; ikk++){
      if (im) {
        if (im==2 || im==3) {
          hKinematics[ikk][im][ih]->Fill( fillvalues[ikk] , 1.0);
          hWeighted[ikk][im][ih]->Fill( fillvalues[ikk] , evtwt);
        }
        hKinematics[ikk][1][ih]->Fill( fillvalues[ikk] , 1.0); // single prong
        hWeighted[ikk][1][ih]->Fill( fillvalues[ikk] , evtwt ); // single prong

        }
      
      hKinematics[ikk][0][ih]->Fill( fillvalues[ikk] , 1.0); // all had tau decays 
      hWeighted[ikk][0][ih]->Fill( fillvalues[ikk] , evtwt); // all had tau decays 
    }

    //if (jentry>50) break;

  } // end loop over jentry

  out->Write();
  out->Close();

  delete chain; chain=0;
  return 0;
   
}// end main function
