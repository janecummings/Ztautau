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

int get_flavor(SpinTree *spin, int jtau, bool &up, bool &down, bool &qg, bool &gg);
double getCosTheta_direct(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK);
double getCosTheta(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK );
double getPhi(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK );
bool CSFrame(TLorentzVector Z, TVector3 &CSAxis, TVector3 &xAxis, TVector3 &yAxis, bool qDirOK );
void fillAis( SpinTree *evtTree, int idhist, float evtWeight, int flavor );

float GeV = 1000.;
float kPI = TMath::Pi();

bool doFlavor = false;
bool DEBUG = false; 
//const int nflavors = 10;
//TString flavors[nflavors] = {"up_qq","down_qq","none_qq","up_qg","down_qg","none_qg","up_gg","down_gg","none_gg","all"}; 
const int nflavors = 1;
TString flavors[nflavors] = {"all"}; 

const int nkin = 15;
const int nhand = 2;
const int nmodes = 4;

TString kinematics[nkin] = {"Zpt","Zmass","TScostheta","cosa","xfrac","upsilon","taupt","vistaupt","taueta","vistaueta","CosDirect","CScostheta","x1","x2","scale"};
TString modes[nmodes] = {"0","1","3","4"};

TH1F* hKinematics[nkin][nmodes][nflavors][nhand];
TH2F* h2KinLH[nkin][nmodes][nkin];
TH2F* h2KinRH[nkin][nmodes][nkin];

const int ncoef = 5; const int nmoments = 9;
TH1F* hRW[nkin][nmodes][ncoef][nhand];
TH1F* hMoments[nmoments];
TH1F* hMomentsRW[nmoments][ncoef];
TString coef[ncoef] = {"A0","A1","A2","A3","A4"};
TString moments[nmoments] = {"A0","A1","A2","A3","A4","A5","A6","A7","A"};


int kbins[nkin] = {200,200, 120,  120, 120,  120, 100,  100, 120, 120, 120, 120,120,120,500};
double lbins[nkin] = {0,    0, -1.2, -1.2,  0, -1.2,   0,    0,  -6,  -6, -1.2,  0,  0,  0,  0};
double hbins[nkin] = {200,200,  1.2,  1.2, 1.2, 1.2, 200,  200,   6,   6,  1.2, 1.2,1.2,1.2,500};

const int nsamples = 5;
// Reweighting
const int nbinsRW = 9;
double A0[nsamples][nbinsRW] = {{0.1478, 0.1756, 0.0701, 0.105, 0.0722, 0.055, -0.0323, 0.0177, 0.1952},{0.1063, 0.1246, 0.0474, 0.0704, 0.0489, 0.0339, -0.018, 0.0088, 0.1136},{0.0949, 0.1085, 0.0419, 0.0585, 0.0412, 0.0279, -0.0162, 0.0063, 0.0729},{0.0861, 0.0836, 0.0357, 0.0537, 0.0372, 0.026, -0.0185, 0.0059, 0.0674},{0.0803, 0.0839, 0.0352, 0.0499, 0.0356, 0.0223, -0.0151, 0.0063, 0.0543}};
double A1[nsamples][nbinsRW] = {{-0.6413, -0.4885, -0.6794, -0.7322, -0.5131, -0.351, -0.5354, -0.7229, 0.6776},{-0.3152, -0.2828, -0.3924, -0.406, -0.3343, -0.2859, -0.5758, -0.411, 0.2073},{-0.2449, -0.238, -0.3181, -0.3385, -0.2457, -0.183, -0.1596, 0.5618, 0.0927},{-0.3342, -0.2343, -0.3332, -0.3349, -0.2856, -0.1881, -0.1058, -0.5077, -1.0527},{-0.2502, -0.2303, -0.3378, -0.2542, -0.2346, -0.2529, -0.1637, -0.2605, 0.0748}};
double A2[nsamples][nbinsRW] = {{0.0304, 0.0539, 0.0289, 0.0365, 0.051, 0.0846, 0.0936, -0.1148, -0.0141},{0.0328, 0.05, 0.0246, 0.0309, 0.0432, 0.0627, 0.0601, -0.0926, -0.009},{0.0311, 0.0573, 0.0291, 0.0316, 0.0412, 0.0724, 0.0834, -0.0894, -0.0099},{0.0592, 0.0609, 0.0265, 0.0338, 0.0488, 0.0918, 0.0441, -0.0686, -0.0072},{0.0336, 0.046, 0.0327, 0.0334, 0.0932, 0.0704, 0.1023, -0.128, -0.0111}};
double A3[nsamples][nbinsRW] = {{-0.2052, -0.2655, 1.2358, -0.2358, -0.1727, -0.1312, -0.2243, -0.3852, 0.0418},{-0.1281, -0.0909, 2.9354, -0.1349, -0.0871, -0.0642, -0.1011, -0.1611, 0.0153},{-0.1077, -0.098, 0.6622, -0.108, -0.0755, -0.0501, -0.0838, -0.1267, 0.0124},{-0.0917, -0.0925, 0.2587, -0.0895, -0.0718, -0.0498, -0.092, -0.1359, 0.0089},{-0.0908, -0.0713, -0.2843, -0.0818, -0.0639, -0.0366, -0.0461, -0.0932, 0.0071}};
double A4[nsamples][nbinsRW] = {{1.4144, 1.2486, 1.208, 1.2508, 1.211, 1.243, 1.1972, 1.1907, 1.1252},{1.6876, 1.5391, 3.1213, 1.4787, 1.4297, 1.4832, 1.3921, 1.4185, 1.2652},{1.9485, 1.5135, 2.8237, 1.6879, 1.6663, 1.8056, 1.5687, 1.7599, 1.4972},{2.3313, 2.4483, 3.1606, 1.8358, 1.7884, 1.902, 1.7539, 1.7407, 1.6381},{2.5276, 1.6804, 0.8445, 2.0642, 1.7822, 1.9081, 1.9577, 1.9127, 1.7623}};


int main(int argc, char *argv[])
{
  if (argc<2) {
    //cout << " usage : " << endl;
    return 0;
  }

  //TString outFileName = argv[1];
  TString sample = argv[1];
  TString outFileName = "CS_Reweight_test.root";
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
  if (sample == "Np1") isam = 0;
  else if (sample == "Np2") isam = 1;
  else if (sample == "Np3") isam = 2;
  else if (sample == "Np4") isam = 3;
  else if (sample == "Np5") isam = 4;

  //Set Branches
  //SetBranches(chain);
  
  TFile* out = new TFile(outFileName,"UPDATE");

  TString hist_name;
  for(int ihk=0; ihk < nkin; ihk++){
    for( int im = 0; im < nmodes; im++){
      for (int ifl=0; ifl < nflavors; ifl++){
        hist_name = sample + "_" + kinematics[ihk] + "_" + flavors[ifl] + "_" + modes[im];
        hKinematics[ihk][im][ifl][0] = new TH1F( hist_name + "_LH", hist_name+"_LH", kbins[ihk], lbins[ihk], hbins[ihk]);
        hKinematics[ihk][im][ifl][1] = new TH1F( hist_name + "_RH", hist_name+"_RH", kbins[ihk], lbins[ihk], hbins[ihk]);
        hKinematics[ihk][im][ifl][0]->Sumw2();
        hKinematics[ihk][im][ifl][1]->Sumw2();
        //cout<<hKinematics[ihk][im][ifl][0]->GetBinCenter(hKinematics[ihk][im][ifl][0]->GetNbinsX())<<endl;
      }
      for(int ikk=0; ikk<nkin; ikk++){
        hist_name = sample + "_" +  kinematics[ihk] + "_" + kinematics[ikk] + "_" + modes[im];
        h2KinRH[ihk][im][ikk] = new TH2F( hist_name + "_RH", hist_name+"_RH", kbins[ihk], lbins[ihk],hbins[ihk],kbins[ikk],lbins[ikk],hbins[ikk] );
        h2KinLH[ihk][im][ikk] = new TH2F( hist_name + "_LH", hist_name+"_LH", kbins[ihk], lbins[ihk],hbins[ihk],kbins[ikk],lbins[ikk],hbins[ikk] );
        h2KinRH[ihk][im][ikk]->Sumw2();
        h2KinLH[ihk][im][ikk]->Sumw2();
      }
      for(int ick=0;ick<ncoef;ick++){
        hist_name = sample + "_" + kinematics[ihk] + "_" + coef[ick] + "_" + modes[im];
        hRW[ihk][im][ick][0] = new TH1F( hist_name + "_LH", hist_name+"_LH", kbins[ihk], lbins[ihk], hbins[ihk]);
        hRW[ihk][im][ick][1] = new TH1F( hist_name + "_RH", hist_name+"_RH", kbins[ihk], lbins[ihk], hbins[ihk]);
        hRW[ihk][im][ick][0]->Sumw2(); hRW[ihk][im][ick][1]->Sumw2();
      }
    }
  }

  TString h_name;
  TString hrw_name;
  for(int ihm=0; ihm < nmoments; ihm++){
    h_name = sample + "_" + moments[ihm];
    hMoments[ihm] = new TH1F( h_name, h_name, 140, 60, 200);
    hMoments[ihm]->Sumw2();
    for (int ick=0; ick < ncoef; ick++){
      hrw_name = h_name + "_" + coef[ick];
      hMomentsRW[ihm][ick] = new TH1F( hrw_name, hrw_name, 140, 60, 200);
      hMomentsRW[ihm][ick]->Sumw2();

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
  int n_up = 0;
  int n_down = 0;
  int n_total = 0;
  int n_qq = 0;
  int n_qg = 0;
  int n_gg = 0;


  int itau = 0;

  //loop over events
  for(int jentry = 0; jentry < numEvents; jentry++){
    int ientry = spin->LoadTree(jentry);
    if(ientry<0) break;
    nb = spin->GetEntry(jentry); nbytes += nb;

    // Get Event Flavor
    bool up = false; bool down = false; bool qg = false; bool gg = false;
    string flavor; int f;
    if (doFlavor){
      for(int jtau = 0; jtau < spin->mc_n; jtau++){
        if (spin->mc_pdgId->at(jtau) == 23 && spin->mc_status->at(jtau) == 3){
          n_total += get_flavor( spin, jtau, up, down, qg, gg );

          //cout<<n_total<<endl;
          //cout<<up<<" "<<down<<" "<<qg<<" "<<gg<<endl;
          if (up) { n_up++; flavor = "up_"; f = 0;}
          else if (down) { n_down++; flavor = "down_"; f = 1;}
          else {flavor = "none_"; f = 2;}
          if (gg) {n_gg++; flavor+="gg"; f+= 6;}
          else if (qg) {n_qg++; flavor +="qg"; f+=3;}
          else {n_qq++; flavor+="qq";}
          break;
        } //end inside hard Z
      }
    }
    //if (jentry > 1000) break;
    //cout<<jentry<<endl;
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
    /*
    cout<<m<<endl;
    cout<<"tau pt "<<tau_had.Pt()<<endl;
    cout<<"tau eta "<<tau_had.Eta()<<endl;
    cout<<"vis tau pt "<<tau_vis.Pt()<<endl;
    cout<<"vis tau eta "<<tau_vis.Eta()<<endl;
    cout<<"***********************"<<endl;
    */
    //TString kinematics[nkin] = {"Zpt","Zmass","TScostheta","cosa","xfrac","upsilon","taupt","vistaupt","taueta","vistaueta"};
    vector<double> fillvalues(nkin);
    fillvalues[6] = tau_had.Pt();
    fillvalues[8] = tau_had.Eta();
    fillvalues[7] = tau_vis.Pt();
    fillvalues[9] = tau_vis.Eta();
    fillvalues[5] = (pi.Pt() - pi0.Pt())/(pi.Pt() + pi0.Pt());
    fillvalues[4] = ( tau_vis.Pt() / tau_had.Pt() );

    TLorentzVector Zbos = TauSpinnerHelpers::simp2tlv(dp.X);
    fillvalues[0] = Zbos.Pt();
    fillvalues[1] = Zbos.M();

    TauSpinner::Particle tau_minus, tau_plus; 
    TLorentzVector minus_tlv, plus_tlv;
    if (dp.tau2.pdgid() == 15 ) {
      tau_minus = TauSpinnerHelpers::simp2pcle(dp.tau2); minus_tlv = TauSpinnerHelpers::simp2tlv(dp.tau2);
      tau_plus = TauSpinnerHelpers::simp2pcle(dp.tau); plus_tlv = TauSpinnerHelpers::simp2tlv(dp.tau);
    }
    else {
      tau_minus = TauSpinnerHelpers::simp2pcle(dp.tau); minus_tlv = TauSpinnerHelpers::simp2tlv(dp.tau);
      tau_plus = TauSpinnerHelpers::simp2pcle(dp.tau2); plus_tlv = TauSpinnerHelpers::simp2tlv(dp.tau2);
    }

    TLorentzVector tlv_tt = minus_tlv + plus_tlv;
    double mtt = tlv_tt.M();
    double evtwt[ncoef] = {1.0,1.0,1.0,1.0,1.0} ;
    int wti = -1;
    if (isam >= 0){ 
      if (mtt >= 60 && mtt < 80) wti = 0;
      else if (mtt >= 80 && mtt < 85) wti = 1;
      else if (mtt >= 85 && mtt < 90) wti = 2;
      else if (mtt >= 90 && mtt < 95) wti = 3;
      else if (mtt >= 95 && mtt < 100) wti = 4;
      else if (mtt >= 100 && mtt < 110) wti = 5;
      else if (mtt >= 110 && mtt < 120) wti = 6;
      else if (mtt >= 120 && mtt < 140) wti = 7;
      else if (mtt >= 140) wti = 8;

      if (wti >= 0) {
        evtwt[0] = (20./3)*A0[isam][wti]; (evtwt[1]=A1[isam][wti];evtwt[2]=A2[isam][wti];evtwt[3]=A3[isam][wti];evtwt[4]=A4[isam][wti];
      }
    }

    double cosdirect = getCosTheta_direct(Zbos, minus_tlv , plus_tlv , 1);
    fillvalues[10] = cosdirect;
    double cosCS = getCosTheta( Zbos, minus_tlv, plus_tlv, 1);
    fillvalues[11] = cosCS;

    //hard scatter x1,x2
    fillvalues[12] = spin->mcevt_pdf_x1->at(0);
    fillvalues[13] = spin->mcevt_pdf_x2->at(0);
    fillvalues[14] = spin->mcevt_pdf_scale->at(0);

    double S = dp.X.e()*dp.X.e() - dp.X.px()*dp.X.px() - dp.X.py()*dp.X.py() - dp.X.pz()*dp.X.pz();
    //cout<<TauSpinner::calculate_costhe_and_WID(S, tau_minus, tau_plus )[0]<<endl;

    //tau rest frame
    
    TauSpinner::Particle P_QQ = TauSpinner::Particle ( lf_tau.px() + lf_nu_tau.px(), lf_tau.py() + lf_nu_tau.py(), lf_tau.pz() + lf_nu_tau.pz(), lf_tau.e() + lf_nu_tau.e(), 0 );
    TauSpinner::Particle lf_tau_nu = TauSpinnerHelpers::simp2pcle(tau_daughters[0]); // the neutrino
    // move tau, lep, and tau neutrino to Z rest frame
    lf_tau.boostToRestFrame(P_QQ);
    lf_nu_tau.boostToRestFrame(P_QQ);
    lf_tau_nu.boostToRestFrame(P_QQ);

    // in z rest frame, cos theta calculation a la TauSpinner
    TauSpinner::Particle P_B1(0, 0, 1, 1, 0);
    TauSpinner::Particle P_B2(0, 0,-1, 1, 0);
    P_B1.boostToRestFrame(P_QQ);
    P_B2.boostToRestFrame(P_QQ);
    tau_minus.boostToRestFrame(P_QQ); tau_plus.boostToRestFrame(P_QQ);

    double costheta1 = (tau_plus.px()*P_B1.px()    +tau_plus.py()*P_B1.py()    +tau_plus.pz()*P_B1.pz()    ) /
                   sqrt(tau_plus.px()*tau_plus.px()+tau_plus.py()*tau_plus.py()+tau_plus.pz()*tau_plus.pz()) /
                   sqrt(P_B1.px()    *P_B1.px()    +P_B1.py()    *P_B1.py()    +P_B1.pz()    *P_B1.pz()    );

    double costheta2 = (tau_minus.px()*P_B2.px()    +tau_minus.py()*P_B2.py()    +tau_minus.pz()*P_B2.pz()    ) /
                   sqrt(tau_minus.px()*tau_minus.px()+tau_minus.py()*tau_minus.py()+tau_minus.pz()*tau_minus.pz()) /
                   sqrt(P_B2.px()    *P_B2.px()    +P_B2.py()    *P_B2.py()    +P_B2.pz()    *P_B2.pz()    );
    double sintheta1 = sqrt(1-costheta1*costheta1);
    double sintheta2 = sqrt(1-costheta2*costheta2);

    // Cosine of hard scattering
    double costhe = (costheta1*sintheta2 + costheta2*sintheta1) / (sintheta1 + sintheta2);
    fillvalues[2] = costhe;

    // move tau and tau neutrino so tau is on z axis; boost tau neutrino into the tau rest frame
    double phi = lf_tau.getAnglePhi();
    lf_tau.rotateXY(-1*phi);
    double theta = lf_tau.getAngleTheta();
    lf_tau.rotateXZ(M_PI - theta);
    lf_tau_nu.rotateXY(-1*phi); lf_tau_nu.rotateXZ(M_PI - theta); lf_tau_nu.boostAlongZ(-1*lf_tau.pz(), lf_tau.e());
    // get angle between tau neutrino and tau direction ( z axis )
    TLorentzVector tl_tau = TLorentzVector( lf_tau.px(), lf_tau.py(), lf_tau.pz(), lf_tau.e()  );
    TLorentzVector tl_tau_nu = TLorentzVector( lf_tau_nu.px(), lf_tau_nu.py(), lf_tau_nu.pz(), lf_tau_nu.e()  );
    double cosa = cos(M_PI - tl_tau_nu.Angle(tl_tau.Vect()));
    fillvalues[3] = cosa;

    int ih = hand;
    if (ih == -1) ih = 0;
    itau++;

    //TH1F* hKinematics[nkin][nmodes][nflavors][nhand];
    //TH2F* h2KinLH[nkin][nmodes][nkin];

    fillAis( spin, 0, 1., -1);
    for (int ick=0; ick<ncoef; ick++){
      cout<<ick<<" ";
      fillAis( spin, 0, evtwt[ick], ick);
    }

    continue;

    if (!doFlavor) f=0;
    for (int ikk=0; ikk < nkin; ikk++){
      if (im) {
        hKinematics[ikk][im][f][ih]->Fill( fillvalues[ikk] );
        hKinematics[ikk][1][f][ih]->Fill( fillvalues[ikk] ); // single prong
        if (doFlavor) {
          hKinematics[ikk][im][9][ih]->Fill( fillvalues[ikk] );
          hKinematics[ikk][1][9][ih]->Fill( fillvalues[ikk] );
        }
      }
      hKinematics[ikk][0][f][ih]->Fill( fillvalues[ikk] ); // all had tau decays 
      if (doFlavor) hKinematics[ikk][0][9][ih]->Fill( fillvalues[ikk] );

      for(int ick=0;ick<ncoef;ick++){
        if (im) {
          hRW[ikk][im][ick][ih]->Fill( fillvalues[ikk], evtwt[ick]);
          hRW[ikk][1][ick][ih]->Fill( fillvalues[ikk], evtwt[ick]);
        }
        hRW[ikk][0][ick][ih]->Fill( fillvalues[ikk], evtwt[ick]);
      }
      
      for (int ikj=0; ikj < nkin; ikj++){
        if (ih == 0) {
          if(im) { 
            h2KinLH[ikk][im][ikj]->Fill( fillvalues[ikk],fillvalues[ikj]  );
            if (im != 1) h2KinLH[ikk][1][ikj]->Fill( fillvalues[ikk],fillvalues[ikj]  );
          }
          h2KinLH[ikk][0][ikj]->Fill( fillvalues[ikk],fillvalues[ikj]  );
        }
        else if (ih == 1) {
          if (im) {
            h2KinRH[ikk][im][ikj]->Fill( fillvalues[ikk], fillvalues[ikj]);
            if (im != 1) h2KinRH[ikk][1][ikj]->Fill( fillvalues[ikk], fillvalues[ikj]);
          }
          h2KinRH[ikk][0][ikj]->Fill( fillvalues[ikk], fillvalues[ikj]);
        
       }
     }
      //cout<<kinematics[ikk]<<" "<<fillvalues[ikk]<<endl;
    }


    //if (jentry>5) break;

  } // end loop over jentry
  if (doFlavor){
    cout<<"itau: "<<itau<<endl;
    cout<<"n_total: "<<n_total<<endl;
    cout<<"n_up: "<<n_up<<endl;
    cout<<"n_down: "<<n_down<<endl; 
    cout<<"n_gg: "<<n_gg<<endl;
    cout<<"n_qg: "<<n_qg<<endl;
    cout<<"n_qq: "<<n_qq<<endl;
  }

  // for(int ihk=0; ihk < nkin; ihk++){
  //   for( int im = 0; im < nmodes; im++){
  //     for (int ifl=0; ifl < nflavors; ifl++){
  //       hKinematics[ihk][im][ifl][0]->Write();
  //       hKinematics[ihk][im][ifl][1]->Write();
  //     }
  //     for(int ikk=0; ikk<nkin; ikk++){
  //       h2KinRH[ihk][im][ikk]->Write();
  //       h2KinLH[ihk][im][ikk]->Write();
  //     }
  //     for(int ick=0; ick<ncoef;ick++){


  //     }
  //   }
  // }  

  out->Write();
  out->Close();

  delete chain; chain=0;
  return 0;
   
}// end main function

int get_flavor(SpinTree *spin, int jtau, bool &up, bool &down, bool &qg, bool &gg ) { 
  int ng = 0;
  int pdg = 0;
  if (spin->mc_parent_index->at(jtau).size() != 2) {
    cout<<"Size of Z parent vector unexpected: "<<spin->mc_parent_index->at(jtau).size()<<endl;
    return 0;
  }
  for (int iq=0; iq < spin->mc_parent_index->at(jtau).size(); iq++){
    int qi = spin->mc_parent_index->at(jtau).at(iq);
    //cout<<qi<<" "<<spin->mc_pdgId->at(qi)<<endl;
    pdg = fabs(spin->mc_pdgId->at(qi));
    if (pdg == 1 || pdg == 3 || pdg == 5) {
      down = true; 
    }
    if ( pdg == 2 || pdg == 4 || pdg == 6) {
      up = true; 
    }
    if (spin->mc_pdgId->at(qi) == 21) ng++;
  }

  if (up && down) {
    int ii = spin->mc_parent_index->at(jtau).at(0);
    int ij = spin->mc_parent_index->at(jtau).at(1);
    //cout<<spin->mc_pdgId->at(ii)<<" "<<spin->mc_pdgId->at(ij)<<endl;
  }
  //else if (up) n_up++;
  //else if (down) n_down++;
  if (ng == 2) gg = true;//n_gg++;
  else if (ng == 1) qg = true;//n_qg++;
  //else if (ng == 0) n_qq++; 

  if ( ng == 2){
    int g_up = 0; int g_down = 0;
    int pi = spin->mc_parent_index->at(jtau).at(0);
    //cout<<spin->mc_pdgId->at(pi)<<endl;
    for (int ip = 0; ip < spin->mc_child_index->at(pi).size();ip++){
      int ipp = spin->mc_child_index->at(pi).at(ip);
      //cout<<spin->mc_pdgId->at(ipp)<<" ";
      int ipg =  fabs(spin->mc_pdgId->at(ipp));
      if ( ipg == 1 || ipg == 3 || ipg == 5) g_down++;
      else if (ipg == 2 || ipg == 4 || ipg == 6) g_up++;
    }
  
    if ( g_up > g_down) up = true;
    if ( g_down > g_up) down = true;
  }
  if (up && down) { up = false; down = false;}
  //cout<<up<<" "<<down<<" "<<pdg<<endl;
  return 1;
  //if (up and !down) n_up++;
  //if (down and !up) n_down++;
  //cout<<endl;
  //cout<<ng<<" "<<up<<" "<<down<<" "<<n_up<<" "<<n_down<<endl;
}

double getCosTheta_direct(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK = true){
//--------------------------------------------------                                                                                                
    double cosTheta;
    double sign  = fabs(Z.Z())/Z.Z();
    if (!qDirOK)
       sign *= -1;
    double part1 = 2 / ( Z.M()*sqrt(Z.M()*Z.M() + Z.Pt()*Z.Pt() ));
    double part2 = (muM.Plus()*muP.Minus() - muP.Plus()*muM.Minus())/2;
    cosTheta = part1*part2;
    // add sign according to boost direktion in z
    cosTheta *= sign;

    return cosTheta;
}
double getCosTheta(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK = true){
//--------------------------------------------------                                                                                                
    TVector3 CSAxis, yAxis, xAxis;
    TLorentzVector boostedMu = muM;
    boostedMu.Boost(-Z.BoostVector());

    if (!CSFrame(Z, CSAxis, xAxis, yAxis, qDirOK))
       return -999;
    return cos(boostedMu.Angle(CSAxis));
}
bool CSFrame(TLorentzVector Z, TVector3 &CSAxis, TVector3 &xAxis, TVector3 &yAxis, bool qDirOK = true){
//--------------------------------------------------                                                                                                
     double ProtonMass = 938.272; // MeV
     double BeamEnergy = 4000000; // MeV

     double sign  = fabs(Z.Z())/Z.Z();
     bool isGood = true;
     TLorentzVector p1, p2;

     if (qDirOK){
         p1.SetPxPyPzE(0, 0, sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)); // quark
         p2.SetPxPyPzE(0, 0, -1*sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)); // antiquark
     } else {
         p1.SetPxPyPzE(0, 0, -1*sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)); // quark
         p2.SetPxPyPzE(0, 0, sign*BeamEnergy, TMath::Sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)); // antiquark
     }
     p1.Boost(-Z.BoostVector());
     p2.Boost(-Z.BoostVector());
     CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit();
     yAxis = (p1.Vect().Unit()).Cross((p2.Vect().Unit()));
     yAxis = yAxis.Unit();
     xAxis = yAxis.Cross(CSAxis);
     xAxis = xAxis.Unit();

     return isGood;
}

void fillAis( SpinTree *evtTree, int idhist, float evtWeight, int ick ) {

//----------------------------------------------------                                                                                                               

  // Fill TLorentz vectors
  TLorentzVector tlv_epos, tlv_eneg, tlv_ee;
  TLorentzVector tlv_beam_pos, tlv_beam_neg;

  float ProtonMass = 0.938; // GeV
  float BeamEnergy =  4000; // GeV
  tlv_beam_pos.SetPxPyPzE(0, 0,  BeamEnergy, BeamEnergy); 
  tlv_beam_neg.SetPxPyPzE(0, 0,- BeamEnergy, BeamEnergy); 

  int refStatus = 3; //changed from 1 to 3
  for(int ipart = 0; ipart < evtTree->mc_n; ipart++){
    if(DEBUG)
      std::cout << "mc_status=" << evtTree->mc_status->at(ipart) << "     mc_pdgId=" << evtTree->mc_pdgId->at(ipart)
                << "    pt=" << evtTree->mc_pt->at(ipart)        << "     eta="      <<  evtTree->mc_eta->at(ipart) << std::endl;
    if( evtTree->mc_status->at(ipart) == 23) evtTree->mc_status->at(ipart) == 3; // changed status from 1 to 3
    //pdgd == -11 means charge=1
    if( evtTree->mc_status->at(ipart) == refStatus &&  evtTree->mc_pdgId->at(ipart) == -15){ // changed 11 to 15
      tlv_epos.SetPtEtaPhiM(evtTree->mc_pt->at(ipart)/1000., evtTree->mc_eta->at(ipart), evtTree->mc_phi->at(ipart), 1.777);
    } else if ( evtTree->mc_status->at(ipart) == refStatus &&  evtTree->mc_pdgId->at(ipart) == 15){
      tlv_eneg.SetPtEtaPhiM(evtTree->mc_pt->at(ipart)/1000., evtTree->mc_eta->at(ipart), evtTree->mc_phi->at(ipart), 1.777);
    }
  }

  // combain into pair
  tlv_ee = tlv_epos + tlv_eneg;

  //components of the cos(th*) calculation
  float Lplus, Lminus, Pplus, Pminus;
  Lplus  = tlv_eneg.E()+tlv_eneg.Pz();
  Lminus = tlv_eneg.E()-tlv_eneg.Pz();
  Pplus  = tlv_epos.E()+tlv_epos.Pz();
  Pminus = tlv_epos.E()-tlv_epos.Pz();
            
  //does the cos(th*) calculation
  float cosThetaCS;
  cosThetaCS  = (Lplus*Pminus - Lminus*Pplus);
  cosThetaCS *= fabs(tlv_ee.Pz());
  cosThetaCS /= (tlv_ee.Mag()*tlv_ee.Pz());
  cosThetaCS /= sqrt(tlv_ee.Mag2() + tlv_ee.Pt()*tlv_ee.Pt() );

  if( DEBUG) {  
    std::cout << "tlv_eneg.E()/GeV=: " << tlv_eneg.E()/GeV << std::endl;
    std::cout << "tlv_epos.E()/GeV=: " << tlv_epos.E()/GeV << std::endl;
    std::cout << "tlv_eneg.Pt()/GeV=: " << tlv_eneg.Pt()/GeV << std::endl;
    std::cout << "tlv_epos.Pt()/GeV=: " << tlv_epos.Pt()/GeV << std::endl;
    std::cout << "tlv_eneg.Pseudorapidity()/GeV=: " << tlv_eneg.PseudoRapidity()/GeV << std::endl;
    std::cout << "tlv_epos.PseudoRapidity()/GeV=: " << tlv_epos.PseudoRapidity()/GeV << std::endl;
    std::cout << "tlv_ee.E()/GeV=: " << tlv_ee.E()/GeV << std::endl;
    std::cout << "tlv_ee.Pt()/GeV=: " << tlv_ee.Pt()/GeV << std::endl;
    std::cout << "tlv_ee.Pseudorapidity()/GeV=: " << tlv_ee.PseudoRapidity()/GeV << std::endl;
    std::cout << "tlv_ee.Mag()/GeV=: " << tlv_ee.Mag()/GeV << std::endl;
    std::cout << "cosThetaCS=: " << cosThetaCS << std::endl;
  }

  // Turn off 80 - 100 GeV Mass Window  
  //int accept = 1;
  //if( tlv_ee.Mag()  < 80.0 || tlv_ee.Mag()  > 100.0 ) accept=0;
  //if( accept == 0) return;


  float  phiCS =  getPhi(tlv_ee, tlv_eneg, tlv_epos, 1);

  float cosPhiCS   = cos(phiCS);
  float sinPhiCS   = sin(phiCS);
  float cos2PhiCS  = cos(2*phiCS);
  float sin2PhiCS  = sin(2*phiCS);
  float sinThetaCS = 0.0;
  if( (1-cosThetaCS*cosThetaCS) > 0 ) sinThetaCS = sqrt(1-cosThetaCS*cosThetaCS);
  float sin2ThetaCS = 2 * sinThetaCS*cosThetaCS;

  float wt_A0 = 0.5 * (1 - 3*cosThetaCS*cosThetaCS);
  float wt_A1 = 2*cosThetaCS * sinThetaCS * cosPhiCS;
  float wt_A2 = sinThetaCS * sinThetaCS * cos2PhiCS;
  float wt_A3 = sinThetaCS * cosPhiCS;
  float wt_A4 = cosThetaCS;
  float wt_A5 = sinThetaCS*sinThetaCS*sin2PhiCS;
  float wt_A6 = sin2ThetaCS*sinPhiCS;
  float wt_A7 = sinThetaCS*sinPhiCS;

  cout<<tlv_ee.M()<<" "<<evtWeight<<endl;

  if (ick>=0){
    hMomentsRW[8][ick]->Fill( tlv_ee.M(), evtWeight );
    hMomentsRW[0][ick]->Fill( tlv_ee.M(), wt_A0*evtWeight );  
    hMomentsRW[1][ick]->Fill( tlv_ee.M(), wt_A1*evtWeight );
    hMomentsRW[2][ick]->Fill( tlv_ee.M(), wt_A2*evtWeight );
    hMomentsRW[3][ick]->Fill( tlv_ee.M(), wt_A3*evtWeight );
    hMomentsRW[4][ick]->Fill( tlv_ee.M(), wt_A4*evtWeight );
    hMomentsRW[5][ick]->Fill( tlv_ee.M(), wt_A5*evtWeight );
    hMomentsRW[6][ick]->Fill( tlv_ee.M(), wt_A6*evtWeight );
    hMomentsRW[7][ick]->Fill( tlv_ee.M(), wt_A7*evtWeight );
  }
  else{
    hMoments[0]->Fill( tlv_ee.M(), wt_A0 );
    hMoments[1]->Fill( tlv_ee.M(), wt_A1 );
    hMoments[2]->Fill( tlv_ee.M(), wt_A2 );
    hMoments[3]->Fill( tlv_ee.M(), wt_A3 );
    hMoments[4]->Fill( tlv_ee.M(), wt_A4 );
    hMoments[5]->Fill( tlv_ee.M(), wt_A5 );
    hMoments[6]->Fill( tlv_ee.M(), wt_A6 );
    hMoments[7]->Fill( tlv_ee.M(), wt_A7 );
    hMoments[8]->Fill( tlv_ee.M());
  }
    //kinematics
    //hKinematics[0][f]->Fill( tlv_ee.Pt());
    //hKinematics[1][f]->Fill( tlv_ee.M());
    //hKinematics[2][f]->Fill( getCosTheta_direct(tlv_ee, tlv_eneg, tlv_epos, 1)); 
}

double getPhi(TLorentzVector Z, TLorentzVector muM, TLorentzVector muP, bool qDirOK = true){
//--------------------------------------------------                                                                                                
    TLorentzVector boostedMu = muM;
    TVector3 CSAxis, xAxis, yAxis;

     boostedMu.Boost(-Z.BoostVector());

     CSFrame(Z, CSAxis, xAxis, yAxis, qDirOK);
     double phi = atan2((boostedMu.Vect()*yAxis),(boostedMu.Vect()*xAxis));
     if(phi<0) phi = phi + 2*kPI;

    return phi;
}
