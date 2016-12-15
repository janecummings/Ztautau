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

#include "PythiaTree.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinnerTool.h"

using namespace std;
using namespace TauSpinner;


float GeV = 1000.;
float kPI = TMath::Pi();

bool DEBUG = false;

//const int nkin = 12;
const int nZkin = 2;
const int nkin = 14;
const int nreco = 7;
const int nres = 6;
const int nhand = 2;
const int nmodes = 6;

//TString kinematics[nkin] = {"Zpt","Zmass","TScostheta","cosa","xfrac","upsilon","taupt","vistaupt","taueta","vistaueta","CosDirect","CScostheta"};
TString Zkin[nZkin] = {"Zmass", "Zpt"};
int zbins[nZkin] = {140,120};
double zlbins[nZkin] = {60,0};
double zhbins[nZkin] = {200,120};

TString kinematics[nkin] = {"upsilon", "vistaupt","taupt","leppt","pi0pt","pipt", "xfrac", "costh","xcosth","upsilon20","upsilon60","upsilon40", "upsilonx0","upsilonx1"};
int kbins[nkin] = {120,  100, 100, 100, 80,80,120, 120, 120, 120,120,120,120,120};
double lbins[nkin] = { -1.2, 0,   0,   0,  0,0, 0, -1.2, -1.2 , -1.2,-1.2,-1.2,-1.2,-1.2};
double hbins[nkin] = { 1.2, 100, 100, 100, 80,80,1.2, 1.2, 1.2, 1.2,1.2,1.2,1.2,1.2};


TString reco[nreco] = {"upsilon","upsilon20","upsilon40","upsilon60","taupt", "trackpt", "pi0pt"};
int rkbins[nreco] = {150, 150, 150, 150, 100, 60, 80};
double rlbins[nreco] = {-1., -1,  -1,  -1,   0, 0, 0};
double rhbins[nreco] = {2 , 2, 2, 2, 100, 60, 80};

TString resolution[nres] = {"track","calo", "ups", "trackpt","calopt", "upspro"};
int sbins[nres] = {80,200, 200, 80,  100,120 };
double slbins[nres] = {0.6, 0, 0, 0, 0, -1.2};
double shbins[nres] = {1.4, 2, 2, 80,100, 1.2};


TString modes[nmodes] = {"0","1","2","3","4","5"};
// 0 = All Taus
// 1 = Single prong taus = 1 charged pion
// 2 = 3 prong taus = 3 charged pions
// 3 = pi nu decays
// 4 = pi pi0 nu decays
// 5 = other type 1 decays 

TH1F* hZkin[nZkin][nhand];
TH1F* hKinematics[nkin][nmodes][nhand];
TH1F* hReco[nreco][nmodes][nhand];
TH1F* hResolution[nres][nmodes][nhand];
TH2F* hUpsilon[nreco][nmodes][nhand];


int main(int argc, char *argv[])
{

 if (argc<2) {
    //cout << " usage : " << endl;
    return 0;
  }

  //TString outFileName = argv[1];
  TString sample = argv[1];
  TString outFileName = "PythiaTestXX.root";
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

  //Set Branches
  //SetBranches(chain);

  
  TFile* out = new TFile(outFileName,"UPDATE");

  //Z kinematics
  TString z_hist_name;
  for(int ihk=0; ihk < nZkin; ihk++){
    z_hist_name =  "Zkin_" + Zkin[ihk] + "_0"; 
    hZkin[ihk][0] = new TH1F( z_hist_name + "_LH", z_hist_name+"_LH", zbins[ihk], zlbins[ihk], zhbins[ihk]);
    hZkin[ihk][1] = new TH1F( z_hist_name + "_RH", z_hist_name+"_RH", zbins[ihk], zlbins[ihk], zhbins[ihk]);
    hZkin[ihk][0]->Sumw2();
    hZkin[ihk][1]->Sumw2();
      //cout<<hKinematics[ihk][im][0]->GetBinCenter(hKinematics[ihk][im][0]->GetNbinsX())<<endl;
  }


  // True Kinematic Plots
  TString hist_name;
  for(int ihk=0; ihk < nkin; ihk++){
    for( int im = 0; im < nmodes; im++){
      hist_name =  "True_" + kinematics[ihk] + "_" + modes[im];
      hKinematics[ihk][im][0] = new TH1F( hist_name + "_LH", hist_name+"_LH", kbins[ihk], lbins[ihk], hbins[ihk]);
      hKinematics[ihk][im][1] = new TH1F( hist_name + "_RH", hist_name+"_RH", kbins[ihk], lbins[ihk], hbins[ihk]);
      hKinematics[ihk][im][0]->Sumw2();
      hKinematics[ihk][im][1]->Sumw2();
      //cout<<hKinematics[ihk][im][0]->GetBinCenter(hKinematics[ihk][im][0]->GetNbinsX())<<endl;
    }
  }
  // Reco Kinematic Plots
  TString reco_hist_name;
  for(int ihk=0; ihk < nreco; ihk++){
    for( int im = 0; im < nmodes; im++){
      reco_hist_name = "Reco_" + reco[ihk] + "_" + modes[im];
      hReco[ihk][im][0] = new TH1F( reco_hist_name + "_LH", reco_hist_name+"_LH", rkbins[ihk], rlbins[ihk], rhbins[ihk]);
      hReco[ihk][im][1] = new TH1F( reco_hist_name + "_RH", reco_hist_name+"_RH", rkbins[ihk], rlbins[ihk], rhbins[ihk]);
      hReco[ihk][im][0]->Sumw2();
      hReco[ihk][im][1]->Sumw2();
      //cout<<hReco[ihk][im][0]->GetBinCenter(hReco[ihk][im][0]->GetNbinsX())<<endl;
    }
  }

  //Resolution Plots
  TString res_hist_name;
  for(int ihk=0; ihk < nres; ihk++){
    for( int im = 0; im < nmodes; im++){
      res_hist_name = "Res_" + resolution[ihk] + "_" + modes[im];
      hResolution[ihk][im][0] = new TH1F( res_hist_name + "_LH", res_hist_name+"_LH", sbins[ihk], slbins[ihk], shbins[ihk]);
      hResolution[ihk][im][1] = new TH1F( res_hist_name + "_RH", res_hist_name+"_RH", sbins[ihk], slbins[ihk], shbins[ihk]);
      hResolution[ihk][im][0]->Sumw2();
      hResolution[ihk][im][1]->Sumw2();
      //cout<<hResolution[ihk][im][0]->GetBinCenter(hResolution[ihk][im][0]->GetNbinsX())<<endl;
    }
  }

  //TH2
  TString th2_hist_name;
  for(int ihk=0; ihk < nreco; ihk++){
    for( int im = 0; im < nmodes; im++){
      th2_hist_name = "upsilon_" + reco[ihk] + "_" + modes[im];

      hUpsilon[ihk][im][0] = new TH2F( th2_hist_name + "_LH", th2_hist_name+"_LH", kbins[0], lbins[0], hbins[0], rkbins[ihk], rlbins[ihk], rhbins[ihk]);
      hUpsilon[ihk][im][1] = new TH2F( th2_hist_name + "_RH", th2_hist_name+"_RH", kbins[0], lbins[0], hbins[0], rkbins[ihk], rlbins[ihk], rhbins[ihk]);
      hUpsilon[ihk][im][0]->Sumw2();
      hUpsilon[ihk][im][1]->Sumw2();
      //cout<<hReco[ihk][im][0]->GetBinCenter(hReco[ihk][im][0]->GetNbinsX())<<endl;
    }
  }  

  // Start TauSpinner  
  gROOT->ProcessLine("#include <vector>");
  PythiaTree* spin = new PythiaTree(); 
  spin->Init(chain);  
  TauSpinnerHelpers::D3PD_MC* d3pd_mc = new TauSpinnerHelpers::D3PD_MC();
  //following HSG4 ihplementation but need Ipol = 1, not 2 (pol = <0> for higgs)
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

  int nBADMODE = 0;
  int nTAU = 0; // number of true taus
  int nMATCHED = 0; //number of reco taus matched to 
  int nXMATCHED = 0; // number of true taus matched to >1 reco tau
  int nMATCHEDSP = 0;
  int nSP = 0; // number type 1
  int nNP = 0; // number type 2
  int nPINU = 0; // number type 3
  int nRHO = 0; // number type 4
  int nPINU_X = 0;
  int nRHO_X = 0;
  int nSXP = 0; // number type 5
  int nXX = 0; //others 
  int nKX = 0; // number type XX that are Knu or Kpi0nu



  //loop over events
  for(int jentry = 0; jentry < numEvents; jentry++){
    int ientry = spin->LoadTree(jentry);
    if(ientry<0) break;
    nb = spin->GetEntry(jentry); nbytes += nb;

    if (jentry % 1000 == 0) cout<<jentry<<endl;

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

    if (m==-1){
      nBADMODE++;
      //continue; // skip weird bad modes
    }

    nTAU++; // count all true taus

    //cout<<"**************************"<<endl;
    //cout<<m<<endl;
    // mode filter
    // 1 - e nu nu, e nu nu gamma
    // 2 - mu nu nu, mu nu nu gamma 
    // 3 - pi nu
    // 5 - pi pi0 pi0 nu, pi pi pi nu
    // 6 - K nu
    // 4 - pi pi0 nu
    // 7 - K pi0 nu, pi K0 nu
    // 8 - pi pi pi pi0 nu
    // 9 - pi pi0 pi0 pi0 nu
    int ncharged = 0;
    int nPi0 = 0;
    int nPi = 0;
    int nK = 0;
    int nK0 = 0;
    int nDP = 0;

    for( int p=0; p < tau_daughters.size(); p++){
      int pdg = fabs(tau_daughters[p].pdgid());
      if (pdg == 16 || pdg == 22) continue;
      nDP++;
      if (pdg == 211) nPi++;
      if (pdg == 111) nPi0++;
      if (pdg == 130 || pdg == 310 || pdg == 311) nK0++;
      if (pdg == 321) nK++;
        //else if (fabs(tau_daughters[p].pdgid()) == 321) ncharged++;
    }
    
    //MODES
    int im = 0;
    if (nPi == 1) {
      if (nPi0 == 0 && nDP == 1) im = 3;
      if (nPi0 == 1 && nDP == 2) im = 4;
      if (nPi0 > 1) im = 5;
    } 
    if (nPi == 3) im = 2;


    // COUNT MODES
    if (im == 5 || im == 3 || im == 4) nSP++;
    if (im == 5) nSXP++;
    if (im == 3) nPINU++;
    if (im == 4) nRHO++;
    if (im == 2) nNP++;

    if ( im == 3 && nDP > 1) nPINU_X++;
    if (im == 4 && nDP > 2) nRHO_X++;
    if (nK0 || nK) nKX++;
    if (im == 0) nXX++;

    //HAND
    int ih = hand;
    if (ih == -1) ih = 0;


    // Z kinematics

    vector<double> fillZvalues(nZkin);
    TLorentzVector Zbos = TauSpinnerHelpers::simp2tlv(dp.X);
    fillZvalues[0] = Zbos.M();
    fillZvalues[1] = Zbos.Pt();

        // Fill True Z kinematics
    for (int ikk=0; ikk < nZkin; ikk++){
      hZkin[ikk][ih]->Fill( fillZvalues[ikk] ); // all had tau decays 
    }



    //Fill Truth Kinematics vector
    //{"upsilon","vistaupt","taupt","leppt","pi0pt","pipt", "xfrac"};
    vector<double> fillkinematics(nkin);
    vector<TauSpinner::SimpleParticle>::iterator il;
    //loop over daughters and leave out neutrinos
    float lep = 0;
    for (il=nu_tau_daughters.begin(); il!=nu_tau_daughters.end();il++){
      if ( fabs(il->pdgid()) ==  11 ||  fabs(il->pdgid()) ==  13 ) {
        lep = TLorentzVector(il->px(),il->py(),il->pz(),il->e()).Pt();
        break;
      }
    }


    // visible true tau
    TLorentzVector tau_vis;
    TLorentzVector tau_pi0;
    TLorentzVector tau_pi;
    vector<TauSpinner::SimpleParticle>::iterator id;
    bool gotvis = false;
    //loop over daughters and leave out neutrinos
    for (id=tau_daughters.begin(); id!=tau_daughters.end();id++){
      if ( fabs(id->pdgid()) == 16 || fabs(id->pdgid()) == 22 ) continue;
      gotvis = true;
      tau_vis+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
      if (fabs(id->pdgid()) == 111 ) tau_pi0+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
      else if (fabs(id->pdgid()) == 211 ) tau_pi+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
    }
    if(!gotvis) cout<<"what happened?"<<endl;

    float pipt = tau_pi.Pt();
    float pi0pt = tau_pi0.Pt();
    float upsilon = (pipt - pi0pt)/(pipt+pi0pt);
    float taupt = tau_vis.Pt();
    float xfrac = tau_vis.Pt()/tau_had.Pt();
    fillkinematics[0] = upsilon;
    fillkinematics[1] = taupt;
    fillkinematics[2] = tau_had.Pt();
    fillkinematics[3] = lep;
    fillkinematics[4] = pi0pt;
    fillkinematics[5] = pipt;
    fillkinematics[6] = xfrac;

    // cos theta observable
    TauSpinner::Particle lf_tau_nu = TauSpinnerHelpers::simp2pcle(tau_daughters[0]); 
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
    fillkinematics[7] = cosa;
    fillkinematics[8] = cosa;

    //upsilon studies
    fillkinematics[9] = upsilon;
    fillkinematics[10] = upsilon;
    fillkinematics[11] = upsilon;
    fillkinematics[12] = upsilon;
    fillkinematics[13] = upsilon;
    // Fill True taus
    float wtx = 1;
    for (int ikk=0; ikk < nkin; ikk++){
      if (ikk == 8) wtx = xfrac; 
      else wtx = 1;
      if (ikk == 9 && taupt > 20) continue;
      if (ikk == 10 && (taupt < 20 || taupt>40)) continue;
      if (ikk == 11 && taupt<40) continue;
      if (ikk == 12 && xfrac > 0.4) continue;
      if (ikk == 13 && xfrac < 0.4) continue; 

      if (im) {
        hKinematics[ikk][im][ih]->Fill( fillkinematics[ikk], wtx ); // fill modes
        if (im>2) hKinematics[ikk][1][ih]->Fill( fillkinematics[ikk],wtx ); // single prong
      }
      hKinematics[ikk][0][ih]->Fill( fillkinematics[ikk],wtx ); // all had tau decays 
    }

    // Reconstructed taus

    int nmatched = 0;

    TLorentzVector reco_tau;
    int mtau = -1; // matched tau index
    for (int itau=0; itau<spin->tau_n; itau++){
//      if (spin->tau_numTrack->at(itau)>1) continue; // only single prong reco taus
      reco_tau.SetPtEtaPhiM( spin->tau_pt->at(itau)/1000.0, spin->tau_eta->at(itau),spin->tau_phi->at(itau),spin->tau_m->at(itau));
      if ( reco_tau.DeltaR( tau_vis) < 0.2 ) {
        nmatched++;
        if (nmatched>1) continue;
        mtau = itau;
      }
    }
    if (nmatched>1) {
      nXMATCHED++;
      cout<<"Nmatched is "<<nmatched<<endl;
    }
    if (nmatched != 1) continue;
    nMATCHED++;

    // require sp reco taus 
    if (spin->tau_numTrack->at(mtau) != 1) continue; // only single prong reco taus
    nMATCHEDSP++;

    // fill reco {"upsilon","upsilon20","upsilon40","upsilon60","taupt", "trackpt", "pi0pt"};
    vector<double> fillreco(nreco);
    float rtaupt = spin->tau_pt->at(mtau)/1000.0;
    float rtrack = spin->tau_leadTrkPt->at(mtau)/1000.0;
    float rupsilon = 2*rtrack/rtaupt - 1;

    fillreco[0] = rupsilon;
    fillreco[1] = rupsilon;
    fillreco[2] = rupsilon;
    fillreco[3] = rupsilon;
    fillreco[4] = rtaupt;
    fillreco[5] = rtrack;
    fillreco[6] = rtaupt - rtrack;

    // measure resolution {"track","calo", "ups", "trackpt","calopt", "upspro"};

    vector<double> fillres(nres);
    vector<double> wtres(nres);
    fillres[0] = rtrack/pipt;
    wtres[0] = 1.;
    fillres[1] = rtaupt/taupt;
    wtres[1] = 1.;
    fillres[2] = rupsilon/upsilon;
    wtres[2] = 1.;
    fillres[3] = pipt; 
    wtres[3] = fabs(rtrack - pipt)/pipt; // this is a weight
    fillres[4] = taupt;
    wtres[4] = fabs(rtaupt - taupt)/taupt;
    fillres[5] = upsilon;
    wtres[5] = fabs(rupsilon-upsilon);



    //Fill Reco Taus
    //{"upsilon","upsilon20","upsilon40","upsilon60","taupt", "trackpt", "pi0pt"};
    for (int ikk=0; ikk < nreco; ikk++){
      if (ikk == 1 && taupt > 20) continue;
      if (ikk == 2 && (taupt < 20 || taupt > 40)) continue;
      if (ikk == 3 && taupt<40) continue;

      if (im) {
        hReco[ikk][im][ih]->Fill( fillreco[ikk] ); // fill modes
        if (im>2) hReco[ikk][1][ih]->Fill( fillreco[ikk] ); // single prong
      }
      hReco[ikk][0][ih]->Fill( fillreco[ikk] ); // all had tau decays 
    }
    //Fill Resolution
    for (int ikk=0; ikk < nres; ikk++){
      if (im) {
        hResolution[ikk][im][ih]->Fill( fillres[ikk] , wtres[ikk] ); // fill modes
        if (im>2) hResolution[ikk][1][ih]->Fill( fillres[ikk], wtres[ikk] ); // single prong
      }
      hResolution[ikk][0][ih]->Fill( fillres[ikk] , wtres[ikk] ); // all had tau decays 
    }
    // fill th2 
    //Fill Reco Taus
    //{"upsilon","upsilon20","upsilon40","upsilon60","taupt", "trackpt", "pi0pt"};
    for (int ikk=0; ikk < nreco; ikk++){
      if (ikk == 1 && taupt > 20) continue;
      if (ikk == 2 && (taupt < 20 || taupt > 40)) continue;
      if (ikk == 3 && taupt<40) continue;

      if (im) {
        hUpsilon[ikk][im][ih]->Fill( upsilon, fillreco[ikk] ); // fill modes
        if (im>2) hUpsilon[ikk][1][ih]->Fill( upsilon, fillreco[ikk] ); // single prong
      }
      hUpsilon[ikk][0][ih]->Fill( upsilon, fillreco[ikk] ); // all had tau decays 
    }


    //if (jentry > 1000) break;
    //long nbreak = 3000000;
    long nbreak = 800000;
    if (nMATCHED % nbreak == 0) {
      out->Write();
      cout<<"----------WRITING FILE-------->"<<endl;
    }
  } // end loop over jentry

  //out->Write();
  out->Write();
  out->Close();
  delete chain; chain=0;
  // Make a log file
  TString outputtxt = "pythialog";
  outputtxt+=sample;
  outputtxt+=".txt";
  ofstream outfile(outputtxt);


  outfile<<"nTAU  "<<nTAU<<endl;
  outfile<<"nSP  "<<nSP<<endl;
  outfile<<"nNP  "<<nNP<<endl;
  outfile<<"nPINU  "<<nPINU<<endl;
  outfile<<"nRHO  "<<nRHO<<endl;
  outfile<<"nPINU_X  "<<nPINU_X<<endl;
  outfile<<"nRHO_X  "<<nRHO_X<<endl;
  outfile<<"nSXP  "<<nSXP<<endl;
  outfile<<"nBADMODE  "<<nBADMODE<<endl;
  outfile<<"nXX  "<<nXX<<endl;
  outfile<<"nKX  "<<nKX<<endl;
  outfile<<"nMATCHED  "<<nMATCHED<<endl;
  outfile<<"nXMATCHED  "<<nXMATCHED<<endl;
  outfile<<"nMATCHEDSP "<<nMATCHEDSP<<endl;


  return 0;
   
}// end main function

