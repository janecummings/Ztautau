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
#include "TauSpinner/Particle.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinnerTool.h"

using namespace std;
using namespace TauSpinner;

float GeV = 1000.;
float kPI = TMath::Pi();
bool DEBUG = false; 

const int nmoments = 8;
TH1F* hMoments[nmoments];
TString moments[nmoments] = {"A4_ZpT", "A_ZpT", "A4_Zmass", "A_Zmass",
                             "A4_ZpT_AZ", "A_ZpT_AZ", "A4_Zmass_AZ", "A_Zmass_AZ"};

int mbins[nmoments] = {200,200,140,140,200,200,140,140};
float mlbins[nmoments] = {0,0,60,60,0,0,60,60};
float mhbins[nmoments] = {200,200,200,200,200,200,200,200};

const int nkin = 5;
const int nhands = 3;
TH1F* hKinematics[nkin][nhands];
TString kinematics[nkin] = {"ZpT", "Zmass", "CosTheta", "ZpT_AZ", "CosTheta_AZ"};
TString hands[nhands] = {"LH","RH","all"};
int kbins[nkin] = {200,140,100,200,100};
float klbins[nkin] = {0,60,-1,0,-1};
float khbins[nkin] = {200,200,1,200,1};

const int nlep = 7;
const int nmodes = 4;
TH1F* hLeptons[nlep][nmodes][nhands];

TString modes[nmodes] = {"0","1","3","4"};
TString leptons[nlep] = {"taupt","leppt","taueta","lepeta","vistaupt","xfrac","upsilon"};
int pbins[nlep] =       {200,     200,    120,      120,    200,       120,    120};
int plbins[nlep] =      {0,       0,      -6,       -6,     0,         0,      -1.2};
int phbins[nlep] =      {200,     200,     6,        6,     200,       1.2,     1.2};


/// weighted distributions
TH1F* hMoments_Zpt[nmoments];
TH1F* hMoments_Zmass[nmoments];
TH1F* hKinematics_Zpt[nkin][nhands];
TH1F* hKinematics_Zmass[nkin][nhands];
TH1F* hLeptons_Zpt[nlep][nmoments][nhands];
TH1F* hLeptons_Zmass[nlep][nmoments][nhands];


const int nsamples = 5;


int main(int argc, char *argv[])
{
  if (argc<2) {
    //cout << " usage : " << endl;
    return 0;
  }

  //TString outFileName = argv[1];
  TString sample = argv[1];
  TString outFileName = "CSAlpgen_Reweight_Fix3.root";
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

  TFile* out = new TFile(outFileName,"UPDATE");
  bool doWeight = true;

  TString h_mname;
  TString h_mwname;
  for(int ihm=0; ihm < nmoments; ihm++){
    h_mname = sample + "_" + moments[ihm];
    hMoments[ihm] = new TH1F( h_mname, h_mname, mbins[ihm], mlbins[ihm], mhbins[ihm]);
    hMoments[ihm]->Sumw2();
    if (doWeight){
      h_mwname = h_mname + "_rwZpt";
      hMoments_Zpt[ihm] = new TH1F( h_mwname, h_mwname, mbins[ihm], mlbins[ihm], mhbins[ihm]);
      hMoments_Zpt[ihm]->Sumw2();
      h_mwname = h_mname + "_rwZmass";
      hMoments_Zmass[ihm] = new TH1F( h_mwname, h_mwname, mbins[ihm], mlbins[ihm], mhbins[ihm]);
      hMoments_Zmass[ihm]->Sumw2();

    }

  }

  TString h_kname;
  TString h_kwname;
  for(int ihk=0; ihk < nkin; ihk++){
    for(int ihh=0;ihh<nhands;ihh++){
      h_kname = sample + "_" + kinematics[ihk] + "_" + hands[ihh];
      hKinematics[ihk][ihh] = new TH1F( h_kname, h_kname, kbins[ihk], klbins[ihk], khbins[ihk]);
      hKinematics[ihk][ihh]->Sumw2();
      if (doWeight){
        h_kwname = h_kname + "_rwZpt";
        hKinematics_Zpt[ihk][ihh] = new TH1F( h_kwname, h_kwname, kbins[ihk], klbins[ihk], khbins[ihk]);
        hKinematics_Zpt[ihk][ihh]->Sumw2();
        h_kwname = h_kname + "_rwZmass";
        hKinematics_Zmass[ihk][ihh] = new TH1F( h_kwname, h_kwname, kbins[ihk], klbins[ihk], khbins[ihk]);
        hKinematics_Zmass[ihk][ihh]->Sumw2();

      }
    }
  }

  TString h_pname;
  TString h_pwname;
  for(int ihp=0; ihp < nlep; ihp++){
    for(int ihm=0;ihm<nmodes;ihm++){
      for(int ihh=0;ihh<nhands;ihh++){
        h_pname = sample + "_" + leptons[ihp] + "_" + modes[ihm] + "_" + hands[ihh];
        hLeptons[ihp][ihm][ihh] = new TH1F( h_pname, h_pname, pbins[ihp], plbins[ihp], phbins[ihp]);
        hLeptons[ihp][ihm][ihh]->Sumw2();
        if (doWeight){
          h_pwname = h_pname + "_rwZpt";
          hLeptons_Zpt[ihp][ihm][ihh] = new TH1F( h_pwname, h_pwname, pbins[ihp], plbins[ihp], phbins[ihp]);
          hLeptons_Zpt[ihp][ihm][ihh]->Sumw2();
          h_pwname = h_pname + "_rwZmass";
          hLeptons_Zmass[ihp][ihm][ihh] = new TH1F( h_pwname, h_pwname, pbins[ihp], plbins[ihp], phbins[ihp]);
          hLeptons_Zmass[ihp][ihm][ihh]->Sumw2();
        }
      }
    }
  }


  gROOT->ProcessLine("#include <vector>");
  SpinTree* spin = new SpinTree(); 
  spin->Init(chain);
  
  int itau = 0;
  
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
  TFile* rwf = new TFile("PythiaWeights.root");
  TH1F *pwts_Zpt, *pwts_Zmass, *awts_Zpt, *awts_Zmass;

  if (doWeight){
  
    pwts_Zpt = (TH1F*)rwf->Get("Pythia_A4_ZpT");
    pwts_Zmass = (TH1F*)rwf->Get("Pythia_A4_Zmass");
    awts_Zpt = (TH1F*)rwf->Get("Alpgen_A4_ZpT_fix");
    awts_Zmass = (TH1F*)rwf->Get("Alpgen_A4_Zmass_fix");
  }
  
  out->cd();

  double hand,weight;

  //loop over events
  for(int jentry = 0; jentry < numEvents; jentry++){
    int ientry = spin->LoadTree(jentry);
    if(ientry<0) break;
    nb = spin->GetEntry(jentry); nbytes += nb;

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
    tool.read_event(*d3pd_mc);
    //get weight -- polari gets filled so can used getTauSpin from tau_reweight_lib
    weight = tool.get_spin_weight();
    hand = tool.get_helicity();
    TauSpinnerHelpers::DecayParticles dp = tool.get_decay_particles();

    TLorentzVector Z = TauSpinnerHelpers::simp2tlv(dp.X);
    TLorentzVector tlv_eneg = TauSpinnerHelpers::simp2tlv(dp.tau);
    TLorentzVector tlv_epos = TauSpinnerHelpers::simp2tlv(dp.tau2);
    TLorentzVector tlv_ee = tlv_epos + tlv_eneg;

    TLorentzVector tlv_beam_pos, tlv_beam_neg;
    float ProtonMass = 0.938; // GeV
    float BeamEnergy =  4000; // GeV
    tlv_beam_pos.SetPxPyPzE(0, 0,  BeamEnergy, BeamEnergy); 
    tlv_beam_neg.SetPxPyPzE(0, 0,- BeamEnergy, BeamEnergy); 

    //components of the cos(th*) calculation
    float Lplus, Lminus, Pplus, Pminus;
    Lplus  = tlv_eneg.E()+tlv_eneg.Pz();
    Lminus = tlv_eneg.E()-tlv_eneg.Pz();
    Pplus  = tlv_epos.E()+tlv_epos.Pz();
    Pminus = tlv_epos.E()-tlv_epos.Pz();

    //does the cos(th*) calculation
    
    TauSpinner::Particle P_tautau( tlv_epos.Px()+tlv_eneg.Px(), tlv_epos.Py()+tlv_eneg.Py(), tlv_epos.Pz()+tlv_eneg.Pz(), tlv_epos.E()+tlv_eneg.E(), 0 );
    float cosThetaCS;
    cosThetaCS  = (Lplus*Pminus - Lminus*Pplus);
    cosThetaCS *= fabs(P_tautau.pz());
    cosThetaCS /= (P_tautau.recalculated_mass()*P_tautau.pz());
    cosThetaCS /= sqrt(P_tautau.recalculated_mass()*P_tautau.recalculated_mass() 
           + ( P_tautau.px()*P_tautau.px() + P_tautau.py()*P_tautau.py() ) );
    float wt_A4 = cosThetaCS;

    if (cosThetaCS != cosThetaCS) {
      cout<<"UH OH!!!"<<cosThetaCS<<endl;
      cout<<" "<<(P_tautau.recalculated_mass()*P_tautau.pz())<<endl;
      cout<<" "<<sqrt(P_tautau.recalculated_mass()*P_tautau.recalculated_mass() 
           + ( P_tautau.px()*P_tautau.px() + P_tautau.py()*P_tautau.py() ) )<<endl; 
      continue;
    }


    hMoments[0]->Fill( tlv_ee.Pt(), wt_A4 );
    hMoments[1]->Fill( tlv_ee.Pt());
    hMoments[2]->Fill( tlv_ee.M(), wt_A4 );
    hMoments[3]->Fill( tlv_ee.M());

    if( tlv_ee.Mag()  > 80.0 && tlv_ee.Mag()  < 100.0 ) {
      hMoments[4]->Fill( tlv_ee.Pt(), wt_A4 );
      hMoments[5]->Fill( tlv_ee.Pt());
      hMoments[6]->Fill( tlv_ee.M(), wt_A4 );
      hMoments[7]->Fill( tlv_ee.M());
    }

    float wt_zpt = 1.;
    float wt_zmass = 1.;
    float PA4,AA4;

    if (doWeight){
      int zbin = pwts_Zpt->FindBin(tlv_ee.Pt());
      if (zbin>10) zbin = 10;
      PA4 = pwts_Zpt->GetBinContent(zbin);
      AA4 = awts_Zpt->GetBinContent(zbin);
      wt_zpt = (1 + wt_A4*wt_A4 + PA4*wt_A4) / (1 + wt_A4*wt_A4 + AA4*wt_A4);
      //cout<<tlv_ee.Pt()<<" "<<zbin<<" "<<PA4<<" "<<AA4<<" "<<wt_zpt<<endl;

      int mbin = pwts_Zmass->FindBin(tlv_ee.M());
      if (mbin>36) mbin = 36;
      PA4 = pwts_Zmass->GetBinContent(mbin);
      AA4 = awts_Zmass->GetBinContent(mbin);
      wt_zmass = (1 + wt_A4*wt_A4 + PA4*wt_A4) / (1 + wt_A4*wt_A4 + AA4*wt_A4);
      //cout<<tlv_ee.M()<<" "<<mbin<<" "<<PA4<<" "<<AA4<<" "<<wt_zmass<<endl;
      //cout<<"*******"<<endl;


      hMoments_Zpt[0]->Fill( tlv_ee.Pt(), wt_A4*wt_zpt );
      hMoments_Zpt[1]->Fill( tlv_ee.Pt(), wt_zpt );
      hMoments_Zpt[2]->Fill( tlv_ee.M(), wt_A4*wt_zpt );
      hMoments_Zpt[3]->Fill( tlv_ee.M(), wt_zpt );

      hMoments_Zmass[0]->Fill( tlv_ee.Pt(), wt_A4*wt_zmass );
      hMoments_Zmass[1]->Fill( tlv_ee.Pt(), wt_zmass );
      hMoments_Zmass[2]->Fill( tlv_ee.M(), wt_A4*wt_zmass );
      hMoments_Zmass[3]->Fill( tlv_ee.M(), wt_zmass );

      if( tlv_ee.Mag()  > 80.0 && tlv_ee.Mag()  < 100.0 ) {
        hMoments_Zpt[4]->Fill( tlv_ee.Pt(), wt_A4*wt_zpt );
        hMoments_Zpt[5]->Fill( tlv_ee.Pt(), wt_zpt );
        hMoments_Zpt[6]->Fill( tlv_ee.M(), wt_A4*wt_zpt );
        hMoments_Zpt[7]->Fill( tlv_ee.M(), wt_zpt );

        hMoments_Zmass[4]->Fill( tlv_ee.Pt(), wt_A4*wt_zmass );
        hMoments_Zmass[5]->Fill( tlv_ee.Pt(), wt_zmass );
        hMoments_Zmass[6]->Fill( tlv_ee.M(), wt_A4*wt_zmass );
        hMoments_Zmass[7]->Fill( tlv_ee.M(), wt_zmass );
      }
    }

    float fkin[nkin];
    fkin[0] = tlv_ee.Pt();
    fkin[1] = tlv_ee.M();
    fkin[2] = wt_A4;
    fkin[3] = tlv_ee.Pt();
    fkin[4] = wt_A4;
    int ih;
    if(hand==-1) ih = 0;
    else if (hand==1) ih = 1;

    for(int ik=0; ik<nkin;ik++){
      hKinematics[ik][ih]->Fill(fkin[ik]);
      hKinematics[ik][2]->Fill(fkin[ik]);
      if (doWeight){
        hKinematics_Zpt[ik][ih]->Fill(fkin[ik], wt_zpt);
        hKinematics_Zpt[ik][2]->Fill(fkin[ik], wt_zpt);
        hKinematics_Zmass[ik][ih]->Fill(fkin[ik], wt_zmass);
        hKinematics_Zmass[ik][2]->Fill(fkin[ik], wt_zmass);
      }
      if (ik>2 && tlv_ee.Mag()  > 80.0 && tlv_ee.Mag()  < 100.0 ){
        hKinematics[ik][ih]->Fill(fkin[ik]);
        hKinematics[ik][2]->Fill(fkin[ik]);
        if (doWeight){
          hKinematics_Zpt[ik][ih]->Fill(fkin[ik], wt_zpt);
          hKinematics_Zpt[ik][2]->Fill(fkin[ik], wt_zpt);
          hKinematics_Zmass[ik][ih]->Fill(fkin[ik], wt_zmass);
          hKinematics_Zmass[ik][2]->Fill(fkin[ik], wt_zmass);
        }
      }
    }


    TLorentzVector tau_had, tau_lep;
    vector<TauSpinner::SimpleParticle> tau_had_daughters;
    vector<TauSpinner::SimpleParticle> tau_lep_daughters;
    int mode = TauSpinnerHelpers::get_decay_mode(dp.tau,dp.tau_daughters);
    int mode2 = TauSpinnerHelpers::get_decay_mode(dp.tau2,dp.tau2_daughters);

    int m = mode;
    tau_had = TauSpinnerHelpers::simp2tlv(dp.tau);
    tau_had_daughters = dp.tau_daughters;
    tau_lep_daughters = dp.tau2_daughters;
    int m2 = 0; 
    if (mode == 1 || mode == 2) { 
      m2 = mode; m = mode2;
      tau_had = TauSpinnerHelpers::simp2tlv(dp.tau2);
      tau_lep = TauSpinnerHelpers::simp2tlv(dp.tau);
      tau_had_daughters = dp.tau2_daughters;
      tau_lep_daughters = dp.tau_daughters;
    }
    else if (mode2 == 1 || mode2 == 2 ) {
      m2 = mode2;
      tau_lep = TauSpinnerHelpers::simp2tlv(dp.tau2);
    }

    // Determine hadronic decay mode 
    int ncharged = 0;
    if (m==5){
      for( int p=0; p < tau_had_daughters.size(); p++){
        if (fabs(tau_had_daughters[p].pdgid()) == 211) ncharged++;
        else if (fabs(tau_had_daughters[p].pdgid()) == 321) ncharged++;
      }
    }
    int im = 0;
    if (ncharged == 1) im = 1; // single prong
    if (m==3 || m==6) im = 2; // single charged pion
    if (m==4 || m == 7) im = 3; // 1p1n

    vector<TauSpinner::SimpleParticle>::iterator id;
    TLorentzVector tau_vis,pi0,pi;
    //loop over daughters and leave out neutrinos
    for (id=tau_had_daughters.begin(); id!=tau_had_daughters.end();id++){
      if ( fabs(id->pdgid()) == 16 || fabs(id->pdgid()) == 22 ) continue;
      tau_vis+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
      if( fabs(id->pdgid()) == 211 || fabs(id->pdgid())==321){
        pi += TauSpinnerHelpers::simp2tlv(*id );
      }
      if( fabs(id->pdgid()) == 111 || fabs(id->pdgid()) == 130 || fabs(id->pdgid()) == 310){
        pi0 += TauSpinnerHelpers::simp2tlv(*id);
      }
    }

    //get lep if there is a lep in the event 
    TLorentzVector lep;
    if (m2){
      for (id=tau_lep_daughters.begin(); id!=tau_lep_daughters.end();id++){
        if ( fabs(id->pdgid()) == 16 || fabs(id->pdgid()) == 22 || fabs(id->pdgid()) == 14 || fabs(id->pdgid()) == 12 ) continue;
        lep+=TLorentzVector(id->px(),id->py(),id->pz(),id->e());
      }
    }

    float flep[nlep];
    flep[0] = tau_had.Pt();
    if (m2) {
      flep[1] = lep.Pt();
      flep[3] = lep.Eta();
    }
    flep[2] = tau_had.Eta();
    flep[4] = tau_vis.Pt();
    flep[5] = (tau_vis.Pt()/tau_had.Pt());
    flep[6] = (pi.Pt() - pi0.Pt())/(pi.Pt() + pi0.Pt());

    for(int ip=0; ip<nlep;ip++){
      if ((ip == 1 || ip == 3) && !m2) continue;
      hLeptons[ip][0][ih]->Fill(flep[ip]); // fill LH/RH all modes
      hLeptons[ip][0][2]->Fill(flep[ip]); // fill all hands all modes
      if (doWeight){
        hLeptons_Zpt[ip][0][ih]->Fill(flep[ip], wt_zpt); // fill LH/RH all modes
        hLeptons_Zpt[ip][0][2]->Fill(flep[ip], wt_zpt); // fill all hands all modes
        hLeptons_Zmass[ip][0][ih]->Fill(flep[ip], wt_zmass); // fill LH/RH all modes
        hLeptons_Zmass[ip][0][2]->Fill(flep[ip], wt_zmass); // fill all hands all modes
      }
      if (im){
        hLeptons[ip][1][ih]->Fill(flep[ip]); // fill LH/RH all single prong
        hLeptons[ip][1][2]->Fill(flep[ip]); // fill all hands all single prong
        if (doWeight){
          hLeptons_Zpt[ip][1][ih]->Fill(flep[ip],wt_zpt); // fill LH/RH all single prong
          hLeptons_Zpt[ip][1][2]->Fill(flep[ip],wt_zpt); // fill all hands all single prong
          hLeptons_Zmass[ip][1][ih]->Fill(flep[ip],wt_zmass); // fill LH/RH all single prong
          hLeptons_Zmass[ip][1][2]->Fill(flep[ip],wt_zmass); // fill all hands all single prong
        }
        if (im>1){
          hLeptons[ip][im][ih]->Fill(flep[ip]); // fill LH/RH mode
          hLeptons[ip][im][2]->Fill(flep[ip]); // fill all hands mode
          if (doWeight){
            hLeptons_Zpt[ip][im][ih]->Fill(flep[ip], wt_zpt); // fill LH/RH mode
            hLeptons_Zpt[ip][im][2]->Fill(flep[ip], wt_zpt); // fill all hands mode
            hLeptons_Zmass[ip][im][ih]->Fill(flep[ip], wt_zmass); // fill LH/RH mode
            hLeptons_Zmass[ip][im][2]->Fill(flep[ip], wt_zmass); // fill all hands mode
          }
        }
      }
    }
    
    //if (jentry > 500 ) break;

    if (jentry % 1000000 == 0 && jentry > 10) {
      out->Write();
      //break;
      cout<<"Wrote "<<itau<<"events for "<<sample<<endl;
    }
    //fillAis( spin, 0, 1., -1);
    //if (jentry % 10000 == 0) cout<<jentry<<endl;
    itau++;


  } // end loop over jentry
  cout<<itau<<endl;

  out->Write();
  out->Close();

  delete chain; chain=0;
  return 0;
   
}// end main function
