// -*- C++ -*-
//
// Package:    Generic_Analizer
// Class:      Generic_Analizer
// 
/**\class Generic_Analizer Generic_Analizer.cc FastTiming/Generic_Analizer/plugins/Generic_Analizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luca Pernie
//         Created:  Fri, 06 Feb 2015 18:30:16 GMT
// $Id$
//
//
// system include files
#include <memory>
#include "TRandom3.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "FastTiming/Generic_Analizer/interface/FastTool.h"
#include "FastTiming/Generic_Analizer/plugins/Generic_Analizer.hh"

using namespace std;
using namespace edm;
using namespace reco;
// ------------------------------------------------------------------------------------------
Generic_Analizer::Generic_Analizer(const edm::ParameterSet& iConfig) {
  fPFCands      = iConfig.getUntrackedParameter<edm::InputTag>("candName");
  SimVtx_       = iConfig.getUntrackedParameter<edm::InputTag>("SimVtx");
  Gamma_        = iConfig.getUntrackedParameter<edm::InputTag>("Gamma");
  Jet_          = iConfig.getUntrackedParameter<edm::InputTag>("Jet");
  GenJet_       = iConfig.getUntrackedParameter<edm::InputTag>("GenJet");
  GenPar_       = iConfig.getUntrackedParameter<edm::InputTag>("GenPar");
  JetCHS_       = iConfig.getUntrackedParameter<edm::InputTag>("JetCHS");
  RecoVtx_      = iConfig.getUntrackedParameter<edm::InputTag>("RecoVtx");
  ak5PFRho_     = iConfig.getUntrackedParameter<edm::InputTag>("ak5PFRho");
  EB_LAYER_     = iConfig.getUntrackedParameter<double>("EB_LAYER",7.5);
  EE_LAYER_     = iConfig.getUntrackedParameter<double>("EE_LAYER",3.5);
  smearing_     = iConfig.getUntrackedParameter<double>("smearing", 0.03);
  isNOPUHgg_    = iConfig.getUntrackedParameter<bool>("isNOPUHgg");
  WannaFitT0Vtx_= iConfig.getUntrackedParameter<bool>("WannaFitT0Vtx");
  OneGamma_     = iConfig.getUntrackedParameter<bool>("OneGamma_");
  Use_MinEne_   = iConfig.getUntrackedParameter<bool>("Use_MinEne");
  Use_R9_       = iConfig.getUntrackedParameter<bool>("Use_R9");
  DoSumEt_      = iConfig.getUntrackedParameter<bool>("DoSumEt");
  DoMass_       = iConfig.getUntrackedParameter<bool>("DoMass");
  isHgg_        = iConfig.getUntrackedParameter<bool>("isHgg");
  SubT0TOF_     = iConfig.getUntrackedParameter<double>("SubT0TOF");
  OutName_      = iConfig.getUntrackedParameter<string>("OutName");
  //Selection
  MinPt_Gen     = iConfig.getUntrackedParameter<double>("MinPt_Gen",30);
  MinPt_Reco    = iConfig.getUntrackedParameter<double>("MinPt_Reco",20);
  MinPt_RecoPu  = iConfig.getUntrackedParameter<double>("MinPt_RecoPu",25);
  MinDR_asso    = iConfig.getUntrackedParameter<double>("MinDR_asso",0.1);
  MinDR_pu      = iConfig.getUntrackedParameter<double>("MinDR_pu",0.6);
  debug = false;
#ifdef DEBUG
  debug = true;
#endif
    cout<<"Time: you choose a SUBTRACTION level "<<SubT0TOF_<<" Where 0 is no Subtraction (your time is perfect), 1 need to sub T0, 2 need to sub TOF, 3 need to sub. T0 and TOF, 4 need to ad back TOF(0,0,0) and sub T0 and TOF."<<endl;
    if(debug) cout<<"DEBUG mode selected."<<endl;
    outfile = new TFile(OutName_.c_str(),"RECREATE");
    outfile->cd();
    h_EventFlow          = new TH1F("h_EventFlow", "", 5, -0.5, 4.5);
    h_EventFlow->GetXaxis()->SetBinLabel(1,"Events"); h_EventFlow->GetXaxis()->SetBinLabel(2,"MC Jets"); h_EventFlow->GetXaxis()->SetBinLabel(3,"RECO Jets"); h_EventFlow->GetXaxis()->SetBinLabel(4,"MC Phot"); h_EventFlow->GetXaxis()->SetBinLabel(5,"RECO Phot");
    h_T0                 = new TH1F("h_T0", "", 100, -0.5, 0.5);
    h_SumEt_cut          = new TH1F("h_SumEt_cut", "", 1000, 0., 2000.);
    h_SumEt_15cut        = new TH1F("h_SumEt_15cut", "", 1000, 0., 2000.);
    h_SumEt_30cut        = new TH1F("h_SumEt_30cut", "", 1000, 0., 2000.);
    h_SumEt_50cut        = new TH1F("h_SumEt_50cut", "", 1000, 0., 2000.);
    h_SumEt_500cut       = new TH1F("h_SumEt_500cut", "", 1000, 0., 2000.);
    h_BestTime_Fir_RemovalSumEt_zoom = new TH1F("h_BestTime_Fir_RemovalSumEt_zoom","", 100, -0.02, 0.2); h_BestTime_Fir_RemovalSumEt_zoom->GetXaxis()->SetTitle("Time [ns]"); 
    h_BestTime_Fir_RemovalSumEt      = new TH1F("h_BestTime_Fir_RemovalSumEt","", 100, -0.3, 0.3); h_BestTime_Fir_RemovalSumEt->GetXaxis()->SetTitle("Time [ns]"); 
    h_SumEt              = new TH1F("h_SumEt", "", 1000, 0., 2000.);
    h_TOT_SumEt_cut      = new TH1F("h_TOT_SumEt_cut", "", 1000, 0., 2000.);
    h_TOT_SumEt_15cut    = new TH1F("h_TOT_SumEt_15cut", "", 1000, 0., 2000.);
    h_TOT_SumEt_30cut    = new TH1F("h_TOT_SumEt_30cut", "", 1000, 0., 2000.);
    h_TOT_SumEt_50cut    = new TH1F("h_TOT_SumEt_50cut", "", 1000, 0., 2000.);
    h_TOT_SumEt_500cut   = new TH1F("h_TOT_SumEt_500cut", "", 1000, 0., 2000.);
    h_TOT_SumEt          = new TH1F("h_TOT_SumEt", "", 1000, 0., 2000.);
    h_Time               = new TH1F("h_Time", "", 100, -0.5, 0.5);
    h_Time_we            = new TH1F("h_Time_we", "", 100, -0.5, 0.5);
    h_TimeSmeared        = new TH1F("h_TimeSmeared", "", 100, -0.5, 0.5);
    h_TimeSmeared_we     = new TH1F("h_TimeSmeared_we", "", 100, -0.5, 0.5);
    h_GoodJet_t          = new TH1F("h_GoodJet_t", "", 100, -0.5, 0.5);
    h_GoodJet_tEB        = new TH1F("h_GoodJet_tEB", "", 100, -0.5, 0.5);
    h_GoodJet_tEE        = new TH1F("h_GoodJet_tEE", "", 100, -0.5, 0.5);
    h_GoodJet_tEB2       = new TH1F("h_GoodJet_tEB2", "", 100, -0.01, 0.03);
    h_GoodJet_tEE2       = new TH1F("h_GoodJet_tEE2", "", 100, -0.01, 0.03);
    h_BadJet_t           = new TH1F("h_BadJet_t", "", 100, -0.5, 0.5);
    h_GoodGamma_t        = new TH1F("h_GoodGamma_t", "", 100, -0.5, 0.5);
    h_GoodGamma_tEB      = new TH1F("h_GoodGamma_tEB", "", 100, -0.5, 0.5);
    h_GoodGamma_tEE      = new TH1F("h_GoodGamma_tEE", "", 100, -0.5, 0.5);
    h_GoodGamma_tEB2     = new TH1F("h_GoodGamma_tEB2", "", 100, -0.01, 0.03);
    h_GoodGamma_tEE2     = new TH1F("h_GoodGamma_tEE2", "", 100, -0.01, 0.03);
    h_Jet_DR             = new TH1F("h_Jet_DR", "", 100, 0., 2.);
    h_Phot_DR            = new TH1F("h_Phot_DR", "", 100, 0., 1.);
    h_PtGenJet           = new TH1F("h_PtGenJet", "", 50, 0., 500.);
    h_EtaGenJet          = new TH1F("h_EtaGenJet", "", 50, 0., 3.);
    h_EffEta_phot1       = new TH1F("h_EffEta_phot1", "", 10, 0., 3.);
    h_EffEta_phot2       = new TH1F("h_EffEta_phot2", "", 10, 0., 3.);
    h_EffEta_phot3       = new TH1F("h_EffEta_phot3", "", 10, 0., 3.);
    h_NEffEta_phot1      = new TH1F("h_NEffEta_phot1", "", 10, 0., 3.);
    h_NEffEta_phot2      = new TH1F("h_NEffEta_phot2", "", 10, 0., 3.);
    h_NEffEta_phot3      = new TH1F("h_NEffEta_phot3", "", 10, 0., 3.);
    h_EffEta_jet1        = new TH1F("h_EffEta_jet1", "", 10, 0., 3.);
    h_EffEta_jet2        = new TH1F("h_EffEta_jet2", "", 10, 0., 3.);
    h_EffEta_jet3        = new TH1F("h_EffEta_jet3", "", 10, 0., 3.);
    h_NEffEta_jet1       = new TH1F("h_NEffEta_jet1", "", 10, 0., 3.);
    h_NEffEta_jet2       = new TH1F("h_NEffEta_jet2", "", 10, 0., 3.);
    h_NEffEta_jet3       = new TH1F("h_NEffEta_jet3", "", 10, 0., 3.);
    h_EffEta_phot        = new TH1F("h_EffEta_phot", "", 10, 0., 3.);
    h_EffEtaTot_phot     = new TH1F("h_EffEtaTot_phot", "", 10, 0., 3.);
    h_EffEta_jet         = new TH1F("h_EffEta_jet", "", 10, 0., 3.);
    h_EffPt_jet          = new TH1F("h_EffPt_jet", "", 100, 0., 500.);
    h_EffEtaTot_jet      = new TH1F("h_EffEtaTot_jet", "", 10, 0., 3.);
    h_EffPtTot_jet       = new TH1F("h_EffPtTot_jet", "", 100, 0., 500.);
    h_BadGamma_t         = new TH1F("h_BadGamma_t", "", 100, -0.5, 0.5);
    h_Rh0                = new TH1F("h_Rh0", "", 100, 0., 100.);
    h_NVtx               = new TH1F("h_NVtx", "", 70, 0., 140.);
    h_HiggsMass          = new TH1F("h_HiggsMass", "", 40, 100., 160.);
    h_HiggsPFMass        = new TH1F("h_HiggsPFMass", "", 40, 100., 160.);
    Ereso1               = new TH1F("Ereso1", "", 50, 0.75, 1.25);
    Ereso2               = new TH1F("Ereso2", "", 50, 0.75, 1.25);
    Ereso1_Mit           = new TH1F("Ereso1_Mit", "", 50, 0.75, 1.25);
    Ereso2_Mit           = new TH1F("Ereso2_Mit", "", 50, 0.75, 1.25);
    h_HiggsPFMass_Mit    = new TH1F("h_HiggsPFMass_Mit", "", 40, 100., 160.);
    h_PU_frac_1          = new TH1F("h_PU_frac_1", "", 100, 0.5, 1.5);
    h_PU_frac_2          = new TH1F("h_PU_frac_2", "", 100, 0.5, 1.5);
    h_HiggsPFMass_Vtx    = new TH1F("h_HiggsPFMass_Vtx", "", 40, 100., 160.);
    h_HiggsMass_MC       = new TH1F("h_HiggsMass_MC", "", 100, 120., 130.);
    Associated = 0.; Associated_EB = 0.;  Associated_EE = 0.; Associated_time = 0.;  Associated_time_EB = 0.;  Associated_time_EE = 0.; Associated_tot = 0.; Associated_tot_EB = 0.; Associated_tot_EE = 0.;
    h_Association        = new TH1F("h_Association", "", 3, -0.5, 2.5);
    h_Association_EB     = new TH1F("h_Association_EB", "", 3, -0.5, 2.5);
    h_Association_EE     = new TH1F("h_Association_EE", "", 3, -0.5, 2.5);
    h_PtGammaAssoEB      = new TH1F("h_PtGammaAssoEB", "", 100, -0.5, 200);
    h_EGammaAssoEB       = new TH1F("h_EGammaAssoEB", "", 100, -0.5, 200);
    h_DRAsso             = new TH1F("h_DRAsso", "", 100, 0., 0.5);
    h_TimeGammaNOTAsso   = new TH1F("h_TimeGammaNOTAsso", "", 100, -0.6, 1.);
    h_TimeGammaNOTAssoW  = new TH1F("h_TimeGammaNOTAssoW", "", 100, -0.6, 1.);
    h_TimeGammaAssoEB_sme15= new TH1F("h_TimeGammaAssoEB_sme15", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEE_sme15= new TH1F("h_TimeGammaAssoEE_sme15", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEB_sme30= new TH1F("h_TimeGammaAssoEB_sme30", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEE_sme30= new TH1F("h_TimeGammaAssoEE_sme30", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEB_sme50= new TH1F("h_TimeGammaAssoEB_sme50", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEE_sme50= new TH1F("h_TimeGammaAssoEE_sme50", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEB_sme500= new TH1F("h_TimeGammaAssoEB_sme500", "", 100, -1.5, 1.5);
    h_TimeGammaAssoEE_sme500= new TH1F("h_TimeGammaAssoEE_sme500", "", 100, -1.5, 1.5);
    h_TimeGammaAssoEB    = new TH1F("h_TimeGammaAssoEB", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEB_L  = new TH1F("h_TimeGammaAssoEB_L", "", 100, -5, 5);
    h_TimeGammaAssoEB_L2 = new TH1F("h_TimeGammaAssoEB_L2", "", 100, -50, 50);
    h_PtGammaAssoEE      = new TH1F("h_PtGammaAssoEE", "", 100, -0.5, 200);
    h_EGammaAssoEE       = new TH1F("h_EGammaAssoEE", "", 100, -0.5, 200);
    h_TimeGammaAssoEE    = new TH1F("h_TimeGammaAssoEE", "", 100, -0.5, 0.5);
    h_TimeGammaAssoEE_L  = new TH1F("h_TimeGammaAssoEE_L", "", 100, -5, 5);
    h_TimeGammaAssoEE_L2 = new TH1F("h_TimeGammaAssoEE_L2", "", 100, -50, 50);
    h_NclustAsso_EB1     = new TH1F("h_NclustAsso_EB1", "", 8, -0.5, 7.5 ); 
    h_NclustAsso_EB2     = new TH1F("h_NclustAsso_EB2", "", 8, -0.5, 7.5 ); 
    h_NclustAsso_EE1     = new TH1F("h_NclustAsso_EE1", "", 8, -0.5, 7.5 ); 
    h_NclustAsso_EE2     = new TH1F("h_NclustAsso_EE2", "", 8, -0.5, 7.5 ); 
    h_DR_vs_Time_EB      = new TH2F("h_DR_vs_Time_EB", "", 5, 0., 0.0855, 10000, -0.2, 100. ); 
    h_DR_vs_Time_EB2     = new TH2F("h_DR_vs_Time_EB2", "", 100, 0., 0.0855, 100, -0.1, 1. ); 
    h_DR_vs_Time_EB_reb  = new TH2F("h_DR_vs_Time_EB_reb", "", 5, 0., 0.0855, 1100, -100., 1000. ); 
    h_DR_vs_Time_EB_reb2 = new TH2F("h_DR_vs_Time_EB_reb2", "", 10, 0., 0.0855, 1000, 0., 1. ); 
    h_DR_vs_Time_EB_b    = new TH2F("h_DR_vs_Time_EB_b", "", 5, 0., 0.0855, 10000, -0.2, 100. ); 
    h_DR_vs_Time_EB_reb_b= new TH2F("h_DR_vs_Time_EB_reb_b", "", 5, 0., 0.0855, 1100, -100., 1000. ); 
    h_DR_vs_Time_EB_reb2_b= new TH2F("h_DR_vs_Time_EB_reb2_b", "", 10, 0., 0.0855, 1000, 0., 1. ); 
    h_DR_vs_Time_L_EB    = new TH2F("h_DR_vs_Time_L_EB", "", 50, 0., 0.05, 1000, -1.1, 2.1 ); 
    h_DR_vs_Time_L_EB_b  = new TH2F("h_DR_vs_Time_L_EB_b", "", 50, 0., 0.05, 1000, -1.1, 2.1 ); 
    h_energyForDR_EB     = new TH1F("h_energyForDR_EB", "", 100, -0.4, 3.5 ); 
    h_timeForDR_EB       = new TH1F("h_timeForDR_EB", "", 1000, -5., 5 ); 
    h_timevsEne_EB       = new TH2F("h_timevsEne_EB", "", 100, -2., 2, 200, -1., 5 ); h_timevsEne_EB->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EB->GetYaxis()->SetTitle("Energy");
    h_timevsEne_EB2      = new TH2F("h_timevsEne_EB2", "", 1000, -10., 100, 200, -1., 5 ); h_timevsEne_EB2->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EB2->GetYaxis()->SetTitle("Energy");
    h_energyForDR_EB_b   = new TH1F("h_energyForDR_EB_b", "", 100, -0.4, 3.5 ); 
    h_timeForDR_EB_b     = new TH1F("h_timeForDR_EB_b", "", 1000, -5., 5 ); 
    h_timevsEne_EB_b     = new TH2F("h_timevsEne_EB_b", "", 100, -2., 2, 200, -1., 5 ); h_timevsEne_EB->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EB->GetYaxis()->SetTitle("Energy");
    h_timevsEne_EB2_b    = new TH2F("h_timevsEne_EB2_b", "", 1000, -10., 100, 200, -1., 5 ); h_timevsEne_EB2->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EB2->GetYaxis()->SetTitle("Energy");
    h_MEANRMS_vs_DR_EB   = new TH1F("h_MEANRMS_vs_DR_EB", "", 5, 0., 0.0855 ); h_MEANRMS_vs_DR_EB->GetXaxis()->SetTitle("#Delta R"); h_MEANRMS_vs_DR_EB->GetYaxis()->SetTitle("Mean Time (with RMS)");
    h_MEANRMS_vs_DR_EB_c = new TH1F("h_MEANRMS_vs_DR_EB_c", "", 5, 0., 0.0855 ); h_MEANRMS_vs_DR_EB_c->GetXaxis()->SetTitle("#Delta R"); h_MEANRMS_vs_DR_EB_c->GetYaxis()->SetTitle("Mean Time (with RMS)");
    h_DR_vs_Time_EE      = new TH2F("h_DR_vs_Time_EE", "", 5, 0., 0.0855, 10000, -0.2, 100. ); 
    h_DR_vs_Time_EE2     = new TH2F("h_DR_vs_Time_EE2", "", 100, 0., 0.0855, 100, -0.1, 1. ); 
    h_DR_vs_Time_EE_reb  = new TH2F("h_DR_vs_Time_EE_reb", "", 5, 0., 0.0855, 1100, -100., 1000. ); 
    h_DR_vs_Time_EE_reb2 = new TH2F("h_DR_vs_Time_EE_reb2", "", 10, 0., 0.0855, 1000, 0., 1. ); 
    h_DR_vs_Time_EE_b    = new TH2F("h_DR_vs_Time_EE_b", "", 5, 0., 0.0855, 10000, -0.2, 100. ); 
    h_DR_vs_Time_EE_reb_b= new TH2F("h_DR_vs_Time_EE_reb_b", "", 5, 0., 0.0855, 1100, -100., 1000. ); 
    h_DR_vs_Time_EE_reb2_b= new TH2F("h_DR_vs_Time_EE_reb2_b", "", 10, 0., 0.0855, 1000, 0., 1. ); 
    h_DR_vs_Time_L_EE    = new TH2F("h_DR_vs_Time_L_EE", "", 50, 0., 0.05, 1000, -1.1, 2.1 ); 
    h_DR_vs_Time_L_EE_b  = new TH2F("h_DR_vs_Time_L_EE_b", "", 50, 0., 0.05, 1000, -1.1, 2.1 ); 
    h_energyForDR_EE     = new TH1F("h_energyForDR_EE", "", 100, -0.4, 3.5 ); 
    h_timeForDR_EE       = new TH1F("h_timeForDR_EE", "", 1000, -5., 5 ); 
    h_timevsEne_EE       = new TH2F("h_timevsEne_EE", "", 100, -2., 2, 200, -1., 5 ); h_timevsEne_EE->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EE->GetYaxis()->SetTitle("Energy");
    h_timevsEne_EE2      = new TH2F("h_timevsEne_EE2", "", 1000, -10., 100, 200, -1., 5 ); h_timevsEne_EE2->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EE2->GetYaxis()->SetTitle("Energy");
    h_energyForDR_EE_b   = new TH1F("h_energyForDR_EE_b", "", 100, -0.4, 3.5 ); 
    h_timeForDR_EE_b     = new TH1F("h_timeForDR_EE_b", "", 1000, -5., 5 ); 
    h_timevsEne_EE_b     = new TH2F("h_timevsEne_EE_b", "", 100, -2., 2, 200, -1., 5 ); h_timevsEne_EE->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EE->GetYaxis()->SetTitle("Energy");
    h_timevsEne_EE2_b    = new TH2F("h_timevsEne_EE2_b", "", 1000, -10., 100, 200, -1., 5 ); h_timevsEne_EE2->GetXaxis()->SetTitle("Time [ns]"); h_timevsEne_EE2->GetYaxis()->SetTitle("Energy");
    h_MEANRMS_vs_DR_EE   = new TH1F("h_MEANRMS_vs_DR_EE", "", 5, 0., 0.0855 ); h_MEANRMS_vs_DR_EE->GetXaxis()->SetTitle("#Delta R"); h_MEANRMS_vs_DR_EE->GetYaxis()->SetTitle("Mean Time (with RMS)");
    h_MEANRMS_vs_DR_EE_c = new TH1F("h_MEANRMS_vs_DR_EE_c", "", 5, 0., 0.0855 ); h_MEANRMS_vs_DR_EE_c->GetXaxis()->SetTitle("#Delta R"); h_MEANRMS_vs_DR_EE_c->GetYaxis()->SetTitle("Mean Time (with RMS)");
    h_R91                = new TH1F("h_R91","", 100, 0., 2.); h_R91->GetXaxis()->SetTitle("R9"); 
    h_R92                = new TH1F("h_R92","", 100, 0., 2.); h_R92->GetXaxis()->SetTitle("R9"); 
    if( WannaFitT0Vtx_ ){
      Tree_Vtx = new TTree("Tree_Vtx","TTree for vertexing");
      Tree_Vtx->Branch("VtxDet_T1",&VtxDet_T1,"VtxDet_T1/F");
      Tree_Vtx->Branch("VtxDet_VBFT1",&VtxDet_VBFT1,"VtxDet_VBFT1/F");
      Tree_Vtx->Branch("VtxDet_GT1",&VtxDet_GT1,"VtxDet_GT1/F");
      Tree_Vtx->Branch("VtxDet_VBFGT1",&VtxDet_VBFGT1,"VtxDet_VBFGT1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_X1",&VtxDet_PosXtal_X1,"VtxDet_PosXtal_X1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_Y1",&VtxDet_PosXtal_Y1,"VtxDet_PosXtal_Y1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_Z1",&VtxDet_PosXtal_Z1,"VtxDet_PosXtal_Z1/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_X1",&VtxDet_PosVBF_X1,"VtxDet_PosVBF_X1/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_Y1",&VtxDet_PosVBF_Y1,"VtxDet_PosVBF_Y1/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_Z1",&VtxDet_PosVBF_Z1,"VtxDet_PosVBF_Z1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCX1",&VtxDet_PosXtal_MCX1,"VtxDet_PosXtal_MCX1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCY1",&VtxDet_PosXtal_MCY1,"VtxDet_PosXtal_MCY1/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCZ1",&VtxDet_PosXtal_MCZ1,"VtxDet_PosXtal_MCZ1/F");
      Tree_Vtx->Branch("VtxDet_PtRecoJet_1",&VtxDet_PtRecoJet_1,"VtxDet_PtRecoJet_1/F");
      Tree_Vtx->Branch("VtxDet_PtMCJet_1",&VtxDet_PtMCJet_1,"VtxDet_PtMCJet_1/F");
      Tree_Vtx->Branch("VtxDet_time",&VtxDet_time,"VtxDet_time/F");
      Tree_Vtx->Branch("VtxDet_T2",&VtxDet_T2,"VtxDet_T2/F");
      Tree_Vtx->Branch("VtxDet_VBFT2",&VtxDet_VBFT2,"VtxDet_VBFT2/F");
      Tree_Vtx->Branch("VtxDet_GT2",&VtxDet_GT2,"VtxDet_GT2/F");
      Tree_Vtx->Branch("VtxDet_VBFGT2",&VtxDet_VBFGT2,"VtxDet_VBFGT2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_X2",&VtxDet_PosXtal_X2,"VtxDet_PosXtal_X2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_Y2",&VtxDet_PosXtal_Y2,"VtxDet_PosXtal_Y2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_Z2",&VtxDet_PosXtal_Z2,"VtxDet_PosXtal_Z2/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_X2",&VtxDet_PosVBF_X2,"VtxDet_PosVBF_X2/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_Y2",&VtxDet_PosVBF_Y2,"VtxDet_PosVBF_Y2/F");
      Tree_Vtx->Branch("VtxDet_PosVBF_Z2",&VtxDet_PosVBF_Z2,"VtxDet_PosVBF_Z2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCX2",&VtxDet_PosXtal_MCX2,"VtxDet_PosXtal_MCX2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCY2",&VtxDet_PosXtal_MCY2,"VtxDet_PosXtal_MCY2/F");
      Tree_Vtx->Branch("VtxDet_PosXtal_MCZ2",&VtxDet_PosXtal_MCZ2,"VtxDet_PosXtal_MCZ2/F");
      Tree_Vtx->Branch("VtxDet_PtRecoJet_2",&VtxDet_PtRecoJet_2,"VtxDet_PtRecoJet_2/F");
      Tree_Vtx->Branch("VtxDet_PtMCJet_2",&VtxDet_PtMCJet_2,"VtxDet_PtMCJet_2/F");
      Tree_Vtx->Branch("vzMC",&vzMC,"vzMC/F");
      Tree_Vtx->Branch("vyMC",&vyMC,"vyMC/F");
      Tree_Vtx->Branch("vxMC",&vxMC,"vxMC/F");
      Tree_Vtx->Branch("vzRECO",&vzRECO,"vzRECO/F");
      Tree_Vtx->Branch("vyRECO",&vyRECO,"vyRECO/F");
      Tree_Vtx->Branch("vxRECO",&vxRECO,"vxRECO/F");
    }
  ranGaus_ = new TRandom(0);
  FTool_ = new FastTool( );
}
// ------------------------------------------------------------------------------------------
Generic_Analizer::~Generic_Analizer(){
    outfile->Write();
    outfile->Close();
  delete ranGaus_;
  delete FTool_;
}
// ------------------------------------------------------------------------------------------
void Generic_Analizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //-------------------------------------------------LOAD Collections--------------------------------------------------------------
  // PFCandidate
  edm::Handle<reco::PFCandidateCollection> hPFProduct;
  iEvent.getByLabel(fPFCands, hPFProduct);
  const reco::PFCandidateCollection *PFCol = 0;
  reco::PFCandidateCollection tmpColl;
  if (hPFProduct.isValid())
   PFCol = hPFProduct.product();
  else {
    edm::Handle< std::vector<edm::FwdPtr<reco::PFCandidate> > > hPFFwdProduct; 
    iEvent.getByLabel(fPFCands, hPFFwdProduct);

    for (auto & fwdCand : *hPFFwdProduct) {
      tmpColl.push_back(*fwdCand);
    }
    PFCol = & tmpColl;
  }
  // Geometry
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  geometry_ = geoHandle.product();
  // RecHit barrel
  edm::Handle<EcalRecHitCollection> recHitsEB;
  iEvent.getByLabel(  edm::InputTag("ecalDetailedTimeRecHit","EcalRecHitsEB"), recHitsEB );
  if ( ! recHitsEB.isValid() ) {
    edm::LogWarning("EBRecoSummary") << "recHitsEB not found"; 
  }
  // RecHit endcap
  edm::Handle<EcalRecHitCollection> recHitsEE;
  iEvent.getByLabel( edm::InputTag("ecalDetailedTimeRecHit","EcalRecHitsEK"), recHitsEE );
  if ( ! recHitsEE.isValid() ) {
    edm::LogWarning("EERecoSummary") << "recHitsEE not found"; 
  }
  // Photons
  edm::Handle<std::vector<reco::Photon>> PhotonsProd;
  iEvent.getByLabel(Gamma_, PhotonsProd);
  const reco::PhotonCollection *Photons = 0;
  if ( ! PhotonsProd.isValid() ) {
    edm::LogWarning("PhotonsSummary") << "Photons not found";
  }
  Photons = PhotonsProd.product();
  // Jet
  edm::Handle<std::vector<reco::PFJet>> JetProd;
  iEvent.getByLabel(Jet_, JetProd);
  const reco::PFJetCollection *Jets = 0;
  if ( ! JetProd.isValid() ) {
    edm::LogWarning("JetSummary") << "Jets not found";
  }
  Jets = JetProd.product();
  // RecoVertex
  edm::Handle<std::vector<reco::Vertex>> RecoVtx;
  iEvent.getByLabel(RecoVtx_, RecoVtx);
  if ( ! RecoVtx.isValid() ) {
    edm::LogWarning("RecoVtxSummary") << "RecoVtx not found";
  }
  std::vector<reco::Vertex>::const_iterator recoVtxFirst = RecoVtx->begin();
  GlobalPoint Vtx_reco( recoVtxFirst->position().x(), recoVtxFirst->position().y(), recoVtxFirst->position().z() );
  h_NVtx->Fill( RecoVtx->size() );
  // GenJet
  edm::Handle<std::vector<reco::GenJet>> GenJetProd;
  iEvent.getByLabel(GenJet_, GenJetProd);
  const reco::GenJetCollection *GenJets = 0;
  if ( ! GenJetProd.isValid() ) {
    edm::LogWarning("GenJetSummary") << "GenJets not found";
  }
  GenJets = GenJetProd.product();
  // GenParticle
  edm::Handle<std::vector<reco::GenParticle>> GenParProd;
  iEvent.getByLabel(GenPar_, GenParProd);
  const reco::GenParticleCollection *GenPars = 0;
  if ( ! GenParProd.isValid() ) {
    edm::LogWarning("GenParSummary") << "GenPars not found";
  }
  GenPars = GenParProd.product();
  // SimVertex
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel(SimVtx_, SimVtx);
  if ( ! SimVtx.isValid() ) {
    edm::LogWarning("SimVtxSummary") << "SimVtx not found";
  }
  FTool_->Inizialization( SimVtx, EB_LAYER_, EE_LAYER_ );
  GlobalPoint Vtx_sim( FTool_->GiveVtxX(), FTool_->GiveVtxY(), FTool_->GiveVtxZ() );
  float T0_Vtx_MC = FTool_->GiveT0();
  h_T0->Fill( T0_Vtx_MC );

  //----------------------------------------Let's START--------------------------------------------------
  h_EventFlow->Fill(0);
  //Good Jets list
  vector<const reco::PFJet*> GoodJetList, BadJetList; GoodJetList.clear(); BadJetList.clear();
  float maxGenJetPt=0., maxGenJetEta=0;
  for (auto& pfGenJet : *GenJets){
    if( pfGenJet.p4().Pt() > maxGenJetPt ) { maxGenJetPt = pfGenJet.p4().Pt(); maxGenJetEta = fabs(pfGenJet.p4().Eta()); }
    if( pfGenJet.p4().Pt() < MinPt_Gen ) continue;
    h_EventFlow->Fill(1);
    float DR_min = MinDR_asso, DR_minH=99.;
    GlobalPoint PosGenJet( pfGenJet.p4().X(), pfGenJet.p4().Y(), pfGenJet.p4().Z() );
    h_EffEtaTot_jet->Fill( fabs(PosGenJet.eta()) );
    h_EffPtTot_jet->Fill( pfGenJet.p4().Pt() );
    const reco::PFJet* GoodJet = 0;
    for (auto& pfJet : *Jets){
      if( pfJet.p4().Pt() < MinPt_Reco ) continue;
	GlobalPoint PosJet( pfJet.p4().X(), pfJet.p4().Y(), pfJet.p4().Z() );
	float DR = DeltaR( PosJet, PosGenJet );
	if( DR<DR_min ){
	  DR_min = DR;
	  GoodJet = &pfJet;
	  h_EventFlow->Fill(2);
	}
	if( DR<DR_minH ){
	  DR_minH = DR;
	}
    }
    h_Jet_DR->Fill( DR_minH );
    if( GoodJet!=0 ) GoodJetList.push_back( GoodJet );
  }
  h_PtGenJet->Fill( maxGenJetPt );
  h_EtaGenJet->Fill( maxGenJetEta );
  //Bad Jets list
  for (auto& pfJet : *Jets){
    if( sqrt( pow(pfJet.p4().Px(),2)+pow(pfJet.p4().Py(),2) ) < MinPt_RecoPu ) continue;
    GlobalPoint PosJet( pfJet.p4().X(),  pfJet.p4().Y(), pfJet.p4().Z() );
    bool isPU = true;
    for (auto& pfGenJet : *GenJets){
	GlobalPoint PosGenJet( pfGenJet.p4().X(), pfGenJet.p4().Y(), pfGenJet.p4().Z() );
	float DR = DeltaR( PosJet, PosGenJet );
	if( DR<MinDR_pu ){
	  isPU = false;
	}
    }
    if( isPU ){
	for (auto& GenPar : *GenPars){
	  GlobalPoint PosGenP( GenPar.p4().X(), GenPar.p4().Y(), GenPar.p4().Z() );
	  float DR = DeltaR( PosJet, PosGenP );
	  if( DR<MinDR_pu ){
	    isPU = false;
	  }
	}
    }
    if( isPU ){
	const reco::PFJet* thisJet;
	thisJet = &pfJet;
	BadJetList.push_back( thisJet );
    }
  }
  //Plot with good and bad jets
  for(int i=0; i<int(GoodJetList.size()); i++){
    vector<float> TimeXYZ = GetTimeFromJet( GoodJetList[i], recHitsEB, recHitsEE, geometry_, FTool_ );
    float time = TimeXYZ[0];
    GlobalPoint thisPos( TimeXYZ[1], TimeXYZ[2], TimeXYZ[3] );
    GlobalPoint Zero( 0., 0., 0. );
    if( SubT0TOF_==1 )  time-=T0_Vtx_MC; 
    if( SubT0TOF_==2 ){ time-=computeTOF( Vtx_sim,thisPos ); }
    if( SubT0TOF_==3 ){ time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim,thisPos ); }
    if( SubT0TOF_==4 ){ time+=computeTOF( Zero, thisPos ); time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    h_GoodJet_t->Fill( time );
    if( fabs(GoodJetList[i]->p4().eta())<1.47 ){ h_GoodJet_tEB->Fill( time ); h_GoodJet_tEB2->Fill( time );}
    if( fabs(GoodJetList[i]->p4().eta())>1.50 ){ h_GoodJet_tEE->Fill( time ); h_GoodJet_tEE2->Fill( time );}
    h_EffEta_jet->Fill( fabs(GoodJetList[i]->p4().eta()) );
    h_EffPt_jet->Fill( GoodJetList[i]->pt() );
    if( time>-0.5 && time < 0.5 )  h_EffEta_jet1->Fill( fabs(GoodJetList[i]->p4().eta()) );
    if( time> 0.0 && time < 0.1 )  h_EffEta_jet2->Fill( fabs(GoodJetList[i]->p4().eta()) );
    if( time>-200 && time < 200. ) h_EffEta_jet3->Fill( fabs(GoodJetList[i]->p4().eta()) );
  }
  for(int i=0; i<int(BadJetList.size()); i++){
    vector<float> TimeXYZ = GetTimeFromJet( BadJetList[i], recHitsEB, recHitsEE, geometry_, FTool_ );
    float time = TimeXYZ[0];
    GlobalPoint thisPos( TimeXYZ[1], TimeXYZ[2], TimeXYZ[3] );
    GlobalPoint Zero( 0., 0., 0. );
    if( SubT0TOF_==1 )  time-=T0_Vtx_MC;
    if( SubT0TOF_==2 ){ time-=computeTOF( Vtx_sim,thisPos ); }
    if( SubT0TOF_==3 ){ time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim,thisPos ); }
    if( SubT0TOF_==4 ){ time+=computeTOF( Zero, thisPos ); time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    h_BadJet_t->Fill( time );
    if( time>-0.5 && time < 0.5 )  h_NEffEta_jet1->Fill( fabs(BadJetList[i]->p4().eta()) );
    if( time> 0.0 && time < 0.1 )  h_NEffEta_jet2->Fill( fabs(BadJetList[i]->p4().eta()) );
    if( time>-200 && time < 200. ) h_NEffEta_jet3->Fill( fabs(BadJetList[i]->p4().eta()) );
  }

  //Good Photons list
  vector<const reco::Photon*> BadPhotList; BadPhotList.clear();
  vector<const reco::PFJet*>  GoodPhotList;GoodPhotList.clear();
  for (auto& GenPar : *GenPars){
    //if( GenPar.p4().Pt() < MinPt_Gen || GenPar.pdgId()!=22 || GenPar.status()!=1 ) continue;
    if( GenPar.pdgId()!=22 || GenPar.status()!=1 ) continue;
    h_EventFlow->Fill(3);
    float DR_min = MinDR_asso, DR_minH=99.;
    const reco::PFJet* GoodPhot = 0;
    GlobalPoint PosGen( GenPar.p4().X(), GenPar.p4().Y(), GenPar.p4().Z() );
    h_EffEtaTot_phot->Fill( fabs(GenPar.eta()) );
    for (auto& pfJet : *Jets){
      if( pfJet.p4().Pt() < MinPt_Reco ) continue;
	GlobalPoint PosJet( pfJet.p4().X(), pfJet.p4().Y(), pfJet.p4().Z() );
	float DR = DeltaR( PosJet, PosGen );
	if( DR<DR_min ){
	  DR_min = DR;
	  GoodPhot = &pfJet;
	  h_EventFlow->Fill(4);
	}
	if( DR<DR_minH ){
	  DR_minH = DR;
	}
    }
    h_Phot_DR->Fill( DR_minH ); 
    if( GoodPhot!=0 ) GoodPhotList.push_back( GoodPhot );
  }
  //Bad Gamma
  for (auto& photon : *Photons){
    if( photon.p4().Pt() < MinPt_RecoPu ) continue;
    GlobalPoint PosPhot( photon.p4().X(),  photon.p4().Y(), photon.p4().Z() );
    bool isPU = true;
    for (auto& GenPar : *GenPars){
	GlobalPoint PosGen( GenPar.p4().X(), GenPar.p4().Y(), GenPar.p4().Z() );
	float DR = DeltaR( PosPhot, PosGen );
	if( DR<MinDR_pu ){
	  isPU = false;
	}
    }
    if( isPU ){
	const reco::Photon* thisPhoto;
	thisPhoto = &photon;
	BadPhotList.push_back( thisPhoto );
    }
  }
  //Plots with Good and Bad photons
  for(int i=0; i<int(GoodPhotList.size()); i++){
    vector<float> TimeXYZ = GetTimeFromJet( GoodPhotList[i], recHitsEB, recHitsEE, geometry_, FTool_ );
    float time = TimeXYZ[0];
    GlobalPoint thisPos( TimeXYZ[1], TimeXYZ[2], TimeXYZ[3] );
    GlobalPoint Zero( 0.,0.,0. );
    if( SubT0TOF_==1 )  time-=T0_Vtx_MC;
    if( SubT0TOF_==2 ){ time-=computeTOF( Vtx_sim, thisPos ); }
    if( SubT0TOF_==3 ){ time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    if( SubT0TOF_==4 ){ time+=computeTOF( Zero, thisPos ); time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    h_GoodGamma_t->Fill( time );
    if( fabs(GoodPhotList[i]->p4().eta())<1.47 ){ h_GoodGamma_tEB->Fill( time ); h_GoodGamma_tEB2->Fill( time );}
    if( fabs(GoodPhotList[i]->p4().eta())>1.50 ){ h_GoodGamma_tEE->Fill( time ); h_GoodGamma_tEE2->Fill( time );}
    h_EffEta_phot->Fill( fabs(GoodPhotList[i]->p4().eta()) );
    if( time>-0.5 && time < 0.5 )  h_EffEta_phot1->Fill( fabs(GoodPhotList[i]->p4().eta()) );
    if( time> 0.0 && time < 0.1 )  h_EffEta_phot2->Fill( fabs(GoodPhotList[i]->p4().eta()) );
    if( time>-200 && time < 200. ) h_EffEta_phot3->Fill( fabs(GoodPhotList[i]->p4().eta()) );
  }
  for(int i=0; i<int(BadPhotList.size()); i++){
    vector<float> TimeXYZ = GetTimeFromGamma( BadPhotList[i], recHitsEB, recHitsEE, geometry_, FTool_ );
    float time = TimeXYZ[0];
    GlobalPoint thisPos( TimeXYZ[1], TimeXYZ[2], TimeXYZ[3] );
    GlobalPoint Zero( 0.,0.,0. );
    if( SubT0TOF_==1 )  time-=T0_Vtx_MC;
    if( SubT0TOF_==2 ){ time-=computeTOF( Vtx_sim, thisPos ); }
    if( SubT0TOF_==3 ){ time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    if( SubT0TOF_==4 ){ time+=computeTOF( Zero, thisPos ); time-=T0_Vtx_MC; time-=computeTOF( Vtx_sim, thisPos ); }
    h_BadGamma_t->Fill( time );
    if( time>-0.5 && time < 0.5 )  h_NEffEta_phot1->Fill( fabs(BadPhotList[i]->p4().eta()) );
    if( time> 0.0 && time < 0.1 )  h_NEffEta_phot2->Fill( fabs(BadPhotList[i]->p4().eta()) );
    if( time>-200 && time < 200. ) h_NEffEta_phot3->Fill( fabs(BadPhotList[i]->p4().eta()) );
  }
  //Hgg MC Photons
  TLorentzVector Gamma1, Gamma2; Gamma1.SetPtEtaPhiE( -1., -1., -1., -1. ); Gamma2.SetPtEtaPhiE( -1., -1., -1., -1. );
  TLorentzVector VBF1, VBF2; VBF1.SetPtEtaPhiE( -1., -1., -1., -1. ); VBF2.SetPtEtaPhiE( -1., -1., -1., -1. );
  bool MC_pres=false;
  if( isHgg_ ){
    bool firstnotfound = true;
    //VBG Jets
    float tmpPt = 20.;
    for (auto& pfGenJet : *GenJets){
	if( pfGenJet.p4().Pt()<20. ) continue;
	if( pfGenJet.p4().Pt() > tmpPt ){
	  tmpPt = pfGenJet.p4().Pt();
	  VBF1.SetPtEtaPhiE( pfGenJet.p4().Pt(), pfGenJet.p4().Eta(), pfGenJet.p4().Phi(), pfGenJet.p4().E() );
	}
    }
    tmpPt = 20.;
    for (auto& pfGenJet : *GenJets){
	if( pfGenJet.p4().Pt()<20. || pfGenJet.p4().Pt() == VBF1.Pt() ) continue;
	if( pfGenJet.p4().Pt() > tmpPt ){
	  tmpPt = pfGenJet.p4().Pt();
	  VBF2.SetPtEtaPhiE( pfGenJet.p4().Pt(), pfGenJet.p4().Eta(), pfGenJet.p4().Phi(), pfGenJet.p4().E() );
	}
    }
    //GAMMAS FROM H
    for (auto& GenPar : *GenPars){
	if( GenPar.pdgId()==22 && GenPar.mother()->pdgId()==25 && firstnotfound ){
	  Gamma1.SetPtEtaPhiE( GenPar.pt(), GenPar.p4().Eta(), GenPar.p4().Phi(), GenPar.p4().E() );
	  firstnotfound = false;
	  continue;
	}
	if( GenPar.pdgId()==22 && GenPar.mother()->pdgId()==25 && GenPar.p4().Eta() != Gamma1.Eta() ){
	  Gamma2.SetPtEtaPhiE( GenPar.pt(), GenPar.p4().Eta(), GenPar.p4().Phi(), GenPar.p4().E() );
	}
    }
    MC_pres = Gamma1.Pt()>15. && Gamma2.Pt()>15. && fabs(Gamma1.Eta())<2.5 && fabs(Gamma2.Eta())<2.5;
    if( MC_pres ){
	h_HiggsMass_MC->Fill( (Gamma1+Gamma2).M() );
	Associated_tot+=2;
	if( fabs(Gamma1.Eta()) < 1.48) Associated_tot_EB++;
	if( fabs(Gamma1.Eta()) > 1.48) Associated_tot_EE++;
	if( fabs(Gamma2.Eta()) < 1.48) Associated_tot_EB++;
	if( fabs(Gamma2.Eta()) > 1.48) Associated_tot_EE++;
    }
  }
  //SumEt Plot
  if( DoSumEt_ && isHgg_ ){
    float SumEt_tot=0., SumEt_cutted_tot=0.,SumEt_15cutted_tot=0., SumEt_30cutted_tot=0., SumEt_50cutted_tot=0., SumEt_500cutted_tot=0.;
    float TOT_SumEt_tot=0., TOT_SumEt_cutted_tot=0., TOT_SumEt_15cutted_tot=0., TOT_SumEt_30cutted_tot=0., TOT_SumEt_50cutted_tot=0., TOT_SumEt_500cutted_tot=0.;
    bool blabla=false;
    bool blabla2=false;
    std::vector<TLorentzVector> Associated_v1, Associated_v2; Associated_v1.clear(); Associated_v2.clear();
    std::vector<TLorentzVector> Associated_v1_noTime, Associated_v2_noTime; Associated_v1_noTime.clear(); Associated_v2_noTime.clear();
    std::vector<float> AssociatedTime_v1_noTime, AssociatedTime_v2_noTime; AssociatedTime_v1_noTime.clear(); AssociatedTime_v2_noTime.clear();
    std::vector<float> AssociatedDR_v1_noTime, AssociatedDR_v2_noTime; AssociatedDR_v1_noTime.clear(); AssociatedDR_v2_noTime.clear();
    if(blabla || blabla2) cout<<"START EVENT:-------------"<<endl;
    if( blabla2 && MC_pres ) cout<<" The MC Gamma1 are: Pt "<<Gamma1.Pt()<<" Eta "<<Gamma1.Eta()<<" Phi "<<Gamma1.Phi()<<endl;
    if( blabla2 && MC_pres ) cout<<" The MC Gamma2 are: Pt "<<Gamma2.Pt()<<" Eta "<<Gamma2.Eta()<<" Phi "<<Gamma2.Phi()<<endl;
    vector<reco::PFClusterRef> ClusetAlreadyUsed; ClusetAlreadyUsed.clear();
    for (auto& pfCand : *PFCol)
    {
	std::vector<EcalRecHit> V_seeds; V_seeds.clear(); std::vector<reco::PFClusterRef> V_cluster; V_cluster.clear(); std::vector<GlobalPoint> V_seedsPos; V_seedsPos.clear();
	if( pfCand.particleId() != reco::PFCandidate::gamma ) continue; 
	if(blabla) cout<<"--START PFCAND: "<<endl;
	const EcalRecHit* seedHit = 0;
	GlobalPoint seedHitPos;
	float SumEt_pf = 0, SumEt_cutted_pf = 0, SumEt_15cutted_pf =0., SumEt_30cutted_pf =0., SumEt_50cutted_pf =0., SumEt_500cutted_pf =0.;
	TLorentzVector TotSum;
	TotSum.SetPtEtaPhiE( pfCand.energy()/cosh( pfCand.eta() ), pfCand.eta(), pfCand.phi(), pfCand.energy() );
	float Total_SumEt_pf = TotSum.Et(), Total_SumEt_cutted_pf = TotSum.Et(), Total_SumEt_15cutted_pf = TotSum.Et(), Total_SumEt_30cutted_pf = TotSum.Et(), Total_SumEt_50cutted_pf = TotSum.Et(), Total_SumEt_500cutted_pf = TotSum.Et();

	if(blabla) cout<<"  -> For this PFCand we loop on pfCand.elementsInBlocks: "<<pfCand.particleId()<<" "<<pfCand.eta()<<" "<<pfCand.phi()<<" E: "<<pfCand.energy()<<" "<<TotSum.Et()<<endl;
	reco::PFClusterRef thiscluster;
	for (auto& blockPair : pfCand.elementsInBlocks()) // that's a misnomer - it gives a pair<PFBlock, unsigned>
	{
	  unsigned int pos = blockPair.second;
	  const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	  {
	    if (blockElement.type() != 4) continue;

	    reco::PFClusterRef cluster = blockElement.clusterRef();
	    if (cluster.isAvailable()) {

		// now we're close!
		DetId seedID = cluster->seed();
		thiscluster = cluster;
		if (cluster->layer() == PFLayer::ECAL_BARREL)
		{
		  //SumET
		  TLorentzVector TSum;
		  TSum.SetPtEtaPhiE( cluster->energy()/cosh( cluster->eta() ), cluster->eta(), cluster->phi(), cluster->energy() );
		  SumEt_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_15cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_30cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_50cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_500cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  if(blabla) cout<<"  -> PFLayer::ECAL_BARREL: "<<cluster->eta()<<" "<<cluster->phi()<<" -> "<<TSum.Et()<<" -> "<<SumEt_cutted_pf<<" E "<<cluster->energy()<<endl;
		  for (auto& rhEB : *recHitsEB)
		  {
		    if (rhEB.id() == seedID){
			seedHit = &rhEB;
			EBDetId idEB(rhEB.id());
			const CaloCellGeometry* cell=geometry_->getGeometry(idEB);
			seedHitPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EB_LAYER_ );
		    }
		  }
		}
		else if (cluster->layer() == PFLayer::ECAL_ENDCAP)
		{
		  //SumET
		  TLorentzVector TSum;
		  TSum.SetPtEtaPhiE( cluster->energy()/cosh( cluster->eta() ), cluster->eta(), cluster->phi(), cluster->energy() );
		  SumEt_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_15cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_30cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  SumEt_50cutted_pf = (pfCand.energy()/cosh( pfCand.eta()));
		  if(blabla) cout<<"  -> PFLayer::ECAL_ENDCAP: "<<cluster->eta()<<" "<<cluster->phi()<<" -> "<<TSum.Et()<<" -> "<<SumEt_cutted_pf<<" E "<<cluster->energy()<<endl;
		  for (auto& rhEE : *recHitsEE)
		  {
		    if (rhEE.id() == seedID)
		    {
			seedHit = &rhEE;
			EKDetId idEE(rhEE.id());
			const CaloCellGeometry* cell=geometry_->getGeometry(idEE);
			seedHitPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EE_LAYER_ );
		    }
		  } 
		}
		if(seedHit){ V_seeds.push_back( *seedHit ); V_cluster.push_back( cluster ); V_seedsPos.push_back( seedHitPos ); if(blabla) cout<<"   -> ADDING one Cluster "<<endl; }
	    }
	  }
	}//All Blocks
	if(blabla) cout<<" -> End of All Blocks"<<endl;
	//Now Select the bigger cluster into the SC
	float BestTime = -1, BestEne = -1, Emin=0;
	TLorentzVector ClustTL;
	for( int nClu=0; nClu<int(V_cluster.size()); nClu++ ){
	  if( V_cluster[nClu]->energy() > Emin ){
	    Emin = V_cluster[nClu]->energy();
	    BestTime = V_seeds[nClu].time();
	    if( SubT0TOF_== 1 ) BestTime-=T0_Vtx_MC;
	    GlobalPoint thisPos( V_seedsPos[nClu].x(), V_seedsPos[nClu].y(), V_seedsPos[nClu].z() );
	    GlobalPoint Zero( 0.,0.,0. );
	    if( SubT0TOF_== 2 ){ BestTime-=computeTOF( Vtx_sim, thisPos ); }
	    if( SubT0TOF_== 3 ){ BestTime-=T0_Vtx_MC; BestTime-=computeTOF( Vtx_sim, thisPos ); }
	    if( SubT0TOF_== 4 ){ BestTime+=computeTOF( Zero, thisPos ); BestTime-=T0_Vtx_MC; BestTime-=computeTOF( Vtx_sim, thisPos ); }
	    BestEne = V_seeds[nClu].energy();
	    ClustTL.SetPtEtaPhiE( V_cluster[nClu]->energy()/cosh(V_cluster[nClu]->eta()) , V_cluster[nClu]->eta(), V_cluster[nClu]->phi(), V_cluster[nClu]->energy() );
	  }
	}
	//If there is a Seed
	h_BestTime_Fir_RemovalSumEt->Fill(BestTime);
	h_BestTime_Fir_RemovalSumEt_zoom->Fill(BestTime);
	if( BestTime!=-1 ){
	  if(blabla) cout<<" -> Trying to remove Time"<<endl;
	  if( (BestTime<-0.01 || BestTime>0.12) ){
	    if(blabla) cout<<" -> DONE!!"<<endl;
	    SumEt_cutted_pf = 0.;
	    Total_SumEt_cutted_pf = 0.;
	  }
	  //Association with a time selection
	  else{
	    if( MC_pres ){
		TLorentzVector PFCan; PFCan.SetPtEtaPhiE( pfCand.energy()/cosh( pfCand.eta()), pfCand.eta(), pfCand.phi(), pfCand.energy() );
		bool isAsso1 = ClustTL.DeltaR( Gamma1 )<0.1 ? true : false; //### ClustTL, best / PFCan
		bool isAsso2 = ClustTL.DeltaR( Gamma2 )<0.1 ? true : false;
		if( isAsso1 ) Associated_v1.push_back( PFCan );
		if( isAsso2 ) Associated_v2.push_back( PFCan );
	    }
	  }
	  //Association w/o a time selection
	  if( MC_pres ){
	    TLorentzVector PFCan; PFCan.SetPtEtaPhiE( pfCand.energy()/cosh( pfCand.eta()), pfCand.eta(), pfCand.phi(), pfCand.energy() );
	    float DR1 = ClustTL.DeltaR( Gamma1 ), DR2 = ClustTL.DeltaR( Gamma2 );
	    bool isAsso1 = DR1<0.5 ? true : false;
	    bool isAsso2 = DR2<0.5 ? true : false;
	    if( isAsso1 ){
		Associated_v1_noTime.push_back( PFCan ); AssociatedTime_v1_noTime.push_back( BestTime ); AssociatedDR_v1_noTime.push_back( DR1 );
		if( blabla2 && MC_pres && DR1<0.1 ) cout<<" ASSO 1: dr: "<<DR1<<" Pt: "<<PFCan.Pt()<<" eta: "<<PFCan.Eta()<<" phi "<<PFCan.Phi()<<endl;
	    }
	    if( isAsso2 ){
		Associated_v2_noTime.push_back( PFCan ); AssociatedTime_v2_noTime.push_back( BestTime ); AssociatedDR_v2_noTime.push_back( DR2 );
		if( blabla2 && MC_pres && DR2<0.1 ) cout<<" ASSO 2: dr: "<<DR1<<" Pt: "<<PFCan.Pt()<<" eta: "<<PFCan.Eta()<<" phi "<<PFCan.Phi()<<endl;
	    }
	    if( DR1>0.1 && DR2>0.1 ){ h_TimeGammaNOTAsso->Fill( BestTime ); h_TimeGammaNOTAssoW->Fill( BestTime, pfCand.energy() ); }
	  }
	  float TimeSmear_15 = SmearTime(BestTime, 0.015, ranGaus_);
	  float TimeSmear_30 = SmearTime(BestTime, 0.03, ranGaus_);
	  float TimeSmear_50 = SmearTime(BestTime, 0.05, ranGaus_);
	  float TimeSmear_500 = SmearTime(BestTime, 0.5, ranGaus_);
	  if( (TimeSmear_15<-0.05 || TimeSmear_15>0.13) ){
	    SumEt_15cutted_pf = 0.;
	    Total_SumEt_15cutted_pf =0.;
	  }
	  if( (TimeSmear_30<-0.07 || TimeSmear_30>0.15) ){
	    SumEt_30cutted_pf = 0.;
	    Total_SumEt_30cutted_pf =0.;
	  }
	  if( (TimeSmear_50<-0.09 || TimeSmear_50>0.16) ){
	    SumEt_50cutted_pf = 0.;
	    Total_SumEt_50cutted_pf =0.;
	  }
	  if( (TimeSmear_500<-1. || TimeSmear_500>1.) ){
	    SumEt_500cutted_pf = 0.;
	    Total_SumEt_500cutted_pf =0.;
	  }
	  h_Time->Fill( BestTime );
	  h_Time_we->Fill( BestTime, BestEne );
	  h_TimeSmeared->Fill( TimeSmear_30 );
	  h_TimeSmeared_we->Fill( TimeSmear_30, BestEne );
	}
	if(blabla) cout<<" -> Now add it to SumEt "<<SumEt_cutted_tot<<" "<<SumEt_cutted_pf<<endl;
	SumEt_tot += SumEt_pf;
	SumEt_cutted_tot += SumEt_cutted_pf;
	SumEt_15cutted_tot += SumEt_15cutted_pf;
	SumEt_30cutted_tot += SumEt_30cutted_pf;
	SumEt_50cutted_tot += SumEt_50cutted_pf;
	SumEt_500cutted_tot += SumEt_500cutted_pf;
	TOT_SumEt_tot += Total_SumEt_pf;
	TOT_SumEt_cutted_tot += Total_SumEt_cutted_pf;
	TOT_SumEt_15cutted_tot += Total_SumEt_15cutted_pf;
	TOT_SumEt_30cutted_tot += Total_SumEt_30cutted_pf;
	TOT_SumEt_50cutted_tot += Total_SumEt_50cutted_pf;
	TOT_SumEt_500cutted_tot += Total_SumEt_500cutted_pf;
    }//END Loop PFCand

    //Best Associated Gamma1
    float bestReso = 99.;
    int BestIndex = -1;
    for(int i=0; i<int(Associated_v1.size()); i++ ){ 
	float resol = fabs(Associated_v1[i].Pt()-Gamma1.Pt())/Gamma1.Pt();
	if( resol < bestReso ){
	  bestReso = resol;
	  BestIndex = i;
	}
    }
    if( BestIndex!=-1 ){
	Associated++;
	if( fabs(Associated_v1[BestIndex].Eta())<1.48 ) Associated_EB++;
	if( fabs(Associated_v1[BestIndex].Eta())>1.48 ) Associated_EE++;
    }
    //Best Associated Gamma2
    bestReso = 99.; BestIndex = -1;
    for(int i=0; i<int(Associated_v2.size()); i++ ){ 
	float resol = fabs(Associated_v2[i].Pt()-Gamma2.Pt())/Gamma2.Pt();
	if( resol < bestReso ){
	  bestReso = resol;
	  BestIndex = i;
	}
    }
    if( BestIndex!=-1 ){
	Associated++;
	if( fabs(Associated_v2[BestIndex].Eta())<1.48 ) Associated_EB++;
	if( fabs(Associated_v2[BestIndex].Eta())>1.48 ) Associated_EE++;
    }
    //Best Associated Gamma1 No Time 
    bestReso = 99.; BestIndex = -1;
    float DR_max = 0.5;
    for(int i=0; i<int(Associated_v1_noTime.size()); i++ ){
	float resol = fabs(Associated_v1_noTime[i].Pt()-Gamma1.Pt())/Gamma1.Pt();
	if( resol < bestReso && AssociatedDR_v1_noTime[i]<0.1 ){
	  bestReso = resol;
	  BestIndex = i;
	}
	if(  AssociatedDR_v1_noTime[i]<DR_max ){
	  DR_max = AssociatedDR_v1_noTime[i];
	}
    }
    h_DRAsso->Fill( DR_max );
    if( BestIndex!=-1 ){
	Associated_time++;
	if( fabs(Associated_v1_noTime[BestIndex].Eta())<1.48 ){
	  Associated_time_EB++;
	  h_TimeGammaAssoEB->Fill( AssociatedTime_v1_noTime[BestIndex] ); h_TimeGammaAssoEB_L->Fill( AssociatedTime_v1_noTime[BestIndex] ); h_TimeGammaAssoEB_L2->Fill( AssociatedTime_v1_noTime[BestIndex] );
	  h_PtGammaAssoEB->Fill( Associated_v1_noTime[BestIndex].Pt() ); h_EGammaAssoEB->Fill( Associated_v1_noTime[BestIndex].E() );
	  h_TimeGammaAssoEB_sme15->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.015, ranGaus_) );
	  h_TimeGammaAssoEB_sme30->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.03, ranGaus_) );
	  h_TimeGammaAssoEB_sme50->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.05, ranGaus_) );
	  h_TimeGammaAssoEB_sme500->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.5, ranGaus_) );
	}
	if( fabs(Associated_v1_noTime[BestIndex].Eta())>1.48 ){
	  Associated_time_EE++;
	  h_TimeGammaAssoEE->Fill( AssociatedTime_v1_noTime[BestIndex] ); h_TimeGammaAssoEE_L->Fill( AssociatedTime_v1_noTime[BestIndex] ); h_TimeGammaAssoEE_L2->Fill( AssociatedTime_v1_noTime[BestIndex] );
	  h_PtGammaAssoEE->Fill( Associated_v1_noTime[BestIndex].Pt() ); h_EGammaAssoEE->Fill(Associated_v1_noTime[BestIndex].E() );
	  h_TimeGammaAssoEE_sme15->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.015, ranGaus_) );
	  h_TimeGammaAssoEE_sme30->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.03, ranGaus_) );
	  h_TimeGammaAssoEE_sme50->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.05, ranGaus_) );
	  h_TimeGammaAssoEE_sme500->Fill(SmearTime(AssociatedTime_v1_noTime[BestIndex], 0.5, ranGaus_) );
	}
    }
    //Best Associated Gamma2 No Time 
    bestReso = 99.; BestIndex = -1;
    DR_max = 0.5;
    for(int i=0; i<int(Associated_v2_noTime.size()); i++ ){
	float resol = fabs(Associated_v2_noTime[i].Pt()-Gamma2.Pt())/Gamma2.Pt();
	if( resol < bestReso && AssociatedDR_v2_noTime[i]<0.1 ){
	  bestReso = resol;
	  BestIndex = i;
	}
	if(  AssociatedDR_v2_noTime[i]<DR_max ){
	  DR_max = AssociatedDR_v2_noTime[i];
	}
    }
    h_DRAsso->Fill( DR_max );
    if( BestIndex!=-1 ){
	Associated_time++;
	if( fabs(Associated_v2_noTime[BestIndex].Eta())<1.48 ){
	  Associated_time_EB++;
	  h_TimeGammaAssoEB->Fill( AssociatedTime_v2_noTime[BestIndex] ); h_TimeGammaAssoEB_L->Fill( AssociatedTime_v2_noTime[BestIndex] ); h_TimeGammaAssoEB_L2->Fill( AssociatedTime_v2_noTime[BestIndex] );
	  h_PtGammaAssoEB->Fill( Associated_v2_noTime[BestIndex].Pt() ); h_EGammaAssoEB->Fill( Associated_v2_noTime[BestIndex].E() );
	  h_TimeGammaAssoEB_sme15->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.015, ranGaus_) );
	  h_TimeGammaAssoEB_sme30->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.03, ranGaus_) );
	  h_TimeGammaAssoEB_sme50->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.05, ranGaus_) );
	  h_TimeGammaAssoEB_sme500->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.5, ranGaus_) );
	}
	if( fabs(Associated_v2_noTime[BestIndex].Eta())>1.48 ){
	  Associated_time_EE++;
	  h_TimeGammaAssoEE->Fill( AssociatedTime_v2_noTime[BestIndex] ); h_TimeGammaAssoEE_L->Fill( AssociatedTime_v2_noTime[BestIndex] ); h_TimeGammaAssoEE_L2->Fill( AssociatedTime_v2_noTime[BestIndex] );
	  h_PtGammaAssoEE->Fill( Associated_v2_noTime[BestIndex].Pt() ); h_EGammaAssoEE->Fill(Associated_v2_noTime[BestIndex].E() );
	  h_TimeGammaAssoEE_sme15->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.015, ranGaus_) );
	  h_TimeGammaAssoEE_sme30->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.03, ranGaus_) );
	  h_TimeGammaAssoEE_sme50->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.05, ranGaus_) );
	  h_TimeGammaAssoEE_sme500->Fill(SmearTime(AssociatedTime_v2_noTime[BestIndex], 0.5, ranGaus_) );
	}
    }

    h_SumEt->Fill(SumEt_tot);
    h_SumEt_cut->Fill(SumEt_cutted_tot);
    h_SumEt_15cut->Fill(SumEt_15cutted_tot);
    h_SumEt_30cut->Fill(SumEt_30cutted_tot);
    h_SumEt_50cut->Fill(SumEt_50cutted_tot);
    h_SumEt_500cut->Fill(SumEt_500cutted_tot);
    h_TOT_SumEt->Fill(TOT_SumEt_tot);
    h_TOT_SumEt_cut->Fill(TOT_SumEt_cutted_tot);
    h_TOT_SumEt_15cut->Fill(TOT_SumEt_15cutted_tot);
    h_TOT_SumEt_30cut->Fill(TOT_SumEt_30cutted_tot);
    h_TOT_SumEt_50cut->Fill(TOT_SumEt_50cutted_tot);
    h_TOT_SumEt_500cut->Fill(TOT_SumEt_500cutted_tot);
  }// if DoSumEt_
  // Mass and Energy resolution cleaning
  if( DoMass_ && isHgg_ ){
    //Hgg Mass----
    // Associate MC gammas with PFCand == gamma
    TLorentzVector gamma_reco1, gamma_reco2; gamma_reco1.SetPtEtaPhiE( -1., -1., -1., -1. ); gamma_reco2.SetPtEtaPhiE( -1., -1., -1., -1. );
    reco::PFCandidate Pfgamma1_re, Pfgamma2_re;
    const reco::PFJet* VBF1_re=0, *VBF2_re=0;
    float PU_frac_1(1.), PU_frac_2(1.);
    float R91(-1.), R92(-1.);
    std::vector<float> v_R9; v_R9.clear();
    bool foundVBF1=false, foundVBF2=false;
    if( MC_pres ){
	//Assocation VBFJets
	float DR_min = MinDR_asso;
	for (auto& pfJet : *Jets){
	  if( pfJet.p4().Pt() < MinPt_Reco ) continue;
	  GlobalPoint PosJet( pfJet.p4().X(), pfJet.p4().Y(), pfJet.p4().Z() );
	  float DR = DeltaR( PosJet, VBF1 );
	  if( DR<DR_min ){
	    DR_min = DR;
	    VBF1_re = &pfJet;
	    foundVBF1 = true;
	  }
	}
	DR_min = MinDR_asso;
	for (auto& pfJet : *Jets){
	  if( pfJet.p4().Pt() < MinPt_Reco ) continue;
	  GlobalPoint PosJet( pfJet.p4().X(), pfJet.p4().Y(), pfJet.p4().Z() );
	  float DR = DeltaR( PosJet, VBF2 );
	  if( DR<DR_min ){
	    DR_min = DR;
	    VBF2_re = &pfJet;
	    foundVBF2 = true;
	  }
	}
	//Association Gammas
	float miniDr = 0.1; 
	for (auto& pfcan : *PFCol){
	  float ErreNove = FillLateralDevel( pfcan, recHitsEB, recHitsEE, false ); //For Sig + Bkg
	  v_R9.push_back( ErreNove );
	  if( pfcan.pt()<15. && pfcan.particleId() != reco::PFCandidate::gamma ) continue;
	  TLorentzVector PosPhot; PosPhot.SetPtEtaPhiE( pfcan.pt(), pfcan.p4().Eta(), pfcan.p4().Phi(), pfcan.p4().E() );
	  float DR = PosPhot.DeltaR( Gamma1 );
	  if( DR<miniDr ){ miniDr = DR; gamma_reco1.SetPtEtaPhiE( pfcan.pt(), pfcan.p4().Eta(), pfcan.p4().Phi(), pfcan.p4().E() ); Pfgamma1_re = pfcan; PU_frac_1 = Compute_PUfrac( pfcan, recHitsEB, recHitsEE ); R91 = ErreNove;}
	}
	miniDr = 0.1;
	int NumPf=0;
	for (auto& pfcan : *PFCol){
	  if( pfcan.pt()<15.  && pfcan.particleId() != reco::PFCandidate::gamma ) continue;
	  TLorentzVector PosPhot; PosPhot.SetPtEtaPhiE( pfcan.pt(), pfcan.p4().Eta(), pfcan.p4().Phi(), pfcan.p4().E() );
	  float DR = PosPhot.DeltaR( Gamma2 );
	  if( DR<miniDr && pfcan.pt() != gamma_reco1.Pt() ){ miniDr = DR; gamma_reco2.SetPtEtaPhiE( pfcan.pt(), pfcan.p4().Eta(), pfcan.p4().Phi(), pfcan.p4().E() ); Pfgamma2_re = pfcan; PU_frac_2 = Compute_PUfrac( pfcan, recHitsEB, recHitsEE ); R92 = v_R9[NumPf];}
	  NumPf++;
	}
    }
    h_R91->Fill( R91 );
    h_R92->Fill( R92 );
    float Min_R9 = Use_R9_ ? 0.6 : -19191.;
    bool isR9( R91>Min_R9 && R92>Min_R9 );
    bool isRecon( gamma_reco1.Pt()>15. && gamma_reco2.Pt()>15 && fabs(gamma_reco1.Pt()-Gamma1.Pt())/Gamma1.Pt()<0.35 && fabs(gamma_reco2.Pt()-Gamma2.Pt())/Gamma2.Pt()<0.35 && isR9 );
    if( isRecon ){
	//Lateral Time Development
	if( fabs(gamma_reco1.Pt()-Gamma1.Pt())/Gamma1.Pt()<0.1 && fabs(gamma_reco2.Pt()-Gamma2.Pt())/Gamma2.Pt()<0.1 ){
	  R91 = FillLateralDevel( Pfgamma1_re, recHitsEB, recHitsEE, true ); //For Sig
	  R92 = FillLateralDevel( Pfgamma2_re, recHitsEB, recHitsEE, true ); //For Sig
	}
	if( WannaFitT0Vtx_ ){
	  //Std method
	  vector<float> Seed_info1 = fabs(Pfgamma1_re.p4().eta())<1.5 ? GetSeedFromSC( true, Pfgamma1_re, recHitsEB, recHitsEE, false, true ) : GetSeedFromSC( false, Pfgamma1_re, recHitsEB, recHitsEE, false, true ); //Seed_PtEtaPhiTime_CluPeEtaPhiE (def. -999.)
	  vector<float> Seed_info2 = fabs(Pfgamma2_re.p4().eta())<1.5 ? GetSeedFromSC( true, Pfgamma2_re, recHitsEB, recHitsEE, false, true ) : GetSeedFromSC( false, Pfgamma2_re, recHitsEB, recHitsEE, false, true ); //Seed_PtEtaPhiTime_CluPeEtaPhiE (def. -999.)
	  TLorentzVector Seed1_tl; Seed1_tl.SetPtEtaPhiM( Seed_info1[0], Seed_info1[1], Seed_info1[2], 0. );
	  TLorentzVector Seed2_tl; Seed2_tl.SetPtEtaPhiM( Seed_info2[0], Seed_info2[1], Seed_info2[2], 0. );
	  if( Seed_info1[3]!=-999 && Seed_info2[3]!=-999 ){
	    vector<float> null; null.push_back(-999.); null.push_back(-999.); null.push_back(-999.);
	    vector<float> VBF1_info = foundVBF1 ? GetTimeFromJet( VBF1_re, recHitsEB, recHitsEE, geometry_, FTool_ ) : null;
	    vector<float> VBF2_info = foundVBF2 ? GetTimeFromJet( VBF2_re, recHitsEB, recHitsEE, geometry_, FTool_ ) : null;
	    GlobalPoint PosTot(  Seed1_tl.X()-Vtx_sim.x(),  Seed1_tl.Y()-Vtx_sim.y() ,  Seed1_tl.Z()-Vtx_sim.z() );
	    float tof1    = T0_Vtx_MC +  sqrt( pow(PosTot.x(),2) + pow(PosTot.y(),2) + pow(PosTot.z(),2) )/(LIGHT_SPEED);
	    //float tofVBF1 = T0_Vtx_MC +  sqrt( pow(foundVBF1 ? VBF1_re->p4().X() : -999.,2) + pow(foundVBF1 ? VBF1_re->p4().Y() : -999.,2) + pow(foundVBF1 ? VBF1_re->p4().Z() : -999.,2) )/(LIGHT_SPEED);
	    float tofVBF1 = foundVBF1 ? (T0_Vtx_MC +  sqrt( pow( VBF1_info[1],2) + pow( VBF1_info[2],2) + pow( VBF1_info[3],2) )/(LIGHT_SPEED)) : -999.; 
	    GlobalPoint PosTot2(  Seed2_tl.X()-Vtx_sim.x(),  Seed2_tl.Y()-Vtx_sim.y() ,  Seed2_tl.Z()-Vtx_sim.z() );
	    float tof2    = T0_Vtx_MC +  sqrt( pow(PosTot2.x(),2) + pow(PosTot2.y(),2) + pow(PosTot2.z(),2) )/(LIGHT_SPEED);
	    //float tofVBF2 = T0_Vtx_MC +  sqrt( pow(foundVBF2 ? VBF2_re->p4().X() : -999.,2) + pow(foundVBF2 ? VBF2_re->p4().Y() : -999.,2) + pow(foundVBF2 ? VBF2_re->p4().Z() : -999.,2) )/(LIGHT_SPEED);
	    float tofVBF2 = foundVBF2 ? (T0_Vtx_MC +  sqrt( pow( VBF2_info[1],2) + pow( VBF2_info[2],2) + pow( VBF2_info[3],2) )/(LIGHT_SPEED)) : -999.;
	    VtxDet_VBFT1       = foundVBF1 ? VBF1_info[0] : -999.; 
	    VtxDet_VBFT2       = foundVBF2 ? VBF2_info[0] : -999.;
	    VtxDet_T1          = Seed_info1[3];
	    VtxDet_T2          = Seed_info2[3];
	    VtxDet_GT1         = tof1;
	    VtxDet_GT2         = tof2;
	    VtxDet_VBFGT1      = foundVBF1 ? tofVBF1 : -999.;
	    VtxDet_VBFGT2      = foundVBF1 ? tofVBF2 : -999.;
	    VtxDet_time        = T0_Vtx_MC;
	    VtxDet_PosXtal_X1  = Seed1_tl.X();
	    VtxDet_PosXtal_Y1  = Seed1_tl.Y();
	    VtxDet_PosXtal_Z1  = Seed1_tl.Z();
	    VtxDet_PosXtal_X2  = Seed2_tl.X();
	    VtxDet_PosXtal_Y2  = Seed2_tl.Y();
	    VtxDet_PosXtal_Z2  = Seed2_tl.Z();
	    VtxDet_PosVBF_X1   = foundVBF1 ? VBF1_info[1] : -999.;
	    VtxDet_PosVBF_Y1   = foundVBF1 ? VBF1_info[2] : -999.;
	    VtxDet_PosVBF_Z1   = foundVBF1 ? VBF1_info[3] : -999.;
	    VtxDet_PosVBF_X2   = foundVBF2 ? VBF2_info[1] : -999.;
	    VtxDet_PosVBF_Y2   = foundVBF2 ? VBF2_info[2] : -999.;
	    VtxDet_PosVBF_Z2   = foundVBF2 ? VBF2_info[3] : -999.;
	    VtxDet_PosXtal_MCX1= Gamma1.X();
	    VtxDet_PosXtal_MCY1= Gamma1.Y();
	    VtxDet_PosXtal_MCZ1= Gamma1.Z();
	    VtxDet_PosXtal_MCX2= Gamma2.X();
	    VtxDet_PosXtal_MCY2= Gamma2.Y();
	    VtxDet_PosXtal_MCZ2= Gamma2.Z();
	    VtxDet_PtRecoJet_1 = Pfgamma1_re.p4().Pt();
	    VtxDet_PtRecoJet_2 = Pfgamma2_re.p4().Pt();;
	    VtxDet_PtMCJet_1   = Gamma1.Pt();
	    VtxDet_PtMCJet_2   = Gamma2.Pt();
	    vzMC               = Vtx_sim.z();
	    vyMC               = Vtx_sim.y();
	    vxMC               = Vtx_sim.x();
	    vzRECO             = Vtx_reco.z();
	    vyRECO             = Vtx_reco.y();
	    vxRECO             = Vtx_reco.x();
	    Tree_Vtx->Fill();
	  }
	}//WannaFitT0Vtx
	//Correction for right vertex and PFMAss
	TVector3 gamma_reco1_poscorr; gamma_reco1_poscorr.SetXYZ( gamma_reco1.X() + Vtx_reco.x(), gamma_reco1.Y() + Vtx_reco.y(), gamma_reco1.Z() + Vtx_reco.z() );
	TVector3 gamma_reco2_poscorr; gamma_reco2_poscorr.SetXYZ( gamma_reco2.X() + Vtx_reco.x(), gamma_reco2.Y() + Vtx_reco.y(), gamma_reco2.Z() + Vtx_reco.z() );
	gamma_reco1_poscorr.SetXYZ( gamma_reco1_poscorr.X() - Vtx_sim.x(), gamma_reco1_poscorr.Y() - Vtx_sim.y(), gamma_reco1_poscorr.Z() - Vtx_sim.z() );
	gamma_reco2_poscorr.SetXYZ( gamma_reco2_poscorr.X() - Vtx_sim.x(), gamma_reco2_poscorr.Y() - Vtx_sim.y(), gamma_reco2_poscorr.Z() - Vtx_sim.z() );
	//gamma_reco1_poscorr.SetXYZ( gamma_reco1.X() - Vtx_sim.x(), gamma_reco1.Y() - Vtx_sim.y(), gamma_reco1.Z() - Vtx_sim.z() );
	//gamma_reco2_poscorr.SetXYZ( gamma_reco2.X() - Vtx_sim.x(), gamma_reco2.Y() - Vtx_sim.y(), gamma_reco2.Z() - Vtx_sim.z() );
	////TLorentzVector gamma_reco1_corr; gamma_reco1_corr.SetPtEtaPhiE( gamma_reco1.E()/cosh(gamma_reco1_poscorr.Eta()), gamma_reco1_poscorr.Eta(), gamma_reco1_poscorr.Phi(), gamma_reco1.E() );
	////TLorentzVector gamma_reco2_corr; gamma_reco2_corr.SetPtEtaPhiE( gamma_reco2.E()/cosh(gamma_reco2_poscorr.Eta()), gamma_reco2_poscorr.Eta(), gamma_reco2_poscorr.Phi(), gamma_reco2.E() );
	TLorentzVector gamma_reco1_corr; gamma_reco1_corr.SetPtEtaPhiM( gamma_reco1.E()/cosh(gamma_reco1_poscorr.Eta()), gamma_reco1_poscorr.Eta(), gamma_reco1_poscorr.Phi(), 0.);
	TLorentzVector gamma_reco2_corr; gamma_reco2_corr.SetPtEtaPhiM( gamma_reco2.E()/cosh(gamma_reco2_poscorr.Eta()), gamma_reco2_poscorr.Eta(), gamma_reco2_poscorr.Phi(), 0.);
	h_HiggsPFMass_Vtx->Fill( (gamma_reco1_corr+gamma_reco2_corr).M() );
	h_HiggsPFMass->Fill( (gamma_reco1+gamma_reco2).M() );
	Ereso1->Fill( gamma_reco1.E()/Gamma1.E() );
	Ereso2->Fill( gamma_reco2.E()/Gamma2.E() );
	//Now Mitigation for the PU faction
	gamma_reco1.SetPtEtaPhiM( gamma_reco1.Pt()*PU_frac_1, gamma_reco1.Eta(), gamma_reco1.Phi(), 0. );
	gamma_reco2.SetPtEtaPhiM( gamma_reco2.Pt()*PU_frac_2, gamma_reco2.Eta(), gamma_reco2.Phi(), 0. );
	h_PU_frac_1->Fill( PU_frac_1 );
	h_PU_frac_2->Fill( PU_frac_2 );
	h_HiggsPFMass_Mit->Fill( (gamma_reco1+gamma_reco2).M() );
	Ereso1_Mit->Fill( gamma_reco1.E()/Gamma1.E() );
	Ereso2_Mit->Fill( gamma_reco2.E()/Gamma2.E() );
    }//isRecon
    // Associate with Gamma
    gamma_reco1.SetPtEtaPhiE( -1., -1., -1., -1. ); gamma_reco2.SetPtEtaPhiE( -1., -1., -1., -1. );
    if( Gamma1.Pt()>15. && Gamma2.Pt()>15. ){
	float miniDr = 0.1; 
	for(auto& photon : *Photons){
	  if( photon.pt()<15. ) continue;
	  TLorentzVector PosPhot; PosPhot.SetPtEtaPhiE( photon.pt(), photon.p4().Eta(), photon.p4().Phi(), photon.p4().E() );
	  float DR = PosPhot.DeltaR( Gamma1 );
	  if( DR<miniDr ){ miniDr = DR; gamma_reco1.SetPtEtaPhiE( photon.pt(), photon.p4().Eta(), photon.p4().Phi(), photon.p4().E() ); }
	}
	miniDr = 0.1;
	for(auto& photon : *Photons){
	  if( photon.pt()<15. ) continue;
	  TLorentzVector PosPhot; PosPhot.SetPtEtaPhiE( photon.pt(), photon.p4().Eta(), photon.p4().Phi(), photon.p4().E() );
	  float DR = PosPhot.DeltaR( Gamma2 );
	  if( DR<miniDr && photon.pt() != gamma_reco1.Pt() ){ miniDr = DR; gamma_reco2.SetPtEtaPhiE( photon.pt(), photon.p4().Eta(), photon.p4().Phi(), photon.p4().E() ); }
	}
    }
    isRecon = gamma_reco1.Pt()>15. && gamma_reco2.Pt()>15 && fabs(gamma_reco1.Pt()-Gamma1.Pt())/Gamma1.Pt()<0.35 && fabs(gamma_reco2.Pt()-Gamma2.Pt())/Gamma2.Pt()<0.35;
    if( isRecon ) h_HiggsMass->Fill( (gamma_reco1+gamma_reco2).M() );
  }//if DoMass_

}//End Analyzer

//------------------------------------------------------------------------------------------------
float Generic_Analizer::Compute_PUfrac( reco::PFCandidate pfcan, edm::Handle<edm::SortedCollection<EcalRecHit> >& theBarrelEcalRecHits, edm::Handle<edm::SortedCollection<EcalRecHit> >& theEndcapEcalRecHits )
{
  float Min_ene = Use_MinEne_ ? 0.2 : -19191.;
  reco::PFCandidate MY_cand( pfcan );
  std::vector<EcalRecHit> V_seeds; V_seeds.clear(); std::vector<reco::PFClusterRef> V_cluster; V_cluster.clear(); std::vector<GlobalPoint> V_PFPos; V_PFPos.clear();
  const EcalRecHit* seedHit = 0;
  for (auto& blockPair : MY_cand.elementsInBlocks()){
    unsigned int pos = blockPair.second;
    const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
    {
	if (blockElement.type() != 4) continue;
	reco::PFClusterRef cluster = blockElement.clusterRef();
	if (cluster.isAvailable())
	{
	  GlobalPoint PFPos(-99., -99., -99);
	  DetId seedID = cluster->seed();
	  if (cluster->layer() == PFLayer::ECAL_BARREL){
	    for (auto& rhEB : *theBarrelEcalRecHits){
		if (rhEB.id() == seedID){
		  seedHit = &rhEB;
		  EBDetId IdXtal( seedID );
		  const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
		  PFPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EB_LAYER_ );
		}
	    }
	  }
	  else if (cluster->layer() == PFLayer::ECAL_ENDCAP){
	    for (auto& rhEE : *theEndcapEcalRecHits){
		if (rhEE.id() == seedID){
		  seedHit = &rhEE;
		  EKDetId IdXtal( seedID );
		  const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
		  PFPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EE_LAYER_ );
		}
	    }
	  }
	  if(seedHit){ V_seeds.push_back( *seedHit ); V_cluster.push_back( cluster ); V_PFPos.push_back( PFPos ); }
	}
    }
  }//All Blocks
  float Cut_MIN[6]={0., 0., 0., 0., 0., 0.};
  //float Cut_MAX[6]={0.11, 0.39, 0.61, 0.74, 0.78, 1.};
  //float Cut_MAX[6]={0.2, 1.4, 20., 50., 66., 1000.}; //Quant 1-1 with nois, still good performenaces
  float Cut_MAX[6]={0.18, 1.5, 8., 14., 12., 1000.}; //Quant 1-1 w/o noise, GREAT performances
  //float Cut_MAX[6]={0.15, 1.8, 7., 9., 12., 1000.}; //Quant 1-1 w/o noise, R9 cut
  //float Cut_MAX[6]={0.21, 4., 28., 50., 70., 1000.}; //Super loose with nois, DR GOOD
  //float Cut_MAX[6]={0.21, 2., 9., 9., 6., 1000.}; //Super loose w/o nois, DR GOOD
  //float Cut_MAX[6]={0.25, 1.5, 10., 1000., 1000., 1000.}; //Quant 1-1 w/o noise, R9 cut, DR precise
  //For each Cluster perform the cleaning:
  float TotalE=0., AfterCut=0.;
  std::vector<uint32_t> Used_id; Used_id.clear();
  for( int nClu=0; nClu<int(V_cluster.size()); nClu++ ){
    TotalE += V_cluster[nClu]->energy();
    float SeedTime = V_seeds[nClu].time();
    //First: look if the seed is ok
    if( SeedTime<Cut_MIN[0] || SeedTime>Cut_MAX[0] ) continue;
    //Second: if so, clean the recHits
    else{
	//AfterCut += V_cluster[nClu]->energy();
	std::vector< std::pair<DetId, float> > RecHitsFrac =  V_cluster[nClu]->hitsAndFractions();
	for( int i=0; i<int(RecHitsFrac.size()); i++ ){
	  bool found = true;
	  DetId My_id( RecHitsFrac[i].first );
	  if(std::find(Used_id.begin(), Used_id.end(), My_id.rawId()) != Used_id.end()) continue;
	  Used_id.push_back( My_id.rawId() );
	  EcalRecHitCollection::const_iterator My_rec = theBarrelEcalRecHits->find( My_id );
	  if( My_rec==theBarrelEcalRecHits->end() ){
	    My_rec = theEndcapEcalRecHits->find( My_id );
	    if( My_rec==theEndcapEcalRecHits->end() ){ found = false; }
	  }
	  if( !found ) continue;
	  else{
	    GlobalPoint RecPos;
	    if( fabs(pfcan.p4().Eta())<1.476 ){
		EBDetId IdXtal( My_id );
		const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
		RecPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EB_LAYER_ );
	    }
	    else{
		EKDetId IdXtal( My_id );
		const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
		RecPos = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EE_LAYER_ );
	    }
	    float time = My_rec->time();
	    float deltaEr = DeltaR( V_PFPos[nClu].eta(), RecPos.eta(), V_PFPos[nClu].phi(), RecPos.phi() );
	    //float deltaEr = DeltaR( V_cluster[nClu]->eta(), RecPos.eta(), V_cluster[nClu]->phi(), RecPos.phi() );
	    float DR_step = 0.017;
	    if( My_rec->energy() > Min_ene ){ if( time>Cut_MIN[int(deltaEr/DR_step)] && time<Cut_MAX[int(deltaEr/DR_step)] ) AfterCut += My_rec->energy(); }
	    else AfterCut += My_rec->energy();
	  }
	}
    }//else
  }//all Clusters
  float PU_frac = (TotalE==0.) ? 0. : float(AfterCut/TotalE);
  return PU_frac;
}
// ------------------------------------------------------------------------------------------
float Generic_Analizer::FillLateralDevel( reco::PFCandidate Gamma, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEB, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEE, bool isSig ){

  float erNine(-1.);
  float Min_ene = Use_MinEne_ ? 0.2 : -19191.;
  float Min_R9 = Use_R9_ ? 0.6 : -19191.;
  if( isSig ){
    bool isEB( fabs(Gamma.p4().Eta())<1.5 );  // bool isEB( fabs(Gamma.Eta())<1.5 );
    if( isEB ){
	//Find the Seed Info (Seed of more Ene SC)
	vector<float> Seed_info = GetSeedFromSC( isEB, Gamma, recHitsEB, recHitsEE, true, false ); //Seed_PtEtaPhiTime_Clu_PtEtaPhiE (def. -999.)
	erNine = Seed_info[8];
	TLorentzVector Gamma_posit; Gamma_posit.SetPtEtaPhiE( Seed_info[4], Seed_info[5], Seed_info[6], Seed_info[7]);
	//If found open a cone
	if( Seed_info[0]!=-999. ){
	  TLorentzVector Seed_pos; Seed_pos.SetPtEtaPhiM( Seed_info[0], Seed_info[1], Seed_info[2], 0. );
	  float seedtime = Seed_info[3];
	  if( Seed_info[0]>Min_ene && erNine>Min_R9 ){
	    h_DR_vs_Time_EB->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EB2->Fill( 0.0000001, seedtime ); h_DR_vs_Time_L_EB->Fill( 0.0000001, seedtime );
	    h_DR_vs_Time_EB_reb->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EB_reb2->Fill( 0.0000001, seedtime );
	  }
	  for( auto& rech : *recHitsEB ){
	    if( rech.energy()<=0 ) continue;
	    EBDetId IdXtal( rech.id() );
	    const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	    GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    TLorentzVector PosRechit_tl; PosRechit_tl.SetPtEtaPhiM( rech.energy()/cosh(PosRechit.eta() ), PosRechit.eta(), PosRechit.phi(), 0. );
	    float DR = PosRechit_tl.DeltaR( Seed_pos );
	    //float DR = PosRechit_tl.DeltaR( Gamma_posit );
	    if( DR<0.0855 && rech.energy()>Min_ene && erNine>Min_R9 ){
		h_DR_vs_Time_EB->Fill( DR, rech.time() ); h_DR_vs_Time_EB2->Fill( DR, rech.time() ); h_DR_vs_Time_EB_reb->Fill( DR, rech.time() ); h_DR_vs_Time_EB_reb2->Fill( DR, rech.time() );
		h_DR_vs_Time_L_EB->Fill( DR, rech.time() );
		h_energyForDR_EB->Fill( rech.energy() );
		h_timeForDR_EB->Fill( rech.time() );
		h_timevsEne_EB->Fill( rech.time(), rech.energy() );
		h_timevsEne_EB2->Fill( rech.time(), rech.energy() );
	    }
	  }// end for
	}
    }//isEB
    else{
	//Find the Seed Info (Seed of more Ene SC)
	vector<float> Seed_info = GetSeedFromSC( isEB, Gamma, recHitsEB, recHitsEE, true, false ); //Seed_PtEtaPhiTime_CluPtEtaPhiE (def. -999.)
	erNine = Seed_info[8];
	TLorentzVector Gamma_posit; Gamma_posit.SetPtEtaPhiE( Seed_info[4], Seed_info[5], Seed_info[6], Seed_info[7]);
	//If found open a cone
	if( Seed_info[0]!=-999. ){
	  TLorentzVector Seed_pos; Seed_pos.SetPtEtaPhiM( Seed_info[0], Seed_info[1], Seed_info[2], 0. );
	  float seedtime = Seed_info[3];
	  if( Seed_info[0]>Min_ene && erNine>Min_R9 ){
	    h_DR_vs_Time_EE->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EE2->Fill( 0.0000001, seedtime ); h_DR_vs_Time_L_EE->Fill( 0.0000001, seedtime );
	    h_DR_vs_Time_EE_reb->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EE_reb2->Fill( 0.0000001, seedtime );
	  }
	  for( auto& rech : *recHitsEE ){
	    EKDetId IdXtal( rech.id() );
	    if( rech.energy()<=0 ) continue;
	    const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	    GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    TLorentzVector PosRechit_tl; PosRechit_tl.SetPtEtaPhiM( rech.energy()/cosh(PosRechit.eta() ), PosRechit.eta(), PosRechit.phi(), 0. );
	    float DR = PosRechit_tl.DeltaR( Seed_pos );
	    //float DR = PosRechit_tl.DeltaR( Gamma_posit );
	    if( DR<0.0855 && rech.energy()>Min_ene && erNine>Min_R9 ){
		h_DR_vs_Time_EE->Fill( DR, rech.time() ); h_DR_vs_Time_EE2->Fill( DR, rech.time() ); h_DR_vs_Time_EE_reb->Fill( DR, rech.time() ); h_DR_vs_Time_EE_reb2->Fill( DR, rech.time() ); 
		h_DR_vs_Time_L_EE->Fill( DR, rech.time() );
		h_energyForDR_EE->Fill( rech.energy() );
		h_timeForDR_EE->Fill( rech.time() );
		h_timevsEne_EE->Fill( rech.time(), rech.energy() );
		h_timevsEne_EE2->Fill( rech.time(), rech.energy() );
	    }
	  }// end for
	}
    }//isEE
  }
  else{
    bool isEB( fabs(Gamma.p4().Eta())<1.5 );  // bool isEB( fabs(Gamma.Eta())<1.5 );
    if( isEB ){
	//Find the Seed Info (Seed of more Ene SC)
	vector<float> Seed_info = GetSeedFromSC( isEB, Gamma, recHitsEB, recHitsEE, true, false ); //Seed_PtEtaPhiTime_CluPtEtaPhiE (def. -999.)
	erNine = Seed_info[8];
	TLorentzVector Gamma_posit; Gamma_posit.SetPtEtaPhiE( Seed_info[4], Seed_info[5], Seed_info[6], Seed_info[7]);
	//If found open a cone
	if( Seed_info[0]!=-999. ){
	  TLorentzVector Seed_pos; Seed_pos.SetPtEtaPhiM( Seed_info[0], Seed_info[1], Seed_info[2], 0. );
	  float seedtime = Seed_info[3];
	  if( Seed_info[0]>Min_ene && erNine>Min_R9 ){
	    h_DR_vs_Time_EB_b->Fill( 0.0000001, seedtime ); h_DR_vs_Time_L_EB_b->Fill( 0.0000001, seedtime );
	    h_DR_vs_Time_EB_reb_b->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EB_reb2_b->Fill( 0.0000001, seedtime );
	  }
	  for( auto& rech : *recHitsEB ){
	    if( rech.energy()<=0 ) continue;
	    EBDetId IdXtal( rech.id() );
	    const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	    GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    TLorentzVector PosRechit_tl; PosRechit_tl.SetPtEtaPhiM( rech.energy()/cosh(PosRechit.eta() ), PosRechit.eta(), PosRechit.phi(), 0. );
	    float DR = PosRechit_tl.DeltaR( Seed_pos );
	    //float DR = PosRechit_tl.DeltaR( Gamma_posit );
	    if( DR<0.0855 && rech.energy()>Min_ene && erNine>Min_R9 ){
		h_DR_vs_Time_EB_b->Fill( DR, rech.time() ); h_DR_vs_Time_EB_reb_b->Fill( DR, rech.time() ); h_DR_vs_Time_EB_reb2_b->Fill( DR, rech.time() );
		h_DR_vs_Time_L_EB_b->Fill( DR, rech.time() );
		h_energyForDR_EB_b->Fill( rech.energy() );
		h_timeForDR_EB_b->Fill( rech.time() );
		h_timevsEne_EB_b->Fill( rech.time(), rech.energy() );
		h_timevsEne_EB2_b->Fill( rech.time(), rech.energy() );
	    }
	  }// end for
	}
    }//isEB
    else{
	//Find the Seed Info (Seed of more Ene SC)
	vector<float> Seed_info = GetSeedFromSC( isEB, Gamma, recHitsEB, recHitsEE, true, false ); //Seed_PtEtaPhiTime_CluPtEtaPhiE (def. -999.)
	erNine = Seed_info[8];
	TLorentzVector Gamma_posit; Gamma_posit.SetPtEtaPhiE( Seed_info[4], Seed_info[5], Seed_info[6], Seed_info[7]);
	//If found open a cone

	if( Seed_info[0]!=-999. ){
	  TLorentzVector Seed_pos; Seed_pos.SetPtEtaPhiM( Seed_info[0], Seed_info[1], Seed_info[2], 0. );
	  float seedtime = Seed_info[3];
	  if( Seed_info[0]>Min_ene && erNine>Min_R9 ){
	    h_DR_vs_Time_EE_b->Fill( 0.0000001, seedtime ); h_DR_vs_Time_L_EE_b->Fill( 0.0000001, seedtime );
	    h_DR_vs_Time_EE_reb_b->Fill( 0.0000001, seedtime ); h_DR_vs_Time_EE_reb2_b->Fill( 0.0000001, seedtime );
	  }
	  for( auto& rech : *recHitsEE ){
	    if( rech.energy()<=0 ) continue;
	    EKDetId IdXtal( rech.id() );
	    const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	    GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    TLorentzVector PosRechit_tl; PosRechit_tl.SetPtEtaPhiM( rech.energy()/cosh(PosRechit.eta() ), PosRechit.eta(), PosRechit.phi(), 0. );
	    float DR = PosRechit_tl.DeltaR( Seed_pos );
	    //float DR = PosRechit_tl.DeltaR( Gamma_posit );
	    if( DR<0.0855 && rech.energy()>Min_ene && erNine>Min_R9 ){
		h_DR_vs_Time_EE_b->Fill( DR, rech.time() ); h_DR_vs_Time_EE_reb_b->Fill( DR, rech.time() ); h_DR_vs_Time_EE_reb2_b->Fill( DR, rech.time() ); 
		h_DR_vs_Time_L_EE_b->Fill( DR, rech.time() );
		h_energyForDR_EE_b->Fill( rech.energy() );
		h_timeForDR_EE_b->Fill( rech.time() );
		h_timevsEne_EE_b->Fill( rech.time(), rech.energy() );
		h_timevsEne_EE2_b->Fill( rech.time(), rech.energy() );
	    }
	  }// end for
	}
    }//isEE
  }
  return erNine;
}

// ------------------------------------------------------------------------------------------
std::vector<float> Generic_Analizer::GetSeedFromSC( bool isEB, reco::PFCandidate Gamma, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEB, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEE, bool atZero, bool HFill ) //SeedPtEtaPhiTime_CluPeEtaPhiE
{
  float erNine(-1.);
  vector<float> SeedInfo; SeedInfo.clear();
  float S_Pt=-999., S_Eta=-999., S_Phi=-999., S_Time=-999., C_Pt=-999., C_Eta=-999., C_Phi=-999., C_E=-999.;
  std::vector<EcalRecHit> V_seeds; V_seeds.clear(); std::vector<reco::PFClusterRef> V_cluster; V_cluster.clear();
  PFCandidate MY_cand( Gamma );
  const EcalRecHit* seedHit = 0;
  size_t nEcalClusters = 0;
  for (auto& blockPair : MY_cand.elementsInBlocks()){
    unsigned int pos = blockPair.second;
    const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
    {
	if (blockElement.type() != 4) continue;
	reco::PFClusterRef cluster = blockElement.clusterRef();
	if ( cluster.isAvailable() ) {
	  ++nEcalClusters;
	  DetId seedID = cluster->seed();
	  if (cluster->layer() == PFLayer::ECAL_BARREL){
	    for (auto& rhEB : *recHitsEB){
		if (rhEB.id() == seedID){
		  seedHit = &rhEB;
		}
	    }
	  }
	  else if (cluster->layer() == PFLayer::ECAL_ENDCAP){
	    for (auto& rhEE : *recHitsEE){
		if (rhEE.id() == seedID){
		  seedHit = &rhEE;
		}
	    }
	  }
	  if(seedHit){ V_seeds.push_back( *seedHit ); V_cluster.push_back( cluster ); }
	}
    }
  }//All Blocks
  if(HFill){
    if( fabs(MY_cand.p4().eta())<1.0 )       h_NclustAsso_EB1->Fill(nEcalClusters);
    else if( fabs(MY_cand.p4().eta())<1.5 )  h_NclustAsso_EB2->Fill(nEcalClusters);
    else if( fabs(MY_cand.p4().eta())<2.4 )  h_NclustAsso_EE1->Fill(nEcalClusters);
    else                                     h_NclustAsso_EE2->Fill(nEcalClusters);
  }
  //Now Select the more Ene cluster into the SC
  float Emin=0;
  for( int nClu=0; nClu<int(V_cluster.size()); nClu++ ){
    if( V_cluster[nClu]->energy() > Emin ){
	Emin = V_cluster[nClu]->energy(); //Take seed from most energetic cluster
	S_Time = V_seeds[nClu].time();
	float NineDiv = 0.;
	if( isEB ){
	  EBDetId IdXtal(  V_seeds[nClu].id() );
	  const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	  GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( atZero ? 0: EB_LAYER_ );
	  S_Eta = PosRechit.eta();
	  S_Phi = PosRechit.phi();
	  S_Pt  = V_seeds[nClu].energy()/cosh(  PosRechit.eta() );
	  C_Pt  = V_cluster[nClu]->energy()/cosh(V_cluster[nClu]->eta());
	  C_Eta = V_cluster[nClu]->eta();
	  C_Phi = V_cluster[nClu]->phi();
	  C_E   = V_cluster[nClu]->energy();
	  if(Use_R9_){
	    for (auto& rhEB : *recHitsEB){
		float iEta = EBDetId(rhEB.id()).ieta(), iPhi = EBDetId(rhEB.id()).iphi();
		if ( fabs( iEta-EBDetId(V_seeds[nClu].id()).ieta() ) < 2 && fabs( iPhi-EBDetId(V_seeds[nClu].id()).iphi() ) < 2 ){
		  NineDiv+=rhEB.energy();
		}
	    }
	  }
	}
	else{
	  EKDetId IdXtal(  V_seeds[nClu].id() );
	  const CaloCellGeometry* cell=geometry_->getGeometry(IdXtal);
	  GlobalPoint PosRechit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( atZero ? 0: EE_LAYER_ );
	  S_Eta = PosRechit.eta();
	  S_Phi = PosRechit.phi();
	  S_Pt  = V_seeds[nClu].energy()/cosh(  PosRechit.eta() );
	  C_Pt  = V_cluster[nClu]->energy()/cosh(V_cluster[nClu]->eta());
	  C_Eta = V_cluster[nClu]->eta();
	  C_Phi = V_cluster[nClu]->phi();
	  C_E   = V_cluster[nClu]->energy();
	  if(Use_R9_){
	    for (auto& rhEE : *recHitsEE){
		float iX = EKDetId(rhEE.id()).ix(), iY = EKDetId(rhEE.id()).iy(), iZ = EKDetId(rhEE.id()).zside();
		if ( fabs( iX-EKDetId(V_seeds[nClu].id()).ix() ) < 2 && fabs( iY-EKDetId(V_seeds[nClu].id()).iy() ) < 2 && iZ == EKDetId(V_seeds[nClu].id()).zside() ){
		  NineDiv+=rhEE.energy();
		}
	    }
	  }
	}
	erNine = Use_R9_ ? float(V_seeds[nClu].energy()/NineDiv) : 0.;
    }
  }//Loop More E cluster
  SeedInfo.push_back( S_Pt ); SeedInfo.push_back( S_Eta ); SeedInfo.push_back( S_Phi ); SeedInfo.push_back( S_Time );  SeedInfo.push_back( C_Pt );  SeedInfo.push_back( C_Eta ); SeedInfo.push_back( C_Phi ); SeedInfo.push_back( C_E ); SeedInfo.push_back( erNine );
  return SeedInfo;
}

// ------------------------------------------------------------------------------------------
void Generic_Analizer::beginJob() {
}
// ------------------------------------------------------------------------------------------
void Generic_Analizer::endJob() {
  outfile->cd();
  h_EventFlow->Write();
  h_T0->Write();
  h_BestTime_Fir_RemovalSumEt_zoom->Write();
  h_BestTime_Fir_RemovalSumEt->Write();
  h_SumEt->Write();
  h_SumEt_cut->Write();
  h_SumEt_15cut->Write();
  h_SumEt_30cut->Write();
  h_SumEt_50cut->Write();
  h_SumEt_500cut->Write();
  h_TOT_SumEt->Write();
  h_TOT_SumEt_cut->Write();
  h_TOT_SumEt_15cut->Write();
  h_TOT_SumEt_30cut->Write();
  h_TOT_SumEt_50cut->Write();
  h_TOT_SumEt_500cut->Write();
  h_Time->Write();
  h_Time_we->Write();
  h_TimeSmeared->Write();
  h_TimeSmeared_we->Write();
  h_GoodJet_t->Write();
  h_GoodJet_tEB->Write();
  h_GoodJet_tEE->Write();
  h_GoodJet_tEB2->Write();
  h_GoodJet_tEE2->Write();
  h_BadJet_t->Write();
  h_GoodGamma_t->Write();
  h_GoodGamma_tEB->Write();
  h_GoodGamma_tEE->Write();
  h_GoodGamma_tEB2->Write();
  h_GoodGamma_tEE2->Write();
  h_Phot_DR->Write();
  h_Jet_DR->Write();
  h_PtGenJet->Write();
  h_EtaGenJet->Write();
  h_EffEta_phot1->Write();
  h_EffEta_phot2->Write();
  h_EffEta_phot3->Write();
  h_NEffEta_phot1->Write();
  h_NEffEta_phot2->Write();
  h_NEffEta_phot3->Write();
  h_EffEta_jet1->Write();
  h_EffEta_jet2->Write();
  h_EffEta_jet3->Write();
  h_NEffEta_jet1->Write();
  h_NEffEta_jet2->Write();
  h_NEffEta_jet3->Write();
  h_EffEtaTot_phot->Write();
  h_EffEta_phot->Write();
  h_EffEta_jet->Write();
  h_EffPt_jet->Write();
  h_EffEtaTot_jet->Write();
  h_EffPtTot_jet->Write();
  h_BadGamma_t->Write();
  h_Rh0->Write();
  h_NVtx->Write();
  h_HiggsMass->Write();
  h_HiggsPFMass->Write();
  Ereso1->Write();
  Ereso2->Write();
  Ereso1_Mit->Write();
  Ereso2_Mit->Write();
  h_HiggsPFMass_Mit->Write();
  h_PU_frac_1->Write();
  h_PU_frac_2->Write();
  h_HiggsPFMass_Vtx->Write();
  h_HiggsMass_MC->Write();
  h_Association->SetBinContent( 1, Associated_tot);
  h_Association->SetBinContent( 2, Associated);
  h_Association->SetBinContent( 3, Associated_time);
  h_Association->Write();
  h_Association_EB->SetBinContent( 1, Associated_tot_EB);
  h_Association_EB->SetBinContent( 2, Associated_EB);
  h_Association_EB->SetBinContent( 3, Associated_time_EB);
  h_Association_EB->Write();
  h_Association_EE->SetBinContent( 1, Associated_tot_EE);
  h_Association_EE->SetBinContent( 2, Associated_EE);
  h_Association_EE->SetBinContent( 3, Associated_time_EE);
  h_Association_EE->Write();
  h_PtGammaAssoEB->Write();
  h_EGammaAssoEB->Write();
  h_DRAsso->Write();
  h_TimeGammaNOTAsso->Write();
  h_TimeGammaNOTAssoW->Write();
  h_TimeGammaAssoEB_sme15->Write();
  h_TimeGammaAssoEE_sme15->Write();
  h_TimeGammaAssoEB_sme30->Write();
  h_TimeGammaAssoEE_sme30->Write();
  h_TimeGammaAssoEB_sme50->Write();
  h_TimeGammaAssoEE_sme50->Write();
  h_TimeGammaAssoEB_sme500->Write();
  h_TimeGammaAssoEE_sme500->Write();
  h_TimeGammaAssoEB->Write();
  h_TimeGammaAssoEB_L->Write();
  h_TimeGammaAssoEB_L2->Write();
  h_PtGammaAssoEE->Write();
  h_EGammaAssoEE->Write();
  h_TimeGammaAssoEE->Write();
  h_TimeGammaAssoEE_L->Write();
  h_TimeGammaAssoEE_L2->Write();
  h_NclustAsso_EB1->Write();
  h_NclustAsso_EB2->Write();
  h_NclustAsso_EE1->Write();
  h_NclustAsso_EE2->Write();
  //
  h_DR_vs_Time_EB->Write();
  h_DR_vs_Time_EB2->Write();
  h_DR_vs_Time_EB_reb->Write();
  h_DR_vs_Time_EB_reb2->Write();
  h_DR_vs_Time_EB_b->Write();
  h_DR_vs_Time_EB_reb_b->Write();
  h_DR_vs_Time_EB_reb2_b->Write();
  h_DR_vs_Time_L_EB->Write();
  h_DR_vs_Time_L_EB_b->Write();
  h_energyForDR_EB->Write();
  h_timeForDR_EB->Write();
  h_timevsEne_EB->Write();
  h_timevsEne_EB2->Write();
  h_energyForDR_EB_b->Write();
  h_timeForDR_EB_b->Write();
  h_timevsEne_EB_b->Write();
  h_timevsEne_EB2_b->Write();
  for(int i=0; i<h_DR_vs_Time_EB->GetNbinsX(); i++){ //_reb
    TH1D *h_Time_vs_DR_RMS = h_DR_vs_Time_EB->ProjectionY( "_py", i+1, i+1 );
    ostringstream Bin; Bin << i;
    TString Name = "Time_dr_EB_" + Bin.str();
    h_Time_vs_DR_RMS->SetName( Name.Data() ); h_Time_vs_DR_RMS->SetNameTitle( Name.Data(),  Name.Data() );
    if( h_Time_vs_DR_RMS->Integral(1,h_Time_vs_DR_RMS->GetNbinsX())>200 ){
	double quant[2], value[2];
	value[0]=(h_Time_vs_DR_RMS->Integral()*99./100)/h_Time_vs_DR_RMS->Integral();
	value[1]=(h_Time_vs_DR_RMS->Integral()*0.5/100)/h_Time_vs_DR_RMS->Integral();
	h_Time_vs_DR_RMS->GetQuantiles(2, quant, value);
	float Max = quant[0]; float Min = quant[1];
	h_Time_vs_DR_RMS->GetXaxis()->SetRangeUser(Min, Max);
	h_MEANRMS_vs_DR_EB->SetBinContent( i+1, h_Time_vs_DR_RMS->GetMean() );
	h_MEANRMS_vs_DR_EB->SetBinError( i+1, h_Time_vs_DR_RMS->GetRMS() );
	h_MEANRMS_vs_DR_EB_c->SetBinContent( i+1, Min+(Max-Min)/2. );
	h_MEANRMS_vs_DR_EB_c->SetBinError( i+1, (Max-Min)/2. );
	h_Time_vs_DR_RMS->Write();
    }
  }
  h_MEANRMS_vs_DR_EB->SetMarkerColor(1); h_MEANRMS_vs_DR_EB->SetMarkerStyle(20); h_MEANRMS_vs_DR_EB->SetMarkerSize(1);
  h_MEANRMS_vs_DR_EB->Write();
  h_MEANRMS_vs_DR_EB_c->SetMarkerColor(1); h_MEANRMS_vs_DR_EB_c->SetMarkerStyle(20); h_MEANRMS_vs_DR_EB_c->SetMarkerSize(1);
  h_MEANRMS_vs_DR_EB_c->Write();
  //EE
  h_DR_vs_Time_EE->Write();
  h_DR_vs_Time_EE2->Write();
  h_DR_vs_Time_EE_reb->Write();
  h_DR_vs_Time_EE_reb2->Write();
  h_DR_vs_Time_EE_b->Write();
  h_DR_vs_Time_EE_reb_b->Write();
  h_DR_vs_Time_EE_reb2_b->Write();
  h_DR_vs_Time_L_EE->Write();
  h_DR_vs_Time_L_EE_b->Write();
  h_energyForDR_EE->Write();
  h_timeForDR_EE->Write();
  h_timevsEne_EE->Write();
  h_timevsEne_EE2->Write();
  h_energyForDR_EE_b->Write();
  h_timeForDR_EE_b->Write();
  h_timevsEne_EE_b->Write();
  h_timevsEne_EE2_b->Write();
  for(int i=0; i<h_DR_vs_Time_EE->GetNbinsX(); i++){//No _reb
    TH1D *h_Time_vs_DR_RMS = h_DR_vs_Time_EE->ProjectionY( "_py", i+1, i+1 );
    ostringstream Bin; Bin << i;
    TString Name = "Time_dr_EE_" + Bin.str();
    h_Time_vs_DR_RMS->SetName( Name.Data() ); h_Time_vs_DR_RMS->SetNameTitle( Name.Data(),  Name.Data() );
    if( h_Time_vs_DR_RMS->Integral(1,h_Time_vs_DR_RMS->GetNbinsX())>200 ){
	double quant[2], value[2];
	value[0]=(h_Time_vs_DR_RMS->Integral()*99./100)/h_Time_vs_DR_RMS->Integral();
	value[1]=(h_Time_vs_DR_RMS->Integral()*0.5/100)/h_Time_vs_DR_RMS->Integral();
	h_Time_vs_DR_RMS->GetQuantiles(2, quant, value);
	float Max = quant[0]; float Min = quant[1];
	h_Time_vs_DR_RMS->GetXaxis()->SetRangeUser(Min, Max);
	h_MEANRMS_vs_DR_EE->SetBinContent( i+1, h_Time_vs_DR_RMS->GetMean() );
	h_MEANRMS_vs_DR_EE->SetBinError( i+1, h_Time_vs_DR_RMS->GetRMS() );
	h_MEANRMS_vs_DR_EE_c->SetBinContent( i+1, Min+(Max-Min)/2. );
	h_MEANRMS_vs_DR_EE_c->SetBinError( i+1, (Max-Min)/2. );
	h_Time_vs_DR_RMS->Write();
    }
  }
  h_MEANRMS_vs_DR_EE->SetMarkerColor(1); h_MEANRMS_vs_DR_EE->SetMarkerStyle(20); h_MEANRMS_vs_DR_EE->SetMarkerSize(1);
  h_MEANRMS_vs_DR_EE->Write();
  h_MEANRMS_vs_DR_EE_c->SetMarkerColor(1); h_MEANRMS_vs_DR_EE_c->SetMarkerStyle(20); h_MEANRMS_vs_DR_EE_c->SetMarkerSize(1);
  h_MEANRMS_vs_DR_EE_c->Write();
  h_R91->Write();
  h_R92->Write();
  if( WannaFitT0Vtx_ ){
    Tree_Vtx->Write();
  }
}

// ------------------------------------------------------------------------------------------
std::vector<DetId> Generic_Analizer::getPFJetRecHitsDR(reco::PFCandidate pfCa, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEB, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEE, const edm::EventSetup& iSetup)
{
  std::vector<DetId> PFrechitsId;
  double etajet = pfCa.p4().eta();
  double phijet = pfCa.p4().phi();

  for (auto& rhEB : *recHitsEB){
    EcalRecHit myhit = (rhEB);
    // thisamp is the EB amplitude of the current rechit
    double thisamp  = myhit.energy () ;

    if (thisamp<1.) continue;

    edm::ESHandle<CaloGeometry> pGeometry ;
    iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
    const CaloGeometry * theGeometry = pGeometry.product () ;

    GlobalPoint pos = theGeometry->getPosition((myhit).detid());
    double etarechit =  pos.eta();
    double phirechit =  pos.phi();

    double dr = sqrt(Delta_phi(phijet,phirechit)*Delta_phi(phijet,phirechit)+Delta_eta(etajet,etarechit)*Delta_eta(etajet,etarechit));
    if(dr>0.5)continue;
    PFrechitsId.push_back(myhit.id());
  }

  for (auto& rhEE : *recHitsEE){
    EcalRecHit myhit = (rhEE);
    // thisamp is the EB amplitude of the current rechit
    double thisamp  = myhit.energy () ;

    if (thisamp<1.) continue;

    edm::ESHandle<CaloGeometry> pGeometry ;
    iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
    const CaloGeometry * theGeometry = pGeometry.product () ;

    GlobalPoint pos = theGeometry->getPosition((myhit).detid());
    double etarechit =  pos.eta();
    double phirechit =  pos.phi();
    double dr = sqrt(Delta_phi(phijet,phirechit)*Delta_phi(phijet,phirechit)+Delta_eta(etajet,etarechit)*Delta_eta(etajet,etarechit));
    if(dr>0.5)continue;
    PFrechitsId.push_back(myhit.id());
  }
  return PFrechitsId;
}

// ------------------------------------------------------------------------------------------
void Generic_Analizer::beginRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void Generic_Analizer::endRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void Generic_Analizer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void Generic_Analizer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}


//define this as a plug-in
DEFINE_FWK_MODULE(Generic_Analizer);
