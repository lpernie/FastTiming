#ifndef FastTiming_Generic_Analizer_Generic_Analizer_h_
#define FastTiming_Generic_Analizer_Generic_Analizer_h_
// system include files
#include <memory>
#include "TRandom3.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FastTiming/Generic_Analizer/interface/FastTool.h"

#define DEBUG
#define PI 3.14159
#define LIGHT_SPEED 29.9792458 //[cm/ns]
#define MAXRECHIT 80000
#define MAXPF 1000
using namespace std;
// ------------------------------------------------------------------------------------------
class Generic_Analizer : public edm::EDAnalyzer {

  public:
    explicit Generic_Analizer(const edm::ParameterSet&);
    ~Generic_Analizer();

    //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    typedef math::XYZTLorentzVector LorentzVector;

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void endRun(edm::Run&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    float FillLateralDevel( reco::PFCandidate, edm::Handle<edm::SortedCollection<EcalRecHit> >&, edm::Handle<edm::SortedCollection<EcalRecHit> >&, bool isSig );
    std::vector<float> GetSeedFromSC( bool isEB, reco::PFCandidate Gamma, edm::Handle<edm::SortedCollection<EcalRecHit> >&, edm::Handle<edm::SortedCollection<EcalRecHit> >&, bool atZero, bool HFill );
    float Compute_PUfrac( reco::PFCandidate pfcan, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEB, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEE );
    std::vector<DetId> getPFJetRecHitsDR(reco::PFCandidate pfCa, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEB, edm::Handle<edm::SortedCollection<EcalRecHit> >& recHitsEE, const edm::EventSetup& iSetup);
    edm::InputTag fPFCands;
    edm::InputTag SimVtx_;
    edm::InputTag Gamma_;
    edm::InputTag Jet_;
    edm::InputTag GenJet_;
    edm::InputTag GenPar_;
    edm::InputTag JetCHS_;
    edm::InputTag RecoVtx_;
    edm::InputTag ak5PFRho_;
    double EB_LAYER_;
    double EE_LAYER_;
    double smearing_;
    TRandom *ranGaus_;
    FastTool *FTool_;
    const CaloGeometry* geometry_;
    bool  debug;
    //Selec
    float MinPt_Gen;
    float MinPt_Reco;
    float MinPt_RecoPu;
    float MinDR_asso;
    float MinDR_pu;
    //
    float Associated;
    float Associated_EB;
    float Associated_EE;
    float Associated_time;
    float Associated_time_EB;
    float Associated_time_EE;
    float Associated_tot;
    float Associated_tot_EB;
    float Associated_tot_EE;
    bool isNOPUHgg_;
    bool WannaFitT0Vtx_;
    bool OneGamma_;
    bool Use_MinEne_;
    bool Use_R9_;
    bool DoSumEt_;
    bool DoMass_;
    bool isHgg_;
    double SubT0TOF_;
    string OutName_;
    TFile*outfile;
    TH1F *h_EventFlow;
    TH1F *h_T0;
    TH1F *h_BestTime_Fir_RemovalSumEt_zoom;
    TH1F *h_BestTime_Fir_RemovalSumEt;
    TH1F *h_SumEt;
    TH1F *h_SumEt_cut;
    TH1F *h_SumEt_15cut;
    TH1F *h_SumEt_30cut;
    TH1F *h_SumEt_50cut;
    TH1F *h_SumEt_500cut;
    TH1F *h_TOT_SumEt;
    TH1F *h_TOT_SumEt_cut;
    TH1F *h_TOT_SumEt_15cut;
    TH1F *h_TOT_SumEt_30cut;
    TH1F *h_TOT_SumEt_50cut;
    TH1F *h_TOT_SumEt_500cut;
    TH1F *h_HiggsMass;
    TH1F *h_HiggsPFMass;
    TH1F *Ereso1;
    TH1F *Ereso2;
    TH1F *Ereso1_Mit;
    TH1F *Ereso2_Mit;
    TH1F *h_HiggsPFMass_Mit;
    TH1F *h_PU_frac_1;
    TH1F *h_PU_frac_2;
    TH1F *h_HiggsPFMass_Vtx;
    TH1F *h_HiggsMass_MC;
    TH1F *h_Association;
    TH1F *h_Association_EB;
    TH1F *h_Association_EE;
    TH1F *h_PtGammaAssoEB;
    TH1F *h_EGammaAssoEB;
    TH1F *h_DRAsso;
    TH1F *h_TimeGammaNOTAsso;
    TH1F *h_TimeGammaNOTAssoW;
    TH1F *h_TimeGammaAssoEB_sme15;
    TH1F *h_TimeGammaAssoEE_sme15;
    TH1F *h_TimeGammaAssoEB_sme30;
    TH1F *h_TimeGammaAssoEE_sme30;
    TH1F *h_TimeGammaAssoEB_sme50;
    TH1F *h_TimeGammaAssoEE_sme50;
    TH1F *h_TimeGammaAssoEB_sme500;
    TH1F *h_TimeGammaAssoEE_sme500;
    TH1F *h_TimeGammaAssoEB;
    TH1F *h_TimeGammaAssoEB_L;
    TH1F *h_TimeGammaAssoEB_L2;
    TH1F *h_PtGammaAssoEE;
    TH1F *h_EGammaAssoEE;
    TH1F *h_TimeGammaAssoEE;
    TH1F *h_TimeGammaAssoEE_L;
    TH1F *h_TimeGammaAssoEE_L2;
    TH1F *h_Time;
    TH1F *h_Time_we;
    TH1F *h_TimeSmeared;
    TH1F *h_TimeSmeared_we;
    TH1F *h_GoodJet_t;
    TH1F *h_GoodJet_tEB;
    TH1F *h_GoodJet_tEE;
    TH1F *h_GoodJet_tEB2;
    TH1F *h_GoodJet_tEE2;
    TH1F *h_BadJet_t;
    TH1F *h_GoodGamma_t;
    TH1F *h_GoodGamma_tEB;
    TH1F *h_GoodGamma_tEE;
    TH1F *h_GoodGamma_tEB2;
    TH1F *h_GoodGamma_tEE2;
    TH1F *h_Phot_DR;
    TH1F *h_Jet_DR;
    TH1F *h_PtGenJet;
    TH1F *h_EtaGenJet;
    TH1F *h_EffEta_phot1;
    TH1F *h_EffEta_phot2;
    TH1F *h_EffEta_phot3;
    TH1F *h_NEffEta_phot1;
    TH1F *h_NEffEta_phot2;
    TH1F *h_NEffEta_phot3;
    TH1F *h_EffEta_jet1;
    TH1F *h_EffEta_jet2;
    TH1F *h_EffEta_jet3;
    TH1F *h_NEffEta_jet1;
    TH1F *h_NEffEta_jet2;
    TH1F *h_NEffEta_jet3;
    TH1F *h_EffEta_phot;
    TH1F *h_EffEtaTot_phot;
    TH1F *h_EffEtaTot_jet;
    TH1F *h_EffPtTot_jet;
    TH1F *h_EffEta_jet;
    TH1F *h_EffPt_jet;
    TH1F *h_BadGamma_t;
    TH1F *h_Rh0;
    TH1F *h_NVtx;
    TH1F *h_NclustAsso_EB1;
    TH1F *h_NclustAsso_EB2;
    TH1F *h_NclustAsso_EE1;
    TH1F *h_NclustAsso_EE2;
    TH2F *h_DR_vs_Time_EB;
    TH2F *h_DR_vs_Time_EB2;
    TH2F *h_DR_vs_Time_EB_reb;
    TH2F *h_DR_vs_Time_EB_reb2;
    TH2F *h_DR_vs_Time_EB_b;
    TH2F *h_DR_vs_Time_EB_reb_b;
    TH2F *h_DR_vs_Time_EB_reb2_b;
    TH2F *h_DR_vs_Time_L_EB;
    TH2F *h_DR_vs_Time_L_EB_b;
    TH1F *h_energyForDR_EB;
    TH1F *h_timeForDR_EB;
    TH2F *h_timevsEne_EB;
    TH2F *h_timevsEne_EB2;
    TH1F *h_energyForDR_EB_b;
    TH1F *h_timeForDR_EB_b;
    TH2F *h_timevsEne_EB_b;
    TH2F *h_timevsEne_EB2_b;
    TH1F *h_MEANRMS_vs_DR_EB;
    TH1F *h_MEANRMS_vs_DR_EB_c;
    TH2F *h_DR_vs_Time_EE;
    TH2F *h_DR_vs_Time_EE2;
    TH2F *h_DR_vs_Time_EE_reb;
    TH2F *h_DR_vs_Time_EE_reb2;
    TH2F *h_DR_vs_Time_EE_b;
    TH2F *h_DR_vs_Time_EE_reb_b;
    TH2F *h_DR_vs_Time_EE_reb2_b;
    TH2F *h_DR_vs_Time_L_EE;
    TH2F *h_DR_vs_Time_L_EE_b;
    TH1F *h_energyForDR_EE;
    TH1F *h_timeForDR_EE;
    TH2F *h_timevsEne_EE;
    TH2F *h_timevsEne_EE2;
    TH1F *h_energyForDR_EE_b;
    TH1F *h_timeForDR_EE_b;
    TH2F *h_timevsEne_EE_b;
    TH2F *h_timevsEne_EE2_b;
    TH1F *h_MEANRMS_vs_DR_EE;
    TH1F *h_MEANRMS_vs_DR_EE_c;
    TH1F *h_R91;
    TH1F *h_R92;
    TTree * Tree_Vtx;
    Float_t VtxDet_T1;          
    Float_t VtxDet_VBFT1;          
    Float_t VtxDet_GT1;          
    Float_t VtxDet_VBFGT1;          
    Float_t VtxDet_T2;          
    Float_t VtxDet_VBFT2;          
    Float_t VtxDet_GT2;          
    Float_t VtxDet_VBFGT2;          
    Float_t VtxDet_time;        
    Float_t VtxDet_PosXtal_X1;  
    Float_t VtxDet_PosXtal_Y1;  
    Float_t VtxDet_PosXtal_Z1;  
    Float_t VtxDet_PosXtal_X2;  
    Float_t VtxDet_PosXtal_Y2;  
    Float_t VtxDet_PosXtal_Z2;  
    Float_t VtxDet_PosVBF_X1;  
    Float_t VtxDet_PosVBF_Y1;  
    Float_t VtxDet_PosVBF_Z1;  
    Float_t VtxDet_PosVBF_X2;  
    Float_t VtxDet_PosVBF_Y2;  
    Float_t VtxDet_PosVBF_Z2;  
    Float_t VtxDet_PosXtal_MCX1;  
    Float_t VtxDet_PosXtal_MCY1;  
    Float_t VtxDet_PosXtal_MCZ1;  
    Float_t VtxDet_PosXtal_MCX2;  
    Float_t VtxDet_PosXtal_MCY2;  
    Float_t VtxDet_PosXtal_MCZ2;  
    Float_t VtxDet_PtRecoJet_1; 
    Float_t VtxDet_PtRecoJet_2; 
    Float_t VtxDet_PtMCJet_1;   
    Float_t VtxDet_PtMCJet_2;   
    Float_t vzMC;               
    Float_t vyMC;               
    Float_t vxMC;               
    Float_t vzRECO;               
    Float_t vyRECO;               
    Float_t vxRECO;               
};
#endif
