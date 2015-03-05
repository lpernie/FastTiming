#include <stdlib.h>
#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <TMinuit.h>

#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TProfile.h"
#include "TChain.h"
#include "TMath.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TFitter.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "tdrStyle.C"

//#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#define NSTEPMAX 10000
#define PI 3.14159
using namespace RooFit;
using namespace std;

void chi2Hist(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  Double_t chi2(0.);
  TVector3 Vtx;     Vtx.SetXYZ( par[12], par[13], par[10]);
  TVector3 posExp1; posExp1.SetXYZ(par[1], par[2], par[3] );
  TVector3 posExp2; posExp2.SetXYZ(par[6], par[7], par[8] );
  posExp1 -= Vtx; posExp2 -= Vtx;
  float distExp1 = sqrt( pow(posExp1.X(),2) + pow(posExp1.Y(),2) + pow(posExp1.Z(),2) );
  float distExp2 = sqrt( pow(posExp2.X(),2) + pow(posExp2.Y(),2) + pow(posExp2.Z(),2) );
  float Texp1 = par[11] + (distExp1)/29.9792458;  //Time c [cm/ns]
  float Texp2 = par[11] + (distExp2)/29.9792458;  //Time c [cm/ns]
  chi2 = pow(par[0]-Texp1, 2)/pow(par[4],2) + pow(par[5]-Texp2, 2)/pow(par[9],2);

  f=chi2;
}

vector<float> ComputeVertex( float Gtime1, float Gtime2, float sigma_T1, float sigma_T2, TVector3 posExp1, TVector3 posExp2, bool T0Free, float VtxDet_timem, float Vx, float Vy);
float Makemin( float a, float b);
float Makemax( float a, float b);

// NameFile     = INPUT FILE
// T0Free       = T0 parameter have to be free to variate?
// OnlyEB       = 0: EB and EE together - 1: only EB - 2: only EE
// OutPutFolder = Name of the output folder
//.x VertexDeterminator_minuit.C+("../GenericAnalyzer_Higgs_140.root", false, 0, "VTX_Output_T0FIX_EBEE_140PU", "std") 
//void VertexDeterminator_minuit( TString NameFile, bool T0Free = true, int OnlyEB = 0, TString OutPutFolder = "VTX_Output", TString sele = "std","0" ){
void VertexDeterminator_minuit( TString NameFile, bool T0Free , int OnlyEB , TString OutPutFolder, TString sele = "std" ){

  //Initial Stuffs
  if(T0Free) cout<<"----> START VertexDeterminator. T0 is Free"<<endl;
  else       cout<<"----> START VertexDeterminator. T0 is Fix"<<endl;
  setTDRStyle();
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);

  //Get Input files
  TString Comm = "mkdir -p " + OutPutFolder;
  system( Comm.Data() );
  TFile* File = TFile::Open( NameFile.Data() );
  if( !File ) cout<<"WARNING: File "<<NameFile.Data()<<" do not exist"<<endl;
  TString TreeName = "Tree_Vtx";
  TTree* Tree = (TTree*) File->Get( TreeName.Data() );
  if( !Tree ) cout<<"WARNING: Tree "<<TreeName.Data()<<" do not exist"<<endl;
  myc1->cd();
  gStyle->SetPalette(1);

  //Get Variables
  Float_t VtxDet_T1, VtxDet_T2, VtxDet_GT1, VtxDet_GT2, VtxDet_time;
  Float_t VtxDet_PosXtal_X1, VtxDet_PosXtal_Y1, VtxDet_PosXtal_Z1, VtxDet_PosXtal_MCX1, VtxDet_PosXtal_MCY1, VtxDet_PosXtal_MCZ1, VtxDet_PtRecoJet_1, VtxDet_PtMCJet_1;
  Float_t VtxDet_PosXtal_X2, VtxDet_PosXtal_Y2, VtxDet_PosXtal_Z2, VtxDet_PosXtal_MCX2, VtxDet_PosXtal_MCY2, VtxDet_PosXtal_MCZ2, VtxDet_PtRecoJet_2, VtxDet_PtMCJet_2;
  Float_t vzMC, vxMC, vyMC, vzRECO, vxRECO, vyRECO;
  Tree->SetBranchAddress( "VtxDet_T1", &VtxDet_T1);
  Tree->SetBranchAddress( "VtxDet_GT1", &VtxDet_GT1);
  Tree->SetBranchAddress( "VtxDet_T2", &VtxDet_T2);
  Tree->SetBranchAddress( "VtxDet_GT2", &VtxDet_GT2);
  Tree->SetBranchAddress( "VtxDet_time", &VtxDet_time);
  Tree->SetBranchAddress( "VtxDet_PosXtal_X1", &VtxDet_PosXtal_X1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_Y1", &VtxDet_PosXtal_Y1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_Z1", &VtxDet_PosXtal_Z1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_X2", &VtxDet_PosXtal_X2);
  Tree->SetBranchAddress( "VtxDet_PosXtal_Y2", &VtxDet_PosXtal_Y2);
  Tree->SetBranchAddress( "VtxDet_PosXtal_Z2", &VtxDet_PosXtal_Z2);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCX1", &VtxDet_PosXtal_MCX1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCY1", &VtxDet_PosXtal_MCY1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCZ1", &VtxDet_PosXtal_MCZ1);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCX2", &VtxDet_PosXtal_MCX2);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCY2", &VtxDet_PosXtal_MCY2);
  Tree->SetBranchAddress( "VtxDet_PosXtal_MCZ2", &VtxDet_PosXtal_MCZ2);
  Tree->SetBranchAddress( "VtxDet_PtRecoJet_1",&VtxDet_PtRecoJet_1);
  Tree->SetBranchAddress( "VtxDet_PtRecoJet_2",&VtxDet_PtRecoJet_2);
  Tree->SetBranchAddress( "VtxDet_PtMCJet_1",  &VtxDet_PtMCJet_1);
  Tree->SetBranchAddress( "VtxDet_PtMCJet_2",  &VtxDet_PtMCJet_2);
  Tree->SetBranchAddress( "vzMC", &vzMC);
  Tree->SetBranchAddress( "vxMC", &vxMC);
  Tree->SetBranchAddress( "vyMC", &vyMC);
  Tree->SetBranchAddress( "vzRECO", &vzRECO);
  Tree->SetBranchAddress( "vxRECO", &vxRECO);
  Tree->SetBranchAddress( "vyRECO", &vyRECO);

  //Histos
  TH1F *ev_flow            = new TH1F("ev_flow","Event Flow", 8, -0.5, 7.5);
  ev_flow->GetXaxis()->SetBinLabel(1,"Total Ev."); ev_flow->GetXaxis()->SetBinLabel(2,"Effective Seeds"); ev_flow->GetXaxis()->SetBinLabel(3,"Resolution 10%");
  ev_flow->GetXaxis()->SetBinLabel(4,"Pt > 5 GeV"); ev_flow->GetXaxis()->SetBinLabel(5,"#Delta #eta"); ev_flow->GetXaxis()->SetBinLabel(6,"Eta Cuts"); ev_flow->GetXaxis()->SetBinLabel(7,"Vtx recon"); ev_flow->GetXaxis()->SetBinLabel(8,"EtaPhi recon");
  TH1F *hdr                = new TH1F("hdr","#Delta R #gamma_1 #gamma_2", 100, 0., 6.); hdr->GetXaxis()->SetTitle("#Delta R");
  TH1F *hdeta              = new TH1F("hdeta","#Delta #eta #gamma_1 #gamma_2", 100, 0., 6.); hdeta->GetXaxis()->SetTitle("#Delta #eta");
  TH1F *T1                 = new TH1F("T1","First Cluster Time", 100, -0.6, 0.6); T1->GetXaxis()->SetTitle("Time [ns]");
  TH1F *T2                 = new TH1F("T2","Second Cluster Time", 100, -0.6, 0.6); T2->GetXaxis()->SetTitle("Time [ns]");
  TH1F *XtalsEta           = new TH1F("XtalsEta","#eta Xtals", 100, -3, 3.); XtalsEta->GetXaxis()->SetTitle("#eta");
  TH1F *XtalsPhi           = new TH1F("XtalsPhi","#phi Xtals", 100, -3.14, 3.14); XtalsPhi->GetXaxis()->SetTitle("#phi");
  TH1F *isEBEE             = new TH1F("isEBEE","", 4, -0.5, 3.5); isEBEE->GetXaxis()->SetBinLabel(1,"ALL"); isEBEE->GetXaxis()->SetBinLabel(2,"EB"); isEBEE->GetXaxis()->SetBinLabel(3,"EE"); isEBEE->GetXaxis()->SetBinLabel(4,"MIX");
  TH2F *Correl_Resol_Vtx   = new TH2F("Correl_Resol_Vtx","", 100, 0.0, 0.04, 600, -8., 8.); Correl_Resol_Vtx->GetXaxis()->SetTitle("Time Resolution"); Correl_Resol_Vtx->GetYaxis()->SetTitle("Vtx Resolution");
  TH2F *Correl_Resol_T0    = new TH2F("Correl_Resol_T0","", 100, 0.0, 0.04, 600, -0.4, 0.4); Correl_Resol_T0->GetXaxis()->SetTitle("Time Resolution"); Correl_Resol_T0->GetYaxis()->SetTitle("T0 Resolution");
  TH1F *VertexReso_z       = new TH1F("VertexReso_z","", 100, -0.05, 0.05); VertexReso_z->GetXaxis()->SetTitle("vtx Z_MC - Z_RECO");
  TH1F *VertexReso_y       = new TH1F("VertexReso_y","", 100, -0.05, 0.05); VertexReso_y->GetXaxis()->SetTitle("vtx Y_MC - Y_RECO");
  TH1F *VertexReso_x       = new TH1F("VertexReso_x","", 100, -0.05, 0.05); VertexReso_x->GetXaxis()->SetTitle("vtx X_MC - X_RECO");
  TH1F *h_DeltaX1          = new TH1F("h_DeltaX1","", 100, -50, 50); h_DeltaX1->GetXaxis()->SetTitle("X1_MC - X1_RECO");
  TH1F *h_DeltaY1          = new TH1F("h_DeltaY1","", 100, -50, 50); h_DeltaY1->GetXaxis()->SetTitle("Y1_MC - Y1_RECO");
  TH1F *h_DeltaZ1          = new TH1F("h_DeltaZ1","", 100, -50, 50); h_DeltaZ1->GetXaxis()->SetTitle("Z1_MC - Z1_RECO");
  TH1F *h_DeltaEta1        = new TH1F("h_DeltaEta1","", 100, -0.2, 0.2); h_DeltaEta1->GetXaxis()->SetTitle("Eta1_MC - Eta1_RECO");
  TH1F *h_DeltaPhi1        = new TH1F("h_DeltaPhi1","", 100, -0.2, 0.2); h_DeltaPhi1->GetXaxis()->SetTitle("Phi1_MC - Phi1_RECO");
  TH1F *h_DeltaX2          = new TH1F("h_DeltaX2","", 100, -50, 50); h_DeltaX2->GetXaxis()->SetTitle("X2_MC - X2_RECO");
  TH1F *h_DeltaY2          = new TH1F("h_DeltaY2","", 100, -50, 50); h_DeltaY2->GetXaxis()->SetTitle("Y2_MC - Y2_RECO");
  TH1F *h_DeltaZ2          = new TH1F("h_DeltaZ2","", 100, -50, 50); h_DeltaZ2->GetXaxis()->SetTitle("Z2_MC - Z2_RECO");
  TH1F *h_DeltaEta2        = new TH1F("h_DeltaEta2","", 100, -0.2, 0.2); h_DeltaEta2->GetXaxis()->SetTitle("Eta2_MC - Eta2_RECO");
  TH1F *h_DeltaPhi2        = new TH1F("h_DeltaPhi2","", 100, -0.2, 0.2); h_DeltaPhi2->GetXaxis()->SetTitle("Phi2_MC - Phi2_RECO");
  TH2F *h_DeltaEta_vs_Reso = new TH2F("h_DeltaEta_vs_Reso","", 100, -0, 0.2, 100, -8., 8.); h_DeltaEta_vs_Reso->GetXaxis()->SetTitle("#Delta Eta"); h_DeltaEta_vs_Reso->GetYaxis()->SetTitle("#Delta Vtx [cm]");
  TH2F *h_DeltaPhi_vs_Reso = new TH2F("h_DeltaPhi_vs_Reso","", 100, -0, 0.2, 100, -8., 8.); h_DeltaPhi_vs_Reso->GetXaxis()->SetTitle("#Delta #Phi"); h_DeltaPhi_vs_Reso->GetYaxis()->SetTitle("#Delta Vtx [cm]");
  TH2F *h_DeltaR_vs_Reso   = new TH2F("h_DeltaR_vs_Reso","", 100, 0., 0.2, 100, -8., 8.); h_DeltaR_vs_Reso->GetXaxis()->SetTitle("#Delta R"); h_DeltaR_vs_Reso->GetYaxis()->SetTitle("#Delta Vtx [cm]");
  TH2F *h_DETA12_vs_Reso   = new TH2F("h_DETA12_vs_Reso","", 100, 0., 5., 100, -8., 8.); h_DETA12_vs_Reso->GetXaxis()->SetTitle("#Delta Eta 1-2"); h_DETA12_vs_Reso->GetYaxis()->SetTitle("#Delta Vtx [cm]");


  //Detector Smearing
  int NumSmearing=18;
  //int NumSmearing=1;
  for(int resInd=0; resInd<NumSmearing; resInd++ ){

    float MinSmearing = 0.0019;
    cout<<" Resolution Number "<<resInd<<endl;
    stringstream SmearInd;
    SmearInd << resInd;

    TH1F *Final_Vtx    = new TH1F("Final_Vtx","Vtx Reco", 600, -8., 8.); Final_Vtx->GetXaxis()->SetTitle("Z [mm]");
    TH1F *DeltaVtx     = new TH1F("DeltaVtx","Vtx Reco - Vtx MC", 600, -8., 8.); DeltaVtx->GetXaxis()->SetTitle("Z [mm]");
    TH1F *DeltaT0      = new TH1F("DeltaT0","T0 Reco - Vtx T0", 100, -0.2, 0.2); DeltaT0->GetXaxis()->SetTitle("T0 [ns]");
    TH1F *hBest_t0     = new TH1F("hBest_t0","Time 0", 600, -0.1, -0.1); hBest_t0->GetXaxis()->SetTitle("Time [ns]");
    TH2F *VtxMc_VtxRec = new TH2F("VtxMc_VtxRec","", 100, -15., 15., 100, -15., 15.); VtxMc_VtxRec->GetXaxis()->SetTitle("Z_MC [cm]"); VtxMc_VtxRec->GetYaxis()->SetTitle("Z_Reco [cm]");
    TH1F *Pull_Vtx     = new TH1F("Pull_Vtx","Pull Vertex Z", 100, -3., 3.); Pull_Vtx->GetXaxis()->SetTitle("(VtxMC-VtxReco)/Error");
    TH1F *Pull_T0      = new TH1F("Pull_T0","Pull T0", 100, -3., 3.); Pull_T0->GetXaxis()->SetTitle("(T0MC-T0Reco)/Error");

    //Loop
    Long64_t nentries = Tree->GetEntriesFast();
    cout<<" Starting Loop on "<<nentries<<" events."<<endl;
    for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
	//Beginning
	if( nentries>500 && iEntry%100==0 ) cout<<" You are at the "<<iEntry<<"/"<<nentries<<" event"<<endl;
	if( nentries<500 && iEntry%10==0  ) cout<<" You are at the "<<iEntry<<"/"<<nentries<<" event"<<endl;
	Tree->GetEntry(iEntry);
	if( resInd==0 ) ev_flow->Fill(0.);
	VtxDet_T1-=VtxDet_time; VtxDet_T2-=VtxDet_time;
	if( VtxDet_T1<-900. || VtxDet_PosXtal_X1<-900. || VtxDet_PosXtal_Y1<-900. || VtxDet_PosXtal_Z1<-900. || VtxDet_PtRecoJet_1<-900. || VtxDet_PtMCJet_1<-900. ) continue;
	if( VtxDet_T2<-900. || VtxDet_PosXtal_X2<-900. || VtxDet_PosXtal_Y2<-900. || VtxDet_PosXtal_Z2<-900. || VtxDet_PtRecoJet_2<-900. || VtxDet_PtMCJet_2<-900. ) continue;
	if( resInd==0 ) ev_flow->Fill(1.);
	//Resolution
	if( (VtxDet_PtRecoJet_1/VtxDet_PtMCJet_1<0.8 || VtxDet_PtRecoJet_1/VtxDet_PtMCJet_1>1.2) || (VtxDet_PtRecoJet_2/VtxDet_PtMCJet_2<0.8 || VtxDet_PtRecoJet_2/VtxDet_PtMCJet_2>1.2) ) continue;
	if( resInd==0 ) ev_flow->Fill(2.);
	//PT Cut
	TVector3 pos1; pos1.SetXYZ(VtxDet_PosXtal_X1, VtxDet_PosXtal_Y1, VtxDet_PosXtal_Z1);
	TVector3 pos2; pos2.SetXYZ(VtxDet_PosXtal_X2, VtxDet_PosXtal_Y2, VtxDet_PosXtal_Z2);
	TVector3 pos1MC; pos1MC.SetXYZ(VtxDet_PosXtal_MCX1, VtxDet_PosXtal_MCY1, VtxDet_PosXtal_MCZ1);
	TVector3 pos2MC; pos2MC.SetXYZ(VtxDet_PosXtal_MCX2, VtxDet_PosXtal_MCY2, VtxDet_PosXtal_MCZ2);
	float PtCut = 15;
	if(sele=="0") PtCut = 0;
	if( VtxDet_PtRecoJet_1 < PtCut || VtxDet_PtRecoJet_2 < PtCut ) continue;
	if( resInd==0 ) ev_flow->Fill(3.);
	//DEta Cut
	double deltaPhi = pos1.Phi() - pos2.Phi() ;
	double deltaEta  = pos1.Eta() - pos2.Eta();
	if (deltaPhi > PI)  deltaPhi -= 2.*PI;
	if (deltaPhi < -PI) deltaPhi += 2.*PI;
	double deltaR = std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
	if( resInd==0 ) hdr->Fill(deltaR);
	if( resInd==0 ) hdeta->Fill( fabs(deltaEta) );
	if( fabs(deltaEta) < (OnlyEB==2 ? 1.6 : 0.6) ) continue;
	if( resInd==0 ) ev_flow->Fill(4.);
	//All or EB or EE
	if( resInd==0 ){
	  T1->Fill(VtxDet_T1); T2->Fill(VtxDet_T2);
	  XtalsEta->Fill( pos1.Eta() ); XtalsEta->Fill( pos2.Eta() ); 
	  XtalsPhi->Fill( pos1.Phi() ); XtalsPhi->Fill( pos2.Phi() );
	  isEBEE->Fill(0); 
	  if( fabs(pos1.Eta())<1.5 && fabs(pos2.Eta())<1.5  )      isEBEE->Fill(1); 
	  else if( fabs(pos1.Eta())>1.5 && fabs(pos2.Eta())>1.5  ) isEBEE->Fill(2); 
	  else                                                     isEBEE->Fill(3); 
	}
	if( OnlyEB==1 && (fabs(pos1.Eta())>1.5 || fabs(pos2.Eta())>1.5) ) continue;
	if( OnlyEB==2 && fabs(pos1.Eta())<1.5 && fabs(pos2.Eta())<1.5 ) continue;
	//if( OnlyEB==2 && (fabs(pos1.Eta())<1.5 || fabs(pos2.Eta())<1.5 ) ) continue;
	if( resInd==0 ) ev_flow->Fill(5.);
	//vertex Reco
	if( resInd==0 ) VertexReso_z->Fill( vzRECO-vzMC );
	if( resInd==0 ) VertexReso_y->Fill( vyRECO-vyMC );
	if( resInd==0 ) VertexReso_x->Fill( vxRECO-vxMC );
	float Vtx_cut = 0.02;
	if(sele=="0") Vtx_cut = 99.;
	if( fabs(vzRECO-vzMC)>Vtx_cut ) continue;
	if( resInd==0 ) ev_flow->Fill(6.);
	//Etahi reconstruction
	if( resInd==0 ){
	  h_DeltaX1->Fill( VtxDet_PosXtal_MCX1-VtxDet_PosXtal_X1);
	  h_DeltaY1->Fill( VtxDet_PosXtal_MCY1-VtxDet_PosXtal_Y1);
	  h_DeltaZ1->Fill( VtxDet_PosXtal_MCZ1-VtxDet_PosXtal_Z1);
	  h_DeltaEta1->Fill( pos1MC.Eta()-pos1.Eta() );
	  h_DeltaPhi1->Fill( pos1MC.Phi()-pos1.Phi() );
	  h_DeltaX2->Fill( VtxDet_PosXtal_MCX2-VtxDet_PosXtal_X2);
	  h_DeltaY2->Fill( VtxDet_PosXtal_MCY2-VtxDet_PosXtal_Y2);
	  h_DeltaZ2->Fill( VtxDet_PosXtal_MCZ2-VtxDet_PosXtal_Z2);
	  h_DeltaEta2->Fill( pos2MC.Eta()-pos2.Eta() );
	  h_DeltaPhi2->Fill( pos2MC.Phi()-pos2.Phi() );
	}
	float MaxEta = Makemax( fabs(pos1MC.Eta()-pos1.Eta()), fabs(pos2MC.Eta()-pos2.Eta()) );
	float MaxPhi = Makemax( fabs(pos1MC.Phi()-pos1.Phi()), fabs(pos2MC.Phi()-pos2.Phi()) );
	bool passEtaPhiReco( MaxEta < 0.05 && MaxPhi < 0.025 );
	if( !passEtaPhiReco && sele!="0" ) continue;
	if( resInd==0 ) ev_flow->Fill(7.);

	//Random Smearing
	TRandom *ran = new TRandom(0);
	float effSigma = MinSmearing + ((MinSmearing*float(resInd*100.))/100.);
	float SmearedTime1 = ran->Gaus(VtxDet_T1, effSigma);
	float SmearedTime2 = ran->Gaus(VtxDet_T2, effSigma);
	if( resInd==0 ){
	  SmearedTime1 = VtxDet_T1;
	  SmearedTime2 = VtxDet_T2;
	  effSigma = MinSmearing;
	}
	delete ran;

	//Compute Eta_TOF
	float Gtime1 = SmearedTime1 + VtxDet_GT1, Gtime2 = SmearedTime2 + VtxDet_GT2;

	double sigma_T1( sqrt(effSigma*effSigma + MinSmearing*MinSmearing) );
	double sigma_T2( sqrt(effSigma*effSigma + MinSmearing*MinSmearing) );

	//Vertexing
	vector<float> Results = ComputeVertex( Gtime1, Gtime2, sigma_T1, sigma_T2, pos1, pos2, T0Free, VtxDet_time, vxMC, vyMC);
	float Best_Vz = Results[0], Vtx_err = Results[1];
	float Best_t0 = Results[2], T0_err = Results[3];

	//Pull
	Pull_Vtx->Fill( (Best_Vz-vzMC)/(Vtx_err) );
	Pull_T0->Fill( (Best_t0-VtxDet_time)/(T0_err) );

	//Final Plots
	Final_Vtx->Fill(Best_Vz);
	DeltaVtx->Fill(Best_Vz-vzMC);
	DeltaT0->Fill(Best_t0-VtxDet_time);
	VtxMc_VtxRec->Fill(vzMC, Best_Vz);
	hBest_t0->Fill(Best_t0);
	if( resInd==0 ) h_DeltaEta_vs_Reso->Fill( Makemax( fabs(pos1MC.Eta()-pos1.Eta()), fabs(pos2MC.Eta()-pos2.Eta()) ), Best_Vz-vzMC );
	if( resInd==0 ) h_DeltaPhi_vs_Reso->Fill( Makemax( fabs(pos1MC.Phi()-pos1.Phi()), fabs(pos2MC.Phi()-pos2.Phi()) ), Best_Vz-vzMC );
	if( resInd==0 ) h_DeltaR_vs_Reso->Fill( Makemax( fabs(pos1MC.DeltaR(pos1)), fabs(pos2MC.DeltaR(pos2)) ), Best_Vz-vzMC );
	if( resInd==0 ) h_DETA12_vs_Reso->Fill( fabs(pos1.Eta() - pos2.Eta() ), Best_Vz-vzMC );
	Correl_Resol_Vtx->Fill(sigma_T1, Best_Vz-vzMC );
	Correl_Resol_T0->Fill(sigma_T1, Best_t0-VtxDet_time );
    }//All events
    string smr = SmearInd.str();
    TString Hname = OutPutFolder + "/Final_Vertex_" + TString(smr) + ".png";
    Final_Vtx->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/Delta_Vertex_" + TString(smr) + ".png";
    DeltaVtx->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/Delta_T0_" + TString(smr) + ".png";
    DeltaT0->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/Time0_" + TString(smr) + ".png";
    hBest_t0->Draw(); gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/VtxMc_VtxRec_" + TString(smr) + ".png";
    VtxMc_VtxRec->Draw("colz"); gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/Pull_Vtx_" + TString(smr) + ".png";
    Pull_Vtx->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );

    Hname = OutPutFolder + "/Pull_T0_" + TString(smr) + ".png";
    Pull_T0->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );

    delete Final_Vtx;
    delete DeltaVtx;
    delete DeltaT0;
    delete hBest_t0;
    delete VtxMc_VtxRec;
    delete Pull_Vtx;
    delete Pull_T0;

  }//Resol
  TString Hname = OutPutFolder + "/Event_flow.png";
  ev_flow->SetMinimum(0.);
  ev_flow->Draw(); myc1->SaveAs( Hname.Data() ); 

  Hname = OutPutFolder + "/DR.png";
  hdr->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/DEta.png";
  hdeta->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/Original_T1.png";
  T1->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/Original_T2.png";
  T2->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/XtalsEta.png";
  XtalsEta->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/isEBEE.png";
  isEBEE->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/XtalsPhi.png";
  XtalsPhi->Draw(); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/VertexReso_z.png";
  VertexReso_z->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/VertexReso_y.png";
  VertexReso_y->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/VertexReso_x.png";
  VertexReso_x->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);

  Hname = OutPutFolder + "/DeltaX1.png";
  h_DeltaX1->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaY1.png";
  h_DeltaY1->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaZ1.png";
  h_DeltaZ1->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaEta1.png";
  h_DeltaEta1->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaPhi1.png";
  h_DeltaPhi1->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaX2.png";
  h_DeltaX2->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaY2.png";
  h_DeltaY2->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaZ2.png";
  h_DeltaZ2->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaEta2.png";
  h_DeltaEta2->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);
  Hname = OutPutFolder + "/DeltaPhi2.png";
  h_DeltaPhi2->Draw(); gStyle->SetOptStat(111111); myc1->SaveAs( Hname.Data() );  gStyle->SetOptStat(0);

  Hname = OutPutFolder + "/DeltaEta_vs_Reso.png";
  h_DeltaEta_vs_Reso->Draw("colz"); myc1->SaveAs( Hname.Data() );
  Hname = OutPutFolder + "/DeltaPhi_vs_Reso.png";
  h_DeltaPhi_vs_Reso->Draw("colz"); myc1->SaveAs( Hname.Data() );
  Hname = OutPutFolder + "/DeltaR_vs_Reso.png";
  h_DeltaR_vs_Reso->Draw("colz"); myc1->SaveAs( Hname.Data() );
  Hname = OutPutFolder + "/DeltaEta1_2_vs_Reso.png";
  h_DETA12_vs_Reso->Draw("colz"); myc1->SaveAs( Hname.Data() );

  bool XminXmax = false;
  bool XminXmax0 = false;
  float QuantMin = T0Free ? 15. : 2.5, QuantMax = T0Free ? 85. : 97.5;

  Hname = OutPutFolder + "/Correl_Resol_Vtx.png";
  Correl_Resol_Vtx->Draw("colz"); gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );
  TH1D *VtxRMS = new TH1D("VtxRMS", "", 100, 0.001, 0.04);
  int index(0);
  for(int i=0; i<Correl_Resol_Vtx->GetNbinsX(); i++){
cout<<"AAAA!!1 "<<i<<endl;
    if( Correl_Resol_Vtx->ProjectionY(" ",i+1,i+1)->Integral()>0. ){
cout<<"AAAA!!2 "<<i<<endl;
	double quant[2], value[2];
	TH1D *h1 = Correl_Resol_Vtx->ProjectionY(" ",i+1,i+1);
	//for(int NN=0; NN<h1->GetNbinsX(); NN++){ if( h1->GetBinContent(NN+1)>1 ) h1->SetBinContent(NN+1, h1->GetBinContent(NN+1) ); }
	value[0]=(h1->Integral()*QuantMax/100)/h1->Integral();
	value[1]=(h1->Integral()*QuantMin/100)/h1->Integral();
	h1->GetQuantiles(2, quant, value);
	float Max = quant[0], Min = quant[1];
	if(OnlyEB==1) h1->Rebin(4);
	if(OnlyEB==2){ h1->Rebin(2); }//if(!XminXmax0 || (XminXmax0 && i>25) ) h1->Rebin(2); }
	h1->GetXaxis()->SetRangeUser(Min, Max);
//	if(i==6  && OnlyEB==2 )  h1->GetXaxis()->SetRangeUser(-0.3, 0.3);
//	if(i==6  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-1., 1.);
//	if(i==10  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-1., 1.);
//	if(i==15  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-2., 2.);
//	if(i==19  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-2., 2.);
//	if(i==24  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-2., 2.);
//	if(i==28  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-3., 2.);
//	if(i==33  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-2.5, 2.5);
//	if(i==6  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-0.5, 0.5);
//	if(i==10  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-0.65, 0.65);
//	if(i==15  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-0.8, 0.8);
//	if(i==19  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-1., 1.);
//	if(i==24  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-1., 1.);
//	if(i==28  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-1.5, 1.);
//	if(i==33  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-1.5, 1.5);
//	if(i==38  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-1.8, 1.8);
//	if(i==43  && OnlyEB==2 && XminXmax0) h1->GetXaxis()->SetRangeUser(-2., 2.);

	TF1 *MyRms = new TF1("MyRms","gaus", Min+0.01, Max-0.01);
	MyRms->SetParName(0,"Normalizzazione"); MyRms->SetParameters(0,h1->Integral());
	MyRms->SetParName(1,"Mean");            MyRms->SetParameters(1, h1->GetBinCenter( h1->GetMaximumBin() ) );
	MyRms->SetParName(2,"Sigma");           MyRms->SetParameters(2, 0.02 + (index*1*0.01) );
	//MyRms->SetParName(2,"Sigma");           MyRms->SetParameters(2, 0.02 + (index*0.005) );
	h1->Fit("MyRms");
	index++;

	//	RooRealVar x("x","X", Min+0.01, Max-0.01, "");
	//	RooDataHist dh("dh","Resolution",RooArgList(x),h1);
	//	RooRealVar mean("mean","mean", 0., -0.1,0.1,"");
	//	RooRealVar sigma("sigma","#sigma", 27*Correl_Resol_Vtx->GetBinCenter(i+1), 0., 3. ,"");
	//	mean.setRange(-0.1, 0.1);
	//	sigma.setRange(0., 3.);
	//	RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);
	//	RooRealVar Nsig("Nsig","Yield",1000.,0.,1.e7);
	//	Nsig.setVal( h1->Integral() );
	//	Nsig.setRange(h1->Integral()/50., 1.5*h1->Integral());
	//	RooAddPdf model1("model","only_gaus",RooArgList(gaus),RooArgList(Nsig));
	//	RooAbsPdf* model=0; model = &model1;
	//	RooNLLVar nll("nll","log likelihood var",*model,dh,Extended());
	//	RooMinuit m(nll);
	//	m.setVerbose(kFALSE);
	//	m.migrad();
	//	RooFitResult* res = m.save();

cout<<"AAAA!!3 "<<i<<endl;
	VtxRMS->SetBinContent( i+1, Makemin( h1->GetRMS(), MyRms->GetParameter(2)) );
	VtxRMS->SetBinError( i+1, h1->GetRMSError() );
	stringstream Ind; Ind << i; string Indst = Ind.str();cout<<"AAAAOOOO "<<Indst<<endl; h1->Draw(); Hname = OutPutFolder + "/VTXRMS_" + TString(Indst)  + ".png"; gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );
	delete h1;
	//delete MyRms;
    }
  }
  //TF1 *Mypol1 = new TF1("Mypol1","[0]+[1]*x",0.0, 0.04);
  TF1 *Mypol1 = new TF1("Mypol1","[0]*x",0.0, 0.04);
  //Mypol1->SetParName(0,"int.");
  Mypol1->SetParName(0,"coeff.");
  VtxRMS->Fit("Mypol1");
  VtxRMS->GetXaxis()->SetTitle("Time Resolution [ns]"); VtxRMS->GetYaxis()->SetTitle("Vertex Resolution [cm]");
  VtxRMS->SetMinimum(0.); VtxRMS->SetMarkerStyle(34); VtxRMS->SetMarkerColor(kRed);VtxRMS->SetLineColor(kRed); VtxRMS->SetFillColor(kRed); VtxRMS->Draw("P");
  TLatex lat; char line[300];
  lat.SetNDC(); lat.SetTextSize(0.040); lat.SetTextColor(1);
  float xmin(0.45), yhi(0.80), ypass(0.05);
  sprintf(line,"Int: %.2f [cm]", Mypol1->GetParameter(0) );
  lat.DrawLatex(xmin,yhi-ypass, line);
  sprintf(line,"Coeff: %.2f", Mypol1->GetParameter(1) );
  lat.DrawLatex(xmin,yhi-2.*ypass, line );
  Hname = OutPutFolder + "/Correl_Resol_Vtx_prfRMS.png";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );
  Hname = OutPutFolder + "/Correl_Resol_Vtx_prfRMS.C";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );

  Hname = OutPutFolder + "/Correl_Resol_T0.png";
  Correl_Resol_T0->Draw("colz"); gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );
  TH1D *T0RMS = new TH1D("T0RMS", "", 100, 0.001, 0.04);
  for(int i=0; i<Correl_Resol_T0->GetNbinsX(); i++){
    if( Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->Integral()>0. ){
	double quant[2], value[2];
	value[0]=(Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->Integral()*QuantMax/100)/Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->Integral();
	value[1]=(Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->Integral()*QuantMin/100)/Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->Integral();
	Correl_Resol_T0->ProjectionY(" ",i+1,i+1)->GetQuantiles(2, quant, value);
	float Max = quant[0], Min = quant[1];
	TH1D *h1 = Correl_Resol_T0->ProjectionY(" ",i+1,i+1);
	/*if(i>37)*/
	if( OnlyEB==1 ) h1->Rebin(4);
	if( OnlyEB==2 && i>27 ) h1->Rebin(4);
	if( OnlyEB==2 && i>37 ) h1->Rebin(2);
	if(i<7  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-0.005, 0.02);
	else if(i<11  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-0.01, 0.03);
	else if(i<25  && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-0.015, 0.035);
	else if(i<39 && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-0.05, 0.05);
	else if(i<53 && OnlyEB==2 && XminXmax) h1->GetXaxis()->SetRangeUser(-0.08, 0.08);
	else          h1->GetXaxis()->SetRangeUser(Min, Max);
	TF1 *MyRms = new TF1("MyRms","gaus", Min+0.22, Max-0.2);
	MyRms->SetParName(0,"Normalizzazione"); //MyRms->SetParameters(0,h1->Integral());
	MyRms->SetParName(1,"Mean");            MyRms->SetParameters(1, 0.);
	MyRms->SetParName(2,"Sigma");
	h1->Fit("MyRms");
	T0RMS->SetBinContent( i+1, Makemin( h1->GetRMS(), MyRms->GetParameter(2)) );
	T0RMS->SetBinError( i+1, h1->GetRMSError() );
	stringstream Ind; Ind << i; string Indstr = Ind.str(); h1->Draw(); Hname = OutPutFolder + "/T0RMS_" + TString(Indstr)  + ".png"; gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );
	delete h1;
    }
  }
  //TF1 *Mypol2 = new TF1("Mypol2","[0]+[1]*x",0.0, 0.04);
  TF1 *Mypol2 = new TF1("Mypol2","[0]*x",0.0, 0.04);
  //Mypol2->SetParName(0,"int.");
  Mypol2->SetParName(0,"coeff.");// Mypol1->SetParameter(1, 8.9);
  T0RMS->Fit("Mypol2");
  T0RMS->GetXaxis()->SetTitle("Time Resolution [ns]"); T0RMS->GetYaxis()->SetTitle("T0 Resolution [ns]");
  T0RMS->SetMinimum(0.0); T0RMS->SetMarkerStyle(34); T0RMS->SetMarkerColor(kRed);T0RMS->SetLineColor(kRed); T0RMS->SetFillColor(kRed); T0RMS->Draw("P");
  TLatex lat2; char line2[300];
  lat2.SetNDC(); lat2.SetTextSize(0.040); lat2.SetTextColor(1);
  sprintf(line2,"Int: %.2f [ns]", Mypol2->GetParameter(0) );
  lat2.DrawLatex(xmin,yhi-ypass, line2);
  sprintf(line2,"Coeff: %.2f", Mypol2->GetParameter(1) );
  lat2.DrawLatex(xmin,yhi-2.*ypass, line2 );
  Hname = OutPutFolder + "/Correl_Resol_T0_prfRMS.png";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );

  //Deleting Stuffs
  delete myc1;
  cout<<"----> END of VertexDeterminator... Thanks for choosing VertexDeterminator_minuit.C!"<<endl;
}//VertexDeterminator.C


vector<float> ComputeVertex( float Gtime1, float Gtime2, float sigma_T1, float sigma_T2, TVector3 posExp1, TVector3 posExp2, bool T0Free, float VtxDet_time, float Vx, float Vy){

  // initialize minuit
  TMinuit aMinuit(12);
  int ierflg;
  aMinuit.SetFCN(chi2Hist);
  aMinuit.mnparm(0, "T1Meas", Gtime1,      0., Gtime1,      Gtime1,      ierflg); 
  aMinuit.mnparm(1, "PosX1",  posExp1.X(), 0., posExp1.X(), posExp1.X(), ierflg); 
  aMinuit.mnparm(2, "PosY1",  posExp1.Y(), 0., posExp1.Y(), posExp1.Y(), ierflg); 
  aMinuit.mnparm(3, "PosZ1",  posExp1.Z(), 0., posExp1.Z(), posExp1.Z(), ierflg); 
  aMinuit.mnparm(4, "T1Unc",  sigma_T1,    0., sigma_T1,    sigma_T1,    ierflg); 
  aMinuit.mnparm(5, "T2Meas", Gtime2,      0., Gtime2,      Gtime2,      ierflg); 
  aMinuit.mnparm(6, "PosX2",  posExp2.X(), 0., posExp2.X(), posExp2.X(), ierflg); 
  aMinuit.mnparm(7, "PosY2",  posExp2.Y(), 0., posExp2.Y(), posExp2.Y(), ierflg); 
  aMinuit.mnparm(8, "PosZ2",  posExp2.Z(), 0., posExp2.Z(), posExp2.Z(), ierflg); 
  aMinuit.mnparm(9, "T2Unc",  sigma_T2,    0., sigma_T2,    sigma_T2,    ierflg); 
  aMinuit.mnparm(10, "Vz",    0.,          0.050, -15.0, 15.0, ierflg); 
  if(T0Free){
    aMinuit.mnparm(11, "T0",    0.0,       0.001, -0.3, 0.3, ierflg); 
  }
  else {
    aMinuit.mnparm(11, "T0",  VtxDet_time, 0., VtxDet_time, VtxDet_time, ierflg); 
    aMinuit.FixParameter(11);
  }
  aMinuit.mnparm(12, "Vx",    Vx,          0., Vx, Vx, ierflg);
  aMinuit.mnparm(13, "Vy",    Vy,          0., Vy, Vy, ierflg);
  aMinuit.FixParameter(0); aMinuit.FixParameter(1); aMinuit.FixParameter(2); aMinuit.FixParameter(3); aMinuit.FixParameter(4); aMinuit.FixParameter(5);
  aMinuit.FixParameter(6); aMinuit.FixParameter(7); aMinuit.FixParameter(8); aMinuit.FixParameter(9); aMinuit.FixParameter(2); aMinuit.FixParameter(13);
  Double_t arglis[10];
  arglis[0] = 2000000; // Max calls
  arglis[1] = 0.00001; // Tolerance
  aMinuit.SetPrintLevel(-1);
  aMinuit.mnexcm("MIGRAD", arglis, 2, ierflg);
  aMinuit.mnexcm("MINOS", arglis, 2, ierflg);
  //aMinuit.mnexcm("HESSE", arglis, 2, ierflg);

  //Get Param
  double VZ_best, VZ_best_err;
  double T0_best, T0_best_err;
  aMinuit.GetParameter(10, VZ_best, VZ_best_err);
  aMinuit.GetParameter(11, T0_best, T0_best_err);

  //Write resuts
  vector<float> results;
  results.push_back( VZ_best );
  results.push_back( VZ_best_err );
  results.push_back( T0_best );
  results.push_back( T0_best_err );
  return results;

}//ComputeVertex

float Makemin( float a, float b){
  if(a<=b) return a;
  else     return b;
}

float Makemax( float a, float b){
  if(a>=b) return a;
  else     return b;
}
