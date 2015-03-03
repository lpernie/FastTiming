#include <stdlib.h>
#include <string>
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>

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
#include "TLatex.h"

#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
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

using namespace RooFit;

//.x SumEtPlot.C+("Efficiencies","effi")
//.x SumEtPlot.C+("Efficiencies","effi","Hgg")
//.x SumEtPlot.C+("HiggsMass","MassH2")
//.x SumEtPlot.C+("SumEtPlots","SumEtPlots")
void SumEtPlot( TString OutPutFolder="plots", TString mode="SumEtPlots", TString mode2="" ){

  //Initial Stuffs
  std::cout<<"----> START SumEtPlot.C"<<std::endl;
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
  using namespace std;
  gStyle->SetOptStat(0);
  //Get Input files
  if(OutPutFolder != ""){
    if(mode=="effi" && mode2!="Hgg") OutPutFolder = OutPutFolder + "_QCD";
    if(mode=="effi" && mode2=="Hgg") OutPutFolder = OutPutFolder + "_Hgg";
    TString Comm = "mkdir -p " + OutPutFolder;
    system( Comm.Data() );
  }
  if( mode == "SumEtPlots" ){
    //Open Files
    TFile* File140 = TFile::Open( "../../../Hgg_140PU.root" );
    if( !File140 ) cout<<"WARNING: File GenericAnalyzer_Higgs_140.root do not exist"<<endl;
    TFile* File0 = TFile::Open( "../../../Hgg_noPU.root" );
    if( !File0 ) cout<<"WARNING: File GenericAnalyzer_Higgs_0.root do not exist"<<endl;
    myc1->cd();
    gStyle->SetPalette(1);
    //Open Histos
    TH1F *h140_SumEt        = (TH1F*) File140->Get( "h_SumEt" );
    TH1F *h140_SumEt_cut    = (TH1F*) File140->Get( "h_SumEt_cut" );
    TH1F *h140_SumEt_cut_15 = (TH1F*) File140->Get( "h_SumEt_15cut" );
    TH1F *h140_SumEt_cut_30 = (TH1F*) File140->Get( "h_SumEt_30cut" );
    TH1F *h140_SumEt_cut_50 = (TH1F*) File140->Get( "h_SumEt_50cut" );
    TH1F *h140_SumEt_cut_500= (TH1F*) File140->Get( "h_SumEt_500cut" );
    TH1F *h0_SumEt          = (TH1F*) File0->Get( "h_SumEt" );
    TH1F *h0_SumEt_cut      = (TH1F*) File0->Get( "h_SumEt_cut" );
    TH1F *h0_SumEt_cut_15   = (TH1F*) File0->Get( "h_SumEt_15cut" );
    TH1F *h0_SumEt_cut_30   = (TH1F*) File0->Get( "h_SumEt_30cut" );
    TH1F *h0_SumEt_cut_50   = (TH1F*) File0->Get( "h_SumEt_50cut" );
    TH1F *h0_SumEt_cut_500  = (TH1F*) File0->Get( "h_SumEt_500cut" );
    h0_SumEt_cut->GetXaxis()->SetTitle("#Sigma Et [GeV]"); h0_SumEt_cut->SetTitle("CMS Simulation Preliminary");
    //Rebin and Normalize
    h140_SumEt->Rebin(8); h140_SumEt_cut->Rebin(8); h0_SumEt->Rebin(8); h0_SumEt_cut->Rebin(8);
    h140_SumEt->Scale( 1./h140_SumEt->Integral() );  h140_SumEt_cut->Scale( 1./h140_SumEt_cut->Integral() ); h0_SumEt->Scale( 1./h0_SumEt->Integral() );  h0_SumEt_cut->Scale(1./h0_SumEt_cut->Integral() );
    //Draw
    h0_SumEt_cut->SetLineColor(kBlue); h0_SumEt_cut->Draw(); h0_SumEt->SetLineColor(kRed); h0_SumEt->Draw("same");
    h140_SumEt_cut->SetLineColor(kGreen); h140_SumEt_cut->Draw("same"); h140_SumEt->SetLineColor(kBlack); h140_SumEt->Draw("same");
    TLegend* leg = new TLegend(.38, 0.7, .86, .85);
    leg->SetFillColor(0);
    leg->SetTextSize(0.031);
    leg->AddEntry(h0_SumEt, "Total #Sigma Et (no PU)", "L");
    leg->AddEntry(h0_SumEt_cut, "Total #Sigma Et, time cut (no PU)", "L");
    leg->AddEntry(h140_SumEt, "Total #Sigma Et (140PU)", "L");
    leg->AddEntry(h140_SumEt_cut, "Total #Sigma Et, time cut (140PU)", "L");
    leg->Draw("same");
    TLatex *tex = new TLatex(0.18,0.87,"VBF H #rightarrow #gamma #gamma sample");
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
    TString name = OutPutFolder + "/SumEt.png";
    myc1->SaveAs( name.Data() );
    name = OutPutFolder + "/SumEt.C";
    myc1->SaveAs( name.Data() );
    //Rebin and Normalize SMEARING 15
    h140_SumEt_cut_15->Rebin(8); h0_SumEt_cut_15->Rebin(8);
    h140_SumEt_cut_15->Scale( 1./h140_SumEt_cut_15->Integral() ); h0_SumEt_cut_15->Scale( 1./h0_SumEt_cut_15->Integral() );
    //Draw
    h0_SumEt_cut_15->SetLineColor(kBlue); h0_SumEt_cut_15->Draw(); h0_SumEt->SetLineColor(kRed); h0_SumEt->Draw("same");
    h140_SumEt_cut_15->SetLineColor(kGreen); h140_SumEt_cut_15->Draw("same"); h140_SumEt->SetLineColor(kBlack); h140_SumEt->Draw("same");
    TLegend* leg2a = new TLegend(.38, 0.7, .86, .85);
    leg2a->SetFillColor(0);
    leg2a->SetTextSize(0.031);
    leg2a->AddEntry(h0_SumEt, "Total #Sigma Et (no PU)", "L");
    leg2a->AddEntry(h0_SumEt_cut_15, "Total #Sigma Et, time cut (no PU)", "L");
    leg2a->AddEntry(h140_SumEt, "Total #Sigma Et (140PU)", "L");
    leg2a->AddEntry(h140_SumEt_cut_15, "Total #Sigma Et, time cut (140PU)", "L");
    leg2a->Draw("same");
    TLatex *tex2a = new TLatex(0.18,0.87,"VBF H #rightarrow #gamma #gamma sample (15 ps smear.)");
    tex2a->SetNDC();
    tex2a->SetTextSize(0.04);
    tex2a->SetLineWidth(2);
    tex2a->Draw();
    name = OutPutFolder + "/SumEt_15.png";
    myc1->SaveAs( name.Data() );
    name = OutPutFolder + "/SumEt_15.C";
    myc1->SaveAs( name.Data() );
    //Rebin and Normalize SMEARING 30
    h140_SumEt_cut_30->Rebin(8); h0_SumEt_cut_30->Rebin(8);
    h140_SumEt_cut_30->Scale( 1./h140_SumEt_cut_30->Integral() ); h0_SumEt_cut_30->Scale( 1./h0_SumEt_cut_30->Integral() );
    //Draw
    h0_SumEt_cut_30->SetLineColor(kBlue); h0_SumEt_cut_30->Draw(); h0_SumEt->SetLineColor(kRed); h0_SumEt->Draw("same");
    h140_SumEt_cut_30->SetLineColor(kGreen); h140_SumEt_cut_30->Draw("same"); h140_SumEt->SetLineColor(kBlack); h140_SumEt->Draw("same");
    TLegend* leg2b = new TLegend(.38, 0.7, .86, .85);
    leg2b->SetFillColor(0);
    leg2b->SetTextSize(0.031);
    leg2b->AddEntry(h0_SumEt, "Total #Sigma Et (no PU)", "L");
    leg2b->AddEntry(h0_SumEt_cut_30, "Total #Sigma Et, time cut (no PU)", "L");
    leg2b->AddEntry(h140_SumEt, "Total #Sigma Et (140PU)", "L");
    leg2b->AddEntry(h140_SumEt_cut_30, "Total #Sigma Et, time cut (140PU)", "L");
    leg2b->Draw("same");
    TLatex *tex2b = new TLatex(0.18,0.87,"VBF H #rightarrow #gamma #gamma sample (30 ps smear.)");
    tex2b->SetNDC();
    tex2b->SetTextSize(0.04);
    tex2b->SetLineWidth(2);
    tex2b->Draw();
    name = OutPutFolder + "/SumEt_30.png";
    myc1->SaveAs( name.Data() );
    name = OutPutFolder + "/SumEt_30.C";
    myc1->SaveAs( name.Data() );
    //Rebin and Normalize SMEARING 50
    h140_SumEt_cut_50->Rebin(8); h0_SumEt_cut_50->Rebin(8);
    h140_SumEt_cut_50->Scale( 1./h140_SumEt_cut_50->Integral() ); h0_SumEt_cut_50->Scale( 1./h0_SumEt_cut_50->Integral() );
    //Draw
    h0_SumEt_cut_50->SetLineColor(kBlue); h0_SumEt_cut_50->Draw(); h0_SumEt->SetLineColor(kRed); h0_SumEt->Draw("same");
    h140_SumEt_cut_50->SetLineColor(kGreen); h140_SumEt_cut_50->Draw("same"); h140_SumEt->SetLineColor(kBlack); h140_SumEt->Draw("same");
    TLegend* leg2c = new TLegend(.38, 0.7, .86, .85);
    leg2c->SetFillColor(0);
    leg2c->SetTextSize(0.031);
    leg2c->AddEntry(h0_SumEt, "Total #Sigma Et (no PU)", "L");
    leg2c->AddEntry(h0_SumEt_cut_50, "Total #Sigma Et, time cut (no PU)", "L");
    leg2c->AddEntry(h140_SumEt, "Total #Sigma Et (140PU)", "L");
    leg2c->AddEntry(h140_SumEt_cut_50, "Total #Sigma Et, time cut (140PU)", "L");
    leg2c->Draw("same");
    TLatex *tex2c = new TLatex(0.18,0.87,"VBF H #rightarrow #gamma #gamma sample (50 ps smear.)");
    tex2c->SetNDC();
    tex2c->SetTextSize(0.04);
    tex2c->SetLineWidth(2);
    tex2c->Draw();
    name = OutPutFolder + "/SumEt_50.png";
    myc1->SaveAs( name.Data() );
    name = OutPutFolder + "/SumEt_50.C";
    myc1->SaveAs( name.Data() );
    //Rebin and Normalize SMEARING 500
    h140_SumEt_cut_500->Rebin(8); h0_SumEt_cut_500->Rebin(8);
    h140_SumEt_cut_500->Scale( 1./h140_SumEt_cut_500->Integral() ); h0_SumEt_cut_500->Scale( 1./h0_SumEt_cut_500->Integral() );
    //Draw
    h0_SumEt_cut_500->SetLineColor(kBlue); h0_SumEt_cut_500->Draw(); h0_SumEt->SetLineColor(kRed); h0_SumEt->Draw("same");
    h140_SumEt_cut_500->SetLineColor(kGreen); h140_SumEt_cut_500->Draw("same"); h140_SumEt->SetLineColor(kBlack); h140_SumEt->Draw("same");
    TLegend* leg2d = new TLegend(.38, 0.7, .86, .85);
    leg2d->SetFillColor(0);
    leg2d->SetTextSize(0.031);
    leg2d->AddEntry(h0_SumEt, "Total #Sigma Et (no PU)", "L");
    leg2d->AddEntry(h0_SumEt_cut_500, "Total #Sigma Et, time cut (no PU)", "L");
    leg2d->AddEntry(h140_SumEt, "Total #Sigma Et (140PU)", "L");
    leg2d->AddEntry(h140_SumEt_cut_500, "Total #Sigma Et, time cut (140PU)", "L");
    leg2d->Draw("same");
    TLatex *tex2d = new TLatex(0.18,0.87,"VBF H #rightarrow #gamma #gamma sample (500 ps smear.)");
    tex2d->SetNDC();
    tex2d->SetTextSize(0.04);
    tex2d->SetLineWidth(2);
    tex2d->Draw();
    name = OutPutFolder + "/SumEt_500.png";
    myc1->SaveAs( name.Data() );
    name = OutPutFolder + "/SumEt_500.C";
    myc1->SaveAs( name.Data() );
  }
  else if( mode == "MassH1" ){
    //Open Files
    TFile* File140_noclean = TFile::Open( "OutAnalizerFile_Higgs_NOClean_140.root" );
    if( !File140_noclean ) cout<<"WARNING: OutAnalizerFile_Higgs_NOClean_140.root do not exist"<<endl;
    TFile* File0_noclean = TFile::Open( "OutAnalizerFile_Higgs_NOClean_noPU.root" );
    if( !File0_noclean ) cout<<"WARNING: File OutAnalizerFile_Higgs_NOClean_noPU.roott do not exist"<<endl;
    TFile* File140_clean = TFile::Open( "OutAnalizerFile_Higgs_Clean_140.root" );
    if( !File140_clean ) cout<<"WARNING: File OutAnalizerFile_Higgs_Clean_140.root do not exist"<<endl;
    TFile* File0_clean = TFile::Open( "OutAnalizerFile_Higgs_Clean_noPU.root" );
    if( !File0_clean ) cout<<"WARNING: File OutAnalizerFile_Higgs_Clean_noPU.root do not exist"<<endl;    
    //Open GAMMA Histos
    TH1F *h_HiggsMass140_noclean = (TH1F*) File140_noclean->Get( "h_HiggsMass" ); 
    TH1F *h_HiggsMass140_clean   = (TH1F*) File140_clean->Get( "h_HiggsMass" );
    TH1F *h_HiggsMass0_noclean   = (TH1F*) File0_noclean->Get( "h_HiggsMass" ); 
    TH1F *h_HiggsMass0_clean     = (TH1F*) File0_clean->Get( "h_HiggsMass" );  
    h_HiggsMass140_noclean->Scale( 1./h_HiggsMass140_noclean->Integral() ); h_HiggsMass140_noclean->SetLineColor( kRed );
    h_HiggsMass140_noclean->Draw();
    h_HiggsMass140_clean->Scale( 1./h_HiggsMass140_clean->Integral() ); h_HiggsMass140_clean->SetLineColor( kBlue );
    h_HiggsMass140_clean->Draw("same"); 
    TString name = OutPutFolder + "/HiggsMass_140.png"; myc1->SaveAs( name.Data() );
    h_HiggsMass0_noclean->Scale( 1./h_HiggsMass0_noclean->Integral() ); h_HiggsMass0_noclean->SetLineColor( kRed );
    h_HiggsMass0_noclean->Draw();
    h_HiggsMass0_clean->Scale( 1./h_HiggsMass0_clean->Integral() ); h_HiggsMass0_clean->SetLineColor( kBlue );
    h_HiggsMass0_clean->Draw("same"); 
    name = OutPutFolder + "/HiggsMass_0.png"; myc1->SaveAs( name.Data() );
    //Open PFGAMMA Histos
    TH1F *h_HiggsPFMass140_noclean = (TH1F*) File140_noclean->Get( "h_HiggsPFMass" );
    TH1F *h_HiggsPFMass140_clean   = (TH1F*) File140_clean->Get( "h_HiggsPFMass" );
    TH1F *h_HiggsPFMass0_noclean   = (TH1F*) File0_noclean->Get( "h_HiggsPFMass" );
    TH1F *h_HiggsPFMass0_clean     = (TH1F*) File0_clean->Get( "h_HiggsPFMass" );
    myc1->cd();
    gStyle->SetPalette(1);
    //Fit
    bool UseBkg=true;
    for(int nToFit=0; nToFit<4; nToFit++){
	cout<<"------------------------- FIT: "<<nToFit<<" ------------------------------"<<endl;
	TH1F *histoToFit;
	if(nToFit==0) histoToFit = h_HiggsPFMass140_noclean;
	if(nToFit==1) histoToFit = h_HiggsPFMass140_clean;
	if(nToFit==2) histoToFit = h_HiggsPFMass0_noclean;
	if(nToFit==3) histoToFit = h_HiggsPFMass0_clean;
	histoToFit->Rebin(2);
	RooRealVar x("x","#gamma#gamma invariant mass", 101, 159, "GeV/c^2");
	RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x), histoToFit);
	RooRealVar mass("mass","#gamma#gamma peak position", 127., 110., 128.,"GeV/c^{2}"); mass.setRange( 110., 140. ); mass.setVal( 125. );
	RooRealVar sigma("sigma","#gamma#gamma sigma", 3., 0.01, 5.,""); sigma.setRange( 0.01, 4. ); sigma.setVal( 3. );
	RooRealVar alpha("alpha","#gamma#gamma alpha", 1., 0.01, 20,""); alpha.setRange( 0.01, 20. ); alpha.setVal( 1. );
	RooRealVar enne("enne","#gamma#gamma n", 3, 0, 10,""); enne.setRange( 0., 10. ); enne.setVal( 3. );
	RooCBShape CriBa("CriBa","Crystal Ball",x, mass, sigma, alpha, enne);
	RooRealVar Nsig("Nsig","Total yields", histoToFit->Integral(), histoToFit->Integral()*0.7, histoToFit->Integral()*1.3);
	Nsig.setRange( histoToFit->Integral()*0.2, histoToFit->Integral()*100 ); Nsig.setVal( histoToFit->Integral() );
	RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
	RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
	RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
	RooArgList cbpars(cb0,cb1,cb2);
	RooChebychev bkg("bkg","bkg model", x, cbpars );
	RooRealVar Nbkg("Nbkg","Bkg yealds", histoToFit->Integral()*0.1, 0., histoToFit->Integral()*0.5);
	Nbkg.setRange( 0., histoToFit->Integral()*0.5 ); Nsig.setVal( histoToFit->Integral()*0.1 );
	RooAbsPdf* model=0;
	RooAddPdf model1("model","CB only", RooArgList(CriBa) ,RooArgList(Nsig));
	RooAddPdf model2("model","CB + bkg",RooArgList(CriBa,bkg),RooArgList(Nsig,Nbkg));
	if(!UseBkg) model = &model1;
	else        model = &model2;
	RooNLLVar nll("nll","log likelihood var", *model, dh );
	RooMinuit m(nll);
	m.setVerbose(kFALSE);
	m.migrad();
	//RooFitResult* res = m.save() ;
	//RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
	//int ndof = histoToFit->GetNbinsX() - res->floatParsFinal().getSize();
	RooPlot*  xframe = x.frame( histoToFit->GetNbinsX() );
	if(nToFit==0) xframe->SetTitle( "Higgs Mass 140 PU (Not Cleaned)" );
	if(nToFit==1) xframe->SetTitle( "Higgs Mass 140 PU (Cleaned)" );
	if(nToFit==2) xframe->SetTitle( "Higgs Mass 0 PU (Not Cleaned)" );
	if(nToFit==3) xframe->SetTitle( "Higgs Mass 0 PU (Cleaned)" );
	dh.plotOn(xframe);
	model->plotOn(xframe,Components(CriBa),LineStyle(kDashed), LineColor(kRed));
	model->plotOn(xframe);
	myc1->cd(); xframe->Draw();
	TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030); lat.SetTextColor(1);
	char line[300];
	float xmin(0.55), yhi(0.80), ypass(0.05);
	sprintf(line,"m_{#gamma#gamma}: %.3f #pm %.3f GeV", mass.getVal(), mass.getError() );
	lat.DrawLatex(xmin,yhi, line);
	sprintf(line,"Sigma: %.3f #pm %.3f GeV", sigma.getVal(), sigma.getError() );
	lat.DrawLatex(xmin,yhi-ypass, line);
	if(nToFit==0) name = OutPutFolder + "/FIT_HiggsPFMass_140_NotCleaned.png";
	if(nToFit==1) name = OutPutFolder + "/FIT_HiggsPFMass_140_Cleaned.png";
	if(nToFit==2) name = OutPutFolder + "/FIT_HiggsPFMass_0_NotCleaned.png";
	if(nToFit==3) name = OutPutFolder + "/FIT_HiggsPFMass_0_Cleaned.png";
	myc1->SaveAs( name.Data() );
    }
    //Draw
    h_HiggsPFMass140_clean->Scale( 1./h_HiggsPFMass140_clean->Integral() ); h_HiggsPFMass140_clean->SetLineColor( kBlue );
    h_HiggsPFMass140_clean->Draw();
    h_HiggsPFMass140_noclean->Scale( 1./h_HiggsPFMass140_noclean->Integral() ); h_HiggsPFMass140_noclean->SetLineColor( kRed );
    h_HiggsPFMass140_noclean->Draw("same");
    name = OutPutFolder + "/HiggsPFMass_140.png"; myc1->SaveAs( name.Data() );
    h_HiggsPFMass0_clean->Scale( 1./h_HiggsPFMass0_clean->Integral() ); h_HiggsPFMass0_clean->SetLineColor( kBlue );
    h_HiggsPFMass0_clean->Draw();
    h_HiggsPFMass0_noclean->Scale( 1./h_HiggsPFMass0_noclean->Integral() ); h_HiggsPFMass0_noclean->SetLineColor( kRed );
    h_HiggsPFMass0_noclean->Draw("same");
    name = OutPutFolder + "/HiggsPFMass_0.png"; myc1->SaveAs( name.Data() );
  }//Mass1
  else if( mode == "MassH2" ){
    //Files Opening
    TFile* File140 = TFile::Open( "time_S_Higgs_140PU_Tot.root" );
    if( !File140 ) cout<<"WARNING: GenericAnalyzer_Higgs_140_HiggsMassCleaned.root do not exist"<<endl;
    TFile* File0   = TFile::Open( "time_S_Higgs_noPU_Tot.root" );
    if( !File0 )   cout<<"WARNING: GenericAnalyzer_Higgs_0_HiggsMassCleaned.root do not exist"<<endl;
    TH1F *h_HiggsMass140_noclean = (TH1F*) File140->Get( "h_HiggsPFMass" );
    TH1F *h_HiggsMass140_clean   = (TH1F*) File140->Get( "h_HiggsPFMass_Mit" );
    TH1F *h_HiggsMass0_noclean   = (TH1F*) File0->Get( "h_HiggsPFMass" );
    TH1F *h_HiggsMass0_clean     = (TH1F*) File0->Get( "h_HiggsPFMass_Mit" );   
    TH1F *h_Ereso1_0             = (TH1F*) File0->Get( "Ereso1" );   
    TH1F *h_Ereso1_Mit_0         = (TH1F*) File0->Get( "Ereso1_Mit" );   
    TH1F *h_Ereso1_140           = (TH1F*) File140->Get( "Ereso1" );   
    TH1F *h_Ereso1_Mit_140       = (TH1F*) File140->Get( "Ereso1_Mit" );   
    float eff=1.;// = h_HiggsMass0_clean->GetEntries()/h_HiggsMass0_noclean->GetEntries();
    TString name;
    myc1->cd();
    gStyle->SetPalette(1);
    //Fit
    bool UseBkg=true;
    for(int nToFit=0; nToFit<8; nToFit++){
	cout<<"------------------------- FIT: "<<nToFit<<" ------------------------------"<<endl;
	TH1F *histoToFit;
	if(nToFit==0) histoToFit = h_HiggsMass140_noclean;
	if(nToFit==1) histoToFit = h_HiggsMass140_clean;
	if(nToFit==2) histoToFit = h_HiggsMass0_noclean;
	if(nToFit==3) histoToFit = h_HiggsMass0_clean;
	if(nToFit==4) histoToFit = h_Ereso1_0;
	if(nToFit==5) histoToFit = h_Ereso1_Mit_0;
	if(nToFit==6) histoToFit = h_Ereso1_140;
	if(nToFit==7) histoToFit = h_Ereso1_Mit_140;
	if( nToFit < 4 ) histoToFit->Rebin(2);
	float x_min=101., x_max=159;
	if( nToFit==4 ||  nToFit==5 ){ x_min=0.95; x_max=1.05; }
	if( nToFit==6 ||  nToFit==7 ){ x_min=0.972; x_max=1.1; }
	RooRealVar x("x","#gamma#gamma invariant mass", x_min, x_max, "GeV/c^2");
	RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x), histoToFit);
	float m_init = 125., m_min = 110., m_max = 128.;
	if( nToFit==4 || nToFit==5 ) { m_init = 1.; m_min = 0.96; m_max = 1.04; }
	if( nToFit==6 || nToFit==7 ) { m_init = 1.02; m_min = 0.96; m_max = 1.04; }
	RooRealVar mass("mass","#gamma#gamma peak position", m_init, m_min, m_max,"GeV/c^{2}");
	mass.setRange( nToFit<4 ? 110. : 0.96, nToFit<4 ? 140. : 1.04 ); mass.setVal(  nToFit<4 ? 125. : 1. );
	RooRealVar sigma("sigma","#gamma#gamma sigma", 3., 0.01, 5.,""); sigma.setRange( 0.01, 5. ); sigma.setVal( 3. );
	RooRealVar alpha("alpha","#gamma#gamma alpha", 1., 0.01, 20,""); alpha.setRange( 0.01, 20. ); alpha.setVal( 1. );
	RooRealVar enne("enne","#gamma#gamma n", 3, 0, 10,""); enne.setRange( 0., 10. ); enne.setVal( 3. );
	RooGaussian gaus("gaus","Core Gaussian",x, mass, sigma);
	RooCBShape CriBa("CriBa","Crystal Ball",x, mass, sigma, alpha, enne);
	RooRealVar Nsig("Nsig","Total yields", histoToFit->Integral(), histoToFit->Integral()*0.7, histoToFit->Integral()*1.3);
	Nsig.setRange( histoToFit->Integral()*0.2, histoToFit->Integral()*100 ); Nsig.setVal( histoToFit->Integral() );
	RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
	RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
	RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
	RooArgList cbpars(cb0,cb1,cb2);
	RooChebychev bkg("bkg","bkg model", x, cbpars );
	RooRealVar Nbkg("Nbkg","Bkg yealds", histoToFit->Integral()*0.1, 0., histoToFit->Integral()*0.5);
	Nbkg.setRange( 0., histoToFit->Integral()*0.5 ); Nsig.setVal( histoToFit->Integral()*0.1 );
	RooAbsPdf* model=0;
	RooAddPdf model0("model","Gaus only", RooArgList(gaus,bkg) ,RooArgList(Nsig,Nbkg));
	RooAddPdf model1("model","CB only", RooArgList(CriBa) ,RooArgList(Nsig));
	RooAddPdf model2("model","CB + bkg",RooArgList(CriBa,bkg),RooArgList(Nsig,Nbkg));
	if( nToFit > 3 ) model = &model0;
	else{
	  if(!UseBkg) model = &model1;
	  else        model = &model2;
	}
	RooNLLVar nll("nll","log likelihood var", *model, dh );
	RooMinuit m(nll);
	m.setVerbose(kFALSE);
	m.migrad();
	//RooFitResult* res = m.save() ;
	//RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
	//int ndof = histoToFit->GetNbinsX() - res->floatParsFinal().getSize();
	float Gauss_area =  Nsig.getVal();
	float Bkg_area   =  Nbkg.getVal();
	cout<<"nToFit: "<<Gauss_area<<" "<<Bkg_area<<" Int h "<<histoToFit->Integral("width")<<endl;
	RooPlot*  xframe = x.frame( histoToFit->GetNbinsX() );
	if(nToFit==0) xframe->SetTitle( "Higgs Mass 140 PU (Not Cleaned)" );
	if(nToFit==1) xframe->SetTitle( "Higgs Mass 140 PU (Cleaned)" );
	if(nToFit==2) xframe->SetTitle( "Higgs Mass 0 PU (Not Cleaned)" );
	if(nToFit==3) xframe->SetTitle( "Higgs Mass 0 PU (Cleaned)" );
	if(nToFit==4) xframe->SetTitle( "E Reso 0 PU (Not Cleaned)" );
	if(nToFit==5) xframe->SetTitle( "E Reso 0 PU (Cleaned)" );
	if(nToFit==6) xframe->SetTitle( "E Reso 140 PU (Not Cleaned)" );
	if(nToFit==7) xframe->SetTitle( "E Reso 140 PU (Cleaned)" );
	dh.plotOn(xframe);
	if( nToFit < 4 ){ model->plotOn(xframe,Components(CriBa),LineStyle(kDashed), LineColor(kRed)); }
	else             { model->plotOn(xframe,Components(gaus),LineStyle(kDashed), LineColor(kRed));  }
	model->plotOn(xframe);
	myc1->cd(); xframe->Draw();
	TLatex lat; lat.SetNDC(); lat.SetTextSize(0.030); lat.SetTextColor(1);
	char line[300];
	float xmin(0.55), yhi(0.80), ypass(0.05);
	sprintf(line,"m_{#gamma#gamma}: %.3f #pm %.3f GeV", mass.getVal(), mass.getError() );
	lat.DrawLatex(xmin,yhi, line);
	sprintf(line,"Sigma: %.3f #pm %.3f GeV", sigma.getVal(), sigma.getError() );
	lat.DrawLatex(xmin,yhi-ypass, line);
	//sprintf(line,"Gaus Int.: %.3f", Gauss_area );
	//lat.DrawLatex(xmin,yhi-ypass*2, line);
	//sprintf(line,"Total Int.: %.3f", Gauss_area+Bkg_area );
	//lat.DrawLatex(xmin,yhi-ypass*3, line);
	if(nToFit==0) name = OutPutFolder + "/FIT_HiggsPFMass_140_NotCleaned.png";
	if(nToFit==1) name = OutPutFolder + "/FIT_HiggsPFMass_140_Cleaned.png";
	if(nToFit==2) name = OutPutFolder + "/FIT_HiggsPFMass_0_NotCleaned.png";
	if(nToFit==3) name = OutPutFolder + "/FIT_HiggsPFMass_0_Cleaned.png";
	if(nToFit==4) name = OutPutFolder + "/FIT_Ereso_0_NotCleaned.png";
	if(nToFit==5) name = OutPutFolder + "/FIT_Ereso_0_Cleaned.png";
	if(nToFit==6) name = OutPutFolder + "/FIT_Ereso_140_NotCleaned.png";
	if(nToFit==7) name = OutPutFolder + "/FIT_Ereso_140_Cleaned.png";
	myc1->SaveAs( name.Data() );
	if(nToFit==2) eff /= Nsig.getVal();
	if(nToFit==3) eff *= Nsig.getVal();
    }
    //Draw
    h_HiggsMass140_clean->Scale( 1./h_HiggsMass140_clean->Integral() ); h_HiggsMass140_clean->SetLineColor( kBlue );
    h_HiggsMass140_clean->Draw();
    h_HiggsMass140_noclean->Scale( 1./h_HiggsMass140_noclean->Integral() ); h_HiggsMass140_noclean->SetLineColor( kRed );
    h_HiggsMass140_noclean->Draw("same");
    name = OutPutFolder + "/HiggsPFMass_140.png"; myc1->SaveAs( name.Data() );
    h_HiggsMass0_clean->Scale( 1./h_HiggsMass0_clean->Integral() ); h_HiggsMass0_clean->SetLineColor( kBlue );
    h_HiggsMass0_clean->Draw();
    h_HiggsMass0_noclean->Scale( 1./h_HiggsMass0_noclean->Integral() ); h_HiggsMass0_noclean->SetLineColor( kRed );
    h_HiggsMass0_noclean->Draw("same");
    name = OutPutFolder + "/HiggsPFMass_0.png"; myc1->SaveAs( name.Data() );
    h_Ereso1_Mit_140->Scale( 1./h_Ereso1_Mit_140->Integral() ); h_Ereso1_Mit_140->SetLineColor( kBlue );
    h_Ereso1_Mit_140->Draw();
    h_Ereso1_140->Scale( 1./h_Ereso1_140->Integral() ); h_Ereso1_140->SetLineColor( kRed );
    h_Ereso1_140->Draw("same");
    name = OutPutFolder + "/Ereso_140.png"; myc1->SaveAs( name.Data() );
    h_Ereso1_Mit_0->Scale( 1./h_Ereso1_Mit_0->Integral() ); h_Ereso1_Mit_0->SetLineColor( kBlue );
    h_Ereso1_Mit_0->Draw();
    h_Ereso1_0->Scale( 1./h_Ereso1_0->Integral() ); h_Ereso1_0->SetLineColor( kRed );
    h_Ereso1_0->Draw("same");
    name = OutPutFolder + "/Ereso_0.png"; myc1->SaveAs( name.Data() );
    cout<<"EFFICIENCY on SIGNAL: "<<eff<<endl;
  }//Mass2
  else if( mode == "effi" ){
    //Files Opening
    TString Name140 = "time_S_QCD_140PU.root";
    TString Name70  = "time_S_QCD_70PU.root";
    TString Name0  = "time_S_QCD_noPU.root";
    if( mode2 == "Hgg" ){
	Name140 = "time_S_Higgs_140PU.root";
	Name70  = "time_S_Higgs_140PU.root";
	Name0   = "time_S_Higgs_noPU.root";
    }
    TFile* File140 = TFile::Open( Name140.Data() );
    if( !File140 ) cout<<"WARNING: "<<Name140.Data()<<" do not exist"<<endl;
    TFile* File70  = TFile::Open( Name70.Data() );
    if( !File70 )  cout<<"WARNING: "<<Name70.Data()<<" do not exist"<<endl;
    TFile* File0   = TFile::Open( Name0.Data() );
    if( !File70 )   cout<<"WARNING: "<<Name0.Data()<<" do not exist"<<endl;
    TH1F *h_GoodJet_140        = (TH1F*) File140->Get( "h_GoodJet_t" );      h_GoodJet_140->SetLineColor(kRed);           h_GoodJet_140->Scale(1./h_GoodJet_140->Integral());
    TH1F *h_GoodJet_140EB      = (TH1F*) File140->Get( "h_GoodJet_tEB" );    h_GoodJet_140EB->SetLineColor(kRed);         h_GoodJet_140EB->Scale(1./h_GoodJet_140EB->Integral());
    TH1F *h_GoodJet_140EE      = (TH1F*) File140->Get( "h_GoodJet_tEE" );    h_GoodJet_140EE->SetLineColor(kRed);         h_GoodJet_140EE->Scale(1./h_GoodJet_140EE->Integral());
    TH1F *h_BadJet_140         = (TH1F*) File140->Get( "h_BadJet_t" );       h_BadJet_140->SetLineColor(kRed);            h_BadJet_140->Scale(1./h_BadJet_140->Integral());
    TH1F *h_GoodGamma_140      = (TH1F*) File140->Get( "h_GoodGamma_t" );    h_GoodGamma_140->SetLineColor(kRed);         h_GoodGamma_140->Scale(1./h_GoodGamma_140->Integral());
    TH1F *h_GoodGamma_140EB    = (TH1F*) File140->Get( "h_GoodGamma_tEB" );  h_GoodGamma_140EB->SetLineColor(kRed);       h_GoodGamma_140EB->Scale(1./h_GoodGamma_140EB->Integral());
    TH1F *h_GoodGamma_140EE    = (TH1F*) File140->Get( "h_GoodGamma_tEE" );  h_GoodGamma_140EE->SetLineColor(kRed);       h_GoodGamma_140EE->Scale(1./h_GoodGamma_140EE->Integral());
    TH1F *h_BadGamma_140       = (TH1F*) File140->Get( "h_BadGamma_t" );     h_BadGamma_140->SetLineColor(kRed);          h_BadGamma_140->Scale(1./h_BadGamma_140->Integral());
    TH1F *h_PtGenJet140        = (TH1F*) File140->Get( "h_PtGenJet" );       h_PtGenJet140->SetLineColor(kRed);           h_PtGenJet140->Scale(1./h_PtGenJet140->Integral()); h_PtGenJet140->Rebin(2);
    TH1F *h_EtaGenJet140       = (TH1F*) File140->Get( "h_EtaGenJet" );      h_EtaGenJet140->SetLineColor(kRed);          h_EtaGenJet140->Scale(1./h_EtaGenJet140->Integral()); h_EtaGenJet140->Rebin(2);
    TH1F *h_EffEtaTot_phot_140 = (TH1F*) File140->Get( "h_EffEtaTot_phot" ); h_EffEtaTot_phot_140->SetLineColor(kRed);    
    TH1F *h_EffEtaTot_jet_140  = (TH1F*) File140->Get( "h_EffEtaTot_jet" );  h_EffEtaTot_jet_140->SetLineColor(kRed);
    TH1F *h_EffEta1_phot_140   = (TH1F*) File140->Get( "h_EffEta_phot1" );   h_EffEta1_phot_140->SetLineColor(kRed);    
    TH1F *h_EffEta1_jet_140    = (TH1F*) File140->Get( "h_EffEta_jet1" );    h_EffEta1_jet_140->SetLineColor(kRed);
    TH1F *h_EffEta2_phot_140   = (TH1F*) File140->Get( "h_EffEta_phot2" );   h_EffEta2_phot_140->SetLineColor(kRed);    
    TH1F *h_EffEta2_jet_140    = (TH1F*) File140->Get( "h_EffEta_jet2" );    h_EffEta2_jet_140->SetLineColor(kRed);
    TH1F *h_EffEta3_phot_140   = (TH1F*) File140->Get( "h_EffEta_phot3" );   h_EffEta3_phot_140->SetLineColor(kRed);    
    TH1F *h_EffEta3_jet_140    = (TH1F*) File140->Get( "h_EffEta_jet3" );    h_EffEta3_jet_140->SetLineColor(kRed);
    TH1F *h_GoodJet_70         = (TH1F*) File70->Get( "h_GoodJet_t" );       h_GoodJet_70->SetLineColor(kGreen+1);        h_GoodJet_70->Scale(1./h_GoodJet_70->Integral());
    TH1F *h_GoodJet_70EB       = (TH1F*) File70->Get( "h_GoodJet_tEB" );     h_GoodJet_70EB->SetLineColor(kGreen+1);      h_GoodJet_70EB->Scale(1./h_GoodJet_70EB->Integral());
    TH1F *h_GoodJet_70EE       = (TH1F*) File70->Get( "h_GoodJet_tEE" );     h_GoodJet_70EE->SetLineColor(kGreen+1);      h_GoodJet_70EE->Scale(1./h_GoodJet_70EE->Integral());
    TH1F *h_BadJet_70          = (TH1F*) File70->Get( "h_BadJet_t" );        h_BadJet_70->SetLineColor(kGreen+1);         h_BadJet_70->Scale(1./h_BadJet_70->Integral());
    TH1F *h_GoodGamma_70       = (TH1F*) File70->Get( "h_GoodGamma_t" );     h_GoodGamma_70->SetLineColor(kGreen+1);      h_GoodGamma_70->Scale(1./h_GoodGamma_70->Integral());
    TH1F *h_GoodGamma_70EB     = (TH1F*) File70->Get( "h_GoodGamma_tEB" );   h_GoodGamma_70EB->SetLineColor(kGreen+1);    h_GoodGamma_70EB->Scale(1./h_GoodGamma_70EB->Integral());
    TH1F *h_GoodGamma_70EE     = (TH1F*) File70->Get( "h_GoodGamma_tEE" );   h_GoodGamma_70EE->SetLineColor(kGreen+1);    h_GoodGamma_70EE->Scale(1./h_GoodGamma_70EE->Integral());
    TH1F *h_BadGamma_70        = (TH1F*) File70->Get( "h_BadGamma_t" );      h_BadGamma_70->SetLineColor(kGreen+1);       h_BadGamma_70->Scale(1./h_BadGamma_70->Integral());
    TH1F *h_PtGenJet70         = (TH1F*) File70->Get( "h_PtGenJet" );        h_PtGenJet70->SetLineColor(kGreen+1);        h_PtGenJet70->Scale(1./h_PtGenJet70->Integral()); h_PtGenJet70->Rebin(2);
    TH1F *h_EtaGenJet70        = (TH1F*) File70->Get( "h_EtaGenJet" );       h_EtaGenJet70->SetLineColor(kGreen+1);       h_EtaGenJet70->Scale(1./h_EtaGenJet70->Integral()); h_EtaGenJet70->Rebin(2);
    TH1F *h_EffEtaTot_phot_70  = (TH1F*) File70->Get( "h_EffEtaTot_phot" );  h_EffEtaTot_phot_70->SetLineColor(kGreen+1);
    TH1F *h_EffEtaTot_jet_70   = (TH1F*) File70->Get( "h_EffEtaTot_jet" );   h_EffEtaTot_jet_70->SetLineColor(kGreen+1);
    TH1F *h_EffEta1_phot_70    = (TH1F*) File70->Get( "h_EffEta_phot1" );    h_EffEta1_phot_70->SetLineColor(kGreen+1);    
    TH1F *h_EffEta1_jet_70     = (TH1F*) File70->Get( "h_EffEta_jet1" );     h_EffEta1_jet_70->SetLineColor(kGreen+1);
    TH1F *h_EffEta2_phot_70    = (TH1F*) File70->Get( "h_EffEta_phot2" );    h_EffEta2_phot_70->SetLineColor(kGreen+1);    
    TH1F *h_EffEta2_jet_70     = (TH1F*) File70->Get( "h_EffEta_jet2" );     h_EffEta2_jet_70->SetLineColor(kGreen+1);
    TH1F *h_EffEta3_phot_70    = (TH1F*) File70->Get( "h_EffEta_phot3" );    h_EffEta3_phot_70->SetLineColor(kGreen+1);    
    TH1F *h_EffEta3_jet_70     = (TH1F*) File70->Get( "h_EffEta_jet3" );     h_EffEta3_jet_70->SetLineColor(kGreen+1);
    TH1F *h_GoodJet_0          = (TH1F*) File0->Get( "h_GoodJet_t" );        h_GoodJet_0->SetLineColor(kBlue);            h_GoodJet_0->Scale(1./h_GoodJet_0->Integral());
    TH1F *h_GoodJet_0EB        = (TH1F*) File0->Get( "h_GoodJet_tEB" );      h_GoodJet_0EB->SetLineColor(kBlue);          h_GoodJet_0EB->Scale(1./h_GoodJet_0EB->Integral());
    TH1F *h_GoodJet_0EE        = (TH1F*) File0->Get( "h_GoodJet_tEE" );      h_GoodJet_0EE->SetLineColor(kBlue);          h_GoodJet_0EE->Scale(1./h_GoodJet_0EE->Integral());
    TH1F *h_BadJet_0           = (TH1F*) File0->Get( "h_BadJet_t" );         h_BadJet_0->SetLineColor(kBlue);             h_BadJet_0->Scale(1./h_BadJet_0->Integral());
    TH1F *h_GoodGamma_0        = (TH1F*) File0->Get( "h_GoodGamma_t" );      h_GoodGamma_0->SetLineColor(kBlue);          h_GoodGamma_0->Scale(1./h_GoodGamma_0->Integral());
    TH1F *h_GoodGamma_0EB      = (TH1F*) File0->Get( "h_GoodGamma_tEB" );    h_GoodGamma_0EB->SetLineColor(kBlue);        h_GoodGamma_0EB->Scale(1./h_GoodGamma_0EB->Integral());
    TH1F *h_GoodGamma_0EE      = (TH1F*) File0->Get( "h_GoodGamma_tEE" );    h_GoodGamma_0EE->SetLineColor(kBlue);        h_GoodGamma_0EE->Scale(1./h_GoodGamma_0EE->Integral());
    TH1F *h_BadGamma_0         = (TH1F*) File0->Get( "h_BadGamma_t" );       h_BadGamma_0->SetLineColor(kBlue);           h_BadGamma_0->Scale(1./h_BadGamma_0->Integral());
    TH1F *h_PtGenJet0          = (TH1F*) File0->Get( "h_PtGenJet" );         h_PtGenJet0->SetLineColor(kBlue);            h_PtGenJet0->Scale(1./h_PtGenJet0->Integral()); h_PtGenJet0->Rebin(2);
    TH1F *h_EtaGenJet0         = (TH1F*) File0->Get( "h_EtaGenJet" );        h_EtaGenJet0->SetLineColor(kBlue);           h_EtaGenJet0->Scale(1./h_EtaGenJet0->Integral()); h_EtaGenJet0->Rebin(2);
    TH1F *h_EffEtaTot_phot_0   = (TH1F*) File0->Get( "h_EffEtaTot_phot" );   h_EffEtaTot_phot_0->SetLineColor(kBlue);
    TH1F *h_EffEtaTot_jet_0    = (TH1F*) File0->Get( "h_EffEtaTot_jet" );    h_EffEtaTot_jet_0->SetLineColor(kBlue);
    TH1F *h_EffEta1_phot_0     = (TH1F*) File0->Get( "h_EffEta_phot1" );     h_EffEta1_phot_0->SetLineColor(kBlue);    
    TH1F *h_EffEta1_jet_0      = (TH1F*) File0->Get( "h_EffEta_jet1" );      h_EffEta1_jet_0->SetLineColor(kBlue);
    TH1F *h_EffEta2_phot_0     = (TH1F*) File0->Get( "h_EffEta_phot2" );     h_EffEta2_phot_0->SetLineColor(kBlue);    
    TH1F *h_EffEta2_jet_0      = (TH1F*) File0->Get( "h_EffEta_jet2" );      h_EffEta2_jet_0->SetLineColor(kBlue);
    TH1F *h_EffEta3_phot_0     = (TH1F*) File0->Get( "h_EffEta_phot3" );     h_EffEta3_phot_0->SetLineColor(kBlue);    
    TH1F *h_EffEta3_jet_0      = (TH1F*) File0->Get( "h_EffEta_jet3" );      h_EffEta3_jet_0->SetLineColor(kBlue);
    //Draw
    TString name  = OutPutFolder + "/Time_GoodJet.png"; myc1->cd();
    h_GoodJet_0->Draw(); h_GoodJet_70->Draw("same"); h_GoodJet_140->Draw("same"); myc1->SaveAs( name.Data() );
     name  = OutPutFolder + "/Time_GoodJet_EB.png"; myc1->cd();
    h_GoodJet_0EB->Draw(); h_GoodJet_70EB->Draw("same"); h_GoodJet_140EB->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Time_GoodJet_EE.png"; myc1->cd();
    h_GoodJet_0EE->Draw(); h_GoodJet_70EE->Draw("same"); h_GoodJet_140EE->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Time_GoodGamma.png";
    h_GoodGamma_0->Draw(); h_GoodGamma_70->Draw("same"); h_GoodGamma_140->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Time_GoodGamma_EB.png";
    h_GoodGamma_0EB->Draw(); h_GoodGamma_70EB->Draw("same"); h_GoodGamma_140EB->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Time_GoodGamma_EE.png";
    h_GoodGamma_0EE->Draw(); h_GoodGamma_70EE->Draw("same"); h_GoodGamma_140EE->Draw("same"); myc1->SaveAs( name.Data() );
    if( mode2 != "Hgg" ){
	name  = OutPutFolder + "/PtGenJet.png";
	h_PtGenJet0->Draw(); h_PtGenJet70->Draw("same"); h_PtGenJet140->Draw("same"); myc1->SaveAs( name.Data() );
	name  = OutPutFolder + "/EtaGenJet.png";
	h_EtaGenJet0->Draw(); h_EtaGenJet70->Draw("same"); h_EtaGenJet140->Draw("same"); myc1->SaveAs( name.Data() );
    }
    name  = OutPutFolder + "/Nreco_Jet_vs_eta.png"; h_EffEta2_jet_0->SetMinimum(0);
    h_EffEta2_jet_0->Draw(); h_EffEta2_jet_70->Draw("same"); h_EffEta2_jet_140->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Nreco_Gamma_vs_eta.png"; h_EffEta2_phot_140->SetMinimum(0); h_EffEta2_phot_0->SetMinimum(0);
    if( mode2 != "Hgg" ){ h_EffEta2_phot_140->Draw(); h_EffEta2_phot_70->Draw("same"); h_EffEta2_phot_0->Draw("same"); myc1->SaveAs( name.Data() ); }
    if( mode2 == "Hgg" ){ h_EffEta2_phot_0->Draw(); h_EffEta2_phot_70->Draw("same"); h_EffEta2_phot_140->Draw("same"); myc1->SaveAs( name.Data() ); }
    name  = OutPutFolder + "/Ntot_Gamma_vs_eta.png"; h_EffEtaTot_phot_140->SetMinimum(0); h_EffEtaTot_phot_0->SetMinimum(0);
    if( mode2 != "Hgg" ){ h_EffEtaTot_phot_140->Draw(); h_EffEtaTot_phot_70->Draw("same"); h_EffEtaTot_phot_0->Draw("same"); myc1->SaveAs( name.Data() );}
    if( mode2 == "Hgg" ){ h_EffEtaTot_phot_0->Draw(); h_EffEtaTot_phot_70->Draw("same"); h_EffEtaTot_phot_140->Draw("same"); myc1->SaveAs( name.Data() );}
    name  = OutPutFolder + "/Ntot_Jet_vs_eta.png"; h_EffEtaTot_jet_0->SetMinimum(0);
    h_EffEtaTot_jet_0->Draw(); h_EffEtaTot_jet_70->Draw("same"); h_EffEtaTot_jet_140->Draw("same"); myc1->SaveAs( name.Data() );
    //
    float int0 = h_EffEtaTot_jet_0->Integral(), int70 = h_EffEtaTot_jet_70->Integral(), int140 = h_EffEtaTot_jet_140->Integral();
    name  = OutPutFolder + "/CorrectAsso_Jet.png";
    h_EffEtaTot_jet_0->Scale(1./int0); h_EffEtaTot_jet_70->Scale(1./int70); h_EffEtaTot_jet_140->Scale(1./int140);
    h_EffEtaTot_jet_0->SetMinimum(0); h_EffEtaTot_jet_0->Draw(); h_EffEtaTot_jet_70->Draw("same"); h_EffEtaTot_jet_140->Draw("same"); myc1->SaveAs( name.Data() );
    h_EffEtaTot_jet_0->Scale(int0); h_EffEtaTot_jet_70->Scale(int70); h_EffEtaTot_jet_140->Scale(int140);
    name  = OutPutFolder + "/CorrectAsso_Phot.png";
    int0 = h_EffEtaTot_phot_0->Integral(); int70 = h_EffEtaTot_phot_70->Integral(); int140 = h_EffEtaTot_phot_140->Integral();
    h_EffEtaTot_phot_0->Scale(1./int0); h_EffEtaTot_phot_70->Scale(1./int70); h_EffEtaTot_phot_140->Scale(1./int140);
    h_EffEtaTot_phot_0->SetMinimum(0); h_EffEtaTot_phot_0->Draw(); h_EffEtaTot_phot_70->Draw("same"); h_EffEtaTot_phot_140->Draw("same"); myc1->SaveAs( name.Data() );
    h_EffEtaTot_phot_0->Scale(int0); h_EffEtaTot_phot_70->Scale(int70); h_EffEtaTot_phot_140->Scale(int140);
    //
    for( int i=0; i<h_EffEtaTot_phot_140->GetNbinsX(); i++ ){ h_EffEtaTot_phot_140->SetBinContent( i+1, h_EffEta2_phot_140->GetBinContent(i+1)/h_EffEtaTot_phot_140->GetBinContent(i+1) ); }
    for( int i=0; i<h_EffEtaTot_jet_140->GetNbinsX(); i++ ){  h_EffEtaTot_jet_140->SetBinContent(  i+1, h_EffEta2_jet_140->GetBinContent(i+1) / h_EffEtaTot_jet_140->GetBinContent(i+1) ); }
    for( int i=0; i<h_EffEtaTot_phot_70->GetNbinsX(); i++ ){ h_EffEtaTot_phot_70->SetBinContent( i+1, h_EffEta2_phot_70->GetBinContent(i+1)/h_EffEtaTot_phot_70->GetBinContent(i+1) ); }
    for( int i=0; i<h_EffEtaTot_jet_70->GetNbinsX(); i++ ){  h_EffEtaTot_jet_70->SetBinContent(  i+1, h_EffEta2_jet_70->GetBinContent(i+1) / h_EffEtaTot_jet_70->GetBinContent(i+1) ); }
    for( int i=0; i<h_EffEtaTot_phot_0->GetNbinsX(); i++ ){ h_EffEtaTot_phot_0->SetBinContent( i+1, h_EffEta2_phot_0->GetBinContent(i+1)/h_EffEtaTot_phot_0->GetBinContent(i+1) ); }
    for( int i=0; i<h_EffEtaTot_jet_0->GetNbinsX(); i++ ){  h_EffEtaTot_jet_0->SetBinContent(  i+1, h_EffEta2_jet_0->GetBinContent(i+1) / h_EffEtaTot_jet_0->GetBinContent(i+1) ); }
    name  = OutPutFolder + "/Effi_Gamma.png"; 
    h_EffEtaTot_phot_0->SetMinimum(0.7); h_EffEtaTot_phot_0->SetMaximum(1.02);
    h_EffEtaTot_phot_0->Draw(); h_EffEtaTot_phot_70->Draw("same"); h_EffEtaTot_phot_140->Draw("same"); myc1->SaveAs( name.Data() );
    name  = OutPutFolder + "/Effi_Jet.png";
    h_EffEtaTot_jet_0->SetMinimum(0.); h_EffEtaTot_jet_0->SetMaximum(1.2);
    h_EffEtaTot_jet_0->Draw(); h_EffEtaTot_jet_70->Draw("same"); h_EffEtaTot_jet_140->Draw("same"); myc1->SaveAs( name.Data() );

  }  
  delete myc1;
  cout<<"----> END of SumEtPlot... Thanks for choosing SumEtPlot.C!"<<endl;
}//SumEtPlot.C
