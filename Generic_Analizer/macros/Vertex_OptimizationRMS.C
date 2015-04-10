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

#define NSTEPMAX 10000
#define PI 3.14159
using namespace std;

float Makemin( float a, float b);
float Makemax( float a, float b);

//.x Vertex_OptimizationRMS.C+( "VTX_Output_EBEE_T0Free_noVBF", "OutFile.root", true, 0 )
void Vertex_OptimizationRMS( TString OutPutFolder="VTX_Output_EBEE_T0Free_noVBF", TString NameFile="OutFile.root", bool T0Free=true, int OnlyEB=0 ){

  //Initial Stuffs
  setTDRStyle();
  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);

  //Get Input files
  TString filename = OutPutFolder + '/' + NameFile;
  TFile* File = TFile::Open( filename.Data() );
  if( !File ) cout<<"WARNING: File "<<NameFile.Data()<<" do not exist"<<endl;
  myc1->cd(); gStyle->SetPalette(1);

  float QuantMin = T0Free ? 10. : 2.5, QuantMax = T0Free ? 90. : 97.5;
  TString Hname;
  bool XminXmax = false;
  bool XminXmax0 = false;
  //Reso VTX
  TH2F *Correl_Resol_Vtx = (TH2F*) File->Get("Correl_Resol_Vtx");
  TH1D *VtxRMS = new TH1D("VtxRMS", "", 100, 0.001, 0.05);
  int index(0);
  for(int i=0; i<Correl_Resol_Vtx->GetNbinsX(); i++){
    if( Correl_Resol_Vtx->ProjectionY(" ",i+1,i+1)->Integral()>0. ){
	double quant[2], value[2];
	TH1D *h1 = Correl_Resol_Vtx->ProjectionY(" ",i+1,i+1);
	value[0]=(h1->Integral()*QuantMax/100)/h1->Integral();
	value[1]=(h1->Integral()*QuantMin/100)/h1->Integral();
	h1->GetQuantiles(2, quant, value);
	float Max = quant[0], Min = quant[1];
	h1->GetXaxis()->SetRangeUser(Min, Max);
	if( i>25 && T0Free ) h1->Rebin(2);
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
	h1->Fit("MyRms");
	index++;
	//VtxRMS->SetBinContent( i+1, Makemin( h1->GetRMS(), MyRms->GetParameter(2)) );
	VtxRMS->SetBinContent( i+1,MyRms->GetParameter(2) );
	VtxRMS->SetBinError( i+1, h1->GetRMSError() );
	stringstream Ind; Ind << i; string Indst = Ind.str(); h1->Draw(); Hname = OutPutFolder + "/VTXRMS_" + TString(Indst)  + ".png"; gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );
	delete h1;
    }
  }
  TF1 *Mypol1 = new TF1("Mypol1","[0]*x",0.0, 0.05);
  Mypol1->SetParName(0,"coeff.");
  VtxRMS->Fit("Mypol1");
  VtxRMS->GetXaxis()->SetTitle("Time Resolution [ns]"); VtxRMS->GetYaxis()->SetTitle("Vertex Resolution [cm]");
  VtxRMS->SetMinimum(0.); VtxRMS->SetMarkerStyle(34); VtxRMS->SetMarkerColor(kRed);VtxRMS->SetLineColor(kRed); VtxRMS->SetFillColor(kRed); VtxRMS->Draw("P");
  TLatex lat; char line[300];
  lat.SetNDC(); lat.SetTextSize(0.040); lat.SetTextColor(1);
  float xmin(0.45), yhi(0.80), ypass(0.05);
  sprintf(line,"Int: %.2f [cm]", Mypol1->GetParameter(0) );
  lat.DrawLatex(xmin,yhi-ypass, line);
  sprintf(line,"Coeff: %.4f", Mypol1->GetParameter(1) );
  lat.DrawLatex(xmin,yhi-2.*ypass, line );
  Hname = OutPutFolder + "/Correl_Resol_Vtx_prfRMS.png";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );
  Hname = OutPutFolder + "/Correl_Resol_Vtx_prfRMS.C";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );
  //RECO T0
  TH2F *Correl_Resol_T0 = (TH2F*) File->Get("Correl_Resol_T0");
  TH1D *T0RMS = new TH1D("T0RMS", "", 100, 0.001, 0.05);
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
	//T0RMS->SetBinContent( i+1, Makemin( h1->GetRMS(), MyRms->GetParameter(2)) );
	T0RMS->SetBinContent( i+1, MyRms->GetParameter(2) );
	T0RMS->SetBinError( i+1, h1->GetRMSError() );
	stringstream Ind; Ind << i; string Indstr = Ind.str(); h1->Draw(); Hname = OutPutFolder + "/T0RMS_" + TString(Indstr)  + ".png"; gStyle->SetOptStat(1111); myc1->SaveAs( Hname.Data() );
	delete h1;
    }
  }
  TF1 *Mypol2 = new TF1("Mypol2","[0]*x",0.0, 0.05);
  Mypol2->SetParName(0,"coeff.");// Mypol1->SetParameter(1, 8.9);
  T0RMS->Fit("Mypol2");
  T0RMS->GetXaxis()->SetTitle("Time Resolution [ns]"); T0RMS->GetYaxis()->SetTitle("T0 Resolution [ns]");
  T0RMS->SetMinimum(0.0); T0RMS->SetMarkerStyle(34); T0RMS->SetMarkerColor(kRed);T0RMS->SetLineColor(kRed); T0RMS->SetFillColor(kRed); T0RMS->Draw("P");
  TLatex lat2; char line2[300];
  lat2.SetNDC(); lat2.SetTextSize(0.040); lat2.SetTextColor(1);
  sprintf(line2,"Int: %.2f [ns]", Mypol2->GetParameter(0) );
  lat2.DrawLatex(xmin,yhi-ypass, line2);
  sprintf(line2,"Coeff: %.4f", Mypol2->GetParameter(1) );
  lat2.DrawLatex(xmin,yhi-2.*ypass, line2 );
  Hname = OutPutFolder + "/Correl_Resol_T0_prfRMS.png";
  gStyle->SetOptStat(0); myc1->SaveAs( Hname.Data() );

  //Deleting Stuffs
  delete myc1;
  cout<<"----> END of VertexDeterminator... Thanks for choosing VertexDeterminator_minuit.C!"<<endl;
}//VertexDeterminator.C

float Makemin( float a, float b){
  if(a<=b) return a;
  else     return b;
}

float Makemax( float a, float b){
  if(a>=b) return a;
  else     return b;
}
