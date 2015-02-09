#ifndef FastTool_h
#define FastTool_h

#include <vector>
#include <iostream>
#include "TFile.h"
#include <memory>
#include "TRandom3.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include <cstring>
#include <sstream>
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FastTiming/Generic_Analizer/interface/FastTool.h"
#define PI 3.14159
#define LIGHT_SPEED 29.9792458 //[cm/ns]

class FastTool{

  public:
    FastTool(){ std::cout<<"Hi, I'm FastTool and I will help you in your Analysis."<<std::endl; }
    ~FastTool(){ std::cout<<"End of FastTool... See you soon."<<std::endl; }
    
    void Inizialization( edm::Handle<edm::SimVertexContainer> SimVtx );
    float GiveT0(){   return T0_;   };
    float GiveVtxX(){ return VtxX_; };
    float GiveVtxY(){ return VtxY_; };
    float GiveVtxZ(){ return VtxZ_; };
  private:
    float T0_;
    float VtxX_;
    float VtxY_;
    float VtxZ_;
};

float computeTOF( GlobalPoint Vtx, EBDetId IdXtal, const CaloGeometry* geometry, float EB_LAYER  ){
  const CaloCellGeometry* cell=geometry->getGeometry(IdXtal);
  GlobalPoint Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EB_LAYER );
  GlobalPoint Dist( Xtal.x()-Vtx.x(),  Xtal.y()-Vtx.y(),  Xtal.z()-Vtx.z() );
  float TOF = sqrt( Dist.x()*Dist.x() + Dist.y()*Dist.y() + Dist.z()*Dist.z() ) / LIGHT_SPEED;
  return TOF;
}
float computeTOF( GlobalPoint Vtx, EKDetId IdXtal, const CaloGeometry* geometry, float EE_LAYER ){
  const CaloCellGeometry* cell=geometry->getGeometry(IdXtal);
  GlobalPoint Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( EE_LAYER );
  GlobalPoint Dist( Xtal.x()-Vtx.x(),  Xtal.y()-Vtx.y(),  Xtal.z()-Vtx.z() );
  float TOF = sqrt( Dist.x()*Dist.x() + Dist.y()*Dist.y() + Dist.z()*Dist.z() ) / LIGHT_SPEED;
  return TOF;
}

float computeTOF( GlobalPoint Vtx, GlobalPoint Xtal ){
  GlobalPoint Dist( Xtal.x()-Vtx.x(),  Xtal.y()-Vtx.y(),  Xtal.z()-Vtx.z() );
  float TOF = sqrt( Dist.x()*Dist.x() + Dist.y()*Dist.y() + Dist.z()*Dist.z() ) / LIGHT_SPEED;
  return TOF;
}

float computeTOF( GlobalPoint Vtx, TLorentzVector Xtal ){
  GlobalPoint Dist( Xtal.X()-Vtx.x(),  Xtal.Y()-Vtx.y(),  Xtal.Z()-Vtx.z() );
  float TOF = sqrt( Dist.x()*Dist.x() + Dist.y()*Dist.y() + Dist.z()*Dist.z() ) / LIGHT_SPEED;
  return TOF;
}

float SmearTime(float time, float smearing, TRandom *ranGaus){
  float smearedTime = ranGaus->Gaus(time, smearing);
  return smearedTime;
}

float DeltaR( GlobalPoint a, GlobalPoint b ){
  float deltaEta = a.eta() - b.eta();
  float deltaPhi = a.phi() - b.phi();
  if (deltaPhi > PI)  deltaPhi -= 2.*PI;
  if (deltaPhi < -PI) deltaPhi += 2.*PI;
  return( sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi  ) );
}

double DeltaR( double eta1, double eta2, double phi1, double phi2 ){
  float deltaPhi = eta1 - eta2;
  float deltaEta = phi1 - phi2;
  if (deltaPhi > PI)  deltaPhi -= 2.*PI;
  if (deltaPhi < -PI) deltaPhi += 2.*PI;
  return sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
}

float Delta_phi(float phi1, float phi2) {
  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}
float Delta_eta(float eta1, float eta2) {
  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}
#endif
