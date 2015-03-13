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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
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
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#define PI 3.14159
#define LIGHT_SPEED 29.9792458 //[cm/ns]

using namespace std;
using namespace edm;
using namespace reco;

class FastTool{

  public:
    FastTool(){ std::cout<<"Hi, I'm FastTool and I will help you in your Analysis."<<std::endl; }
    ~FastTool(){ std::cout<<"End of FastTool... See you soon."<<std::endl; }
    
    void Inizialization( edm::Handle<edm::SimVertexContainer> SimVtx, float EB_LAYER, float EE_LAYER );
    float GiveT0(){     return T0_;   };
    float GiveVtxX(){   return VtxX_; };
    float GiveVtxY(){   return VtxY_; };
    float GiveVtxZ(){   return VtxZ_; };
    float GetEBLayer(){ return EB_LAYER_; };
    float GetEELayer(){ return EE_LAYER_; };
  private:
    float EB_LAYER_;
    float EE_LAYER_;
    float T0_;
    float VtxX_;
    float VtxY_;
    float VtxZ_;
};
//[0]T [1]X [2]Y [3]Z
vector<float> GetTimeFromJet( const reco::PFJet *PFJets, edm::Handle<edm::SortedCollection<EcalRecHit> >& theBarrelEcalRecHits, edm::Handle<edm::SortedCollection<EcalRecHit> >& theEndcapEcalRecHits, const CaloGeometry* geometry, FastTool *FTool_)
{
  vector<float> TXYZ; TXYZ.clear();
  float finalTime=-999., finalX=-999., finalY=-999., finalZ=-999., minE=0.;
  std::vector <reco::PFCandidatePtr> PFCand = PFJets->getPFConstituents();
  for(unsigned int i=0; i<PFCand.size(); i++){
    PFCandidate MY_cand(PFCand[i]);
    std::vector<EcalRecHit> V_seeds; V_seeds.clear(); std::vector<reco::PFClusterRef> V_cluster; V_cluster.clear(); vector<GlobalPoint> V_seedPos; V_seedPos.clear();
    const EcalRecHit* seedHit = 0;
    GlobalPoint Xtal;
    double maxClusterEnergy = 0.;
    for (auto& blockPair : MY_cand.elementsInBlocks()){
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	{
	  if (blockElement.type() != 4)
	    continue;
	  reco::PFClusterRef cluster = blockElement.clusterRef();
	  if (cluster.isAvailable())
	  {
	    if (cluster->energy() > maxClusterEnergy)
		maxClusterEnergy = cluster->energy();
	    DetId seedID = cluster->seed();
	    if (cluster->layer() == PFLayer::ECAL_BARREL){
		for (auto& rhEB : *theBarrelEcalRecHits){
		  if (rhEB.id() == seedID){
		    seedHit = &rhEB;
		    EBDetId idEB(rhEB.id());
		    const CaloCellGeometry* cell=geometry->getGeometry(idEB);
		    Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( FTool_->GetEBLayer() );
		  }
		}
	    }
	    else if (cluster->layer() == PFLayer::ECAL_ENDCAP){
		for (auto& rhEE : *theEndcapEcalRecHits){
		  if (rhEE.id() == seedID){
		    seedHit = &rhEE;
		    EKDetId idEE(rhEE.id());
		    const CaloCellGeometry* cell=geometry->getGeometry(idEE);
		    Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( FTool_->GetEELayer() );
		  }
		}
	    }
	    if(seedHit){ V_seeds.push_back( *seedHit ); V_cluster.push_back( cluster ); V_seedPos.push_back( Xtal ); }
	  }
	}
    }//All Blocks
    //Now Select the bigger cluster into the SC
    float BestTime = -1, BestEne = -1, Emin=0;
    float BestX = -999., BestY = -999., BestZ = -999.;
    for( int nClu=0; nClu<int(V_cluster.size()); nClu++ ){
	if( V_cluster[nClu]->energy() > Emin ){
	  Emin = V_cluster[nClu]->energy(); //Take seed from most energetic cluster
	  BestTime = V_seeds[nClu].time(); BestEne = V_seeds[nClu].energy();
	  BestX = V_seedPos[nClu].x(); BestY = V_seedPos[nClu].y(); BestZ = V_seedPos[nClu].z();
	}
    }
    //cout<<" I took: "<<BestTime<<" clu size: "<<int(V_cluster.size())<<endl;
    if( BestTime==-1 ) continue;
    if(BestEne>minE){
	//cout<<"  And to PF I giveit !: "<<BestTime<<endl;
	minE = BestEne;
	finalTime = BestTime;
	finalX = BestX;
	finalY = BestY;
	finalZ = BestZ;
    }
  }//End PFCand Loop
  //cout<<"   final time "<<finalTime<<endl;
 
  TXYZ.push_back(finalTime); TXYZ.push_back(finalX); TXYZ.push_back(finalY); TXYZ.push_back(finalZ);
  return TXYZ;
}


float GetTimeFromJet( const reco::PFJet *PFJets, edm::Handle<edm::SortedCollection<EcalRecHit> >& theBarrelEcalRecHits, edm::Handle<edm::SortedCollection<EcalRecHit> >& theEndcapEcalRecHits  )
{
  float finalTime=-999., minE=0.;
  std::vector <reco::PFCandidatePtr> PFCand = PFJets->getPFConstituents();
  for(unsigned int i=0; i<PFCand.size(); i++){
    PFCandidate MY_cand(PFCand[i]);
    std::vector<EcalRecHit> V_seeds; V_seeds.clear(); std::vector<reco::PFClusterRef> V_cluster; V_cluster.clear();
    const EcalRecHit* seedHit = 0;
    double maxClusterEnergy = 0.;
    for (auto& blockPair : MY_cand.elementsInBlocks()){
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	{
	  if (blockElement.type() != 4)
	    continue;
	  reco::PFClusterRef cluster = blockElement.clusterRef();
	  if (cluster.isAvailable())
	  {
	    if (cluster->energy() > maxClusterEnergy)
		maxClusterEnergy = cluster->energy();
	    DetId seedID = cluster->seed();
	    if (cluster->layer() == PFLayer::ECAL_BARREL){
		for (auto& rhEB : *theBarrelEcalRecHits){
		  if (rhEB.id() == seedID){
		    seedHit = &rhEB;
		  }
		}
	    }
	    else if (cluster->layer() == PFLayer::ECAL_ENDCAP){
		for (auto& rhEE : *theEndcapEcalRecHits){
		  if (rhEE.id() == seedID){
		    seedHit = &rhEE;
		  }
		}
	    }
	    if(seedHit){ V_seeds.push_back( *seedHit ); V_cluster.push_back( cluster ); }
	  }
	}
    }//All Blocks
    //Now Select the bigger cluster into the SC
    float BestTime = -1, BestEne = -1, Emin=0;
    for( int nClu=0; nClu<int(V_cluster.size()); nClu++ ){
	if( V_cluster[nClu]->energy() > Emin ){
	  Emin = V_cluster[nClu]->energy(); //Take seed from most energetic cluster
	  BestTime = V_seeds[nClu].time(); BestEne = V_seeds[nClu].energy();
	  //cout<<"Looking into the clusters of this PFCand: time: "<<V_seeds[nClu].time()<<" ene seed: "<<V_seeds[nClu].energy()<<" eene clu: "<<V_cluster[nClu]->energy()<<endl;
	}
    }
    //cout<<" I took: "<<BestTime<<" clu size: "<<int(V_cluster.size())<<endl;
    if( BestTime==-1 ) continue;
    if(BestEne>minE){
	//cout<<"  And to PF I giveit !: "<<BestTime<<endl;
	minE = BestEne;
	finalTime = BestTime;
    }
  }//End PFCand Loop
  //cout<<"   final time "<<finalTime<<endl;
  return finalTime;
}

vector<float> GetTimeFromGamma( const reco::Photon* photon, edm::Handle<edm::SortedCollection<EcalRecHit> >& theBarrelEcalRecHits, edm::Handle<edm::SortedCollection<EcalRecHit> >& theEndcapEcalRecHits, const CaloGeometry* geometry, FastTool *FTool_ )
{
  vector<float> T_XYZ; T_XYZ.clear();
  float posX=-99., posY=-99., posZ=-99.;
  float ThisTime=-99.;
  const EcalRecHit* seedHit = 0;
  DetId seedId = photon->superCluster()->seed()->seed();
  GlobalPoint Xtal;
  if( photon->isEB() ){
    for (auto& rhEB : *theBarrelEcalRecHits){
	if (rhEB.id() == seedId ){
	  seedHit = &rhEB;
	  EBDetId idEB(rhEB.id());
	  const CaloCellGeometry* cell=geometry->getGeometry(idEB);
	  Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( FTool_->GetEBLayer() );
	}
    }
  }
  else{
    for (auto& rhEE : *theEndcapEcalRecHits){
	if (rhEE.id() == seedId ){
	  seedHit = &rhEE;
	  EKDetId idEE(rhEE.id());
	  const CaloCellGeometry* cell=geometry->getGeometry(idEE);
	  Xtal = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( FTool_->GetEELayer() );
	}
    }
  }
  if(seedHit){
    ThisTime = seedHit->time();
    posX = Xtal.x(); posY = Xtal.y(); posZ = Xtal.z();
  }
  T_XYZ.push_back(ThisTime); T_XYZ.push_back(posX); T_XYZ.push_back(posY); T_XYZ.push_back(posZ);
  return T_XYZ;
}


float GetTimeFromGamma( const reco::Photon* photon, edm::Handle<edm::SortedCollection<EcalRecHit> >& theBarrelEcalRecHits, edm::Handle<edm::SortedCollection<EcalRecHit> >& theEndcapEcalRecHits  )
{
  float ThisTime=-99;
  const EcalRecHit* seedHit = 0;
  DetId seedId = photon->superCluster()->seed()->seed();
  if( photon->isEB() ){
    for (auto& rhEB : *theBarrelEcalRecHits){
	if (rhEB.id() == seedId ){
	  seedHit = &rhEB;
	}
    }
  }
  else{
    for (auto& rhEE : *theEndcapEcalRecHits){
	if (rhEE.id() == seedId ){
	  seedHit = &rhEE;
	}
    }
  }
  if(seedHit){
    ThisTime = seedHit->time();
  }
  return ThisTime;
}

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

float DeltaR( GlobalPoint a, TLorentzVector b ){
  float deltaEta = a.eta() - b.Eta();
  float deltaPhi = a.phi() - b.Phi();
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
