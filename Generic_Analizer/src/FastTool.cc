#include "TFile.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FastTiming/Generic_Analizer/interface/FastTool.h"

void FastTool::Inizialization( edm::Handle<edm::SimVertexContainer> SimVtx  ){

  edm::SimVertexContainer::const_iterator simVtxFirst = SimVtx->begin();
  GlobalPoint Vtx_sim( simVtxFirst->position().x(), simVtxFirst->position().y(), simVtxFirst->position().z() );
  float T0_Vtx_MC = simVtxFirst->position().t() * pow(10,9);
 
  T0_   = T0_Vtx_MC;
  VtxX_ = 0.;
  VtxY_ = 0.;
  VtxZ_ = 0.; 
}
