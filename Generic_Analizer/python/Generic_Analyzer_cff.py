import FWCore.ParameterSet.Config as cms

GenericAnalizer  = cms.EDAnalyzer("Generic_Analizer",
    candName     = cms.untracked.InputTag('particleFlow'),
    maxPt        = cms.untracked.double(20.),
    SimVtx       = cms.untracked.InputTag('g4SimHits','',''),
    Gamma        = cms.untracked.InputTag('photons','','RECO'),
    RecoVtx      = cms.untracked.InputTag('offlinePrimaryVertices','','RECO'),
    Jet          = cms.untracked.InputTag('ak5PFJets','','RECO'),
    GenJet       = cms.untracked.InputTag('ak5GenJets','',''),
    GenPar       = cms.untracked.InputTag('genParticles','',''),
    JetCHS       = cms.untracked.InputTag('ak5PFJetsCHS','','RECO'),
    ak5PFRho     = cms.untracked.InputTag('ak5PFJets','rho','RECO'),
    isNOPUHgg    = cms.untracked.bool(False),
    WannaFitT0Vtx= cms.untracked.bool(False),
    OneGamma_    = cms.untracked.bool(False),
    Use_MinEne   = cms.untracked.bool(True),
    Use_R9       = cms.untracked.bool(False),
    DoSumEt      = cms.untracked.bool(True),
    DoMass       = cms.untracked.bool(True),
    OutName      = cms.untracked.string("Outfile.root")
)
