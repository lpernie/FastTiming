import FWCore.ParameterSet.Config as cms

process = cms.Process('Analyzer')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaper_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 50

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

#Generic Analyzer
process.load("FastTiming.Generic_Analizer.Generic_Analyzer_cff")
process.GenericAnalizer.DoSumEt       = cms.untracked.bool(True)
process.GenericAnalizer.DoMass        = cms.untracked.bool(False)
process.GenericAnalizer.WannaFitT0Vtx = cms.untracked.bool(False)
process.GenericAnalizer.isHgg         = cms.untracked.bool(True)
process.GenericAnalizer.OutName = cms.untracked.string("Hgg_noPU.root")
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.ak5PFJets_cfi')

#process.AK5PFJets.src         = cms.InputTag('timeClean', 'timeCleanedPFCandidates')
#process.AK5caPFJetsPruned.src = cms.InputTag('timeClean', 'timeCleanedPFCandidates')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames  = cms.untracked.vstring(
   #'file:root://eoscms//eos/cms/store/caf/user/lpernie/Upgrade/SHASHLIK_Time/SingleGammaPt35_SHASHTIME_noPU_step3.root'
   'file:Production/Hgg_0PU/step3.root'
   #'file:/afs/cern.ch/work/l/lpernie/ECALpro/gitHubUpgrade/Shashlik/CMSSW_6_2_0_SLHC16/src/SingeGamma_0PU_Shash/Arabella_Twiki/step3.root'
                            )
)

process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

import os
cmssw_base = os.environ['CMSSW_BASE']

process.OurPath = cms.Sequence(
                               process.GenericAnalizer
                              )
process.p = cms.Path(process.OurPath)
process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  outputCommands = cms.untracked.vstring('drop *', 'keep recoPFCandidates_*_*_*'),                                                                                                                      
                                  fileName       = cms.untracked.string ("BBB")                                                                                                                    
                                 )
# process.outpath  = cms.EndPath(process.output)      
