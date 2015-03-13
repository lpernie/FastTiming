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
process.GenericAnalizer.DoSumEt       = cms.untracked.bool(False)
process.GenericAnalizer.DoMass        = cms.untracked.bool(False)
process.GenericAnalizer.WannaFitT0Vtx = cms.untracked.bool(False)
process.GenericAnalizer.isHgg         = cms.untracked.bool(False)
process.GenericAnalizer.OutName = cms.untracked.string("QCD_noPU.root")
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.ak5PFJets_cfi')

#process.AK5PFJets.src         = cms.InputTag('timeClean', 'timeCleanedPFCandidates')
#process.AK5caPFJetsPruned.src = cms.InputTag('timeClean', 'timeCleanedPFCandidates')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
                            fileNames  = cms.untracked.vstring(
   #'file:/afs/cern.ch/user/s/spigazzi/work/public/xLucaP/RAW-RECO_noPU.root'
   #"root://xrootd-cms.infn.it//store/user/lpernie/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step1_140PU_V1/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step3_140PU_V1/7d29be1aeee7600c1e4046dde551424b/step3_985_1_R7l.root",
   #"root://xrootd-cms.infn.it//store/user/lpernie/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step1_140PU_V1/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step3_140PU_V1/7d29be1aeee7600c1e4046dde551424b/step3_986_1_yRl.root",
   #"root://xrootd-cms.infn.it//store/user/lpernie/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step1_140PU_V1/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step3_140PU_V1/7d29be1aeee7600c1e4046dde551424b/step3_992_1_A6n.root",
   #"root://xrootd-cms.infn.it//store/user/lpernie/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step1_140PU_V1/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step3_140PU_V1/7d29be1aeee7600c1e4046dde551424b/step3_994_1_CUV.root",
   #"root://xrootd-cms.infn.it//store/user/lpernie/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step1_140PU_V1/SingleGamma_Pt35_FastTiming_Shashlik_PU140_TuneZ2starpythia6_cff_py_GEN_step3_140PU_V1/7d29be1aeee7600c1e4046dde551424b/step3_998_4_Uvv.root"
   #"root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_RECO_PU140inTime/QCD_14TeV_\upgradePLS3_RECO_140PUinTime_20.root"
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_1.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_10.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_100.root",  
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_101.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_102.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_103.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_104.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_105.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_106.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_107.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_108.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_109.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_11.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_110.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_111.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_112.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_113.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_114.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_115.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_116.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_117.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_118.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_119.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_12.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_120.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_121.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_122.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_123.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_124.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_125.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_126.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_127.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_128.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_129.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_13.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_130.root", 
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_131.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_132.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_133.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_134.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_135.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_136.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_137.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_138.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_139.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_14.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_140.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_141.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_142.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_143.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_144.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_145.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_146.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_147.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_148.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_149.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_15.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_150.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_151.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_152.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_153.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_154.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_155.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_156.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_157.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_158.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_159.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_16.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_160.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_161.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_162.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_163.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_164.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_165.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_166.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_167.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_168.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_169.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_17.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_170.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_171.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_172.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_173.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_174.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_175.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_176.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_177.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_178.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_179.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_18.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_180.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_181.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_182.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_183.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_184.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_185.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_186.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_187.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_188.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_189.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_19.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_190.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_191.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_192.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_193.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_194.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_195.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_196.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_197.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_198.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_199.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_2.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_20.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_200.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_21.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_22.root", 
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_23.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_24.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_25.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_26.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_28.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_29.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_3.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_30.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_31.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_32.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_33.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_34.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_35.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_36.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_37.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_38.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_39.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_4.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_40.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_41.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_42.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_43.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_44.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_45.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_46.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_47.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_48.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_49.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_5.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_50.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_51.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_52.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_53.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_54.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_55.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_56.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_57.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_58.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_59.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_6.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_60.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_61.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_62.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_63.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_64.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_65.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_66.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_67.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_68.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_69.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_7.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_70.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_71.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_72.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_73.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_74.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_75.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_76.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_77.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_78.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_79.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_8.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_80.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_81.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_82.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_83.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_84.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_85.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_86.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_87.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_88.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_89.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_9.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_90.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_91.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_92.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_93.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_94.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_95.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_96.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_97.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_98.root",
   "root://eoscms//eos/cms/store/caf/user/amartell/fastTiming/QCD_Pt15_500_14TeV_Tune4C_RECO_PU140inTime/QCD_14TeV_upgradePLS3_RECO_140PUinTime_99.root"
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
