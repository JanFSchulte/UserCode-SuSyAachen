from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.setName_("TEST")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring() )
#process.source.fileNames.append('file:/user/edelhoff/mcData/CMSSW_336patch1/diLepton/LM1AndMimic_v1.baseCuts.SUSY_LM1_FastSim_PAT.EDM.root')
process.source.fileNames.append('file:/user/edelhoff/mcData/CMSSW_336patch1/diLepton/LM1AndMimic_v1.baseCuts.SUSY_LM1_mimic1_FastSim_PAT.EDM.root')
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents.input = -1

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
    fileName = cms.string( "counterTest.root" ),
    closeFileFast = cms.untracked.bool(True) ) 
)

process.load('SuSyAachen.Skimming.diLeptonFilter_cff')
electronSrc = "baseCutsIsoElectrons"
muonSrc = "baseCutsIsoMuons"
tauSrc = "baseCutsIsoTaus"

process.SSEE.primarySrc = electronSrc
process.SSEE.secondarySrc = electronSrc

process.SSEMu.primarySrc = electronSrc
process.SSEMu.secondarySrc = muonSrc

process.SSMuMu.primarySrc = muonSrc
process.SSMuMu.secondarySrc = muonSrc

process.SSETau.primarySrc = electronSrc
process.SSETau.secondarySrc = tauSrc

process.SSMuTau.primarySrc = muonSrc
process.SSMuTau.secondarySrc = tauSrc

process.SSTauTau.primarySrc = tauSrc
process.SSTauTau.secondarySrc = tauSrc

from SuSyAachen.Histograms.triggerResultsCounter_cfi import triggerResultsCounter, makeFilterPaths
process.myCounter = triggerResultsCounter.clone()

process.p = cms.Path( process.seqSSDiLeptons)


process.out.outputCommands = cms.untracked.vstring('keep *')
process.out.splitLevel = cms.untracked.int32(99)  
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')
process.out.fileName = 'edmFile.root'

process.dump = cms.EDAnalyzer('EventContentAnalyzer') 

process.ende = cms.EndPath(process.out* process.myCounter)
makeFilterPaths(process)
     

