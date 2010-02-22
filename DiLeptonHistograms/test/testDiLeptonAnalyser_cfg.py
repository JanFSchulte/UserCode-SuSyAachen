from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.setName_("TEST")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring() )
process.source.fileNames.append(
'file:/tmp/edelhoff/LM1_7TeV_342/LM1_7TeV_342_file_0.root'
#'file:/tmp/edelhoff/SUSYPAT.root'
)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents.input = 100

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
    fileName = cms.string( "diLeptonTest.root" ),
    closeFileFast = cms.untracked.bool(True) ) 
)

from SuSyAachen.Skimming.defaults.GenLeptons_cff import GenLeptons
GenLeptons(process)

from SuSyAachen.DiLeptonHistograms.MCDiLeptonHistograms_cfi import MCDiLeptonHistograms
process.promptMCHistograms = MCDiLeptonHistograms.clone()

from SuSyAachen.DiLeptonHistograms.DiLeptonHistograms_cfi import DiLeptonAnalysis
process.myDiLeptonAnalysis = DiLeptonAnalysis.clone()

process.p = cms.Path(process.seqGenLeptons + process.promptMCHistograms + process.myDiLeptonAnalysis)

#process.out.outputCommands = cms.untracked.vstring('drop *')
#process.out.splitLevel = cms.untracked.int32(99)  
#process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
#process.out.dropMetaData = cms.untracked.string('DROPPED')
#process.out.fileName = 'edmFile.root'

#process.dump = cms.EDAnalyzer('EventContentAnalyzer') 
     

