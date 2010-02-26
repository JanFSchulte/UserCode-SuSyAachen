import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/n/nmohr/scratch0/CMSSW_3_4_1/src/SUSYBSMAnalysis/SusyCAF/test/SUSYPAT.root'
    )
)

process.TFileService = cms.Service('TFileService', fileName = cms.string('MuonTnPOutput.root'))

process.demo = cms.EDAnalyzer('MuonTnPTreeWriter',
mcInfoAvailable = cms.untracked.bool(True),

mcSource = cms.InputTag("genParticles"),
tagSource = cms.InputTag("cleanLayer1Muons"),
probeSource = cms.InputTag("generalTracks"),
passProbeSource = cms.InputTag("cleanLayer1Muons"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),

)


process.p = cms.Path(process.demo)
