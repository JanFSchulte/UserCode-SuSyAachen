import FWCore.ParameterSet.Config as cms

DiLeptonAnalysis = cms.EDAnalyzer('DiLeptonHistograms',

debug = cms.untracked.bool(False),
mcInfo = cms.untracked.bool(True),
tauIsPrompt = cms.untracked.bool(True),
treeInfo = cms.untracked.bool(False),
effInfo = cms.untracked.bool(False),

mcSource = cms.InputTag("genParticles"),
beamSpotSource = cms.InputTag("offlineBeamSpot"),
primaryVertexSource = cms.InputTag("offlinePrimaryVertices"),
muonSource = cms.InputTag("cleanLayer1Muons"),
electronSource = cms.InputTag("cleanLayer1Electrons"),
tauSource = cms.InputTag("cleanLayer1Taus"),
triggerSource = cms.InputTag("TriggerResults","","HLT"),
metSource = cms.InputTag("layer1METsAK5"),
jetSource = cms.InputTag("cleanLayer1JetsAK5"),

CSA_weighted = cms.untracked.bool(False),

acc_MuonPt = cms.untracked.double(5.), 
acc_MuonEta = cms.untracked.double(2.), 

acc_ElectronPt = cms.untracked.double(5.), 
acc_ElectronEta = cms.untracked.double(2.) ,

user_bJetAlgo = cms.untracked.string("trackCountingHighPurBJetTags"),
user_bTagDiscriminator = cms.untracked.double(3.),
#to be removed
trackSource = cms.InputTag("generalTracks"),
#jetMcSource = cms.InputTag(""),

maxJetsForAlphaT = cms.uint32(10)

)
