import FWCore.ParameterSet.Config as cms

DiLeptonAnalysis = cms.EDFilter('DiLeptonHistograms',

debug = cms.untracked.bool(False),
mcInfo = cms.untracked.bool(True),

mcSource = cms.InputTag("genParticles"),
backmapSource = cms.InputTag("genParticles"),
muonSource = cms.InputTag("selectedLayer1Muons"),
electronSource = cms.InputTag("selectedLayer1Electrons"),
tauSource = cms.InputTag("selectedLayer1Taus"),
metSource = cms.InputTag("selectedLayer1METs"),
jetSource = cms.InputTag("selectedLayer1Jets"),

CSA_weighted = cms.untracked.bool(False),

acc_MuonPt = cms.untracked.double(10.), 
acc_MuonEta = cms.untracked.double(2.), 

acc_ElectronPt = cms.untracked.double(10.), 
acc_ElectronEta = cms.untracked.double(2.) ,

#to be removed
user_bJetAlgo = cms.untracked.string("jetProbabilityBJetTags"),
user_bTagDiscriminator = cms.untracked.double(0.4),
beamSpotSource = cms.InputTag("offlineBeamSpot"),
trackSource = cms.InputTag("generalTracks")
#jetMcSource = cms.InputTag(""),

)
