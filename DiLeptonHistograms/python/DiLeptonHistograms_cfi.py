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
acc_MuonEta = cms.untracked.double(2.5), 

acc_ElectronPt = cms.untracked.double(5.), 
acc_ElectronEta = cms.untracked.double(2.5) ,

user_bJetAlgo = cms.untracked.string("trackCountingHighEffBJetTags"),
user_bTagDiscriminator = cms.untracked.double(1.7),
#to be removed
trackSource = cms.InputTag("generalTracks"),
#jetMcSource = cms.InputTag(""),

maxJetsForAlphaT = cms.uint32(10),
fakeRates =  cms.PSet(
        electrons = cms.VPSet(
            cms.PSet(
                etaMin = cms.double(-2.4), etaMax = cms.double(2.4),
                ptMin = cms.double(0), ptMax = cms.double(10999999),
                
                weight = cms.double(0.5)
            ),
            cms.PSet(
                etaMin = cms.double(-1.3),  etaMax = cms.double(1.3),
                ptMin = cms.double(0), ptMax = cms.double(10999999),
    
                weight = cms.double(0.2)
            ),
        ),

        muons =  cms.VPSet(
            cms.PSet(
                etaMin = cms.double(-1.3),  etaMax = cms.double(1.3),
                ptMin = cms.double(0), ptMax = cms.double(10999999),
    
                weight = cms.double(0.2)
            ),                   
        ),
        taus = cms.VPSet(
                                     cms.PSet(
                etaMin = cms.double(-1.3),  etaMax = cms.double(1.3),
                ptMin = cms.double(0), ptMax = cms.double(10999999),
    
                weight = cms.double(0.2)
            ),
        )
),
#fakeRates = cms.bool(False),

)
