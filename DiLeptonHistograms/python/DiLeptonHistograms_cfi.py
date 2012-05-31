import FWCore.ParameterSet.Config as cms

from SuSyAachen.DiLeptonHistograms.fakeRates.electronCenter_cff import fakes as electronCenterFakeRates
from SuSyAachen.DiLeptonHistograms.fakeRates.muonCenter_cff import fakes as muonCenterFakeRates
from SuSyAachen.DiLeptonHistograms.fakeRates.tauCenter_cff import fakes as tauCenterFakeRates
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

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
muonLooseSource = cms.InputTag("cleanLayer1Muons"),
electronLooseSource = cms.InputTag("cleanLayer1Electrons"),
tauLooseSource = cms.InputTag("cleanLayer1Taus"),
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
fakeRates = cms.bool(False),
fakeRatesMc = cms.bool(False),
vertexWeights = vertexWeightPars,
isolationDefinitions = isolationDefinitions,

)

DiLeptonAnalysisInclFake = DiLeptonAnalysis.clone(
    fakeRates =  cms.PSet(
        electrons = electronCenterFakeRates,
        muons =  muonCenterFakeRates,
        taus = tauCenterFakeRates,
    ),
    fakeRatesMc =  cms.PSet(
        electrons = electronCenterFakeRates,
        muons =  muonCenterFakeRates,
        taus = tauCenterFakeRates,
    ),                                          
)
