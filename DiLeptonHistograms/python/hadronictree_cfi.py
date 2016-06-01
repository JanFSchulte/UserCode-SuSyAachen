import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars

mctreewriter = cms.EDAnalyzer('HadronicTree',
    genSrc = cms.InputTag("genParticles"),
    genJets = cms.InputTag("ak4genJetsNoNu"),
    jets = cms.InputTag("triggerMatchedPatJetsPF"),
    met = cms.InputTag("patMETsPF"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    pfCandidates = cms.InputTag("particleFlow"),
    triggers = cms.vstring("HLT_Jet15U"),
    vertexWeights = vertexWeightPars
)
