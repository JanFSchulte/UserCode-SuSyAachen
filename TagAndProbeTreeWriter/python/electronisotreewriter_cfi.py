import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

electronisotreewriter = cms.EDAnalyzer('ElectronIsoTreeWriter',
                                       
src = cms.InputTag("cleanLayer1Electrons"),
jets = cms.InputTag("triggerMatchedPatJetsPF"),
met = cms.InputTag("patMETsPF"),
vertices = cms.InputTag("offlinePrimaryVertices"),

useTauExtensions = cms.bool(False),
useMcInfo = cms.bool(True),
vertexWeights = vertexWeightPars,
isolationDefinitions = isolationDefinitions, 
)

electronisotreewriterWithSecondLepton = electronisotreewriter.clone(
   secondLeptonElectronSrc = cms.InputTag("selectedPatElectrons"),
   secondLeptonMuonSrc = cms.InputTag("selectedPatMuons"),
   secondLeptonTauSrc = cms.InputTag("selectedPatTaus"),
   secondLeptonMinDeltaR = cms.double(0.01)
)
