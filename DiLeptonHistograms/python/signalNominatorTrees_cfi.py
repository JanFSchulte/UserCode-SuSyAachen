import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition


signalNominatorTrees = cms.EDAnalyzer("signalNominatorTrees",
   genParticles = cms.InputTag("prunedGenParticles"),
   #~ LHEInfo = cms.InputTag("externalLHEProducer"),   	   	   	   
   pdgIdDefinition = defaultPdgIdDefinition,   
)

