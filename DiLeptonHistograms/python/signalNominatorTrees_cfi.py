import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition


signalNominatorTrees = cms.EDAnalyzer("signalNominatorTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   jets = cms.InputTag("qualityJets"),	   	 	  
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),	
   #~ LHEInfo = cms.InputTag("externalLHEProducer"),		         	      
   LHEInfo = cms.InputTag("source"),		         	      
   rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"), 	   	   
   pdgIdDefinition = defaultPdgIdDefinition,	
)

