import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition


signalNominatorTrees = cms.EDAnalyzer("signalNominatorTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
   jets = cms.InputTag("qualityJets"),	   	   
   genJets = cms.InputTag("slimmedGenJets"),		   	   
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"),	  	 	  
   met = cms.InputTag("slimmedMETs"),  	 	  
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),	
   LHEInfo = cms.InputTag("externalLHEProducer"),		         	      
   rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),	
)

