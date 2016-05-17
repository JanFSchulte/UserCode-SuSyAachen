# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
MinimalDiLeptonTreesFromMiniAOD = cms.EDAnalyzer("MinimalDiLeptonTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   jets = cms.InputTag("qualityJets"),	   	      	   
   bJets = cms.InputTag("qualityBJets"),  	
   met = cms.InputTag("slimmedMETs","","Analysis"),      	    
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),	
   LHEInfo = cms.InputTag("externalLHEProducer"),
   vertexWeights = vertexWeightPars,  	   	   	   
   pdgIdDefinition = defaultPdgIdDefinition,         
                          
)
