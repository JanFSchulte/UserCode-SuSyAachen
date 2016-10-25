# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions
from SuSyAachen.DiLeptonHistograms.triggerDefinitionMiniAOD_cff import defaultTriggerDefinition as triggerDefinitions
DiLeptonTreesFromMiniAODNoTaus = cms.EDAnalyzer("DiLeptonTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("qualityJets"),	   	   
   genJets = cms.InputTag("slimmedGenJets"),	   	   
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"),	
   #~ met = cms.InputTag("slimmedMETs","","Analysis"),  	
   met = cms.InputTag("slimmedMETs"),  	 	     	    
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),	
   LHEInfo = cms.InputTag("externalLHEProducer"),		         	      
   rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),	   
   susyVars = cms.VPSet(),
   pdfWeightTags = cms.VInputTag(),
   vertexWeights = vertexWeightPars,
   vertexWeightsUp = vertexWeightParsUp,
   vertexWeightsDown = vertexWeightParsDown,   	   	   	   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
   triggerDefinitions = triggerDefinitions,
   writeID = cms.untracked.bool(False),
   writeTrigger = cms.untracked.bool(True),
   triggerNames=cms.vstring(													#	1e34	7e33
						"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",			# 	0 		1
						"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",			# 	1 		1
						"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",				    # 	1 		1
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",				# 	1 		1
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",				# 	1 		1
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",				
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",				
						"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",				
						"HLT_Mu27_TkMu8_v",										#	50 		1
						"HLT_Mu30_TkMu11_v",									#	1 		1
						"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",	# 	0 		1
						"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",	# 	1 		1
						"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",	
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",	# 	1 		1
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",	
						"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",		# 	0 		1				
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",		# 	1 		1			
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",		# 	1 		1			
						"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",		# 	1 		1			
						"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",		# 	1 		1			
						"HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",					# 	1 		1
						
						"HLT_PFHT125_v",
						"HLT_PFHT200_v",
						"HLT_PFHT250_v",
						"HLT_PFHT300_v",
						"HLT_PFHT350_v", 
						"HLT_PFHT400_v",
						"HLT_PFHT475_v",
						"HLT_PFHT600_v",
						"HLT_PFHT650_v",
						"HLT_PFHT800_v",
						"HLT_PFHT900_v",
						"HLT_MET200_v",
						"HLT_MET250_v",
						"HLT_MET300_v",
						"HLT_MET600_v",
						"HLT_MET700_v",
						
						"HLT_PFMET300_v",
						"HLT_PFMET400_v",
						"HLT_PFMET500_v",
						"HLT_PFMET600_v",
						"HLT_PFMET170_NotCleaned_v",
						"HLT_PFMET170_HBHECleaned_v",
						"HLT_PFMET170_BeamHaloCleaned_v",
						"HLT_PFMET170_HBHE_BeamHaloCleaned_v",
						"HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned_v",
						"HLT_PFMET170_JetIdCleaned_v",
						"HLT_PFMET170_NoiseCleaned_v",
						       
						"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_v",		# 	1 		1		Emergency 0
						"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v",		# 	1 		1
						"HLT_DoubleMu8_Mass8_PFHT250_v",						# 	1 		1		Emergency 0
						"HLT_DoubleMu8_Mass8_PFHT300_v",						# 	1 		1
						"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_v",		# 	1 		1		Emergency 0
						"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v",		# 	1 		1
	), 
)






