# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
#~ from SuSyAachen.DiLeptonHistograms.LeptonFullSimScaleFactorMap_cfi import LeptonFullSimScaleFactorMap as LeptonFullSimScaleFactorMapPars
from SuSyAachen.DiLeptonHistograms.btagEffMap_cfi import bTagEffMap as bTagEffMapPars
from SuSyAachen.DiLeptonHistograms.BTagCalibration_cfi import BTagCalibration as BTagCalibrationPars
from SuSyAachen.DiLeptonHistograms.BTagCalibrationReader_cfi import BTagCalibrationReader as BTagCalibrationReaderPars
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import isolationDefinitions
from SuSyAachen.DiLeptonHistograms.triggerDefinitionMiniAOD_cff import defaultTriggerDefinition as triggerDefinitions
DiLeptonTreesFromMiniAODNoTausMC = cms.EDAnalyzer("DiLeptonTreesFromMiniAODMC",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("qualityJets"),          
   genJets = cms.InputTag("slimmedGenJets"),          
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"), 
   #metNoCleaning = cms.InputTag("slimmedMETs","","Analysis"),    
   #met = cms.InputTag("slimmedMETs","","Analysis"),   
   met = cms.InputTag("slimmedMETs"),                  
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),   
   LHEInfo = cms.InputTag("externalLHEProducer"),                       
   rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),    
   idMapSource = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
   storeMetFilters = cms.untracked.bool(True),
   badPFMuonFilter = cms.InputTag("badPFMuonFilter"),
   badChargedCandidateFilter = cms.InputTag("badChargedCandidateFilter"),
   cloneGlobalMuonFilter = cms.InputTag("cloneGlobalMuonFilter"),
   badGlobalMuonFilter = cms.InputTag("badGlobalMuonFilter"),     
   susyVars = cms.VPSet(),
   pdfWeightTags = cms.VInputTag(),
   bTagEfficiencies = bTagEffMapPars,
   BTagCalibration = BTagCalibrationPars,
   BTagCalibrationReader = BTagCalibrationReaderPars,
   #~ LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars,
   vertexWeights = vertexWeightPars,
   vertexWeightsUp = vertexWeightParsUp,
   vertexWeightsDown = vertexWeightParsDown,                   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
   triggerDefinitions = triggerDefinitions,
   writeID = cms.untracked.bool(False),
   writeTrigger = cms.untracked.bool(False),
   doMETUncert = cms.untracked.bool(False), 
   triggerNames=cms.untracked.vstring(                                     #  1e34  7e33
                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",   
                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",   
                  
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",           #  1     1
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",            #  1     1
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",           
                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",            
                  "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",          
                  "HLT_Mu27_TkMu8_v",                             #  50       1
                  "HLT_Mu30_TkMu11_v",                         #  1     1
                  
                  "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",  #  0     1
                  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",   #  1     1
                  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",   
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",  #  1     1
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  
                  "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",      #  0     1           
                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",      #  1     1        
                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",      #  1     1        
                  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     #  1     1        
                  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",     #  1     1        
                  "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",             #  1     1
                  "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v",             #  1     1
                  
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
   ), 
   metFilterNames=cms.untracked.vstring(                                      
                  "Flag_HBHENoiseFilter", 
                  "Flag_HBHENoiseIsoFilter",
                  "Flag_globalTightHalo2016Filter",   
                  "Flag_goodVertices", 
                  "Flag_EcalDeadCellTriggerPrimitiveFilter",   
                  "Flag_eeBadScFilter",                     

   ), 
)






