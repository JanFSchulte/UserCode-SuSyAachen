# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.LeptonFullSimScaleFactorMap_cfi import LeptonFullSimScaleFactorMap as LeptonFullSimScaleFactorMapPars
from SuSyAachen.DiLeptonHistograms.btagEffMap_cfi import bTagEffMap as bTagEffMapPars
from SuSyAachen.DiLeptonHistograms.BTagCalibration_cfi import BTagCalibration as BTagCalibrationPars
from SuSyAachen.DiLeptonHistograms.BTagCalibrationReader_cfi import BTagCalibrationReader as BTagCalibrationReaderPars
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import isolationDefinitions
from SuSyAachen.DiLeptonHistograms.triggerDefinitionMiniAOD_cff import defaultTriggerDefinition as triggerDefinitions
DiLeptonTreesFromMiniAODNoTaus = cms.EDAnalyzer("DiLeptonTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
   jets = cms.InputTag("qualityJets"),          
   genJets = cms.InputTag("slimmedGenJets"),          
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"),    
   met = cms.InputTag("slimmedMETsModifiedMET", "", "Analysis"),               
   met_normal = cms.InputTag("slimmedMETs"),               
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("isolatedTracks"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),   
   LHEInfo = cms.InputTag("externalLHEProducer"),                       
   rho = cms.InputTag("fixedGridRhoFastjetAll"),    
   idMapSource = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
   storeMetFilters = cms.untracked.bool(True),
   pdfWeightTags = cms.VInputTag(),
   bTagEfficiencies = bTagEffMapPars,
   BTagCalibration = BTagCalibrationPars,
   BTagCalibrationReader = BTagCalibrationReaderPars,
   LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars,
   vertexWeights = vertexWeightPars,
   vertexWeightsUp = vertexWeightParsUp,
   vertexWeightsDown = vertexWeightParsDown,                   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
   triggerDefinitions = triggerDefinitions,
   writeID = cms.untracked.bool(False),
   writeTrigger = cms.untracked.bool(True),
   doMETUncert = cms.untracked.bool(False),
   
   eeTriggerNames=cms.untracked.vstring( 
                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                  "HLT_DoubleEle33_CaloIdL_MW_v",
                  "HLT_DoubleEle25_CaloIdL_MW_v",
   ),
   
   emTriggerNames=cms.untracked.vstring(
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",       
                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",       
                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",            
                  "HLT_Mu27_Ele37_CaloIdL_MW_v",       
                  "HLT_Mu37_Ele27_CaloIdL_MW_v",  
   ),
   
   mmTriggerNames=cms.untracked.vstring(                
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", # from Run2017C
                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", # only Run2017B
                  "HLT_Mu37_TkMu27_v",        
   ),
   
   htTriggerNames=cms.untracked.vstring(
                  "HLT_PFHT180_v",
                  "HLT_PFHT250_v",
                  "HLT_PFHT370_v",
                  "HLT_PFHT430_v",
                  "HLT_PFHT510_v", 
                  "HLT_PFHT590_v",
                  "HLT_PFHT680_v",
                  "HLT_PFHT780_v",
                  "HLT_PFHT890_v",
                  "HLT_PFHT1050_v" 
   ),        
          
   metTriggerNames=cms.untracked.vstring(
                  "HLT_PFMET120_PFMHT120_IDTight_v",
                  "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
                  
   ),
   metFilterNames=cms.untracked.vstring(                                      
                  "Flag_goodVertices",             
                  "Flag_globalSuperTightHalo2016Filter",  
                  "Flag_HBHENoiseFilter",
                  "Flag_HBHENoiseIsoFilter",
                  "Flag_EcalDeadCellTriggerPrimitiveFilter",
                  "Flag_BadPFMuonFilter",
                  "Flag_BadChargedCandidateFilter",
                  "Flag_eeBadScFilter",
                  "Flag_ecalBadCalibFilter"         

   ), 
)
