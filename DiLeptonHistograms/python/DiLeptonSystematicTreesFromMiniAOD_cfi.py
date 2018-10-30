# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
#~ from SuSyAachen.DiLeptonHistograms.LeptonFastSimScaleFactorMap_cfi import LeptonFastSimScaleFactorMap as LeptonFastSimScaleFactorMapPars
#~ from SuSyAachen.DiLeptonHistograms.LeptonFullSimScaleFactorMap_cfi import LeptonFullSimScaleFactorMap as LeptonFullSimScaleFactorMapPars
from SuSyAachen.DiLeptonHistograms.btagEffMap_cfi import bTagEffMap as bTagEffMapPars
#from SuSyAachen.DiLeptonHistograms.btagEffMapFastSim_cfi import bTagEffMap as bTagEffMapPars
from SuSyAachen.DiLeptonHistograms.BTagCalibration_cfi import BTagCalibration as BTagCalibrationPars
from SuSyAachen.DiLeptonHistograms.BTagCalibrationReader_cfi import BTagCalibrationReader as BTagCalibrationReaderPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignal_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignalUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignalDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import isolationDefinitions

DiLeptonSystematicTreesFromMiniAODNoTaus = cms.EDAnalyzer("DiLeptonSystematicTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("qualityJets"),          
   genJets = cms.InputTag("slimmedGenJets"),             
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"), 
   met = cms.InputTag("slimmedMETs","","Analysis"),        
   #met = cms.InputTag("slimmedMETs"),          
   #met = cms.InputTag("slimmedMETsMuClean","","Analysis"),  
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   pfCands = cms.InputTag("packedPFCandidates"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),   
   LHEInfo = cms.InputTag("source"),                        
   rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),    
   pdfWeightTags = cms.VInputTag(),
   bTagEfficiencies = bTagEffMapPars,
   BTagCalibration = BTagCalibrationPars,
   BTagCalibrationReader = BTagCalibrationReaderPars,
   #~ LeptonFastSimScaleFactors = LeptonFastSimScaleFactorMapPars,
   #~ LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars,
   vertexWeights = vertexWeightPars,
   vertexWeightsUp = vertexWeightParsUp,
   vertexWeightsDown = vertexWeightParsDown,                   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
   writeID = cms.untracked.bool(False), 
   doMETUncert = cms.untracked.bool(True),
   
   storeMetFilters = cms.untracked.bool(True),
   metFilterNames=cms.untracked.vstring(                                      
                  "Flag_goodVertices",              
                  "Flag_HBHENoiseFilter",
                  "Flag_HBHENoiseIsoFilter",
                  "Flag_EcalDeadCellTriggerPrimitiveFilter",
                  "Flag_BadPFMuonFilter",
                  "Flag_BadChargedCandidateFilter",
                  "Flag_ecalBadCalibFilter"                      

   ), 
)
