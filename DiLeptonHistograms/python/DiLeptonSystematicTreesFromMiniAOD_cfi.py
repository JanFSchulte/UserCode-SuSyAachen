# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
#~ from SuSyAachen.DiLeptonHistograms.LeptonFastSimScaleFactorMap_cfi import LeptonFastSimScaleFactorMap as LeptonFastSimScaleFactorMapPars
#~ from SuSyAachen.DiLeptonHistograms.LeptonFullSimScaleFactorMap_cfi import LeptonFullSimScaleFactorMap as LeptonFullSimScaleFactorMapPars
from SuSyAachen.DiLeptonHistograms.btagEffMap_cfi import *
from SuSyAachen.DiLeptonHistograms.BTagCalibration_cfi import *
from SuSyAachen.DiLeptonHistograms.BTagCalibrationReader_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignal_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignalUp_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeightsSignalDown_cfi import *
from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import isolationDefinitions
from SuSyAachen.DiLeptonHistograms.metFilterLists_cfi import *
from SuSyAachen.DiLeptonHistograms.triggerLists_cfi import *

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
   bTagEfficiencies = bTagEffMapParsFastSim2017,
   BTagCalibration = BTagCalibrationPars2017,
   BTagCalibrationReader = BTagCalibrationReaderPars2017,
   #~ LeptonFastSimScaleFactors = LeptonFastSimScaleFactorMapPars,
   #~ LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars,
   vertexWeights = vertexWeightsPars2017,
   vertexWeightsUp = vertexWeightsParsUp2017,
   vertexWeightsDown = vertexWeightsParsDown2017,                         
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
   writeID = cms.untracked.bool(False), 
   doMETUncert = cms.untracked.bool(True),
   
   storeMetFilters = cms.untracked.bool(True),
   metFilterNames=metFilterNamesFastSim2017,
)
