# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.LeptonFullSimScaleFactorMap_cfi import *
from SuSyAachen.DiLeptonHistograms.btagEffMap_cfi import *
from SuSyAachen.DiLeptonHistograms.BTagCalibration_cfi import *
from SuSyAachen.DiLeptonHistograms.BTagCalibrationReader_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import *
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import *
from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import *
from SuSyAachen.DiLeptonHistograms.metFilterLists_cfi import *
from SuSyAachen.DiLeptonHistograms.triggerLists_cfi import *

DiLeptonTreesFromMiniAODNoTaus = cms.EDAnalyzer("DiLeptonTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   looseElectrons = cms.InputTag("LooseElectrons"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   looseMuons = cms.InputTag("LooseMuons"),
   fatJets = cms.InputTag("updatedPatJetsAK8"),
   jets = cms.InputTag("qualityJets"),          
   genJets = cms.InputTag("slimmedGenJets"),          
   looseBJets = cms.InputTag("looseBJets"),
   bJets = cms.InputTag("qualityBJets"),
   bJets35 = cms.InputTag("qualityBJets35"),    
   met = cms.InputTag("slimmedMETsModifiedMET", "", "Analysis"),                             
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   isoTracks = cms.InputTag("isolatedTracks"),
   genParticles = cms.InputTag("prunedGenParticles"),
   pdfInfo = cms.InputTag("generator"),   
   LHEInfo = cms.InputTag("externalLHEProducer"),                       
   rho = cms.InputTag("fixedGridRhoFastjetAll"),    
   storeMetFilters = cms.untracked.bool(True),
   pdfWeightTags = cms.VInputTag(),
   bTagEfficiencies = bTagEffMapPars2017,
   BTagCalibration = BTagCalibrationPars2017,
   BTagCalibrationReader = BTagCalibrationReaderPars2017,
   LeptonFullSimScaleFactors = LeptonFullSimScaleFactorMapPars2017,
   vertexWeights = vertexWeightsPars2017,
   vertexWeightsUp = vertexWeightsParsUp2017,
   vertexWeightsDown = vertexWeightsParsDown2017,                   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions2017,
   writeID = cms.untracked.bool(False),
   writeTrigger = cms.untracked.bool(True),
   doMETUncert = cms.untracked.bool(False),
   
   eeTriggerNames=eeTriggerNames2017,
   emTriggerNames=emTriggerNames2017,
   mmTriggerNames=mmTriggerNames2017,  
   htTriggerNames=htTriggerNames2017,
   metTriggerNames=metTriggerNames2017,
   
   metFilterNames=metFilterNamesData2017,
)
