# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.DiLeptonHistograms.efficiencies.electronEffPSet_cff import electronCenterEfficiencies as electronEfficiency
from SuSyAachen.DiLeptonHistograms.efficiencies.muonEffPSet_cff import muonCenterEfficiencies as muonEfficiency
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions
from SuSyAachen.DiLeptonHistograms.triggerDefinitionMiniAOD_cff import defaultTriggerDefinition as triggerDefinitions
DiLeptonTreesFromMiniAODNoTaus = cms.EDAnalyzer("DiLeptonTreesFromMiniAOD",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("qualityJets"),	   	   
   genJets = cms.InputTag("slimmedGenJets"),	   	   
   #~ genJets = cms.InputTag("selectedPatJetsAK4PFPuppi","genJets"),	   	     	   
   jetsPuppi = cms.InputTag("qualityJetsPuppi"),  	   
   bJets = cms.InputTag("qualityBJets"),
   bJetsPuppi = cms.InputTag("qualityBJetsPuppi"),  	
   met = cms.InputTag("slimmedMETs","","Analysis"),  	
   #~ met = cms.InputTag("slimmedMETsCHS"),  	
   metPuppi = cms.InputTag("slimmedMETsPuppi"),  	     	    
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
   writeID = cms.untracked.bool(True),
   writeTrigger = cms.untracked.bool(True),
   triggerNames=cms.vstring(
						"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
						"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
						"HLT_Mu27_TkMu8_v",
						"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",					
						"HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
						"HLT_PFHT200_v",
						"HLT_PFHT250_v",
						"HLT_PFHT300_v",
						"HLT_PFHT350_v", 
						"HLT_PFHT400_v",
						"HLT_PFHT475_v",
						"HLT_PFHT600_v",
						"HLT_PFHT650_v",
						"HLT_PFHT800_v",
						"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v",
						"HLT_DoubleMu8_Mass8_PFHT300_v",
						"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v",
	), 
     
   NOelectronCorrections = cms.VPSet(
    cms.PSet(
    absEta = cms.double(0.4),
    correction = cms.double(1.00594),
    ),
    cms.PSet(
    absEta = cms.double(0.8),
    correction = cms.double(1.0082),
    ),
    cms.PSet(
    absEta = cms.double(1.2),
    correction = cms.double(1.01092),
    ),
    cms.PSet(
    absEta = cms.double(1.5),
    correction = cms.double(1.01659),
    ),
    cms.PSet(
    absEta = cms.double(2.),
    correction = cms.double(1.01271),
    ),
    cms.PSet(
    absEta = cms.double(2.5),
    correction = cms.double(1.02727),
    ),

    ),                                  

# NOT USED RIGHT NOW
   fakeRates =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  cms.VPSet(),
        taus = cms.VPSet()
   ),
   efficiencies =  cms.PSet(
        electrons = electronEfficiency,
        muons =  muonEfficiency,
        taus = cms.VPSet()
   )
#        cms.VPSet(
#            cms.PSet(
#                etaMin = cms.double(-2.4), etaMax = cms.double(-1.3),
#                ptMin = cms.double(0), ptMax = cms.double(10),
#                
#                weight = cms.double(0.5)
#            ),
#            cms.PSet(
#                etaMin = cms.double(-1.3),  etaMax = cms.double(1.3),
#                ptMin = cms.double(0), ptMax = cms.double(10),
#    
#                weight = cms.double(0.2)
#            ),
#        )
#    )

)

DiLeptonTreesFromMiniAOD = DiLeptonTreesFromMiniAODNoTaus.clone(
    taus = cms.InputTag("triggerMatchedPatTausPF"),
    tauId = cms.string("byTaNCfrHalfPercent"), 
    )

DiLeptonTreesFromMiniAODForLMPoints = DiLeptonTreesFromMiniAOD.clone(
          pdfWeightTags = cms.VInputTag(
#            "susyScanPdfWeights:cteq66",
#            "susyScanPdfWeights:MSTW2008nlo68cl",
#            "susyScanPdfWeights:NNPDF20"
            "susyScanPdfWeights:cteq66",
            "susyScanPdfWeights:MRST2006nnlo",
            "susyScanPdfWeights:NNPDF10"
            )
          )

mSugraVars =  [
    cms.PSet(var = cms.InputTag("susycafscan","susyScanA0"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susycafscan","susyScanM0"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susycafscan","susyScanM12"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susycafscan","susyScantanbeta"), type = cms.string("float")),
      #       cms.PSet(var = cms.InputTag("susyScanRun"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susycafscan","susyScanMu"), type = cms.string("int"))
    ]

xSecVars = [
    cms.PSet(var = cms.InputTag("susycafscan","susyScanLOXSection"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susycafscan","susyScanGenFilterEfficiency"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susyScanNLOCrossSection"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale2"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale05"), type = cms.string("float")),
    cms.PSet(var = cms.InputTag("susyScankFactor"), type = cms.string("float")),
]

DiLeptonTreesFromMiniAODmSugra = DiLeptonTreesFromMiniAOD.clone(
   susyVars = cms.VPSet( mSugraVars + xSecVars),
   pdfWeightTags = cms.VInputTag(
        "susyScanPdfWeights:cteq66",
#        "susyScanPdfWeights:MSTW2008nlo68cl",
#        "susyScanPdfWeights:NNPDF20"
#        "susyScanPdfWeights:MRST2006nnlo",
#        "susyScanPdfWeights:NNPDF10"
       )
)

simplifedVars = cms.VPSet(
        cms.PSet(var = cms.InputTag("susycafscan","susyScanMLSP"), type = cms.string("float")),
            cms.PSet(var = cms.InputTag("susycafscan","susyScanMGL"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("seqSimplifiedPars","susyScanLower"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("seqSimplifiedPars","susyScanHigher"), type = cms.string("float")),
            cms.PSet(var = cms.InputTag("susycafscan","susyScanXChi"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanCrossSection"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("susyScanNLOCrossSection"), type = cms.string("float")),
        #                  cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale2"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale05"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("susyScankFactor"), type = cms.string("float")),
            #       cms.PSet(var = cms.InputTag("susyScanRun"), type = cms.string("float")),
        #    cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanMu"), type = cms.string("int"))
            )

DiLeptonTreesFromMiniAODSimplified = DiLeptonTreesFromMiniAOD.clone(
#    vertexWeights = vertexWeightPars.clone(doWeight = False)
    susyVars = simplifedVars,
    pdfWeightTags = cms.VInputTag(
        "susyScanPdfWeights:cteq66",
#        "susyScanPdfWeights:MSTW2008nlo68cl",
#        "susyScanPdfWeights:NNPDF20"
#    "susyScanPdfWeights:MRST2006nnlo",
#    "susyScanPdfWeights:NNPDF10"
    )
    )
DiLeptonTreesFromMiniAODSimplified.vertexWeights.doWeight=False
