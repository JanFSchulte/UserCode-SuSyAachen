import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.pdgIdDefinition_cff import defaultPdgIdDefinition
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockA_cfi import vertexWeightsBlockA as vertexWeightParsBlockA
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockB_cfi import vertexWeightsBlockB as vertexWeightParsBlockB
from SuSyAachen.DiLeptonHistograms.vertexWeightsUp_cfi import vertexWeightsUp as vertexWeightParsUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockAUp_cfi import vertexWeightsBlockAUp as vertexWeightParsBlockAUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockBUp_cfi import vertexWeightsBlockBUp as vertexWeightParsBlockBUp
from SuSyAachen.DiLeptonHistograms.vertexWeightsDown_cfi import vertexWeightsDown as vertexWeightParsDown
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockADown_cfi import vertexWeightsBlockADown as vertexWeightParsBlockADown
from SuSyAachen.DiLeptonHistograms.vertexWeightsBlockBDown_cfi import vertexWeightsBlockBDown as vertexWeightParsBlockBDown
from SuSyAachen.DiLeptonHistograms.efficiencies.electronEffPSet_cff import electronCenterEfficiencies as electronEfficiency
from SuSyAachen.DiLeptonHistograms.efficiencies.muonEffPSet_cff import muonCenterEfficiencies as muonEfficiency
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions

MultiLeptonTreesNoTaus = cms.EDAnalyzer("MultiLeptonTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("triggerMatchedPatJetsPF"),	   	   
   bJets = cms.InputTag("triggerMatchedPatJetsPF"),
   met = cms.InputTag("patMETsPF"),  	   
   tcmet = cms.InputTag("patMETsTC"), 
   type1met = cms.InputTag("pfType1CorrectedMet"),
   caloMet = cms.InputTag("patMETsTC"),
   genMetTrue = cms.InputTag("genMetTrue"),	   	   
   vertices = cms.InputTag("offlinePrimaryVertices"),
   tauId = cms.string("byTaNCfrHalfPercent"),
   pfCands = cms.InputTag("particleFlow"),
   genParticles = cms.InputTag("genParticles"),
   rho = cms.InputTag("kt6PFJetsForIsolation", "rho"),	   
   susyVars = cms.VPSet(),
   pdfWeightTags = cms.VInputTag("susyScanPdfWeights:CT10","susyScanPdfWeights:MSTW2008nlo68cl","susyScanPdfWeights:NNPDF10"),
   pdfInfo = cms.InputTag("generator"),		   
   metUncertaintyInputs = cms.VInputTag(cms.InputTag("patPFMet"),cms.InputTag("patPFMetElectronEnDown"),cms.InputTag("patPFMetElectronEnUp"),cms.InputTag("patPFMetJetEnDown"),cms.InputTag("patPFMetJetEnUp"),cms.InputTag("patPFMetJetResDown"),cms.InputTag("patPFMetJetResUp"),cms.InputTag("patPFMetMuonEnDown"),cms.InputTag("patPFMetMuonEnUp"),cms.InputTag("patPFMetTauEnDown"),cms.InputTag("patPFMetTauEnUp"),cms.InputTag("patPFMetUnclusteredEnDown"),cms.InputTag("patPFMetUnclusteredEnUp")),	   
   vertexWeights = vertexWeightPars,
   vertexWeightsBlockA = vertexWeightParsBlockA,
   vertexWeightsBlockB = vertexWeightParsBlockB,
   vertexWeightsUp = vertexWeightParsUp,
   vertexWeightsBlockAUp = vertexWeightParsBlockAUp,
   vertexWeightsBlockBUp = vertexWeightParsBlockBUp,
   vertexWeightsDown = vertexWeightParsDown,
   vertexWeightsBlockADown = vertexWeightParsBlockADown,
   vertexWeightsBlockBDown = vertexWeightParsBlockBDown,	   	   	   	   
   pdgIdDefinition = defaultPdgIdDefinition,
   isolationDefinitions = isolationDefinitions,
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

MultiLeptonTrees = MultiLeptonTreesNoTaus.clone(
    taus = cms.InputTag("triggerMatchedPatTausPF"),
    tauId = cms.string("byTaNCfrHalfPercent"), 
    )

MultiLeptonTreesForLMPoints = MultiLeptonTrees.clone(
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

MultiLeptonTreesmSugra = MultiLeptonTrees.clone(
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

MultiLeptonTreesSimplified = MultiLeptonTrees.clone(
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
MultiLeptonTreesSimplified.vertexWeights.doWeight=False
