import FWCore.ParameterSet.Config as cms
from SuSyAachen.DiLeptonHistograms.vertexWeights_cfi import vertexWeights as vertexWeightPars

DiLeptonTrees = cms.EDAnalyzer("DiLeptonTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("triggerMatchedPatJetsPF"),
   bJets = cms.InputTag("triggerMatchedPatJetsPF"),
   met = cms.InputTag("patMETsPF"),
   vertices = cms.InputTag("offlinePrimaryVertices"),
   tauId = cms.string("byTaNCfrHalfPercent"),
   susyVars = cms.VPSet(),
   pdfWeightTags = cms.VInputTag(),
   vertexWeights = vertexWeightPars,

# NOT USED RIGHT NOW
   fakeRates =  cms.PSet(
        elsectrons = cms.VPSet(),
        muons =  cms.VPSet(),
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

DiLeptonTreesmSugra = DiLeptonTrees.clone(
   susyVars = cms.VPSet(
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanA0"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanM0"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanM12"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScantanbeta"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanCrossSection"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("susyScanNLOCrossSection"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale2"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("susyScanNLOCrossSectionScale05"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("susyScankFactor"), type = cms.string("float")),
#       cms.PSet(var = cms.InputTag("susyScanRun"), type = cms.string("float")),
       cms.PSet(var = cms.InputTag("seqSUSYPARS","susyScanMu"), type = cms.string("int"))
       ),
   pdfWeightTags = cms.VInputTag(
        "susyScanPdfWeights:cteq66",
        "susyScanPdfWeights:MRST2006nnlo",
        "susyScanPdfWeights:NNPDF10"
       )
)
