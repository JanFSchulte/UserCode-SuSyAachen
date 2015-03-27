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

from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronGSFScaleFactors as electronGSFScaleFactors
from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronIDScaleFactors as electronIDScaleFactors

from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronLowerGSFScaleFactors as electronGSFScaleFactorsDown
from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronLowerIDScaleFactors as electronIDScaleFactorsDown

from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronUpperGSFScaleFactors as electronGSFScaleFactorsUp
from SuSyAachen.DiLeptonHistograms.efficiencies.electronScaleFactorsPSet_cff import electronUpperIDScaleFactors as electronIDScaleFactorsUp

from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonTrackingScaleFactors as muonTrackingScaleFactors
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonIDScaleFactors as muonIDScaleFactors
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonIsolationScaleFactors as muonIsolationScaleFactors

from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonLowerTrackingScaleFactors as muonTrackingScaleFactorsDown
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonLowerIDScaleFactors as muonIDScaleFactorsDown
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonLowerIsolationScaleFactors as muonIsolationScaleFactorsDown

from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonUpperTrackingScaleFactors as muonTrackingScaleFactorsUp
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonUpperIDScaleFactors as muonIDScaleFactorsUp
from SuSyAachen.DiLeptonHistograms.efficiencies.muonScaleFactorsPSet_cff import muonUpperIsolationScaleFactors as muonIsolationScaleFactorsUp

from SuSyAachen.DiLeptonHistograms.efficiencies.fastSimPSet_cff import eeScaleFactors as eeScaleFactorsFastSim
from SuSyAachen.DiLeptonHistograms.efficiencies.fastSimPSet_cff import emScaleFactors as emScaleFactorsFastSim
from SuSyAachen.DiLeptonHistograms.efficiencies.fastSimPSet_cff import mmScaleFactors as mmScaleFactorsFastSim

DiLeptonSystematicTreesNoTaus = cms.EDAnalyzer("DiLeptonSystematicTrees",
   electrons = cms.InputTag("triggerMatchedPatElectronsPF"),
   muons = cms.InputTag("triggerMatchedPatMuonsPF"),
#   taus = cms.InputTag("triggerMatchedPatTausPF"),
   jets = cms.InputTag("triggerMatchedPatJetsPF"),
   jetsEnUp = cms.InputTag("shiftedJetsUp"),
   jetsEnDown = cms.InputTag("shiftedJetsDown"),
   jetsSmeared = cms.InputTag("smearedJets"),
   jetsSmearedUp = cms.InputTag("smearedJetsUp"),	
   jetsSmearedDown = cms.InputTag("smearedJetsDown"),		   	   	   	   	   	   	   	   	   	   
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
   pdfWeightTags = cms.VInputTag(),
   #pdfWeightTags = cms.VInputTag(),	   
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
   ),
   IDScaleFactors =  cms.PSet(
        electrons = electronIDScaleFactors,
        muons =  muonIDScaleFactors,
        taus = cms.VPSet()
   ),
   GSFScaleFactors =  cms.PSet(
        electrons = electronGSFScaleFactors,
        muons =  cms.VPSet(),
        taus = cms.VPSet()
   ),		   
   trackingScaleFactors =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonTrackingScaleFactors,
        taus = cms.VPSet()
   ),
   isolationScaleFactors =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonIsolationScaleFactors,
        taus = cms.VPSet()
   ),	
   IDScaleFactorsDown =  cms.PSet(
        electrons = electronIDScaleFactorsDown,
        muons =  muonIDScaleFactorsDown,
        taus = cms.VPSet()
   ),
   GSFScaleFactorsDown =  cms.PSet(
        electrons = electronGSFScaleFactorsDown,
        muons =  cms.VPSet(),
        taus = cms.VPSet()
   ),		   
   trackingScaleFactorsDown =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonTrackingScaleFactorsDown,
        taus = cms.VPSet()
   ),
   isolationScaleFactorsDown =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonIsolationScaleFactorsDown,
        taus = cms.VPSet()
   ),	
   IDScaleFactorsUp =  cms.PSet(
        electrons = electronIDScaleFactorsUp,
        muons =  muonIDScaleFactorsUp,
        taus = cms.VPSet()
   ),
   GSFScaleFactorsUp =  cms.PSet(
        electrons = electronGSFScaleFactorsUp,
        muons =  cms.VPSet(),
        taus = cms.VPSet()
   ),		   
   trackingScaleFactorsUp =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonTrackingScaleFactorsUp,
        taus = cms.VPSet()
   ),
   isolationScaleFactorsUp =  cms.PSet(
        electrons = cms.VPSet(),
        muons =  muonIsolationScaleFactorsUp,
        taus = cms.VPSet()
   ),
		   			   		   		   		   		   
   fastSimScaleFactors =  cms.PSet(
        dielectrons = eeScaleFactorsFastSim,
        electronmuons =  emScaleFactorsFastSim,
        dimuons = mmScaleFactorsFastSim,
	ditaus = cms.VPSet(),
	electrontaus = cms.VPSet(),
	muontaus = cms.VPSet(),		
   ),						       
			    
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

DiLeptonSystematicTrees = DiLeptonSystematicTreesNoTaus.clone(
    taus = cms.InputTag("triggerMatchedPatTausPF"),
    tauId = cms.string("byTaNCfrHalfPercent"), 
    )

DiLeptonSystematicTreesForLMPoints = DiLeptonSystematicTrees.clone(
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

DiLeptonSystematicTreesmSugra = DiLeptonSystematicTrees.clone(
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

DiLeptonSystematicTreesSimplified = DiLeptonSystematicTrees.clone(
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
DiLeptonSystematicTreesSimplified.vertexWeights.doWeight=False
