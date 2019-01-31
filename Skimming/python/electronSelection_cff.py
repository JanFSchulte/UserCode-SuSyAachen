import FWCore.ParameterSet.Config as cms

filterElectrons = cms.bool(False)

basicElectrons = cms.EDFilter("PATElectronSelector", filter = filterElectrons,
  src = cms.InputTag("cleanLayer1Electrons"),
  cut = cms.string('abs( eta ) <= 2. & pt >= 10')#GeV
)

qualityElectrons = cms.EDFilter("PATElectronSelector", filter = filterElectrons,
  src = cms.InputTag("basicElectrons"),
  cut = cms.string('electronID("eidTight")==1.0')
)

d0Electrons = cms.EDFilter("PATElectronD0Selector", filter = filterElectrons,
  src = cms.InputTag("qualityElectrons"),
  d0Min = cms.double(0.2),
  beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
)


mvaIdElectrons = cms.EDFilter("PATElectronMVAIDSelector", filter = filterElectrons,
  src = cms.InputTag("d0Electrons"),
  idMapSource  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues") #offlinePrimeryVertices
)

isoElectrons = cms.EDFilter("PATElectronSelector", filter = filterElectrons,
  src = cms.InputTag("cleanElectrons"),
  #cut = cms.string('hcalIsoDeposit.candEnergy < 999 &  ecalIsoDeposit.candEnergy < 999')
  cut = cms.string('(trackIso + ecalIso + hcalIso) / pt < 0.4')
)

from SuSyAachen.DiLeptonHistograms.isolationFunctor_cfi import *
effectiveAreaIsoElectrons = cms.EDFilter("PATElectronEffectiveAreaSelector", filter = cms.bool(True),
                                           src = cms.InputTag("cleanElectrons"),
                                           rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
                                           isolationDefinitions = isolationDefinitions,
                                           isoMin = cms.double(-1.),
                                           isoMax = cms.double(0.09),                                         
                                           )
                                           
isoElectrons = cms.EDFilter("PATElectronIsolationSelector", filter = cms.bool(True),
                                           src = cms.InputTag("cleanElectrons"),
                                           isolationDefinitions = isolationDefinitions,
                                           method = cms.string("miniIsoEA"),                                           
                                           isoMin = cms.double(-1.),
                                           isoMax = cms.double(0.1),                                         
                                           )                                           
                                           
noMatchedConversionsElectrons = cms.EDFilter("PATElectronMatchedConversionSelector", filter = cms.bool(True),
                                             src = cms.InputTag("cleanElectrons"),
                                             conversionsSource = cms.InputTag("allConversions"),
                                             beamspotSource = cms.InputTag("offlineBeamSpot"),
                                             )

pfElectronProducer = cms.EDProducer('PfElectronProducer',
  src = cms.InputTag("basicElectrons"),
)

seqElectrons = cms.Sequence(basicElectrons * qualityElectrons * d0Electrons * isoElectrons)
