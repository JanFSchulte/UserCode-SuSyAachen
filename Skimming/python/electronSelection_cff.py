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

cleanElectrons = cms.EDFilter("PATElectronSelector", filter = filterElectrons,
  src = cms.InputTag("d0Electrons"),
  cut = cms.string('')
)

isoElectrons = cms.EDFilter("PATElectronSelector", filter = filterElectrons,
  src = cms.InputTag("cleanElectrons"),
  #cut = cms.string('hcalIsoDeposit.candEnergy < 999 &  ecalIsoDeposit.candEnergy < 999')
  cut = cms.string('(trackIso + ecalIso + hcalIso) / pt < 0.2')
)

seqElectrons = cms.Sequence(basicElectrons + qualityElectrons + d0Electrons + cleanElectrons +
                        isoElectrons)
