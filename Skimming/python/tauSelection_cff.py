import FWCore.ParameterSet.Config as cms

filterTaus = cms.bool(False)

basicTaus = cms.EDFilter("PATTauSelector", filter = filterTaus,
  src = cms.InputTag("cleanLayer1Taus"),
  cut = cms.string('abs( eta ) <= 2.5 & pt >= 7')#GeV
 )

qualityTaus = cms.EDFilter("PATTauSelector", filter = filterTaus,
  src = cms.InputTag("basicTaus"),
  cut = cms.string('leadPFChargedHadrCand.trackRef.numberOfValidHits >= 11' 
                  +'&(signalPFChargedHadrCands.size = 1 | signalPFChargedHadrCands.size = 3)'
                  +'& abs(charge) = 1')
)#leadTrack.numberOfValidHits >= 11. & 

cleanTaus = cms.EDFilter("PATTauSelector", filter = filterTaus,
  src = cms.InputTag("qualityTaus"),
  cut = cms.string('')
)

isoTaus = cms.EDFilter("PATTauSelector", filter = filterTaus,
  src = cms.InputTag("cleanTaus"),
  #cut = cms.string('hcalIsoDeposit.candEnergy < 999 &  ecalIsoDeposit.candEnergy < 999')
  #cut = cms.string('(trackIso + ecalIso + hcalIso) / pt < 0.2')
  #cut = cms.string('isolationTracksPtSum <= 1 & isolationTracksPtSum/pt <= 0.05') 
  cut = cms.string('isolationPFChargedHadrCandsPtSum <= 1 & isolationPFChargedHadrCandsPtSum/pt <= 0.05')
)

seqTaus = cms.Sequence(basicTaus + qualityTaus + cleanTaus +
                        isoTaus)
