import FWCore.ParameterSet.Config as cms

filterMuons = cms.bool(False)

basicMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                          src = cms.InputTag("cleanLayer1Muons"),
                          cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
                          )

globalMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                       src = cms.InputTag("basicMuons"),
                       cut = cms.string('muonID( "AllGlobalMuons" )')
                       )

qualityMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                            src = cms.InputTag("globalMuons"),
                            cut = cms.string('globalTrack.normalizedChi2 <= 10. & track.numberOfValidHits >= 11.')
                            )

d0Muons = cms.EDFilter("PATMuonD0Selector", filter = filterMuons,
                       src = cms.InputTag("qualityMuons"),
                       d0Min = cms.double(0.2),
                       beamSpotSource  = cms.InputTag("offlineBeamSpot") #offlinePrimeryVertices
                       )

cleanMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                            src = cms.InputTag("d0Muons"),
                            cut = cms.string('')
                            )

bJetMuonProducer = cms.EDProducer('bJetMuonProducer',
    src = cms.InputTag("basicMuons"),
    jetSrc = cms.InputTag("basicJets"),
    dRJetLepton = cms.double(0.2),
    dPhiOppositeJetLepton = cms.double(2.7),
    user_bJetAlgo = cms.string("trackCountingHighPurBJetTags"),
    user_bTagDiscriminator = cms.double(3.),
)
from SuSyAachen.TagAndProbeTreeWriter.isolationFunctor_cfi import isolationDefinitions
isoMuons = cms.EDFilter("PATMuonIsolationSelector", filter = cms.bool(True),
                                           src = cms.InputTag("cleanMuons"),
                                           isolationDefinitions = isolationDefinitions,
                                           method = cms.string("miniIsoEA"),
                                           isoMin = cms.double(-1.),
                                           isoMax = cms.double(0.1),                                         
                                           )   
seqMuons = cms.Sequence(basicMuons * globalMuons * qualityMuons * d0Muons * cleanMuons *
                        isoMuons)
