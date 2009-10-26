import FWCore.ParameterSet.Config as cms

filterMuons = cms.bool(False)

etaMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                        src = cms.InputTag("cleanLayer1Muons"),
                        cut = cms.string('abs( eta ) <= 2.5')
                        )

ptMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                       src = cms.InputTag("etaMuons"),
                       cut = cms.string('pt >= 10')#GeV                       
                       )

globalMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                       src = cms.InputTag("ptMuons"),
                       cut = cms.string('muonID( "AllGlobalMuons" )')
                       )

qualityMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                            src = cms.InputTag("globalMuons"),
                            cut = cms.string('globalTrack.normalizedChi2 < 10. & track.numberOfValidHits >= 11.')
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

isoMuons = cms.EDFilter("PATMuonSelector", filter = filterMuons,
                           src = cms.InputTag("cleanMuons"),
                           #cut = cms.string('hcalIsoDeposit.candEnergy < 999 &  ecalIsoDeposit.candEnergy < 999')
                           cut = cms.string('(trackIso + ecalIso + hcalIso) / pt < 0.2')
                           )

seqMuons = cms.Sequence(etaMuons + ptMuons + globalMuons + qualityMuons + d0Muons + cleanMuons +
                        isoMuons)
