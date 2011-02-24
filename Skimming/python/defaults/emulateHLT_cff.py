import FWCore.ParameterSet.Config as cms

def emulateHLT(process):
    process.load("RecoJets.JetProducers.ak5CaloJets_cfi")
    process.load("PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff")
    process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
#    process.load("RecoJets.Configuration.JetIDProducers_cff")
    process.patJets.tagInfoSources = cms.VInputTag()
    process.patJets.discriminatorSources = cms.VInputTag()
    process.patJets.addJetID = False
    
    process.hltJet30Ht = cms.EDFilter("PATJetSelector", 
        filter = cms.bool(False),
        src = cms.InputTag("patJets"),
        cut = cms.string('pt >= 30 & emEnergyFraction >= 1e-6 & n90 >= 2')
    )
    
    process.seqemulateHLT = cms.Sequence(process.ak5CaloJets *
        process.ak5JetTracksAssociatorAtVertex *# process.ak5JetID *# process.btagging *
        process.makePatJets + process.hltJet30Ht
    )