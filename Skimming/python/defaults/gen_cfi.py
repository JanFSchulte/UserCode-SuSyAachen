import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("GenParticleSelector", 
        filter = cms.bool(True),
        src = cms.InputTag("genParticles"),
        cut = cms.string('abs( eta ) <= 5.0 & pt >= 5')#GeV
)

signatureFilter = cms.EDFilter("SignatureFilter",
    src = cms.InputTag("genParticles"),
    IDs = cms.vint32( 23 ),
    status = cms.vint32(3),
    daughters = cms.VPSet(
        cms.PSet(
            IDs = cms.vint32( 15 ),
            status = cms.vint32(3)
        ),
        cms.PSet(
            IDs = cms.vint32( -15 ),
            status = cms.vint32(3)
        ),
    )
)