import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter(
    "BadPFMuonFilterMiniAOD",
    PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
    src  = cms.InputTag("muons"),   # Collection to test
    #~ muons = src,
    taggingMode   = cms.bool(False),
    debug         = cms.bool(False),
    algo          = cms.int32(14),
    minDZ         = cms.double(0.1),              # dz threshold on PF muons to consider; this is not used
    minMuPt       = cms.double(100),               # pt threshold on PF muons 
    minTrkPtError  = cms.double(0.5),               # threshold on inner track pt Error
)
