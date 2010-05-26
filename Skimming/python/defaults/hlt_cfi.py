import FWCore.ParameterSet.Config as cms

#import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

defaultSelector = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("HLT"),
    HLTPaths = cms.vstring(),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True)    # throw exception on unknown path names
)
