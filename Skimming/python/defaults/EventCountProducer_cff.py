import FWCore.ParameterSet.Config as cms
def EventCountProducer(process): 
    process.seqEventCountProducer = cms.EDProducer("EventCountProducer")
