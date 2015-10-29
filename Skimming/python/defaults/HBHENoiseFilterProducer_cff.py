import FWCore.ParameterSet.Config as cms

def HBHENoiseFilterProducer(process):


	process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
	process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
	process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(False)
	process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
	
	process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
		inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
		reverseDecision = cms.bool(False)
	)
	
	process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
		inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
		reverseDecision = cms.bool(False)
	)
	
	#~ process.seqHBHENoiseFilterProducer = cms.Sequence(process.HBHENoiseFilterResultProducer*
												#~ process.ApplyBaselineHBHENoiseFilter*
												#~ process.ApplyBaselineHBHEIsoNoiseFilter)
	process.seqHBHENoiseFilterProducer = cms.Sequence(process.HBHENoiseFilterResultProducer)
	process.seqHBHENoiseFilterPath = cms.Path(process.seqHBHENoiseFilterProducer)								
