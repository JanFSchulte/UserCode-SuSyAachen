import FWCore.ParameterSet.Config as cms

def prefireWeightProducer17(process):
        from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
        process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
            DataEra = cms.string("2017BtoF"),
            UseJetEMPt = cms.bool(False),
            PrefiringRateSystematicUncty = cms.double(0.2),
            SkipWarnings = False)


        
        process.seqprefireWeightProducer17 = cms.Sequence(process.prefiringweight)
        process.prefireWeightPath = cms.Path(process.seqprefireWeightProducer17)


def prefireWeightProducer16(process):
        from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
        process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
            DataEra = cms.string("2016BtoH"),
            UseJetEMPt = cms.bool(False),
            PrefiringRateSystematicUncty = cms.double(0.2),
            SkipWarnings = False)


        process.seqprefireWeightProducer16 = cms.Sequence(process.prefiringweight)
        process.prefireWeightPath = cms.Path(process.seqprefireWeightProducer16)
