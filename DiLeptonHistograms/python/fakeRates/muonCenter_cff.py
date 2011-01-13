import FWCore.ParameterSet.Config as cms

fakes =  cms.VPSet(
        cms.PSet(
            weight = cms.double(0.052095),
            etaMin = cms.double(-9999.000000),
            ptMin = cms.double(5.000000),
            ptMax = cms.double(50.000000),
            etaMax = cms.double(-1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.037500),
            etaMin = cms.double(-9999.000000),
            ptMin = cms.double(50.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(-1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.031072),
            etaMin = cms.double(-1.000000),
            ptMin = cms.double(5.000000),
            ptMax = cms.double(50.000000),
            etaMax = cms.double(1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.031802),
            etaMin = cms.double(-1.000000),
            ptMin = cms.double(50.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.054115),
            etaMin = cms.double(1.000000),
            ptMin = cms.double(5.000000),
            ptMax = cms.double(50.000000),
            etaMax = cms.double(9999.000000),
        ),
        cms.PSet(
            weight = cms.double(0.054795),
            etaMin = cms.double(1.000000),
            ptMin = cms.double(50.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(9999.000000),
        ),
)
