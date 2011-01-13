import FWCore.ParameterSet.Config as cms

fakes =  cms.VPSet(
        cms.PSet(
            weight = cms.double(0.033307),
            etaMin = cms.double(-9999.000000),
            ptMin = cms.double(10.000000),
            ptMax = cms.double(20.000000),
            etaMax = cms.double(-1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.033195),
            etaMin = cms.double(-9999.000000),
            ptMin = cms.double(20.000000),
            ptMax = cms.double(60.000000),
            etaMax = cms.double(-1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.123727),
            etaMin = cms.double(-9999.000000),
            ptMin = cms.double(60.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(-1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.017574),
            etaMin = cms.double(-1.000000),
            ptMin = cms.double(10.000000),
            ptMax = cms.double(20.000000),
            etaMax = cms.double(1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.015595),
            etaMin = cms.double(-1.000000),
            ptMin = cms.double(20.000000),
            ptMax = cms.double(60.000000),
            etaMax = cms.double(1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.096984),
            etaMin = cms.double(-1.000000),
            ptMin = cms.double(60.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(1.000000),
        ),
        cms.PSet(
            weight = cms.double(0.029024),
            etaMin = cms.double(1.000000),
            ptMin = cms.double(10.000000),
            ptMax = cms.double(20.000000),
            etaMax = cms.double(9999.000000),
        ),
        cms.PSet(
            weight = cms.double(0.028736),
            etaMin = cms.double(1.000000),
            ptMin = cms.double(20.000000),
            ptMax = cms.double(60.000000),
            etaMax = cms.double(9999.000000),
        ),
        cms.PSet(
            weight = cms.double(0.172440),
            etaMin = cms.double(1.000000),
            ptMin = cms.double(60.000000),
            ptMax = cms.double(99995.000000),
            etaMax = cms.double(9999.000000),
        ),
)
