import FWCore.ParameterSet.Config as cms
isolationDefinitions = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
    candSource = cms.InputTag("packedPFCandidates"),
    effAreaElectronEta = cms.vdouble(  0.0,     1.0,    1.479,  2.0,    2.2,    2.3,    2.4),
    effAreaElectronValue = cms.vdouble(0.1440,  0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654),
    effAreaMuonEta = cms.vdouble ( 0.0,     0.8,    1.3,    2.0,    2.2),
    effAreaMuonValue = cms.vdouble(0.0566,  0.0562, 0.0363, 0.0119, 0.0064),
    ) 

isolationDefinitions2017 = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
    candSource = cms.InputTag("packedPFCandidates"),
    effAreaElectronEta = cms.vdouble(  0.0,     1.0,    1.479,  2.0,    2.2,    2.3,    2.4),
    effAreaElectronValue = cms.vdouble(0.1440,  0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654),
    effAreaMuonEta = cms.vdouble ( 0.0,     0.8,    1.3,    2.0,    2.2),
    effAreaMuonValue = cms.vdouble(0.0566,  0.0562, 0.0363, 0.0119, 0.0064),
    ) 


