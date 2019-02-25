import FWCore.ParameterSet.Config as cms

isolationDefinitions2016 = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    candSource = cms.InputTag("packedPFCandidates"),
    effAreaElectronEta = cms.vdouble(  0.0,     1.0,    1.479,  2.0,    2.2,    2.3,    2.4),
    effAreaElectronValue = cms.vdouble(0.1752,  0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687),
    effAreaMuonEta = cms.vdouble(  0.0,     0.8,    1.3,    2.0,    2.2),
    effAreaMuonValue = cms.vdouble(0.0735,  0.0619, 0.0465, 0.0433, 0.0577),
    )     
    
isolationDefinitions2017 = cms.PSet(
    rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
    candSource = cms.InputTag("packedPFCandidates"),
    effAreaElectronEta = cms.vdouble(  0.0,     1.0,    1.479,  2.0,    2.2,    2.3,    2.4),
    effAreaElectronValue = cms.vdouble(0.1440,  0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654),
    effAreaMuonEta = cms.vdouble ( 0.0,     0.8,    1.3,    2.0,    2.2),
    effAreaMuonValue = cms.vdouble(0.0566,  0.0562, 0.0363, 0.0119, 0.0064),
    ) 
    
isolationDefinitions2018 = cms.PSet( # not changed for 2018 yet, using 2017 values
    rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
    candSource = cms.InputTag("packedPFCandidates"),
    effAreaElectronEta = cms.vdouble(  0.0,     1.0,    1.479,  2.0,    2.2,    2.3,    2.4),
    effAreaElectronValue = cms.vdouble(0.1440,  0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654),
    effAreaMuonEta = cms.vdouble ( 0.0,     0.8,    1.3,    2.0,    2.2),
    effAreaMuonValue = cms.vdouble(0.0566,  0.0562, 0.0363, 0.0119, 0.0064),
    ) 


