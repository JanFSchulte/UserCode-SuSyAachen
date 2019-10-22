import FWCore.ParameterSet.Config as cms

# order 0: x:pt, y:eta
# order 1: x:eta, y:pt
# order 10: x:pt, y:abseta
# order 11: x:abseta, y:pt

LeptonFullSimScaleFactorMapPars2016 = cms.PSet(
    electronPtThreshold = cms.double(500.0),
    muonPtThreshold = cms.double(120.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ElectronScaleFactors_Run2016.root'), histName = cms.string('Run2016_MVATightTightIP2D3D'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ElectronScaleFactors_Run2016.root'), histName = cms.string('Run2016_Mini'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ElectronScaleFactors_Run2016.root'), histName = cms.string('Run2016_ConvIHit0'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root'), histName = cms.string('EGamma_SF2D'), order = cms.int32(1) ),
    ),
    muonScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonID.root'), histName = cms.string('SF'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonMiniIso.root'), histName = cms.string('SF'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonIP2D.root'), histName = cms.string('SF'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonSIP3D.root'), histName = cms.string('SF'), order = cms.int32(10) ),
    )
)

LeptonFullSimScaleFactorMapPars2017 = cms.PSet(
    electronPtThreshold = cms.double(500.0),
    muonPtThreshold = cms.double(120.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ElectronScaleFactors_Run2017.root'), histName = cms.string('Run2017_MVATightTightIP2D3D'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ElectronScaleFactors_Run2017.root'), histName = cms.string('Run2017_MVAVLooseTightIP2DMini'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ElectronScaleFactors_Run2017.root'), histName = cms.string('Run2017_ConvIHit0'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/EGM2D_runBCDEF_passingRECO.root'), histName = cms.string('EGamma_SF2D'), order = cms.int32(1) ),
    ),
    muonScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ScaleFactorMuonID.root'), histName = cms.string('NUM_MediumPromptID_DEN_genTracks_pt_abseta'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ScaleFactorMuonMiniIso.root'), histName = cms.string('TnP_MC_NUM_MiniIso02Cut_DEN_MediumCutidPromptCut_PAR_pt_eta'), order = cms.int32(10) ),
    )
)

LeptonFullSimScaleFactorMapPars2018 = cms.PSet(
    electronPtThreshold = cms.double(500.0),
    muonPtThreshold = cms.double(120.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/ElectronScaleFactors_Run2018.root'), histName = cms.string('Run2018_MVATightTightIP2D3D'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/ElectronScaleFactors_Run2018.root'), histName = cms.string('Run2018_Mini'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/ElectronScaleFactors_Run2018.root'), histName = cms.string('Run2018_ConvIHit0'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/EGM2D_updatedAll.root'), histName = cms.string('EGamma_SF2D'), order = cms.int32(1) ),
    ),
    muonScaleFactors = cms.VPSet( # no SF for 2018 until UL, so using 2017 SFs
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/ScaleFactorMuonID.root'), histName = cms.string('NUM_MediumPromptID_DEN_genTracks_pt_abseta'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/ScaleFactorMuonMiniIso.root'), histName = cms.string('TnP_MC_NUM_MiniIso02Cut_DEN_MediumCutidPromptCut_PAR_pt_eta'), order = cms.int32(10) ),
    )
)

