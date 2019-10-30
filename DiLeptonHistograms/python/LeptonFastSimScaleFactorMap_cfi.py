import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as cms

# order 0: x:pt, y:eta
# order 1: x:eta, y:pt
# order 10: x:pt, y:abseta
# order 11: x:abseta, y:pt

LeptonFastSimScaleFactorMapPars2016 = cms.PSet(
    electronPtThreshold = cms.double(200.0),
    muonPtThreshold = cms.double(200.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_el_inhit_eq0.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_el_mini01.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_el_tight2d3d.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
    ),
    muonScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_mu_mediumID.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_mu_mediumID_mini02.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/sf_mu_mediumID_tightIP2D.root'), histName = cms.string('histo2D'), order = cms.int32(10) ),
    )
)

LeptonFastSimScaleFactorMapPars2017 = cms.PSet(
    electronPtThreshold = cms.double(500.0),
    muonPtThreshold = cms.double(500.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/detailed_ele_full_fast_sf_17.root'), histName = cms.string('MVATightTightIP2D3D_sf'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/detailed_ele_full_fast_sf_17.root'), histName = cms.string('MVAVLooseTightIP2DMini_sf'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/detailed_ele_full_fast_sf_17.root'), histName = cms.string('ConvIHit0_sf'), order = cms.int32(1) ),
    ),
    muonScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/detailed_mu_full_fast_sf_17.root'), histName = cms.string('miniIso02_MediumPrompt_sf'), order = cms.int32(10) ),
    )
)

LeptonFastSimScaleFactorMapPars2018 = cms.PSet(
    electronPtThreshold = cms.double(500.0),
    muonPtThreshold = cms.double(500.0),
    electronScaleFactors = cms.VPSet(
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/detailed_ele_full_fast_sf_18.root'), histName = cms.string('MVATightTightIP2D3D_sf'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/detailed_ele_full_fast_sf_18.root'), histName = cms.string('MVAVLooseTightIP2DMini_sf'), order = cms.int32(1) ),
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/detailed_ele_full_fast_sf_18.root'), histName = cms.string('ConvIHit0_sf'), order = cms.int32(1) ),
    ),
    muonScaleFactors = cms.VPSet( # no SF for 2018 until UL, so using 2017 SFs
        cms.PSet( fileName = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/detailed_mu_full_fast_sf_18.root'), histName = cms.string('miniIso02_MediumPrompt_sf'), order = cms.int32(10) ),
    )
)
