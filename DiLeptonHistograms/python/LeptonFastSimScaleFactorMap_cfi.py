import FWCore.ParameterSet.Config as cms

LeptonFastSimScaleFactorMap = cms.PSet(
    FastSimScaleFactorFile_mu = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonFastSimScaleFactors.root'),
    FastSimScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronFastSimScaleFactors.root'),
    FastSimScaleFactorHisto_mu = cms.string('histo3D'),
    FastSimScaleFactorHisto_ele = cms.string('histo3D')
)
