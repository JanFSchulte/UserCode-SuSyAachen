import FWCore.ParameterSet.Config as cms

LeptonFastSimScaleFactorMap = cms.PSet(
    FastSimScaleFactorFile_mu_Id = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonFastSimScaleFactorsId.root'),
    FastSimScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonFastSimScaleFactorsIso.root'),
    FastSimScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonFastSimScaleFactorsImpact.root'),
    FastSimScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronFastSimScaleFactors.root'),
    
    FastSimScaleFactorHisto_mu_Id = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Iso = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Impact = cms.string('histo2D'),
    FastSimScaleFactorHisto_ele = cms.string('histo2D')
)
