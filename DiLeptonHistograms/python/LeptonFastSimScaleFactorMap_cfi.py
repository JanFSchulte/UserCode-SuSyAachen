import FWCore.ParameterSet.Config as cms

LeptonFastSimScaleFactorMap = cms.PSet(
    FastSimScaleFactorFile_mu_Id = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorMuonID.root'),
    FastSimScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorMuonIso.root'),
    FastSimScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorMuonIP2D.root'),
    FastSimScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorMuonSIP3D.root'),
    
    FastSimScaleFactorFile_ele_Id = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorElectronID.root'),
    FastSimScaleFactorFile_ele_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorElectronIso.root'),
    FastSimScaleFactorFile_ele_ConvHit = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/FastSimScaleFactorElectronIso.root'),
    
    FastSimScaleFactorHisto_mu_Id = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Iso = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Impact = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_SIP3D = cms.string('histo2D'),
    
    FastSimScaleFactorHisto_ele_Id = cms.string('histo2D'),
    FastSimScaleFactorHisto_ele_Iso = cms.string('histo2D'),
    FastSimScaleFactorHisto_ele_ConvHit = cms.string('histo2D')
)
