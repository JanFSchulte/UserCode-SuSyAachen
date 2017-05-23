import FWCore.ParameterSet.Config as cms

LeptonFastSimScaleFactorMap = cms.PSet(
    FastSimScaleFactorFile_mu_Id = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorMuonID.root'),
    FastSimScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorMuonIso.root'),
    FastSimScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorMuonIP2D.root'),
    FastSimScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorMuonSIP3D.root'),
    
    FastSimScaleFactorFile_ele_Id = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorElectronID.root'),
    FastSimScaleFactorFile_ele_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorElectronIso.root'),
    FastSimScaleFactorFile_ele_ConvHit = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/FastSimScaleFactorElectronIso.root'),
    
    FastSimScaleFactorHisto_mu_Id = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Iso = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_Impact = cms.string('histo2D'),
    FastSimScaleFactorHisto_mu_SIP3D = cms.string('histo2D'),
    
    FastSimScaleFactorHisto_ele_Id = cms.string('histo2D'),
    FastSimScaleFactorHisto_ele_Iso = cms.string('histo2D'),
    FastSimScaleFactorHisto_ele_ConvHit = cms.string('histo2D')
)
