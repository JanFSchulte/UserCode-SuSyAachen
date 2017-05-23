import FWCore.ParameterSet.Config as cms

LeptonFullSimScaleFactorMap = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonID.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonMiniIso.root'),
    dataMCScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonIP2D.root'),
    #~ dataMCScaleFactorFile_mu_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonTrackScaleFactors.root'),
    dataMCScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonSIP3D.root'),
    
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorsElectrons.root'),
    dataMCScaleFactorFile_ele_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/TrackScaleFactorsElectrons.root'),
    
    dataMCScaleFactorHisto_mu_ID = cms.string('SF'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('SF'),
    dataMCScaleFactorHisto_mu_Impact = cms.string('SF'),
    #~ dataMCScaleFactorHisto_mu_Track = cms.string('mutrksfptg10'),
    dataMCScaleFactorHisto_mu_SIP3D  = cms.string('SF'),
    
    dataMCScaleFactorHisto_ele_ID = cms.string('GsfElectronToMVATightTightIP2DSIP3D4'),
    dataMCScaleFactorHisto_ele_Iso = cms.string('MVAVLooseElectronToMini'),
    dataMCScaleFactorHisto_ele_ConvHit = cms.string('MVATightElectronToConvVetoIHit0'),
    dataMCScaleFactorHisto_ele_Track = cms.string('EGamma_SF2D')
)
