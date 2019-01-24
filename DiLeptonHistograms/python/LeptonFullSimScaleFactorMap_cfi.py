import FWCore.ParameterSet.Config as cms

LeptonFullSimScaleFactorMapPars2017 = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ScaleFactorMuonID.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ScaleFactorMuonMiniIso.root'),
    #dataMCScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonIP2D.root'),
    #dataMCScaleFactorFile_mu_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonTrackScaleFactors.root'),
    #dataMCScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/ScaleFactorMuonSIP3D.root'),
    
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/ElectronScaleFactors_Run2017.root'),
    #dataMCScaleFactorFile_ele_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/TrackScaleFactorsElectrons.root'),
    
    dataMCScaleFactorHisto_mu_ID = cms.string('NUM_MediumID_DEN_genTracks_pt_abseta'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta'),
    #dataMCScaleFactorHisto_mu_Impact = cms.string('SF'),
    #dataMCScaleFactorHisto_mu_Track = cms.string('mutrksfptg10'),
    #dataMCScaleFactorHisto_mu_SIP3D  = cms.string('SF'),
    
    dataMCScaleFactorHisto_ele_ID = cms.string('Run2017_MVATightTightIP2D3D'),
    dataMCScaleFactorHisto_ele_Iso = cms.string('Run2017_MVAVLooseTightIP2DMini'),
    dataMCScaleFactorHisto_ele_ConvHit = cms.string('Run2017_ConvIHit0'),
    #dataMCScaleFactorHisto_ele_Track = cms.string('EGamma_SF2D')
)

LeptonFullSimScaleFactorMapPars2016 = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonID.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonMiniIso.root'),
    dataMCScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonIP2D.root'),
    #~ dataMCScaleFactorFile_mu_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonTrackScaleFactors.root'),
    dataMCScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorMuonSIP3D.root'),
    
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/ScaleFactorsElectrons.root'),
    dataMCScaleFactorFile_ele_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/TrackScaleFactorsElectrons.root'),
    
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
