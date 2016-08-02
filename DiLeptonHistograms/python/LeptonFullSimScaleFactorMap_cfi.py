import FWCore.ParameterSet.Config as cms

LeptonFullSimScaleFactorMap = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonIDScaleFactors.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonIsoScaleFactors.root'),
    dataMCScaleFactorFile_mu_Impact = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonImpactParameterScaleFactors.root'),
    dataMCScaleFactorFile_mu_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonTrackScaleFactors.root'),
    #~ dataMCScaleFactorFile_mu_SIP3D = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonSIP3DScaleFactors.root'),
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronScaleFactors.root'),
    dataMCScaleFactorFile_ele_Track = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronTrackScaleFactors.root'),
    dataMCScaleFactorHisto_mu_ID = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass'),
    dataMCScaleFactorHisto_mu_Impact = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass'),
    dataMCScaleFactorHisto_mu_Track = cms.string('mutrksfptg10'),
    #~ dataMCScaleFactorHisto_mu_SIP3D = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass'),
    dataMCScaleFactorHisto_ele_ID = cms.string('GsfElectronToTight2D3D'),
    dataMCScaleFactorHisto_ele_Iso = cms.string('MVAVLooseElectronToMini'),
    dataMCScaleFactorHisto_ele_ConvHit = cms.string('MVATightElectronToConvIHit0'),
    dataMCScaleFactorHisto_ele_Track = cms.string('EGamma_SF2D')
)
