import FWCore.ParameterSet.Config as cms

LeptonScaleFactorMap = cms.PSet(
    dataMCScaleFactorFile_mu_ID = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonIDScaleFactors.root'),
    dataMCScaleFactorFile_mu_Iso = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonIsoScaleFactors.root'),
    dataMCScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronScaleFactors.root'),
    dataMCScaleFactorHisto_mu_ID = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass'),
    dataMCScaleFactorHisto_mu_Iso = cms.string('pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass'),
    dataMCScaleFactorHisto_ele_ID = cms.string('MVATight_and_TightIP2D_and_TightIP3D'),
    dataMCScaleFactorHisto_ele_Iso = cms.string('MiniIso0p1_vs_AbsEta'),
    FastSimScaleFactorFile_mu = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/muonFastSimScaleFactors.root'),
    FastSimScaleFactorFile_ele = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/electronFastSimScaleFactors.root'),
    FastSimScaleFactorHisto_mu = cms.string('histo3D'),
    FastSimScaleFactorHisto_ele = cms.string('histo3D'),
    
    useFastSim = cms.bool(False),
)
