import FWCore.ParameterSet.Config as cms


bTagEffMapPars2016 = cms.PSet(
    fullSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/btageff__ttbar_powheg_pythia8_25ns_Moriond17_deepCSV.root'),
    bEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    fastSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2016/btageff__SMS-T1tttt_2016_80X_deepCSV.root'),
    bEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    
    useFastSim = cms.bool(False),
)

bTagEffMapPars2017 = cms.PSet(
    fullSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/btageff__ttbar_amc_94X_deepCSV.root'),
    bEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    fastSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2017/btageff__SMS-T1tttt_2017_94X_deepCSV.root'),
    bEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    
    useFastSim = cms.bool(False),
)

bTagEffMapPars2018 = cms.PSet(
    fullSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/btageff__ttbar_amc_102X_deepCSV.root'),
    bEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    fastSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/2018/btageff__DeepCSV_SMS_T2tt_fastsim_Autumn18.root'),
    bEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    
    useFastSim = cms.bool(False),
)


bTagEffMapParsFastSim2016 = bTagEffMapPars2016.clone(
    useFastSim = cms.bool(True)
)

bTagEffMapParsFastSim2017 = bTagEffMapPars2017.clone(
    useFastSim = cms.bool(True)
)

bTagEffMapParsFastSim2018 = bTagEffMapPars2018.clone(
    useFastSim = cms.bool(True)
)

