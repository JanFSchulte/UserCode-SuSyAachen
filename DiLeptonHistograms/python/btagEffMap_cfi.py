import FWCore.ParameterSet.Config as cms

bTagEffMap = cms.PSet(
    fullSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/btageff__ttbar_powheg_pythia8_25ns.root'),
    bEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFullSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    fastSimFile = cms.string('${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/btageff__SMS-T1bbbb-T1qqqq_fastsim.root'),
    bEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_b'),
    cEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_c'),
    lightEffFastSimName = cms.string('h2_BTaggingEff_csv_med_Eff_udsg'),
    
    useFastSim = cms.bool(True),
)
