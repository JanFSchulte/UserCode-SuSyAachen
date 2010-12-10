import FWCore.ParameterSet.Config as cms
def SUSYPAT(process):  
    from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands, removeMCDependence
    process.out = cms.OutputModule("PoolOutputModule",
      fileName = cms.untracked.string('dummy.root'),
      outputCommands = cms.untracked.vstring( 'drop *')
    )
#    process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

    mcInfo = True
    hltName = 'HLT'
    jetCorrections = ['L2Relative', 'L3Absolute']
    mcVersion = ''
    jetTypes = ['AK5PF']
    doValidation = False
    doExtensiveMatching = False
    doSusyTopProjection = True

    addKeep = ['drop *_selectedPatJetsPF__PAT',
                        'drop *_selectedPatElectronsPF_*_*',
                        'drop *_selectedPatMuonsPF_*_*',
                        'drop *_selectedPatTausPF_*_*',
                        'drop *_cleanPatJetsAK5Calo__PAT',
                        'drop *_cleanPatJetsIC5Calo__PAT',
                        'drop *_cleanPatJetsAK5PF__PAT',
                        'drop *_cleanPatElectrons_*_*',
                        'drop *_cleanPatMuons_*_*',
                        'drop *_cleanPatTaus_*_*',
                        'drop *_cleanPatPhotons_*_*',
                        'keep *_selectedPat*TriggerMatchPF_*_*',
                        'keep *_cleanPat*TriggerMatch*_*_*',
                        'keep *_patTaus_*_*',
                        'keep *_patMETsTypeIPF_*_' ]

    
    electronHLTMatches = ['HLT_Ele15_LW_L1R','HLT_Ele15_SW_L1R','HLT_Ele15_SW_CaloEleId_L1R','HLT_Ele15_SW_EleId_L1R','HLT_Ele17_SW_TightEleId_L1R','HLT_Ele17_SW_TighterEleId_L1R_v1','HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1','HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2','HLT_Ele22_SW_TighterEleId_L1R_v1','HLT_Ele22_SW_TighterEleId_L1R_v2','HLT_Ele22_SW_TighterEleId_L1R_v3']
    muonHLTMatches = ['HLT_Mu15_v1', 'HLT_Mu_9', 'HLT_Mu_11']
    photonHLTMatches = ['HLT_Photon30']
    tauHLTMatches = ['HLT_Jet15U', 'HLT_Jet30U', 'HLT_Jet50U']
    jetHLTMatches = ['HLT_Jet15U', 'HLT_Jet30U', 'HLT_Jet50U']
    addDefaultSUSYPAT(process,mcInfo,hltName,jetCorrections,mcVersion,jetTypes,doValidation,doExtensiveMatching,doSusyTopProjection,electronHLTMatches,muonHLTMatches,tauHLTMatches,jetHLTMatches,photonHLTMatches)
    #addDefaultSUSYPAT(process,True,'HLT','Summer09_7TeV_ReReco332') #no up-to-date JetMET corrections for FastSim these were recomendate somewhere...
    #SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process ) 

    del process.out
    process.seqSUSYPAT = process.susyPatDefaultSequence
                         
