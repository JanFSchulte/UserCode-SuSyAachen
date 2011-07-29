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
    jetCorrections = ['L1FastJet','L2Relative', 'L3Absolute']
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

    
    addDefaultSUSYPAT(process,mcInfo,hltName,jetCorrections,mcVersion,jetTypes,doValidation,doExtensiveMatching,doSusyTopProjection)
    #addDefaultSUSYPAT(process,True,'HLT','Summer09_7TeV_ReReco332') #no up-to-date JetMET corrections for FastSim these were recomendate somewhere...
    #SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process ) 

    del process.out
    process.seqSUSYPAT = process.susyPatDefaultSequence
                         
