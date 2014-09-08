import FWCore.ParameterSet.Config as cms
def SUSYPATmcSignal(process):  
    from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands, removeMCDependence
    process.out = cms.OutputModule("PoolOutputModule",
      fileName = cms.untracked.string('dummy.root'),
      outputCommands = cms.untracked.vstring( 'drop *')
    )

    mcInfo = True
    hltName = 'HLT'
    jetCorrections = ['L1FastJet', "L2Relative", "L3Absolute"]
    mcVersion = ''
    jetTypes = ['AK5PF']
    doValidation = False
    doExtensiveMatching = False
    doSusyTopProjection = True
    doType1MetCorrection=False

    #Apply SUSYPAT
    addDefaultSUSYPAT(process,mcInfo,hltName,jetCorrections, mcVersion, jetTypes, doValidation, doExtensiveMatching, doSusyTopProjection,doType1MetCorrection,useMETFilters=False)
    SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )
    ############################## END SUSYPAT specifics ####################################

    ### AACHEN specific, need better place for this #####
    #from SuSyAachen.Configuration.AachenSUSYPAT_cff import reduceEventsize, additionalTaus
    #additionalTaus(process,postfix="PF")
    #reduceEventsize(process)
    
    del process.out
    process.seqSUSYPATmcSignal = process.susyPatDefaultSequence

    
