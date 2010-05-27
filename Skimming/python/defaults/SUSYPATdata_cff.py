import FWCore.ParameterSet.Config as cms
def SUSYPATdata(process):  
    from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands, removeMCDependence
    process.out = cms.OutputModule("PoolOutputModule",
      fileName = cms.untracked.string('dummy.root'),
      outputCommands = cms.untracked.vstring( 'drop *')
    )
#    process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    addDefaultSUSYPAT(process,False,'HLT','Spring10','35x',['IC5Calo','AK5PF'])
    #addDefaultSUSYPAT(process,True,'HLT','Summer09_7TeV_ReReco332') #no up-to-date JetMET corrections for FastSim these were recomendate somewhere...
    #SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process ) 

    del process.out
    process.seqSUSYPATdata = process.susyPatDefaultSequence
                         
