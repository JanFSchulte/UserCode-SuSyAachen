import FWCore.ParameterSet.Config as cms
def SUSYPATdata2012(process):  
    process.out = cms.OutputModule("PoolOutputModule",
      fileName = cms.untracked.string('dummy.root'),
      outputCommands = cms.untracked.vstring( 'drop *')
    )


    ### AACHEN specific, need better place for this #####
    from SuSyAachen.Configuration.SusyPAT_data53PromptReco_cfg import susyPatDefaultSequence
    
    del process.out
    process.seqSUSYPATdata = process.susyPatDefaultSequence

    
