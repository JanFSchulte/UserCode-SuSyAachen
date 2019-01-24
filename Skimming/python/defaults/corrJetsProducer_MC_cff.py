# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

def corrJetsProducer_MC(process):


        usePrivateSQlite=True
        
        if usePrivateSQlite:
                era="Fall17_17Nov2017_V32_94X_MC"
                from CondCore.CondDB.CondDB_cfi import CondDB
                CondDBJECFile = CondDB.clone(connect = cms.string('sqlite_file:'+era+'.db'))
                process.jec = cms.ESSource("PoolDBESSource",
                        CondDBJECFile,
                        timetype = cms.string('runnumber'),
                        toGet = cms.VPSet(
                        cms.PSet(
                                record = cms.string('JetCorrectionsRecord'),
                                tag    = cms.string('JetCorrectorParametersCollection_'+era+'_AK4PFchs'),
                                label  = cms.untracked.string('AK4PFchs')
                                ),
                        ## here you add as many jet types as you need
                        ## note that the tag name is specific for the particular sqlite file 
                        ), 
                )
                ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
                process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

        updateJetCollection(
           process,
           jetSource = cms.InputTag('slimmedJets'),
           labelName = 'UpdatedJEC',
           jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
        )
        
        process.seqcorrJetsProducerMC = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)
        process.seqcorrJetsPathMC = cms.Path(process.seqcorrJetsProducerMC)


def corrJetsProducer_MC17(process):


        usePrivateSQlite=True
        
        if usePrivateSQlite:
                era="Fall17_17Nov2017_V32_94X_MC"
                from CondCore.CondDB.CondDB_cfi import CondDB
                CondDBJECFile = CondDB.clone(connect = cms.string('sqlite_file:'+era+'.db'))
                process.jec = cms.ESSource("PoolDBESSource",
                        CondDBJECFile,
                        timetype = cms.string('runnumber'),
                        toGet = cms.VPSet(
                        cms.PSet(
                                record = cms.string('JetCorrectionsRecord'),
                                tag    = cms.string('JetCorrectorParametersCollection_'+era+'_AK4PFchs'),
                                label  = cms.untracked.string('AK4PFchs')
                                ),
                        ## here you add as many jet types as you need
                        ## note that the tag name is specific for the particular sqlite file 
                        ), 
                )
                ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
                process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

        from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

        updateJetCollection(
           process,
           jetSource = cms.InputTag('slimmedJets'),
           labelName = 'UpdatedJEC',
           jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
        )
        
        process.seqcorrJetsProducer_MC17 = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)
        process.seqcorrJetsPath_MC = cms.Path(process.seqcorrJetsProducer_MC17)
