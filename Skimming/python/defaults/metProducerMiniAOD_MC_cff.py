# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms



def metProducerMiniAOD_MC17(process):

        #Run corrected MET maker

        #configurable options =======================================================================
        runOnData=False #data/MC switch
        usePrivateSQlite=True #use external JECs (sqlite file)
        redoPuppi=False # rebuild puppiMET
        #===================================================================

        if usePrivateSQlite:
                import os
                era="Fall17_17Nov2017_V32_94X_MC"
                from CondCore.CondDB.CondDB_cfi import CondDB
                CondDBJECFile = CondDB.clone(connect = cms.string('sqlite_file:'+era+'.db'))
                process.jec = cms.ESSource("PoolDBESSource",
                                       CondDBJECFile,
                                       timetype = cms.string('runnumber'),
                                       toGet =  cms.VPSet(
                    cms.PSet(
                        record = cms.string("JetCorrectionsRecord"),
                        tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                        label= cms.untracked.string("AK4PF")
                        ),
                    cms.PSet(
                        record = cms.string("JetCorrectionsRecord"),
                        tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                        label= cms.untracked.string("AK4PFchs")
                        ),
                    )
                               )
                process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

        ### =====================================================================================================


        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD


        #default configuration for miniAOD reprocessing, change the isData flag to run on data
        #for a full met computation, remove the pfCandColl input
        runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               )

        # temporary fix for prefire issue, ignore jec for certain forward jets
        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

        runMetCorAndUncFromMiniAOD (
                process,
                isData = False, # false for MC
                fixEE2017 = True,
                fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
                postfix = "ModifiedMET"
        )

        if redoPuppi:
                from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
                makePuppiesFromMiniAOD( process );
                
                runMetCorAndUncFromMiniAOD(process,
                                                                 isData=runOnData,
                                                                 pfCandColl=cms.InputTag("puppiForMET"),
                                                                 recoMetFromPFCs=True,
                                                                 reclusterJets=True,
                                                                 jetFlavor="AK4PFPuppi",
                                                                 postfix="Puppi"
                                                                 )

        # end Run corrected MET maker
        
        process.seqmetProducerMiniAOD_MC17 = cms.Sequence(process.fullPatMetSequence*process.fullPatMetSequenceModifiedMET)
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MC17)


def metProducerMiniAOD_MC(process):

        #Run corrected MET maker

        #configurable options =======================================================================
        runOnData=False #data/MC switch
        usePrivateSQlite=True #use external JECs (sqlite file)
        redoPuppi=False # rebuild puppiMET
        #===================================================================

        if usePrivateSQlite:
                import os
                era="Fall17_17Nov2017_V32_94X_MC"
                from CondCore.CondDB.CondDB_cfi import CondDB
                CondDBJECFile = CondDB.clone(connect = cms.string('sqlite_file:'+era+'.db'))
                process.jec = cms.ESSource("PoolDBESSource",
                                       CondDBJECFile,
                                       timetype = cms.string('runnumber'),
                                       toGet =  cms.VPSet(
                    cms.PSet(
                        record = cms.string("JetCorrectionsRecord"),
                        tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                        label= cms.untracked.string("AK4PF")
                        ),
                    cms.PSet(
                        record = cms.string("JetCorrectionsRecord"),
                        tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                        label= cms.untracked.string("AK4PFchs")
                        ),
                    )
                               )
                process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

        ### =====================================================================================================


        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD


        #default configuration for miniAOD reprocessing, change the isData flag to run on data
        #for a full met computation, remove the pfCandColl input
        runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               )

        # temporary fix for prefire issue, ignore jec for certain forward jets
        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

        runMetCorAndUncFromMiniAOD (
                process,
                isData = False, # false for MC
                fixEE2017 = True,
                fixEE2017Params = {'userawPt': True, 'PtThreshold':50.0, 'MinEtaThreshold':2.65, 'MaxEtaThreshold': 3.139} ,
                postfix = "ModifiedMET"
        )

        if redoPuppi:
                from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
                makePuppiesFromMiniAOD( process );
                
                runMetCorAndUncFromMiniAOD(process,
                                                                 isData=runOnData,
                                                                 pfCandColl=cms.InputTag("puppiForMET"),
                                                                 recoMetFromPFCs=True,
                                                                 reclusterJets=True,
                                                                 jetFlavor="AK4PFPuppi",
                                                                 postfix="Puppi"
                                                                 )

        # end Run corrected MET maker
        
        process.seqmetProducerMiniAOD_MC = cms.Sequence(process.fullPatMetSequenceModifiedMET*process.fullPatMetSequence)
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MC)
