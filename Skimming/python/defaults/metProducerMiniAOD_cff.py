# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from CondCore.CondDB.CondDB_cfi import *
from SuSyAachen.DiLeptonHistograms.jecToUse_cfi import *

def metProducerMiniAOD(process, year, runOnData=True, usePrivateSQlite=False, era=None):
        if usePrivateSQlite:
                import os
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
                
        from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
        #default configuration for miniAOD reprocessing, change the isData flag to run on data
        #for a full met computation, remove the pfCandColl input
        runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               )
        if year == "2017":
                runMetCorAndUncFromMiniAOD(
                        process,
                        isData = runOnData, # false for MC
                        fixEE2017 = True,
                        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
                        postfix = "ModifiedMET"
                )
        
        if year == "2017":
                return cms.Sequence(process.fullPatMetSequence*process.fullPatMetSequenceModifiedMET)
        else:
                return cms.Sequence(process.fullPatMetSequence)
        
def metProducerMiniAOD16(process):
        process.seqmetProducerMiniAOD16 = metProducerMiniAOD(process, "2016", runOnData=True, usePrivateSQlite=dataUseDB["2016"], era=dataDB["2016"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD16)
        
def metProducerMiniAOD17(process):
        process.seqmetProducerMiniAOD17 = metProducerMiniAOD(process, "2017", runOnData=True, usePrivateSQlite=dataUseDB["2017"], era=dataDB["2017"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD17)
        
def metProducerMiniAOD18(process):
        process.seqmetProducerMiniAOD18 = metProducerMiniAOD(process, "2018", runOnData=True, usePrivateSQlite=dataUseDB["2018"], era=dataDB["2018"])
        process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD18)

       


