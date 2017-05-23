# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
from CondCore.DBCommon.CondDBSetup_cfi import *

def metProducerMiniAOD_MCSignal(process):

    ### from Vince: https://github.com/cmstas/NtupleMaker/blob/CMS3_V07-04-08/test/DataProduction2015_NoFilter_RECO_cfg.py#L87-L168 
	#Run corrected MET maker

	#configurable options =======================================================================
	runOnData=False #data/MC switch
	usePrivateSQlite=True #use external JECs (sqlite file)
	useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false
	applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or for analysis checks 	and developments (not the official recommendation!).
	redoPuppi=False # rebuild puppiMET
	#===================================================================

	if usePrivateSQlite:
		import os
		era="Spring16_25nsFastSimMC_V1"
		from CondCore.CondDB.CondDB_cfi import *
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


	### ---------------------------------------------------------------------------
	### Removing the HF from the MET computation
	### ---------------------------------------------------------------------------
	#~ if not useHFCandidates:
			#~ process.noHFCands = cms.EDFilter("CandPtrSelector",
    	                                 #~ src=cms.InputTag("packedPFCandidates"),
    	                                 #~ cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
    	                                 #~ )

	#jets are rebuilt from those candidates by the tools, no need to do anything else
	### =================================================================================

	from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD


	#default configuration for miniAOD reprocessing, change the isData flag to run on data
	#for a full met computation, remove the pfCandColl input
	runMetCorAndUncFromMiniAOD(process,
    	                       isData=runOnData,
    	                       )

	#~ if not useHFCandidates:
			#~ runMetCorAndUncFromMiniAOD(process,
    	                           #~ isData=runOnData,
    	                           #~ pfCandColl=cms.InputTag("noHFCands"),
    	                           #~ postfix="NoHF"
    	                           #~ )
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

	
	### -------------------------------------------------------------------
	### the lines below remove the L2L3 residual corrections when processing data
	### -------------------------------------------------------------------
	if not applyResiduals:
			process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
			process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
			process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
			process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
			process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
			process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")

			#~ if not useHFCandidates:
				  #~ process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
				  #~ process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
				  #~ process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
				  #~ process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
				  #~ process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
				  #~ process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
	### ------------------------------------------------------------------

	# end Run corrected MET maker


	
	process.seqmetProducerMiniAOD_MCSignal = cms.Sequence()
	process.seqmetProducerMiniAODPath = cms.Path(process.seqmetProducerMiniAOD_MCSignal)
