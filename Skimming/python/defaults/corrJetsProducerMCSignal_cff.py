# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

def corrJetsProducerMCSignal(process):


	usePrivateSQlite=True
	
	if usePrivateSQlite:
		era="Spring16_25nsFastSimMC_V1"
		#~ process.load("CondCore.DBCommon.CondDBCommon_cfi")
		#~ from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
		from CondCore.CondDB.CondDB_cfi import *
		CondDBJECFile = CondDB.clone(connect = cms.string('sqlite_file:'+era+'.db'))
		process.jec = cms.ESSource("PoolDBESSource",
		      CondDBJECFile,
		      #~ DBParameters = cms.PSet(
		        #~ messageLevel = cms.untracked.int32(0)
		        #~ ),		      
		      timetype = cms.string('runnumber'),
		      toGet = cms.VPSet(
		      cms.PSet(
		            record = cms.string('JetCorrectionsRecord'),
		            tag    = cms.string('JetCorrectorParametersCollection_'+era+'_AK4PFchs'),
		            #~ tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV6_MC_AK4PFchs'),
		            label  = cms.untracked.string('AK4PFchs')
		            ),
		
		
			## here you add as many jet types as you need
			## note that the tag name is specific for the particular sqlite file 
			 ), 
		)
		## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
		process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')




	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
	process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
  	src = cms.InputTag("slimmedJets"),
  	levels = ['L1FastJet', 
    	    'L2Relative', 
    	    'L3Absolute'],
 	 payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
	
	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets	
	process.patJetsReapplyJEC = updatedPatJets.clone(
  	jetSource = cms.InputTag("slimmedJets"),
  	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  	)


	
	process.seqcorrJetsProducerMCSignal = cms.Sequence(process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC )
	process.seqcorrJetsPathMC = cms.Path(process.seqcorrJetsProducerMCSignal)
