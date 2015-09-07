import FWCore.ParameterSet.Config as cms

def corrJetsProducer(process):


	process.load("CondCore.DBCommon.CondDBCommon_cfi")
	from CondCore.DBCommon.CondDBSetup_cfi import *	
	process.jec = cms.ESSource("PoolDBESSource",
	      DBParameters = cms.PSet(
	        messageLevel = cms.untracked.int32(0)
	        ),
	      timetype = cms.string('runnumber'),
	      toGet = cms.VPSet(
	      cms.PSet(
	            record = cms.string('JetCorrectionsRecord'),
	            tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV4_DATA_AK4PFchs'),
	            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V3_MC_AK5PF'),
	            label  = cms.untracked.string('AK4PFchs')
	            ),


      	## here you add as many jet types as you need
      	## note that the tag name is specific for the particular sqlite file 
     	 ), 
      	connect = cms.string('sqlite:Summer15_50nsV4_DATA.db')
    	 # uncomment above tag lines and this comment to use MC JEC
    	 # connect = cms.string('sqlite:Summer12_V7_MC.db')
	)
	## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
	process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')




	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
	process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
  	src = cms.InputTag("slimmedJets"),
  	levels = ['L1FastJet', 
    	    'L2Relative', 
    	    'L3Absolute',
    	    'L2L3Residuals'],
 	 payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
	
	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated	
	process.patJetsReapplyJEC = process.patJetsUpdated.clone(
  	jetSource = cms.InputTag("slimmedJets"),
  	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  	)


	
	process.seqCorrJetsProducer = cms.Sequence(process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )
	process.seqCorrJetsPath = cms.Path(process.seqCorrJetsProducer)
