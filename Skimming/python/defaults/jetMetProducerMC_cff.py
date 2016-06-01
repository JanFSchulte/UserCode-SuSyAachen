import FWCore.ParameterSet.Config as cms

def clean_met_(met):
    del met.t01Variation
    del met.t1Uncertainties
    del met.t1SmearedVarsAndUncs
    del met.tXYUncForRaw
    del met.tXYUncForT1
    del met.tXYUncForT01
    del met.tXYUncForT1Smear
    del met.tXYUncForT01Smear

def get_cone_size(algo):
    import re
    cone_size = re.search('(\d+)', algo)
    if cone_size is None:
        raise ValueError('Cannot extract cone size from algorithm name')

    return int(cone_size.group(1))

def useJECFromDB(process, db):
    process.load("CondCore.DBCommon.CondDBCommon_cfi")

    print("Using database %r for JECs" % db)

    process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(0)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(),

            connect = cms.string('sqlite:%s' % db)
         
            )

    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

def checkForTag(db_file, tag):
    import sqlite3

    db_file = db_file.replace('sqlite:', '')

    connection = sqlite3.connect(db_file)
    
    res = connection.execute('select TAG_NAME from IOV where TAG_NAME=?', tag).fetchall()

    return len(res) != 0

def appendJECToDB(process, payload, prefix):

    for set in process.jec.toGet:
        if set.label == payload:
            return

    tag = 'JetCorrectorParametersCollection_%s_%s' % (prefix, payload)
    if not checkForTag(process.jec.connect.value(), (tag,)):
        print("WARNING: The JEC payload %r is not present in the database you want to use. Corrections for this payload will be loaded from the Global Tag" % payload)
        return

    process.jec.toGet += [cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(tag),
            label  = cms.untracked.string(payload)
            )]

def jetMetProducerMC(process):
	readJECFromDB=False 
	jec_database=None 
	jec_db_prefix=None
	isMC = True

	process.load('CommonTools.PileupAlgos.Puppi_cff')
    
	process.puppi.candName = cms.InputTag('PackedPFCandidates')
	process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
	process.puppi.useExistingWeights = True

	if readJECFromDB and (not jec_database or not jec_db_prefix):
		raise LogicError("You must specify the parameters jec_database and jec_db_prefix when reading JEC from DB")
	
	if readJECFromDB:
		useJECFromDB(process, jec_database)
		
	def get_cone_size(algo):
		import re
		cone_size = re.search('(\d+)$', algo)
		if cone_size is None:
			raise ValueError('Cannot extract cone size from algorithm name')
	
		return int(cone_size.group(1))

	def get_jec_payload(algo, pu_method):
	
		# FIXME: Until PUPPI and SK payloads are in the GT, use CHS corrections
		jec_payloads = {
				'Puppi': 'AK%dPFchs',
				'CHS': 'AK%dPFchs',
				'SK': 'AK%dPFchs',
				'': 'AK%dPF',
				}
	
	
		cone_size = get_cone_size(algo)
	
		if not pu_method in jec_payloads:
			print('WARNING: JEC payload not found for method %r. Using default one.' % pu_method)
			return 'None'
	
		return jec_payloads[pu_method] % cone_size

	def get_jec_levels(pu_method):
	
		jec_levels = {
				'Puppi': ['L2Relative', 'L3Absolute'],
				'CHS': ['L1FastJet', 'L2Relative', 'L3Absolute'],
				'SK': ['L2Relative', 'L3Absolute'],
				'': ['L1FastJet', 'L2Relative', 'L3Absolute'],
				}
	
	
		if not pu_method in jec_levels:
			print('WARNING: JEC levels not found for method %r. Using default ones.' % pu_method)
			return ['None']
	
		return jec_levels[pu_method]
	
	jetsCollections = {
			#~ 'AK1': {
				#~ 'algo': 'ak1',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK2': {
				#~ 'algo': 'ak2',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK3': {
				#~ 'algo': 'ak3',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	
			'AK4': {
				'algo': 'ak4',
				'pu_methods': ['Puppi', 'CHS', 'SK', ''],
				'pu_jet_id': False,
				'qg_tagger': False,
				},
	
			#~ 'AK5': {
				#~ 'algo': 'ak5',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK6': {
				#~ 'algo': 'ak6',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK7': {
				#~ 'algo': 'ak7',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK8': {
				#~ 'algo': 'ak8',
				#~ 'pu_methods': ['Puppi', 'CHS', 'SK', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK9': {
				#~ 'algo': 'ak9',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'AK10': {
				#~ 'algo': 'ak10',
				#~ 'pu_methods': ['Puppi', 'CHS', ''],
				#~ 'pu_jet_id': False,
				#~ },
	#~ 
			#~ 'CA10': {
				#~ 'algo': 'ca10',
				#~ 'pu_methods': ['Puppi', 'CHS', 'SK', ''],
				#~ 'pu_jet_id': False,
				#~ },
	
			}


	# Jet corrections
	process.load('JetMETCorrections.Configuration.JetCorrectorsAllAlgos_cff')
	
	from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
	from PhysicsTools.PatAlgos.tools.helpers import loadWithPostfix, applyPostfix
	
	process.load('RecoJets.JetProducers.QGTagger_cfi')
	
	for name, params in jetsCollections.items():
		for index, pu_method in enumerate(params['pu_methods']):
			# Add the jet collection
	
			jec_payload = get_jec_payload(params['algo'], pu_method)
			jec_levels = get_jec_levels(pu_method)
	
			if readJECFromDB:
				appendJECToDB(process, jec_payload, jec_db_prefix)
	
			jetToolbox(process, params['algo'], 'dummy', 'out', runOnMC=isMC, PUMethod = pu_method, JETCorrPayload = jec_payload, JETCorrLevels = jec_levels )
	
	
			algo = params['algo'].upper()
			jetCollection = '%sPFJets%s' % (params['algo'], pu_method)
			postfix = '%sPF%s' % (algo, pu_method)
	
			# FIXME: PU Jet id is not working with puppi jets or SK jets
			if params['pu_jet_id'] and pu_method != 'Puppi' and pu_method != 'SK':
	
				# PU jet Id
				loadWithPostfix(process, 'RecoJets.JetProducers.pileupjetidproducer_cfi', postfix)
				applyPostfix(process, "pileupJetIdEvaluator", postfix).jets = cms.InputTag(jetCollection)
				applyPostfix(process, "pileupJetIdCalculator", postfix).jets = cms.InputTag(jetCollection)
				applyPostfix(process, "pileupJetIdEvaluator", postfix).rho = cms.InputTag("fixedGridRhoFastjetAll")
				applyPostfix(process, "pileupJetIdEvaluator", postfix).vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
				applyPostfix(process, "pileupJetIdCalculator", postfix).rho = cms.InputTag("fixedGridRhoFastjetAll")
				applyPostfix(process, "pileupJetIdCalculator", postfix).vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
	
				# Add informations as userdata: easily accessible
				applyPostfix(process, 'patJets', postfix).userData.userFloats.src += ['pileupJetIdEvaluator%s:fullDiscriminant' % postfix]
				applyPostfix(process, 'patJets', postfix).userData.userInts.src += ['pileupJetIdEvaluator%s:cutbasedId' % postfix, 'pileupJetIdEvaluator%s:fullId' % postfix]
	
			# Quark / gluon discriminator
			# FIXME: Puppi needs some love
			# FIXME: So does SK
			if 'qg_tagger' in params and params['qg_tagger'] and pu_method != 'Puppi' and pu_method != 'SK':
	
				taggerPayload = 'QGL_%sPF%s' % (algo, pu_method.lower())
	
				setattr(process, 'QGTagger%s' % postfix, process.QGTagger.clone(
						srcJets = cms.InputTag(jetCollection),
						jetsLabel = cms.string(taggerPayload)
					))
	
				applyPostfix(process, "patJets", postfix).userData.userFloats.src += ['QGTagger%s:qgLikelihood' % postfix]
	  
	# Create METs from CHS and PUPPI
	from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
	
	## Gen MET
	
	### Copied from https://github.com/cms-sw/cmssw/blob/2b75137e278b50fc967f95929388d430ef64710b/RecoMET/Configuration/python/GenMETParticles_cff.py#L37
	process.genParticlesForMETAllVisible = cms.EDProducer(
			"InputGenJetsParticleSelector",
			src = cms.InputTag("packedGenParticles"),
			partonicFinalState = cms.bool(False),
			excludeResonances = cms.bool(False),
			excludeFromResonancePids = cms.vuint32(),
			tausAsJets = cms.bool(False),
	
			ignoreParticleIDs = cms.vuint32(
				1000022,
				1000012, 1000014, 1000016,
				2000012, 2000014, 2000016,
				1000039, 5100039,
				4000012, 4000014, 4000016,
				9900012, 9900014, 9900016,
				39, 12, 14, 16
				)
			)
	process.load('RecoMET.METProducers.genMetTrue_cfi')
	
	## Raw PF METs
	process.load('RecoMET.METProducers.PFMET_cfi')
	
	process.pfMet.src = cms.InputTag('packedPFCandidates')
	addMETCollection(process, labelName='patPFMet', metSource='pfMet') # RAW MET
	process.patPFMet.addGenMET = isMC
	
	process.pfMetCHS = process.pfMet.clone()
	process.pfMetCHS.src = cms.InputTag("chs")
	process.pfMetCHS.alias = cms.string('pfMetCHS')
	addMETCollection(process, labelName='patPFMetCHS', metSource='pfMetCHS') # RAW CHS MET
	process.patPFMetCHS.addGenMET = isMC
	
	process.pfMetPuppi = process.pfMet.clone()
	process.pfMetPuppi.src = cms.InputTag("puppi")
	process.pfMetPuppi.alias = cms.string('pfMetPuppi')
	addMETCollection(process, labelName='patPFMetPuppi', metSource='pfMetPuppi') # RAW puppi MET
	process.patPFMetPuppi.addGenMET = isMC
	
	## Type 1 corrections
	from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1
	from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1
	
	### Standard
	if not hasattr(process, 'ak4PFJets'):
		print("WARNING: No AK4 jets produced. Type 1 corrections for MET are not available.")
	else:
		process.corrPfMetType1 = corrPfMetType1.clone(
			src = 'ak4PFJets',
			jetCorrLabel = 'ak4PFL1FastL2L3Corrector',
			offsetCorrLabel = 'ak4PFL1FastjetCorrector'
		)
		process.pfMetT1 = pfMetT1.clone(
			src = 'pfMet',
			srcCorrections = [ cms.InputTag("corrPfMetType1","type1") ]
		)
	
		addMETCollection(process, labelName='patMET', metSource='pfMetT1') # T1 MET
		process.patMET.addGenMET = isMC
	
	### PUPPI
	if not hasattr(process, 'ak4PFJetsPuppi'):
		print("WARNING: No AK4 puppi jets produced. Type 1 corrections for puppi MET are not available.")
	else:
		process.corrPfMetType1Puppi = corrPfMetType1.clone(
			src = 'ak4PFJetsPuppi',
			jetCorrLabel = 'ak4PFCHSL2L3Corrector', #FIXME: Use PUPPI corrections when available?
		)
		del process.corrPfMetType1Puppi.offsetCorrLabel # no L1 for PUPPI jets
		process.pfMetT1Puppi = pfMetT1.clone(
			src = 'pfMetPuppi',
			srcCorrections = [ cms.InputTag("corrPfMetType1Puppi","type1") ]
		)
	
		addMETCollection(process, labelName='patMETPuppi', metSource='pfMetT1Puppi') # T1 puppi MET
		process.patMETPuppi.addGenMET = isMC
	
	### CHS
	if not hasattr(process, 'ak4PFJetsCHS'):
		print("WARNING: No AK4 CHS jets produced. Type 1 corrections for CHS MET are not available.")
	else:
		process.corrPfMetType1CHS = corrPfMetType1.clone(
			src = 'ak4PFJetsCHS',
			jetCorrLabel = 'ak4PFCHSL1FastL2L3Corrector',
			offsetCorrLabel = 'ak4PFCHSL1FastjetCorrector'
		)
		process.pfMetT1CHS = pfMetT1.clone(
			src = 'pfMetCHS',
			srcCorrections = [ cms.InputTag("corrPfMetType1CHS","type1") ]
		)
	
		addMETCollection(process, labelName='patMETCHS', metSource='pfMetT1CHS') # T1 CHS MET
		process.patMETCHS.addGenMET = isMC
	
	
	## Slimmed METs
	from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
	#### CaloMET is not available in MiniAOD
	del slimmedMETs.caloMET
	
	### Standard
	process.slimmedMETs = slimmedMETs.clone()
	if hasattr(process, 'patMET'):
		# Create MET from Type 1 PF collection
		process.patMET.addGenMET = isMC
		process.slimmedMETs.src = cms.InputTag("patMET")
		process.slimmedMETs.rawUncertainties = cms.InputTag("patPFMet") # only central value
	else:
		# Create MET from RAW PF collection
		process.patPFMet.addGenMET = isMC
		process.slimmedMETs.src = cms.InputTag("patPFMet")
		del process.slimmedMETs.rawUncertainties # not available
	
	clean_met_(process.slimmedMETs)
	
	### PUPPI
	process.slimmedMETsPuppi = slimmedMETs.clone()
	if hasattr(process, "patMETPuppi"):
		# Create MET from Type 1 PF collection
		process.patMETPuppi.addGenMET = isMC
		process.slimmedMETsPuppi.src = cms.InputTag("patMETPuppi")
		process.slimmedMETsPuppi.rawUncertainties = cms.InputTag("patPFMetPuppi") # only central value
	else:
		# Create MET from RAW PF collection
		process.patPFMetPuppi.addGenMET = isMC
		process.slimmedMETsPuppi.src = cms.InputTag("patPFMetPuppi")
		del process.slimmedMETsPuppi.rawUncertainties # not available
	
	clean_met_(process.slimmedMETsPuppi)
	
	### CHS
	process.slimmedMETsCHS = slimmedMETs.clone()
	if hasattr(process, "patMETCHS"):
		# Create MET from Type 1 PF collection
		process.patMETCHS.addGenMET = isMC
		process.slimmedMETsCHS.src = cms.InputTag("patMETCHS")
		process.slimmedMETsCHS.rawUncertainties = cms.InputTag("patPFMetCHS") # only central value
	else:
		# Create MET from RAW PF collection
		process.patPFMetCHS.addGenMET = isMC
		process.slimmedMETsCHS.src = cms.InputTag("patPFMetCHS")
		del process.slimmedMETsCHS.rawUncertainties # not available
	
	clean_met_(process.slimmedMETsCHS)
	
	
	process.seqjetMetProducerMC = cms.Sequence()
	process.seqjetMetProducerMCPath = cms.Path(process.seqjetMetProducerMC)
	


	

