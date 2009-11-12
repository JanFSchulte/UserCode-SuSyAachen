import FWCore.ParameterSet.Config as cms

#-- PAT standard config -------------------------------------------------------
from PhysicsTools.PatAlgos.patSequences_cff import *

#-- Changes for electron and photon ID ----------------------------------------
# Turn off photon-electron cleaning (i.e., flag only)
cleanLayer1Photons.checkOverlaps.electrons.requireNoOverlaps = False
# Remove embedding of superClusters, will keep entire superCluster collection
allLayer1Electrons.embedSuperCluster = False
allLayer1Photons.embedSuperCluster   = False

#-- Tuning of Monte Carlo matching --------------------------------------------
# Also match with leptons of opposite charge
electronMatch.checkCharge = False
muonMatch.checkCharge     = False
tauMatch.checkCharge      = False

#-- Jet plus tracks -----------------------------------------------------------
from JetMETCorrections.Configuration.JetCorrectionsRecord_cfi import *

# ESSources delivering ZSP correctors
ZSPJetCorrectorIcone5 = cms.ESSource("ZSPJetCorrectionService",tagName=cms.string('ZSP_CMSSW219_Iterative_Cone_05'),label=cms.string('ZSPJetCorrectorIcone5'))
ZSPJetCorrectorSiscone5 = cms.ESSource("ZSPJetCorrectionService",tagName=cms.string('ZSP_CMSSW219_Iterative_Cone_05'),label=cms.string('ZSPJetCorrectorSiscone5'))
ZSPJetCorrectorAntiKt5 = cms.ESSource("ZSPJetCorrectionService",tagName=cms.string('ZSP_CMSSW219_Iterative_Cone_05'),label=cms.string('ZSPJetCorrectorAntiKt5'))
# EDProducers using ZSP correctors
ZSPJetCorJetIcone5 = cms.EDProducer("CaloJetCorrectionProducer",src=cms.InputTag("iterativeCone5CaloJets"),correctors=cms.vstring('ZSPJetCorrectorIcone5'),alias=cms.untracked.string('ZSPJetCorJetIcone5'))
ZSPJetCorJetSiscone5=cms.EDProducer( "CaloJetCorrectionProducer",src=cms.InputTag("sisCone5CaloJets"),correctors=cms.vstring('ZSPJetCorrectorSiscone5'),alias=cms.untracked.string('ZSPJetCorJetSiscone5'))
ZSPJetCorJetAntiKt5=cms.EDProducer( "CaloJetCorrectionProducer",src=cms.InputTag("antikt5CaloJets"),correctors=cms.vstring('ZSPJetCorrectorAntiKt5'),alias=cms.untracked.string('ZSPJetCorJetAntiKt5'))
# Jet-track association for IC5, SC5 and AK5
from RecoJets.Configuration.RecoJetAssociations_cff  import *
from RecoJets.JetAssociationProducers.iterativeCone5JTA_cff  import *

ZSPiterativeCone5JetTracksAssociatorAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone() 
ZSPiterativeCone5JetTracksAssociatorAtVertex.jets = cms.InputTag("ZSPJetCorJetIcone5")
ZSPiterativeCone5JetTracksAssociatorAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
ZSPiterativeCone5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ZSPJetCorJetIcone5")
ZSPiterativeCone5JetExtender = iterativeCone5JetExtender.clone() 
ZSPiterativeCone5JetExtender.jets = cms.InputTag("ZSPJetCorJetIcone5")
ZSPiterativeCone5JetExtender.jet2TracksAtCALO = cms.InputTag("ZSPiterativeCone5JetTracksAssociatorAtCaloFace")
ZSPiterativeCone5JetExtender.jet2TracksAtVX = cms.InputTag("ZSPiterativeCone5JetTracksAssociatorAtVertex")

ZSPSisCone5JetTracksAssociatorAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone()
ZSPSisCone5JetTracksAssociatorAtVertex.jets = cms.InputTag("ZSPJetCorJetSiscone5")
ZSPSisCone5JetTracksAssociatorAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
ZSPSisCone5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ZSPJetCorJetSiscone5")
ZSPSisCone5JetExtender = iterativeCone5JetExtender.clone()
ZSPSisCone5JetExtender.jets = cms.InputTag("ZSPJetCorJetSiscone5")
ZSPSisCone5JetExtender.jet2TracksAtCALO = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtCaloFace")
ZSPSisCone5JetExtender.jet2TracksAtVX = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtVertex")

ZSPAntiKt5JetTracksAssociatorAtVertex = iterativeCone5JetTracksAssociatorAtVertex.clone()
ZSPAntiKt5JetTracksAssociatorAtVertex.jets = cms.InputTag("ZSPJetCorJetAntiKt5")
ZSPAntiKt5JetTracksAssociatorAtCaloFace = iterativeCone5JetTracksAssociatorAtCaloFace.clone()
ZSPAntiKt5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ZSPJetCorJetAntiKt5")
ZSPAntiKt5JetExtender = iterativeCone5JetExtender.clone()
ZSPAntiKt5JetExtender.jets = cms.InputTag("ZSPJetCorJetAntiKt5")
ZSPAntiKt5JetExtender.jet2TracksAtCALO = cms.InputTag("ZSPAntiKt5JetTracksAssociatorAtCaloFace")
ZSPAntiKt5JetExtender.jet2TracksAtVX = cms.InputTag("ZSPAntiKt5JetTracksAssociatorAtVertex")

# ESSources delivering JPT correctors
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cfi import *
JetPlusTrackZSPCorrectorIcone5 = cms.ESSource( "JetPlusTrackCorrectionService", cms.PSet(JPTZSPCorrectorICone5), label = cms.string('JetPlusTrackZSPCorrectorIcone5'),)
JetPlusTrackZSPCorrectorIcone5.JetTrackCollectionAtVertex = cms.InputTag("ZSPiterativeCone5JetTracksAssociatorAtVertex")
JetPlusTrackZSPCorrectorIcone5.JetTrackCollectionAtCalo = cms.InputTag("ZSPiterativeCone5JetTracksAssociatorAtCaloFace")
JetPlusTrackZSPCorrectorIcone5.SplitMergeP = cms.int32(0)
JetPlusTrackZSPCorrectorIcone5.eIDValueMap = cms.InputTag("eidTight")

JetPlusTrackZSPCorrectorSiscone5 = cms.ESSource( "JetPlusTrackCorrectionService", cms.PSet(JPTZSPCorrectorICone5), label = cms.string('JetPlusTrackZSPCorrectorSiscone5'),)
JetPlusTrackZSPCorrectorSiscone5.JetTrackCollectionAtVertex = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtVertex")
JetPlusTrackZSPCorrectorSiscone5.JetTrackCollectionAtCalo = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtCaloFace")
JetPlusTrackZSPCorrectorSiscone5.SplitMergeP = cms.int32(1)
JetPlusTrackZSPCorrectorSiscone5.eIDValueMap = cms.InputTag("eidTight")

JetPlusTrackZSPCorrectorAntiKt5 = cms.ESSource( "JetPlusTrackCorrectionService", cms.PSet(JPTZSPCorrectorICone5), label = cms.string('JetPlusTrackZSPCorrectorAntiKt5'),)
JetPlusTrackZSPCorrectorAntiKt5.JetTrackCollectionAtVertex = cms.InputTag("ZSPAntiKt5JetTracksAssociatorAtVertex")
JetPlusTrackZSPCorrectorAntiKt5.JetTrackCollectionAtCalo = cms.InputTag("ZSPAntiKt5JetTracksAssociatorAtCaloFace")
JetPlusTrackZSPCorrectorAntiKt5.SplitMergeP = cms.int32(2)
JetPlusTrackZSPCorrectorAntiKt5.eIDValueMap = cms.InputTag("eidTight")

# EDProducers using JPT correctors
JetPlusTrackZSPCorJetIcone5 = cms.EDProducer( "CaloJetCorrectionProducer", src = cms.InputTag("ZSPJetCorJetIcone5"), correctors = cms.vstring('JetPlusTrackZSPCorrectorIcone5'), alias = cms.untracked.string('JetPlusTrackZSPCorJetIcone5'))
JetPlusTrackZSPCorJetSiscone5 = cms.EDProducer( "CaloJetCorrectionProducer", src = cms.InputTag("ZSPJetCorJetSiscone5"), correctors = cms.vstring('JetPlusTrackZSPCorrectorSiscone5'), alias = cms.untracked.string('JetPlusTrackZSPCorJetSiscone5'))
JetPlusTrackZSPCorJetAntiKt5 = cms.EDProducer( "CaloJetCorrectionProducer", src = cms.InputTag("ZSPJetCorJetAntiKt5"), correctors = cms.vstring('JetPlusTrackZSPCorrectorAntiKt5'), alias = cms.untracked.string('JetPlusTrackZSPCorJetAntiKt5'))

# Sequences
JetPlusTrackCorrectionsIcone5 = cms.Sequence(ZSPJetCorJetIcone5*ZSPiterativeCone5JetTracksAssociatorAtVertex*ZSPiterativeCone5JetTracksAssociatorAtCaloFace*ZSPiterativeCone5JetExtender*JetPlusTrackZSPCorJetIcone5)
JetPlusTrackCorrectionsSisCone5 = cms.Sequence(ZSPJetCorJetSiscone5*ZSPSisCone5JetTracksAssociatorAtVertex*ZSPSisCone5JetTracksAssociatorAtCaloFace*ZSPSisCone5JetExtender*JetPlusTrackZSPCorJetSiscone5)
JetPlusTrackCorrectionsAntiKt5 = cms.Sequence(ZSPJetCorJetAntiKt5*ZSPAntiKt5JetTracksAssociatorAtVertex*ZSPAntiKt5JetTracksAssociatorAtCaloFace*ZSPAntiKt5JetExtender*JetPlusTrackZSPCorJetAntiKt5)
JetPlusTrackCorrections = cms.Sequence(JetPlusTrackCorrectionsIcone5*JetPlusTrackCorrectionsSisCone5*JetPlusTrackCorrectionsAntiKt5)

#-- Extra Jet/MET collections -------------------------------------------------
# Add a few jet collections...
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.metTools import *

# Add latest HcalNoiseSummary
from RecoMET.METProducers.hcalnoiseinfoproducer_cfi import *
hcalnoise.refillRefVectors = True
hcalnoise.hcalNoiseRBXCollName = "hcalnoise" # This has changed in 33X

#-- Track Jets ----------------------------------------------------------------
# Select tracks for track jets
from PhysicsTools.RecoAlgos.TrackWithVertexSelector_cfi  import *
trackWithVertexSelector.src              = cms.InputTag("generalTracks")
trackWithVertexSelector.ptMax            = cms.double(500.0) 
trackWithVertexSelector.normalizedChi2   = cms.double(100.0)
trackWithVertexSelector.vertexTag        = cms.InputTag("offlinePrimaryVertices")
trackWithVertexSelector.copyTrajectories = cms.untracked.bool(False)
trackWithVertexSelector.vtxFallback      = cms.bool(False)
trackWithVertexSelector.useVtx           = cms.bool(False)
from SimGeneral.HepPDTESSource.pythiapdt_cfi  import *
tracksForJets = cms.EDProducer("ConcreteChargedCandidateProducer",
                                       src = cms.InputTag("trackWithVertexSelector"),
                                       particleType = cms.string('pi+')
                                       )
# Add jet collections
from RecoJets.JetProducers.SISConeJetParameters_cfi import *
from RecoJets.JetProducers.FastjetParameters_cfi import *
SISCone5TrackJets = cms.EDProducer("SISConeJetProducer",
                                           SISConeJetParameters,
                                           FastjetNoPU,
                                           src = cms.InputTag("tracksForJets"),
                                           jetType = cms.untracked.string('BasicJet'),
                                           alias = cms.untracked.string('SISCone5TrackJets'),
                                           coneRadius = cms.double(0.5),
                                           jetPtMin = cms.double(0.3),
                                           inputEMin = cms.double(0.0),
                                           inputEtMin = cms.double(0.2),
                                           )
addTrackJets = cms.Sequence(  trackWithVertexSelector
                                    * tracksForJets
                                    * SISCone5TrackJets )

from PhysicsTools.PatAlgos.tools.trigTools import *

def addSUSYCollections( process ):
    #CaloJets
    addJetCollection(process, cms.InputTag('antikt5CaloJets'),
                     'AK5',
                     doJTA            = True,
                     doBTagging       = True,
                     jetCorrLabel     = ('AK5','Calo'),
                     doType1MET       = True,
                     genJetCollection = cms.InputTag("antikt5GenJets")
                     )
    addJetCollection(process,cms.InputTag('sisCone5CaloJets'),
                     'SC5',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = ('SC5','Calo'),
                     doType1MET   = True,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection=cms.InputTag("sisCone5GenJets")
                     )
    # PF jets
    addJetCollection(process,cms.InputTag('antikt5PFJets'),
                     'AK5PF',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = True,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection=cms.InputTag("antikt5GenJets")
                     )
    addJetCollection(process,cms.InputTag('sisCone5PFJets'),
                     'SC5PF',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = True,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection=cms.InputTag("sisCone5GenJets")
                     )
	# JPT jets
    addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                     'AK5JPT',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection = cms.InputTag("antikt5GenJets")
                     )
    addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetSiscone5'),
                     'SC5JPT',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection = cms.InputTag("sisCone5GenJets")
                     )
    
    # Add tcMET and PFMET
    addTcMET(process,'TC')
    addPfMET(process,'PF')
    
    addJetCollection(process,cms.InputTag('SISCone5TrackJets'),
                     'SC5Track',
                     doJTA        = False,
                     doBTagging   = True,
                     jetCorrLabel = None,
                     doType1MET   = False,
                     doL1Cleaning = True,
                     doL1Counters = True,
                     genJetCollection = cms.InputTag("sisCone5GenJets")
                     )

    #-- Tune contents of jet collections  -----------------------------------------
    for jetName in ( '', 'AK5', 'SC5', 'AK5PF', 'SC5PF', 'AK5JPT', 'SC5JPT', 'SC5Track' ):
        module = getattr(process,'allLayer1Jets'+jetName)
        module.addTagInfos = False    # Remove tag infos
        module.addJetID    = True     # Add JetID variables
        module.embedGenJetMatch = False # Only keep reference, since we anyway keep the genJet collections

    #-- Trigger matching ----------------------------------------------------------
    switchOnTrigger( process )
    process.patTriggerSequence.remove( process.patTriggerMatcher )
    process.patTriggerEvent.patTriggerMatches  = ()

    # Rename default jet collection for uniformity
    process.cleanLayer1JetsIC5 = process.cleanLayer1Jets
    process.layer1METsIC5      = process.layer1METs
    
    # Modify subsequent modules
    process.cleanLayer1Hemispheres.patJets = process.cleanLayer1JetsIC5.label()
    process.countLayer1Jets.src            = process.cleanLayer1JetsIC5.label()
    
    # Modify counters' input
    process.allLayer1Summary.candidates.remove(cms.InputTag('layer1METs'))
    process.allLayer1Summary.candidates.append(cms.InputTag('layer1METsIC5'))
    process.cleanLayer1Summary.candidates.remove(cms.InputTag('cleanLayer1Jets'))
    process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1JetsIC5'))
    # Add new jet collections to counters (MET done automatically)
    for jets in ( 'AK5', 'SC5', 'AK5PF', 'SC5PF', 'AK5JPT', 'SC5JPT', 'SC5Track' ):
        process.allLayer1Summary.candidates.append(cms.InputTag('allLayer1Jets'+jets))
        process.selectedLayer1Summary.candidates.append(cms.InputTag('selectedLayer1Jets'+jets))
        process.cleanLayer1Summary.candidates.append(cms.InputTag('cleanLayer1Jets'+jets))
    
    # Full path
    process.seqSUSY_pattuple = cms.Sequence( process.hcalnoise*process.addTrackJets*process.JetPlusTrackCorrections
                                     * process.patDefaultSequence
                                     * process.patTrigger*process.patTriggerEvent )


def getSUSY_pattuple_outputCommands( name ):
    return [ # PAT Objects
        'keep *_cleanLayer1Photons_*_*',
        'keep *_cleanLayer1Electrons_*_*',
        'keep *_cleanLayer1Muons_*_*',
        'keep *_cleanLayer1Taus_*_*',
        'keep *_cleanLayer1Jets*_*_*',       # All Jets
        'keep *_layer1METs*_*_*',            # All METs
        'keep *_cleanLayer1Hemispheres_*_*',
        'keep *_cleanLayer1PFParticles_*_*',
        # Generator information
        'keep GenEventInfoProduct_generator_*_*',
        'keep GenRunInfoProduct_generator_*_*',
        # Generator particles/jets/MET
        'keep recoGenParticles_genParticles_*_*',
        'keep recoGenJets_iterativeCone5GenJets_*_*',
        'keep recoGenJets_sisCone5GenJets_*_*',
        'keep recoGenJets_antikt5GenJets_*_*',
        'keep recoGenMETs_*_*_*',
        # Trigger information
        'keep edmTriggerResults_TriggerResults_*_HLT',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep L1GlobalTriggerObjectMapRecord_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_*_*_*',
        # Others
        'keep *_muon*METValueMapProducer_*_*',   # Muon corrections to MET
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep *_towerMaker_*_*',                 # Keep CaloTowers for cross-cleaning
        'keep recoTracks_generalTracks_*_*',
        'keep recoSuperClusters_corrected*_*_*',
        'keep recoConversions_conversions_*_*',
        'keep recoTracks_*onversions_*_*',
        'keep HcalNoiseSummary_*_*_'+name, # Only keep the one we create
        'keep recoPFCandidates_particleFlow_*_*'
        ] 


