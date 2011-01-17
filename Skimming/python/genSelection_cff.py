import FWCore.ParameterSet.Config as cms
import SuSyAachen.Skimming.patMatchedSelectors_cfi as patMatchers

filterGenElectrons = cms.bool(False)
filterGenMuons = cms.bool(False)
filterGenTaus = cms.bool(False)

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

# list of odgIds : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimGeneral/HepPDTESSource/data/pythiaparticle.tbl?view=markup
#1000000 to 3000000 are SuSy

#------------ Electrons
electronGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenElectrons,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep (pdgId = {e+} &status=1)| (pdgId = {e-}  &status=1)",
#    "drop (pdgId = {e+} & status = 2) | (pdgId = {e-} & status = 2)"
    )
)
import SuSyAachen.Skimming.electronSelection_cff as electronSelection
electronBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenElectrons,
   src = cms.InputTag("electronGenParticles"),
   cut = electronSelection.basicElectrons.cut
)

anyElectronMatchedParticles = patMatchers.patMatchedElectronSelector.clone( filter = filterGenElectrons,
)

basicElectronMatchedParticles = patMatchers.patMatchedElectronSelector.clone( filter = filterGenElectrons,
   src = cms.InputTag("basicElectrons"),
)

isoElectronMatchedParticles = patMatchers.patMatchedElectronSelector.clone( filter = filterGenElectrons,
   src = cms.InputTag("isoElectrons"),
)

electronMatchedGenParticles = cms.EDProducer('GenParticleElectronProducer',
    src = cms.InputTag("basicElectrons"),
)

#------------ Muons
muonGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenMuons,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep (pdgId = {mu+} & status=1) | (pdgId = {mu-} & status=1)",
#    "drop (pdgId = {mu+} & status = 2) | (pdgId = {mu-} & status = 2)"
    )
)

import SuSyAachen.Skimming.muonSelection_cff as muonSelection
muonBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenMuons,
   src = cms.InputTag("muonGenParticles"),
   cut = muonSelection.basicMuons.cut
)

anyMuonMatchedParticles = patMatchers.patMatchedMuonSelector.clone( filter = filterGenMuons,
)

basicMuonMatchedParticles = patMatchers.patMatchedMuonSelector.clone( filter = filterGenMuons,
   src = cms.InputTag("basicMuons"),
)

isoMuonMatchedParticles = patMatchers.patMatchedMuonSelector.clone( filter = filterGenMuons,
   src = cms.InputTag("isoMuons"),
)

muonMatchedGenParticles = cms.EDProducer('GenParticleMuonProducer',
    src = cms.InputTag("basicMuons"),
)

#------------ Taus
genTausWithHistory = cms.EDFilter( "GenParticlePruner", filter = filterGenTaus,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep++ pdgId = {tau+} | pdgId = {tau-}",
    #NOTE the taus will decay one needs to select the decay mode and then drop the daughters
    )
)

hadronicGenTaus = cms.EDFilter("GenDaughterExcluder", filter = filterGenTaus,
    src = cms.InputTag("genTausWithHistory"),
    daughterIds = cms.vint32( 11,-11, 13,-13) #e+ e- mu+ mu-
)

tauGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenTaus,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep pdgId = {tau+} | pdgId = {tau-}",
    "drop  status = 3"
    )
)

hadronicTauGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenTaus,
    src = cms.InputTag("hadronicGenTaus"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep pdgId = {tau+} | pdgId = {tau-}",
    "drop  status = 3"
    )
)

# require generated tau to decay hadronically
from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets
hadronicGenTauJets = cms.EDFilter("TauGenJetDecayModeSelector",
   src = cms.InputTag("tauGenJets"),
   select = cms.vstring('oneProng0Pi0', 'oneProng1Pi0', 'oneProng2Pi0', 'oneProngOther',
                        'threeProng0Pi0', 'threeProng1Pi0', 'threeProngOther', 'rare'),
   filter = cms.bool(False)
)


nutauGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenTaus,
  src = cms.InputTag("genParticles"),
  select = cms.vstring(
    "drop  *", # this is the default
    "keep ( pdgId = {nu_tau} & status = 1)  | (pdgId = {nu_taubar} & status = 1)",
  )
)

import SuSyAachen.Skimming.tauSelection_cff as tauSelection
tauBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenTaus,
   src = cms.InputTag("tauGenParticles"),
   cut = cms.string("abs(eta) <= 2.5 & pt >= 10")#tauSelection.basicTaus.cut                       
)

anyTauMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
)

anyTauJetMatchedParticles = patMatchers.patJetMatchedTauSelector.clone( filter = filterGenTaus,
)

basicTauJetMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
   src = cms.InputTag("basicTaus"),
)

isoTauJetMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
   src = cms.InputTag("isoTaus"),
)

susyCrossSectionProducer = cms.EDProducer( "SusyXSecProducer", src = cms.InputTag("genParticles"),
    susyVars = cms.VPSet(
    cms.PSet(var = cms.string("susyScanM0")),
    cms.PSet(var = cms.string("susyScanM12")),
    cms.PSet(var = cms.string("susyScantanbeta"))
    ),
    crossSections = cms.VPSet(
    cms.PSet(susyScanM0 = cms.double(270.), susyScanM12 = cms.double(450.), susyScantanbeta = cms.double(3.), susyScanA0 = cms.double(0.), susyScanMu = cms.double(-1.), ng = cms.double(0.224611), ns =  cms.double(0.0101205), nn = cms.double(53.1719), ll = cms.double(3.5998e-06), sb = cms.double(0.00592), ss = cms.double(0.0386), tb = cms.double(0.024632), bb = cms.double(0.001203), gg = cms.double(66.1), sg = cms.double(2.3)),
    )
)

# Produce PDF weights (maximum is 3)
pdfWeights = cms.EDProducer("PdfWeightProducer",
    # Fix POWHEG if buggy (this PDF set will also appear on output, 
    #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
    GenTag = cms.untracked.InputTag("genParticles"),
    PdfInfoTag = cms.untracked.InputTag("generator"),
    PdfSetNames = cms.untracked.vstring(
        "cteq66.LHgrid",
        #"MRST2006nnlo.LHgrid",
        #"MRST2007lomod.LHgrid"
    )
)

#------------ Sequences
seqGenParticles = cms.Sequence( electronGenParticles + muonGenParticles + tauGenParticles +nutauGenParticles + (tauGenJets*hadronicGenTauJets)
                                + (genTausWithHistory * hadronicGenTaus * hadronicTauGenParticles) )
seqBasicGenParticles = cms.Sequence( muonBasicGenParticles + electronBasicGenParticles + tauBasicGenParticles)

seqMatchedParticles = cms.Sequence( isoElectronMatchedParticles + isoMuonMatchedParticles + isoTauJetMatchedParticles 
                                    + anyElectronMatchedParticles + anyMuonMatchedParticles+ anyTauMatchedParticles +anyTauJetMatchedParticles
                                    + basicElectronMatchedParticles + basicMuonMatchedParticles+ basicTauJetMatchedParticles)
                                    
