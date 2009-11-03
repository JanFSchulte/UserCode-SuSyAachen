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
    src = cms.InputTag("hadronicGenTaus"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep pdgId = {tau+} | pdgId = {tau-}",
    "drop  status = 3"
    )
)

import SuSyAachen.Skimming.tauSelection_cff as tauSelection
tauBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenTaus,
   src = cms.InputTag("tauGenParticles"),
   cut = tauSelection.basicTaus.cut                       
)

anyTauMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
)

basicTauMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
   src = cms.InputTag("basicTaus"),
)

isoTauMatchedParticles = patMatchers.patMatchedTauSelector.clone( filter = filterGenTaus,
   src = cms.InputTag("isoTaus"),
)
#------------ Sequences
seqGenParticles = cms.Sequence( electronGenParticles + muonGenParticles + 
                                (genTausWithHistory * hadronicGenTaus * tauGenParticles) )
seqBasicGenParticles = cms.Sequence( muonBasicGenParticles + electronBasicGenParticles + tauBasicGenParticles)

seqMatchedParticles = cms.Sequence( isoElectronMatchedParticles + isoMuonMatchedParticles + isoTauMatchedParticles 
                                    + anyElectronMatchedParticles + anyMuonMatchedParticles+ anyTauMatchedParticles
                                    + basicElectronMatchedParticles + basicMuonMatchedParticles+ basicTauMatchedParticles)
                                    
