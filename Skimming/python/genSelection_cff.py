import FWCore.ParameterSet.Config as cms

filterGenParticles = cms.bool(False)

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

# list of odgIds : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimGeneral/HepPDTESSource/data/pythiaparticle.tbl?view=markup
#1000000 to 3000000 are SuSy

#------------ Electrons
electronGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenParticles,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep (pdgId = {e+} &status=1)| (pdgId = {e-}  &status=1)",
#    "drop (pdgId = {e+} & status = 2) | (pdgId = {e-} & status = 2)"
    )
)
import SuSyAachen.Skimming.electronSelection_cff as electronSelection
electronBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenParticles,
   src = cms.InputTag("electronGenParticles"),
   cut = electronSelection.basicElectrons.cut
)

#------------ Muons
muonGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenParticles,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep (pdgId = {mu+} & status=1) | (pdgId = {mu-} & status=1)",
#    "drop (pdgId = {mu+} & status = 2) | (pdgId = {mu-} & status = 2)"
    )
)

import SuSyAachen.Skimming.muonSelection_cff as muonSelection
muonBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenParticles,
   src = cms.InputTag("muonGenParticles"),
   cut = muonSelection.basicMuons.cut
)

#------------ Taus
genTausWithHistory = cms.EDFilter( "GenParticlePruner", filter = filterGenParticles,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep++ pdgId = {tau+} | pdgId = {tau-}",
    #NOTE the taus will decay one needs to select the decay mode and then drop the daughters
    )
)

hadronicGenTaus = cms.EDFilter("GenDaughterExcluder", filter = filterGenParticles,
    src = cms.InputTag("genTausWithHistory"),
    daughterIds = cms.vint32( 11,-11, 13,-13) #e+ e- mu+ mu-
)

tauGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenParticles,
    src = cms.InputTag("hadronicGenTaus"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep pdgId = {tau+} | pdgId = {tau-}",
    "drop  status = 3"
    )
)

import SuSyAachen.Skimming.tauSelection_cff as tauSelection
tauBasicGenParticles = cms.EDFilter("GenParticleSelector", filter = filterGenParticles,
   src = cms.InputTag("tauGenParticles"),
   cut = tauSelection.basicTaus.cut                       
)

#------------ Sequences
seqGenParticles = cms.Sequence( electronGenParticles + muonGenParticles + 
                                genTausWithHistory + hadronicGenTaus + tauGenParticles )

seqBasicGenParticles = cms.Sequence( muonBasicGenParticles + electronBasicGenParticles + tauBasicGenParticles)
