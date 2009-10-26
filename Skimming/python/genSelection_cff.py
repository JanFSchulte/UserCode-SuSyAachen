import FWCore.ParameterSet.Config as cms

filterGenParticles = cms.bool(False)

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

# list of odgIds : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimGeneral/HepPDTESSource/data/pythiaparticle.tbl?view=markup

#1000000 to 3000000
tauGenParticles = cms.EDFilter( "GenParticlePruner", filter = filterGenParticles,
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *", # this is the default
    "keep pdgId = {tau+} | pdgId = {tau-}",
    "drop (pdgId = {tau+} & status = 2) | (pdgId = {tau-} & status = 2)"
    )
)

seqGenParticles = cms.Sequence( tauGenParticles )
