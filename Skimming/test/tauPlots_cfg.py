import FWCore.ParameterSet.Config as cms

process = cms.Process('Analysis')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
    )
process.MessageLogger = cms.Service('MessageLogger',
                                    tauPlots = cms.untracked.PSet(
    INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
    ),
    FwkReport = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    reportEvery = cms.untracked.int32(1000),
    limit = cms.untracked.int32(1000)
    ),
    default = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(100)
    ),
    Root_NoDictionary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(0)
    ),
    FwkJob = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(0)
    ),
    FwkSummary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(10000000)
    ),
    threshold = cms.untracked.string('INFO')
    ),
                                    destinations = cms.untracked.vstring('tauPlots', 'cout')
                                    )

process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(
    'file:/user/edelhoff/mcData/CMSSW_313/LM0-313-SUSYPAT-V00-04-07.root')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

###########------ Histos -------#############
import SuSyAachen.Histograms.leptonCounter_cfi 

#---- All reco particles
process.anyRecoLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "cleanLayer1Electrons",
  muonSource = "cleanLayer1Muons",
  tauSource = "cleanLayer1Taus",
) 

#---- isolated reco particles
process.isoRecoLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  electronSource = "isoElectrons",
  muonSource = "isoMuons",
  tauSource = "isoTaus",
) 

#---- All gen particles
process.anyGenLeptonExclusiveCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  method = 'exclusive',
  electronSource = "electronGenParticles",
  muonSource = "muonGenParticles",
  tauSource = "tauGenParticles",
)

process.anyGenLeptonInclusiveCounter = process.anyGenLeptonExclusiveCounter.clone(
  method = "inclusive"
)

#---- Cut Gen Particles
process.basicGenLeptonCounter = SuSyAachen.Histograms.leptonCounter_cfi.leptonCounter.clone(
  method = 'exclusive',
  electronSource = "electronBasicGenParticles",
  muonSource = "muonBasicGenParticles",
  tauSource = "tauBasicGenParticles",
)

process.seqLeptonCounter = cms.Sequence( 
  process.anyRecoLeptonCounter 
  + process.isoRecoLeptonCounter 
  + process.anyGenLeptonExclusiveCounter + process.anyGenLeptonInclusiveCounter 
  + process.basicGenLeptonCounter
)
###########------ Tools -------#############
process.dump = cms.EDAnalyzer('EventContentAnalyzer')

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("tauGenParticles"),#genParticles"),#genTausWithHistory"),#
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),    
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(False),
    status = cms.untracked.vint32( 1, 2,3 )
  )   

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("hadronicGenTaus"),#genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.TFileService = cms.Service('TFileService', fileName = cms.string('tauPlots.root'))

###########------ Selection Paths -------#############
process.load('SuSyAachen.Skimming.electronSelection_cff')
process.load('SuSyAachen.Skimming.muonSelection_cff')
process.load('SuSyAachen.Skimming.tauSelection_cff')
process.load('SuSyAachen.Skimming.genSelection_cff')

process.GenParticles = cms.Path(
    process.seqGenParticles
    + process.seqBasicGenParticles
)

process.Leptons = cms.Path(
    process.seqMuons
    + process.seqElectrons
    + process.seqTaus
#    + process.dump
#    + seqSelectJetMET
    )

process.Histograms = cms.Path(
    process.seqLeptonCounter
#    process.printTree
#    + process.printDecay
)

###########------ Output -------#############

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('tauPAT.root'),
                               # save only events passing the full path
                               #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('GenParticles') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('keep *' ) 
                               )
process.outpath = cms.EndPath(process.out)

process.schedule = cms.Schedule( process.GenParticles, process.Leptons,
                                 process.Histograms
#, process.outpath 
)

