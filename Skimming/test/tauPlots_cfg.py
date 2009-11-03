import FWCore.ParameterSet.Config as cms

process = cms.Process('Analysis')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
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
#    'file:/user/edelhoff/mcData/CMSSW_313/LM0-313-SUSYPAT-V00-04-07.root'
    'file:/user/edelhoff/mcData/CMSSW_314/LM0-314-SUSYPAT-V00-04-11.root'
#    'file:/user/edelhoff/mcData/CMSSW_314/RelValZTT_314_SUSYPAT_V00-04-12_1000ev.root'
)
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))


###########------ Tools -------#############
process.dump = cms.EDAnalyzer('EventContentAnalyzer')

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genTausWithHistory"),#tauGenParticles"),#genParticles"),#
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),    
    printStatus = cms.untracked.bool(True),
    printIndex = cms.untracked.bool(False),
    status = cms.untracked.vint32( 1, 2 )
  )   

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticles"),#hadronicGenTaus"),#
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.TFileService = cms.Service('TFileService', fileName = cms.string('tauPlots.root'))

###########------ Selection Paths -------#############
process.load('SuSyAachen.Skimming.electronSelection_cff')
process.load('SuSyAachen.Skimming.muonSelection_cff')
process.load('SuSyAachen.Skimming.tauSelection_cff')
process.load('SuSyAachen.Skimming.jetMETSelection_cff')

process.load('SuSyAachen.Skimming.genSelection_cff')

process.GenParticles = cms.Path(
    process.seqGenParticles
    * process.seqBasicGenParticles
    * process.seqMatchedParticles
)

process.Leptons = cms.Path(
    process.seqMuons
    + process.seqElectrons
    + process.seqTaus
    + process.seqSelectJetMET
#    + process.dump
    )
###########------ Histos -------#############
process.load('SuSyAachen.Histograms.leptonCounter_cff')

process.Histograms = cms.Path(
    process.seqLeptonCounter
#    + process.printTree
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
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep patTaus_*_*_*',
                                                                      'keep recoGenParticles_*_*_*') 
                               )
process.outpath = cms.EndPath(process.out)

process.schedule = cms.Schedule( process.Leptons, process.GenParticles,
                                 process.Histograms
, process.outpath 
)

